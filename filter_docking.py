#!/usr/bin/env python3

from pathlib import Path
import re

import numpy as np
import pandas as pd


DOCKING_DIR = Path("docking_results")
RECEPTOR_DIR = Path("receptors")
OUTPUT_FILE = "results/filtered_docking_results.csv"

Y_RESID = 511
T_RESID = 550
DISTANCE_CUTOFF = 5.0
EXCLUDED_ATOMS = {"CA"}

SCORE_PATTERN = re.compile(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)")


def parse_atoms(lines):
    atoms = []

    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue

        try:
            atoms.append({
                "atom": line[12:16].strip(),
                "chain": line[21].strip(),
                "resid": int(line[22:26]),
                "coord": np.array([
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ]),
            })
        except ValueError:
            continue

    return atoms


def split_poses(lines):
    poses = []
    current = []

    for line in lines:
        if line.startswith("MODEL"):
            current = [line]
        elif line.startswith("ENDMDL"):
            current.append(line)
            poses.append(current)
            current = []
        elif current:
            current.append(line)

    return poses or [lines]


def pose_score(lines):
    for line in lines:
        match = SCORE_PATTERN.search(line)
        if match:
            return float(match.group(1))

    return None


def residue_coords(receptor, resid):
    with receptor.open() as f:
        atoms = parse_atoms(f)

    coords = {}

    for atom in atoms:
        if atom["resid"] == resid and atom["atom"] not in EXCLUDED_ATOMS:
            coords.setdefault(atom["chain"], []).append(atom["coord"])

    return {
        chain: np.array(xyz)
        for chain, xyz in coords.items()
    }


def min_distance(a, b):
    distances = a[:, None, :] - b[None, :, :]
    return float(np.sqrt(np.sum(distances**2, axis=2)).min())


def find_valid_pose(docking_file, y511, t550):
    with docking_file.open() as f:
        poses = split_poses(f.readlines())

    chains = sorted(set(y511) & set(t550))

    for pose_number, pose in enumerate(poses, start=1):
        score = pose_score(pose)
        atoms = parse_atoms(pose)

        if score is None or not atoms:
            continue

        ligand = np.array([atom["coord"] for atom in atoms])

        for chain in chains:
            d_y511 = min_distance(y511[chain], ligand)
            d_t550 = min_distance(t550[chain], ligand)

            if d_y511 <= DISTANCE_CUTOFF and d_t550 <= DISTANCE_CUTOFF:
                return {
                    "score": score,
                    "pose": pose_number,
                    "chain": chain,
                    "distance_y511": d_y511,
                    "distance_t550": d_t550,
                }

    return None


def main():
    rows = []

    for model_dir in sorted(DOCKING_DIR.iterdir()):
        if not model_dir.is_dir():
            continue

        receptor = RECEPTOR_DIR / f"{model_dir.name}.pdbqt"

        if not receptor.exists():
            print(f"Missing receptor: {receptor}")
            continue

        y511 = residue_coords(receptor, Y_RESID)
        t550 = residue_coords(receptor, T_RESID)

        for ligand_dir in sorted(model_dir.iterdir()):
            if not ligand_dir.is_dir():
                continue

            for seed_dir in sorted(ligand_dir.iterdir()):
                if not seed_dir.is_dir():
                    continue

                docking_files = sorted(seed_dir.glob("*.pdbqt"))

                for docking_file in docking_files:
                    result = find_valid_pose(docking_file, y511, t550)

                    if result is None:
                        continue

                    rows.append({
                        "model": model_dir.name,
                        "ligand": ligand_dir.name,
                        "seed": seed_dir.name.removeprefix("seed_"),
                        "score": result["score"],
                        "pose": result["pose"],
                        "chain": result["chain"],
                        "distance_y511": result["distance_y511"],
                        "distance_t550": result["distance_t550"],
                    })

                    break

    df = pd.DataFrame(rows)
    df.to_csv(OUTPUT_FILE, index=False)

    print(f"Saved {len(df)} valid seeds to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
