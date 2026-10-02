#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd


DOCKING_DIR = Path("docking_results")
RECEPTOR_DIR = Path("receptors")

INPUT_FILE = Path("filtered_docking_results.csv")
ATOM_FILE = Path("head_tail_atoms.csv")
OUTPUT_FILE = Path("orientation_results.csv")

Y_RESID = 511
T_RESID = 550
EXCLUDED_ATOMS = {"N", "CA", "C", "O"}


def atom_range(text):
    atoms = []

    for part in str(text).split(","):
        part = part.strip()

        if "-" in part:
            start, end = map(int, part.split("-"))
            atoms.extend(range(start, end + 1))
        else:
            atoms.append(int(part))

    return set(atoms)


def parse_atoms(lines):
    atoms = []

    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue

        try:
            atoms.append({
                "serial": int(line[6:11]),
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


def residue_centroid(receptor, resid, chain):
    with receptor.open() as f:
        atoms = parse_atoms(f)

    coords = [
        atom["coord"]
        for atom in atoms
        if (
            atom["resid"] == resid
            and atom["chain"] == chain
            and atom["atom"] not in EXCLUDED_ATOMS
        )
    ]

    if not coords:
        raise ValueError(f"Residue {resid} chain {chain} not found in {receptor}")

    return np.mean(coords, axis=0)


def ligand_centroid(atoms, serials):
    coords = [
        atom["coord"]
        for atom in atoms
        if atom["serial"] in serials
    ]

    if not coords:
        raise ValueError("No matching ligand atoms found")

    return np.mean(coords, axis=0)


def cosine_similarity(a, b):
    norm = np.linalg.norm(a) * np.linalg.norm(b)

    if norm == 0:
        return np.nan

    return float(np.dot(a, b) / norm)


def get_pose(docking_file, pose_number):
    with docking_file.open() as f:
        poses = split_poses(f.readlines())

    return poses[int(pose_number) - 1]


def main():
    results = pd.read_csv(INPUT_FILE)
    atom_table = pd.read_csv(ATOM_FILE)

    atom_map = {
        row.ligand: (
            atom_range(row.head_atoms),
            atom_range(row.tail_atoms),
        )
        for row in atom_table.itertuples()
    }

    rows = []

    for row in results.itertuples():
        if row.ligand not in atom_map:
            raise ValueError(f"No head/tail definition for {row.ligand}")

        receptor = RECEPTOR_DIR / f"{row.model}.pdbqt"
        docking_file = (
            DOCKING_DIR
            / row.model
            / row.ligand
            / f"seed_{row.seed}"
            / f"{row.ligand}_{row.seed}.pdbqt"
        )

        pose = get_pose(docking_file, row.pose)
        ligand_atoms = parse_atoms(pose)

        head_serials, tail_serials = atom_map[row.ligand]

        head = ligand_centroid(ligand_atoms, head_serials)
        tail = ligand_centroid(ligand_atoms, tail_serials)

        y511 = residue_centroid(receptor, Y_RESID, row.chain)
        t550 = residue_centroid(receptor, T_RESID, row.chain)

        ligand_vector = head - tail
        receptor_vector = t550 - y511

        cosine = cosine_similarity(ligand_vector, receptor_vector)

        rows.append({
            "model": row.model,
            "ligand": row.ligand,
            "seed": row.seed,
            "score": row.score,
            "chain": row.chain,
            "cosine": cosine,
            "orientation": "head-up" if cosine > 0 else "head-down",
        })

    output = pd.DataFrame(rows)
    output.to_csv(OUTPUT_FILE, index=False)

    print(f"Saved {len(output)} poses to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
