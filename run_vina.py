#!/usr/bin/env python3

from pathlib import Path
import os
import subprocess

VINA = "./vina_1.2.7_linux_x86_64"

RECEPTOR_DIR = Path("receptors")
LIGAND_DIR = Path("ligands")
OUTPUT_DIR = Path("docking_results")

CENTER = (83, 95, 95)
BOX_SIZE = (64, 64, 64)

EXHAUSTIVENESS = 8
NUM_MODES = 10


def run_docking(receptor, ligand, seed, cpus):
    receptor_name = receptor.stem
    ligand_name = ligand.stem

    out_dir = OUTPUT_DIR / receptor_name / ligand_name / f"seed_{seed}"
    out_dir.mkdir(parents=True, exist_ok=True)

    output_file = out_dir / f"{ligand_name}_{seed}.pdbqt"
    log_file = out_dir / f"{ligand_name}_{seed}.log"
    metadata_file = out_dir / "meta.tsv"

    command = [
        VINA,
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--center_x", str(CENTER[0]),
        "--center_y", str(CENTER[1]),
        "--center_z", str(CENTER[2]),
        "--size_x", str(BOX_SIZE[0]),
        "--size_y", str(BOX_SIZE[1]),
        "--size_z", str(BOX_SIZE[2]),
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", str(NUM_MODES),
        "--seed", str(seed),
        "--cpu", str(cpus),
        "--out", str(output_file),
    ]

    with log_file.open("w") as log:
        subprocess.run(
            command,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )

    write_metadata(
        metadata_file,
        seed,
        receptor,
        ligand,
        output_file,
        log_file,
    )


def write_metadata(metadata_file, seed, receptor, ligand, output_file, log_file):
    job_id = os.environ.get("SLURM_ARRAY_JOB_ID", "")
    task_id = os.environ.get("SLURM_ARRAY_TASK_ID", "")

    with metadata_file.open("w") as f:
        f.write(
            "seed\tjob_id\tarray_task_id\treceptor\tligand\tout_pdbqt\tlog_file\n"
        )
        f.write(
            f"{seed}\t{job_id}\t{task_id}\t"
            f"{receptor}\t{ligand}\t{output_file}\t{log_file}\n"
        )


def main():
    task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", "1"))
    seed = 100000 + task_id

    cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))

    receptors = sorted(RECEPTOR_DIR.glob("*.pdbqt"))
    ligands = sorted(LIGAND_DIR.glob("*.pdbqt"))

    if not receptors:
        raise FileNotFoundError(f"No receptor PDBQT files found in {RECEPTOR_DIR}")

    if not ligands:
        raise FileNotFoundError(f"No ligand PDBQT files found in {LIGAND_DIR}")

    for receptor in receptors:
        for ligand in ligands:
            run_docking(
                receptor=receptor,
                ligand=ligand,
                seed=seed,
                cpus=cpus,
            )


if __name__ == "__main__":
    main()
