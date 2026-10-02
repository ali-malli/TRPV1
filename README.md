# TRPV1 Ensemble Docking Captures Ligand Pungency

This repository contains the workflow used to analyze molecular dynamics-derived TRPV1 conformations, perform ensemble docking, identify poses satisfying binding-site constraints, characterize ligand orientation, calculate population-weighted Boltzmann ensemble docking scores, and compare these scores with experimental pungency and TRPV1 EC50 data.

## Requirements

- Python 3
- NumPy
- pandas
- SciPy
- AutoDock Vina 1.2.7

Install the Python dependencies with:

```bash
pip install numpy pandas scipy
```

## Usage

Run the scripts from the repository root in the following order.

### 1. Run molecular docking

```bash
python run_vina.py
```

This generates the raw docking outputs in `docking_results/`.

### 2. Filter docking poses

```bash
python filter_docking.py
```

Output:

```text
results/filtered_docking_results.csv
```

### 3. Analyze ligand orientations

```bash
python orientation_analysis.py
```

Output:

```text
results/orientation_results.csv
```

### 4. Calculate ensemble docking scores

```bash
python calculate_ensemble_scores.py
```

Output:

```text
results/ensemble_docking_scores.csv
```

### 5. Calculate correlations with experimental data

```bash
python pungency_analysis.py
```

This reports correlations with ln(SHU) and −ln(EC50) measurements.

## Data

Input data are provided in `data/`, ligand structures in `ligands/`, receptor structures in `receptors/`, and processed outputs in `results/`.
The AutoDock Vina executable was obtained from: https://github.com/ccsb-scripps/AutoDock-Vina

## License

See `LICENSE`.

## Citation

If you use this repository or associated data, please cite the accompanying manuscript:

> Malli A. et al. *TRPV1 Ensemble Docking Captures Ligand Pungency.* Manuscript in preparation.

The manuscript citation will be updated upon publication.
