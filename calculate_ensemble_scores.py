#!/usr/bin/env python3

import numpy as np
import pandas as pd


SCORES_FILE = "results/filtered_docking_results.csv"
POPULATIONS_FILE = "data/cluster_populations.csv"
OUTPUT_FILE = "results/ensemble_docking_scores.csv"

R = 0.0019872041  # kcal mol-1 K-1
T = 298.15


def boltzmann_average(scores):
    scores = np.asarray(scores, dtype=float)

    weights = np.exp(-scores / (R * T))
    return np.sum(scores * weights) / np.sum(weights)


def main():
    scores = pd.read_csv(SCORES_FILE)
    populations = pd.read_csv(POPULATIONS_FILE)

    population = populations.set_index("model")["population"]

    rows = []

    for ligand, ligand_data in scores.groupby("ligand"):
        ensemble_score = 0.0
        model_scores = {}

        for model, model_data in ligand_data.groupby("model"):
            score = boltzmann_average(model_data["score"])
            model_scores[model] = score
            ensemble_score += population[model] * score

        rows.append({
            "ligand": ligand,
            "model1_boltzmann": model_scores.get("TRPV1_model1", np.nan),
            "model2_boltzmann": model_scores.get("TRPV1_model2", np.nan),
            "model3_boltzmann": model_scores.get("TRPV1_model3", np.nan),
            "model4_boltzmann": model_scores.get("TRPV1_model4", np.nan),
            "ensemble_score": ensemble_score,
        })

    results = pd.DataFrame(rows)
    results.to_csv(OUTPUT_FILE, index=False)

    print(results.to_string(index=False))
    print(f"\nSaved to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
