#!/usr/bin/env python3

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


SCORES_FILE = "ensemble_docking_scores.csv"
SHU_FILE = "pungency_data.csv"
EC50_FILE = "ec50_data.csv"


def metrics(x, y):
    r, p_pearson = pearsonr(x, y)
    rho, p_spearman = spearmanr(x, y)

    return {
        "n": len(x),
        "pearson_r": r,
        "pearson_p": p_pearson,
        "spearman_rho": rho,
        "spearman_p": p_spearman,
        "r2": r ** 2,
    }


def print_metrics(name, values):
    print(f"\n{name}")
    print("-" * len(name))

    for key, value in values.items():
        if key == "n":
            print(f"{key:15s}: {value}")
        else:
            print(f"{key:15s}: {value:.6f}")


def main():
    scores = pd.read_csv(SCORES_FILE)

    shu = pd.read_csv(SHU_FILE)
    shu["ln_shu"] = np.log(shu["shu"])

    shu_data = scores.merge(shu, on="ligand")
    print_metrics(
        "ln(SHU)",
        metrics(shu_data["ensemble_score"], shu_data["ln_shu"]),
    )

    ec50 = pd.read_csv(EC50_FILE)

    for assay, data in ec50.groupby("assay"):
        merged = scores.merge(data, on="ligand")

        print_metrics(
            f"-ln(EC50): {assay}",
            metrics(
                merged["ensemble_score"],
                merged["neg_ln_ec50"],
            ),
        )


if __name__ == "__main__":
    main()
