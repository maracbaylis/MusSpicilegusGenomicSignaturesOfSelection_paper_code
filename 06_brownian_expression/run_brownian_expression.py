#!/usr/bin/env python3
"""
Compute conditional Brownian-motion expression residuals for each focal species.

Input: CSV containing one row per gene and species-mean voom expression columns:
    gene_name,musculus,castaneus,domesticus,spicilegus,spretus

The covariance matrix below is the five-taxon matrix used in the supplied
analysis script, with tip order:
    musculus, castaneus, domesticus, spicilegus, spretus

For each gene and focal species, the script estimates the ancestral/intercept
mean and Brownian rate from the remaining species, computes the conditional
expectation and variance for the focal species, and reports a standardized
residual z-score.

Both analytic two-sided normal P-values and empirical P-values based on
200,000 draws from |N(0,1)| are reported. BH correction is applied within
each focal species.
"""

from __future__ import annotations

import argparse
import numpy as np
import pandas as pd
from scipy.stats import norm

TIP_ORDER = ["musculus", "castaneus", "domesticus", "spicilegus", "spretus"]

COVARIANCE = np.array(
    [
        [2.0, 1.5, 1.0, 0.0, 0.0],
        [1.5, 2.0, 1.0, 0.0, 0.0],
        [1.0, 1.0, 2.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 2.0, 1.0],
        [0.0, 0.0, 0.0, 1.0, 2.0],
    ],
    dtype=float,
)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("species_means_csv")
    p.add_argument("output_csv")
    p.add_argument("--gene-column", default="gene_name")
    p.add_argument("--empirical-draws", type=int, default=200000)
    p.add_argument("--seed", type=int, default=1)
    return p.parse_args()


def bh_adjust(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    m = len(ranked)
    q = np.empty(m, dtype=float)
    prev = 1.0
    for i in range(m - 1, -1, -1):
        prev = min(prev, ranked[i] * m / (i + 1))
        q[i] = prev
    out = np.empty(m, dtype=float)
    out[order] = q
    return out


def empirical_p_from_abs_z(zvals: np.ndarray, n_draws: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    null_abs_z = np.abs(rng.standard_normal(n_draws))
    null_abs_z.sort()
    observed = np.abs(np.asarray(zvals, dtype=float))
    exceed = n_draws - np.searchsorted(null_abs_z, observed, side="left")
    return (exceed + 1) / (n_draws + 1)


def conditional_bm(y: np.ndarray, focal_idx: int) -> dict[str, float]:
    other_idx = [i for i in range(len(y)) if i != focal_idx]

    y_other = y[other_idx]
    c22 = COVARIANCE[np.ix_(other_idx, other_idx)]
    c12 = COVARIANCE[np.ix_([focal_idx], other_idx)]
    c21 = COVARIANCE[np.ix_(other_idx, [focal_idx])]
    c11 = float(COVARIANCE[focal_idx, focal_idx])

    inv = np.linalg.inv(c22)
    one = np.ones((len(other_idx), 1))
    y_col = y_other.reshape(-1, 1)

    mu_hat = float((one.T @ inv @ y_col) / (one.T @ inv @ one))
    centered = y_other - mu_hat

    predicted = float(mu_hat + (c12 @ inv @ centered.reshape(-1, 1)).item())
    unit_conditional_variance = max(
        float(c11 - (c12 @ inv @ c21).item()),
        1e-12,
    )

    sigma2_hat = max(
        float(centered.T @ inv @ centered) / len(y_other),
        1e-12,
    )
    conditional_variance = unit_conditional_variance * sigma2_hat
    conditional_sd = float(np.sqrt(conditional_variance))

    residual = float(y[focal_idx] - predicted)
    z = residual / conditional_sd
    p_analytic = float(2.0 * norm.sf(abs(z)))

    return {
        "observed": float(y[focal_idx]),
        "predicted": predicted,
        "residual": residual,
        "unit_conditional_variance": unit_conditional_variance,
        "sigma2_hat_from_other_species": sigma2_hat,
        "conditional_variance": conditional_variance,
        "conditional_sd": conditional_sd,
        "z": z,
        "p_analytic": p_analytic,
    }


def main():
    args = parse_args()
    means = pd.read_csv(args.species_means_csv)

    required = [args.gene_column, *TIP_ORDER]
    missing = [col for col in required if col not in means.columns]
    if missing:
        raise SystemExit("Missing required columns: " + ", ".join(missing))

    rows = []
    for record in means[required].itertuples(index=False, name=None):
        gene_name = record[0]
        y = np.asarray(record[1:], dtype=float)
        if not np.all(np.isfinite(y)):
            continue

        for focal_idx, species in enumerate(TIP_ORDER):
            rows.append(
                {
                    "gene_name": gene_name,
                    "species": species,
                    **conditional_bm(y, focal_idx),
                }
            )

    out = pd.DataFrame(rows)

    out["q_analytic"] = (
        out.groupby("species", group_keys=False)["p_analytic"]
        .transform(lambda s: bh_adjust(s.to_numpy()))
    )

    # Use a distinct seed per species but the same number of null draws.
    empirical_parts = []
    for species_idx, species in enumerate(TIP_ORDER):
        part = out[out["species"] == species].copy()
        part["p_empirical"] = empirical_p_from_abs_z(
            part["z"].to_numpy(),
            n_draws=args.empirical_draws,
            seed=args.seed + species_idx,
        )
        part["q_empirical"] = bh_adjust(part["p_empirical"].to_numpy())
        empirical_parts.append(part)

    out = pd.concat(empirical_parts, ignore_index=True)
    out["direction"] = np.where(out["residual"] >= 0, "higher_than_expected", "lower_than_expected")
    out = out.sort_values(["species", "q_analytic", "gene_name"]).reset_index(drop=True)

    out.to_csv(args.output_csv, index=False)

    print(f"Genes analyzed: {out['gene_name'].nunique()}")
    print(f"Species: {', '.join(TIP_ORDER)}")
    print(f"Rows written: {len(out)}")
    print(f"Output: {args.output_csv}")


if __name__ == "__main__":
    main()
