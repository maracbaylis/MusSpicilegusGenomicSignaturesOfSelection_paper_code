#!/usr/bin/env python3
"""
Apply the manuscript's stringent focal-lineage Brownian-expression filters.

Default candidate criteria:
    q < 0.01
    |z| >= 3
    focal |z| - next-most-extreme species |z| >= 2

By default the filter uses q_analytic because the supplied historical filter
script used its `q` field from the analytic Brownian residual table. Use
`--q-column q_empirical` to apply the same criteria to empirical BH q-values.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("residuals_csv")
    p.add_argument("output_csv")
    p.add_argument("--foreground", default="spicilegus")
    p.add_argument("--q-column", default="q_analytic",
                   choices=["q_analytic", "q_empirical"])
    p.add_argument("--q-threshold", type=float, default=0.01)
    p.add_argument("--abs-z-threshold", type=float, default=3.0)
    p.add_argument("--margin", type=float, default=2.0)
    return p.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.residuals_csv)

    required = {"gene_name", "species", "z", args.q_column}
    missing = sorted(required - set(df.columns))
    if missing:
        raise SystemExit("Missing required columns: " + ", ".join(missing))

    species = sorted(df["species"].dropna().unique())
    if args.foreground not in species:
        raise SystemExit(f"Foreground {args.foreground!r} not present in residual table")

    z_wide = df.pivot(index="gene_name", columns="species", values="z")
    q_wide = df.pivot(index="gene_name", columns="species", values=args.q_column)

    focal = pd.DataFrame(index=z_wide.index)
    focal["gene_name"] = focal.index
    focal["foreground"] = args.foreground
    focal["focal_q"] = q_wide[args.foreground]
    focal["focal_z"] = z_wide[args.foreground]
    focal["focal_abs_z"] = focal["focal_z"].abs()

    other_species = [s for s in species if s != args.foreground]
    other_abs = z_wide[other_species].abs()
    focal["next_most_extreme_abs_z"] = other_abs.max(axis=1)
    focal["next_most_extreme_species"] = other_abs.idxmax(axis=1)
    focal["abs_z_margin_over_next_species"] = (
        focal["focal_abs_z"] - focal["next_most_extreme_abs_z"]
    )

    keep = (
        (focal["focal_q"] < args.q_threshold)
        & (focal["focal_abs_z"] >= args.abs_z_threshold)
        & (focal["abs_z_margin_over_next_species"] >= args.margin)
    )

    hits = focal.loc[keep].copy()
    hits["direction"] = hits["focal_z"].map(
        lambda z: "higher_than_expected" if z >= 0 else "lower_than_expected"
    )
    hits["q_column_used"] = args.q_column
    hits = hits.sort_values(
        ["focal_q", "focal_abs_z", "abs_z_margin_over_next_species"],
        ascending=[True, False, False],
    ).reset_index(drop=True)

    Path(args.output_csv).parent.mkdir(parents=True, exist_ok=True)
    hits.to_csv(args.output_csv, index=False)

    print(f"Foreground: {args.foreground}")
    print(f"q column: {args.q_column}")
    print(f"Candidate genes: {len(hits)}")
    print(f"Output: {args.output_csv}")


if __name__ == "__main__":
    main()
