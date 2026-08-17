#!/usr/bin/env python3
"""Run codeml likelihood-ratio tests and BH correction by foreground."""

from __future__ import annotations

import math
from pathlib import Path
import os

import pandas as pd


INPUT_DIR = Path(os.environ.get("MSPIC_INPUT_DIR", "data/inputs")).resolve()
OUTPUT_DIR = Path(os.environ.get("MSPIC_PROJECT_ROOT", ".")).resolve() / 'paml_lrt_bh_results'
NULL_FILE = INPUT_DIR / "PAML_dataframe_null_all.csv"
FOREGROUNDS = {
    "caroli": INPUT_DIR / "PAML_dataframe_model2_nonnull_caroli.csv",
    "wsbeij": INPUT_DIR / "PAML_dataframe_model2_nonnull_wsbeij.csv",
    "spretus": INPUT_DIR / "PAML_dataframe_model2_nonnull_spretus.csv",
    "spicilegus": INPUT_DIR / "PAML_dataframe_model2_nonnull_spicilegus.csv",
    "casteij": INPUT_DIR / "PAML_dataframe_model2_nonnull_casteij.csv",
    "pwkphj": INPUT_DIR / "PAML_dataframe_model2_nonnull_pwkphj.csv",
}
ALPHA = 0.05


def chi2_sf_df1(x: float) -> float:
    """Survival function for chi-square with 1 degree of freedom."""
    if pd.isna(x) or x <= 0:
        return 1.0
    return math.erfc(math.sqrt(x / 2.0))


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    """Return BH-adjusted q-values for a Series of p-values."""
    p = p_values.astype(float)
    q = pd.Series(index=p.index, dtype=float)
    valid = p.dropna().sort_values()
    n = len(valid)
    if n == 0:
        return q

    adjusted = valid * n / pd.Series(range(1, n + 1), index=valid.index)
    adjusted = adjusted.iloc[::-1].cummin().iloc[::-1].clip(upper=1.0)
    q.loc[adjusted.index] = adjusted
    return q


def main() -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)

    null_df = pd.read_csv(NULL_FILE).rename(
        columns={"lnL": "lnL_null", "Omega": "omega_null"}
    )

    summaries = []
    significant_frames = []

    for foreground, path in FOREGROUNDS.items():
        alt_df = pd.read_csv(path).rename(
            columns={
                "lnL": "lnL_alt",
                "Value1": "value1_alt",
                "Value2": "omega_foreground",
            }
        )
        merged = alt_df.merge(null_df, on="Filename", how="left", indicator=True)
        merged["foreground"] = foreground
        merged["matched_null"] = merged["_merge"].eq("both")
        merged["delta_lnL"] = merged["lnL_alt"] - merged["lnL_null"]
        merged["lrt"] = 2 * merged["delta_lnL"]
        merged.loc[merged["lrt"] < 0, "lrt"] = 0
        merged["p_value"] = merged["lrt"].map(chi2_sf_df1)
        merged.loc[
            ~merged["matched_null"],
            [
                "delta_lnL",
                "lrt",
                "p_value",
            ],
        ] = pd.NA
        merged["q_value_bh"] = benjamini_hochberg(merged["p_value"])
        merged["positive_selection_signal"] = (
            merged["matched_null"]
            & merged["q_value_bh"].le(ALPHA)
            & merged["omega_foreground"].gt(1)
            & merged["delta_lnL"].gt(0)
        )

        columns = [
            "foreground",
            "Filename",
            "lnL_null",
            "lnL_alt",
            "delta_lnL",
            "lrt",
            "p_value",
            "q_value_bh",
            "omega_null",
            "value1_alt",
            "omega_foreground",
            "matched_null",
            "positive_selection_signal",
        ]
        merged[columns].to_csv(OUTPUT_DIR / f"{foreground}_lrt_bh.csv", index=False)

        sig = merged.loc[merged["positive_selection_signal"], columns].copy()
        significant_frames.append(sig)
        summaries.append(
            {
                "foreground": foreground,
                "foreground_rows": len(alt_df),
                "matched_null_rows": int(merged["matched_null"].sum()),
                "unmatched_foreground_rows": int((~merged["matched_null"]).sum()),
                "q_lt_0.05_rows": int(merged["q_value_bh"].le(ALPHA).sum()),
                "positive_selection_signal_rows": len(sig),
                "omega_foreground_gt_1_rows": int(merged["omega_foreground"].gt(1).sum()),
            }
        )

    summary = pd.DataFrame(summaries)
    summary.to_csv(OUTPUT_DIR / "summary_counts.csv", index=False)

    all_significant = pd.concat(significant_frames, ignore_index=True)
    all_significant.sort_values(
        ["q_value_bh", "foreground", "Filename"],
        inplace=True,
        ignore_index=True,
    )
    all_significant.to_csv(
        OUTPUT_DIR / "all_significant_positive_selection.csv", index=False
    )

    print(summary.to_string(index=False))
    print()
    print(f"Wrote results to {OUTPUT_DIR}")
    print(
        "Combined positive-selection candidates:",
        OUTPUT_DIR / "all_significant_positive_selection.csv",
    )


if __name__ == "__main__":
    main()
