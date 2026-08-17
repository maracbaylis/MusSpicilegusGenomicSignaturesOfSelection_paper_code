#!/usr/bin/env python3
"""Compute empirical Enformer significance from matched background loci."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def bh_fdr(p_values: pd.Series) -> pd.Series:
    p = p_values.astype(float).to_numpy()
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    out = np.empty(n, dtype=float)
    out[order] = adjusted
    return pd.Series(out, index=p_values.index)


def empirical_upper_tail(observed: pd.Series, null_values: np.ndarray) -> tuple[pd.Series, pd.Series]:
    null_values = np.asarray(null_values, dtype=float)
    if len(null_values) == 0:
        nan = pd.Series(np.nan, index=observed.index)
        return nan, nan
    obs = observed.astype(float).to_numpy()
    p = np.array([(1 + np.sum(null_values >= x)) / (1 + len(null_values)) for x in obs])
    percentile = np.array([100 * np.mean(null_values <= x) for x in obs])
    return pd.Series(p, index=observed.index), pd.Series(percentile, index=observed.index)


def ccre_empirical_p(fg_path: Path, bg_path: Path) -> pd.DataFrame:
    fg = pd.read_csv(fg_path)
    bg = pd.read_csv(bg_path)
    metric = "all_mouse_tracks_abs_delta_sum"
    p, percentile = empirical_upper_tail(fg[metric], bg[metric].to_numpy())
    out = fg.copy()
    out["background_n"] = len(bg)
    out["empirical_p_upper_tail"] = p
    out["background_percentile"] = percentile
    out["empirical_q_bh"] = bh_fdr(out["empirical_p_upper_tail"])
    return out.sort_values(["empirical_q_bh", "empirical_p_upper_tail", metric], ascending=[True, True, False])


def target_empirical_p(fg_path: Path, bg_path: Path) -> pd.DataFrame:
    fg = pd.read_csv(fg_path)
    bg = pd.read_csv(bg_path)
    metric = "abs_delta_sum"
    rows = []
    for channel, fg_group in fg.groupby("mouse_output_channel", sort=False):
        bg_null = bg.loc[bg["mouse_output_channel"] == channel, metric].to_numpy()
        p, percentile = empirical_upper_tail(fg_group[metric], bg_null)
        tmp = fg_group.copy()
        tmp["background_n_same_channel"] = len(bg_null)
        tmp["empirical_p_upper_tail"] = p
        tmp["background_percentile_same_channel"] = percentile
        rows.append(tmp)
    out = pd.concat(rows, ignore_index=True)
    out["empirical_q_bh"] = bh_fdr(out["empirical_p_upper_tail"])
    return out.sort_values(["empirical_q_bh", "empirical_p_upper_tail", metric], ascending=[True, True, False])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--foreground-ccre", required=True, type=Path)
    parser.add_argument("--background-ccre", required=True, type=Path)
    parser.add_argument("--foreground-target", required=True, type=Path)
    parser.add_argument("--background-target", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    ccre = ccre_empirical_p(args.foreground_ccre, args.background_ccre)
    target = target_empirical_p(args.foreground_target, args.background_target)

    ccre_out = args.outdir / "foreground_vs_background_ccre_empirical_p.csv"
    target_out = args.outdir / "foreground_vs_background_strict_neural_target_empirical_p.csv"
    ccre.to_csv(ccre_out, index=False)
    target.to_csv(target_out, index=False)

    (args.outdir / "foreground_vs_background_summary.txt").write_text(
        "\n".join(
            [
                "Empirical Enformer significance summary",
                f"foreground_ccre_rows: {len(ccre)}",
                f"background_ccre_n: {int(ccre['background_n'].iloc[0]) if len(ccre) else 0}",
                f"foreground_target_rows: {len(target)}",
                f"target_tests_q_lt_0_10: {int((target['empirical_q_bh'] < 0.10).sum())}",
                f"target_tests_q_lt_0_05: {int((target['empirical_q_bh'] < 0.05).sum())}",
                f"ccre_tests_q_lt_0_10: {int((ccre['empirical_q_bh'] < 0.10).sum())}",
                f"ccre_tests_q_lt_0_05: {int((ccre['empirical_q_bh'] < 0.05).sum())}",
                "",
            ]
        )
    )

    print(f"Wrote {ccre_out}")
    print(f"Wrote {target_out}")


if __name__ == "__main__":
    main()
