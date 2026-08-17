#!/usr/bin/env python3
"""Plot left cerebral cortex and frontal cortex discovery rates side by side."""

from __future__ import annotations

import csv
import math
import os
from pathlib import Path


LEFT_CEREBRAL = Path("msa_qc_l_cerebral_cortex/results/msa_qc_discovery_rate_summary.csv")
FRONTAL = Path("frontal_cortex_fasta_qc_early_strict_filtered_percent_of_total.tsv")
OUT_PNG = Path("left_cerebral_and_frontal_cortex_discovery_rate.png")
OUT_PDF = Path("left_cerebral_and_frontal_cortex_discovery_rate.pdf")

os.environ.setdefault("MPLCONFIGDIR", str(Path(".plot_cache").resolve()))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SPECIES_ORDER = [
    "M. spretus",
    "M. spicilegus",
    "M. m. domesticus",
    "M. m. musculus",
    "M. m. castaneus",
]
SPECIES_COLORS = {
    "M. spicilegus": "#08519c",
    "M. spretus": "#5da5da",
    "M. m. musculus": "#31a354",
    "M. m. domesticus": "#006d2c",
    "M. m. castaneus": "#a1d99b",
}
SPECIES_MARKERS = {
    "M. spretus": "o",
    "M. spicilegus": "s",
    "M. m. domesticus": "^",
    "M. m. musculus": "D",
    "M. m. castaneus": "v",
}
LEFT_TIMEPOINTS = ["4d", "10d", "14d", "25d", "36d", "2m", "8-10m", "18-20m"]
FRONTAL_TIMEPOINTS = ["10d", "14d", "25d", "36d", "2m"]
TIMEPOINT_LABELS = {
    "4d": "4 days",
    "10d": "10 days",
    "14d": "14 days",
    "25d": "25 days",
    "36d": "36 days",
    "2m": "2 months",
    "8-10m": "8-10 months",
    "18-20m": "18-20 months",
}
YMAX = 0.18


def read_dicts(path: Path, delimiter: str) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def plot_panel(
    ax,
    lookup: dict[tuple[str, str], float],
    err_lookup: dict[tuple[str, str], float],
    data_timepoints: list[str],
    title: str,
    axis_timepoints: list[str] = LEFT_TIMEPOINTS,
) -> None:
    x_lookup = {timepoint: idx for idx, timepoint in enumerate(axis_timepoints)}
    for species in SPECIES_ORDER:
        x = [x_lookup[timepoint] for timepoint in data_timepoints]
        values = [lookup[(species, timepoint)] for timepoint in data_timepoints]
        errors = [err_lookup[(species, timepoint)] for timepoint in data_timepoints]
        ax.errorbar(
            x,
            values,
            yerr=errors,
            color=SPECIES_COLORS[species],
            marker=SPECIES_MARKERS[species],
            linewidth=2.2,
            markersize=5.0,
            markeredgecolor="white",
            markeredgewidth=0.7,
            elinewidth=0.9,
            capsize=2.5,
            capthick=0.9,
        )
    ax.set_xticks(list(range(len(axis_timepoints))))
    ax.set_xticklabels([TIMEPOINT_LABELS[timepoint] for timepoint in axis_timepoints], rotation=35, ha="right")
    ax.set_xlim(-0.3, len(axis_timepoints) - 0.7)
    ax.set_xlabel("Timepoint")
    ax.set_title(title)
    ax.set_ylim(0, YMAX)
    ax.margins(x=0.03)


def binomial_error_percent(pass_count: object, total_count: object) -> float:
    total = int(total_count)
    if total <= 0:
        return 0.0
    prop = int(pass_count) / total
    return 1.96 * math.sqrt(prop * (1.0 - prop) / total) * 100.0


def main() -> None:
    left_rows = read_dicts(LEFT_CEREBRAL, ",")
    frontal_rows = read_dicts(FRONTAL, "\t")
    left_lookup = {
        (row["species"], row["timepoint"]): float(row["discovery_rate_percent"])
        for row in left_rows
    }
    left_err_lookup = {
        (row["species"], row["timepoint"]): binomial_error_percent(
            row["qc_pass_significant_fastas"], row["total_tested"]
        )
        for row in left_rows
    }
    frontal_lookup = {
        (row["species"], row["timepoint"]): float(row["early_strict_pass_percent_total"])
        for row in frontal_rows
    }
    frontal_err_lookup = {
        (row["species"], row["timepoint"]): binomial_error_percent(row["EARLY_STRICT_PASS"], row["total_tested"])
        for row in frontal_rows
    }

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.titlesize": 12,
            "axes.labelsize": 10.5,
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    fig, (ax_left, ax_frontal) = plt.subplots(1, 2, figsize=(13.0, 5.4), sharey=True)
    plot_panel(ax_left, left_lookup, left_err_lookup, LEFT_TIMEPOINTS, "Left cerebral cortex")
    plot_panel(ax_frontal, frontal_lookup, frontal_err_lookup, FRONTAL_TIMEPOINTS, "Frontal cortex")
    ax_left.set_ylabel("Discovery rate of accelerated CREs (%)")
    fig.tight_layout(w_pad=2.5)
    fig.savefig(OUT_PNG, dpi=300)
    fig.savefig(OUT_PDF)
    plt.close(fig)
    print(OUT_PNG)
    print(OUT_PDF)


if __name__ == "__main__":
    main()
