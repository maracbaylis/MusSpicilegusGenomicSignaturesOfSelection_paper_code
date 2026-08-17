#!/usr/bin/env python3
"""Plot mixed strict QC with spicilegus indel-aware evidence."""

import os
from pathlib import Path


CACHE_DIR = Path(".plot_cache")
CACHE_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(CACHE_DIR.resolve()))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


SUMMARY = Path("all_species_spic_indel_aware_fasta_qc_by_timepoint.tsv")
TOTALS = Path("frontal_cortex_foreground_totals.tsv")
OUT_TABLE = Path("all_species_spic_indel_aware_fasta_qc_percent_of_total.tsv")
OUT_PLOT = Path("all_species_spic_indel_aware_pass_percent_of_total_by_timepoint.png")

TIMEPOINTS = ["10d", "14d", "25d", "36d", "2m"]
SPECIES_LABELS = {
    "mus castaneous foreground": "M. m. castaneus",
    "mus domesticus foreground": "M. m. domesticus",
    "mus musculus foreground": "M. m. musculus",
    "spicilegus foreground": "M. spicilegus",
    "spretus foreground": "M. spretus",
}


def main() -> None:
    summary = pd.read_csv(SUMMARY, sep="\t")
    totals = pd.read_csv(TOTALS, sep="\t")
    df = summary.merge(totals[["foreground_folder", "timepoint", "total_tested"]], on=["foreground_folder", "timepoint"])
    df["species"] = df["foreground_folder"].map(SPECIES_LABELS)
    df["timepoint"] = pd.Categorical(df["timepoint"], categories=TIMEPOINTS, ordered=True)
    df["species"] = pd.Categorical(df["species"], categories=list(SPECIES_LABELS.values()), ordered=True)
    df = df.sort_values(["species", "timepoint"])
    df["spic_indel_aware_pass_percent_total"] = df["SPIC_INDEL_AWARE_PASS"] / df["total_tested"] * 100
    df.to_csv(OUT_TABLE, sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(10, 5.8))
    for species in SPECIES_LABELS.values():
        sub = df[df["species"] == species].set_index("timepoint").reindex(TIMEPOINTS)
        ax.plot(
            TIMEPOINTS,
            sub["spic_indel_aware_pass_percent_total"],
            marker="o",
            linewidth=2,
            label=species,
        )
    ax.set_xlabel("Timepoint")
    ax.set_ylabel("QC PASS (% of total tested)")
    ax.set_title("Frontal Cortex QC-confirmed Foreground Hits, Spicilegus Indel-aware")
    ax.grid(axis="y", color="#dddddd", linewidth=0.7)
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(OUT_PLOT, dpi=300)
    plt.close(fig)
    print(OUT_TABLE)
    print(OUT_PLOT)


if __name__ == "__main__":
    main()
