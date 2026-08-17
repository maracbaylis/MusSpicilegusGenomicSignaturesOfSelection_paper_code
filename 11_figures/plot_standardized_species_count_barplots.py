#!/usr/bin/env python3
from pathlib import Path
import os

ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib_cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


SPECIES_ORDER = ["spretus", "spicilegus", "domesticus", "musculus", "castaneus"]
SPECIES_LABELS = {
    "spretus": "M. spretus",
    "spicilegus": "M. spicilegus",
    "domesticus": "M. m.\ndomesticus",
    "musculus": "M. m.\nmusculus",
    "castaneus": "M. m.\ncastaneus",
}
SPECIES_COLORS = {
    "spicilegus": "#08519c",
    "spretus": "#5da5da",
    "domesticus": "#006d2c",
    "musculus": "#31a354",
    "castaneus": "#a1d99b",
}

AXIS_LABEL_SIZE = 12
TICK_LABEL_SIZE = 11
ANNOTATION_SIZE = 10

FIG1_COUNTS = ROOT / "gene_full_sequences_qc_results/auto_qc_pass_barplot_counts_no_caroli_short_labels.csv"
FIG1_OUTDIR = ROOT / "gene_full_sequences_qc_results"

FIG5_COUNTS = (
    ROOT
    / "expression_evolution_5species/results_dom_altmap_rerun/"
    / "bm_mus_complex_polytomy_tree/bm_mus_complex_polytomy_summary_counts.csv"
)
FIG5_OUTDIR = FIG5_COUNTS.parent
FIG5_THRESHOLD = "q01_absz3_biggest_margin2"


def format_axis(ax, ylabel: str, ylim_pad: float = 0.12) -> None:
    ax.set_ylabel(ylabel, fontsize=AXIS_LABEL_SIZE)
    ax.tick_params(axis="y", labelsize=TICK_LABEL_SIZE, width=1.1, length=4)
    ax.tick_params(axis="x", labelsize=TICK_LABEL_SIZE, width=1.1, length=4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for spine in ax.spines.values():
        spine.set_linewidth(1.1)
        spine.set_color("#222222")
    ymax = max(p.get_height() for p in ax.patches)
    ax.set_ylim(0, ymax * (1 + ylim_pad))


def make_species_count_barplot(df: pd.DataFrame, value_col: str, ylabel: str, figsize: tuple[float, float]) -> plt.Figure:
    df = df.set_index("species").loc[SPECIES_ORDER].reset_index()
    x = range(len(df))

    fig, ax = plt.subplots(figsize=figsize)
    bars = ax.bar(
        x,
        df[value_col],
        width=0.72,
        color=[SPECIES_COLORS[species] for species in df["species"]],
        edgecolor=[SPECIES_COLORS[species] for species in df["species"]],
    )
    ax.set_xticks(list(x))
    ax.set_xticklabels([SPECIES_LABELS[species] for species in df["species"]], rotation=0, fontstyle="italic")
    format_axis(ax, ylabel)

    for bar, value in zip(bars, df[value_col]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + ax.get_ylim()[1] * 0.018,
            f"{int(value):,}",
            ha="center",
            va="bottom",
            fontsize=ANNOTATION_SIZE,
            color="#202020",
        )

    fig.tight_layout()
    return fig


def plot_fig1a() -> None:
    counts = pd.read_csv(FIG1_COUNTS).rename(columns={"group": "species", "auto_qc_pass": "count"})
    fig = make_species_count_barplot(
        counts[["species", "count"]],
        "count",
        "# of genes significant that are\ncandidates for positive selection",
        (8.2, 4.8),
    )
    for stem in [
        "auto_qc_pass_positive_selection_barplot_no_caroli_short_labels_tight_wrapped_no_panel",
        "fig1A_positive_selection_barplot_standardized_species_colors",
    ]:
        fig.savefig(FIG1_OUTDIR / f"{stem}.png", dpi=600)
        fig.savefig(FIG1_OUTDIR / f"{stem}.pdf")
    plt.close(fig)


def plot_fig5a() -> None:
    counts = pd.read_csv(FIG5_COUNTS)
    counts = counts[counts["threshold"] == FIG5_THRESHOLD].rename(columns={"gene_count": "count"})
    fig = make_species_count_barplot(
        counts[["species", "count"]],
        "count",
        "BM foreground outlier genes",
        (8.2, 4.8),
    )
    for stem in [
        "bm_mus_complex_polytomy_q01_absz3_biggest_margin2_counts_barplot",
        "bm_mus_complex_polytomy_q01_absz3_biggest_margin2_counts_barplot_navy",
        "fig5A_bm_expression_count_barplot_standardized_species_colors",
    ]:
        fig.savefig(FIG5_OUTDIR / f"{stem}.png", dpi=600)
        fig.savefig(FIG5_OUTDIR / f"{stem}.pdf")
    plt.close(fig)


def main() -> None:
    plot_fig1a()
    plot_fig5a()


if __name__ == "__main__":
    main()
