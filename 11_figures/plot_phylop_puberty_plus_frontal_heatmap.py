from __future__ import annotations

from pathlib import Path
import os

ROOT = Path(os.environ.get("MSPIC_PROJECT_ROOT", ".")).resolve()
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib_cache"))
os.environ.setdefault("XDG_CACHE_HOME", str(ROOT / ".font_cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd


BASE = ROOT / "brain_region_workbook_phyloP_expression_summary"
PLOT_DIR = BASE / "allen_puberty_expression_energy_tests"
ALL_STRUCTURES = (
    BASE / "plots" / "workbook_phyloP_expression_positive_genes_Allen_all_structures_matrix.csv"
)
WORKBOOK_HITS = BASE / "source_workbook_spicilegus_phyloP_hits.csv"
EXISTING_MATRIX = (
    PLOT_DIR / "phyloP_puberty_circuits_plus_all_regions_all_genes_logscale_matrix.csv"
)

OUT_STEM = PLOT_DIR / "phyloP_puberty_plus_frontal_regions_all_genes_no_gene_labels_condensed_heatmap"
OUT_MATRIX = PLOT_DIR / "phyloP_puberty_plus_frontal_regions_all_genes_no_gene_labels_condensed_matrix.csv"

FRONTAL_ROWS = [
    "Prelimbic area",
    "Infralimbic area",
    "Anterior cingulate area",
    "Orbital area",
    "Frontal pole, cerebral cortex",
    "Secondary motor area",
]

NEGATIVE_CONTROL_ROWS = [
    "Lateral reticular nucleus, parvicellular part",
]


def workbook_gene_order() -> list[str]:
    hits = pd.read_csv(WORKBOOK_HITS)
    genes = (
        hits["Gene Name"]
        .dropna()
        .astype(str)
        .str.strip()
        .replace("", np.nan)
        .dropna()
        .drop_duplicates()
        .tolist()
    )
    return genes


def main() -> None:
    existing = pd.read_csv(EXISTING_MATRIX, index_col=0)
    all_structures = pd.read_csv(ALL_STRUCTURES, index_col=0)

    requested_rows = FRONTAL_ROWS + NEGATIVE_CONTROL_ROWS
    missing = [row for row in requested_rows if row not in all_structures.index]
    if missing:
        raise ValueError(f"Missing requested Allen structures: {missing}")

    genes = [gene for gene in workbook_gene_order() if gene in existing.columns]
    existing = existing.loc[:, genes]
    frontal = all_structures.loc[FRONTAL_ROWS, genes]
    controls = all_structures.loc[NEGATIVE_CONTROL_ROWS, genes]
    controls.index = [f"Negative control: {row}" for row in controls.index]
    combined = pd.concat([existing, frontal, controls], axis=0)
    combined.to_csv(OUT_MATRIX)

    values = combined.astype(float).values

    fig, ax = plt.subplots(figsize=(13.8, 7.75))
    im = ax.imshow(
        values,
        aspect="auto",
        interpolation="nearest",
        cmap="Blues",
        norm=Normalize(vmin=0, vmax=30),
    )

    ax.set_xticks([])
    ax.set_xlabel(
        r"Genes linked to accelerated evolution on $\it{cis}$-regulatory elements",
        fontsize=10,
        labelpad=9,
    )
    ax.set_yticks(np.arange(len(combined.index)))
    ax.set_yticklabels(combined.index, fontsize=8.0)

    for i, label in enumerate(ax.get_yticklabels()):
        if i < 5:
            label.set_fontweight("bold")

    for y in [4.5, len(existing.index) - 0.5, len(existing.index) + len(FRONTAL_ROWS) - 0.5]:
        ax.axhline(y, color="black", linewidth=0.65, alpha=0.45)

    ax.set_title(
        r"Allen Mouse Brain Atlas expression in 56 day old C57BL/6J males for genes linked to accelerated "
        "\n"
        r"$\it{cis}$-regulatory evolution in $\it{M.\ spicilegus}$",
        fontsize=11.5,
        fontweight="normal",
        pad=12,
    )

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis="y", length=0)

    cbar = fig.colorbar(im, ax=ax, fraction=0.022, pad=0.018)
    cbar.set_ticks([0, 5, 10, 20, 30])
    cbar.set_label("Allen ISH expression energy", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    fig.subplots_adjust(left=0.305, right=0.91, top=0.84, bottom=0.095)
    fig.savefig(OUT_STEM.with_suffix(".png"), dpi=450)
    fig.savefig(OUT_STEM.with_suffix(".pdf"))
    plt.close(fig)

    print(OUT_STEM.with_suffix(".png"))
    print(OUT_STEM.with_suffix(".pdf"))
    print(OUT_MATRIX)
    print(f"rows: {len(combined.index)} genes: {len(combined.columns)}")


if __name__ == "__main__":
    main()
