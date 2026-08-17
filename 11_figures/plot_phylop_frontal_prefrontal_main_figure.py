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
ALL_STRUCTURES = BASE / "plots" / "workbook_phyloP_expression_positive_genes_Allen_all_structures_matrix.csv"
WORKBOOK_HITS = BASE / "source_workbook_spicilegus_phyloP_hits.csv"
MAIN_MATRIX = PLOT_DIR / "phyloP_puberty_plus_frontal_regions_all_genes_no_gene_labels_condensed_matrix.csv"

OUT_STEM = PLOT_DIR / "phyloP_frontal_prefrontal_cortex_main_figure_original_order_heatmap"
OUT_MATRIX = PLOT_DIR / "phyloP_frontal_prefrontal_cortex_main_figure_original_order_matrix.csv"

ROWS = [
    "Prelimbic area",
    "Infralimbic area",
    "Anterior cingulate area",
    "Orbital area",
    "Frontal pole, cerebral cortex",
    "Secondary motor area",
    "Lateral reticular nucleus, parvicellular part",
]

LABELS = [
    "Prelimbic area",
    "Infralimbic area",
    "Anterior cingulate area",
    "Orbital area",
    "Frontal pole, cerebral cortex",
    "Secondary motor area",
    "Negative control:\nLateral reticular nucleus",
]


def workbook_gene_order() -> list[str]:
    hits = pd.read_csv(WORKBOOK_HITS)
    return (
        hits["Gene Name"]
        .dropna()
        .astype(str)
        .str.strip()
        .replace("", np.nan)
        .dropna()
        .drop_duplicates()
        .tolist()
    )


def main() -> None:
    all_structures = pd.read_csv(ALL_STRUCTURES, index_col=0)
    main_matrix = pd.read_csv(MAIN_MATRIX, index_col=0)
    missing = [row for row in ROWS if row not in all_structures.index]
    if missing:
        raise ValueError(f"Missing requested Allen structures: {missing}")

    genes = [
        gene
        for gene in workbook_gene_order()
        if gene in all_structures.columns and gene in main_matrix.columns
    ]
    matrix = all_structures.loc[ROWS, genes].copy()
    matrix.index = LABELS
    matrix.to_csv(OUT_MATRIX, index_label="row")

    values = np.log1p(matrix.astype(float).values)
    fig, ax = plt.subplots(figsize=(13.4, 4.1))
    im = ax.imshow(
        values,
        aspect="auto",
        interpolation="nearest",
        cmap="Blues",
        norm=Normalize(vmin=0, vmax=np.log1p(30)),
    )

    ax.set_xticks([])
    ax.set_xlabel(
        r"Genes linked to accelerated evolution on $\it{cis}$-regulatory elements",
        fontsize=10,
        labelpad=9,
    )
    ax.set_yticks(np.arange(len(matrix.index)))
    ax.set_yticklabels(matrix.index, fontsize=8.8)
    ax.axhline(5.5, color="black", linewidth=0.65, alpha=0.45)

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis="y", length=0)

    cbar = fig.colorbar(im, ax=ax, fraction=0.026, pad=0.018)
    raw_ticks = [0, 1, 3, 10, 30]
    cbar.set_ticks(np.log1p(raw_ticks))
    cbar.set_ticklabels([str(t) for t in raw_ticks])
    cbar.set_label("Allen ISH expression energy (log-scaled display)", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    fig.subplots_adjust(left=0.285, right=0.91, top=0.92, bottom=0.18)
    fig.savefig(OUT_STEM.with_suffix(".png"), dpi=450)
    fig.savefig(OUT_STEM.with_suffix(".pdf"))
    plt.close(fig)

    print(OUT_STEM.with_suffix(".png"))
    print(OUT_STEM.with_suffix(".pdf"))
    print(OUT_MATRIX)
    print(f"rows: {len(matrix.index)} genes: {len(matrix.columns)}")


if __name__ == "__main__":
    main()
