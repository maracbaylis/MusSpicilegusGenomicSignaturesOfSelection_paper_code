from pathlib import Path
import os

ROOT = Path(os.environ.get("MSPIC_PROJECT_ROOT", ".")).resolve()
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib_cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import pandas as pd


BASE = ROOT / "expression_evolution_5species/results_dom_altmap_rerun"
BM_DIR = BASE / "bm_mus_complex_polytomy_tree"
OUTDIR = BM_DIR / "expression_profile_plots" / "selected_gene_profiles"
OUTDIR.mkdir(parents=True, exist_ok=True)

MEANS_FILE = BASE / "species_mean_voom_expression.csv"
RESID_FILE = BM_DIR / "bm_residuals_mus_complex_polytomy_all_species.csv"

SELECTED_GENES = ["Cdkl5", "Nectin3", "Dlgap1", "Sema3b", "Grn"]
SPECIES_ORDER = ["spretus", "spicilegus", "domesticus", "musculus", "castaneus"]
COMPACT_SPECIES_DISPLAY = [
    "M.\nspretus",
    "M.\nspicilegus",
    "M. m.\ndomesticus",
    "M. m.\nmusculus",
    "M. m.\ncastaneus",
]
SPECIES_COLORS = {
    "spicilegus": "#08519c",
    "spretus": "#5da5da",
    "domesticus": "#006d2c",
    "musculus": "#31a354",
    "castaneus": "#a1d99b",
}


def load_plot_data() -> pd.DataFrame:
    means = pd.read_csv(MEANS_FILE)
    resid = pd.read_csv(RESID_FILE)

    means["gene_name"] = means["gene_name"].astype(str).str.strip()
    resid["gene_name"] = resid["gene_name"].astype(str).str.strip()

    spic = resid[resid["species"] == "spicilegus"].copy()
    spic = spic.rename(
        columns={
            "observed": "spicilegus_observed",
            "predicted": "bm_predicted_spicilegus",
            "residual": "bm_residual_spicilegus",
            "z": "bm_z_spicilegus",
            "q_empirical": "bm_q_empirical_spicilegus",
            "direction": "bm_direction_spicilegus",
        }
    )

    plot_df = means.merge(
        spic[
            [
                "gene_name",
                "bm_predicted_spicilegus",
                "bm_residual_spicilegus",
                "bm_z_spicilegus",
                "bm_q_empirical_spicilegus",
                "bm_direction_spicilegus",
            ]
        ],
        on="gene_name",
        how="left",
    )
    plot_df = plot_df[plot_df["gene_name"].isin(SELECTED_GENES)].copy()
    plot_df["gene_name"] = pd.Categorical(plot_df["gene_name"], SELECTED_GENES, ordered=True)
    plot_df = plot_df.sort_values("gene_name")

    found = set(plot_df["gene_name"].astype(str))
    missing = [gene for gene in SELECTED_GENES if gene not in found]
    if missing:
        raise ValueError(f"Selected genes missing from input tables: {', '.join(missing)}")

    return plot_df


def draw_selected_profiles(plot_df: pd.DataFrame, show_q: bool, suffix: str) -> None:
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        plt.style.use("seaborn-whitegrid")

    fig, axes = plt.subplots(1, len(SELECTED_GENES), figsize=(20.5, 4.0), sharey=False)
    for idx, (ax, row) in enumerate(zip(axes, plot_df.itertuples(index=False))):
        obs = [getattr(row, species) for species in SPECIES_ORDER]
        spic_idx = SPECIES_ORDER.index("spicilegus")

        ax.plot(range(len(SPECIES_ORDER)), obs, color="#64748b", linewidth=1.4, alpha=0.9)
        for i, species in enumerate(SPECIES_ORDER):
            ax.scatter(i, getattr(row, species), s=54, color=SPECIES_COLORS[species], zorder=3)
        ax.scatter(
            spic_idx,
            row.bm_predicted_spicilegus,
            marker="D",
            s=60,
            color="#111827",
            zorder=4,
        )
        ax.vlines(
            spic_idx,
            ymin=min(row.bm_predicted_spicilegus, row.spicilegus),
            ymax=max(row.bm_predicted_spicilegus, row.spicilegus),
            color="#111827",
            linewidth=1.1,
            linestyle="--",
        )
        ax.set_xticks(range(len(SPECIES_ORDER)))
        ax.set_xticklabels(COMPACT_SPECIES_DISPLAY, rotation=0, fontsize=8, fontstyle="italic")
        ax.tick_params(axis="y", labelsize=9)
        ax.grid(False)
        ax.set_title(row.gene_name, fontweight="bold", fontsize=12)
        if show_q:
            ax.text(
                0.98,
                0.98,
                f"q={row.bm_q_empirical_spicilegus:.3g}",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=8,
                bbox=dict(boxstyle="round,pad=0.14", facecolor="white", edgecolor="none", alpha=0.75),
            )
        if idx == 0:
            ax.set_ylabel("Species mean normalized log-expression", fontsize=10)

    handles = [
        mlines.Line2D([], [], color="#64748b", linewidth=1.4, label="Observed species profile"),
        mlines.Line2D([], [], color="#111827", marker="D", linestyle="None", markersize=5, label="Predicted M. spicilegus"),
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, 1.03),
        fontsize=10,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    stem = f"selected_profiles_Cdkl5_Nectin3_Dlgap1_Sema3b_Grn{suffix}"
    fig.savefig(OUTDIR / f"{stem}.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUTDIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    plot_df = load_plot_data()
    plot_df.to_csv(OUTDIR / "selected_profiles_Cdkl5_Nectin3_Dlgap1_Sema3b_Grn_summary.csv", index=False)
    draw_selected_profiles(plot_df, show_q=False, suffix="_no_q")
    draw_selected_profiles(plot_df, show_q=True, suffix="_with_q")
    print(OUTDIR)
    print(len(plot_df))


if __name__ == "__main__":
    main()
