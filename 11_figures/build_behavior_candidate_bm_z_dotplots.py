from pathlib import Path
import os

ROOT = Path(os.environ.get("MSPIC_PROJECT_ROOT", ".")).resolve()
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib_cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE = ROOT / "expression_evolution_5species/results_dom_altmap_rerun"
BM_ALL_FILE = BASE / "bm_mus_complex_polytomy_tree/bm_residuals_mus_complex_polytomy_all_species.csv"
WORKBOOK = Path(os.environ.get("MSPIC_INPUT_DIR", "data/inputs")).resolve() / 'M. Spicilegus as foreground FDR significant results -- PAML, phyloP, BM  .xlsx'
OUTDIR = BASE / "behavior_candidate_summary"
OUTDIR.mkdir(parents=True, exist_ok=True)


PLOT_GROUPS = {
    "cooperative_mound_building": {
        "title": "Candidate expression shifts linked to cooperative mound-building behavior",
        "groups": {
            "Stress / HPA-axis regulation": ["Fkbp5", "Crhr1"],
            "Monoamine and arousal signaling": ["Htr1a", "Maoa", "Maob", "Adrb1"],
            "Excitation-inhibition balance and cognitive flexibility": ["Npas4", "Npas3", "Nrxn2", "Slc10a4", "Gabre"],
            "Neuropeptide / social-behavior signaling": ["Cck", "Tacr1"],
        },
    },
    "delayed_dispersal_puberty": {
        "title": "Candidate expression shifts linked to delayed dispersal and puberty",
        "groups": {
            "Reproductive-axis and steroid signaling": ["Kiss1", "Pttg1ip", "Prokr1", "Hsd17b2"],
            "Energy balance and growth state": ["Mc4r", "Thada", "Brs3"],
            "Circadian / seasonal timing": ["Per3"],
            "Neuropeptide processing and pituitary signaling": ["Pcsk1n", "Prlhr", "Cpe", "Chga"],
        },
    },
}


def load_finalized_source_membership() -> dict[str, set[str]]:
    membership: dict[str, set[str]] = {}
    sheet_to_source = {
        "M. spicilegus PAML Hits": "PAML",
        "M. spicilegus PhyloP Hits": "phyloP",
        "M. spicilegus BM Hits": "BM",
    }
    try:
        sheets = pd.read_excel(WORKBOOK, sheet_name=None)
        for sheet_name, source in sheet_to_source.items():
            df = sheets[sheet_name]
            gene_col = [col for col in df.columns if "Gene Name" in str(col)][0]
            for gene in df[gene_col].dropna().astype(str).str.strip():
                membership.setdefault(gene, set()).add(source)
    except ImportError:
        fallback = {
            "PAML": ["Gabre", "Per3", "Pcsk1n", "Slc17a6", "Tmem108", "Sytl4", "Sdc3", "Cspg5"],
            "phyloP": [
                "Adrb1",
                "Cdkl5",
                "Chl1",
                "Grin2d",
                "Hsd17b2",
                "Maob",
                "Mbd5",
                "Nectin3",
                "Npas3",
                "Npas4",
                "Nrxn2",
                "Prokr1",
                "Pttg1ip",
                "Rai1",
                "Rims1",
                "Sema6b",
                "Slc10a4",
                "Slc2a4",
                "Syndig1",
                "Tacr1",
            ],
            "BM": [
                "Abl2",
                "Begain",
                "Brs3",
                "Cck",
                "Cdkl5",
                "Chga",
                "Crhr1",
                "Disc1",
                "Dlgap1",
                "Fkbp5",
                "Htr1a",
                "Kiss1",
                "Maoa",
                "Mc4r",
                "Mbd5",
                "Nectin3",
                "Nrp2",
                "Per3",
                "Plxna4",
                "Plxnb3",
                "Prlhr",
                "Sema3a",
                "Sema3b",
                "Sema6c",
                "Shank3",
                "Shisa9",
                "Thada",
                "Wnt5a",
            ],
        }
        for source, genes in fallback.items():
            for gene in genes:
                membership.setdefault(gene, set()).add(source)
    return membership


def load_spic_bm_residuals() -> pd.DataFrame:
    df = pd.read_csv(BM_ALL_FILE)
    df = df[df["species"].eq("spicilegus")].copy()
    return df[
        [
            "gene_name",
            "observed",
            "predicted",
            "residual",
            "z",
            "p",
            "q",
            "direction",
            "is_biggest_abs_z",
            "abs_z_margin_over_next_species",
        ]
    ]


def build_plot_table(plot_key: str, membership: dict[str, set[str]], bm: pd.DataFrame) -> pd.DataFrame:
    groups = PLOT_GROUPS[plot_key]["groups"]
    genes = [gene for group_genes in groups.values() for gene in group_genes]
    bm_by_gene = bm.set_index("gene_name")
    rows = []
    for module, module_genes in groups.items():
        for gene in module_genes:
            row = {"candidate_module": module}
            if gene not in bm_by_gene.index:
                row.update(
                    {
                        "gene_name": gene,
                        "sources_in_finalized_workbook": ";".join(sorted(membership.get(gene, []))),
                        "present_in_bm_background": False,
                    }
                )
                rows.append(row)
                continue
            bm_row = bm_by_gene.loc[gene].to_dict()
            row.update(bm_row)
            row["gene_name"] = gene
            row["sources_in_finalized_workbook"] = ";".join(sorted(membership.get(gene, [])))
            row["present_in_bm_background"] = True
            row["is_strict_spic_bm_hit"] = bool(
                row["q"] < 0.01
                and abs(row["z"]) > 3
                and row["is_biggest_abs_z"]
                and row["abs_z_margin_over_next_species"] >= 2
            )
            rows.append(row)
    return pd.DataFrame(rows)


def plot_candidate_table(plot_key: str, table: pd.DataFrame) -> None:
    plot_df = table[table["present_in_bm_background"]].copy()
    groups = PLOT_GROUPS[plot_key]["groups"]
    order = [gene for genes in groups.values() for gene in genes if gene in set(plot_df["gene_name"])]
    plot_df["gene_name"] = pd.Categorical(plot_df["gene_name"], order, ordered=True)
    plot_df = plot_df.sort_values("gene_name").reset_index(drop=True)

    point_color = "#111827"
    line_color = "#6b7280"
    neutral = "#334155"
    band = "#f8fafc"

    y = np.arange(len(plot_df))
    fig_h = max(5.5, 0.42 * len(plot_df) + 2.0)
    fig, ax = plt.subplots(figsize=(8.6, fig_h))

    module_spans = []
    start = 0
    for module, genes in groups.items():
        found = [gene for gene in genes if gene in set(plot_df["gene_name"].astype(str))]
        if not found:
            continue
        end = start + len(found) - 1
        module_spans.append((module, start, end))
        start = end + 1

    for idx, (module, start, end) in enumerate(module_spans):
        if idx % 2 == 0:
            ax.axhspan(start - 0.5, end + 0.5, color=band, zorder=0)
        ax.text(
            min_x := 0,
            start - 0.35,
            module,
            fontsize=8.2,
            fontweight="bold",
            color="#475569",
            ha="left",
            va="bottom",
            transform=ax.get_yaxis_transform(),
        )
    ax.axvline(0, color=neutral, linewidth=1.0, zorder=1)

    for idx, row in plot_df.iterrows():
        ax.plot([0, row["z"]], [idx, idx], color=line_color, linewidth=1.5, alpha=0.55, zorder=2)
        ax.scatter(
            row["z"],
            idx,
            s=58,
            color=point_color,
            edgecolor="white",
            linewidth=0.7,
            zorder=3,
        )

    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["gene_name"], fontsize=9, fontstyle="italic")
    ax.invert_yaxis()
    min_x = min(-9.5, float(np.floor(plot_df["z"].min() - 1)))
    max_x = max(5.5, float(np.ceil(plot_df["z"].max() + 1)))
    ax.set_xlim(min_x, max_x)
    ax.set_xlabel(r"$\it{M.}$ spicilegus BM residual z-score", fontsize=11)
    ax.set_title(PLOT_GROUPS[plot_key]["title"], fontweight="bold", fontsize=12.5, pad=14)
    ax.text(
        min_x + (0 - min_x) * 0.52,
        -0.85,
        "Lower expression\nthan BM prediction",
        color=neutral,
        fontsize=9,
        ha="center",
        va="center",
        fontweight="bold",
    )
    ax.text(
        0 + (max_x - 0) * 0.52,
        -0.85,
        "Higher expression\nthan BM prediction",
        color=neutral,
        fontsize=9,
        ha="center",
        va="center",
        fontweight="bold",
    )
    ax.grid(axis="x", visible=False)
    ax.grid(axis="y", visible=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#475569")
    ax.spines["bottom"].set_color("#475569")
    ax.tick_params(axis="x", labelsize=9)
    ax.tick_params(axis="y", length=0)

    fig.tight_layout()
    fig.savefig(OUTDIR / f"{plot_key}_candidate_bm_z_dotplot.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / f"{plot_key}_candidate_bm_z_dotplot.pdf", bbox_inches="tight")
    fig.savefig(OUTDIR / f"{plot_key}_candidate_bm_z_dotplot_publication.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUTDIR / f"{plot_key}_candidate_bm_z_dotplot_publication.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    membership = load_finalized_source_membership()
    bm = load_spic_bm_residuals()
    combined = []
    for plot_key in PLOT_GROUPS:
        table = build_plot_table(plot_key, membership, bm)
        table.insert(0, "candidate_plot", plot_key)
        table.to_csv(OUTDIR / f"{plot_key}_candidate_bm_z_table.csv", index=False)
        plot_candidate_table(plot_key, table)
        combined.append(table)
    pd.concat(combined, ignore_index=True).to_csv(OUTDIR / "behavior_candidate_bm_z_tables_combined.csv", index=False)
    print(OUTDIR)


if __name__ == "__main__":
    main()
