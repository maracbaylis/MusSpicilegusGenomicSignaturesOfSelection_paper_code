#!/usr/bin/env python3
"""Compare M. spicilegus candidate genes to external behavior/plasticity gene sets."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import os
import re

import pandas as pd


ROOT = Path(__file__).resolve().parent
EXTERNAL = ROOT / "processed" / "external_gene_sets_long.csv"
OUT = ROOT / "overlaps"
DEFAULT_CANDIDATES = Path(
    os.environ.get(
        "MSPIC_FINAL_HITS_WORKBOOK",
        Path(os.environ.get("MSPIC_INPUT_DIR", "data/inputs"))
        / "M. Spicilegus as foreground FDR significant results -- PAML, phyloP, BM  .xlsx",
    )
).resolve()


def norm_gene(value):
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none"}:
        return None
    text = re.sub(r"\s+", "", text)
    return text


def find_column(df, candidates):
    normalized = {str(c).strip().lower(): c for c in df.columns}
    for candidate in candidates:
        key = candidate.strip().lower()
        if key in normalized:
            return normalized[key]
    raise KeyError(f"Could not find any of columns {candidates} in {list(df.columns)}")


def read_candidate_sets(path: Path):
    sheets = pd.read_excel(path, sheet_name=None)
    candidate_sets = {}
    gene_display = {}
    candidate_long_rows = []

    sheet_specs = [
        ("M. spicilegus PAML Hits", "PAML_positive_selection", ["Gene Name", "Gene Name "]),
        ("M. spicilegus PhyloP Hits", "PhyloP_accel_cCRE_linked", ["Gene Name", "gene_name"]),
        ("M. spicilegus BM Hits", "BM_expression_outlier", ["Gene Name", "gene_name"]),
    ]

    base_sets = {}
    for sheet_name, set_name, possible_gene_cols in sheet_specs:
        df = sheets[sheet_name]
        gene_col = find_column(df, possible_gene_cols)
        genes = []
        for gene in df[gene_col]:
            gene = norm_gene(gene)
            if not gene:
                continue
            upper = gene.upper()
            genes.append(upper)
            gene_display.setdefault(upper, gene)
            candidate_long_rows.append(
                {
                    "candidate_set": set_name,
                    "gene_symbol": gene,
                    "gene_symbol_upper": upper,
                    "source_sheet": sheet_name,
                }
            )
        base_sets[set_name] = set(genes)

    candidate_sets.update(base_sets)
    candidate_sets["Any_candidate_evidence"] = set().union(*base_sets.values())

    evidence_by_gene = defaultdict(set)
    for set_name, genes in base_sets.items():
        for gene in genes:
            evidence_by_gene[gene].add(set_name)

    candidate_sets["Multi_evidence_any_two_or_more"] = {
        gene for gene, evidence in evidence_by_gene.items() if len(evidence) >= 2
    }
    candidate_sets["PhyloP_plus_BM"] = (
        base_sets["PhyloP_accel_cCRE_linked"] & base_sets["BM_expression_outlier"]
    )
    candidate_sets["PAML_plus_PhyloP"] = (
        base_sets["PAML_positive_selection"] & base_sets["PhyloP_accel_cCRE_linked"]
    )
    candidate_sets["PAML_plus_BM"] = (
        base_sets["PAML_positive_selection"] & base_sets["BM_expression_outlier"]
    )
    candidate_sets["All_three_evidence_types"] = (
        base_sets["PAML_positive_selection"]
        & base_sets["PhyloP_accel_cCRE_linked"]
        & base_sets["BM_expression_outlier"]
    )

    return candidate_sets, gene_display, pd.DataFrame(candidate_long_rows)


def compare(candidate_sets, gene_display):
    external = pd.read_csv(EXTERNAL)
    external["gene_symbol_upper"] = external["gene_symbol_upper"].astype(str).str.upper()

    rows = []
    detail_rows = []
    for (dataset, subset), ext_df in external.groupby(["dataset", "subset"], sort=True):
        ext_genes = set(ext_df["gene_symbol_upper"].dropna())
        ext_display = {
            row.gene_symbol_upper: row.gene_symbol
            for row in ext_df[["gene_symbol_upper", "gene_symbol"]].itertuples(index=False)
        }
        source = ext_df["source"].iloc[0]
        for candidate_set, cand_genes in candidate_sets.items():
            overlap = sorted(ext_genes & cand_genes)
            rows.append(
                {
                    "dataset": dataset,
                    "subset": subset,
                    "candidate_set": candidate_set,
                    "external_n": len(ext_genes),
                    "candidate_n": len(cand_genes),
                    "overlap_n": len(overlap),
                    "pct_candidate_overlapped": len(overlap) / len(cand_genes) if cand_genes else 0,
                    "pct_external_overlapped": len(overlap) / len(ext_genes) if ext_genes else 0,
                    "jaccard": len(overlap) / len(ext_genes | cand_genes) if (ext_genes | cand_genes) else 0,
                    "overlap_genes": ";".join(gene_display.get(g, ext_display.get(g, g)) for g in overlap),
                    "source": source,
                }
            )
            for gene in overlap:
                detail_rows.append(
                    {
                        "dataset": dataset,
                        "subset": subset,
                        "candidate_set": candidate_set,
                        "gene_symbol_upper": gene,
                        "candidate_gene_symbol": gene_display.get(gene, gene),
                        "external_gene_symbol": ext_display.get(gene, gene),
                        "source": source,
                    }
                )

    summary = pd.DataFrame(rows).sort_values(
        ["overlap_n", "jaccard", "dataset", "subset", "candidate_set"],
        ascending=[False, False, True, True, True],
    )
    details = pd.DataFrame(detail_rows).sort_values(
        ["dataset", "subset", "candidate_set", "gene_symbol_upper"]
    )
    return summary, details


def write_outputs(summary, details, candidate_long):
    OUT.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUT / "candidate_external_overlap_summary.csv", index=False)
    details.to_csv(OUT / "candidate_external_overlap_genes_long.csv", index=False)
    candidate_long.to_csv(OUT / "spicilegus_candidate_sets_long.csv", index=False)

    wide = summary.pivot_table(
        index=["dataset", "subset"],
        columns="candidate_set",
        values="overlap_n",
        aggfunc="first",
        fill_value=0,
    ).reset_index()
    wide.to_csv(OUT / "candidate_external_overlap_count_matrix.csv", index=False)

    key_sets = [
        "Any_candidate_evidence",
        "Multi_evidence_any_two_or_more",
        "PhyloP_plus_BM",
        "PAML_positive_selection",
        "PhyloP_accel_cCRE_linked",
        "BM_expression_outlier",
    ]
    top = summary[summary["candidate_set"].isin(key_sets)].copy()
    top = top.sort_values(
        ["candidate_set", "overlap_n", "jaccard"],
        ascending=[True, False, False],
    )
    top.to_csv(OUT / "candidate_external_overlap_key_sets.csv", index=False)

    md = OUT / "candidate_external_overlap_README.md"
    md.write_text(
        "# Candidate External Gene-Set Overlaps\n\n"
        "Generated by `external_gene_sets/compare_candidate_overlaps.py`.\n\n"
        "Outputs:\n\n"
        "- `candidate_external_overlap_summary.csv`: one row per external set x candidate set, with counts, percentages, Jaccard, and overlap genes.\n"
        "- `candidate_external_overlap_genes_long.csv`: one row per overlapping gene.\n"
        "- `candidate_external_overlap_count_matrix.csv`: overlap counts as a wide matrix.\n"
        "- `candidate_external_overlap_key_sets.csv`: focused summary for the main candidate sets.\n"
        "- `spicilegus_candidate_sets_long.csv`: normalized candidate gene membership used for the comparison.\n\n"
        "Candidate sets include the three evidence classes from the M. spicilegus workbook plus union/intersection sets: any evidence, multi-evidence, PhyloP+BM, PAML+PhyloP, PAML+BM, and all three.\n",
        encoding="utf-8",
    )


def main():
    candidate_sets, gene_display, candidate_long = read_candidate_sets(DEFAULT_CANDIDATES)
    summary, details = compare(candidate_sets, gene_display)
    write_outputs(summary, details, candidate_long)
    print(f"Wrote overlap outputs to {OUT}")
    print(summary.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
