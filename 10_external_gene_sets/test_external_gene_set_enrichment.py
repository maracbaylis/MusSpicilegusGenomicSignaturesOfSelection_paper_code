#!/usr/bin/env python3
"""Fisher enrichment tests for external behavior/plasticity gene sets."""

from __future__ import annotations

from pathlib import Path
import os
import math
import re

import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "enrichment"
EXTERNAL = ROOT / "processed" / "external_gene_sets_long.csv"
CANDIDATES = Path(
    os.environ.get(
        "MSPIC_FINAL_HITS_WORKBOOK",
        Path(os.environ.get("MSPIC_INPUT_DIR", "data/inputs"))
        / "M. Spicilegus as foreground FDR significant results -- PAML, phyloP, BM  .xlsx",
    )
).resolve()
PHYLOP_BACKGROUND = Path(
    "phyloP_MSA_QC/go_enrichment_spicilegus_accel_173_p_accel_given_total_background/"
    "background_genes_from_valid_p_accel_subtree_given_total_cCREs.txt"
)
BM_BACKGROUND_GLOB = "bm_spicilegus_xlsx_results/*_bm_residuals_all_genes.csv"


def norm_gene(value):
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none"}:
        return None
    return re.sub(r"\s+", "", text).upper()


def bh_adjust(pvals):
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [None] * n
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        true_rank = n - rank + 1
        val = min(prev, pvals[idx] * n / true_rank)
        q[idx] = val
        prev = val
    return q


def log_choose(n, k):
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_logpmf(x, n_total, n_success, n_draws):
    return (
        log_choose(n_success, x)
        + log_choose(n_total - n_success, n_draws - x)
        - log_choose(n_total, n_draws)
    )


def fisher_exact_greater(a, b, c, d):
    """One-sided Fisher exact P(X >= observed) for [[a,b],[c,d]]."""
    n_total = a + b + c + d
    n_success = a + c
    n_draws = a + b
    max_x = min(n_success, n_draws)
    logs = [hypergeom_logpmf(x, n_total, n_success, n_draws) for x in range(a, max_x + 1)]
    max_log = max(logs)
    p = math.exp(max_log) * sum(math.exp(v - max_log) for v in logs)
    if b == 0 or c == 0:
        odds = math.inf if a * d > 0 else float("nan")
    else:
        odds = (a * d) / (b * c)
    return odds, min(1.0, p)


def read_candidate_hits():
    sheets = pd.read_excel(CANDIDATES, sheet_name=None)
    specs = {
        "BM_expression_outlier": ("M. spicilegus BM Hits", "Gene Name"),
        "PhyloP_accel_cCRE_linked": ("M. spicilegus PhyloP Hits", "Gene Name"),
        "PAML_positive_selection": ("M. spicilegus PAML Hits", "Gene Name "),
    }
    sets = {}
    display = {}
    for label, (sheet, col) in specs.items():
        genes = set()
        for value in sheets[sheet][col]:
            gene = norm_gene(value)
            if gene:
                genes.add(gene)
                display.setdefault(gene, str(value).strip())
        sets[label] = genes
    sets["PhyloP_plus_BM"] = sets["PhyloP_accel_cCRE_linked"] & sets["BM_expression_outlier"]
    return sets, display


def read_backgrounds():
    bm = set()
    for path in sorted(Path(".").glob(BM_BACKGROUND_GLOB)):
        df = pd.read_csv(path, usecols=["Genename"])
        bm |= {norm_gene(g) for g in df["Genename"] if norm_gene(g)}

    phylop = {
        norm_gene(line)
        for line in PHYLOP_BACKGROUND.read_text().splitlines()
        if norm_gene(line)
    }

    return {
        "BM_all_tissue_tested_genes": bm,
        "PhyloP_valid_cCRE_linked_gene_background": phylop,
        "BM_and_PhyloP_common_background": bm & phylop,
    }


def read_external_sets():
    external = pd.read_csv(EXTERNAL)
    external["gene_symbol_upper"] = external["gene_symbol_upper"].map(norm_gene)
    external = external.dropna(subset=["gene_symbol_upper"])
    external_sets = {}
    external_display = {}
    for (dataset, subset), df in external.groupby(["dataset", "subset"], sort=True):
        key = (dataset, subset)
        external_sets[key] = set(df["gene_symbol_upper"])
        external_display[key] = {
            row.gene_symbol_upper: row.gene_symbol
            for row in df[["gene_symbol_upper", "gene_symbol"]].itertuples(index=False)
        }
    return external_sets, external_display


def fisher_test(background, candidates, external):
    cand = candidates & background
    ext = external & background
    overlap = cand & ext
    a = len(overlap)
    b = len(cand - ext)
    c = len(ext - cand)
    d = len(background - cand - ext)
    odds, p = fisher_exact_greater(a, b, c, d)
    expected = len(cand) * len(ext) / len(background) if background else float("nan")
    return a, b, c, d, odds, p, expected, overlap, cand, ext


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    candidates, candidate_display = read_candidate_hits()
    backgrounds = read_backgrounds()
    external_sets, external_display = read_external_sets()

    tests = [
        ("BM_expression_outlier", "BM_all_tissue_tested_genes"),
        ("PhyloP_accel_cCRE_linked", "PhyloP_valid_cCRE_linked_gene_background"),
        ("PhyloP_plus_BM", "BM_and_PhyloP_common_background"),
    ]

    rows = []
    overlap_rows = []
    for candidate_label, background_label in tests:
        background = backgrounds[background_label]
        candidate_genes = candidates[candidate_label]
        for (dataset, subset), external_genes in external_sets.items():
            a, b, c, d, odds, p, expected, overlap, cand_bg, ext_bg = fisher_test(
                background, candidate_genes, external_genes
            )
            rows.append(
                {
                    "candidate_set": candidate_label,
                    "background": background_label,
                    "dataset": dataset,
                    "subset": subset,
                    "background_n": len(background),
                    "candidate_n_total": len(candidate_genes),
                    "candidate_n_in_background": len(cand_bg),
                    "external_n_total": len(external_genes),
                    "external_n_in_background": len(ext_bg),
                    "observed_overlap": a,
                    "expected_overlap": expected,
                    "odds_ratio": odds,
                    "fisher_p_greater": p,
                    "overlap_genes": ";".join(sorted(candidate_display.get(g, external_display[(dataset, subset)].get(g, g)) for g in overlap)),
                    "table_a_candidate_and_external": a,
                    "table_b_candidate_not_external": b,
                    "table_c_external_not_candidate": c,
                    "table_d_neither": d,
                }
            )
            for gene in sorted(overlap):
                overlap_rows.append(
                    {
                        "candidate_set": candidate_label,
                        "background": background_label,
                        "dataset": dataset,
                        "subset": subset,
                        "gene_symbol": candidate_display.get(
                            gene, external_display[(dataset, subset)].get(gene, gene)
                        ),
                        "gene_symbol_upper": gene,
                    }
                )

    results = pd.DataFrame(rows)
    results["fisher_q_bh_all_tests"] = bh_adjust(results["fisher_p_greater"].tolist())
    results["fisher_q_bh_within_candidate_set"] = pd.NA
    for candidate_label in results["candidate_set"].unique():
        mask = results["candidate_set"] == candidate_label
        results.loc[mask, "fisher_q_bh_within_candidate_set"] = bh_adjust(
            results.loc[mask, "fisher_p_greater"].tolist()
        )

    results = results.sort_values(
        ["fisher_q_bh_all_tests", "fisher_p_greater", "observed_overlap"],
        ascending=[True, True, False],
    )
    overlaps = pd.DataFrame(overlap_rows).sort_values(
        ["candidate_set", "dataset", "subset", "gene_symbol_upper"]
    )

    results.to_csv(OUT / "external_gene_set_enrichment_fisher.csv", index=False)
    overlaps.to_csv(OUT / "external_gene_set_enrichment_overlap_genes.csv", index=False)

    # A compact table of the rows most likely to matter in manuscript prose.
    keep = results[
        (results["fisher_q_bh_all_tests"] < 0.1)
        | (results["observed_overlap"] >= 5)
        | (results["candidate_set"] == "PhyloP_plus_BM")
    ].copy()
    keep.to_csv(OUT / "external_gene_set_enrichment_key_results.csv", index=False)

    (OUT / "external_gene_set_enrichment_README.md").write_text(
        "# External Gene-Set Enrichment Tests\n\n"
        "Generated by `external_gene_sets/test_external_gene_set_enrichment.py`.\n\n"
        "Fisher's exact tests use candidate-set-specific backgrounds:\n\n"
        "- BM: union of all genes tested in `bm_spicilegus_xlsx_results/*_bm_residuals_all_genes.csv` (34,037 genes).\n"
        "- PhyloP: valid cCRE-linked gene background from `phyloP_MSA_QC/.../background_genes_from_valid_p_accel_subtree_given_total_cCREs.txt` (31,134 genes).\n"
        "- PhyloP+BM: intersection of the BM and PhyloP backgrounds.\n\n"
        "PAML was not given a formal Fisher test here because the available PAML universe is protein-ID based, while external sets are gene-symbol based.\n",
        encoding="utf-8",
    )

    print(results.head(25).to_string(index=False))


if __name__ == "__main__":
    main()
