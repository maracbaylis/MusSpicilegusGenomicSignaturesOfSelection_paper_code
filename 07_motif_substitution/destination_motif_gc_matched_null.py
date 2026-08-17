#!/usr/bin/env python3
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from test_family_transition_gc_matched_null import assign_strata, bh_fdr


HERE = Path(__file__).resolve().parent
OUTDIR = HERE / "analysis"
IN_CONTEXT = OUTDIR / "spic_diff_site_context.tsv"

N_PERM = 20000
SEED = 43

OUT_TESTS = OUTDIR / "gc_matched_destination_motif_tests.tsv"
OUT_FDR_SITES = OUTDIR / "accelerated_CRE_destination_motif_GC_FDR_pass.tsv"
OUT_NOMINAL_SITES = OUTDIR / "accelerated_CRE_destination_motif_GC_nominal_pass.tsv"
OUT_SUMMARY = OUTDIR / "gc_matched_destination_motif_summary.tsv"
OUT_XLSX = OUTDIR / "destination_motif_GC_corrected_results.xlsx"


def genomic_position(cre, pos):
    parts = str(cre).split("__")
    if len(parts) != 2:
        return ""
    coord = parts[1].split("_")
    if len(coord) < 3:
        return ""
    chrom, start = coord[0], int(coord[1])
    return f"{chrom}:{start + int(pos) - 1}"


def permutation_all_categories(df, category_col):
    df = df.reset_index(drop=True).copy()
    categories = sorted(df[category_col].dropna().unique())
    cat_to_idx = {cat: i for i, cat in enumerate(categories)}
    cat_idx = df[category_col].map(cat_to_idx).to_numpy()
    selected = df["is_selected"].astype(int).to_numpy()
    observed = np.bincount(cat_idx[selected == 1], minlength=len(categories))

    rng = np.random.default_rng(SEED)
    counts = np.zeros((N_PERM, len(categories)), dtype=np.int16)
    for _, sub in df.groupby("stratum", sort=False):
        idx = sub.index.to_numpy()
        n_selected = int(sub["is_selected"].sum())
        if n_selected <= 0:
            continue
        if n_selected >= len(idx):
            fixed = np.bincount(cat_idx[idx], minlength=len(categories))
            counts += fixed
            continue
        for i in range(N_PERM):
            chosen = rng.choice(idx, size=n_selected, replace=False)
            counts[i] += np.bincount(cat_idx[chosen], minlength=len(categories)).astype(np.int16)

    expected = counts.mean(axis=0)
    sd = counts.std(axis=0, ddof=1)
    p_greater = (1 + (counts >= observed).sum(axis=0)) / (N_PERM + 1)
    p_less = (1 + (counts <= observed).sum(axis=0)) / (N_PERM + 1)
    p_two = np.minimum(1.0, 2 * np.minimum(p_greater, p_less))
    return categories, observed, expected, sd, p_greater, p_less, p_two


def site_sheet(switches, pass_tests):
    columns = [
        "cCRE",
        "Position",
        "Ancestral motif",
        "Spicilegus motif",
        "Ancestral family",
        "Spicilegus family",
        "Base change",
        "Ancestral motif score",
        "Spicilegus motif score",
        "matched_observed_over_expected",
        "matched_p_greater",
        "matched_p_two_sided",
        "matched_q_greater",
        "matched_q_two_sided",
        "selected_hit",
        "background_hit",
    ]
    if pass_tests.empty:
        return pd.DataFrame(columns=columns)

    selected = switches[switches["set"].eq("selected")].copy()
    selected["Position"] = [genomic_position(cre, pos) for cre, pos in zip(selected["CRE"], selected["pos"])]
    sites = selected.merge(pass_tests, left_on="spic_motif", right_on="destination_motif", how="inner")
    return pd.DataFrame(
        {
            "cCRE": sites["CRE"],
            "Position": sites["Position"],
            "Ancestral motif": sites["other_motif"],
            "Spicilegus motif": sites["spic_motif"],
            "Ancestral family": sites["other_family"],
            "Spicilegus family": sites["spic_family"],
            "Base change": sites["change"],
            "Ancestral motif score": sites["other_rel_score"],
            "Spicilegus motif score": sites["spic_rel_score"],
            "matched_observed_over_expected": sites["matched_observed_over_expected"],
            "matched_p_greater": sites["matched_p_greater"],
            "matched_p_two_sided": sites["matched_p_two_sided"],
            "matched_q_greater": sites["matched_q_greater"],
            "matched_q_two_sided": sites["matched_q_two_sided"],
            "selected_hit": sites["selected_hit"],
            "background_hit": sites["background_hit"],
        }
    ).sort_values(["matched_q_greater", "matched_p_greater", "Spicilegus motif", "cCRE", "Position"])


def main():
    switches = pd.read_csv(IN_CONTEXT, sep="\t")
    switches = switches[switches["motif_switch"] == True].copy()
    switches = assign_strata(switches)

    categories, observed, expected, sd, p_greater, p_less, p_two = permutation_all_categories(switches, "spic_motif")
    selected_total = int((switches["set"] == "selected").sum())
    background_total = int((switches["set"] == "background_matched").sum())

    rows = []
    for i, motif in enumerate(categories):
        selected_hit = int(observed[i])
        background_hit = int(((switches["set"] == "background_matched") & (switches["spic_motif"] == motif)).sum())
        _, p_unadj = fisher_exact(
            [[selected_hit, selected_total - selected_hit], [background_hit, background_total - background_hit]],
            alternative="greater",
        )
        or_unadj = ((selected_hit + 0.5) * (background_total - background_hit + 0.5)) / (
            (selected_total - selected_hit + 0.5) * (background_hit + 0.5)
        )
        family_values = switches.loc[switches["spic_motif"].eq(motif), "spic_family"].dropna().unique()
        rows.append(
            {
                "destination_motif": motif,
                "destination_family": ";".join(sorted(map(str, family_values))),
                "selected_hit": selected_hit,
                "selected_total": selected_total,
                "selected_fraction": selected_hit / selected_total,
                "background_hit": background_hit,
                "background_total": background_total,
                "background_fraction": background_hit / background_total,
                "total_hit": selected_hit + background_hit,
                "unadjusted_or": float(or_unadj),
                "unadjusted_p_greater": float(p_unadj),
                "matched_expected_selected_hit": float(expected[i]),
                "matched_expected_selected_fraction": float(expected[i] / selected_total),
                "matched_observed_over_expected": float(selected_hit / expected[i]) if expected[i] else np.nan,
                "matched_sd": float(sd[i]),
                "matched_p_greater": float(p_greater[i]),
                "matched_p_less": float(p_less[i]),
                "matched_p_two_sided": float(p_two[i]),
            }
        )

    tests = pd.DataFrame(rows)
    tests["matched_q_greater"] = bh_fdr(tests["matched_p_greater"])
    tests["matched_q_two_sided"] = bh_fdr(tests["matched_p_two_sided"])
    tests = tests.sort_values(["matched_q_greater", "matched_p_greater", "matched_observed_over_expected"], ascending=[True, True, False])
    tests.to_csv(OUT_TESTS, sep="\t", index=False)

    fdr_pass = tests[(tests["matched_q_greater"] < 0.05) | (tests["matched_q_two_sided"] < 0.05)].copy()
    nominal_pass = tests[(tests["matched_p_greater"] < 0.05) | (tests["matched_p_two_sided"] < 0.05)].copy()
    fdr_sites = site_sheet(switches, fdr_pass)
    nominal_sites = site_sheet(switches, nominal_pass)
    fdr_sites.to_csv(OUT_FDR_SITES, sep="\t", index=False)
    nominal_sites.to_csv(OUT_NOMINAL_SITES, sep="\t", index=False)

    summary = pd.DataFrame(
        [
            {"metric": "tested_destination_motifs", "value": int(tests.shape[0])},
            {"metric": "FDR_pass_destination_motifs", "value": int(fdr_pass.shape[0])},
            {"metric": "FDR_pass_accelerated_CRE_site_rows", "value": int(fdr_sites.shape[0])},
            {"metric": "nominal_matched_p_pass_destination_motifs", "value": int(nominal_pass.shape[0])},
            {"metric": "nominal_matched_p_pass_accelerated_CRE_site_rows", "value": int(nominal_sites.shape[0])},
            {"metric": "definition", "value": "Destination motif enrichment among motif-switch sites, GC/base-change matched by base class, local GC, and CRE GC."},
        ]
    )
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    try:
        with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as writer:
            summary.to_excel(writer, sheet_name="README", index=False)
            fdr_pass.to_excel(writer, sheet_name="GC_FDR_pass_dest_motifs", index=False)
            fdr_sites.to_excel(writer, sheet_name="GC_FDR_pass_CREs", index=False)
            nominal_pass.to_excel(writer, sheet_name="nominal_p_pass_dest_motifs", index=False)
            nominal_sites.to_excel(writer, sheet_name="nominal_p_pass_CREs", index=False)
            tests.to_excel(writer, sheet_name="all_destination_tests", index=False)
    except ModuleNotFoundError:
        pass

    print(f"Wrote: {OUT_TESTS.name}")
    print(f"Wrote: {OUT_SUMMARY.name}")
    print(summary.to_string(index=False))
    print("\nTop destination motif tests:")
    print(tests.head(25).to_string(index=False))


if __name__ == "__main__":
    main()
