#!/usr/bin/env python3
from collections import Counter, OrderedDict
from pathlib import Path
import os
import re

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from scan_motif_family_switches_msa_sets import load_pssms, best_motif_overlap, clean_name, tf_family, find_taxon
from test_family_transition_gc_matched_null import assign_strata, bh_fdr


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SELECTED_DIR = Path(os.environ.get("MSPIC_INPUT_DIR", "data/inputs")).resolve() / 'selected_MSAs'
BACKGROUND_DIR = Path(os.environ.get("MSPIC_INPUT_DIR", "data/inputs")).resolve() / 'background_MSAs'
OUTDIR = HERE / "analysis"
INPUT_DIR = Path(os.environ.get("MSPIC_INPUT_DIR", "data/inputs")).resolve()
JASPAR_PATH = Path(
    os.environ.get("MSPIC_JASPAR_PATH", INPUT_DIR / "JASPAR.jaspar")
)

N_PERM = 20000
SEED = 41


def parse_fasta(path):
    seqs = OrderedDict()
    name = None
    parts = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(parts)
            name = line[1:]
            parts = []
        else:
            parts.append(line.strip())
    if name is not None:
        seqs[name] = "".join(parts)
    if seqs:
        max_len = max(len(seq) for seq in seqs.values())
        seqs = OrderedDict((name, seq + ("-" * (max_len - len(seq)))) for name, seq in seqs.items())
    return seqs


def spic_alignment_columns(seqs):
    return [i for i, base in enumerate(seqs["mus_spicilegus"]) if base != "-"]


def seq_on_spic_coords(seq, cols):
    out = []
    for col in cols:
        base = seq[col]
        out.append(base if base in "ACGT" else "N")
    return "".join(out)


def consensus_on_spic_coords(seqs, cols):
    others = [name for name in seqs if name != "mus_spicilegus"]
    consensus = []
    support = []
    for col in cols:
        bases = [seqs[name][col] for name in others if seqs[name][col] in "ACGT"]
        if not bases:
            consensus.append("N")
            support.append(0)
            continue
        base, n = Counter(bases).most_common(1)[0]
        consensus.append(base)
        support.append(n)
    return "".join(consensus), support


def gc_fraction(seq):
    bases = [b for b in seq if b in "ACGT"]
    if not bases:
        return float("nan")
    return sum(b in "GC" for b in bases) / len(bases)


def classify_change(change):
    ref, alt = change.split(">")
    if ref in "AT" and alt in "GC":
        return "AT->GC"
    if ref in "GC" and alt in "AT":
        return "GC->AT"
    return "other"


def genomic_position(cre, pos):
    match = re.match(r"(.+)__(chr[^_]+)_(\d+)_(\d+)$", str(cre))
    if not match:
        return ""
    _, chrom, start, _ = match.groups()
    return f"{chrom}:{int(start) + int(pos) - 1}"


def scan_motifs():
    pssms = load_pssms(JASPAR_PATH)
    selected_rows, selected_metrics = analyze_set_robust("selected", SELECTED_DIR, pssms)
    background_rows, background_metrics = analyze_set_robust("background_matched", BACKGROUND_DIR, pssms)

    selected_rows = selected_rows.sort_values(["family_switch", "delta_rel_score", "spic_rel_score"], ascending=[False, False, False])
    background_rows = background_rows.sort_values(["family_switch", "delta_rel_score", "spic_rel_score"], ascending=[False, False, False])
    selected_rows.to_csv(OUTDIR / "selected_motif_switch_sites.tsv", sep="\t", index=False)
    background_rows.to_csv(OUTDIR / "background_matched_motif_switch_sites.tsv", sep="\t", index=False)

    summary = pd.DataFrame([selected_metrics, background_metrics])
    summary["motif_switch_fraction_of_diff_sites"] = summary["motif_switch_sites"] / summary["spic_diff_sites"]
    summary["family_switch_fraction_of_diff_sites"] = summary["family_switch_sites"] / summary["spic_diff_sites"]
    summary["family_switch_fraction_of_motif_switches"] = summary["family_switch_sites"] / summary["motif_switch_sites"]
    summary.to_csv(OUTDIR / "motif_switch_set_summary.tsv", sep="\t", index=False)
    return selected_rows, background_rows, summary


def analyze_set_robust(set_name, msa_dir, pssms):
    rows = []
    metrics = {
        "set": set_name,
        "n_msas": 0,
        "spic_diff_sites": 0,
        "motif_switch_sites": 0,
        "family_switch_sites": 0,
        "same_family_switch_sites": 0,
        "unknown_family_switch_sites": 0,
    }

    for fa_path in sorted(msa_dir.glob("*.fa")):
        metrics["n_msas"] += 1
        cre = fa_path.stem
        seqs = parse_fasta(fa_path)
        cols = spic_alignment_columns(seqs)
        spic_seq = seq_on_spic_coords(seqs["mus_spicilegus"], cols)
        consensus_seq, support = consensus_on_spic_coords(seqs, cols)
        spret_name = find_taxon(seqs, "mus_spretus")
        spret_seq = seq_on_spic_coords(seqs[spret_name], cols) if spret_name else ("N" * len(spic_seq))

        for pos, (spic_base, consensus_base, support_n) in enumerate(zip(spic_seq, consensus_seq, support), start=1):
            if spic_base not in "ACGT" or consensus_base not in "ACGT":
                continue
            if spic_base == consensus_base:
                continue

            metrics["spic_diff_sites"] += 1
            spic_hit = best_motif_overlap(spic_seq, pos, pssms)
            consensus_hit = best_motif_overlap(consensus_seq, pos, pssms)
            if spic_hit is None or consensus_hit is None:
                continue
            if spic_hit["motif_name"] == consensus_hit["motif_name"]:
                continue
            if max(spic_hit["rel_score"], consensus_hit["rel_score"]) < 0.85:
                continue

            metrics["motif_switch_sites"] += 1
            spic_family = tf_family(spic_hit["motif_name"])
            other_family = tf_family(consensus_hit["motif_name"])
            family_switch = spic_family != other_family
            if family_switch:
                metrics["family_switch_sites"] += 1
            else:
                metrics["same_family_switch_sites"] += 1
            if "unknown" in {spic_family, other_family}:
                metrics["unknown_family_switch_sites"] += 1

            rows.append(
                {
                    "set": set_name,
                    "CRE": cre,
                    "pos": pos,
                    "consensus_support": support_n,
                    "change": f"{consensus_base}>{spic_base}",
                    "spret_base": spret_seq[pos - 1],
                    "spic_motif": clean_name(spic_hit["motif_name"]),
                    "spic_family": spic_family,
                    "spic_rel_score": round(spic_hit["rel_score"], 3),
                    "spic_window": f"{spic_hit['start']}-{spic_hit['end']}",
                    "other_motif": clean_name(consensus_hit["motif_name"]),
                    "other_family": other_family,
                    "other_rel_score": round(consensus_hit["rel_score"], 3),
                    "other_window": f"{consensus_hit['start']}-{consensus_hit['end']}",
                    "family_switch": family_switch,
                    "delta_rel_score": round(spic_hit["rel_score"] - consensus_hit["rel_score"], 3),
                }
            )

    return pd.DataFrame(rows), metrics


def build_full_site_table(selected_switches, background_switches):
    switch_lookup = {}
    for set_name, df in [("selected", selected_switches), ("background_matched", background_switches)]:
        for _, row in df.iterrows():
            switch_lookup[(set_name, row["CRE"], int(row["pos"]), row["change"])] = {
                "motif_switch": True,
                "family_switch": bool(row["family_switch"]),
                "spic_motif": row["spic_motif"],
                "other_motif": row["other_motif"],
                "spic_family": row["spic_family"],
                "other_family": row["other_family"],
                "spic_rel_score": row["spic_rel_score"],
                "other_rel_score": row["other_rel_score"],
            }

    rows = []
    comp = []
    for set_name, folder in [("selected", SELECTED_DIR), ("background_matched", BACKGROUND_DIR)]:
        for fa_path in sorted(folder.glob("*.fa")):
            cre = fa_path.stem
            seqs = parse_fasta(fa_path)
            cols = spic_alignment_columns(seqs)
            spic_seq = seq_on_spic_coords(seqs["mus_spicilegus"], cols)
            consensus_seq, support = consensus_on_spic_coords(seqs, cols)
            cre_gc = gc_fraction(consensus_seq)
            comp.append({"set": set_name, "CRE": cre, "GC": cre_gc, "length": len(consensus_seq)})

            for pos, (spic_base, consensus_base, support_n) in enumerate(zip(spic_seq, consensus_seq, support), start=1):
                if spic_base not in "ACGT" or consensus_base not in "ACGT":
                    continue
                if spic_base == consensus_base:
                    continue
                change = f"{consensus_base}>{spic_base}"
                start = max(1, pos - 5)
                end = min(len(consensus_seq), pos + 5)
                local_other = consensus_seq[start - 1 : end]
                info = switch_lookup.get((set_name, cre, pos, change))
                rows.append(
                    {
                        "set": set_name,
                        "CRE": cre,
                        "pos": pos,
                        "change": change,
                        "support": support_n,
                        "motif_switch": info is not None,
                        "family_switch": False if info is None else bool(info["family_switch"]),
                        "spic_motif": "none" if info is None else info["spic_motif"],
                        "other_motif": "none" if info is None else info["other_motif"],
                        "spic_family": "none" if info is None else info["spic_family"],
                        "other_family": "none" if info is None else info["other_family"],
                        "spic_rel_score": np.nan if info is None else info["spic_rel_score"],
                        "other_rel_score": np.nan if info is None else info["other_rel_score"],
                        "local_gc_other": gc_fraction(local_other),
                        "cre_gc": cre_gc,
                    }
                )
    comp_df = pd.DataFrame(comp)
    comp_df.to_csv(OUTDIR / "CRE_composition.tsv", sep="\t", index=False)
    full = pd.DataFrame(rows)
    full["base_class"] = full["change"].map(classify_change)
    full["is_selected"] = full["set"].eq("selected").astype(int)
    full.to_csv(OUTDIR / "spic_diff_site_context.tsv", sep="\t", index=False)
    return full


def motif_lookup(selected_switches):
    df = selected_switches.copy()
    df["Position"] = [genomic_position(cre, pos) for cre, pos in zip(df["CRE"], df["pos"])]
    out = pd.DataFrame(
        {
            "cCRE": df["CRE"],
            "Position": df["Position"],
            "Ancestral motif": df["other_motif"],
            "Spicilegus motif": df["spic_motif"],
            "Ancestral family": df["other_family"],
            "Spicilegus family": df["spic_family"],
            "Base change": df["change"],
            "Ancestral motif score": df["other_rel_score"],
            "Spicilegus motif score": df["spic_rel_score"],
            "Family switch": df["family_switch"],
        }
    )
    out = out.sort_values(["cCRE", "Position", "Ancestral motif", "Spicilegus motif"])
    out.to_csv(OUTDIR / "accelerated_CRE_TF_motif_switch_lookup.tsv", sep="\t", index=False)
    out[out["Family switch"] == True].to_csv(OUTDIR / "accelerated_CRE_TF_family_switch_lookup.tsv", sep="\t", index=False)
    return out


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


def test_motif_to_motif(full):
    switches = full[full["motif_switch"]].copy()
    switches["transition"] = switches["other_motif"].astype(str) + " -> " + switches["spic_motif"].astype(str)
    switches = assign_strata(switches)
    categories, observed, expected, sd, p_greater, p_less, p_two = permutation_all_categories(switches, "transition")

    rows = []
    selected_total = int((switches["set"] == "selected").sum())
    background_total = int((switches["set"] == "background_matched").sum())
    for i, transition in enumerate(categories):
        src, dst = transition.split(" -> ", 1)
        selected_hit = int(observed[i])
        background_hit = int(((switches["set"] == "background_matched") & (switches["transition"] == transition)).sum())
        _, p_unadj = fisher_exact(
            [[selected_hit, selected_total - selected_hit], [background_hit, background_total - background_hit]],
            alternative="greater",
        )
        or_unadj = ((selected_hit + 0.5) * (background_total - background_hit + 0.5)) / (
            (selected_total - selected_hit + 0.5) * (background_hit + 0.5)
        )
        rows.append(
            {
                "ancestral_motif": src,
                "spicilegus_motif": dst,
                "transition": transition,
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
    tests.to_csv(OUTDIR / "gc_matched_motif_to_motif_transition_tests.tsv", sep="\t", index=False)
    return switches, tests


def make_pass_sheets(switches, tests):
    fdr_pass = tests[(tests["matched_q_greater"] < 0.05) | (tests["matched_q_two_sided"] < 0.05)].copy()
    nominal_pass = tests[(tests["matched_p_greater"] < 0.05) | (tests["matched_p_two_sided"] < 0.05)].copy()

    def site_sheet(pass_tests):
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
        sites = selected.merge(pass_tests, on="transition", how="inner")
        out = pd.DataFrame(
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
        )
        return out.sort_values(["matched_q_greater", "matched_p_greater", "cCRE", "Position"])

    fdr_sites = site_sheet(fdr_pass)
    nominal_sites = site_sheet(nominal_pass)
    fdr_sites.to_csv(OUTDIR / "accelerated_CRE_TF_motif_switches_GC_FDR_pass.tsv", sep="\t", index=False)
    nominal_sites.to_csv(OUTDIR / "accelerated_CRE_TF_motif_switches_GC_nominal_pass.tsv", sep="\t", index=False)

    summary = pd.DataFrame(
        [
            {"metric": "selected_MSA_count", "value": len(list(SELECTED_DIR.glob("*.fa")))},
            {"metric": "background_MSA_count", "value": len(list(BACKGROUND_DIR.glob("*.fa")))},
            {"metric": "tested_motif_to_motif_transitions", "value": int(tests.shape[0])},
            {"metric": "FDR_pass_transitions", "value": int(fdr_pass.shape[0])},
            {"metric": "FDR_pass_accelerated_CRE_site_rows", "value": int(fdr_sites.shape[0])},
            {"metric": "nominal_matched_p_pass_transitions", "value": int(nominal_pass.shape[0])},
            {"metric": "nominal_matched_p_pass_accelerated_CRE_site_rows", "value": int(nominal_sites.shape[0])},
            {"metric": "definition", "value": "FDR pass = matched_q_greater < 0.05 or matched_q_two_sided < 0.05 after GC/base-change matched permutation testing."},
        ]
    )
    summary.to_csv(OUTDIR / "gc_matched_motif_to_motif_transition_summary.tsv", sep="\t", index=False)

    try:
        with pd.ExcelWriter(OUTDIR / "accelerated_CRE_TF_motif_switches_passing_GC_correction.xlsx", engine="openpyxl") as writer:
            summary.to_excel(writer, sheet_name="README", index=False)
            fdr_pass.to_excel(writer, sheet_name="GC_FDR_pass_motif_pairs", index=False)
            fdr_sites.to_excel(writer, sheet_name="GC_FDR_pass_CREs", index=False)
            nominal_pass.to_excel(writer, sheet_name="nominal_p_pass_motif_pairs", index=False)
            nominal_sites.to_excel(writer, sheet_name="nominal_p_pass_CREs", index=False)
            tests.to_excel(writer, sheet_name="all_motif_pair_tests", index=False)
    except ModuleNotFoundError:
        pass
    return summary, fdr_pass, nominal_pass


def main():
    OUTDIR.mkdir(exist_ok=True)
    selected_switches, background_switches, set_summary = scan_motifs()
    full = build_full_site_table(selected_switches, background_switches)
    motif_lookup(selected_switches)
    switches, tests = test_motif_to_motif(full)
    summary, fdr_pass, nominal_pass = make_pass_sheets(switches, tests)

    print("Motif-switch summary:")
    print(set_summary.to_string(index=False))
    print("\nGC-matched motif-to-motif summary:")
    print(summary.to_string(index=False))
    print("\nTop motif-pair tests:")
    print(tests.head(20).to_string(index=False))
    if not fdr_pass.empty:
        print("\nFDR-pass motif pairs:")
        print(fdr_pass.to_string(index=False))
    elif not nominal_pass.empty:
        print("\nNo FDR-pass motif pairs. Nominal matched-p motif pairs:")
        print(nominal_pass.to_string(index=False))
    else:
        print("\nNo FDR-pass or nominal matched-p motif pairs.")


if __name__ == "__main__":
    main()
