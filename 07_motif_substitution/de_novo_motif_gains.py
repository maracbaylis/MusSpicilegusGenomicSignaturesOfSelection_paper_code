#!/usr/bin/env python3
from collections import Counter, OrderedDict
from pathlib import Path
import os

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from scan_motif_family_switches_msa_sets import load_pssms, best_motif_overlap, clean_name, tf_family
from test_family_transition_gc_matched_null import assign_strata, bh_fdr


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTDIR = HERE / "analysis"
SELECTED_DIR = Path(os.environ.get("MSPIC_INPUT_DIR", "data/inputs")).resolve() / 'selected_MSAs'
BACKGROUND_DIR = Path(os.environ.get("MSPIC_INPUT_DIR", "data/inputs")).resolve() / 'background_MSAs'
INPUT_DIR = Path(os.environ.get("MSPIC_INPUT_DIR", "data/inputs")).resolve()
JASPAR_PATH = Path(
    os.environ.get("MSPIC_JASPAR_PATH", INPUT_DIR / "JASPAR.jaspar")
)

MIN_REL_SCORE = 0.85
N_PERM = 20000
SEED = 47

OUT_SITES = OUTDIR / "de_novo_motif_gain_sites.tsv"
OUT_DEST_TESTS = OUTDIR / "gc_matched_de_novo_gain_destination_motif_tests.tsv"
OUT_SUMMARY = OUTDIR / "de_novo_motif_gain_summary.tsv"
OUT_FDR_SITES = OUTDIR / "de_novo_gain_destination_motif_GC_FDR_pass.tsv"
OUT_NOMINAL_SITES = OUTDIR / "de_novo_gain_destination_motif_GC_nominal_pass.tsv"
OUT_XLSX = OUTDIR / "de_novo_motif_gain_GC_corrected_results.xlsx"


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
    return "".join(seq[col] if seq[col] in "ACGT" else "N" for col in cols)


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
    parts = str(cre).split("__")
    if len(parts) != 2:
        return ""
    coord = parts[1].split("_")
    if len(coord) < 3:
        return ""
    return f"{coord[0]}:{int(coord[1]) + int(pos) - 1}"


def scan_sites():
    pssms = load_pssms(JASPAR_PATH)
    rows = []
    for set_name, folder in [("selected", SELECTED_DIR), ("background_matched", BACKGROUND_DIR)]:
        for fa_path in sorted(folder.glob("*.fa")):
            cre = fa_path.stem
            seqs = parse_fasta(fa_path)
            cols = spic_alignment_columns(seqs)
            spic_seq = seq_on_spic_coords(seqs["mus_spicilegus"], cols)
            consensus_seq, support = consensus_on_spic_coords(seqs, cols)
            cre_gc = gc_fraction(consensus_seq)
            for pos, (spic_base, consensus_base, support_n) in enumerate(zip(spic_seq, consensus_seq, support), start=1):
                if spic_base not in "ACGT" or consensus_base not in "ACGT" or spic_base == consensus_base:
                    continue
                spic_hit = best_motif_overlap(spic_seq, pos, pssms)
                other_hit = best_motif_overlap(consensus_seq, pos, pssms)
                spic_rel = np.nan if spic_hit is None else float(spic_hit["rel_score"])
                other_rel = np.nan if other_hit is None else float(other_hit["rel_score"])
                spic_has = bool(spic_hit is not None and spic_rel >= MIN_REL_SCORE)
                other_has = bool(other_hit is not None and other_rel >= MIN_REL_SCORE)
                if not (spic_has or other_has):
                    continue
                change = f"{consensus_base}>{spic_base}"
                local = consensus_seq[max(0, pos - 6) : min(len(consensus_seq), pos + 5)]
                state = "unchanged_above_threshold"
                if spic_has and not other_has:
                    state = "de_novo_gain"
                elif other_has and not spic_has:
                    state = "loss"
                elif spic_has and other_has and clean_name(spic_hit["motif_name"]) != clean_name(other_hit["motif_name"]):
                    state = "switch_between_existing_matches"
                rows.append(
                    {
                        "set": set_name,
                        "CRE": cre,
                        "Position": genomic_position(cre, pos),
                        "pos": pos,
                        "change": change,
                        "base_class": classify_change(change),
                        "consensus_support": support_n,
                        "local_gc_other": gc_fraction(local),
                        "cre_gc": cre_gc,
                        "state": state,
                        "spic_motif": "none" if spic_hit is None else clean_name(spic_hit["motif_name"]),
                        "spic_family": "none" if spic_hit is None else tf_family(spic_hit["motif_name"]),
                        "spic_rel_score": spic_rel,
                        "other_motif": "none" if other_hit is None else clean_name(other_hit["motif_name"]),
                        "other_family": "none" if other_hit is None else tf_family(other_hit["motif_name"]),
                        "other_rel_score": other_rel,
                        "spic_has_match": spic_has,
                        "other_has_match": other_has,
                    }
                )
    out = pd.DataFrame(rows)
    out.to_csv(OUT_SITES, sep="\t", index=False)
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
            counts += np.bincount(cat_idx[idx], minlength=len(categories))
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


def destination_tests(sites):
    gains = sites[sites["state"].eq("de_novo_gain")].copy()
    gains["is_selected"] = gains["set"].eq("selected").astype(int)
    gains = assign_strata(gains)
    categories, observed, expected, sd, p_greater, p_less, p_two = permutation_all_categories(gains, "spic_motif")
    selected_total = int(gains["set"].eq("selected").sum())
    bg_total = int(gains["set"].eq("background_matched").sum())
    rows = []
    for i, motif in enumerate(categories):
        selected_hit = int(observed[i])
        bg_hit = int(((gains["set"].eq("background_matched")) & (gains["spic_motif"].eq(motif))).sum())
        _, p_unadj = fisher_exact(
            [[selected_hit, selected_total - selected_hit], [bg_hit, bg_total - bg_hit]],
            alternative="greater",
        )
        or_unadj = ((selected_hit + 0.5) * (bg_total - bg_hit + 0.5)) / (
            (selected_total - selected_hit + 0.5) * (bg_hit + 0.5)
        )
        family = ";".join(sorted(map(str, gains.loc[gains["spic_motif"].eq(motif), "spic_family"].dropna().unique())))
        rows.append(
            {
                "destination_motif": motif,
                "destination_family": family,
                "selected_hit": selected_hit,
                "selected_total": selected_total,
                "selected_fraction": selected_hit / selected_total,
                "background_hit": bg_hit,
                "background_total": bg_total,
                "background_fraction": bg_hit / bg_total,
                "unadjusted_or": float(or_unadj),
                "unadjusted_p_greater": float(p_unadj),
                "matched_expected_selected_hit": float(expected[i]),
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
    tests.to_csv(OUT_DEST_TESTS, sep="\t", index=False)
    return gains, tests


def write_site_sheets(gains, tests):
    fdr = tests[(tests["matched_q_greater"] < 0.05) & (tests["matched_observed_over_expected"] > 1)].copy()
    nominal = tests[(tests["matched_p_greater"] < 0.05) & (tests["matched_observed_over_expected"] > 1)].copy()

    def sheet(pass_tests):
        if pass_tests.empty:
            return pd.DataFrame()
        selected = gains[gains["set"].eq("selected")].copy()
        out = selected.merge(pass_tests, left_on="spic_motif", right_on="destination_motif", how="inner")
        keep = [
            "CRE", "Position", "other_motif", "spic_motif", "other_family", "spic_family", "change",
            "other_rel_score", "spic_rel_score", "matched_observed_over_expected", "matched_p_greater",
            "matched_q_greater", "selected_hit", "background_hit",
        ]
        return out[keep].sort_values(["matched_q_greater", "matched_p_greater", "spic_motif", "CRE", "Position"])

    fdr_sites = sheet(fdr)
    nominal_sites = sheet(nominal)
    fdr_sites.to_csv(OUT_FDR_SITES, sep="\t", index=False)
    nominal_sites.to_csv(OUT_NOMINAL_SITES, sep="\t", index=False)
    return fdr, nominal, fdr_sites, nominal_sites


def main():
    OUTDIR.mkdir(exist_ok=True)
    sites = scan_sites()
    gains, tests = destination_tests(sites)
    fdr, nominal, fdr_sites, nominal_sites = write_site_sheets(gains, tests)
    summary = pd.DataFrame(
        [
            {"metric": "selected_de_novo_gain_sites", "value": int(gains["set"].eq("selected").sum())},
            {"metric": "background_de_novo_gain_sites", "value": int(gains["set"].eq("background_matched").sum())},
            {"metric": "tested_de_novo_destination_motifs", "value": int(tests.shape[0])},
            {"metric": "FDR_pass_destination_motifs", "value": int(fdr.shape[0])},
            {"metric": "FDR_pass_accelerated_site_rows", "value": int(fdr_sites.shape[0])},
            {"metric": "nominal_pass_destination_motifs", "value": int(nominal.shape[0])},
            {"metric": "nominal_pass_accelerated_site_rows", "value": int(nominal_sites.shape[0])},
            {"metric": "de_novo_gain_definition", "value": f"spic motif relative score >= {MIN_REL_SCORE} and consensus/reference best overlapping motif relative score < {MIN_REL_SCORE}"},
        ]
    )
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    try:
        with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as writer:
            summary.to_excel(writer, sheet_name="README", index=False)
            fdr.to_excel(writer, sheet_name="GC_FDR_pass_dest_motifs", index=False)
            fdr_sites.to_excel(writer, sheet_name="GC_FDR_pass_CREs", index=False)
            nominal.to_excel(writer, sheet_name="nominal_p_pass_dest_motifs", index=False)
            nominal_sites.to_excel(writer, sheet_name="nominal_p_pass_CREs", index=False)
            tests.to_excel(writer, sheet_name="all_destination_tests", index=False)
            sites.to_excel(writer, sheet_name="all_above_threshold_sites", index=False)
    except ModuleNotFoundError:
        pass
    print(f"Wrote: {OUT_SITES.name}")
    print(f"Wrote: {OUT_DEST_TESTS.name}")
    print(f"Wrote: {OUT_SUMMARY.name}")
    print(summary.to_string(index=False))
    print("\nTop de novo destination tests:")
    print(tests.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
