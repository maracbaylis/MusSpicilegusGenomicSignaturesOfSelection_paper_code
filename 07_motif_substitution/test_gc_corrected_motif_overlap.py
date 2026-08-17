#!/usr/bin/env python3
from collections import Counter, OrderedDict
from pathlib import Path
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import fisher_exact


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MSA_SEARCH_ROOTS = [
    ROOT / "phyloP_MSA_QC",
    HERE,
    Path(os.environ.get("MSPIC_INPUT_DIR", "data/inputs")).resolve(),
]
IN_SUBS = HERE / "substitution_features.tsv"

OUT_AUGMENTED = HERE / "motif_overlap_substitution_features_gc_context.tsv"
OUT_SUMMARY = HERE / "gc_corrected_motif_overlap_results.tsv"
OUT_PNG = HERE / "fig_gc_corrected_motif_overlap.png"
OUT_PDF = HERE / "fig_gc_corrected_motif_overlap.pdf"

N_PERM = 20000
SEED = 53


def parse_fasta(path):
    seqs = OrderedDict()
    name = None
    parts = []
    for line in path.read_text().splitlines():
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
        seqs = OrderedDict((name, seq + "-" * (max_len - len(seq))) for name, seq in seqs.items())
    return seqs


def spic_alignment_columns(seqs):
    return [i for i, base in enumerate(seqs["mus_spicilegus"]) if base != "-"]


def seq_on_spic_coords(seq, cols):
    return "".join(seq[col] if seq[col] in "ACGT" else "N" for col in cols)


def consensus_on_spic_coords(seqs, cols):
    others = [name for name in seqs if name != "mus_spicilegus"]
    out = []
    for col in cols:
        bases = [seqs[name][col] for name in others if seqs[name][col] in "ACGT"]
        if not bases:
            out.append("N")
        else:
            out.append(Counter(bases).most_common(1)[0][0])
    return "".join(out)


def gc_fraction(seq):
    bases = [base for base in seq if base in "ACGT"]
    if not bases:
        return np.nan
    return sum(base in "GC" for base in bases) / len(bases)


def add_gc_context():
    if OUT_AUGMENTED.exists():
        return pd.read_csv(OUT_AUGMENTED, sep="\t")

    subs = pd.read_csv(IN_SUBS, sep="\t")
    needed = set(subs["CRE"].astype(str).unique())
    fasta_map = {}
    for root in MSA_SEARCH_ROOTS:
        if not root.exists():
            continue
        for fa in root.rglob("*.fa"):
            if fa.stem in needed and fa.stem not in fasta_map:
                fasta_map[fa.stem] = fa

    rows = []
    missing = 0
    for cre, group in subs.groupby("CRE", sort=False):
        fa = fasta_map.get(cre)
        if fa is None:
            missing += len(group)
            continue
        seqs = parse_fasta(fa)
        cols = spic_alignment_columns(seqs)
        consensus = consensus_on_spic_coords(seqs, cols)
        cre_gc = gc_fraction(consensus)
        for row in group.itertuples(index=False):
            pos = int(row.pos)
            if pos < 1 or pos > len(consensus):
                local_gc = np.nan
            else:
                start = max(1, pos - 5)
                end = min(len(consensus), pos + 5)
                local_gc = gc_fraction(consensus[start - 1 : end])
            rows.append(
                {
                    "CRE": cre,
                    "pos": pos,
                    "is_spic": int(row.is_spic),
                    "in_motif": int(row.in_motif),
                    "local_gc": local_gc,
                    "cre_gc": cre_gc,
                }
            )
    out = pd.DataFrame(rows).dropna(subset=["local_gc", "cre_gc"])
    out.to_csv(OUT_AUGMENTED, sep="\t", index=False)
    if missing:
        print(f"Skipped {missing} rows with missing MSA files")
    return out


def qcut_strata(df):
    out = df.copy()
    for col in ["local_gc", "cre_gc"]:
        out[f"{col}_bin"] = pd.qcut(out[col], q=4, duplicates="drop").astype(str)
    out["stratum"] = out["local_gc_bin"].astype(str) + "|cre=" + out["cre_gc_bin"].astype(str)
    return out


def stratified_permutation(df):
    df = qcut_strata(df).reset_index(drop=True)
    in_motif = df["in_motif"].astype(int).to_numpy()
    is_spic = df["is_spic"].astype(int).to_numpy()
    observed = int(((in_motif == 1) & (is_spic == 1)).sum())
    total_spic = int(is_spic.sum())

    rng = np.random.default_rng(SEED)
    counts = np.zeros(N_PERM, dtype=int)
    for _, sub in df.groupby("stratum", sort=False):
        idx = sub.index.to_numpy()
        n_spic = int(sub["is_spic"].sum())
        if n_spic <= 0:
            continue
        if n_spic >= len(idx):
            counts += int(in_motif[idx].sum())
            continue
        for i in range(N_PERM):
            chosen = rng.choice(idx, size=n_spic, replace=False)
            counts[i] += int(in_motif[chosen].sum())

    expected = float(counts.mean())
    p_greater = float((1 + np.sum(counts >= observed)) / (N_PERM + 1))
    p_less = float((1 + np.sum(counts <= observed)) / (N_PERM + 1))
    p_two = float(min(1.0, 2 * min(p_greater, p_less)))
    return observed, total_spic, expected, float(counts.std(ddof=1)), p_greater, p_two


def logistic_regression(df):
    x = pd.DataFrame(
        {
            "const": 1.0,
            "in_motif": df["in_motif"].astype(float),
            "local_gc": df["local_gc"].astype(float),
            "cre_gc": df["cre_gc"].astype(float),
        }
    )
    y = df["is_spic"].astype(float)
    fit = sm.Logit(y, x).fit(disp=False, maxiter=200)
    conf = fit.conf_int()
    return {
        "logit_motif_or": float(np.exp(fit.params["in_motif"])),
        "logit_motif_ci_low": float(np.exp(conf.loc["in_motif", 0])),
        "logit_motif_ci_high": float(np.exp(conf.loc["in_motif", 1])),
        "logit_motif_p": float(fit.pvalues["in_motif"]),
    }


def make_plot(summary):
    raw_non = float(summary["raw_spic_fraction_non_motif"].iloc[0])
    raw_motif = float(summary["raw_spic_fraction_motif"].iloc[0])
    expected = float(summary["matched_expected_spic_in_motif"].iloc[0])
    observed = float(summary["observed_spic_in_motif"].iloc[0])
    oe = float(summary["matched_observed_over_expected"].iloc[0])

    plt.rcParams.update({"font.family": "Arial", "font.size": 9.5, "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 3.1), gridspec_kw={"width_ratios": [1.0, 1.0]})

    colors = ["#B8C99D", "#6F8F45"]
    axes[0].bar([0, 1], [raw_non, raw_motif], color=colors, edgecolor="none", width=0.62)
    axes[0].set_xticks([0, 1], ["Non motif", "Motif"])
    axes[0].set_ylabel(r"Fraction $\it{M.\ spicilegus}$ specific")
    axes[0].set_ylim(0, max(raw_non, raw_motif) * 1.45)
    axes[0].text(0, raw_non + max(raw_non, raw_motif) * 0.06, f"{raw_non*100:.1f}%", ha="center")
    axes[0].text(1, raw_motif + max(raw_non, raw_motif) * 0.06, f"{raw_motif*100:.1f}%", ha="center")
    axes[0].text(0.5, max(raw_non, raw_motif) * 1.28, "**", ha="center", fontsize=14, fontweight="bold")

    axes[1].bar([0, 1], [expected, observed], color=["#DDE5CF", "#6F8F45"], edgecolor="none", width=0.62)
    axes[1].set_xticks([0, 1], ["GC-matched\nexpected", "Observed"])
    axes[1].set_ylabel(r"$\it{M.\ spicilegus}$ substitutions\nin motif bases")
    axes[1].set_ylim(0, max(expected, observed) * 1.25)
    axes[1].text(0, expected + max(expected, observed) * 0.04, f"{expected:.1f}", ha="center")
    axes[1].text(1, observed + max(expected, observed) * 0.04, f"{observed:.0f}", ha="center")
    axes[1].text(0.5, max(expected, observed) * 1.14, "**", ha="center", fontsize=14, fontweight="bold")
    axes[1].text(0.5, -0.24, f"O/E = {oe:.2f}", transform=axes[1].transAxes, ha="center", va="top")

    for ax in axes:
        ax.grid(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.tight_layout(w_pad=2.0)
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(OUT_PDF, bbox_inches="tight")


def main():
    df = add_gc_context()
    n_motif = int(df["in_motif"].sum())
    n_non = int((df["in_motif"] == 0).sum())
    spic_motif = int(((df["in_motif"] == 1) & (df["is_spic"] == 1)).sum())
    spic_non = int(((df["in_motif"] == 0) & (df["is_spic"] == 1)).sum())

    _, raw_p = fisher_exact([[spic_motif, n_motif - spic_motif], [spic_non, n_non - spic_non]], alternative="greater")
    raw_or = ((spic_motif + 0.5) * (n_non - spic_non + 0.5)) / ((n_motif - spic_motif + 0.5) * (spic_non + 0.5))
    observed, total_spic, expected, sd, p_greater, p_two = stratified_permutation(df)
    logit = logistic_regression(df)

    summary = pd.DataFrame(
        [
            {
                "n_substitution_rows": int(df.shape[0]),
                "n_motif_rows": n_motif,
                "n_non_motif_rows": n_non,
                "spic_in_motif": spic_motif,
                "spic_in_non_motif": spic_non,
                "raw_spic_fraction_motif": spic_motif / n_motif,
                "raw_spic_fraction_non_motif": spic_non / n_non,
                "raw_or_motif_vs_non_motif": raw_or,
                "raw_p_greater": raw_p,
                "observed_spic_in_motif": observed,
                "total_spic": total_spic,
                "matched_expected_spic_in_motif": expected,
                "matched_observed_over_expected": observed / expected if expected else np.nan,
                "matched_sd": sd,
                "matched_p_greater": p_greater,
                "matched_p_two_sided": p_two,
                **logit,
                "note": "GC-corrected test matches/permutes M. spicilegus labels within local-GC and CRE-GC strata. Base-change class was not used because non-spic substitution changes are not stored in substitution_features.tsv.",
            }
        ]
    )
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    make_plot(summary)
    print(f"Wrote: {OUT_AUGMENTED}")
    print(f"Wrote: {OUT_SUMMARY}")
    print(f"Wrote: {OUT_PNG}")
    print(f"Wrote: {OUT_PDF}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
