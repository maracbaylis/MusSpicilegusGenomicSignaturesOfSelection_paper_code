#!/usr/bin/env python3
from collections import Counter, OrderedDict
from pathlib import Path
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import motifs
from Bio.Seq import Seq
from scipy.stats import fisher_exact, mannwhitneyu


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
INPUT_DIR = Path(os.environ.get("MSPIC_INPUT_DIR", "data/inputs")).resolve()
JASPAR = Path(os.environ.get("MSPIC_JASPAR_PATH", INPUT_DIR / "JASPAR.jaspar"))
SETS = {"acc": "selected_MSAs", "bg": "background_MSAs_matched"}
MIN_REL_SCORE = 0.85
HIGH_DELTA = 0.15


def parse_fasta(path):
    seqs = OrderedDict()
    name = None
    parts = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(parts)
            name = line[1:].strip()
            parts = []
        else:
            parts.append(line.strip())
    if name is not None:
        seqs[name] = "".join(parts)
    return seqs


def load_pssms(path):
    bg = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
    out = []
    with open(path) as handle:
        for motif in motifs.parse(handle, "jaspar"):
            pwm = motif.counts.normalize(pseudocounts=0.5)
            pssm = pwm.log_odds(bg)
            out.append(
                {
                    "motif_id": motif.matrix_id,
                    "motif_name": motif.name,
                    "length": len(motif),
                    "pssm": pssm,
                    "min": float(pssm.min),
                    "max": float(pssm.max),
                }
            )
    return out


def clean_name(name):
    parts = str(name).split(".")
    return ".".join(parts[2:]) if len(parts) >= 3 else str(name)


def tf_family(name):
    clean = clean_name(name)
    primary = clean.split("::")[0]
    upper = primary.upper()
    if upper.startswith(("FOX",)):
        return "forkhead"
    if upper.startswith(("ATF", "CEBP", "BACH", "FOS", "JUN", "MAF", "DDIT3", "XBP", "CREB")):
        return "bZIP"
    if upper.startswith(("ATOH", "ASCL", "BHLH", "HAND", "HES", "HEY", "MYF", "MYOD", "MYOG", "NEUROD", "NEUROG", "NHLH", "OLIG", "TCF12", "TCF21", "TWIST", "USF", "MAX", "MITF", "MSC", "TAL", "SCL")):
        return "bHLH"
    if upper.startswith(("BARHL", "DLX", "HOX", "ARX", "ALX", "MSX", "POU", "PITX", "LHX", "SIX", "OTX", "EMX", "GBX", "EN", "ISL", "RHOX", "CDX")):
        return "homeobox"
    if upper.startswith(("RXR", "RARA", "RARB", "RARG", "THRA", "THRB", "NR", "PPAR", "ESR", "HNF4", "ROR", "VDR", "PGR", "AR", "GR", "PPARD")):
        return "nuclear receptor"
    if upper.startswith(("ZNF", "ZBTB", "ZFP", "ZIC", "WT1", "GLIS", "KLF", "SP", "RREB", "PLAG", "MAZ", "PATZ", "EGR", "INSM")):
        return "C2H2 zinc finger"
    return "unknown"


def spic_cols(seqs):
    return [i for i, base in enumerate(seqs["mus_spicilegus"]) if base != "-"]


def seq_on_cols(seq, cols):
    return "".join(seq[i] if seq[i] in "ACGT" else "N" for i in cols)


def consensus_on_cols(seqs, cols):
    others = [name for name in seqs if name != "mus_spicilegus"]
    consensus = []
    support = []
    for col in cols:
        bases = [seqs[name][col] for name in others if seqs[name][col] in "ACGT"]
        if not bases:
            consensus.append("N")
            support.append(0)
        else:
            base, n = Counter(bases).most_common(1)[0]
            consensus.append(base)
            support.append(n)
    return "".join(consensus), support


def base_class(ref, alt):
    if ref in "AT" and alt in "GC":
        return "AT_to_GC"
    if ref in "GC" and alt in "AT":
        return "GC_to_AT"
    return "other"


def rel_score_for_kmer(motif, kmer):
    if set(kmer) - set("ACGT"):
        return None
    plus = float(motif["pssm"].calculate(kmer))
    minus = float(motif["pssm"].calculate(str(Seq(kmer).reverse_complement())))
    score = max(plus, minus)
    denom = motif["max"] - motif["min"]
    if denom <= 0:
        return None
    return (score - motif["min"]) / denom


def best_rel_overlap(seq, pos, motif):
    length = motif["length"]
    seq_len = len(seq)
    best = None
    best_start = None
    start_min = max(1, pos - length + 1)
    start_max = min(pos, seq_len - length + 1)
    for start in range(start_min, start_max + 1):
        kmer = seq[start - 1 : start - 1 + length]
        rel = rel_score_for_kmer(motif, kmer)
        if rel is None:
            continue
        if best is None or rel > best:
            best = rel
            best_start = start
    return best, best_start


def scan_mutation_site(consensus_seq, spic_seq, pos, pssms):
    best_abs = None
    best_gain = None
    best_loss = None
    best_spic_any = None
    best_cons_any = None

    for motif in pssms:
        cons_rel, cons_start = best_rel_overlap(consensus_seq, pos, motif)
        spic_rel, spic_start = best_rel_overlap(spic_seq, pos, motif)
        if cons_rel is None or spic_rel is None:
            continue
        delta = spic_rel - cons_rel
        row = {
            "motif_id": motif["motif_id"],
            "motif_name": clean_name(motif["motif_name"]),
            "motif_family": tf_family(motif["motif_name"]),
            "consensus_rel_score": cons_rel,
            "spic_rel_score": spic_rel,
            "delta_rel_score": delta,
            "abs_delta_rel_score": abs(delta),
            "consensus_window": f"{cons_start}-{cons_start + motif['length'] - 1}",
            "spic_window": f"{spic_start}-{spic_start + motif['length'] - 1}",
        }
        if best_abs is None or row["abs_delta_rel_score"] > best_abs["abs_delta_rel_score"]:
            best_abs = row
        if delta > 0 and (best_gain is None or delta > best_gain["delta_rel_score"]):
            best_gain = row
        if delta < 0 and (best_loss is None or delta < best_loss["delta_rel_score"]):
            best_loss = row
        if best_spic_any is None or spic_rel > best_spic_any["spic_rel_score"]:
            best_spic_any = row
        if best_cons_any is None or cons_rel > best_cons_any["consensus_rel_score"]:
            best_cons_any = row

    return best_abs, best_gain, best_loss, best_spic_any, best_cons_any


def collect_sites(pssms):
    rows = []
    for set_name, folder in SETS.items():
        for fasta in sorted((HERE / folder).glob("*.fa")):
            seqs = parse_fasta(fasta)
            cols = spic_cols(seqs)
            spic_seq = seq_on_cols(seqs["mus_spicilegus"], cols)
            consensus_seq, support = consensus_on_cols(seqs, cols)
            for pos, (ref, alt, sup) in enumerate(zip(consensus_seq, spic_seq, support), start=1):
                if ref not in "ACGT" or alt not in "ACGT" or ref == alt:
                    continue
                best_abs, best_gain, best_loss, best_spic_any, best_cons_any = scan_mutation_site(
                    consensus_seq, spic_seq, pos, pssms
                )
                if best_abs is None:
                    continue
                gain = best_gain if best_gain is not None else {}
                loss = best_loss if best_loss is not None else {}
                rows.append(
                    {
                        "set": set_name,
                        "CRE": fasta.stem,
                        "pos": pos,
                        "change": f"{ref}>{alt}",
                        "base_class": base_class(ref, alt),
                        "consensus_support": sup,
                        "best_abs_motif": best_abs["motif_name"],
                        "best_abs_family": best_abs["motif_family"],
                        "best_abs_consensus_rel": best_abs["consensus_rel_score"],
                        "best_abs_spic_rel": best_abs["spic_rel_score"],
                        "best_abs_delta_rel": best_abs["delta_rel_score"],
                        "best_abs_abs_delta_rel": best_abs["abs_delta_rel_score"],
                        "best_gain_motif": gain.get("motif_name"),
                        "best_gain_family": gain.get("motif_family"),
                        "best_gain_delta_rel": gain.get("delta_rel_score", 0.0),
                        "best_gain_spic_rel": gain.get("spic_rel_score", np.nan),
                        "best_loss_motif": loss.get("motif_name"),
                        "best_loss_family": loss.get("motif_family"),
                        "best_loss_delta_rel": loss.get("delta_rel_score", 0.0),
                        "best_loss_consensus_rel": loss.get("consensus_rel_score", np.nan),
                        "best_spic_any_motif": best_spic_any["motif_name"],
                        "best_spic_any_family": best_spic_any["motif_family"],
                        "best_spic_any_rel": best_spic_any["spic_rel_score"],
                        "best_consensus_any_motif": best_cons_any["motif_name"],
                        "best_consensus_any_family": best_cons_any["motif_family"],
                        "best_consensus_any_rel": best_cons_any["consensus_rel_score"],
                    }
                )
    return pd.DataFrame(rows)


def fisher_row(df, group_col, group_a, group_b, outcome_col, label):
    a = df[df[group_col].eq(group_a)]
    b = df[df[group_col].eq(group_b)]
    table = [
        [int(a[outcome_col].sum()), int((~a[outcome_col]).sum())],
        [int(b[outcome_col].sum()), int((~b[outcome_col]).sum())],
    ]
    odds, p = fisher_exact(table, alternative="greater")
    return {
        "test": label,
        "group_a": group_a,
        "group_b": group_b,
        "group_a_hits": table[0][0],
        "group_a_total": len(a),
        "group_a_fraction": table[0][0] / len(a),
        "group_b_hits": table[1][0],
        "group_b_total": len(b),
        "group_b_fraction": table[1][0] / len(b),
        "odds_ratio": odds,
        "p_greater": p,
    }


def summarize(df):
    df = df.copy()
    df["high_abs_delta"] = df["best_abs_abs_delta_rel"] >= HIGH_DELTA
    df["high_gain_delta"] = df["best_gain_delta_rel"] >= HIGH_DELTA
    df["high_loss_delta"] = df["best_loss_delta_rel"] <= -HIGH_DELTA
    df["crosses_gain_threshold"] = (df["best_abs_spic_rel"] >= MIN_REL_SCORE) & (df["best_abs_consensus_rel"] < MIN_REL_SCORE)
    df["crosses_loss_threshold"] = (df["best_abs_consensus_rel"] >= MIN_REL_SCORE) & (df["best_abs_spic_rel"] < MIN_REL_SCORE)
    df["crosses_any_threshold"] = df["crosses_gain_threshold"] | df["crosses_loss_threshold"]
    df.to_csv(HERE / "mutation_motif_score_delta_sites.tsv", sep="\t", index=False)

    rows = []
    for set_name, group in df.groupby("set"):
        for cls, sub in group.groupby("base_class"):
            rows.append(
                {
                    "set": set_name,
                    "base_class": cls,
                    "n_sites": len(sub),
                    "median_abs_delta": sub["best_abs_abs_delta_rel"].median(),
                    "median_gain_delta": sub["best_gain_delta_rel"].median(),
                    "median_loss_delta": sub["best_loss_delta_rel"].median(),
                    "high_abs_delta_fraction": sub["high_abs_delta"].mean(),
                    "crosses_any_threshold_fraction": sub["crosses_any_threshold"].mean(),
                }
            )
    pd.DataFrame(rows).to_csv(HERE / "mutation_motif_score_delta_summary_by_class.tsv", sep="\t", index=False)

    tests = []
    for set_name in ["acc", "bg"]:
        sub = df[df["set"].eq(set_name)]
        atgc = sub[sub["base_class"].eq("AT_to_GC")]
        other = sub[~sub["base_class"].eq("AT_to_GC")]
        for metric in ["best_abs_abs_delta_rel", "best_gain_delta_rel", "best_loss_delta_rel"]:
            stat, p = mannwhitneyu(atgc[metric], other[metric], alternative="two-sided")
            tests.append(
                {
                    "test": f"{set_name}_AT_to_GC_vs_other_{metric}",
                    "n_atgc": len(atgc),
                    "n_other": len(other),
                    "median_atgc": atgc[metric].median(),
                    "median_other": other[metric].median(),
                    "mannwhitney_u": stat,
                    "p_value": p,
                }
            )
        for outcome in ["high_abs_delta", "high_gain_delta", "high_loss_delta", "crosses_any_threshold"]:
            tests.append(fisher_row(sub.assign(is_atgc=sub["base_class"].eq("AT_to_GC")), "is_atgc", True, False, outcome, f"{set_name}_AT_to_GC_enriched_{outcome}"))

    for outcome in ["high_abs_delta", "high_gain_delta", "high_loss_delta", "crosses_any_threshold"]:
        tests.append(fisher_row(df, "set", "acc", "bg", outcome, f"acc_vs_bg_enriched_{outcome}"))

    pd.DataFrame(tests).to_csv(HERE / "mutation_motif_score_delta_tests.tsv", sep="\t", index=False)
    return df, pd.DataFrame(tests)


def plot(df):
    plt.rcParams.update({"font.family": "Arial", "font.size": 10, "pdf.fonttype": 42, "ps.fonttype": 42})
    colors = {"AT_to_GC": "#238b45", "GC_to_AT": "#e6550d", "other": "#969696"}
    high_cutoff = df.loc[df["set"].eq("bg"), "best_abs_abs_delta_rel"].quantile(0.90)
    df = df.copy()
    df["high_abs_bg90"] = df["best_abs_abs_delta_rel"] >= high_cutoff

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.8), constrained_layout=True)
    ax = axes[0]
    labels = []
    data = []
    box_colors = []
    for set_name, set_label in [("acc", "Accelerated"), ("bg", "Background")]:
        for cls in ["AT_to_GC", "GC_to_AT", "other"]:
            sub = df[(df["set"].eq(set_name)) & (df["base_class"].eq(cls))]
            labels.append(f"{set_label}\n{cls.replace('_', '->')}")
            data.append(sub["best_abs_abs_delta_rel"].to_numpy())
            box_colors.append(colors[cls])
    bp = ax.boxplot(data, labels=labels, patch_artist=True, showfliers=False)
    for patch, color in zip(bp["boxes"], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
    ax.set_ylabel("Largest motif-score change per mutation\n(abs delta relative score)")
    ax.tick_params(axis="x", labelrotation=45)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[1]
    frac_rows = (
        df.groupby(["set", "base_class"])["high_abs_bg90"]
        .agg(["sum", "count", "mean"])
        .reset_index()
    )
    x = np.arange(len(frac_rows))
    ax.bar(x, frac_rows["mean"], color=[colors[c] for c in frac_rows["base_class"]])
    ax.set_xticks(
        x,
        [
            f"{'Accelerated' if s == 'acc' else 'Background'}\n{c.replace('_', '->')}\n{int(n)}/{int(t)}"
            for s, c, n, t in zip(frac_rows["set"], frac_rows["base_class"], frac_rows["sum"], frac_rows["count"])
        ],
        rotation=45,
        ha="right",
    )
    ax.set_ylabel(f"Fraction high-impact\n(abs delta >= background 90th pct.; {high_cutoff:.3f})")
    ax.set_ylim(0, max(0.05, frac_rows["mean"].max() * 1.25))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for ext in ["png", "pdf"]:
        plt.savefig(HERE / f"fig_mutation_motif_score_delta.{ext}", dpi=300, bbox_inches="tight")
    print("Wrote: fig_mutation_motif_score_delta.png/pdf")


def main():
    pssms = load_pssms(JASPAR)
    df = collect_sites(pssms)
    df, tests = summarize(df)
    plot(df)
    print("Wrote: mutation_motif_score_delta_sites.tsv")
    print("Wrote: mutation_motif_score_delta_summary_by_class.tsv")
    print("Wrote: mutation_motif_score_delta_tests.tsv")
    print(tests.to_string(index=False))


if __name__ == "__main__":
    main()
