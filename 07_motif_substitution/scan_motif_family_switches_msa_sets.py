#!/usr/bin/env python3
from collections import Counter, OrderedDict
from pathlib import Path
import os
import re

import pandas as pd
from Bio import motifs
from Bio.Seq import Seq


MIN_REL_SCORE = 0.85
SETS = {
    "selected": "selected_MSAs",
    "background_matched": "background_MSAs_matched",
}


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
    return seqs


def spic_alignment_columns(seqs):
    return [i for i, base in enumerate(seqs["mus_spicilegus"]) if base != "-"]


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
        cnt = Counter(bases)
        base, n = cnt.most_common(1)[0]
        consensus.append(base)
        support.append(n)
    return "".join(consensus), support


def seq_on_spic_coords(seq, cols):
    out = []
    for col in cols:
        base = seq[col]
        out.append(base if base in "ACGT" else "N")
    return "".join(out)


def find_taxon(seqs, target):
    for name in seqs:
        if name == target:
            return name
    for name in seqs:
        if target in name:
            return name
    return None


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


def score_kmer(pssm, kmer):
    if set(kmer) - set("ACGT"):
        return None
    plus = float(pssm.calculate(kmer))
    minus = float(pssm.calculate(str(Seq(kmer).reverse_complement())))
    if plus >= minus:
        return plus, "+"
    return minus, "-"


def best_motif_overlap(seq, pos, pssms):
    best = None
    seq_len = len(seq)
    for motif in pssms:
        length = motif["length"]
        best_score = None
        best_start = None
        best_strand = None
        best_kmer = None

        start_min = max(1, pos - length + 1)
        start_max = min(pos, seq_len - length + 1)
        for start in range(start_min, start_max + 1):
            kmer = seq[start - 1 : start - 1 + length]
            scored = score_kmer(motif["pssm"], kmer)
            if scored is None:
                continue
            score, strand = scored
            if best_score is None or score > best_score:
                best_score = score
                best_start = start
                best_strand = strand
                best_kmer = kmer

        if best_score is None:
            continue

        rel = (best_score - motif["min"]) / (motif["max"] - motif["min"]) if motif["max"] > motif["min"] else 0.0
        hit = {
            "motif_id": motif["motif_id"],
            "motif_name": motif["motif_name"],
            "rel_score": rel,
            "start": best_start,
            "end": best_start + length - 1,
            "strand": best_strand,
            "kmer": best_kmer,
        }
        if best is None or hit["rel_score"] > best["rel_score"]:
            best = hit
    return best


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


def analyze_set(set_name, msa_dir, pssms):
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
            if max(spic_hit["rel_score"], consensus_hit["rel_score"]) < MIN_REL_SCORE:
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


def main():
    here = Path(__file__).resolve().parent
    input_dir = Path(os.environ.get("MSPIC_INPUT_DIR", "data/inputs")).resolve()
    jaspar_path = Path(
        os.environ.get("MSPIC_JASPAR_PATH", input_dir / "JASPAR.jaspar")
    )
    pssms = load_pssms(jaspar_path)

    all_rows = []
    all_metrics = []
    for set_name, folder in SETS.items():
        rows, metrics = analyze_set(set_name, here / folder, pssms)
        rows = rows.sort_values(["family_switch", "delta_rel_score", "spic_rel_score"], ascending=[False, False, False])
        rows.to_csv(here / f"{set_name}_motif_switch_sites.tsv", sep="\t", index=False)
        all_rows.append(rows)
        all_metrics.append(metrics)

    summary = pd.DataFrame(all_metrics)
    summary["motif_switch_fraction_of_diff_sites"] = summary["motif_switch_sites"] / summary["spic_diff_sites"]
    summary["family_switch_fraction_of_diff_sites"] = summary["family_switch_sites"] / summary["spic_diff_sites"]
    summary["family_switch_fraction_of_motif_switches"] = summary["family_switch_sites"] / summary["motif_switch_sites"]
    summary.to_csv(here / "motif_switch_set_summary.tsv", sep="\t", index=False)

    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
