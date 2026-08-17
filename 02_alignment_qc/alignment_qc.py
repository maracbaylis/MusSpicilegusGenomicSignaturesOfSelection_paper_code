#!/usr/bin/env python3
"""
QC per-locus nucleotide multiple-sequence alignments for lineage-specific
evolution analyses.

Implements the alignment/QC logic described in Methods section 2:
  - equal alignment lengths
  - identifiable focal foreground
  - callable foreground coverage
  - foreground-private substitutions
  - private substitution rate
  - maximum private substitutions in any 10-bp window
  - longest serial cluster of private substitutions separated by <=2 bp
  - foreground-private gaps
  - manuscript artifact filters

A foreground-private substitution is counted when:
  * foreground base is A/C/G/T
  * at least 5 nonforeground taxa have callable A/C/G/T bases
  * the nonforeground majority base frequency is >= 0.75
  * foreground differs from that majority base

By default, loci are flagged FAIL if:
  private substitution rate > 0.08
  OR private substitution count >= 25
  OR max private substitutions in any 10-bp window >= 7
  OR foreground-private gaps >= 3
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

CALLABLE = set("ACGT")
MISSING = set("NRYKMSWBDHV?-")

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("fasta_dir", type=Path,
                   help="Directory containing one aligned FASTA per locus.")
    p.add_argument("--foreground", required=True,
                   help="Foreground sequence name, e.g. mus_spicilegus.")
    p.add_argument("--out", type=Path, required=True,
                   help="Output CSV.")
    p.add_argument("--min-nonforeground-callable", type=int, default=5)
    p.add_argument("--majority-frequency", type=float, default=0.75)
    p.add_argument("--max-private-rate", type=float, default=0.08)
    p.add_argument("--max-private-count", type=int, default=25)
    p.add_argument("--max-private-in-10bp", type=int, default=7)
    p.add_argument("--max-private-gaps", type=int, default=3)
    return p.parse_args()

def read_fasta(path: Path) -> dict[str, str]:
    records = {}
    header = None
    chunks = []
    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    if header in records:
                        raise ValueError(f"Duplicate FASTA header {header!r} in {path}")
                    records[header] = "".join(chunks).upper()
                header = line[1:].split()[0]
                chunks = []
            else:
                if header is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                chunks.append(line)
    if header is not None:
        if header in records:
            raise ValueError(f"Duplicate FASTA header {header!r} in {path}")
        records[header] = "".join(chunks).upper()
    return records

def longest_cluster(positions: list[int], max_sep: int = 2) -> int:
    if not positions:
        return 0
    best = cur = 1
    for a, b in zip(positions, positions[1:]):
        if b - a <= max_sep:
            cur += 1
        else:
            cur = 1
        best = max(best, cur)
    return best

def max_in_window(positions: list[int], width: int = 10) -> int:
    if not positions:
        return 0
    best = 0
    left = 0
    for right, pos in enumerate(positions):
        while pos - positions[left] >= width:
            left += 1
        best = max(best, right - left + 1)
    return best

def evaluate_alignment(path: Path, args) -> dict:
    seqs = read_fasta(path)
    lengths = [len(s) for s in seqs.values()]
    equal_length = len(set(lengths)) <= 1
    foreground_present = args.foreground in seqs

    row = {
        "locus": path.stem,
        "file": str(path),
        "n_sequences": len(seqs),
        "equal_alignment_length": equal_length,
        "foreground_present": foreground_present,
        "alignment_length": lengths[0] if lengths and equal_length else "",
        "foreground_callable_sites": 0,
        "foreground_callable_fraction": "",
        "private_substitution_count": 0,
        "private_substitution_rate": "",
        "max_private_substitutions_10bp": 0,
        "longest_private_cluster_sep_le_2bp": 0,
        "foreground_private_gap_count": 0,
        "qc_status": "FAIL",
        "qc_reasons": "",
    }

    reasons = []
    if not seqs:
        reasons.append("empty_alignment")
        row["qc_reasons"] = ";".join(reasons)
        return row
    if not equal_length:
        reasons.append("unequal_alignment_lengths")
        row["qc_reasons"] = ";".join(reasons)
        return row
    if not foreground_present:
        reasons.append("foreground_missing")
        row["qc_reasons"] = ";".join(reasons)
        return row

    fg = seqs[args.foreground]
    others = [seq for name, seq in seqs.items() if name != args.foreground]
    L = len(fg)

    callable_fg = 0
    private_positions = []
    private_gaps = 0

    for i in range(L):
        fbase = fg[i]
        other_bases = [seq[i] for seq in others]

        # Private foreground gaps: foreground gap while sufficiently many
        # nonforeground taxa provide callable bases at the same column.
        nonfg_callable = [b for b in other_bases if b in CALLABLE]
        if fbase == "-" and len(nonfg_callable) >= args.min_nonforeground_callable:
            private_gaps += 1

        if fbase not in CALLABLE:
            continue

        callable_fg += 1

        if len(nonfg_callable) < args.min_nonforeground_callable:
            continue

        counts = Counter(nonfg_callable)
        majority_base, majority_count = counts.most_common(1)[0]
        majority_freq = majority_count / len(nonfg_callable)

        if majority_freq >= args.majority_frequency and fbase != majority_base:
            private_positions.append(i)

    private_count = len(private_positions)
    private_rate = private_count / callable_fg if callable_fg else 0.0

    row.update({
        "foreground_callable_sites": callable_fg,
        "foreground_callable_fraction": callable_fg / L if L else 0.0,
        "private_substitution_count": private_count,
        "private_substitution_rate": private_rate,
        "max_private_substitutions_10bp": max_in_window(private_positions, 10),
        "longest_private_cluster_sep_le_2bp": longest_cluster(private_positions, 2),
        "foreground_private_gap_count": private_gaps,
    })

    if callable_fg == 0:
        reasons.append("no_callable_foreground_bases")
    if private_rate > args.max_private_rate:
        reasons.append("private_substitution_rate_gt_0.08")
    if private_count >= args.max_private_count:
        reasons.append("private_substitution_count_ge_25")
    if row["max_private_substitutions_10bp"] >= args.max_private_in_10bp:
        reasons.append("private_substitutions_10bp_ge_7")
    if private_gaps >= args.max_private_gaps:
        reasons.append("foreground_private_gaps_ge_3")

    row["qc_status"] = "PASS" if not reasons else "FAIL"
    row["qc_reasons"] = ";".join(reasons)
    return row

def main():
    args = parse_args()
    paths = sorted(
        list(args.fasta_dir.glob("*.fa")) +
        list(args.fasta_dir.glob("*.fasta"))
    )
    if not paths:
        raise SystemExit(f"No FASTA files found in {args.fasta_dir}")

    rows = [evaluate_alignment(path, args) for path in paths]
    args.out.parent.mkdir(parents=True, exist_ok=True)

    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    n_pass = sum(r["qc_status"] == "PASS" for r in rows)
    print(f"Alignments tested: {len(rows)}")
    print(f"PASS: {n_pass}")
    print(f"FAIL: {len(rows) - n_pass}")
    print(f"Output: {args.out}")

if __name__ == "__main__":
    main()
