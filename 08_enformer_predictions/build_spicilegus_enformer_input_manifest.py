#!/usr/bin/env python3
"""Build an Enformer-ready manifest for spicilegus accelerated cCRE hits.

This prepares the Pollard/SuPreMo-style input layer for allele replacement:
mouse genomic context is handled later from a full mouse reference genome, while
this script packages the accelerated cCRE coordinates, MSA-derived spicilegus
ortholog sequence, and alignment-derived mouse consensus sequence.
"""

from __future__ import annotations

import argparse
import csv
import glob
import re
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_HITS = (
    ROOT
    / "phyloP_MSA_QC"
    / "manual_msa_inspection_tables"
    / "spicilegus_phyloP_accelerated_285_full_stats.csv"
)
DEFAULT_OUTDIR = ROOT / "enformer_spicilegus_phyloP" / "results"
DEFAULT_MSA_DIRS = [
    ROOT / "cCRE_reulatory_regions_fasta_procesing_for_MSA" / "fa_cCRE_brain.8taxa.stitched",
    ROOT / "cCRE_reulatory_regions_fasta_procesing_for_MSA" / "fa_cCRE_brain.8taxa",
    ROOT / "phyloP_MSA_QC" / "more_accelerated_than_expected_given_total_fastas" / "spicilegus",
    ROOT / "phyloP_MSA_QC" / "more_accelerated_than_expected_given_total_PASS_fastas" / "spicilegus",
]
MOUSE_PREFIX = "mus_musculus"
SPIC_HEADER = "mus_spicilegus"


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_id_list(path: Path) -> set[str]:
    ids: set[str] = set()
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            ids.add(line.split(",")[0].split("\t")[0])
    return ids


def parse_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    header: str | None = None
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].split()[0]
                seqs.setdefault(header, [])
            elif header is not None:
                seqs[header].append(line.upper())
    return {name: "".join(parts) for name, parts in seqs.items()}


def wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def ungap(seq: str) -> str:
    return seq.replace("-", "").upper()


def clean_seq(seq: str) -> str:
    return re.sub("[^ACGTN]", "N", seq.upper())


def mouse_consensus(seqs: dict[str, str]) -> str:
    mouse = [seq for name, seq in seqs.items() if name.startswith(MOUSE_PREFIX)]
    if not mouse:
        return ""
    max_len = max(len(seq) for seq in mouse)
    consensus: list[str] = []
    for idx in range(max_len):
        bases = [seq[idx] for seq in mouse if idx < len(seq)]
        counts = Counter(base for base in bases if base in "ACGTN")
        if counts:
            consensus.append(counts.most_common(1)[0][0])
        else:
            consensus.append("-")
    return "".join(consensus)


def find_msa(ccre_id: str, msa_dirs: list[Path]) -> Path | None:
    matches: list[Path] = []
    for msa_dir in msa_dirs:
        matches.extend(Path(p) for p in glob.glob(str(msa_dir / f"{ccre_id}__*.fa")))
    if not matches:
        return None
    # Prefer stitched 8-taxon alignments because they represent the full cCRE
    # interval most directly; otherwise use the shortest path as a stable tie-break.
    return sorted(matches, key=lambda p: ("stitched" not in str(p), len(str(p)), str(p)))[0]


def infer_changes(mouse_aln: str, spic_aln: str) -> tuple[int, int, int, str]:
    substitutions = insertions = deletions = compared = 0
    events: list[str] = []
    mouse_pos = 0
    for m_base, s_base in zip(mouse_aln, spic_aln):
        if m_base != "-":
            mouse_pos += 1
        if m_base == "-" and s_base != "-":
            insertions += 1
            events.append(f"{mouse_pos}:ins:{s_base}")
        elif m_base != "-" and s_base == "-":
            deletions += 1
            events.append(f"{mouse_pos}:del:{m_base}")
        elif m_base in "ACGTN" and s_base in "ACGTN":
            compared += 1
            if m_base != s_base:
                substitutions += 1
                events.append(f"{mouse_pos}:{m_base}>{s_base}")
    return substitutions, insertions, deletions, ";".join(events)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hits", type=Path, default=DEFAULT_HITS)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--explicit-id-list",
        type=Path,
        help="Optional newline/CSV/TSV list of cCRE IDs to restrict the hit set, e.g. the final 283.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        help="Optional number of top rows to keep after sorting by q_p_accel_subtree_given_total.",
    )
    args = parser.parse_args()

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    rows = read_table(args.hits)
    if args.explicit_id_list:
        keep = read_id_list(args.explicit_id_list)
        rows = [row for row in rows if row["cCRE_id"] in keep]

    rows = sorted(
        rows,
        key=lambda row: (
            float(row.get("q_p_accel_subtree_given_total") or 1),
            row.get("chrom", ""),
            int(float(row.get("start") or 0)),
        ),
    )
    if args.limit:
        rows = rows[: args.limit]

    manifest_path = outdir / "spicilegus_phyloP_hits_enformer_manifest.csv"
    variants_path = outdir / "spicilegus_phyloP_hits_alignment_changes.tsv"
    fasta_path = outdir / "spicilegus_phyloP_hits_ccre_alleles.fa"
    missing_path = outdir / "spicilegus_phyloP_hits_missing_msa.csv"
    summary_path = outdir / "spicilegus_phyloP_hits_manifest_summary.txt"

    manifest_fields = [
        "cCRE_id",
        "nearest_gene_name",
        "cCRE_classification",
        "chrom",
        "start",
        "end",
        "midpoint",
        "ccre_len_bp",
        "q_p_accel_subtree_given_total",
        "p_accel_subtree_given_total",
        "msa_path",
        "mouse_alignment_consensus_len",
        "spic_ortholog_len",
        "length_delta_spic_minus_mouse",
        "substitution_count",
        "insertion_count",
        "deletion_count",
        "has_indel",
        "ready_for_full_enformer_window",
    ]

    missing: list[dict[str, str]] = []
    manifest_rows: list[dict[str, object]] = []
    variants_rows: list[dict[str, object]] = []
    fasta_records: list[tuple[str, str]] = []

    for row in rows:
        ccre_id = row["cCRE_id"]
        msa_path = find_msa(ccre_id, DEFAULT_MSA_DIRS)
        chrom = row["chrom"]
        start = int(float(row["start"]))
        end = int(float(row["end"]))
        midpoint = (start + end) // 2
        if msa_path is None:
            missing.append({"cCRE_id": ccre_id, "reason": "no_matching_msa"})
            continue

        seqs = parse_fasta(msa_path)
        spic_aln = seqs.get(SPIC_HEADER, "")
        mouse_aln = mouse_consensus(seqs)
        if not spic_aln or not mouse_aln:
            missing.append({"cCRE_id": ccre_id, "reason": "missing_spic_or_mouse_alignment"})
            continue

        mouse_seq = clean_seq(ungap(mouse_aln))
        spic_seq = clean_seq(ungap(spic_aln))
        subs, ins, dels, events = infer_changes(mouse_aln, spic_aln)
        has_indel = bool(ins or dels or len(mouse_seq) != len(spic_seq))

        manifest_rows.append(
            {
                "cCRE_id": ccre_id,
                "nearest_gene_name": row.get("nearest_gene_name", ""),
                "cCRE_classification": row.get("cCRE_classification", ""),
                "chrom": chrom,
                "start": start,
                "end": end,
                "midpoint": midpoint,
                "ccre_len_bp": end - start,
                "q_p_accel_subtree_given_total": row.get("q_p_accel_subtree_given_total", ""),
                "p_accel_subtree_given_total": row.get("p_accel_subtree_given_total", ""),
                "msa_path": str(msa_path),
                "mouse_alignment_consensus_len": len(mouse_seq),
                "spic_ortholog_len": len(spic_seq),
                "length_delta_spic_minus_mouse": len(spic_seq) - len(mouse_seq),
                "substitution_count": subs,
                "insertion_count": ins,
                "deletion_count": dels,
                "has_indel": has_indel,
                "ready_for_full_enformer_window": "needs_mouse_reference_genome",
            }
        )
        variants_rows.append(
            {
                "cCRE_id": ccre_id,
                "substitution_count": subs,
                "insertion_count": ins,
                "deletion_count": dels,
                "alignment_change_events": events,
            }
        )
        fasta_records.append((f"{ccre_id}|allele=mouse_alignment_consensus", mouse_seq))
        fasta_records.append((f"{ccre_id}|allele=spicilegus_ortholog", spic_seq))

    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_fields)
        writer.writeheader()
        writer.writerows(manifest_rows)

    with variants_path.open("w", newline="") as handle:
        fieldnames = ["cCRE_id", "substitution_count", "insertion_count", "deletion_count", "alignment_change_events"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(variants_rows)

    with fasta_path.open("w") as handle:
        for header, seq in fasta_records:
            handle.write(f">{header}\n{wrap(seq)}\n")

    with missing_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["cCRE_id", "reason"])
        writer.writeheader()
        writer.writerows(missing)

    indel_count = sum(bool(row["has_indel"]) for row in manifest_rows)
    with summary_path.open("w") as handle:
        handle.write("spicilegus phyloP Enformer input manifest summary\n")
        handle.write(f"source_hits: {args.hits}\n")
        handle.write(f"input_rows_after_filters: {len(rows)}\n")
        handle.write(f"manifest_rows_with_msa: {len(manifest_rows)}\n")
        handle.write(f"missing_msa_rows: {len(missing)}\n")
        handle.write(f"rows_with_spic_indels_or_length_delta: {indel_count}\n")
        handle.write(f"substitution_only_rows: {len(manifest_rows) - indel_count}\n")
        handle.write("\nNext step: add a full mm10/mm39 reference FASTA and construct 393216-bp Enformer windows centered on each cCRE midpoint.\n")

    print(f"Wrote {manifest_path}")
    print(f"Wrote {variants_path}")
    print(f"Wrote {fasta_path}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
