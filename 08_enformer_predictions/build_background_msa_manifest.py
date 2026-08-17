#!/usr/bin/env python3
"""Build an Enformer manifest from a directory of background cCRE MSAs."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_spicilegus_enformer_input_manifest import (  # noqa: E402
    SPIC_HEADER,
    clean_seq,
    infer_changes,
    mouse_consensus,
    parse_fasta,
    ungap,
    wrap,
)


FILENAME_RE = re.compile(r"^(?P<ccre_id>[^_]+)__+(?P<chrom>chr[^_]+)_(?P<start>\d+)_(?P<end>\d+)\.fa$")


def read_exclude(path: Path | None) -> set[str]:
    if path is None:
        return set()
    ids: set[str] = set()
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            ids.add(line.split(",")[0].split("\t")[0])
    return ids


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--msa-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--exclude-id-list", type=Path)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    exclude = read_exclude(args.exclude_id_list)

    manifest_path = args.outdir / "background_enformer_manifest.csv"
    variants_path = args.outdir / "background_alignment_changes.tsv"
    fasta_path = args.outdir / "background_ccre_alleles.fa"
    excluded_path = args.outdir / "background_excluded_or_failed.csv"
    summary_path = args.outdir / "background_manifest_summary.txt"

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

    manifest_rows: list[dict[str, object]] = []
    variants_rows: list[dict[str, object]] = []
    fasta_records: list[tuple[str, str]] = []
    excluded_rows: list[dict[str, str]] = []

    for msa_path in sorted(args.msa_dir.glob("*.fa")):
        match = FILENAME_RE.match(msa_path.name)
        if not match:
            excluded_rows.append({"cCRE_id": msa_path.stem, "reason": "filename_parse_failed"})
            continue
        ccre_id = match.group("ccre_id")
        if ccre_id in exclude:
            excluded_rows.append({"cCRE_id": ccre_id, "reason": "excluded_id"})
            continue
        chrom = match.group("chrom")
        start = int(match.group("start"))
        end = int(match.group("end"))

        seqs = parse_fasta(msa_path)
        spic_aln = seqs.get(SPIC_HEADER, "")
        mouse_aln = mouse_consensus(seqs)
        if not spic_aln or not mouse_aln:
            excluded_rows.append({"cCRE_id": ccre_id, "reason": "missing_spic_or_mouse_alignment"})
            continue

        mouse_seq = clean_seq(ungap(mouse_aln))
        spic_seq = clean_seq(ungap(spic_aln))
        subs, ins, dels, events = infer_changes(mouse_aln, spic_aln)
        has_indel = bool(ins or dels or len(mouse_seq) != len(spic_seq))
        midpoint = (start + end) // 2

        manifest_rows.append(
            {
                "cCRE_id": ccre_id,
                "nearest_gene_name": "",
                "cCRE_classification": "background",
                "chrom": chrom,
                "start": start,
                "end": end,
                "midpoint": midpoint,
                "ccre_len_bp": end - start,
                "q_p_accel_subtree_given_total": "",
                "p_accel_subtree_given_total": "",
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
        writer = csv.DictWriter(
            handle,
            fieldnames=["cCRE_id", "substitution_count", "insertion_count", "deletion_count", "alignment_change_events"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(variants_rows)

    with fasta_path.open("w") as handle:
        for header, seq in fasta_records:
            handle.write(f">{header}\n{wrap(seq)}\n")

    with excluded_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["cCRE_id", "reason"])
        writer.writeheader()
        writer.writerows(excluded_rows)

    indel_count = sum(bool(row["has_indel"]) for row in manifest_rows)
    with summary_path.open("w") as handle:
        handle.write("background Enformer input manifest summary\n")
        handle.write(f"source_msa_dir: {args.msa_dir}\n")
        handle.write(f"msa_files: {len(list(args.msa_dir.glob('*.fa')))}\n")
        handle.write(f"manifest_rows: {len(manifest_rows)}\n")
        handle.write(f"excluded_or_failed_rows: {len(excluded_rows)}\n")
        handle.write(f"rows_with_spic_indels_or_length_delta: {indel_count}\n")
        handle.write(f"substitution_only_rows: {len(manifest_rows) - indel_count}\n")

    print(f"Wrote {manifest_path}")
    print(f"Wrote {fasta_path}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
