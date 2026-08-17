#!/usr/bin/env python3
"""Construct full-length Enformer reference and spicilegus-edited windows."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ENFORMER_INPUT_LEN = 393_216


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    header: str | None = None
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].split()[0]
                records.setdefault(header, [])
            elif header is not None:
                records[header].append(line.upper())
    return {header: "".join(parts) for header, parts in records.items()}


def wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


class IndexedFasta:
    def __init__(self, fasta: Path):
        self.fasta = fasta
        self.fai = fasta.with_suffix(fasta.suffix + ".fai")
        if not self.fai.exists():
            self._build_index()
        self.index = self._read_index()
        self.handle = fasta.open("rb")

    def _build_index(self) -> None:
        rows: list[tuple[str, int, int, int, int]] = []
        with self.fasta.open("rb") as handle:
            name = None
            length = 0
            seq_offset = 0
            line_bases = 0
            line_width = 0
            while True:
                offset = handle.tell()
                line = handle.readline()
                if not line:
                    if name is not None:
                        rows.append((name, length, seq_offset, line_bases, line_width))
                    break
                if line.startswith(b">"):
                    if name is not None:
                        rows.append((name, length, seq_offset, line_bases, line_width))
                    name = line[1:].decode().strip().split()[0]
                    length = 0
                    seq_offset = handle.tell()
                    line_bases = 0
                    line_width = 0
                else:
                    stripped = line.rstrip(b"\r\n")
                    if name is None or not stripped:
                        continue
                    if line_bases == 0:
                        line_bases = len(stripped)
                        line_width = len(line)
                    length += len(stripped)
        with self.fai.open("w") as out:
            for row in rows:
                out.write("\t".join(map(str, row)) + "\n")

    def _read_index(self) -> dict[str, tuple[int, int, int, int]]:
        index: dict[str, tuple[int, int, int, int]] = {}
        with self.fai.open() as handle:
            for line in handle:
                name, length, offset, line_bases, line_width = line.rstrip("\n").split("\t")[:5]
                index[name] = (int(length), int(offset), int(line_bases), int(line_width))
        return index

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Fetch 0-based half-open sequence, padding with N outside contig bounds."""
        if chrom not in self.index:
            raise KeyError(f"{chrom} not present in {self.fasta}")
        chrom_len, offset, line_bases, line_width = self.index[chrom]
        left_pad = max(0, -start)
        right_pad = max(0, end - chrom_len)
        query_start = max(0, start)
        query_end = min(chrom_len, end)
        if query_end <= query_start:
            return "N" * (end - start)

        pieces: list[str] = []
        pos = query_start
        while pos < query_end:
            line_idx = pos // line_bases
            in_line = pos % line_bases
            take = min(query_end - pos, line_bases - in_line)
            self.handle.seek(offset + line_idx * line_width + in_line)
            pieces.append(self.handle.read(take).decode().upper())
            pos += take
        return "N" * left_pad + "".join(pieces) + "N" * right_pad


def load_replacement_alleles(path: Path, allele_name: str) -> dict[str, str]:
    records = read_fasta(path)
    alleles: dict[str, str] = {}
    for header, seq in records.items():
        ccre_id, *parts = header.split("|")
        if any(part == f"allele={allele_name}" for part in parts):
            alleles[ccre_id] = seq
    return alleles


def edited_window(
    fasta: IndexedFasta,
    chrom: str,
    ccre_start: int,
    ccre_end: int,
    replacement: str,
    input_len: int,
) -> tuple[str, int, int, int, int]:
    midpoint = (ccre_start + ccre_end) // 2
    ref_window_start = midpoint - input_len // 2
    left_len = ccre_start - ref_window_start
    right_len = input_len - left_len - len(replacement)
    if right_len < 0:
        raise ValueError(
            f"replacement length {len(replacement)} is too long for {chrom}:{ccre_start}-{ccre_end}"
        )
    left = fasta.fetch(chrom, ccre_start - left_len, ccre_start)
    right = fasta.fetch(chrom, ccre_end, ccre_end + right_len)
    edited = left + replacement + right
    return edited, ccre_start - left_len, ccre_end + right_len, left_len, right_len


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--allele-fasta", type=Path, required=True)
    parser.add_argument("--replacement-allele", default="spicilegus_ortholog")
    parser.add_argument("--replacement-label", default="spicilegus_replacement")
    parser.add_argument("--mouse-fasta", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--input-len", type=int, default=ENFORMER_INPUT_LEN)
    parser.add_argument(
        "--coordinate-mode",
        choices=["half-open", "inclusive-end"],
        default="half-open",
        help="Use inclusive-end if the cCRE table end coordinate is 1-based/inclusive.",
    )
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    fasta_dir = args.outdir / "fasta"
    fasta_dir.mkdir(exist_ok=True)

    rows = read_manifest(args.manifest)
    if args.limit:
        rows = rows[: args.limit]
    replacement_alleles = load_replacement_alleles(args.allele_fasta, args.replacement_allele)
    mouse = IndexedFasta(args.mouse_fasta)

    summary_fields = [
        "cCRE_id",
        "chrom",
        "ccre_start",
        "ccre_end_for_window",
        "input_len",
        "reference_window_start",
        "reference_window_end",
        "edited_window_start",
        "edited_window_end",
        "spic_replacement_len",
        "reference_fasta",
        "spic_edited_fasta",
    ]
    summary_rows: list[dict[str, object]] = []

    for row in rows:
        ccre_id = row["cCRE_id"]
        chrom = row["chrom"]
        ccre_start = int(row["start"])
        ccre_end = int(row["end"]) + (1 if args.coordinate_mode == "inclusive-end" else 0)
        midpoint = (ccre_start + ccre_end) // 2
        ref_start = midpoint - args.input_len // 2
        ref_end = ref_start + args.input_len
        ref_seq = mouse.fetch(chrom, ref_start, ref_end)
        replacement = replacement_alleles[ccre_id]
        edit_seq, edit_start, edit_end, _left_len, _right_len = edited_window(
            mouse, chrom, ccre_start, ccre_end, replacement, args.input_len
        )
        if len(ref_seq) != args.input_len or len(edit_seq) != args.input_len:
            raise RuntimeError(f"bad length for {ccre_id}: ref={len(ref_seq)} edit={len(edit_seq)}")

        ref_path = fasta_dir / f"{ccre_id}.mouse_reference.enformer.fa"
        edit_path = fasta_dir / f"{ccre_id}.{args.replacement_label}.enformer.fa"
        ref_header = f"{ccre_id}|allele=mouse_reference|{chrom}:{ref_start}-{ref_end}"
        edit_header = f"{ccre_id}|allele={args.replacement_label}|{chrom}:{edit_start}-{edit_end}|replacement_len={len(replacement)}"
        ref_path.write_text(f">{ref_header}\n{wrap(ref_seq)}\n")
        edit_path.write_text(f">{edit_header}\n{wrap(edit_seq)}\n")

        summary_rows.append(
            {
                "cCRE_id": ccre_id,
                "chrom": chrom,
                "ccre_start": ccre_start,
                "ccre_end_for_window": ccre_end,
                "input_len": args.input_len,
                "reference_window_start": ref_start,
                "reference_window_end": ref_end,
                "edited_window_start": edit_start,
                "edited_window_end": edit_end,
                "spic_replacement_len": len(replacement),
                "reference_fasta": str(ref_path),
                "spic_edited_fasta": str(edit_path),
            }
        )

    summary_path = args.outdir / "enformer_window_manifest.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"Wrote {summary_path}")
    print(f"Wrote {len(summary_rows) * 2} FASTA files in {fasta_dir}")


if __name__ == "__main__":
    main()
