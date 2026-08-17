#!/usr/bin/env python3
"""
Extract reference-coordinate intervals from a MAF file and write block-separated FASTA files suitable as input to stitch_blocks.py.
Input interval coordinates are BED-style: 0-based, half-open.

Important:
  * Target intervals MUST use the same reference assembly/coordinate system as the reference row in the MAF.
  * The script currently requires the reference sequence row to be on the '+' strand. This is appropriate for a mouse-reference MAF; if this assertion fails, inspect the MAF rather than silently transforming coordinates.
"""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class Interval:
    chrom: str
    start: int
    end: int
    locus_id: str


@dataclass
class MafSeq:
    src: str
    start: int
    size: int
    strand: str
    src_size: int
    text: str


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--maf", type=Path, required=True,
                   help="Mouse-reference MAF file.")
    p.add_argument("--intervals", type=Path, required=True,
                   help="BED-like file: chrom start end locus_id.")
    p.add_argument("--species-map", type=Path, required=True,
                   help="TSV: regex_pattern<TAB>canonical_name.")
    p.add_argument("--reference-pattern", required=True,
                   help="Regex uniquely matching the mouse-reference MAF src field.")
    p.add_argument("--outdir", type=Path, required=True)
    return p.parse_args()


def read_intervals(path: Path) -> list[Interval]:
    out = []
    with path.open() as fh:
        for n, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                raise ValueError(
                    f"{path}:{n}: expected >=4 tab-separated fields: "
                    "chrom start end locus_id"
                )
            chrom, start, end, locus_id = fields[:4]
            start_i, end_i = int(start), int(end)
            if start_i < 0 or end_i <= start_i:
                raise ValueError(f"{path}:{n}: invalid interval {start_i}-{end_i}")
            out.append(Interval(chrom, start_i, end_i, locus_id))
    return out


def read_species_map(path: Path) -> list[tuple[re.Pattern, str]]:
    mappings = []
    with path.open() as fh:
        for n, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 2:
                raise ValueError(
                    f"{path}:{n}: expected regex_pattern<TAB>canonical_name"
                )
            mappings.append((re.compile(fields[0], re.I), fields[1]))
    return mappings


def canonical_species(src: str, mappings: list[tuple[re.Pattern, str]]) -> str | None:
    hits = [name for pattern, name in mappings if pattern.search(src)]
    hits = list(dict.fromkeys(hits))
    if len(hits) > 1:
        raise ValueError(f"MAF src {src!r} matches multiple species mappings: {hits}")
    return hits[0] if hits else None


def parse_maf_blocks(path: Path) -> Iterable[list[MafSeq]]:
    block: list[MafSeq] = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                if block:
                    yield block
                    block = []
                continue
            if line.startswith("s "):
                f = line.split()
                if len(f) < 7:
                    raise ValueError(f"Malformed MAF s-line: {line}")
                block.append(
                    MafSeq(
                        src=f[1],
                        start=int(f[2]),
                        size=int(f[3]),
                        strand=f[4],
                        src_size=int(f[5]),
                        text=f[6],
                    )
                )
        if block:
            yield block


def src_chrom(src: str) -> str:
    """
    Return a permissive chromosome token from a MAF src string.
    Because Ensembl source naming can differ among exports, matching is also
    attempted by suffix in chrom_matches().
    """
    return src.split(".")[-1]


def chrom_matches(interval_chrom: str, ref_src: str) -> bool:
    want = interval_chrom
    want_no_chr = want[3:] if want.lower().startswith("chr") else want
    token = src_chrom(ref_src)
    token_no_chr = token[3:] if token.lower().startswith("chr") else token

    return (
        token == want
        or token_no_chr == want_no_chr
        or ref_src.endswith("." + want)
        or ref_src.endswith("." + want_no_chr)
        or ref_src.endswith(":" + want)
        or ref_src.endswith(":" + want_no_chr)
    )


def reference_column_slice(ref: MafSeq, overlap_start: int, overlap_end: int) -> tuple[int, int]:
    """
    Convert a genomic interval on a '+'-strand MAF reference row into alignment
    column bounds [col_start, col_end).
    MAF `start` is 0-based. Gaps in the reference text do not advance genomic
    coordinates.
    """
    if ref.strand != "+":
        raise ValueError(
            f"Reference row {ref.src!r} is on strand {ref.strand!r}; "
            "this extractor requires a '+' mouse-reference row"
        )

    gpos = ref.start
    col_start = None
    col_end = None

    for col, base in enumerate(ref.text):
        if base != "-":
            if gpos == overlap_start and col_start is None:
                col_start = col
            if gpos == overlap_end:
                col_end = col
                break
            gpos += 1

    if col_start is None and overlap_start == ref.start + ref.size:
        col_start = len(ref.text)
    if col_end is None:
        if overlap_end == ref.start + ref.size:
            col_end = len(ref.text)
        else:
            # overlap_end may be encountered immediately after the final base represented in the selected block
            gpos2 = ref.start
            for col, base in enumerate(ref.text):
                if base != "-":
                    gpos2 += 1
                    if gpos2 == overlap_end:
                        col_end = col + 1
                        break

    if col_start is None or col_end is None or col_end <= col_start:
        raise ValueError(
            f"Could not map reference coordinates {overlap_start}-{overlap_end} "
            f"onto MAF row {ref.src}:{ref.start}-{ref.start + ref.size}"
        )
    return col_start, col_end


def write_block(handle, seqs: dict[str, str]) -> None:
    for name, seq in seqs.items():
        handle.write(f">{name}\n{seq}\n")
    handle.write("\n")


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    intervals = read_intervals(args.intervals)
    mappings = read_species_map(args.species_map)
    ref_re = re.compile(args.reference_pattern, re.I)

    by_chrom: dict[str, list[Interval]] = defaultdict(list)
    for iv in intervals:
        by_chrom[iv.chrom].append(iv)
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda x: (x.start, x.end))

    # Open lazily so loci with no overlap do not get empty files
    handles = {}
    block_counts = defaultdict(int)

    try:
        for block in parse_maf_blocks(args.maf):
            refs = [row for row in block if ref_re.search(row.src)]
            if not refs:
                continue
            if len(refs) != 1:
                raise ValueError(
                    "Reference pattern matched multiple rows in one MAF block: "
                    + ", ".join(r.src for r in refs)
                )
            ref = refs[0]
            ref_start = ref.start
            ref_end = ref.start + ref.size

            candidate_intervals = []
            for chrom, ivs in by_chrom.items():
                if chrom_matches(chrom, ref.src):
                    candidate_intervals.extend(ivs)

            for iv in candidate_intervals:
                if iv.end <= ref_start or iv.start >= ref_end:
                    continue

                ov_start = max(iv.start, ref_start)
                ov_end = min(iv.end, ref_end)
                c0, c1 = reference_column_slice(ref, ov_start, ov_end)

                seqs = {}
                for row in block:
                    name = canonical_species(row.src, mappings)
                    if name is None:
                        continue
                    frag = row.text[c0:c1]
                    if name in seqs:
                        raise ValueError(
                            f"Duplicate focal species {name!r} in one MAF block "
                            f"for locus {iv.locus_id}; inspect duplicated mappings."
                        )
                    seqs[name] = frag

                if not seqs:
                    continue

                out_path = args.outdir / f"{iv.locus_id}.fa"
                if iv.locus_id not in handles:
                    handles[iv.locus_id] = out_path.open("w")
                write_block(handles[iv.locus_id], seqs)
                block_counts[iv.locus_id] += 1

    finally:
        for h in handles.values():
            h.close()

    missing = [iv.locus_id for iv in intervals if block_counts[iv.locus_id] == 0]
    print(f"Intervals supplied: {len(intervals)}")
    print(f"Loci with >=1 extracted block: {len(block_counts)}")
    print(f"Loci with no extracted block: {len(missing)}")
    if missing:
        missing_path = args.outdir / "loci_without_maf_overlap.txt"
        missing_path.write_text("\n".join(sorted(set(missing))) + "\n")
        print(f"Missing-locus list: {missing_path}")


if __name__ == "__main__":
    main()
