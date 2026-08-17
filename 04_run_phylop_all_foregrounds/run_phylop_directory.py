#!/usr/bin/env python3
"""
Run phyloP SPH subtree tests across FASTA alignments in a directory.

This runner uses the documented PHAST/phyloP SPH subtree interface and writes
one raw phyloP report per alignment. If your historical phyloP build used
additional options, supply them with --extra-arg and document them in README.md.
"""

from __future__ import annotations
import argparse
import subprocess
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("alignment_dir", type=Path)
    p.add_argument("neutral_model", type=Path)
    p.add_argument("subtree", help="Node name identifying the focal lineage/subtree.")
    p.add_argument("output_dir", type=Path)
    p.add_argument("--phyloP", default="phyloP")
    p.add_argument("--msa-format", default="FASTA")
    p.add_argument(
        "--extra-arg",
        action="append",
        default=[],
        help="Additional phyloP argument; may be supplied multiple times."
    )
    return p.parse_args()


def main():
    args = parse_args()
    aln_dir = args.alignment_dir.resolve()
    model = args.neutral_model.resolve()
    outdir = args.output_dir.resolve()

    if not aln_dir.is_dir():
        raise SystemExit(f"Alignment directory not found: {aln_dir}")
    if not model.is_file():
        raise SystemExit(f"Neutral model not found: {model}")

    outdir.mkdir(parents=True, exist_ok=True)
    alignments = sorted(
        list(aln_dir.glob("*.fa")) +
        list(aln_dir.glob("*.fasta"))
    )
    if not alignments:
        raise SystemExit(f"No FASTA alignments found in {aln_dir}")

    failures = []
    for aln in alignments:
        out = outdir / f"{aln.stem}.phyloP.txt"
        cmd = [
            args.phyloP,
            "--method", "SPH",
            "--mode", "CONACC",
            "--subtree", args.subtree,
            "--msa-format", args.msa_format,
            *args.extra_arg,
            str(model),
            str(aln.resolve()),
        ]
        print("[phyloP]", aln.name)
        with out.open("w") as handle:
            proc = subprocess.run(cmd, stdout=handle, stderr=subprocess.PIPE, text=True)
        if proc.returncode != 0:
            failures.append((aln.name, proc.returncode, proc.stderr.strip()))

    print(f"Alignments processed: {len(alignments)}")
    print(f"Failed phyloP runs: {len(failures)}")
    if failures:
        fail_path = outdir / "failed_phylop_runs.tsv"
        with fail_path.open("w") as handle:
            handle.write("alignment\treturn_code\tstderr\n")
            for name, rc, err in failures:
                handle.write(f"{name}\t{rc}\t{err.replace(chr(9), ' ')}\n")
        print(f"Failure log: {fail_path}")


if __name__ == "__main__":
    main()
