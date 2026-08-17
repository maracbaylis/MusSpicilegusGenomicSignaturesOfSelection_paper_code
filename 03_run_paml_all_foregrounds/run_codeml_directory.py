
"""
Run PAML codeml across all PHYLIP alignments in a directory.

The runner uses a codeml control-file template and replaces only:
  - seqfile
  - treefile
  - outfile

This allows the same script to run either the null one-ratio model or a
foreground two-ratio model by supplying the corresponding template and tree.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run codeml on every .phy alignment in a directory."
    )
    parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing PHYLIP codon alignments (*.phy).",
    )
    parser.add_argument(
        "output_dir",
        type=Path,
        help="Directory for codeml output files.",
    )
    parser.add_argument(
        "tree_file",
        type=Path,
        help="Tree file used by codeml.",
    )
    parser.add_argument(
        "template_ctl",
        type=Path,
        help="codeml control-file template (e.g. null or model2).",
    )
    parser.add_argument(
        "--codeml",
        default="codeml",
        help="Path/name of codeml executable. Default: codeml",
    )
    parser.add_argument(
        "--suffix",
        default=None,
        help=(
            "Output suffix. Defaults to '.null.codeml_out' for a template whose "
            "name contains 'null', otherwise '.model2.codeml_out'."
        ),
    )
    return parser.parse_args()


def output_suffix(template_ctl: Path, explicit_suffix: str | None) -> str:
    if explicit_suffix:
        return explicit_suffix
    if "null" in template_ctl.name.lower():
        return ".null.codeml_out"
    return ".model2.codeml_out"


def write_control_file(
    template_ctl: Path,
    control_path: Path,
    seqfile: Path,
    treefile: Path,
    outfile: Path,
) -> None:
    """
    Copy the template while replacing seqfile, treefile, and outfile entries.
    """
    seen = {"seqfile": False, "treefile": False, "outfile": False}

    with template_ctl.open() as source, control_path.open("w") as dest:
        for line in source:
            stripped = line.strip()

            if stripped.startswith("seqfile"):
                dest.write(f"seqfile = {seqfile.resolve()}\n")
                seen["seqfile"] = True
            elif stripped.startswith("treefile"):
                dest.write(f"treefile = {treefile.resolve()}\n")
                seen["treefile"] = True
            elif stripped.startswith("outfile"):
                dest.write(f"outfile = {outfile.resolve()}\n")
                seen["outfile"] = True
            else:
                dest.write(line)

    missing = [key for key, found in seen.items() if not found]
    if missing:
        raise ValueError(
            f"{template_ctl} is missing required control fields: {', '.join(missing)}"
        )


def main() -> None:
    args = parse_args()

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    tree_file = args.tree_file.resolve()
    template_ctl = args.template_ctl.resolve()

    if not input_dir.is_dir():
        raise SystemExit(f"Input directory not found: {input_dir}")
    if not tree_file.is_file():
        raise SystemExit(f"Tree file not found: {tree_file}")
    if not template_ctl.is_file():
        raise SystemExit(f"Template control file not found: {template_ctl}")

    output_dir.mkdir(parents=True, exist_ok=True)

    alignments = sorted(input_dir.glob("*.phy"))
    if not alignments:
        raise SystemExit(f"No .phy alignments found in {input_dir}")

    suffix = output_suffix(template_ctl, args.suffix)
    control_path = output_dir / "codeml_current.ctl"

    failed = []

    for alignment in alignments:
        gene = alignment.stem
        outfile = output_dir / f"{gene}{suffix}"

        write_control_file(
            template_ctl=template_ctl,
            control_path=control_path,
            seqfile=alignment,
            treefile=tree_file,
            outfile=outfile,
        )

        print(f"[codeml] {gene}")
        result = subprocess.run(
            [args.codeml, str(control_path)],
            check=False,
        )

        if result.returncode != 0:
            failed.append((gene, result.returncode))

    if control_path.exists():
        control_path.unlink()

    print(f"Alignments processed: {len(alignments)}")
    print(f"Failed codeml runs: {len(failed)}")

    if failed:
        failed_path = output_dir / "failed_codeml_runs.tsv"
        with failed_path.open("w") as handle:
            handle.write("gene\treturn_code\n")
            for gene, return_code in failed:
                handle.write(f"{gene}\t{return_code}\n")
        print(f"Failure log: {failed_path}")


if __name__ == "__main__":
    main()

