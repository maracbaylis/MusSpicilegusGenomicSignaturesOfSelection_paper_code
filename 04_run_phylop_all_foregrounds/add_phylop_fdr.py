#!/usr/bin/env python3
"""Add BH/FDR q-value columns to phyloP CSV result files."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def parse_float(value: str) -> float | None:
    value = (value or "").strip()
    if not value:
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    if not math.isfinite(parsed):
        return None
    return min(max(parsed, 0.0), 1.0)


def bh_adjust(values: list[float | None]) -> list[float | None]:
    indexed = [(idx, value) for idx, value in enumerate(values) if value is not None]
    n = len(indexed)
    adjusted: list[float | None] = [None] * len(values)
    if n == 0:
        return adjusted

    indexed.sort(key=lambda item: item[1])
    running_min = 1.0
    for rank_from_end, (idx, pvalue) in enumerate(reversed(indexed), start=1):
        rank = n - rank_from_end + 1
        qvalue = min(running_min, pvalue * n / rank)
        running_min = qvalue
        adjusted[idx] = min(qvalue, 1.0)
    return adjusted


def format_q(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.10g}"


def insert_q_columns(fieldnames: list[str], p_columns: list[str]) -> list[str]:
    out: list[str] = []
    p_set = set(p_columns)
    for name in fieldnames:
        out.append(name)
        if name in p_set:
            out.append(f"q_{name[2:]}")
    return out


def process_csv(path: Path, output_dir: Path) -> dict[str, object]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {path}")
        fieldnames = reader.fieldnames

    p_columns = [name for name in fieldnames if name.startswith("p_")]
    q_by_column: dict[str, list[float | None]] = {}
    p_by_column: dict[str, list[float | None]] = {}

    for p_column in p_columns:
        values = [parse_float(row.get(p_column, "")) for row in rows]
        p_by_column[p_column] = values
        q_by_column[p_column] = bh_adjust(values)

    output_fieldnames = insert_q_columns(fieldnames, p_columns)
    output_path = output_dir / f"{path.stem}_with_fdr.csv"
    output_dir.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fieldnames)
        writer.writeheader()
        for row_idx, row in enumerate(rows):
            out_row: dict[str, str] = {}
            for field in output_fieldnames:
                if field.startswith("q_"):
                    p_column = f"p_{field[2:]}"
                    out_row[field] = format_q(q_by_column[p_column][row_idx])
                else:
                    out_row[field] = row.get(field, "")
            writer.writerow(out_row)

    summary: dict[str, object] = {
        "input_file": str(path),
        "output_file": str(output_path),
        "rows": len(rows),
    }
    for p_column in p_columns:
        values = p_by_column[p_column]
        qvalues = q_by_column[p_column]
        valid = [value for value in values if value is not None]
        significant = [q for q in qvalues if q is not None and q <= 0.05]
        summary[f"{p_column}_n_tested"] = len(valid)
        summary[f"q_{p_column[2:]}_le_0.05"] = len(significant)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    inputs = sorted(args.input_dir.glob("*.csv"))
    if not inputs:
        raise SystemExit(f"No CSV files found in {args.input_dir}")

    summaries = [process_csv(path, args.output_dir) for path in inputs]
    summary_path = args.output_dir / "fdr_summary.csv"
    summary_fields: list[str] = []
    for summary in summaries:
        for key in summary:
            if key not in summary_fields:
                summary_fields.append(key)

    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows(summaries)

    print(f"Wrote {len(summaries)} FDR-adjusted CSVs to {args.output_dir}")
    print(f"Wrote summary to {summary_path}")


if __name__ == "__main__":
    main()
