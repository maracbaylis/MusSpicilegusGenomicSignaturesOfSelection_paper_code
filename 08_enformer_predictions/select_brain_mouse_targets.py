#!/usr/bin/env python3
"""Select brain/neural mouse Enformer target metadata rows."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


DEFAULT_KEYWORDS = [
    "brain",
    "cortex",
    "cortical",
    "forebrain",
    "hindbrain",
    "midbrain",
    "hypothalam",
    "hippocamp",
    "cerebell",
    "neuron",
    "neural",
    "embryonic day",
    "e14",
    "e15",
    "e16",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--keywords", nargs="*", default=DEFAULT_KEYWORDS)
    args = parser.parse_args()

    targets = pd.read_csv(args.targets, sep="\t").reset_index(names="mouse_output_channel")
    pattern = "|".join(args.keywords)
    text = targets.fillna("").astype(str).agg(" ".join, axis=1).str.lower()
    selected = targets[text.str.contains(pattern, regex=True)].copy()
    selected.insert(0, "selection_keywords", ",".join(args.keywords))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.out, index=False)
    print(f"Selected {len(selected)} / {len(targets)} targets")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
