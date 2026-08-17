#!/usr/bin/env python3
"""Run Enformer predictions for paired mouse-reference/spicilegus windows."""

from __future__ import annotations

import argparse
import csv
import os
import ssl
from pathlib import Path

import numpy as np
import pandas as pd
import tensorflow as tf


MODEL_URL = "https://tfhub.dev/deepmind/enformer/1"
OUTPUT_START_IN_INPUT = (393_216 - 114_688) // 2
BIN_SIZE = 128


def read_fasta_one(path: Path) -> str:
    return "".join(line.strip().upper() for line in path.open() if not line.startswith(">"))


def one_hot(seq: str) -> np.ndarray:
    arr = np.zeros((len(seq), 4), dtype=np.float32)
    mapping = {"A": 0, "C": 1, "G": 2, "T": 3}
    for idx, base in enumerate(seq):
        col = mapping.get(base)
        if col is not None:
            arr[idx, col] = 1.0
    return arr


def bins_for_interval(rel_start: int, rel_end: int, output_bins: int) -> list[int]:
    start_bin = int(np.floor((rel_start - OUTPUT_START_IN_INPUT) / BIN_SIZE))
    end_bin = int(np.ceil((rel_end - OUTPUT_START_IN_INPUT) / BIN_SIZE))
    start_bin = max(0, start_bin)
    end_bin = min(output_bins, end_bin)
    return list(range(start_bin, end_bin))


def predict_mouse(model, seq: str) -> np.ndarray:
    inputs = tf.convert_to_tensor(one_hot(seq)[None, :, :])
    outputs = model.predict_on_batch(inputs)
    if isinstance(outputs, dict):
        mouse = outputs["mouse"]
    else:
        mouse = outputs
    return np.asarray(mouse)[0]


def predict_mouse_pair(model, ref_seq: str, edit_seq: str) -> tuple[np.ndarray, np.ndarray]:
    inputs = tf.convert_to_tensor(np.stack([one_hot(ref_seq), one_hot(edit_seq)], axis=0))
    outputs = model.predict_on_batch(inputs)
    if isinstance(outputs, dict):
        mouse = outputs["mouse"]
    else:
        mouse = outputs
    arr = np.asarray(mouse)
    return arr[0], arr[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--window-manifest", type=Path, required=True)
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--model", default=MODEL_URL)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--target-limit", type=int)
    parser.add_argument("--write-bin-summaries", action="store_true")
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from partial CSV outputs in --outdir, skipping completed cCRE IDs.",
    )
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    windows = pd.read_csv(args.window_manifest)
    if args.limit:
        windows = windows.head(args.limit).copy()

    targets = pd.read_csv(args.targets)
    if "mouse_output_channel" not in targets.columns:
        targets = targets.reset_index(names="mouse_output_channel")
    if args.target_limit:
        targets = targets.head(args.target_limit).copy()

    try:
        import certifi

        os.environ.setdefault("SSL_CERT_FILE", certifi.where())
        os.environ.setdefault("REQUESTS_CA_BUNDLE", certifi.where())
        ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=certifi.where())
    except Exception:
        pass

    import tensorflow_hub as hub

    print(f"Loading Enformer model: {args.model}", flush=True)
    model = hub.load(args.model).model
    print("Model loaded", flush=True)

    target_partial = args.outdir / "enformer_local_delta_by_target.partial.csv"
    ccre_partial = args.outdir / "enformer_local_delta_by_ccre.partial.csv"
    long_rows: list[dict[str, object]] = []
    ccre_rows: list[dict[str, object]] = []
    bin_rows: list[dict[str, object]] = []

    completed: set[str] = set()
    if args.resume and target_partial.exists() and ccre_partial.exists():
        prev_targets = pd.read_csv(target_partial)
        prev_ccres = pd.read_csv(ccre_partial)
        if "cCRE_id" in prev_targets and "cCRE_id" in prev_ccres:
            completed = set(prev_ccres["cCRE_id"].astype(str))
            long_rows = prev_targets.to_dict("records")
            ccre_rows = prev_ccres.to_dict("records")
            print(f"Resuming with {len(completed)} completed cCREs", flush=True)

    for row_idx, row in windows.iterrows():
        ccre_id = row["cCRE_id"]
        if ccre_id in completed:
            continue
        print(f"[{row_idx + 1}/{len(windows)}] {ccre_id}", flush=True)
        ref_seq = read_fasta_one(Path(row["reference_fasta"]))
        edit_seq = read_fasta_one(Path(row["spic_edited_fasta"]))
        ref_pred, edit_pred = predict_mouse_pair(model, ref_seq, edit_seq)

        output_bins, output_tracks = ref_pred.shape
        rel_ref_start = int(row["ccre_start"]) - int(row["reference_window_start"])
        rel_ref_end = int(row["ccre_end_for_window"]) - int(row["reference_window_start"])
        rel_edit_start = int(row["ccre_start"]) - int(row["edited_window_start"])
        rel_edit_end = rel_edit_start + int(row["spic_replacement_len"])
        ref_bins = bins_for_interval(rel_ref_start, rel_ref_end, output_bins)
        edit_bins = bins_for_interval(rel_edit_start, rel_edit_end, output_bins)
        union_bins = sorted(set(ref_bins) | set(edit_bins))
        if not union_bins:
            raise RuntimeError(f"No output bins overlap {ccre_id}")

        all_delta = edit_pred[union_bins, :] - ref_pred[union_bins, :]
        ccre_rows.append(
            {
                "cCRE_id": ccre_id,
                "n_output_bins": len(union_bins),
                "all_mouse_tracks_abs_delta_sum": float(np.abs(all_delta).sum()),
                "all_mouse_tracks_delta_sum": float(all_delta.sum()),
                "all_mouse_tracks_max_abs_delta": float(np.abs(all_delta).max()),
            }
        )

        for _target_row_idx, target in targets.iterrows():
            channel = int(target["mouse_output_channel"])
            if channel >= output_tracks:
                continue
            ref_vals = ref_pred[union_bins, channel]
            edit_vals = edit_pred[union_bins, channel]
            delta_vals = edit_vals - ref_vals
            long_rows.append(
                {
                    "cCRE_id": ccre_id,
                    "mouse_output_channel": channel,
                    "target_index_metadata": target.get("index", ""),
                    "identifier": target.get("identifier", ""),
                    "description": target.get("description", ""),
                    "n_output_bins": len(union_bins),
                    "mouse_reference_sum": float(ref_vals.sum()),
                    "spic_replacement_sum": float(edit_vals.sum()),
                    "delta_sum": float(delta_vals.sum()),
                    "abs_delta_sum": float(np.abs(delta_vals).sum()),
                    "max_abs_delta": float(np.abs(delta_vals).max()),
                    "mean_log2_fc": float(
                        np.mean(np.log2((edit_vals + 1e-6) / (ref_vals + 1e-6)))
                    ),
                }
            )
            if args.write_bin_summaries:
                for output_bin, ref_val, edit_val in zip(union_bins, ref_vals, edit_vals):
                    bin_rows.append(
                        {
                            "cCRE_id": ccre_id,
                            "mouse_output_channel": channel,
                            "output_bin": output_bin,
                            "mouse_reference": float(ref_val),
                            "spic_replacement": float(edit_val),
                            "delta": float(edit_val - ref_val),
                        }
                    )

        pd.DataFrame(long_rows).to_csv(target_partial, index=False)
        pd.DataFrame(ccre_rows).to_csv(ccre_partial, index=False)

    pd.DataFrame(long_rows).to_csv(args.outdir / "enformer_local_delta_by_target.csv", index=False)
    pd.DataFrame(ccre_rows).to_csv(args.outdir / "enformer_local_delta_by_ccre.csv", index=False)
    if args.write_bin_summaries:
        pd.DataFrame(bin_rows).to_csv(args.outdir / "enformer_local_delta_by_bin.csv", index=False)
    print(f"Wrote outputs to {args.outdir}", flush=True)


if __name__ == "__main__":
    main()
