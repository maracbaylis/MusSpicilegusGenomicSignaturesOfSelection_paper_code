#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ASSET_DIR="${ROOT}/enformer_spicilegus_phyloP/assets"
mkdir -p "${ASSET_DIR}"

MM10_URL="https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz"
MM10_GZ="${ASSET_DIR}/mm10.fa.gz"
MM10_FA="${ASSET_DIR}/mm10.fa"

if [[ ! -s "${MM10_GZ}" && ! -s "${MM10_FA}" ]]; then
  curl -L "${MM10_URL}" -o "${MM10_GZ}"
fi

if [[ ! -s "${MM10_FA}" ]]; then
  python3 - <<PY
import gzip
from pathlib import Path
src = Path("${MM10_GZ}")
dst = Path("${MM10_FA}")
with gzip.open(src, "rb") as inp, dst.open("wb") as out:
    while True:
        chunk = inp.read(1024 * 1024)
        if not chunk:
            break
        out.write(chunk)
PY
fi

TARGETS_URL="https://raw.githubusercontent.com/calico/basenji/master/manuscripts/cross2020/targets_mouse.txt"
TARGETS_OUT="${ASSET_DIR}/enformer_targets_mouse.txt"

if [[ ! -s "${TARGETS_OUT}" ]]; then
  curl -L "${TARGETS_URL}" -o "${TARGETS_OUT}" || true
fi

echo "Assets directory: ${ASSET_DIR}"
echo "Mouse FASTA: ${MM10_FA}"
echo "Mouse targets: ${TARGETS_OUT}"
