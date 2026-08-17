#!/usr/bin/env python3
import sys, os, glob, textwrap

if len(sys.argv) != 3:
    print("usage: stitch_blocks.py INDIR OUTDIR", file=sys.stderr); sys.exit(1)

INDIR, OUTDIR = sys.argv[1], sys.argv[2]
os.makedirs(OUTDIR, exist_ok=True)

# Desired taxa + order
order = [
    "rattus_norvegicus",
    "mus_pahari",
    "mus_caroli",
    "mus_spretus",
    "mus_spicilegus",
    "mus_musculus_wsbeij",
    "mus_musculus_pwkphj",
    "mus_musculus_casteij",
]

def write_fasta(path, seqs):
    with open(path, "w") as out:
        for name in order:
            out.write(f">{name}\n")
            for line in textwrap.wrap(seqs[name], width=80):
                out.write(line + "\n")

for fa in glob.glob(os.path.join(INDIR, "*.fa")):
    # Parse into blocks (blank line separates blocks)
    blocks = []
    cur = {}
    cur_name = None
    with open(fa) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                if cur:
                    blocks.append(cur); cur = {}; cur_name = None
                continue
            if line.startswith(">"):
                cur_name = line[1:].split()[0]
                cur.setdefault(cur_name, "")
            else:
                if cur_name is None:
                    raise RuntimeError(f"sequence before header in {fa}")
                cur[cur_name] += line
    if cur:
        blocks.append(cur)

    # Determine each block length; sanity check uniformity within block
    blens = []
    for i, b in enumerate(blocks):
        lens = {len(s) for s in b.values()}
        if len(lens) != 1:
            print(f"[warn] {os.path.basename(fa)} block {i+1} has varying lengths: {sorted(lens)}",
                  file=sys.stderr)
        blens.append(max(lens) if lens else 0)

    # Stitch across blocks, padding missing taxa with '-' of block length
    stitched = {name: "" for name in order}
    for b, L in zip(blocks, blens):
        for name in order:
            stitched[name] += b.get(name, "-" * L)

    # Write output
    out_path = os.path.join(OUTDIR, os.path.basename(fa))
    write_fasta(out_path, stitched)

    # Quick checks
    lens = {len(s) for s in stitched.values()}
    if len(lens) != 1:
        print(f"[warn] {os.path.basename(fa)} stitched sequences differ in length: {sorted(lens)}",
              file=sys.stderr)
