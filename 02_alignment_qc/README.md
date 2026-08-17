# 02 — Alignment quality control

## Purpose

Apply alignment-quality filters before lineage-specific PAML and phyloP
analyses. These checks are intended to reduce false signals caused by assembly,
orthology, or local alignment artifacts.

## QC procedure

Each locus is required to have:

- an identifiable foreground sequence;
- equal alignment lengths across retained taxa; and
- callable bases in the foreground sequence.

Columns containing a gap or ambiguous base in the foreground are excluded from
foreground-private substitution summaries.

A foreground-private substitution is counted only when:

1. the foreground base is callable (`A`, `C`, `G`, or `T`);
2. at least five nonforeground taxa have callable bases at that position;
3. the nonforeground taxa support a single majority base at frequency >= 0.75;
   and
4. the foreground base differs from that majority state.

For each alignment, `alignment_qc.py` reports:

- foreground callable coverage;
- foreground-private substitution count;
- foreground-private substitution rate;
- maximum number of foreground-private substitutions in any 10-bp window;
- longest serial cluster of private substitutions separated by no more than
  2 bp; and
- number of foreground-private gaps.

Alignments are excluded when they show any of the manuscript artifact criteria:

```text
private substitution rate > 0.08
private substitution count >= 25
private substitutions in any 10-bp window >= 7
foreground-private gaps >= 3
```

## Running QC

Example for the *M. spicilegus* foreground:

```bash
python alignment_qc.py \
    stitched_alignments/ \
    --foreground mus_spicilegus \
    --out outputs/spicilegus_alignment_qc.csv
```

Run separately for each focal lineage used in a lineage-specific analysis.

## Output

The output CSV contains one row per alignment, including the QC metrics,
`qc_status` (`PASS` or `FAIL`), and a semicolon-delimited `qc_reasons` field.

Only alignments passing the relevant QC filters should be forwarded to the
downstream PAML or phyloP analyses.

## Repository contents

```text
02_alignment_qc/
├── README.md
└── alignment_qc.py
```

## Requirements

- Python >= 3.8
