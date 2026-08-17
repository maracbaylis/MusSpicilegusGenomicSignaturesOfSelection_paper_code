# 08 — Enformer regulatory-effect prediction

## Purpose

Estimate the predicted regulatory effects of *M. spicilegus*-specific sequence
changes in accelerated brain-active cCREs using Enformer and compare those
effects with matched non-accelerated cCREs.

The primary workflow uses paired sequence inputs representing:

1. the inferred ancestral/non-*M. spicilegus* allele in mouse genomic context;
   and
2. the corresponding *M. spicilegus*-derived allele.

Predicted differences are summarized across mouse Enformer output tracks, with
the primary analysis focused on neural tissue/cell-context tracks.

## Repository contents

```text
08_enformer_predictions/
├── README.md
├── download_assets.sh
├── build_spicilegus_enformer_input_manifest.py
├── build_background_msa_manifest.py
├── make_enformer_windows.py
├── select_brain_mouse_targets.py
├── run_enformer_predictions.py
└── compute_empirical_enformer_significance.py
```

## External assets

The workflow requires:

- the UCSC mm10 mouse reference genome; and
- Enformer/Basenji mouse target metadata (`targets_mouse.txt`).

These can be downloaded with:

```bash
bash download_assets.sh
```

The Enformer model itself is loaded from TensorFlow Hub:

```text
https://tfhub.dev/deepmind/enformer/1
```

## 1. Build accelerated-cCRE allele manifest

`build_spicilegus_enformer_input_manifest.py` takes the accelerated cCRE table
and corresponding cCRE multiple-sequence alignments and records:

- cCRE coordinates;
- inferred mouse/ancestral alignment consensus;
- *M. spicilegus* ortholog sequence;
- substitution and indel counts; and
- allele FASTA records used for window construction.

Use the exact accelerated cCRE set reported in the manuscript.

## 2. Build matched-background manifest

`build_background_msa_manifest.py` performs the same allele/variant extraction
for the matched non-accelerated cCRE set.

The background cCRE set should be the same matched set described in the
manuscript, selected to resemble the accelerated cCREs in sequence length, GC
content, and *M. spicilegus*-specific substitution burden.

## 3. Construct Enformer input windows

`make_enformer_windows.py` constructs full-length Enformer inputs centered on
the focal cCRE using the mm10 mouse reference genome.

The Enformer input length is:

```text
393216 bp
```

The reference sequence is taken from mm10. The alternative sequence is
constructed by replacing the focal cCRE interval with the inferred
*M. spicilegus*-derived sequence while leaving the surrounding genomic context
unchanged.

Run the same window-building procedure for accelerated and matched-background
manifests.

## 4. Select neural mouse output tracks

`select_brain_mouse_targets.py` filters the Enformer mouse target metadata to
the neural tissue/cell-context tracks used for the primary analysis.

```bash
python select_brain_mouse_targets.py \
    --targets assets/enformer_targets_mouse.txt \
    --out neural_mouse_targets.csv
```

The exact selected target table should be retained with the analysis outputs so
that the reported set of neural tracks is unambiguous.

## 5. Run Enformer predictions

`run_enformer_predictions.py` loads the DeepMind Enformer model from TensorFlow
Hub and predicts mouse genomic output tracks for paired reference and
*M. spicilegus*-edited windows.

```bash
python run_enformer_predictions.py \
    --window-manifest enformer_windows.csv \
    --targets neural_mouse_targets.csv \
    --outdir enformer_predictions/
```

For each cCRE and output track, the script summarizes the prediction difference
over Enformer output bins overlapping the focal cCRE.

## 6. Empirical significance using matched background loci

`compute_empirical_enformer_significance.py` compares predicted effect sizes
for accelerated cCREs with the empirical distribution from matched
non-accelerated cCREs.

```bash
python compute_empirical_enformer_significance.py \
    --foreground-ccre foreground/enformer_local_delta_by_ccre.csv \
    --background-ccre background/enformer_local_delta_by_ccre.csv \
    --foreground-target foreground/enformer_local_delta_by_target.csv \
    --background-target background/enformer_local_delta_by_target.csv \
    --outdir enformer_significance/
```

At the cCRE level, the primary effect metric is the summed absolute prediction
difference across mouse tracks. Empirical upper-tail P-values are calculated
from the matched-background distribution and adjusted across accelerated cCREs
using the Benjamini-Hochberg procedure.

Track-level empirical tests are calculated against matched-background values
for the corresponding mouse output channel.

## Requirements

- Python >= 3.10
- numpy
- pandas
- TensorFlow 2.21.0
- TensorFlow Hub 0.16.1
- certifi
- DeepMind Enformer (`https://tfhub.dev/deepmind/enformer/1`)
