# 01 — Prepare orthologs and alignments

## Purpose

Prepare the locus-specific multiple-sequence alignments used downstream for the PAML coding-sequence analyses and phyloP cCRE analyses.


## Source multiple alignments

The comparative-genomic source is release-112 of the Ensembl 21-way Murinae
Enredo-Pecan–Ortheus (EPO) multiple alignment.

Ensembl source pages:

- MLSS information:
  https://useast.ensembl.org/info/genome/compara/mlss.html?mlss=9598
- EPO README:
  https://ftp.ensembl.org/pub/release-112/emf/ensembl-compara/multiple_alignments/21_murinae.epo/README.21_murinae.epo


## Focal taxa

The analysis retains these eight taxa:

| Repository name | Ensembl assembly / strain |
|---|---|
| `rattus_norvegicus` | Rat, mRatBN7.2 |
| `mus_pahari` | PAHARI_EIJ_v1.1 |
| `mus_caroli` | CAROLI_EIJ_v1.1 |
| `mus_spretus` | SPRET_EiJ_v1 |
| `mus_spicilegus` | MUSP714 |
| `mus_musculus_wsbeij` | WSB/EiJ, WSB_EiJ_v1 |
| `mus_musculus_pwkphj` | PWK/PhJ, PWK_PhJ_v1 |
| `mus_musculus_casteij` | CAST/EiJ, CAST_EiJ_v1 |

`species_map.tsv` maps Ensembl MAF source names to these repository names.

## cCRE source

Mouse candidate cis-regulatory elements (cCREs) and their genomic
coordinates were obtained from the 2025 release of the ENCODE
Registry/SCREEN:

https://screen.wenglab.org/downloads

The analyses used the Mouse "cCREs by Cell and Tissue Type"
resources corresponding to the tissues and developmental time
points described in the manuscript.

cCRE coordinates were used in the mm10/GRCm38 reference coordinate
system provided by SCREEN and represented as BED-like intervals:

    chrom    start    end    locus_id

Coordinates are 0-based and half-open.


## Reference coordinate systems

Coding-sequence and regulatory analyses were processed in their respective reference coordinate systems. Coding-sequence analyses used GRCm39-reference alignments. cCRE analyses retained the mm10/GRCm38 coordinates of the SCREEN regulatory annotations and used corresponding mouse-reference comparative alignments. No coordinate liftOver between mm10 and GRCm39 was applied.


## Workflow

### 1. Convert Ensembl EMF to MAF

Ensembl Compara EPO alignments distributed in EMF format were
converted to MAF using the Ensembl Compara `emf2maf.pl` utility.

A representative invocation is:

    perl <ensembl-compara>/scripts/dumps/emf2maf.pl \
        < chromosome_alignment.emf \
        > chromosome_alignment.maf

### 2. Prepare BED-like locus files

Target loci for each analysis were represented as tab-delimited genomic intervals

    chrom    start    end    locus_id

Examples:

    chr1    100000    100450    ENSMUSG00000000001
    chr2    500000    500350    EH38E0000001

Coding and cCRE intervals were maintained separately because they use different mouse reference-coordinate systems.

### 3. Extract overlapping MAF blocks

`extract_maf_intervals.py` extracts the columns of every MAF block overlapping
each target interval. It writes one block-separated FASTA file per locus.

Example:

    python extract_maf_intervals.py \
        --maf chromosome_alignment.maf \
        --intervals chromosome_targets.bed \
        --species-map species_map.tsv \
        --reference-pattern '<mouse-reference identifier>' \
        --outdir extracted_blocks/

The reference pattern must uniquely identify the mouse-reference `s` row in
each MAF block. Inspect a representative MAF file before choosing it.

### 4. Stitch blocks belonging to each locus

Run `stitch_blocks.py` on the directory produced above:

    python stitch_blocks.py extracted_blocks/chr1 stitched_alignments/chr1

`stitch_blocks.py` concatenates the extracted blocks and pads focal taxa that
are absent from an individual block with gaps.

### 5. Downstream QC

Coding-sequence QC is performed separately with
`qc_gene_full_sequences.py` and the additional alignment-QC procedures
described in Methods section 2.

## Files in this directory

- `README.md` — provenance and reproducible workflow.
- `extract_maf_intervals.py` — reconstructed interval extractor.
- `species_map.tsv` — mapping from Ensembl MAF source names to the eight focal taxa.
- `stitch_blocks.py` — existing block-stitching script.
- `qc_gene_full_sequences.py` — existing coding-sequence QC utility; conceptually belongs mainly to Methods section 2.
