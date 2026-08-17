#!/usr/bin/env python3
"""Build normalized external gene sets for M. spicilegus candidate overlap tests."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"
OUT = ROOT / "processed"


def norm_gene(value):
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none"}:
        return None
    text = re.sub(r"\s+", "", text)
    return text


def add_rows(rows, dataset, subset, source, genes, evidence="", notes="", species=""):
    seen = set()
    for gene in genes:
        gene = norm_gene(gene)
        if not gene:
            continue
        key = gene.upper()
        if key in seen:
            continue
        seen.add(key)
        rows.append(
            {
                "dataset": dataset,
                "subset": subset,
                "gene_symbol": gene,
                "gene_symbol_upper": key,
                "source": source,
                "evidence": evidence,
                "species_or_mapping": species,
                "notes": notes,
            }
        )


def read_johnson(rows):
    source = "Johnson et al. 2023 Nat Commun, doi:10.1038/s41467-023-40331-9"

    ieg = pd.read_excel(RAW / "Johnson2023_MOESM9.xlsx", sheet_name="IEGs")
    add_rows(
        rows,
        "Johnson2023_cichlid_bower_snRNAseq",
        "IEG_like_genes",
        source + "; Supplementary Data 6, sheet IEGs",
        ieg["human"],
        "Genes exhibiting IEG-like expression",
        species="human ortholog column from mzebra genes",
    )

    deg_primary = pd.read_excel(RAW / "Johnson2023_MOESM10.xlsx", sheet_name="Primary Cluster Results")
    deg_secondary = pd.read_excel(RAW / "Johnson2023_MOESM10.xlsx", sheet_name="Secondary Cluster Results")
    for label, df in [("primary_clusters", deg_primary), ("secondary_clusters", deg_secondary)]:
        for category in sorted(df["category"].dropna().unique()):
            sub = df[df["category"] == category]
            add_rows(
                rows,
                "Johnson2023_cichlid_bower_snRNAseq",
                f"{category}_{label}",
                source + f"; Supplementary Data 7, {label}",
                sub["human"],
                "Linear mixed-effect behavior-associated expression result",
                species="human ortholog column from mzebra genes",
            )

    rg_all = pd.read_excel(RAW / "Johnson2023_MOESM14.xlsx", sheet_name="cond_deg_rg_all")
    rg_sub = pd.read_excel(RAW / "Johnson2023_MOESM14.xlsx", sheet_name="cond_deg_rg_subclusters")
    for label, df in [("radial_glia_all", rg_all), ("radial_glia_subclusters", rg_sub)]:
        for category in sorted(df["category"].dropna().unique()):
            sub = df[df["category"] == category]
            add_rows(
                rows,
                "Johnson2023_cichlid_bower_snRNAseq",
                f"{category}_{label}",
                source + f"; Supplementary Data 11, {label}",
                sub["hgnc"],
                "Building-associated expression in radial glia",
                species="HGNC ortholog column from mzebra genes",
            )

    proneuro = pd.read_csv(RAW / "Johnson2023_MOESM13.csv")
    add_rows(
        rows,
        "Johnson2023_cichlid_bower_snRNAseq",
        "positive_regulation_of_neurogenesis",
        source + "; Supplementary Data 10",
        proneuro["human"],
        "Genes with GO positive regulation of neurogenesis in zebrafish and mice",
        species="human ortholog column from mzebra genes",
    )

    rg_states = pd.read_csv(RAW / "Johnson2023_MOESM15.csv")
    rg_states["rg_state_clean"] = rg_states["rg_state"].astype(str).str.strip()
    for state in sorted(rg_states["rg_state_clean"].dropna().unique()):
        sub = rg_states[rg_states["rg_state_clean"] == state]
        state_label = re.sub(r"[^A-Za-z0-9]+", "_", state).strip("_")
        add_rows(
            rows,
            "Johnson2023_cichlid_bower_snRNAseq",
            f"radial_glia_state_{state_label}",
            source + "; Supplementary Data 12",
            sub["human"],
            "Markers of radial-glial functional states",
            species="human ortholog column from mzebra genes",
        )

    cdg = pd.read_excel(RAW / "Johnson2023_MOESM18.xlsx", sheet_name="CDG and CDG Module")
    add_rows(
        rows,
        "Johnson2023_cichlid_bower_comparative_genomics",
        "castle_divergent_genes_CDG",
        source + "; Supplementary Data 15, sheet CDG and CDG Module",
        cdg["mouse"],
        "Castle-divergent gene list, mouse ortholog column",
        species="mouse ortholog column from mzebra genes",
    )
    add_rows(
        rows,
        "Johnson2023_cichlid_bower_comparative_genomics",
        "castle_divergent_genes_CDG_module",
        source + "; Supplementary Data 15, sheet CDG and CDG Module",
        cdg.loc[cdg["in_CDG_module"] == True, "mouse"],
        "Castle-divergent gene module subset, mouse ortholog column",
        species="mouse ortholog column from mzebra genes",
    )


def read_thompson(rows):
    source = "Thompson et al. 2012 PLoS ONE, doi:10.1371/journal.pone.0035119"
    files = [
        ("Thompson2012_s001.xls", "HVC_breeding_induction"),
        ("Thompson2012_s002.xls", "RA_breeding_induction"),
        ("Thompson2012_s003.xls", "HVC_regression_to_nonbreeding"),
        ("Thompson2012_s004.xls", "RA_regression_to_nonbreeding"),
    ]
    for filename, subset in files:
        df = pd.read_excel(RAW / filename, header=1)
        if "HGNC symbol" not in df.columns:
            continue
        add_rows(
            rows,
            "Thompson2012_seasonal_songbird_plasticity",
            subset,
            f"{source}; {filename}",
            df["HGNC symbol"],
            "Spots varying >1.5-fold with raw p < 0.01 in song-control nucleus",
            species="HGNC-style symbols from zebra finch microarray annotation",
        )


def read_olson(rows):
    source = "Olson et al. 2015 Dev Neurobiol, doi:10.1002/dneu.22286"
    groups = {
        "HVC_higher_in_adults": ["ADAM23", "COL12A1", "COL21A1", "CXCR7", "MPZL1", "PVALB", "SLC8A1"],
        "HVC_higher_in_juveniles": ["CACNG3", "CNTN4", "FLRT2", "KCNT2", "LGMN", "LRRTM2", "NRP1", "SAP30L"],
        "shelf_higher_in_juveniles": ["MAP1B", "NDRG4", "RCAN2", "SEMA6A", "SV2B"],
        "marker_across_ages": ["CADPS2", "CDH6", "KCNA1", "NETO1", "NROB1", "PPM1E", "SCUBE1", "SEMA3E", "SYT4"],
        "area_X_marker_examined": ["TAC1"],
    }
    for subset, genes in groups.items():
        add_rows(
            rows,
            "Olson2015_zebra_finch_song_learning",
            subset,
            source + "; Table 1 and text",
            genes,
            "Known robust molecular markers assayed across song-learning development",
            species="zebra finch marker symbols, matched by orthologous symbol where possible",
        )


def write_outputs(rows):
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows).sort_values(["dataset", "subset", "gene_symbol_upper", "gene_symbol"])
    df.to_csv(OUT / "external_gene_sets_long.csv", index=False)

    summary = (
        df.groupby(["dataset", "subset"], as_index=False)
        .agg(n_genes=("gene_symbol_upper", "nunique"), source=("source", "first"))
        .sort_values(["dataset", "subset"])
    )
    summary.to_csv(OUT / "external_gene_sets_summary.csv", index=False)

    # Convenient GMT-like file for enrichment scripts.
    with (OUT / "external_gene_sets.gmt").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for (dataset, subset), sub in df.groupby(["dataset", "subset"], sort=True):
            genes = sorted(sub["gene_symbol_upper"].unique())
            writer.writerow([f"{dataset}::{subset}", sub["source"].iloc[0], *genes])

    manifest = OUT / "external_gene_sets_manifest.md"
    manifest.write_text(
        "# External Gene Sets\n\n"
        "Generated by `external_gene_sets/build_external_gene_sets.py`.\n\n"
        "Included machine-readable sets:\n\n"
        + "\n".join(
            f"- `{r.dataset}::{r.subset}`: {r.n_genes} genes"
            for r in summary.itertuples(index=False)
        )
        + "\n\nNot yet included as full gene lists:\n\n"
        "- York et al. 2018 PNAS Dataset S1: PNAS direct download returned a Cloudflare challenge page in this environment.\n"
        "- Won et al. 2019 HAR/HGE/HLE target genes: PMC/Nature identify Supplementary Data workbooks, but binary downloads required a proof-of-work/challenge page here; the article reports 1,648 HAR putative target genes.\n"
        "- Weber et al. 2013 Peromyscus burrowing: the paper is QTL/conceptual for this use case, not a clean orthologous gene-list overlap set.\n"
        "- Warren et al. 2010 zebra finch genome: useful background for song-regulated rapid evolution, but not extracted here as a clean one-gene-per-row external set.\n",
        encoding="utf-8",
    )


def main():
    rows = []
    read_johnson(rows)
    read_thompson(rows)
    read_olson(rows)
    write_outputs(rows)


if __name__ == "__main__":
    main()
