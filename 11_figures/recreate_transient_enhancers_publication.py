from __future__ import annotations

from pathlib import Path
import os
from zipfile import ZipFile
from xml.etree import ElementTree as ET
import colorsys

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter, MultipleLocator


INPUT_XLSX = Path(os.environ.get("MSPIC_INPUT_DIR", "data/inputs")).resolve() / 'enhancer_presence_matrix.xlsx'
OUTPUT_PNG = Path(os.environ.get("MSPIC_PROJECT_ROOT", ".")).resolve() / 'transient_enhancers_publication_ready.png'
OUTPUT_PDF = Path(os.environ.get("MSPIC_PROJECT_ROOT", ".")).resolve() / 'transient_enhancers_publication_ready.pdf'
ALT_OUTPUT_PNG = Path(os.environ.get("MSPIC_PROJECT_ROOT", ".")).resolve() / 'transient_enhancers_publication_ready_green.png'
ALT_OUTPUT_PDF = Path(os.environ.get("MSPIC_PROJECT_ROOT", ".")).resolve() / 'transient_enhancers_publication_ready_green.pdf'
PINK_OUTPUT_PNG = Path(os.environ.get("MSPIC_PROJECT_ROOT", ".")).resolve() / 'transient_enhancers_publication_ready_pink.png'
PINK_OUTPUT_PDF = Path(os.environ.get("MSPIC_PROJECT_ROOT", ".")).resolve() / 'transient_enhancers_publication_ready_pink.pdf'

TIMEPOINTS = ["4d", "10d", "14d", "25d", "36d", "2m", "8-10m", "18-20m"]
TIMEPOINT_LABELS = [
    "4 days",
    "10 days",
    "14 days",
    "25 days",
    "36 days",
    "2 months",
    "8-10 months",
    "18-20 months",
]

# These discovery-rate totals come from the original plotted figure.
TIMEPOINT_TOTALS = {
    "4d": 0.0084,
    "10d": 0.0118,
    "14d": 0.0160,
    "25d": 0.0225,
    "36d": 0.0158,
    "2m": 0.0294,
    "8-10m": 0.0133,
    "18-20m": 0.0111,
}

XML_NS = {"main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def _column_letters(cell_ref: str) -> str:
    letters = []
    for char in cell_ref:
        if char.isalpha():
            letters.append(char)
        else:
            break
    return "".join(letters)


def read_xlsx_without_openpyxl(path: Path) -> pd.DataFrame:
    with ZipFile(path) as workbook:
        shared_strings = []
        if "xl/sharedStrings.xml" in workbook.namelist():
            string_root = ET.fromstring(workbook.read("xl/sharedStrings.xml"))
            for string_item in string_root:
                text = "".join(
                    node.text or ""
                    for node in string_item.iter("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t")
                )
                shared_strings.append(text)

        sheet_root = ET.fromstring(workbook.read("xl/worksheets/sheet1.xml"))
        row_nodes = sheet_root.find("main:sheetData", XML_NS)
        if row_nodes is None:
            raise ValueError(f"No worksheet data found in {path}")

        rows = []
        for row in row_nodes:
            values = {}
            for cell in row:
                ref = cell.attrib["r"]
                column = _column_letters(ref)
                value_node = cell.find("main:v", XML_NS)
                value = "" if value_node is None else value_node.text or ""
                if cell.attrib.get("t") == "s" and value:
                    value = shared_strings[int(value)]
                values[column] = value
            rows.append(values)

    header_map = rows[0]
    records = []
    for row in rows[1:]:
        record = {column_name: row.get(column_letter, "") for column_letter, column_name in header_map.items()}
        records.append(record)

    dataframe = pd.DataFrame(records)
    for timepoint in TIMEPOINTS + ["Total_Appearances"]:
        dataframe[timepoint] = pd.to_numeric(dataframe[timepoint], errors="coerce").fillna(0)
    return dataframe


def build_transient_dataframe(dataframe: pd.DataFrame) -> pd.DataFrame:
    transient_df = dataframe.loc[dataframe["Total_Appearances"] < 4].copy()
    transient_df["display_label"] = transient_df["Gene Name"].fillna("").str.strip()
    blank_labels = transient_df["display_label"].eq("")
    transient_df.loc[blank_labels, "display_label"] = transient_df.loc[blank_labels, "Enhancer_ID"]
    return transient_df


def build_color_lookup(
    transient_df: pd.DataFrame,
    hue_start: float = 0.55,
    hue_end: float = 0.64,
) -> dict[str, tuple[float, float, float]]:
    enhancer_ids = sorted(transient_df["Enhancer_ID"].tolist())
    n_colors = len(enhancer_ids)
    hues = np.linspace(hue_start, hue_end, n_colors)
    palette = []
    for idx, hue in enumerate(hues):
        saturation = 0.42 if idx % 2 == 0 else 0.58
        lightness = 0.50 if idx % 3 == 0 else 0.58
        palette.append(colorsys.hls_to_rgb(float(hue), float(lightness), float(saturation)))
    return dict(zip(enhancer_ids, palette))


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "axes.titlesize": 16,
            "axes.titleweight": "normal",
            "axes.labelsize": 13,
            "axes.labelweight": "normal",
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "axes.linewidth": 1.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "savefig.bbox": "tight",
        }
    )


def build_figure(
    transient_df: pd.DataFrame,
    output_png: Path,
    output_pdf: Path,
    hue_start: float = 0.55,
    hue_end: float = 0.64,
) -> None:
    configure_style()
    color_lookup = build_color_lookup(transient_df, hue_start=hue_start, hue_end=hue_end)

    x = np.arange(len(TIMEPOINTS))
    totals = np.array([TIMEPOINT_TOTALS[tp] for tp in TIMEPOINTS])

    fig, ax = plt.subplots(figsize=(8.8, 5.6))

    for idx, timepoint in enumerate(TIMEPOINTS):
        tp_df = transient_df.loc[transient_df[timepoint] > 0].copy()
        tp_df = tp_df.sort_values(["Enhancer_ID", "Gene Name"], kind="stable")
        total = TIMEPOINT_TOTALS[timepoint]
        segment_height = total / len(tp_df)

        bottom = 0.0
        for row in tp_df.itertuples(index=False):
            ax.bar(
                idx,
                segment_height,
                bottom=bottom,
                width=0.62,
                color=color_lookup[row.Enhancer_ID],
                edgecolor=(0, 0, 0, 0.16),
                linewidth=0.32,
            )
            bottom += segment_height

    ax.set_title(
        "Active CREs with significant phyloP acceleration\nin left cerebral cortex",
        pad=14,
    )
    ax.set_ylabel("Discovery rate of accelerated CREs (%)")
    ax.set_xlabel("Developmental timepoint")
    ax.set_xticks(x)
    ax.set_xticklabels(TIMEPOINT_LABELS, rotation=45, ha="right")

    ax.set_ylim(0, 0.031)
    ax.yaxis.set_major_locator(MultipleLocator(0.005))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=False)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#59656F")
    ax.spines["bottom"].set_color("#59656F")
    ax.tick_params(colors="#2F2F2F")

    note = (
        "Stacks are Excel-derived transient CRE identities with consistent colors across timepoints.\n"
        "Bar heights reproduce the original discovery-rate totals from the source figure."
    )
    fig.text(0.015, -0.015, note, ha="left", va="top", fontsize=9, color="#51606C")

    fig.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=600)
    fig.savefig(output_pdf)
    plt.close(fig)


def main() -> None:
    dataframe = read_xlsx_without_openpyxl(INPUT_XLSX)
    transient_df = build_transient_dataframe(dataframe)
    build_figure(transient_df, OUTPUT_PNG, OUTPUT_PDF, hue_start=0.55, hue_end=0.64)
    build_figure(transient_df, ALT_OUTPUT_PNG, ALT_OUTPUT_PDF, hue_start=0.33, hue_end=0.42)
    build_figure(transient_df, PINK_OUTPUT_PNG, PINK_OUTPUT_PDF, hue_start=0.88, hue_end=0.97)
    print(transient_df[["Enhancer_ID", "Gene Name", "Total_Appearances"] + TIMEPOINTS].to_string(index=False))
    print(f"Saved {OUTPUT_PNG}")
    print(f"Saved {OUTPUT_PDF}")
    print(f"Saved {ALT_OUTPUT_PNG}")
    print(f"Saved {ALT_OUTPUT_PDF}")
    print(f"Saved {PINK_OUTPUT_PNG}")
    print(f"Saved {PINK_OUTPUT_PDF}")


if __name__ == "__main__":
    main()
