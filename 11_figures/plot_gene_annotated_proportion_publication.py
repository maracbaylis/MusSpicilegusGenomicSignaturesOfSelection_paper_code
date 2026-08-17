from __future__ import annotations

from pathlib import Path
import os
from zipfile import ZipFile
from xml.etree import ElementTree as ET

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.ticker import MultipleLocator


INPUT_XLSX = Path(os.environ.get("MSPIC_INPUT_DIR", "data/inputs")).resolve() / 'enhancer_presence_matrix.xlsx'
OUTPUT_PNG = Path(os.environ.get("MSPIC_PROJECT_ROOT", ".")).resolve() / 'gene_annotated_proportion_publication.png'
OUTPUT_PDF = Path(os.environ.get("MSPIC_PROJECT_ROOT", ".")).resolve() / 'gene_annotated_proportion_publication.pdf'

TIMEPOINT_ORDER = ["4d", "10d", "14d", "25d", "36d", "2m", "8-10m", "18-20m"]
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
    """Read a simple .xlsx worksheet using only the Python standard library."""
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

    if not rows:
        raise ValueError(f"No rows found in {path}")

    header_map = rows[0]
    records = []
    for row in rows[1:]:
        record = {column_name: row.get(column_letter, "") for column_letter, column_name in header_map.items()}
        records.append(record)

    dataframe = pd.DataFrame(records)
    for timepoint in TIMEPOINT_ORDER + ["Total_Appearances"]:
        if timepoint in dataframe.columns:
            dataframe[timepoint] = pd.to_numeric(dataframe[timepoint], errors="coerce").fillna(0)
    return dataframe


def build_summary_table(dataframe: pd.DataFrame) -> pd.DataFrame:
    summary_rows = []
    for timepoint in TIMEPOINT_ORDER:
        present_mask = dataframe[timepoint] > 0
        present = dataframe.loc[present_mask]
        total_present = len(present)
        annotated = present["Gene Name"].fillna("").astype(str).str.strip().ne("").sum()
        percent = (annotated / total_present * 100) if total_present else 0
        summary_rows.append(
            {
                "timepoint": timepoint,
                "total_present": total_present,
                "annotated": annotated,
                "percent_annotated": percent,
            }
        )
    return pd.DataFrame(summary_rows)


def configure_matplotlib() -> None:
    sns.set_theme(style="whitegrid")
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "axes.titlesize": 16,
            "axes.titleweight": "bold",
            "axes.labelsize": 12,
            "axes.labelweight": "bold",
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.bbox": "tight",
            "savefig.facecolor": "white",
        }
    )


def make_plot(dataframe: pd.DataFrame, summary_df: pd.DataFrame, output_png: Path, output_pdf: Path) -> None:
    configure_matplotlib()

    bar_color = "#295C77"
    edge_color = "#173645"
    highlight_color = "#C9792B"

    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    bars = ax.bar(
        summary_df["timepoint"],
        summary_df["percent_annotated"],
        color=bar_color,
        edgecolor=edge_color,
        linewidth=1.0,
        width=0.72,
    )

    ax.set_ylim(0, 106)
    ax.set_ylabel("Gene-annotated enhancers (%)")
    ax.set_xlabel("Developmental timepoint")
    ax.set_title("All detected enhancers are gene-annotated across timepoints", pad=14)

    ax.yaxis.set_major_locator(MultipleLocator(20))
    ax.grid(axis="y", color="#D7DEE3", linewidth=0.8)
    ax.grid(axis="x", visible=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#55636E")
    ax.spines["bottom"].set_color("#55636E")

    for bar, row in zip(bars, summary_df.itertuples(index=False)):
        x_center = bar.get_x() + bar.get_width() / 2
        ax.text(
            x_center,
            row.percent_annotated + 1.2,
            f"{row.percent_annotated:.0f}%",
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
            color=highlight_color,
        )
        ax.text(
            x_center,
            max(row.percent_annotated - 8.5, 6),
            f"{row.annotated}/{row.total_present}",
            ha="center",
            va="top",
            fontsize=9,
            color="white",
            fontweight="bold",
        )

    total_unique_enhancers = int(dataframe["Enhancer_ID"].nunique())
    annotation = (
        f"Workbook contains {total_unique_enhancers} unique enhancers"
        "\nPercent calculated among enhancers present at each stage"
    )
    fig.text(0.015, -0.02, annotation, ha="left", va="top", fontsize=9.5, color="#4D5B65")

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=600)
    fig.savefig(output_pdf)
    plt.close(fig)


def main() -> None:
    dataframe = read_xlsx_without_openpyxl(INPUT_XLSX)
    summary_df = build_summary_table(dataframe)
    make_plot(dataframe, summary_df, OUTPUT_PNG, OUTPUT_PDF)
    print(summary_df.to_string(index=False))
    print(f"Saved {OUTPUT_PNG}")
    print(f"Saved {OUTPUT_PDF}")


if __name__ == "__main__":
    main()
