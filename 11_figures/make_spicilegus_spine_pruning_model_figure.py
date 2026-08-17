from __future__ import annotations

import html
import math
from pathlib import Path
import os

import pandas as pd


ROOT = Path(os.environ.get("MSPIC_PROJECT_ROOT", ".")).resolve()
OUTDIR = ROOT / "figures"
JOINED = ROOT / "spicilegus_dendritic_spine_pruning_gene_trends.csv"
BM = ROOT / "bm_spicilegus_xlsx_results" / "brain_bm_residuals_all_genes.csv"
BM_Q01_STATS = (
    ROOT
    / "expression_evolution_5species"
    / "results_dom_altmap_rerun"
    / "bm_mus_complex_polytomy_tree"
    / "spicilegus_mus_complex_polytomy_q01_absz3_biggest_margin2_full_stats.csv"
)

OUT_SVG = OUTDIR / "spicilegus_spine_pruning_model_figure.svg"
OUT_DATA = OUTDIR / "spicilegus_spine_pruning_model_figure_source_data.csv"
OUT_PANEL_A = OUTDIR / "spicilegus_spine_pruning_model_panel_A_phenotype.svg"
OUT_PANEL_B = OUTDIR / "spicilegus_spine_pruning_model_panel_B_expression.svg"
OUT_PANEL_C = OUTDIR / "spicilegus_spine_pruning_model_panel_C_BM_residuals.svg"
OUT_PANEL_C_DATA = OUTDIR / "spicilegus_spine_pruning_model_panel_C_BM_residuals_source_data.csv"
OUT_PANEL_D = OUTDIR / "spicilegus_spine_pruning_model_panel_D_mechanism.svg"


W = 1800
H = 1220

COLORS = {
    "ink": "#202124",
    "muted": "#666D75",
    "grid": "#D7DDE3",
    "panel": "#F7F9FB",
    "blue": "#2F6FDB",
    "red": "#C9493D",
    "orange": "#D98B2B",
    "green": "#2F8F68",
    "purple": "#7357B8",
    "teal": "#178A9A",
}


def esc(text: object) -> str:
    return html.escape(str(text), quote=True)


def svg_text(
    x: float,
    y: float,
    text: object,
    size: int = 28,
    weight: int | str = 400,
    color: str = COLORS["ink"],
    anchor: str = "start",
    style: str = "",
) -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" fill="{color}" '
        f'text-anchor="{anchor}" {style}>{esc(text)}</text>'
    )


def rect(x: float, y: float, w: float, h: float, fill: str, stroke: str = "none", sw: float = 1, rx: float = 0) -> str:
    return (
        f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}" rx="{rx:.1f}"/>'
    )


def line(x1: float, y1: float, x2: float, y2: float, stroke: str = COLORS["ink"], sw: float = 2, dash: str | None = None) -> str:
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{stroke}" stroke-width="{sw}"{dash_attr}/>'


def arrow_marker_defs() -> str:
    return """
<defs>
  <marker id="arrowInk" markerWidth="12" markerHeight="12" refX="10" refY="6" orient="auto" markerUnits="strokeWidth">
    <path d="M2,2 L10,6 L2,10 Z" fill="#202124"/>
  </marker>
  <marker id="arrowRed" markerWidth="12" markerHeight="12" refX="10" refY="6" orient="auto" markerUnits="strokeWidth">
    <path d="M2,2 L10,6 L2,10 Z" fill="#C9493D"/>
  </marker>
  <marker id="arrowBlue" markerWidth="12" markerHeight="12" refX="10" refY="6" orient="auto" markerUnits="strokeWidth">
    <path d="M2,2 L10,6 L2,10 Z" fill="#2F6FDB"/>
  </marker>
</defs>
"""


def arrow(x1: float, y1: float, x2: float, y2: float, color: str = COLORS["ink"], sw: float = 3, marker: str = "arrowInk") -> str:
    return (
        f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
        f'stroke="{color}" stroke-width="{sw}" marker-end="url(#{marker})"/>'
    )


def wrapped_text(x: float, y: float, lines: list[str], size: int = 24, weight: int | str = 400, color: str = COLORS["ink"], gap: float = 1.25) -> list[str]:
    return [svg_text(x, y + i * size * gap, s, size=size, weight=weight, color=color) for i, s in enumerate(lines)]


def dendrite(x: float, y: float, spine_count: int, label: str) -> list[str]:
    parts: list[str] = []
    parts.append(svg_text(x, y - 52, label, size=24, weight=700, anchor="middle"))
    parts.append(line(x - 165, y, x + 165, y, stroke="#34495E", sw=8))
    if spine_count <= 1:
        positions = [0]
    else:
        positions = [(-140 + i * (280 / (spine_count - 1))) for i in range(spine_count)]
    for i, dx in enumerate(positions):
        stem = 26 + (i % 3) * 5
        up = i % 2 == 0
        y2 = y - stem if up else y + stem
        parts.append(line(x + dx, y, x + dx + (8 if i % 4 < 2 else -8), y2, stroke="#34495E", sw=4))
        parts.append(f'<circle cx="{x + dx + (8 if i % 4 < 2 else -8):.1f}" cy="{y2:.1f}" r="7" fill="#34495E"/>')
    return parts


def significance_label(q: float) -> str:
    if q < 0.001:
        return "q<0.001"
    if q < 0.01:
        return "q<0.01"
    if q < 0.05:
        return "q<0.05"
    return f"q={q:.2f}"


def significance_stars(q: float) -> str:
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return ""


def bar_panel(
    x: float,
    y: float,
    w: float,
    title: str,
    data: pd.DataFrame,
    value_col: str,
    q_col: str,
    xmin: float,
    xmax: float,
    axis_label: str,
    star_significance: bool = False,
    show_grid: bool = True,
    italicize_species_axis: bool = False,
    row_h: float = 34,
    title_size: int = 25,
    gene_col: str = "Genename",
) -> list[str]:
    parts: list[str] = []
    parts.append(svg_text(x, y, title, size=title_size, weight=700))
    plot_x = x + 132
    plot_y = y + 34
    plot_w = w - 178
    zero_x = plot_x + (0 - xmin) / (xmax - xmin) * plot_w
    if show_grid:
        parts.append(line(plot_x, plot_y - 8, plot_x + plot_w, plot_y - 8, stroke=COLORS["grid"], sw=1))
    parts.append(line(zero_x, plot_y - 18, zero_x, plot_y + row_h * len(data) + 2, stroke="#111111", sw=1.3))
    for tick in range(math.ceil(xmin), math.floor(xmax) + 1):
        tx = plot_x + (tick - xmin) / (xmax - xmin) * plot_w
        if show_grid:
            parts.append(line(tx, plot_y - 12, tx, plot_y + row_h * len(data), stroke=COLORS["grid"], sw=0.7, dash="3 5"))
        parts.append(svg_text(tx, plot_y + row_h * len(data) + 28, tick, size=17, color=COLORS["muted"], anchor="middle"))

    for i, row in enumerate(data.itertuples(index=False)):
        gene = getattr(row, gene_col)
        val = float(getattr(row, value_col))
        qval = float(getattr(row, q_col))
        yy = plot_y + i * row_h
        vx = plot_x + (val - xmin) / (xmax - xmin) * plot_w
        color = COLORS["blue"] if val > 0 else COLORS["red"]
        bx = min(vx, zero_x)
        bw = abs(vx - zero_x)
        parts.append(svg_text(x + 3, yy + 7, gene, size=20, weight=700 if qval < 0.05 else 400))
        parts.append(rect(bx, yy - 15, max(bw, 2), 20, color, rx=2))
        sig = significance_stars(qval) if star_significance else significance_label(qval)
        parts.append(svg_text(plot_x + plot_w + 16, yy + 3, sig, size=22 if star_significance else 17, weight=700 if star_significance else 400, color=COLORS["ink"] if star_significance else COLORS["muted"]))
    axis_x = plot_x + plot_w / 2
    axis_y = plot_y + row_h * len(data) + 58
    if italicize_species_axis:
        parts.append(
            f'<text x="{axis_x:.1f}" y="{axis_y:.1f}" font-family="Arial, Helvetica, sans-serif" '
            f'font-size="18" font-weight="400" fill="{COLORS["muted"]}" text-anchor="middle">'
            f'BM residual z-score in <tspan font-style="italic">M. spicilegus</tspan> brain</text>'
        )
    else:
        parts.append(svg_text(axis_x, axis_y, axis_label, size=18, color=COLORS["muted"], anchor="middle"))
    return parts


def svg_doc(width: int, height: int, body: list[str]) -> str:
    return "\n".join(
        [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
            arrow_marker_defs(),
            rect(0, 0, width, height, "white"),
            *body,
            "</svg>",
        ]
    )


def panel_a(x: float, y: float, w: float, h: float, label: bool = True) -> list[str]:
    parts: list[str] = []
    parts.append(rect(x, y, w, h, COLORS["panel"], stroke="#E2E8EF", sw=1.2, rx=8))
    if label:
        parts.append(svg_text(x + 26, y + 42, "A", size=32, weight=700))
    parts += wrapped_text(
        x + 78,
        y + 39,
        ["Observed phenotype", "higher dendritic spine density in both photoperiods"],
        size=25,
        weight=700,
        gap=1.2,
    )
    parts += dendrite(x + 530, y + 154, 12, "Long day")
    parts += dendrite(x + 980, y + 154, 12, "Short day")
    parts.append(arrow(x + 1210, y + 154, x + 1385, y + 154, color=COLORS["ink"], sw=3))
    parts += wrapped_text(
        x + 1410,
        y + 125,
        ["Species-level", "retention set point"],
        size=25,
        weight=700,
        color=COLORS["green"],
        gap=1.15,
    )
    return parts


def panel_b(x: float, y: float, w: float, h: float, expr: pd.DataFrame, label: bool = True) -> list[str]:
    parts: list[str] = []
    parts.append(rect(x, y, w, h, "white", stroke="#E2E8EF", sw=1.2, rx=8))
    if label:
        parts.append(svg_text(x + 24, y + 43, "B", size=32, weight=700))
    parts += bar_panel(
        x + 78,
        y + 46,
        w - 112,
        "Simple expression shifts",
        expr,
        "spic_minus_other_mean",
        "spic_minus_other_mean_q_bh",
        -6,
        5,
        "Spicilegus - mean(other species)",
    )
    parts.append(svg_text(x + 86, y + 405, "Red bars indicate lower expression/residuals; blue bars indicate higher values.", size=17, color=COLORS["muted"]))
    return parts


def panel_c(
    x: float,
    y: float,
    w: float,
    h: float,
    resid: pd.DataFrame,
    label: bool = True,
    row_h: float = 34,
    title_size: int = 25,
) -> list[str]:
    parts: list[str] = []
    parts.append(rect(x, y, w, h, "white", stroke="#E2E8EF", sw=1.2, rx=8))
    if label:
        parts.append(svg_text(x + 24, y + 43, "C", size=32, weight=700))
    parts += bar_panel(
        x + 78,
        y + 46,
        w - 112,
        "Spine-related genes in q < 0.01 BM set",
        resid,
        "z",
        "q_empirical",
        -8.5,
        7.5,
        "BM residual z-score in M. spicilegus brain",
        star_significance=True,
        show_grid=False,
        italicize_species_axis=True,
        row_h=row_h,
        title_size=title_size,
        gene_col="gene_name",
    )
    parts.append(svg_text(x + w - 24, y + h - 18, "* q < 0.05, ** q < 0.01, *** q < 0.001", size=13, color=COLORS["muted"], anchor="end"))
    return parts


def panel_d(x: float, y: float, w: float, h: float, label: bool = True) -> list[str]:
    parts: list[str] = []
    parts.append(rect(x, y, w, h, COLORS["panel"], stroke="#E2E8EF", sw=1.2, rx=8))
    if label:
        parts.append(svg_text(x + 26, y + 43, "D", size=32, weight=700))
    parts.append(svg_text(x + 78, y + 43, "Mechanistic interpretation", size=26, weight=700))

    box_y = y + 92
    box_h = 90
    box_w = 315
    gap = 58
    boxes = [
        (x + 85, "Reduced semaphorin cue", "Sema3b lower", COLORS["red"]),
        (x + 85 + (box_w + gap), "Reduced complement/microglia pruning", "C4a/C4b, Tyrobp, Grn lower", COLORS["red"]),
        (x + 85 + 2 * (box_w + gap), "Altered ephrin remodeling", "Efna1/Efna5/Efnb1 higher", COLORS["blue"]),
        (x + 85 + 3 * (box_w + gap), "More retained spines", "higher density in LD and SD", COLORS["green"]),
    ]
    for i, (box_x, title, subtitle, color) in enumerate(boxes):
        parts.append(rect(box_x, box_y, box_w, box_h, "white", stroke=color, sw=2.2, rx=7))
        parts.append(svg_text(box_x + 18, box_y + 34, title, size=21, weight=700, color=color))
        parts.append(svg_text(box_x + 18, box_y + 65, subtitle, size=18, color=COLORS["muted"]))
        if i < len(boxes) - 1:
            marker = "arrowBlue" if i == 2 else "arrowRed"
            arrow_color = COLORS["blue"] if i == 2 else COLORS["red"]
            parts.append(arrow(box_x + box_w + 10, box_y + 45, box_x + box_w + gap - 8, box_y + 45, color=arrow_color, sw=3, marker=marker))

    parts.append(
        svg_text(
            x + 85,
            y + 212,
            "Hypothesis: elevated spine density reflects less synapse elimination plus remodeling/stabilization, rather than a global pruning-gene shift.",
            size=20,
            color=COLORS["ink"],
        )
    )
    return parts


def main() -> None:
    OUTDIR.mkdir(exist_ok=True)

    joined = pd.read_csv(JOINED)
    bm_q01 = pd.read_csv(BM_Q01_STATS)

    expression_genes = ["Sema3b", "Shank3", "Efnb2", "Sema4d", "Mertk", "Grn"]
    residual_genes = ["Sema3b", "Grn", "Cdkl5", "Dlgap1", "Nectin3", "Tyrobp"]

    expr = joined[joined["Genename"].isin(expression_genes)].copy()
    expr["panel"] = "B_expression"
    expr["metric"] = "Spicilegus minus other species mean expression"
    expr = expr.sort_values("spic_minus_other_mean")

    resid = bm_q01[bm_q01["gene_name"].isin(residual_genes)].copy()
    resid["panel"] = "C_BM_q01_residual"
    resid["metric"] = "Spicilegus BM residual z-score from q < 0.01 source"
    resid = resid.sort_values("z")
    resid[
        [
            "gene_name",
            "observed",
            "predicted",
            "residual",
            "z",
            "q_empirical",
            "direction",
            "abs_z_margin_over_next_species",
        ]
    ].to_csv(OUT_PANEL_C_DATA, index=False)

    source_rows = []
    for row in expr.itertuples(index=False):
        source_rows.append(
            {
                "panel": "expression",
                "Genename": row.Genename,
                "value": row.spic_minus_other_mean,
                "q_value": row.spic_minus_other_mean_q_bh,
                "interpretation": "lower pruning cue" if row.Genename == "Sema3b" else "candidate remodeling signal",
            }
        )
    for row in resid.itertuples(index=False):
        source_rows.append(
            {
                "panel": "brownian_residual_q01_source",
                "Genename": row.gene_name,
                "value": row.z,
                "q_value": row.q_empirical,
                "interpretation": row.direction,
            }
        )
    pd.DataFrame(source_rows).to_csv(OUT_DATA, index=False)

    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
        arrow_marker_defs(),
        rect(0, 0, W, H, "white"),
        svg_text(70, 62, "Selective reduction of synapse-pruning cues may support elevated spine density in M. spicilegus", size=34, weight=700),
        svg_text(70, 98, "Working model from brain expression and Brownian residual analyses", size=22, color=COLORS["muted"]),
    ]

    # Panel A: phenotype anchor.
    ax, ay, aw, ah = 70, 140, 1660, 250
    svg += panel_a(ax, ay, aw, ah)

    # Panel B/C: data bars.
    bx, by, bw, bh = 70, 430, 805, 445
    cx, cy, cw, ch = 925, 430, 805, 445
    svg += panel_b(bx, by, bw, bh, expr)
    svg += panel_c(cx, cy, cw, ch, resid)

    # Panel D: mechanism model.
    dx, dy, dw, dh = 70, 925, 1660, 235
    svg += panel_d(dx, dy, dw, dh)

    svg.append("</svg>")
    OUT_SVG.write_text("\n".join(svg))

    OUT_PANEL_A.write_text(svg_doc(1660, 250, panel_a(0, 0, 1660, 250)))
    OUT_PANEL_B.write_text(svg_doc(805, 445, panel_b(0, 0, 805, 445, expr)))
    OUT_PANEL_C.write_text(svg_doc(600, 650, panel_c(0, 0, 600, 650, resid, row_h=58, title_size=22)))
    OUT_PANEL_D.write_text(svg_doc(1660, 235, panel_d(0, 0, 1660, 235)))


if __name__ == "__main__":
    main()
