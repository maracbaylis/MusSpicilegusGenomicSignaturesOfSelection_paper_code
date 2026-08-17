#!/usr/bin/env python3
from pathlib import Path
import os

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib_cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle


OUTDIR = ROOT / "phyloP_MSA_QC/go_enrichment_spicilegus_combined_fg_allcCRE_background"

NAVY = "#0b2f5b"
BLUE = "#08519c"
LIGHT_BLUE = "#5da5da"
GREEN = "#006d2c"
MID_GREEN = "#31a354"
LIGHT_GREEN = "#a1d99b"
ORANGE = "#d95f02"
GRAY = "#6b7280"
LIGHT_GRAY = "#f3f6f8"
TEXT = "#1f2933"

MODULES = [
    {
        "title": "Neuromodulatory\ninput",
        "genes": ["ADRB1", "TACR1", "PGRMC1"],
        "body": "Arousal, salience,\nstress/state regulation",
        "color": BLUE,
    },
    {
        "title": "Intracellular signaling\nand calcium pathways",
        "genes": ["PLCB1", "ITGA8", "SLC2A4"],
        "body": "Second messengers,\nCa2+ signaling,\nsynapse stabilization",
        "color": GREEN,
    },
    {
        "title": "Activity-dependent\ntranscription",
        "genes": ["NPAS4", "PLK2"],
        "body": "Activity response,\ninhibitory synapse tuning,\nhomeostatic scaling",
        "color": MID_GREEN,
    },
    {
        "title": "Synaptic remodeling\nand vesicle cycling",
        "genes": ["SYT11"],
        "body": "Vesicle trafficking,\nrelease dynamics,\nsynaptic output",
        "color": LIGHT_BLUE,
    },
]


def rounded_box(ax, xy, width, height, color, face=LIGHT_GRAY, radius=0.08, lw=1.5):
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle=f"round,pad=0.025,rounding_size={radius}",
        linewidth=lw,
        edgecolor=color,
        facecolor=face,
    )
    ax.add_patch(patch)
    return patch


def arrow(ax, start, end, color=GRAY, lw=2.0, mutation_scale=14):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=mutation_scale,
            linewidth=lw,
            color=color,
            shrinkA=4,
            shrinkB=4,
            connectionstyle="arc3,rad=0.0",
        )
    )


def gene_pill(ax, x, y, label, color, width=0.55):
    pill = FancyBboxPatch(
        (x - width / 2, y - 0.09),
        width,
        0.18,
        boxstyle="round,pad=0.02,rounding_size=0.08",
        linewidth=0,
        facecolor=color,
        alpha=0.95,
    )
    ax.add_patch(pill)
    ax.text(x, y, label, ha="center", va="center", fontsize=9.5, color="white", fontweight="bold")


def draw_receptor_icon(ax, x, y, color):
    ax.plot([x - 0.35, x + 0.35], [y, y], color=color, linewidth=3)
    for dx in [-0.18, 0.0, 0.18]:
        ax.add_patch(Rectangle((x + dx - 0.025, y - 0.22), 0.05, 0.44, color=color, alpha=0.9))
    for dx in [-0.25, 0.05, 0.28]:
        ax.add_patch(Circle((x + dx, y + 0.42), 0.055, color=ORANGE, alpha=0.9))
    arrow(ax, (x + 0.1, y + 0.34), (x + 0.03, y + 0.12), color=ORANGE, lw=1.6, mutation_scale=10)


def draw_signaling_icon(ax, x, y, color):
    ax.add_patch(Circle((x - 0.22, y + 0.08), 0.10, edgecolor=color, facecolor="white", linewidth=2))
    ax.add_patch(Circle((x + 0.03, y - 0.06), 0.10, edgecolor=color, facecolor="white", linewidth=2))
    ax.add_patch(Circle((x + 0.28, y + 0.10), 0.10, edgecolor=color, facecolor="white", linewidth=2))
    arrow(ax, (x - 0.12, y + 0.04), (x - 0.06, y - 0.02), color=color, lw=1.5, mutation_scale=9)
    arrow(ax, (x + 0.13, y - 0.02), (x + 0.20, y + 0.04), color=color, lw=1.5, mutation_scale=9)
    ax.text(x, y + 0.36, "Ca2+", ha="center", va="center", fontsize=12, color=color, fontweight="bold")


def draw_nucleus_icon(ax, x, y, color):
    ax.add_patch(Circle((x, y), 0.33, edgecolor=color, facecolor="white", linewidth=2))
    ax.add_patch(Circle((x + 0.07, y + 0.03), 0.08, color=color, alpha=0.22))
    ax.text(x, y, "TF", ha="center", va="center", fontsize=12, color=color, fontweight="bold")
    for i, yy in enumerate([y - 0.43, y - 0.53, y - 0.63]):
        ax.plot([x - 0.32, x + 0.32], [yy, yy], color=color, linewidth=1.6, alpha=0.7 - i * 0.12)


def draw_synapse_icon(ax, x, y, color):
    ax.add_patch(Circle((x - 0.22, y), 0.30, edgecolor=color, facecolor="white", linewidth=2))
    ax.add_patch(Rectangle((x + 0.12, y - 0.30), 0.08, 0.60, color=color, alpha=0.85))
    for dx, dy in [(-0.30, 0.08), (-0.18, 0.16), (-0.08, -0.02)]:
        ax.add_patch(Circle((x + dx, y + dy), 0.045, color=ORANGE, alpha=0.9))
    for dy in [-0.15, 0.0, 0.15]:
        arrow(ax, (x - 0.02, y + dy), (x + 0.12, y + dy), color=ORANGE, lw=1.2, mutation_scale=8)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(16.5, 6.2))
    ax.set_xlim(0, 13.2)
    ax.set_ylim(0, 5.4)
    ax.axis("off")

    ax.text(
        0.25,
        5.65,
        "Accelerated brain cCRE-associated genes converge on memory-related neuronal plasticity",
        ha="left",
        va="center",
        fontsize=17,
        fontweight="bold",
        color=TEXT,
    )
    ax.text(
        0.25,
        5.25,
        "GO:0007613 memory | ADRB1, ITGA8, NPAS4, PGRMC1, PLCB1, PLK2, SLC2A4, SYT11, TACR1",
        ha="left",
        va="center",
        fontsize=10.5,
        color=GRAY,
    )

    # Upstream regulatory context.
    rounded_box(ax, (0.25, 4.25), 1.35, 0.47, ORANGE, face="#fff7ed", radius=0.07, lw=1.2)
    ax.text(0.925, 4.49, "Accelerated\nbrain cCREs", ha="center", va="center", fontsize=11, color=TEXT, fontweight="bold")
    for x in [0.43, 0.58, 0.73, 1.36, 1.21, 1.06]:
        ax.add_patch(Circle((x, 4.18), 0.035, color=ORANGE, alpha=0.75))

    y0 = 1.55
    width = 2.2
    gap = 0.28
    x0 = 0.55
    icon_y = 3.63

    icon_funcs = [draw_receptor_icon, draw_signaling_icon, draw_nucleus_icon, draw_synapse_icon]
    module_centers = []
    for i, module in enumerate(MODULES):
        x = x0 + i * (width + gap)
        module_centers.append((x + width / 2, y0 + 1.35))
        rounded_box(ax, (x, y0), width, 2.55, module["color"], face="#f8fafc", radius=0.09)
        icon_funcs[i](ax, x + width / 2, icon_y, module["color"])
        ax.text(
            x + width / 2,
            y0 + 1.82,
            module["title"],
            ha="center",
            va="center",
            fontsize=10.8,
            fontweight="bold",
            color=TEXT,
            linespacing=1.05,
        )
        gene_y = y0 + 1.22
        for j, gene in enumerate(module["genes"]):
            gene_pill(ax, x + width / 2, gene_y - j * 0.23, gene, module["color"], width=0.66)
        ax.text(
            x + width / 2,
            y0 + 0.38,
            module["body"],
            ha="center",
            va="center",
            fontsize=8.9,
            color=TEXT,
            linespacing=1.15,
        )

    for (sx, sy), (ex, ey) in zip(module_centers[:-1], module_centers[1:]):
        arrow(ax, (sx + width / 2 - 0.03, sy), (ex - width / 2 + 0.03, ey), color=GRAY, lw=2.0)

    # cCRE arrow to pathway.
    arrow(ax, (1.58, 4.38), (2.0, 4.05), color=ORANGE, lw=2.0, mutation_scale=12)

    out_x = x0 + 4 * (width + gap) + 0.05
    rounded_box(ax, (out_x, 1.75), 1.95, 2.15, NAVY, face="#eff6ff", radius=0.09, lw=1.6)
    ax.text(out_x + 0.975, 3.36, "Memory-related\nsynaptic plasticity", ha="center", va="center", fontsize=11.8, fontweight="bold", color=TEXT)
    ax.text(out_x + 0.975, 2.70, "LTP and homeostatic\nregulation", ha="center", va="center", fontsize=9.4, color=TEXT)
    ax.text(out_x + 0.975, 2.18, "Experience-dependent\ncircuit remodeling", ha="center", va="center", fontsize=9.4, color=TEXT)
    arrow(ax, (module_centers[-1][0] + width / 2 - 0.03, module_centers[-1][1]), (out_x + 0.02, module_centers[-1][1]), color=GRAY, lw=2.2)

    ax.text(
        0.55,
        0.70,
        "Core logic: neuronal activity -> neuromodulatory signaling -> Ca2+/second messenger pathways -> activity-dependent transcription -> synaptic remodeling",
        ha="left",
        va="center",
        fontsize=10.2,
        color=TEXT,
    )
    ax.text(
        0.55,
        0.36,
        "Interpretation: accelerated cCREs are associated with genes spanning multiple levels of memory-related neuronal plasticity.",
        ha="left",
        va="center",
        fontsize=10.2,
        color=TEXT,
    )

    for ext in ["png", "pdf", "svg"]:
        fig.savefig(OUTDIR / f"memory_plasticity_circuit_schematic.{ext}", dpi=450, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
