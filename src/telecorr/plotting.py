"""Shared figure style: one categorical palette in a fixed order, thin marks, recessive grid, direct labels.

Utility: figure style of every notebook.

Colours follow the entity, never its rank: a condition keeps its colour in every figure of the repository.
"""
from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt

#: categorical slots in their fixed order (a validated colour-blind-safe palette, light surface)
SLOTS: tuple[str, ...] = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
#: every experimental condition and its colour, used the same way in all notebooks
CONDITION: dict[str, str] = {
    "TOMO-9 correction": SLOTS[0],
    "no correction": SLOTS[1],
    "native (no injected noise)": SLOTS[2],
    "oracle (exact inverse)": SLOTS[3],
    "DD": SLOTS[4],
    "TOMO-9 + DD": SLOTS[5],
    "stored correction (4 Oct)": SLOTS[6],
    # simulation steps 1-7: the same entities keep the same colours
    "noiseless": SLOTS[2],
    "learned from counts": SLOTS[0],
    "ceiling": SLOTS[3],
    "Bob-side correction": SLOTS[6],
    "single start": SLOTS[4],
    "two-qubit correction": SLOTS[7],
}
#: the two test messages, when a figure plots them as separate series
MESSAGE: dict[str, str] = {"product_HT": SLOTS[6], "bell_30_70": SLOTS[4]}
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e4e3df"
SURFACE = "#fcfcfb"


def use_style() -> None:
    """Apply the repository's matplotlib style (call once per notebook)."""
    mpl.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "figure.dpi": 110, "savefig.dpi": 150, "font.size": 10,
        "axes.edgecolor": INK_2, "axes.labelcolor": INK, "xtick.color": INK_2, "ytick.color": INK_2,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "axes.axisbelow": True,
        "lines.linewidth": 2.0, "lines.markersize": 7, "legend.frameon": False,
        "axes.prop_cycle": mpl.cycler(color=list(SLOTS)),
    })


def bar_labels(ax, bars, fmt: str = "{:.3f}", pad: float = 0.005) -> None:
    """Value labels above bars, in text ink (not the series colour)."""
    for b in bars:
        h = b.get_height()
        ax.text(b.get_x() + b.get_width() / 2, h + pad, fmt.format(h), ha="center", va="bottom", fontsize=8.5, color=INK)


def new(figsize=(7.5, 3.6), **kw):
    use_style()
    return plt.subplots(figsize=figsize, **kw)
