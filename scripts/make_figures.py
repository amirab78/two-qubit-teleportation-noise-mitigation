"""Make the figures of the README (``figures/``).

usage:  python scripts/make_figures.py

* ``protocol.png``        the six-qubit protocol with the channel and the correction box;
* ``hardware_result.png`` the main hardware result, recomputed from the raw counts in ``data/`` (about 20 s).
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                   # noqa: E402
import numpy as np                                                # noqa: E402
from matplotlib.patches import FancyBboxPatch                     # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(REPO, "src"))
from telecorr.plotting import CONDITION, INK, INK_2, SURFACE, bar_labels, use_style   # noqa: E402

OUT = os.path.join(REPO, "figures")


def _box(ax, x, y0, y1, text, face, ink="white", w=0.9):
    ax.add_patch(FancyBboxPatch((x - w / 2, y0 - 0.32), w, y1 - y0 + 0.64, boxstyle="round,pad=0.02,rounding_size=0.12",
                                fc=face, ec=face, zorder=3))
    ax.text(x, (y0 + y1) / 2, text, ha="center", va="center", color=ink, fontsize=9, zorder=4, linespacing=1.15)


def protocol() -> str:
    use_style()
    fig, ax = plt.subplots(figsize=(10.5, 3.9))
    ys = {q: 5 - q for q in range(6)}                              # q0 on top
    roles = {0: "q0  message", 1: "q1  message", 2: "q2  channel", 3: "q3  channel", 4: "q4  Bob", 5: "q5  Bob"}
    end = {0: 6.0, 1: 6.45, 2: 6.0, 3: 6.45, 4: 10.6, 5: 10.6}       # Alice's wires end at her measurement
    for q, y in ys.items():
        ax.plot([0.2, end[q]], [y, y], color=INK_2, lw=1.1, zorder=1)
        ax.text(0.1, y, roles[q], ha="right", va="center", fontsize=9.5, color=INK)
    _box(ax, 0.9, ys[1], ys[0], "prepare\nmessage", INK_2)
    for a, b, x in ((2, 4, 2.0), (3, 5, 2.45)):                    # Bell pairs
        ax.plot([x, x], [ys[b], ys[a]], color=INK, lw=1.4, zorder=2)
        for q in (a, b):
            ax.plot(x, ys[q], "o", color=INK, ms=6, zorder=3)
    ax.text(2.22, 5.55, "Bell pairs", ha="center", fontsize=9, color=INK)
    _box(ax, 3.6, ys[3], ys[2], "channel\nnoise", INK)
    _box(ax, 4.8, ys[3], ys[2], "correction\n(6 angles)", CONDITION["learned from counts"])
    for a, b, x in ((0, 2, 6.0), (1, 3, 6.45)):                    # Alice's Bell measurement
        ax.plot([x, x], [ys[b], ys[a]], color=INK, lw=1.4, zorder=2)
        for q in (a, b):
            ax.add_patch(FancyBboxPatch((x - 0.17, ys[q] - 0.2), 0.34, 0.4, boxstyle="round,pad=0.01", fc=SURFACE, ec=INK, zorder=3))
            ax.text(x, ys[q], "M", ha="center", va="center", fontsize=8.5, zorder=4)
    ax.text(6.22, 5.55, "Alice's Bell measurement", ha="center", fontsize=9, color=INK)
    ax.annotate("", xy=(7.85, 0.5), xytext=(6.7, 2.2), arrowprops=dict(arrowstyle="-|>", color=INK_2, lw=1.2, ls="--"))
    ax.text(7.05, 1.75, "4 bits", fontsize=8.5, color=INK_2, rotation=-38)
    _box(ax, 8.3, ys[5], ys[4], "X / Z\ncorrections", INK_2)
    _box(ax, 9.7, ys[5], ys[4], "measure\n(9 bases)", INK_2)
    ax.text(4.8, -0.75, "learned from probe counts only:\nthe receiver never knows the noise", ha="center", fontsize=8.5,
            color=CONDITION["learned from counts"])
    ax.set_xlim(-1.3, 10.8); ax.set_ylim(-1.2, 5.9); ax.axis("off")
    plt.tight_layout()
    path = os.path.join(OUT, "protocol.png")
    fig.savefig(path, dpi=160); plt.close(fig)
    return path


def hardware_result() -> str:
    from telecorr.hardware import analysis as A
    use_style()
    out = A.reproduce()
    r1, r2 = out["r1"], out["r2"]
    F = lambda k: (r1["tomo"].get(k) or r2["tomo"][k])["F"]        # noqa: E731
    S = lambda k: (r1["tomo"].get(k) or r2["tomo"][k])["F_sigma"]  # noqa: E731
    order = [("R1_N1none", "no correction"), ("R1_N1oracle", "oracle (exact inverse)"), ("R2_N1corr", "TOMO-9 correction"),
             ("R1_nativeA", "native (no injected noise)")]
    states = ("product_HT", "bell_30_70")
    fig, ax = plt.subplots(figsize=(7.6, 3.7)); w = 0.2
    for i, (pre, cond) in enumerate(order):
        bars = ax.bar(np.arange(2) + (i - 1.5) * w, [F(f"{pre}_{st}") for st in states], w * 0.92,
                      yerr=[S(f"{pre}_{st}") for st in states], color=CONDITION[cond], label=A.NAMES[pre], capsize=3, ecolor=INK_2)
        bar_labels(ax, bars, pad=0.012)
    ax.set_xticks(range(2), states); ax.set_ylim(0, 1.05); ax.set_ylabel("fidelity F")
    ax.set_title("ibm_kingston, 4 October 2026: an injected error N1, corrected from real counts", fontsize=10)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2, fontsize=8.5)
    plt.tight_layout()
    path = os.path.join(OUT, "hardware_result.png")
    fig.savefig(path, dpi=160); plt.close(fig)
    return path


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for f in (protocol, hardware_result):
        print("written", f())
