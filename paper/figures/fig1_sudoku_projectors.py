r"""fig1_sudoku_projectors.py
============================
Fig. 1 (paper II §3) — A 9x9 Sudoku grid with the row-band ($P_{\mathrm{rb}}$)
and column-stack ($P_{\mathrm{cs}}$) projector overlays. The figure shows
that the Box Band Lemma identity $B_{\mathrm{rb}}^\top E B_{\mathrm{cs}} = 0$
acts on a basis adapted to the 3x3 box decomposition.

Frozen seed: ``RANDOM_SEED = 0`` (only used if the script is asked to
fill in a sample grid; the figure as shipped uses base1).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

RANDOM_SEED = 0
HERE = Path(__file__).resolve().parent

# base1 grid (same as in the certified package).
BASE1 = np.array([
    [5, 3, 4, 6, 7, 8, 9, 1, 2],
    [6, 7, 2, 1, 9, 5, 3, 4, 8],
    [1, 9, 8, 3, 4, 2, 5, 6, 7],
    [8, 5, 9, 7, 6, 1, 4, 2, 3],
    [4, 2, 6, 8, 5, 3, 7, 9, 1],
    [7, 1, 3, 9, 2, 4, 8, 5, 6],
    [9, 6, 1, 5, 3, 7, 2, 8, 4],
    [2, 8, 7, 4, 1, 9, 6, 3, 5],
    [3, 4, 5, 2, 8, 6, 1, 7, 9],
], dtype=np.int64)

# Three row-bands and three column-stacks (each a union of 3 rows / cols).
RB_COLORS = ["#cce5ff", "#ffe4b5", "#d5f5d5"]
CS_HATCH = ["//", "\\\\", "xx"]


def main() -> Path:
    np.random.seed(RANDOM_SEED)
    fig, ax = plt.subplots(figsize=(5.6, 5.6))

    # Row-band background colors (3 horizontal bands).
    for k in range(3):
        ax.add_patch(Rectangle(
            (0, 9 - 3 * (k + 1)), 9, 3,
            facecolor=RB_COLORS[k], edgecolor="none", zorder=0,
        ))

    # Column-stack hatching (3 vertical stacks).
    for k in range(3):
        ax.add_patch(Rectangle(
            (3 * k, 0), 3, 9,
            facecolor="none", edgecolor="0.55",
            hatch=CS_HATCH[k], linewidth=0, zorder=1,
        ))

    # Cell entries.
    for i in range(9):
        for j in range(9):
            ax.text(
                j + 0.5, 9 - i - 0.5, str(int(BASE1[i, j])),
                ha="center", va="center", fontsize=11, color="#222",
                zorder=3,
            )

    # Thin grid lines.
    for k in range(10):
        lw = 1.6 if k % 3 == 0 else 0.5
        ax.plot([k, k], [0, 9], color="black", lw=lw, zorder=2)
        ax.plot([0, 9], [k, k], color="black", lw=lw, zorder=2)

    ax.set_xlim(-0.05, 9.05)
    ax.set_ylim(-0.05, 9.05)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)

    # Legend.
    rb_patches = [Rectangle((0, 0), 1, 1, facecolor=c, edgecolor="black",
                             label=f"row-band $R_{k+1}$")
                  for k, c in enumerate(RB_COLORS)]
    cs_patches = [Rectangle((0, 0), 1, 1, facecolor="white",
                             edgecolor="0.55", hatch=h,
                             label=f"col-stack $C_{k+1}$")
                  for k, h in enumerate(CS_HATCH)]
    ax.legend(handles=rb_patches + cs_patches,
              loc="upper center", bbox_to_anchor=(0.5, -0.02),
              ncol=3, frameon=False, fontsize=8)

    fig.suptitle(
        r"Fig. 1.  Sudoku base$_1$ with row-band and column-stack overlays",
        fontsize=11, y=0.99,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))

    out_pdf = HERE / "fig1_sudoku_projectors.pdf"
    out_png = HERE / "fig1_sudoku_projectors.png"
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f"[OK] {out_pdf.name}, {out_png.name}")
    return out_pdf


if __name__ == "__main__":
    main()
