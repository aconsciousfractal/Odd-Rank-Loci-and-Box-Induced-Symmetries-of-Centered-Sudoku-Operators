"""fig3_D6_cayley.py
=====================
Fig. 3 (paper II §5) — Cayley diagram of the dihedral group $D_6$ of
order $12$ inside $B_3 < S_7$, with generators $r$ (rotation, order 6)
and $s$ (reflection, order 2). Drawn as the standard hexagonal-prism
Cayley graph.

Frozen seed: vertex placement is deterministic; ``RANDOM_SEED = 0``
controls only label jitter (set to 0 in the shipped figure).
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch

RANDOM_SEED = 0
HERE = Path(__file__).resolve().parent


def main() -> Path:
    np.random.seed(RANDOM_SEED)

    # Inner hexagon: identity ring  e, r, r^2, r^3, r^4, r^5.
    # Outer hexagon: reflected ring s, sr, sr^2, sr^3, sr^4, sr^5.
    R_in, R_out = 1.0, 2.0
    inner = [(R_in * math.cos(k * math.pi / 3 + math.pi / 2),
              R_in * math.sin(k * math.pi / 3 + math.pi / 2))
             for k in range(6)]
    outer = [(R_out * math.cos(k * math.pi / 3 + math.pi / 2),
              R_out * math.sin(k * math.pi / 3 + math.pi / 2))
             for k in range(6)]

    inner_lbl = ["$e$", "$r$", "$r^2$", "$r^3$", "$r^4$", "$r^5$"]
    outer_lbl = ["$s$", "$sr$", "$sr^2$", "$sr^3$", "$sr^4$", "$sr^5$"]

    fig, ax = plt.subplots(figsize=(5.6, 5.6))

    # r-edges (red, directed): k -> k+1 on inner; on outer the convention
    # for sr^k * r = sr^{k+1}.
    for k in range(6):
        for ring in (inner, outer):
            x1, y1 = ring[k]
            x2, y2 = ring[(k + 1) % 6]
            ax.add_patch(FancyArrowPatch(
                (x1, y1), (x2, y2),
                arrowstyle="-|>", mutation_scale=10,
                color="#d62728", lw=1.2,
                shrinkA=10, shrinkB=10,
            ))

    # s-edges (blue, undirected): inner k <-> outer k.
    for k in range(6):
        x1, y1 = inner[k]
        x2, y2 = outer[k]
        ax.plot([x1, x2], [y1, y2], color="#1f77b4", lw=1.2)

    # Vertices and labels.
    for (x, y), lbl in zip(inner, inner_lbl):
        ax.plot(x, y, "o", markersize=18, color="white",
                markeredgecolor="black", zorder=4)
        ax.text(x, y, lbl, ha="center", va="center", fontsize=10, zorder=5)
    for (x, y), lbl in zip(outer, outer_lbl):
        ax.plot(x, y, "o", markersize=20, color="white",
                markeredgecolor="black", zorder=4)
        ax.text(x, y, lbl, ha="center", va="center", fontsize=10, zorder=5)

    ax.set_aspect("equal")
    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-2.5, 2.6)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)

    # Legend.
    ax.plot([], [], "-", color="#d62728", lw=1.2, label=r"$r$ (order 6)")
    ax.plot([], [], "-", color="#1f77b4", lw=1.2, label=r"$s$ (order 2)")
    ax.legend(loc="upper right", frameon=False, fontsize=9)

    fig.suptitle(
        r"Fig. 3.  Cayley diagram of $D_6 < B_3 < S_7$"
        r" (the HV1 dihedral $D_6$ of paper II \S5)",
        fontsize=11, y=0.97,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))

    out_pdf = HERE / "fig3_D6_cayley.pdf"
    out_png = HERE / "fig3_D6_cayley.png"
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f"[OK] {out_pdf.name}, {out_png.name}")
    return out_pdf


if __name__ == "__main__":
    main()
