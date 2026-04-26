r"""fig2_rank_distribution.py
=============================
Fig. 2 (paper II §4) — Rank distribution of the family
$\{E_\gamma\}_{\gamma \in S_9}$ on Sudoku base$_1$. Single bar plot,
log-scale on the y-axis, showing 362858 rank-8 vs 22 rank-7 relabelings.

Source: ``certified/base1_odd_rank_22.json`` (repo-root sibling of ``paper/``).
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
# HERE = paper/figures/, HERE.parents[1] = repo root, then /certified/.
CERT = HERE.parents[1] / "certified" / "base1_odd_rank_22.json"


def main() -> Path:
    with CERT.open("r", encoding="utf-8") as f:
        cert = json.load(f)
    dist = cert["outputs"]["rank_distribution"]
    # Sort by rank ascending.
    items = sorted(((int(k), int(v)) for k, v in dist.items()))
    ranks = [r for r, _ in items]
    counts = [c for _, c in items]

    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    bars = ax.bar(
        [str(r) for r in ranks], counts,
        color=["#d62728" if r == 7 else "#1f77b4" for r in ranks],
        edgecolor="black", linewidth=0.7,
    )
    ax.set_yscale("log")
    ax.set_xlabel(r"rank of $E_\gamma$ on $V_{\mathrm{std},9}$")
    ax.set_ylabel(r"$\#\{\gamma \in S_9\}$  (log scale)")
    ax.set_title(
        r"Fig. 2.  Rank distribution on base$_1$"
        r"  ($|S_9| = 362\,880$)",
        fontsize=11,
    )
    for b, c in zip(bars, counts):
        ax.text(b.get_x() + b.get_width() / 2,
                c * 1.25, f"{c:,}",
                ha="center", va="bottom", fontsize=9)
    ax.set_ylim(0.7, max(counts) * 5)
    ax.grid(axis="y", which="both", linestyle=":", linewidth=0.5, alpha=0.6)
    ax.set_axisbelow(True)
    fig.tight_layout()

    out_pdf = HERE / "fig2_rank_distribution.pdf"
    out_png = HERE / "fig2_rank_distribution.png"
    fig.savefig(out_pdf, metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f"[OK] {out_pdf.name}, {out_png.name}")
    return out_pdf


if __name__ == "__main__":
    main()
