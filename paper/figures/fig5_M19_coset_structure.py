r"""fig5_M19_coset_structure.py
================================
Fig. 5 (paper II §8) — Bilateral coset structure of $\Gamma^*_{M_{19}}$
under the Hessian symmetry groups $K_L$ (left) and $K_R$ (right). The
group acts as $K_L \times K_R$ on $\Gamma^*_{M_{19}}$ and partitions
$|\Gamma^*_{M_{19}}| = 4320$ into $|K_L| \cdot |K_R| / |\Stab| = 4320$
labelled cells; we draw a $|K_L| = 72$ row × $|K_R| / |\Stab| = 60$
column heatmap-style schematic.

Source: ``paper/certified/M19_audit.json`` (numbers).
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

RANDOM_SEED = 0
HERE = Path(__file__).resolve().parent
CERT = HERE.parent / "certified" / "M19_audit.json"


def main() -> Path:
    np.random.seed(RANDOM_SEED)

    with CERT.open("r", encoding="utf-8") as f:
        cert = json.load(f)
    out = cert["outputs"]
    KL = int(out.get("K_L_size", 72))
    stab = int(out.get("stab_size", 4320))
    total = stab  # |Gamma*_{M_19}| in the certificate semantics
    cols = max(1, total // KL)  # 60 columns

    # Build a deterministic block matrix where each KL-row is a coset.
    grid = np.zeros((KL, cols))
    for r in range(KL):
        # alternating shading to make the KL rows visually distinct.
        grid[r, :] = (r % 6) * 0.15 + 0.05

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    im = ax.imshow(
        grid, aspect="auto", cmap="Blues", vmin=0, vmax=1,
        interpolation="nearest", origin="lower",
    )
    # KR-block separators every 6 columns (60 = 6 * 10).
    for k in range(0, cols + 1, 6):
        ax.axvline(k - 0.5, color="0.4", lw=0.4)
    # KL-block separators every 6 rows (72 = 6 * 12).
    for k in range(0, KL + 1, 6):
        ax.axhline(k - 0.5, color="0.4", lw=0.4)

    ax.set_xticks([0, cols // 2, cols - 1])
    ax.set_yticks([0, KL // 2, KL - 1])
    ax.set_xlabel(
        rf"$K_R$-coset index  (1 \dots {cols} = $|\Gamma^*_{{M_{{19}}}}| / |K_L| = {cols}$)",
        fontsize=9,
    )
    ax.set_ylabel(rf"$K_L$ row index (1 \dots {KL} = $|K_L|$)", fontsize=9)
    ax.set_title(
        r"Fig. 5.  Bilateral coset structure of $\Gamma^*_{M_{19}}$"
        rf"  ($|K_L|={KL}$, total $={total}$)",
        fontsize=11,
    )
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02,
                 ticks=[0, 1], label="$K_L$ row band")

    fig.tight_layout()
    out_pdf = HERE / "fig5_M19_coset_structure.pdf"
    out_png = HERE / "fig5_M19_coset_structure.png"
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f"[OK] {out_pdf.name}, {out_png.name}")
    return out_pdf


if __name__ == "__main__":
    main()
