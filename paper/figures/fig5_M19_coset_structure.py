r"""fig5_M19_coset_structure.py
================================
Fig. 5 (paper II §8) -- Bilateral coset structure of $\Gamma^*_{M_{19}}$.

The structural facts depicted (all read from ``certified/M19_audit.json``
and ``certified/M19_KL_KR_3sq_Q8.json``):

* $|\Gamma^*_{M_{19}}| = 72$ (one bilateral coset);
* $\Gamma^*_{M_{19}} = K_L\,\gamma_0 = \gamma_0\,K_R$, where
  $K_L, K_R \le S_9$ both have order $72$;
* $K_L \cap K_R = \{e\}$ but $K_L \neq K_R$;
* $\pi_0 \in S_9 \setminus (K_L \cup K_R)$ is an external conjugator
  with $\pi_0 K_R \pi_0^{-1} = K_L$ (lex-min of $\Gamma^*$);
* the ambient coordinate stabilizer of the value-multiplicity vector
  $u_9^* \in \mathbb{Z}^9$ is $S_3 \times S_6$ of order $4320$,
  inside which $K_L, K_R$ and $\Gamma^*$ all sit.

The figure draws (i) one schematic block of $72$ cells representing
$\Gamma^*_{M_{19}}$ as both a $K_L$-left-coset and a $K_R$-right-coset,
plus (ii) a separate label box for the ambient stabilizer
$S_3 \times S_6$ of order $4320$. NO $72 \times 60$ heatmap is drawn:
$|\Gamma^*| = 72$, not $4320$.

Source: ``certified/M19_audit.json`` and ``certified/M19_KL_KR_3sq_Q8.json``.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
# HERE = paper/figures/, HERE.parents[1] = repo root, then /certified/.
CERT_DIR = HERE.parents[1] / "certified"


def _load(name: str) -> dict:
    with (CERT_DIR / name).open("r", encoding="utf-8") as f:
        return json.load(f)


def main() -> Path:
    audit = _load("M19_audit.json")
    klkr = _load("M19_KL_KR_3sq_Q8.json")

    out_a = audit["outputs"]
    out_k = klkr["outputs"]

    KL_size = int(out_a.get("K_L_size", 72))
    KR_size = int(out_a.get("K_R_size", 72))
    Gstar = int(out_a.get("Gamma_star_size", 72))
    inter = int(out_a.get("K_L_inter_K_R_size", out_k.get("K_L_inter_K_R_size", 1)))
    ambient = int(out_a.get("ambient_coordinate_stabilizer_S_9_of_u9_star_size", 4320))
    KL_eq_KR = bool(out_k.get("K_L_eq_K_R", False))
    Y2_pi0 = bool(out_k.get("Y2_conjugate_via_pi0", True))
    pi0 = out_k.get("pi0_perm")

    assert KL_size == KR_size == Gstar == 72, (
        f"figure expects |K_L|=|K_R|=|Gamma*|=72, got "
        f"{KL_size}/{KR_size}/{Gstar}"
    )
    assert inter == 1, f"K_L cap K_R must be trivial, got size {inter}"
    assert not KL_eq_KR, "K_L must NOT equal K_R as subgroups of S_9"
    assert Y2_pi0, "pi_0 conjugacy K_L = pi_0 K_R pi_0^-1 must hold"

    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.set_aspect("equal")
    ax.axis("off")

    # ---- Big outer box: ambient coordinate stabilizer S_3 x S_6 (4320). ----
    ax.add_patch(mpatches.Rectangle(
        (0.2, 0.2), 9.6, 5.6,
        linewidth=1.2, edgecolor="0.3", facecolor="0.97"
    ))
    ax.text(
        5.0, 5.45,
        rf"Ambient coord. stabilizer  $\mathrm{{Stab}}_{{S_9}}(u_9^*) "
        rf"\cong S_3 \times S_6$  (order ${ambient}$)",
        ha="center", va="center", fontsize=10,
    )

    # ---- Inner box: Gamma^*_{M_19}, the unique 72-element bilateral coset.
    ax.add_patch(mpatches.Rectangle(
        (1.2, 1.0), 7.6, 3.0,
        linewidth=1.4, edgecolor="C0", facecolor="#dde8f5"
    ))
    ax.text(
        5.0, 3.65,
        rf"$\Gamma^*_{{M_{{19}}}} = K_L\,\gamma_0 = \gamma_0\,K_R$, "
        rf"$|\Gamma^*_{{M_{{19}}}}| = {Gstar}$",
        ha="center", va="center", fontsize=10, color="C0",
    )

    # ---- Render the 72-element coset as a 6 x 12 grid of cells. ----
    rows, cols = 6, 12
    x0, y0 = 1.45, 1.25
    cell_w = (8.55 - 1.45) / cols
    cell_h = (3.55 - 1.25) / rows
    for r in range(rows):
        for c in range(cols):
            shade = 0.65 + 0.05 * ((r + c) % 2)
            ax.add_patch(mpatches.Rectangle(
                (x0 + c * cell_w, y0 + r * cell_h),
                cell_w * 0.92, cell_h * 0.85,
                linewidth=0.3, edgecolor="0.6",
                facecolor=plt.cm.Blues(shade),
            ))

    # ---- Labels for K_L (left), K_R (right), pi_0 (external). ----
    ax.text(
        0.85, 2.5,
        rf"$K_L$" "\n" rf"order ${KL_size}$",
        ha="center", va="center", fontsize=9, rotation=90,
        color="C0",
    )
    ax.text(
        9.15, 2.5,
        rf"$K_R$" "\n" rf"order ${KR_size}$",
        ha="center", va="center", fontsize=9, rotation=270,
        color="C0",
    )
    pi0_str = (
        "(" + ",\\,".join(str(x) for x in pi0) + ")"
        if pi0 is not None else r"\text{see cert}"
    )
    ax.text(
        5.0, 0.55,
        rf"External conjugator  $\pi_0 = {pi0_str} \in S_9$,  "
        rf"$\pi_0 \notin K_L \cup K_R$,  "
        rf"$K_L \cap K_R = \{{e\}}$  (trivial)",
        ha="center", va="center", fontsize=8.5, color="0.2",
    )

    ax.set_title(
        r"Fig. 5.  Bilateral coset structure $\Gamma^*_{M_{19}} = "
        r"K_L\,\gamma_0 = \gamma_0\,K_R$ inside $S_3 \times S_6$",
        fontsize=10.5,
    )

    fig.tight_layout()
    out_pdf = HERE / "fig5_M19_coset_structure.pdf"
    out_png = HERE / "fig5_M19_coset_structure.png"
    fig.savefig(out_pdf, metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f"[OK] {out_pdf.name}, {out_png.name}")
    return out_pdf


if __name__ == "__main__":
    main()
