# Script index — cryptic legacy names → meaningful purpose

The producer scripts under [`legacy/`](legacy) retain phase-numbered
filenames inherited from the research log. Their purpose, in
plain English, is recorded below.

| Legacy filename | Plain-English purpose | Output(s) under `data/` |
|---|---|---|
| `phase_9_13.py`                       | Exhaustive $S_9$ symbol-relabeling scan of the LS-9 base `base1`; counts ranks of $E_\gamma = \gamma(L) - 5J$ over $V_{\mathrm{std}}$ and lists the 22 rank-7 perms. | `phase_9_13_results.json` |
| `phase_9_13f_certificate.py`          | Per-perm certificate for the 22 base1 odd-rank perms (kernel witness, integer SNF). | `phase_9_13f_certificate.json` |
| `phase_9_14.py`                       | $n=6$ odd-rank census: full-rank counts of $E_\gamma$ over $\mathbb{F}_2$ for Sudoku-6 vs Latin-Square-6, first-row-fixed; produces the two-proportion $z$-statistic. | `phase_9_14_results.json` |
| `phase_9_15.py`                       | $n=7$ box-permutation $B_3 = S_3 \wr S_3$ scan on the laboratory 7×7 Latin square $L_{V_1}$; finds the rank-5 locus $H_{V_1} \subset B_3$ of order 12. | `phase_9_15_results.json` |
| `phase_9_15A_identify_HV1.py`         | Group-theoretic identification of $H_{V_1}$ as $D_6$ (element-order multiset, dihedral relation, $H_{V_1} \subset B_3$). | `phase_9_15A_H_V1_identification.json` |
| `phase_9_18_S9c.py`                   | $S_9$ stabilizer combinatorics for $L_{M_{19}}$; produces the rank-7 candidate set used downstream. | `phase_9_18_S9c_results.json` |
| `phase_9_18_S9x.py`                   | $M_{19}$ rank-7 scan stage 1 (NumPy filter) + stage 2 (Sympy exact rank); enumerates the 100 rank-7 symbol-relabelings of $L_{M_{19}}$. | `phase_9_18_S9x_results.json` |
| `phase_9_18_S9x_bis.py`               | Four-way certificate (rank $\mathbb{Z}$, primitive $V_{\mathrm{std}}$ kernel, integer SNF, $\mathbb{F}_2$ SNF) for the 100 $M_{19}$ rank-7 perms. | `phase_9_18_S9x_bis_results.json` |
| `phase_9_18_S9q_C1.py`                | Coset / quotient analysis: $\Gamma^{\star} = K_L \cdot \gamma_0 = \gamma_0 \cdot K_R$ at $L_{M_{19}}$. | `phase_9_18_S9q_C1.json` |
| `phase_9_18_S9q_beta_lattice.py`      | $\beta$-lattice inclusion $\Lambda_{O_6} \subseteq \Lambda_M$ (Box-Band derived). | `phase_9_18_S9q_beta_lattice.json` |
| `phase_9_18_S9q_gamma1.py`            | $\gamma_1$ stabilizer projection / orbit decomposition at $M_{19}$. | `phase_9_18_S9q_gamma1.json` |
| `phase_9_18_extra_canonical_ranks.py` | Extra canonical-rank distribution check on $N=400$ random Sudoku grids (frozen seed). | `phase_9_18_extra_canonical_ranks.json` |
| `phase13b_exhaustive.py`              | Phase-13 exhaustive $D_4$ / boundary scan (research artifact, not directly cited). | — |
| `task_9_10c_box_subspace.py`          | Box-subspace projector identity test (`P_rb · E · P_cs = 0` symbolic check). | — |
| `check_b3_subgroup.py`                | One-page sanity check that $H_{V_1}$ is closed under $B_3$ composition. | — |

## Why are these names kept?

Two reasons:

1. **Provenance integrity.** The certificate JSONs in `certified/`
   record sha256 hashes of these files under the `scripts/legacy/...`
   path. Renaming them would break the standalone audit unless every
   certificate is rebuilt; we prefer to lock provenance and document
   the names.
2. **Cross-imports.** Several scripts import each other
   (`phase_9_18_S9x.py` imports `phase_9_18_S9c`, `phase_9_18_S9q_*`
   imports `phase_9_13`, `phase_9_18_extra_canonical_ranks` imports
   `phase_9_13`). Renaming would require fixing every import site
   without affecting numerical output.

## What you actually need to run

You almost never need to run the legacy scripts directly. The
canonical entry points are:

- `python certified/build_certified.py` — builds the certificates from
  the already-generated `data/*.json`.
- `python certified/verify_all.py` — verifies everything.
- `python certified/reproduce_all.py` — re-derives certificates byte-equally.

The legacy scripts under `scripts/legacy/` are only needed if you want
to regenerate `data/*.json` from scratch.
