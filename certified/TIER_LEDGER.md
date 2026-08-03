# TIER_LEDGER.md - Claim-by-claim tier discipline

> **Purpose.** Single canonical ledger that maps every claim in the
> paper to (i) its public source artifact/script, (ii) the certified
> JSON in `certified/`, and (iii) the **tier** at which the paper
> exposes it. Every referenced artifact is present in this repository or
> linked to a public companion repository.

> **Tier alphabet (6 levels).**  THEOREM | CERTIFIED_EXHAUSTIVE |
> EXHAUSTIVE | EMPIRICAL | OBSERVED TEMPLATE | CONJECTURE.
> `CERTIFIED_EXHAUSTIVE` means a finite exact-arithmetic computation with
> a checked certificate and reproducible JSON proof object.  Anything in
> the source corpus that lacks an explicit tier
> tag, or is stated at a higher tier than its evidence supports, is
> **demoted** in the row below and a justification line is given in §B.

> **Filter applied.** Every numerical or structural claim in the paper
> must come from one of the rows in §A. Statements not supported at that
> level are publicly bounded or excluded in §B.

---

## A. Section-to-source-to-certificate ledger

| § | Claim (paper-level statement) | Paper label | Tier in paper | Primary source(s) | Producer script(s) | Certified JSON |
|---|---|---|---|---|---|---|
| 2 | Setup: $V_{\mathrm{std},n}$, $E_\gamma = \gamma(L) - m_n J$, $\widehat{E}_{\gamma,L} = P_n^\top E_{\gamma,L} P_n$, root lattice $A_{n-1}$ | (definitions) | — | `paper/main.tex`, `paper/notation.tex` | — | — |
| 3 | Box Band Lemma: $B_{\mathrm{rb}}^\top E B_{\mathrm{cs}} = 0$ on square Sudoku and rectangular $h\times w$ Sudoku; gerechte designs satisfy the corresponding cell-region sum identity; standard boxes are nonzero on a generic cyclic Latin square | Theorem A + Cor. rectangular + Lemma gerechte | THEOREM | `paper/main.tex`, [`box_band_lemma_witness.json`](box_band_lemma_witness.json) | `scripts/legacy/phase_9_13.py`, `scripts/legacy/task_9_10c_box_subspace.py`, `scripts/legacy/phase_9_14.py`, `certified/build_certified.py` | [`box_band_lemma_witness.json`](box_band_lemma_witness.json) |
| 4 | Exhaustive $S_9$ scan of base1: rank distribution $\{8\!:\!362858,\ 7\!:\!22\}$; the 22 odd-rank relabelings are explicitly enumerated and four-way certified (rank, kernel, integer SNF, $\mathbb{F}_2$ SNF) | Proposition P1 | EXHAUSTIVE | `data/phase_9_13_results.json`, `data/phase_9_13f_certificate.json`, `certified/four_way_base1_22.json` | `scripts/legacy/phase13b_exhaustive.py`, `scripts/legacy/phase_9_13f_certificate.py`, `certified/recover_base1_22_perms.py` | [`base1_odd_rank_22.json`](base1_odd_rank_22.json) |
| 5.1 | $n=6$ Roku-Doku: 6.617% full-rank-over-$\mathbb{F}_2$ for the auxiliary matrix $A_{ij}=L_{ij}-L_{i,6}$ in Sudoku-6 vs 6.122% in LS-6, descriptive two-proportion $z = 4.0137$, two-sided $\,\mathrm{erfc}(|z|/\sqrt2)=5.9773\cdot10^{-5}$; this is the determinant-paper rank-parity proxy, not $\rank(P^TEP)$ and not the integer odd-rank locus | Proposition P2 | EXHAUSTIVE | `data/phase_9_14_results.json`, [`n6_odd_rank_census.json`](n6_odd_rank_census.json) | `scripts/legacy/phase_9_14.py` | [`n6_odd_rank_census.json`](n6_odd_rank_census.json) |
| 5.2 | $n=7$, certificate $V_1$: stabilizer $H_{V_1} \cong D_6$ (order 12), embedded inside the signed-pair group $B_3=C_2\wr S_3< S_7$ preserving pairs $\{1,7\},\{2,6\},\{3,5\}$ and fixing $4$; dihedral relation $r^6 = s^2 = (rs)^2 = e$ verified | Certified Proposition B | CERTIFIED_EXHAUSTIVE | `data/phase_9_15_results.json`, `data/phase_9_15A_H_V1_identification.json` | `scripts/legacy/phase_9_15.py`, `scripts/legacy/phase_9_15A_identify_HV1.py`, `scripts/legacy/check_b3_subgroup.py` | [`n7_HV1_D6.json`](n7_HV1_D6.json) |
| 5.3 | Rank parity is *not* used (cite-only bridge to sibling) | (cite) | — | [Determinant Divisibility of Centered Latin Squares](https://github.com/aconsciousfractal/Determinant-Divisibility-of-Centered-Latin-Squares) | — | — |
| 6 | Lift Lemma: $\forall \gamma \in S_9,\ \forall w \in V_{\mathrm{std}}$, $E_\gamma w = 0 \iff \gamma(L) w = 0$ (and dually); kernel of every odd-rank relabeling lifts to $V_{\mathrm{std}}$. Verified by full Sympy nullspace on **all 22** base1 + **all 100** $M_{19}$ rank-7 relabelings + the **72** shared-left-kernel subset (kernel dim 2 in $\mathbb{Q}^9$, dim 1 in $V_{\mathrm{std}}$) | Theorem C | THEOREM | `paper/main.tex`, `lift_lemma_evidence.json`, recovered perm lists in `certified/` | `certified/build_certified.py`, `certified/recover_M19_100_perms.py`, `certified/recover_base1_22_perms.py` | [`lift_lemma_evidence.json`](lift_lemma_evidence.json) |
| 7 | $A_8$ saturation: $\Lambda_M = \Lambda_{O_6} = V_{\mathrm{std},9} \cap \mathbb{Z}^9 = A_8$; rank 8, all elementary divisors equal 1 | Certified Proposition D | CERTIFIED_EXHAUSTIVE | `paper/main.tex`, `A8_saturation.json` | `certified/build_certified.py` | [`A8_saturation.json`](A8_saturation.json) |
| 8 | $M_{19}$ exhaustive audit: 100/362880 rank-7 relabelings; $\Gamma^*_{M_{19}}$ is the 72-element shared-left-kernel / middle-band subset and a bilateral coset; $K_L$ and $K_R$ both have order 72; $K_L \cap K_R = 1$; $K_L \cong K_R \cong 3^2{:}D_4$ with normal regular $3^2$ translation subgroup, $D_4$ point stabilizer, exponent 12, derived series $72\!\to\!18\!\to\!9\!\to\!1$, abelianization $(\mathbb{Z}/2)^2$, trivial centre, 9 conjugacy classes; neither group is contained in $\mathrm{Stab}_{S_9}(w^L_9)$ | Certified Proposition E + Proposition P3 | CERTIFIED_EXHAUSTIVE (affine identification) + EXHAUSTIVE ($S_9$ modular determinant scan + exact rank) | `M19_audit.json`, `M19_KL_KR_affine_D4.json`, `data/phase_9_18_S9x_results.json`, `data/phase_9_18_S9x_bis_results.json` | `scripts/legacy/phase_9_18_S9x.py`, `scripts/legacy/phase_9_18_S9x_bis.py`, `certified/recover_M19_100_perms.py`, `certified/build_certified.py` | [`M19_audit.json`](M19_audit.json), [`M19_KL_KR_affine_D4.json`](M19_KL_KR_affine_D4.json) |
| 8 hypergraph | Middle-band balance criterion and reduced-hypergraph realization at $M_{19}$: $\Gamma^*_{M_{19}}$ is exactly the set of $72$ relabelings $\gamma\in S_9$ for which $\sum_{i=4}^{6}\gamma(M19_{ij})=15$ for every column $j$; equivalently $w^T E_{\gamma,M_{19}}=0$ for $w=(1,1,1,-2,-2,-2,1,1,1)$. The same set is the isomorphism torsor from the reduced M19 middle-band hypergraph to a six-edge magic-triple subhypergraph; $\operatorname{Aut}(\mathcal H_{M19}^{mid})=K_R$, $\operatorname{Aut}(\mathcal T^-_{15})=K_L$, while the full magic-triple hypergraph has automorphism group of order $8$. | Proposition P5 + Certified Proposition hypergraph realization | CERTIFIED_EXHAUSTIVE | direct lemma/certified proposition in paper §8 + `data/phase_9_18_S9x_results.json` | `scripts/legacy/phase_9_18_S9x.py`, `certified/build_certified.py` | [`M19_middle_band_balance.json`](M19_middle_band_balance.json) |
| 9 (i) | Finite dichotomy at $n \in \{7, 9, 11\}$: certified identification $\mathrm{Stab}(M_n) \cong D_4 \times S_{n-4}$, $\mathrm{Stab}(O_n) \cong C_2 \times S_{n-4}$; $\Lambda_{M_n} = \Lambda_{O_n} = A_{n-1}$; orbit-distance $d_{M_n}(\mathrm{orb}(O_n)) = 2$ | Certified Proposition F | CERTIFIED_EXHAUSTIVE | `data/phase_9_18_S9q_gamma1.json`, [`dichotomy_n7_n9_n11.json`](dichotomy_n7_n9_n11.json) | `scripts/legacy/phase_9_18_S9q_gamma1.py` | [`dichotomy_n7_n9_n11.json`](dichotomy_n7_n9_n11.json) |
| 9 (ii) | Stabilizer / orbit-size formulas exhaustively verified for **every** $n$ with $5 \le n \le 13$ | Proposition P4 | EXHAUSTIVE | `data/phase_9_18_S9q_gamma1.json` | `scripts/legacy/phase_9_18_S9q_gamma1.py` | [`stabilizers_5_to_13.json`](stabilizers_5_to_13.json) |
| 9 (iii) | Universal stabilizer and lattice-saturation dichotomy for **every** $n \ge 5$: $\mathrm{Stab}(M_n) \cong D_4\times S_{n-4}$, $\mathrm{Stab}(O_n) \cong C_2\times S_{n-4}$, and $\Lambda_{M_n}=\Lambda_{O_n}=A_{n-1}$ | Theorem H | THEOREM | direct argument in paper §9.3 | (none) | (proof in paper; finite certificates remain independent audits) |
| 10.3 | $\mathrm{Aut}(L_{M_{19}}) = 1$, hence $K_L \neq \mathrm{image}(\mathrm{Aut}(L_{M_{19}}))$; verified by exhaustive scan over the 9! triples $(\alpha, \beta, \gamma)$ | Certified Proposition G | CERTIFIED_EXHAUSTIVE | `data/phase_9_18_S9q_C1.json`, [`Aut_LM19_trivial.json`](Aut_LM19_trivial.json) | `scripts/legacy/phase_9_18_S9q_C1.py` | [`Aut_LM19_trivial.json`](Aut_LM19_trivial.json) |
| 11 | Methods and tier discipline | (cross-cutting) | — | `README.md`, `REPRODUCE.md`, this ledger | `certified/build_certified.py`, `certified/verify_all.py`, `certified/reproduce_all.py` | all 12 certified JSONs |
| 12 | Open problems | Larger-order affine realization for odd $m \ge 5$; all-$n$ orbit-distance beyond the verified range; finite-corpus prevalence questions | CONJECTURE | Discussion section of `paper/main.tex` | — | — |
| App. A | Explicit list of the 22 base1 odd-rank relabelings | (data) | EXHAUSTIVE | `data/phase_9_13_results.json` | `scripts/legacy/phase_9_13.py`, `scripts/legacy/phase_9_13f_certificate.py` | [`base1_odd_rank_22.json`](base1_odd_rank_22.json) |
| App. B | $K_L = 3^2{:}D_4$ affine data and generators (8 explicit permutations of $\{1,\dots,9\}$) | (data) | CERTIFIED_EXHAUSTIVE | `data/phase_9_18_S9x_bis_results.json`, `M19_KL_KR_affine_D4.json` | `scripts/legacy/phase_9_18_S9x_bis.py`, `certified/build_certified.py` | [`M19_KL_KR_affine_D4.json`](M19_KL_KR_affine_D4.json) |
| App. C | Stabilizer / orbit table for $n \in \{5,\dots,13\}$ | (data) | EXHAUSTIVE | `data/phase_9_18_S9q_gamma1.json` | `scripts/legacy/phase_9_18_S9q_gamma1.py` | [`stabilizers_5_to_13.json`](stabilizers_5_to_13.json) |

---

## B. Publicly bounded or excluded statements

These rows preserve the public claim boundary without depending on
unshipped research reports.

| Statement | Public disposition | Reason |
|---|---|---|
| $Q_{239}$ Bonferroni-uniqueness across "the universe of Sudoku bases" | **Finite-corpus proposition, $N=400$** | No proof of universality; the paper cites only the finite corpus count. |
| $f_7$ distribution and general rank-7 frequency bounds | **Finite-corpus proposition** | Measured on the stated $N=400$ corpus only. |
| "Typical" $|\Gamma_B|$ statistics and Sudoku-9 prevalence rates | **Finite-corpus proposition** | Sample-driven, not exhaustive over all bases or $S_9\times S_9$ pairs. |
| Trace-cyclic barrier $\lim M/T\le 1$ and $M/T(n=4)=1.094$ | **Excluded** | Protocol-dependent artifacts; retained only as cautions in paper §11. |
| Universal orbit-distance $d_{M_n}(\mathrm{orb}(O_n))=2$ for every $n\ge5$ | **Open problem** | Exhaustively verified only for $5\le n\le13$; the all-$n$ theorem concerns stabilizers and lattice saturation, not distance. |
| "$3^2{:}Q_8$ appears for every odd $m\ge3$" | **Excluded** | The certified $M_{19}$ group is $3^2{:}D_4$ and no universal Hessian realization is claimed. |
| CWF, writhe, arc length, 600-cell, ZKP, or cryptographic applications | **Out of scope** | Not part of this paper or its 12-certificate package. |

---

## C. Public exclusion boundary

The paper and repository do not claim:

- ❌ "A new theory of Sudoku."
- ❌ "600-cell underlies Sudoku."
- ❌ "Immediate cryptographic applications."
- ❌ "Hessian classifies all Sudoku."
- ❌ "$K_L$ or $K_R$ is contained in $\mathrm{Stab}_{S_9}(w^L_9)$."
- ❌ "47 theorems established."
- ❌ "$\varphi$ is fundamental to Sudoku."
- ❌ Universal orbit-distance $d_{M_n}(\mathrm{orb}(O_n))=2$ outside the verified range $5 \le n \le 13$.

The whitelist is exactly the rows of §A above.  Anything else needs an
explicit ledger amendment before it can enter the manuscript.

---

## D. Verification

This ledger is consistent with:

- the 12 entries of [`MANIFEST.sha256`](MANIFEST.sha256) (one row of §A per certified JSON, plus
  cite-only / definition-only rows);
- [`verify_all.py`](verify_all.py) (claim-level quick checks against each JSON);
- [`reproduce_all.py`](reproduce_all.py) (byte-equal regeneration of all 12 JSONs).

Last refreshed: 2026-08-03 for release `v1.0.1`.
