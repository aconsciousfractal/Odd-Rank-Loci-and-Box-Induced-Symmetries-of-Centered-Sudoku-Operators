# TIER_LEDGER.md — Paper II claim-by-claim tier discipline

> **Purpose.** Single canonical ledger that maps every claim in the
> paper to (i) its primary source document/script, (ii) the certified
> JSON in `paper/certified/`, and (iii) the **tier** at which the paper
> exposes it. Built per task **P0.4** of [`PRE_DRAFT_TASKS.md`](../PRE_DRAFT_TASKS.md);
> mirrors the §-to-source map of [`PAPER_PLAN.md §3`](../PAPER_PLAN.md).

> **Tier alphabet (5 levels).**  THEOREM | EXHAUSTIVE | EMPIRICAL |
> OBSERVED TEMPLATE | CONJECTURE.  Composite: THEOREM (EXHAUSTIVE)
> = analytical theorem whose hypothesis includes a finite exhaustive
> check.  Anything in the source corpus that lacks an explicit tier
> tag, or is stated at a higher tier than its evidence supports, is
> **demoted** in the row below and a justification line is given in §B.

> **Filter applied.** Every numerical or structural claim in the paper
> must come from one of the rows in §A. Claims that appear in the
> source reports but are **not** in §A are either (a) out of perimeter
> (see [`PAPER_PLAN.md §7`](../PAPER_PLAN.md)), or (b) demoted (see §B).

---

## A. Section-to-source-to-certificate ledger

| § | Claim (paper-level statement) | Paper label | Tier in paper | Primary source(s) | Producer script(s) | Certified JSON |
|---|---|---|---|---|---|---|
| 2 | Setup: $V_{\mathrm{std},n}$, $E_\gamma = \gamma(L) - m_n J$, $\widehat{E}_{\gamma,L} = P_n^\top E_{\gamma,L} P_n$, root lattice $A_{n-1}$ | (definitions) | — | [`MATHEMATICAL_FOUNDATIONS.md`](../../../docs/MATHEMATICAL_FOUNDATIONS.md), [`INVARIANT_MAP.md`](../../../docs/INVARIANT_MAP.md), [`paper/notation.tex`](../paper/notation.tex) | `e_transform.py`, `psi_mapping.py` | (none — definitions) |
| 3 | Box Band Lemma: $B_{\mathrm{rb}}^\top E B_{\mathrm{cs}} = 0$ on every Sudoku $L$; nonzero on a generic Latin square | Theorem A | THEOREM | [`PHASE_9_13_REPORT.md`](../../../docs/PHASE_9_13_REPORT.md) §9.13e | `phase_9_13.py`, `task_9_10c_box_subspace.py` | [`box_band_lemma_witness.json`](box_band_lemma_witness.json) |
| 4 | Exhaustive $S_9$ scan of base1: rank distribution $\{8\!:\!362858,\ 7\!:\!22\}$; the 22 odd-rank relabelings are explicitly enumerated and four-way certified (rank, kernel, integer SNF, $\mathbb{F}_2$ SNF) | Proposition P1 | EXHAUSTIVE | [`PHASE_9_13_REPORT.md`](../../../docs/PHASE_9_13_REPORT.md) §9.13c, §9.13f | `phase13b_exhaustive.py`, `phase_9_13f_certificate.py`, `paper/certified/recover_base1_22_perms.py` | [`base1_odd_rank_22.json`](base1_odd_rank_22.json) |
| 5.1 | $n=6$ Roku-Doku: 6.617% odd-rank in Sudoku-6 vs 6.122% in LS-6, two-proportion $z = 4.0137$, $p < 5{\cdot}10^{-5}$ | Proposition P2 | EXHAUSTIVE | [`PHASE_9_14_REPORT.md`](../../../docs/PHASE_9_14_REPORT.md) | `phase_9_14.py` | [`n6_odd_rank_census.json`](n6_odd_rank_census.json) |
| 5.2 | $n=7$, certificate $V_1$: stabilizer $H_{V_1} \cong D_6$ (order 12), embedded inside $B_3 < S_7$; dihedral relation $r^6 = s^2 = (rs)^2 = e$ verified | Theorem B | THEOREM | [`PHASE_9_15_REPORT.md`](../../../docs/PHASE_9_15_REPORT.md), [`PHASE_9_15_AUDIT_FOLLOWUP_REPORT.md`](../../../docs/PHASE_9_15_AUDIT_FOLLOWUP_REPORT.md) | `phase_9_15.py`, `phase_9_15A_identify_HV1.py`, `check_b3_subgroup.py` | [`n7_HV1_D6.json`](n7_HV1_D6.json) |
| 5.3 | Rank parity is *not* used (cite-only bridge to sibling) | (cite) | — | sibling [`paper/parity conjecture/draft_v1.tex`](../../parity%20conjecture/draft_v1.tex) | (none) | (none) |
| 6 | Lift Lemma: $\forall \gamma \in S_9,\ \forall w \in V_{\mathrm{std}}$, $E_\gamma w = 0 \iff \gamma(L) w = 0$ (and dually); kernel of every odd-rank relabeling lifts to $V_{\mathrm{std}}$. Verified by full Sympy nullspace on **all 22** base1 + **all 72** $M_{19}$ relabelings (kernel dim 2 in $\mathbb{Q}^9$, dim 1 in $V_{\mathrm{std}}$) | Theorem C | THEOREM | [`PHASE_9_13_S9_REPORT.md`](../../../docs/PHASE_9_13_S9_REPORT.md) §S9–§S9c, T8.9.1 | `phase_9_18_S9.py`, `phase_9_18_S9c.py`, `phase_9_18_S9d.py`–`phase_9_18_S9i.py` | [`lift_lemma_evidence.json`](lift_lemma_evidence.json) |
| 7 | $A_8$ saturation: $\Lambda_M = \Lambda_{O_6} = V_{\mathrm{std},9} \cap \mathbb{Z}^9 = A_8$; rank 8, all elementary divisors equal 1 | Theorem D | THEOREM | [`PHASE_9_13_S9_REPORT.md`](../../../docs/PHASE_9_13_S9_REPORT.md) §S9r–§S9y, §8.19 (T8.19.1), T8.25.terdecies.1 | `phase_9_18_S9r.py`, `phase_9_18_S9w.py`, `phase_9_18_S9w_S9y_closure_audit.py`, `phase_9_18_S9y_six.py` | [`A8_saturation.json`](A8_saturation.json) |
| 8 | $M_{19}$ exhaustive audit: 100/362880 rank-7 relabelings; $\Gamma^*_{M_{19}}$ is a bilateral coset; $K_L$ and $K_R$ both of order 72; $K_L \cap K_R = 1$; $K_L \cong K_R \cong 3^2{:}Q_8 = \mathrm{SmallGroup}(72,41)$ (exponent 6, derived series $72\!\to\!18\!\to\!9\!\to\!1$, abelianization $\mathbb{Z}/4$, trivial centre, 9 conjugacy classes) | Theorem E + Proposition P3 | THEOREM (Hessian identification) + EXHAUSTIVE ($S_9$ scan) | [`PHASE_9_13_S9_REPORT.md`](../../../docs/PHASE_9_13_S9_REPORT.md) §S9x, §S9x-bis (T8.25.quinquies.2, T8.25.sexies.4), §S9x-ter | `phase_9_18_S9x.py`, `phase_9_18_S9x_bis.py`, `phase_9_18_S9x_ter.py` | [`M19_audit.json`](M19_audit.json), [`M19_KL_KR_3sq_Q8.json`](M19_KL_KR_3sq_Q8.json) |
| 9 (i) | Universal dichotomy at $n \in \{7, 9, 11\}$: analytic identification $\mathrm{Stab}(M_n) \cong D_4 \times S_{n-4}$, $\mathrm{Stab}(O_n) \cong C_2 \times S_{n-4}$; $\Lambda_{M_n} = \Lambda_{O_n} = A_{n-1}$; orbit-distance $d_{M_n}(\mathrm{orb}(O_n)) = 2$ | Theorem F | THEOREM | [`PHASE_9_13_S9_REPORT.md`](../../../docs/PHASE_9_13_S9_REPORT.md) T8.25.duodecies/terdecies/quaterdecies/sexdecies | `phase_9_18_S9q_gamma1.py` | [`dichotomy_n7_n9_n11.json`](dichotomy_n7_n9_n11.json) |
| 9 (ii) | Stabilizer / orbit-size formulas exhaustively verified for **every** $n$ with $5 \le n \le 13$ | Proposition P4 | EXHAUSTIVE | [`PHASE_9_13_S9_REPORT.md`](../../../docs/PHASE_9_13_S9_REPORT.md) §8.27 (γ-1) | `phase_9_18_S9q_gamma1.py` | [`stabilizers_5_to_13.json`](stabilizers_5_to_13.json) |
| 9 (iii) | Universal $\{M_n, O_n\}$ dichotomy for **every** $n \ge 5$ | Conjecture C10 | CONJECTURE | (no source — open) | (none) | (none) |
| 10 | $\mathrm{Aut}(L_{M_{19}}) = 1$, hence $K_L \neq \mathrm{image}(\mathrm{Aut}(L_{M_{19}}))$; verified by exhaustive scan over the 9! triples $(\alpha, \beta, \gamma)$ | Theorem G | THEOREM (EXHAUSTIVE) | [`PHASE_9_13_S9_REPORT.md`](../../../docs/PHASE_9_13_S9_REPORT.md) §8.33 (T8.33.3, Cor. T8.33.4) | `phase_9_18_S9q_C1.py` | [`Aut_LM19_trivial.json`](Aut_LM19_trivial.json) |
| 11 | Methods, tier-discipline protocol | (cross-cutting) | — | [`COUNTING_PROTOCOLS.md`](../../../docs/COUNTING_PROTOCOLS.md), this ledger | (all of the above) | (all of the above) |
| 12 | Open problems | Conjecture C9, C10, P5-realization for odd $m \ge 5$ | CONJECTURE | [`PROJECT_STATE_2026-04-26.md`](../../../docs/PROJECT_STATE_2026-04-26.md) §6 | — | — |
| App. A | Explicit list of the 22 base1 odd-rank relabelings | (data) | EXHAUSTIVE | `data/phase_9_13_results.json` | `phase_9_13.py`, `phase_9_13f_certificate.py` | [`base1_odd_rank_22.json`](base1_odd_rank_22.json) |
| App. B | $K_L = 3^2{:}Q_8$ generators (8 explicit permutations of $\{1,\dots,9\}$) | (data) | THEOREM | `data/phase_9_18_S9x_bis.json` | `phase_9_18_S9x_bis.py` | [`M19_KL_KR_3sq_Q8.json`](M19_KL_KR_3sq_Q8.json) |
| App. C | Stabilizer / orbit table for $n \in \{5,\dots,13\}$ | (data) | EXHAUSTIVE | `data/phase_9_18_S9q_gamma1.json` | `phase_9_18_S9q_gamma1.py` | [`stabilizers_5_to_13.json`](stabilizers_5_to_13.json) |

---

## B. Demoted claims — present in the source corpus, **not** stated as theorem in the paper

These are claims the source reports occasionally cite at theorem-level
(or simply stated without an explicit tier), but which the paper exposes
only as proposition over a finite corpus, or excludes outright.

| Claim (as in source) | Source location | Original tier in source | Tier in paper | Reason for demotion |
|---|---|---|---|---|
| $Q_{239}$ Bonferroni-uniqueness across "the universe of Sudoku bases" | [`PHASE_9_13_S9_REPORT.md`](../../../docs/PHASE_9_13_S9_REPORT.md) §8.* (corpus statistics block) | Theorem-style assertion | **Proposition (finite corpus, $N=400$)** | The uniqueness has been verified only on the $N{=}400$ `sample_sudoku` corpus; no proof of universality. The paper cites the count only as a corpus statistic, not as a structural theorem. |
| $f_7$ distribution / "the rank-7 relabeling count is bounded by …" generality statements | [`PHASE_9_13_S9_REPORT.md`](../../../docs/PHASE_9_13_S9_REPORT.md) corpus statistics block | Theorem-style | **Proposition (finite corpus)** | Distribution measured on $N{=}400$ corpus; paper restricts statement to "in the $N{=}400$ corpus". |
| $\Gamma_B$ size statistics ("typical $|\Gamma_B|$") | [`PHASE_9_13_S9_REPORT.md`](../../../docs/PHASE_9_13_S9_REPORT.md) corpus block | Stated without explicit tier | **Proposition (finite corpus)** | Sample-driven; not exhaustive over $S_9 \times S_9$ pairs. |
| Corpus prevalence rates ("X % of Sudoku-9 grids exhibit feature Y") | various phase 9.x reports | Often presented as fact | **Proposition (finite corpus)** | Rates are over `sample_sudoku` ($N{=}400$); not population statistics. |
| Trace-cyclic barrier $\lim M/T \le 1$ | older trace–cokernel notes | Theorem-style | **Excluded** | Protocol-dependent; explicitly retracted. Mentioned in §11 only as a cautionary example. |
| $M/T(n=4) = 1.094$ | older protocol log | Reported value | **Excluded** | Retracted as artifact; cited in §11 only. |
| Universal $\{M_n, O_n\}$ dichotomy for all $n \ge 5$ | implicit in narrative phrasing across §S9–§8.27 | Sometimes asserted | **Conjecture C10** | Analytic at $n \in \{7, 9, 11\}$; exhaustive at $5 \le n \le 13$; nothing more is proved. |
| "$3^2{:}Q_8$ Hessian appears for every odd $m \ge 3$" | abstract statement of Proposition P5 | Algebraic theorem ($\Leftrightarrow$ on $m$ odd) + empirical Sudoku realization at $n=9$ | **Theorem on the algebraic side, Conjecture C9 on Sudoku realization for $m \ge 5$** | The "$\Leftrightarrow m$ odd" is proven at the abstract group level; realization inside Sudoku of order $n = m^2$ is verified only for $m = 3$. |
| CWF / writhe / arc_length / 600-cell / ZKP / cryptographic-application claims | various older docs | Mixed | **Out of perimeter** | Excluded by [`PAPER_PLAN.md §7`](../PAPER_PLAN.md). |

---

## C. Anti-claims (recap, locked)

The paper's drafting guard-rail is:

- ❌ "A new theory of Sudoku."
- ❌ "600-cell underlies Sudoku."
- ❌ "Immediate cryptographic applications."
- ❌ "Hessian classifies all Sudoku."
- ❌ "47 theorems established."
- ❌ "$\varphi$ is fundamental to Sudoku."
- ❌ Universal $\{M_n, O_n\}$ for all $n$ outside the **Conjecture** statement.

The whitelist is exactly the rows of §A above.  Anything else needs an
explicit ledger amendment before it can enter the manuscript.

---

## D. Verification

This ledger is consistent with:

- the 11 entries of [`MANIFEST.sha256`](MANIFEST.sha256) (one row of §A per certified JSON, plus
  cite-only / definition-only rows);
- [`verify_all.py`](verify_all.py) (claim-level quick checks against each JSON);
- [`reproduce_all.py`](reproduce_all.py) (byte-equal regeneration of all 11 JSONs).

Last refreshed: 2026-04-26 (P0.4 closure).
