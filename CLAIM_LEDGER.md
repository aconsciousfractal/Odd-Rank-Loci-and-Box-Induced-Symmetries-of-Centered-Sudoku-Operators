# Claim Ledger

Companion documentation added 2026-07-08, after publication (no git tags;
release commit e4c827c "Finalize Sudoku paper release package"); it does
not modify the paper. This repository already carries a per-result
certificate system — `certified/README.md` maps each of the 12 canonical
JSONs to its paper result and tier
(`THEOREM | CERTIFIED_EXHAUSTIVE | EXHAUSTIVE | EMPIRICAL | OBSERVED
TEMPLATE | CONJECTURE`), `certified/MANIFEST.sha256` pins them, and
`certified/verify_all.py` re-checks schema, manifest, semantics, and even
greps the paper for forbidden overclaim phrases. **That table is the
authoritative per-result ledger; this file is a thin index over it.**

## Analytic Results (proved in paper, all stated n)

| ID | Statement (scoped) | Source locator | Certificate |
| --- | --- | --- | --- |
| A1 | Box Band Lemma: for Sudoku grids, averaging the centered operator E_γ,L over row-bands and column-stacks gives zero for EVERY relabeling γ; a cyclic Latin-square witness shows the identity is not forced by Latin-squareness alone. Rectangular and gerechte variants as stated. | `thm:box-band`; `cor:rectangular-box-band`; `lem:gerechte-region-sum` | `certified/box_band_lemma_witness.json` (THEOREM tier + counterexample witness) |
| A2 | Lift Lemma: kernel lines on the standard hyperplane admit primitive integral representatives in the root lattice A_{n−1}; for w in V_std, E_γ w = 0 iff γ(L)w = 0 (and dually). | `thm:lift`; `cor:cert-in-A` | `certified/lift_lemma_evidence.json` (THEOREM tier) |
| A3 | General-n stabilizer/lattice dichotomy for the two sparse root-lattice templates {M_n, O_n}, every n ≥ 5, with stabilizers D4 × S_{n−4} vs C2 × S_{n−4}. | `thm:universal-dichotomy` | `certified/dichotomy_n7_n9_n11.json` (certified instances), `certified/stabilizers_5_to_13.json` (5 ≤ n ≤ 13 verification) |

## Certified Finite Results (exhaustive at the stated orders ONLY)

The full list is `certified/README.md`. Highlights:

| ID | Statement (scoped) | Certificate |
| --- | --- | --- |
| F1 | base1 odd-rank census: 362858 rank-8 / 22 rank-7 (exhaustive over S9 relabelings at the stated base). | `base1_odd_rank_22.json` |
| F2 | M19 audit: exhaustive 100/362880 rank-7 at the laboratory base M19. | `M19_audit.json` (`prop:M19-audit`) |
| F3 | Affine coset symmetry K_L ≅ K_R ≅ 3²:D4 AT THE LABORATORY BASE M19; \|Aut(L_M19)\| = 1. | `M19_KL_KR_affine_D4.json`, `Aut_LM19_trivial.json` |
| F4 | Γ*_M19 = exactly the 72 shared-left-kernel / middle-band balanced relabelings = the 72 reduced-hypergraph embeddings into the six-edge magic subhypergraph. | `M19_middle_band_balance.json` |
| F5 | Λ_M = Λ_{O6} = A8 saturation; n=6 auxiliary full-rank F2 census (Sudoku 6.62% vs LS 6.12%); H_{V1} ≅ D6 inside C2 wr S3. | `A8_saturation.json`, `n6_odd_rank_census.json`, `n7_HV1_D6.json` |

## Downstream link — FCIG grammar X05 (added 2026-07-17)

The order-72 group of this repository is the affine symmetry group **3²:D₄ = C₃²⋊D₄ < AGL(2,3)** of the
natural 3×3 (Lo Shu) arrangement of the values 1..9. In the FCIG grammar paper it is the **X05**
exhibit's `G₀`: the census-mined `G₀` is a fixed reference copy, and the certified left/right
stabilizers of **F3** (`K_L`, `K_R`, with `K_L = π₀ K_R π₀⁻¹`, `K_L ≠ K_R`) are its two conjugate
realizations at the base M19. **Claim-level note:** F3 (single base, `K_L ≠ K_R`) is CL3-certified here;
the grammar paper's *cross-grid "single universal `G₀`"* statement is CL1 (computational evidence) until
a certified set-equality over the full census is minted. [LEDGER A-3 / E-2, resolved IDENTITY 2026-07-17.]

## Verification Trail

- `certified/verify_all.py` — ALL CHECKS PASSED (12 certificates, 18
  scripts, 13 sources) on 2026-07-08.
- `certified/reproduce_all.py` — byte-identical re-derivation path (not
  re-run 2026-07-08; verify_all's manifest check covers integrity).
- The verifier's paper-grep step enforces the absence of forbidden
  overclaim phrases in `paper/main.tex` — an in-repo RT-10 check.
