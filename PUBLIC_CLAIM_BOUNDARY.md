# Public Claim Boundary

Companion documentation added 2026-07-08, after publication; it does not
modify the paper.

## Can Say

- The Box Band Lemma, Lift Lemma, and general-n stabilizer/lattice
  dichotomy (n ≥ 5) are proved analytically (ledger A1-A3).
- The finite censuses and symmetry identifications are exhaustive/certified
  at their stated orders and bases, with pinned SHA-256 certificates and a
  standalone verifier (F1-F5).
- The affine coset symmetry 3²:D4 and trivial autotopism group hold AT THE
  LABORATORY BASE M19 (certified).

## Must Not Say

- That finite-order censuses (n = 6, 7, 9, order-9 bases) generalize to
  other orders or bases — each is exhaustive only at its stated scope.
- That the M19 coset symmetry is generic for Sudoku bases — M19 is the
  paper's laboratory base; the certificate is base-specific.
- That any row of tier EMPIRICAL, OBSERVED TEMPLATE, or CONJECTURE in the
  certificate schema is a theorem — the tier field is the authority.
- Anything about relabeling-family structure beyond what the paper and the
  12 certificates state.

## Scope Notes

- The centered operator convention is E_γ,L = γ(L) − ((n+1)/2)J_n
  restricted to the standard hyperplane, per the paper.
- `certified/verify_all.py` includes a forbidden-phrase grep over the
  paper — run it after any wording change.
