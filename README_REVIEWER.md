# Reviewer Guide

Companion documentation added 2026-07-08, after publication; it does not
modify the paper.

## Ten-Minute Path

1. Read the abstract in `paper/main.tex` — note the explicit split between
   "the analytic part" (proved lemmas) and "the certified finite part"
   (exhaustive selected orders).
2. Read `certified/README.md` — one canonical JSON per main result, each
   with an explicit tier (THEOREM / CERTIFIED_EXHAUSTIVE / EXHAUSTIVE /
   EMPIRICAL / OBSERVED TEMPLATE / CONJECTURE).
3. Run the verifier:
   `PYTHONUTF8=1 python certified/verify_all.py`
   (expects "ALL CHECKS PASSED (12 certificates, 18 scripts, 13 sources)";
   verified 2026-07-08, ~1 min).

## Thirty-Minute Path

4. Spot-check one certificate chain: `certified/M19_audit.json` against
   `prop:M19-audit` in the paper and `certified/MANIFEST.sha256`.
5. Run the byte-identical re-derivation for one certificate via
   `certified/reproduce_all.py` (full run is longer; per-certificate
   builders are in `certified/build_certified.py`).
6. Note the verifier's paper-grep step: the repo mechanically checks the
   paper for forbidden overclaim phrases.
7. Read `PUBLIC_CLAIM_BOUNDARY.md` for what may and may not be quoted.

## Main Claims

- Analytic (all stated n): Box Band Lemma, Lift Lemma into A_{n−1},
  general-n {M_n, O_n} stabilizer dichotomy (ledger A1-A3).
- Certified finite (stated orders/bases only): odd-rank censuses, M19
  audit, 3²:D4 affine coset symmetry at M19, A8 saturation, trivial
  Aut(L_M19) (F1-F5).

## Known Limits

- Finite results are exhaustive only at their stated orders and bases;
  no cross-order generalization is claimed.
- Tier fields in the certificate JSONs are the claim-strength authority;
  EMPIRICAL/OBSERVED TEMPLATE/CONJECTURE rows are not theorems.
- `reproduce_all.py` was not fully re-run in the 2026-07-08 retrofit;
  `verify_all.py` (manifest + semantic checks) was run and passed.
