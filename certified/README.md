# Certified results package

This directory contains one canonical JSON per main result in the
paper, plus a manifest and a verifier.

## Files

| File | Paper result | Tier |
|---|---|---|
| `base1_odd_rank_22.json` | Prop P1: 362858 rank-8 / 22 rank-7 on base1 | EXHAUSTIVE |
| `n6_odd_rank_census.json` | Prop P2: $n{=}6$ Sudoku 6.62% vs LS 6.12% | EXHAUSTIVE |
| `n7_HV1_D6.json` | Theorem B: $H_{V_1}\cong D_6$ inside $B_3$ | THEOREM |
| `box_band_lemma_witness.json` | Theorem A: $P_{\mathrm{rb}} E P_{\mathrm{cs}}^\top = 0$ | THEOREM (+ counterexample for generic LS) |
| `lift_lemma_evidence.json` | Theorem C: for $w\in V_{\mathrm{std}}$, $E_\gamma w=0 \iff \gamma(L)w=0$ (and dually) | THEOREM |
| `M19_audit.json` | Prop P3: $M_{19}$ exhaustive 100/362880 rank-7 | EXHAUSTIVE |
| `M19_KL_KR_3sq_Q8.json` | Theorem E: $K_L\cong K_R\cong 3^2{:}Q_8$ | THEOREM |
| `A8_saturation.json` | Theorem D: $\Lambda_M=\Lambda_{O_6}=A_8$ | THEOREM |
| `dichotomy_n7_n9_n11.json` | Theorem F: stabilizer dichotomy at $n\in\{7,9,11\}$ | THEOREM |
| `stabilizers_5_to_13.json` | Prop P4: stabilizer/orbit formulas verified for $5\le n\le 13$ | EXHAUSTIVE |
| `Aut_LM19_trivial.json` | Theorem G: $|\mathrm{Aut}(L_{M_{19}})|=1$ | THEOREM (EXHAUSTIVE) |

## Schema

Each JSON conforms to:

```json
{
  "result_id": "<stable identifier>",
  "tier": "THEOREM | THEOREM_EXHAUSTIVE | EXHAUSTIVE | EMPIRICAL | OBSERVED TEMPLATE | CONJECTURE",
  "paper_label": "Theorem A | Proposition P1 | ...",
  "claim": "<one or two sentences>",
  "inputs": { ... },
  "outputs": { ... },
  "scripts": [
    { "path": "scripts/legacy/<file>.py", "sha256": "<hex64>" },
    ...
  ],
  "source": [
    { "file": "data/<file>.json | certified/<file>.json", "sha256": "<hex64>",
      "extracted_keys": [ ... ] },
    ...
  ],
  "seed": null,
  "produced_utc": "YYYY-MM-DDTHH:MM:SSZ",
  "schema_version": 1
}
```

`scripts[]` is a list (one or more producer scripts under
`scripts/legacy/` or canonical scripts under `certified/`); `source[]`
is a list of the upstream JSON sources under `data/` and recovered
cross-check JSONs under `certified/` from which the certificate is
assembled. The
`MANIFEST.sha256` digest is computed over the canonical serialization
`json.dumps(cert, sort_keys=True, indent=2, ensure_ascii=False)` with
the `produced_utc` field removed.

## Reproducibility entry-points

```
python verify_all.py            # claim-level checks (seconds)
python reproduce_all.py         # full builder rerun + byte-equality
```

`verify_all.py` reads `MANIFEST.sha256`, validates schema, recomputes
content hashes, and re-checks the principal claims (ranks, SNF
invariant factors of the listed witnesses, kernel verifications,
group orders / SmallGroup IDs).

`reproduce_all.py` re-imports `build_certified.py` and invokes every
builder, asserting byte-equality of the canonical JSON output (modulo
the `produced_utc` field). It does **not** re-execute the upstream
producer scripts under `scripts/legacy/` — those have already written
the authoritative source JSONs under `data/`, whose content is pinned
by their own sha256 inside each certificate's `source[]` block. Exit
code 0 iff every certificate reproduces.

## Seed policy

All randomized steps (only the $n\ge 6$ random-walk samplers) read
their seed from `seeds/seedlist.txt`. Recomputation is fully
deterministic.
