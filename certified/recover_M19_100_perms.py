"""recover_M19_100_perms.py
==============================
Re-runs the exhaustive S_9 scan of M_{19} to recover the full list of 100
rank-7 symbol relabelings (the upstream script `phase_9_18_S9x.py` saves
`Gamma_star_perms` of size 72 — the shared left-kernel / middle-band subset — but does not
serialize the full rank-7 list of size 100).

Two-stage protocol:
  1. Exact modular determinant scan over all 9! = 362880 relabelings; flag
      candidates whose projected 8x8 determinant is zero modulo a large prime.
      If the rational determinant is zero, it is zero modulo every prime, so no
      rank-7 relabeling can be missed by this filter.
  2. Exact Sympy rational rank verification on every flagged candidate:
     compute rank of the integer centered operator E_{gamma, M19} = gamma(M19) - 5*J
     by Sympy.Matrix.rank().

Output: `M19_100_perms_recovered.json` next to this script. Then read by
`build_certified.py` and folded into `M19_audit.json` and
`lift_lemma_evidence.json`.

Runtime: ~30-60 seconds.
"""
from __future__ import annotations

import json
import time
from itertools import permutations
from pathlib import Path

import numpy as np
from sympy import Matrix

HERE = Path(__file__).resolve().parent
# HERE = <repo>/certified/  -> REPO_ROOT = <repo>/
REPO_ROOT = HERE.parent
DATA = REPO_ROOT / "data"
MODULUS = 1_000_003


def load_M19() -> np.ndarray:
    with (DATA / "phase_9_18_S9x_results.json").open("r", encoding="utf-8") as f:
        d = json.load(f)
    return np.array(d["M19_grid"], dtype=np.int64)


def exact_rank(L: np.ndarray, gamma: tuple[int, ...]) -> int:
    """Exact rank of E_{gamma, L} = gamma(L) - 5*J via Sympy."""
    g = np.asarray(gamma, dtype=np.int64)
    Lg = g[L - 1]
    E = (Lg - 5).tolist()
    return int(Matrix(E).rank())


def projected_ehat(Lg: np.ndarray) -> np.ndarray:
    """Return P^T E P for P=(e_i-e_9)_{i=1}^8 and E=Lg-5J."""
    E = Lg.astype(np.int64) - 5
    return E[:8, :8] - E[:8, 8:9] - E[8:9, :8] + E[8, 8]


def det_mod(A: np.ndarray, p: int) -> int:
    """Determinant modulo prime p by Gaussian elimination."""
    rows = [[int(x % p) for x in row] for row in A.tolist()]
    n = len(rows)
    det = 1
    for i in range(n):
        pivot = None
        for r in range(i, n):
            if rows[r][i] % p:
                pivot = r
                break
        if pivot is None:
            return 0
        if pivot != i:
            rows[i], rows[pivot] = rows[pivot], rows[i]
            det = (-det) % p
        pivot_value = rows[i][i] % p
        det = (det * pivot_value) % p
        inverse = pow(pivot_value, p - 2, p)
        for r in range(i + 1, n):
            if rows[r][i] % p:
                factor = (rows[r][i] * inverse) % p
                for c in range(i, n):
                    rows[r][c] = (rows[r][c] - factor * rows[i][c]) % p
    return det % p


def main() -> int:
    M19 = load_M19()

    t0 = time.time()
    # Stage 1: modular determinant scan, collect all candidates whose rational
    # projected determinant may be zero.
    candidates: list[tuple[int, ...]] = []
    n_total = 0
    for perm in permutations(range(1, 10)):
        n_total += 1
        g = np.asarray(perm, dtype=np.int64)
        Lg = g[M19 - 1]
        if det_mod(projected_ehat(Lg), MODULUS) == 0:
            candidates.append(perm)
        if n_total % 60480 == 0:
            print(
                f"  stage1 {n_total:>6}/362880  elapsed={time.time()-t0:.1f}s "
                f"candidates={len(candidates)}",
                flush=True,
            )
    t1 = time.time()
    print(f"[stage1] scanned {n_total} in {t1-t0:.1f}s -> {len(candidates)} candidates")

    # Stage 2: exact Sympy verification.
    rank7: list[list[int]] = []
    rank_dist: dict[int, int] = {}
    for perm in candidates:
        r = exact_rank(M19, perm)
        rank_dist[r] = rank_dist.get(r, 0) + 1
        if r == 7:
            rank7.append(list(perm))
    t2 = time.time()
    print(
        f"[stage2] sympy verification on {len(candidates)} candidates "
        f"in {t2-t1:.1f}s"
    )
    print(f"         exact rank distribution among candidates: {rank_dist}")

    # Every non-candidate has determinant nonzero modulo MODULUS, hence nonzero
    # over Z and full projected rank 8 over Q. False positives are harmless and
    # removed by the exact Sympy rank check above.
    full_dist = {8: n_total - len(candidates)}
    for r, c in rank_dist.items():
        full_dist[r] = full_dist.get(r, 0) + c
    print(f"[combined] full rank distribution: {full_dist}")

    elapsed = t2 - t0
    out = {
        "phase": "M19 S_9 recovery (full 100 rank-7 perms, exact)",
        "M19_grid": M19.tolist(),
        "n_total": n_total,
        "rank_distribution": {str(k): v for k, v in full_dist.items() if v},
        "rank7_count": len(rank7),
        "rank7_perms": sorted(rank7),
        "modular_prime": MODULUS,
        "modular_candidate_count": len(candidates),
        "stage1_elapsed_seconds": t1 - t0,
        "stage2_elapsed_seconds": t2 - t1,
        "elapsed_seconds": elapsed,
        "method": "two-stage: exact modular determinant scan + sympy exact rank verification",
    }
    out_path = HERE / "M19_100_perms_recovered.json"
    with out_path.open("w", encoding="utf-8") as f:
        json.dump(out, f, sort_keys=True, indent=2)
    print(f"[OK] wrote {out_path.name}  ({len(rank7)} rank-7 perms)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
