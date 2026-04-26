"""recover_M19_100_perms.py
==============================
Re-runs the exhaustive S_9 scan of M_{19} to recover the full list of 100
rank-7 symbol relabelings (the upstream script `phase_9_18_S9x.py` saves
`Gamma_star_perms` of size 72 — the Hessian-class subset — but does not
serialize the full rank-7 list of size 100).

Two-stage protocol:
  1. Fast NumPy float scan over all 9! = 362880 relabelings; flag candidates
     with rank < 8.
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


def main() -> int:
    M19 = load_M19()

    t0 = time.time()
    # Stage 1: float scan, collect candidates with float rank < 8.
    candidates: list[tuple[int, ...]] = []
    n_total = 0
    for perm in permutations(range(1, 10)):
        n_total += 1
        g = np.asarray(perm, dtype=np.int64)
        Lg = g[M19 - 1]
        E = Lg.astype(np.float64) - 5.0
        r = int(round(np.linalg.matrix_rank(E)))
        if r < 8:
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

    # Sanity: every non-candidate has float-rank 8; we trust those.
    # Total rank distribution (approximate, but float misclassification flagged
    # by stage2 only against candidates).
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
        "stage1_elapsed_seconds": t1 - t0,
        "stage2_elapsed_seconds": t2 - t1,
        "elapsed_seconds": elapsed,
        "method": "two-stage: numpy float scan + sympy exact rank verification",
    }
    out_path = HERE / "M19_100_perms_recovered.json"
    with out_path.open("w", encoding="utf-8") as f:
        json.dump(out, f, sort_keys=True, indent=2)
    print(f"[OK] wrote {out_path.name}  ({len(rank7)} rank-7 perms)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
