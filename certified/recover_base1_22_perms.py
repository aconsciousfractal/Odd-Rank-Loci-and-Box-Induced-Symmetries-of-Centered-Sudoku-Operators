"""recover_base1_22_perms.py
==============================
Re-runs the exhaustive S_9 scan of base1 to recover the full list of 22
odd-rank symbol relabelings (the original phase_9_13.py truncates
``odd_examples`` to 5 for log brevity; the count 22 is preserved).

Output: ``base1_22_perms_recovered.json`` next to this script. The
output is then read by ``build_certified.py`` and folded into
``base1_odd_rank_22.json``.

Runtime: ~2-4 minutes.
"""
from __future__ import annotations

import json
import time
from itertools import permutations
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
# HERE = <repo>/certified/  -> REPO_ROOT = <repo>/
REPO_ROOT = HERE.parent
DATA = REPO_ROOT / "data"

# base1 grid as recorded in data/phase_9_13_results.json -> phase_9_13_c.base_grid
BASE1 = np.array([
    [5, 3, 4, 6, 7, 8, 9, 1, 2],
    [6, 7, 2, 1, 9, 5, 3, 4, 8],
    [1, 9, 8, 3, 4, 2, 5, 6, 7],
    [8, 5, 9, 7, 6, 1, 4, 2, 3],
    [4, 2, 6, 8, 5, 3, 7, 9, 1],
    [7, 1, 3, 9, 2, 4, 8, 5, 6],
    [9, 6, 1, 5, 3, 7, 2, 8, 4],
    [2, 8, 7, 4, 1, 9, 6, 3, 5],
    [3, 4, 5, 2, 8, 6, 1, 7, 9],
], dtype=np.int64)


def main() -> int:
    # Sanity: confirm base1 matches the recorded grid in data JSON.
    with (DATA / "phase_9_13_results.json").open("r", encoding="utf-8") as f:
        d = json.load(f)
    recorded = np.array(d["phase_9_13_c"]["base_grid"], dtype=np.int64)
    assert np.array_equal(recorded, BASE1), "base1 mismatch with recorded grid"

    t0 = time.time()
    odd_perms: list[list[int]] = []
    rank_counts = {7: 0, 8: 0}
    n_total = 0

    # Scan all 9! relabelings.
    for perm in permutations(range(1, 10)):
        n_total += 1
        gamma = np.asarray(perm, dtype=np.int64)
        Lg = gamma[BASE1 - 1]                  # symbol relabeling
        E = Lg.astype(np.float64) - 5.0        # centered operator
        r = int(round(np.linalg.matrix_rank(E)))
        if r == 7:
            odd_perms.append(list(perm))
            rank_counts[7] += 1
        else:
            rank_counts[8] = rank_counts.get(8, 0) + 1
        if n_total % 60480 == 0:
            print(f"  {n_total:>6}/362880  elapsed={time.time()-t0:.1f}s "
                  f"odd={rank_counts[7]}", flush=True)

    elapsed = time.time() - t0
    print(f"\n[OK] scanned {n_total} relabelings in {elapsed:.1f}s")
    print(f"     rank 7: {rank_counts[7]}   rank 8: {rank_counts[8]}")

    out = {
        "phase": "base1 S_9 recovery (full 22 perms)",
        "base_grid": BASE1.tolist(),
        "n_total": n_total,
        "rank_distribution": {str(k): v for k, v in rank_counts.items() if v},
        "odd_rank_count": rank_counts[7],
        "odd_perms": odd_perms,
        "elapsed_seconds": elapsed,
    }
    out_path = HERE / "base1_22_perms_recovered.json"
    with out_path.open("w", encoding="utf-8") as f:
        json.dump(out, f, sort_keys=True, indent=2)
    print(f"[OK] wrote {out_path.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
