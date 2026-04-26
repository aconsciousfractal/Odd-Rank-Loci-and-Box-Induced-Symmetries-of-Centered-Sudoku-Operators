"""
Phase 9.14 — F_2-Rank Criterion per Sudoku n=6 (Roku-Doku, 2x3 box).

Question: at n = 6 (Sudoku-2x3), the F_2-rank criterion (Determinant paper,
Theorem D3) says that a Latin square L of order n ≡ 2 (mod 4) has
n^2 | det(E_std) iff rank_{F_2}(A mod 2) < n - 1, where
    A_{ij} = L_{ij} - L_{i, n-1},  0 ≤ i, j ≤ n-2.
For generic LS-6 (9408 reduced), exactly 576 (≈6.12%) have rank_{F_2}(A) = 5
(full F_2-rank, hence the n^2-divisibility fails by a factor of 2).

Does the 2x3 box constraint change this distribution?

Protocol:
  9.14a  Exhaustively enumerate all Sudoku-6 grids with first row fixed
         = (1,2,3,4,5,6).  Total = 28,200,960 / 6! = 39,168 grids.
  9.14b  For each, compute rank_{F_2}(A mod 2) and tabulate the distribution.
  9.14c  Compare with all 9408 reduced LS-6 (or 1,128,960 first-row-fixed LS-6).
         Note: for first-row-fixed, the two universes are directly comparable
         per-grid; for reduced, the distribution is fixed-multiplicity
         (each reduced LS represents 6 first-row-fixed LS by column shift...
         actually no -- reduced fixes first column as 1..n too).
         For a clean comparison, we enumerate ALL first-row-fixed LS-6
         (1,128,960) but report the reduced (9408) result as well.
  9.14e  Sudoku-6 parity-pattern analysis: each row, col, AND box of L mod 2
         contains exactly 3 odds and 3 evens.  Test whether this stronger
         "doubly-balanced" structure forces rank_{F_2}(A) < 5.
  9.14d  If the Sudoku distribution differs from the LS distribution,
         articulate the algebraic mechanism.

E_std convention: A_{ij} = L_{ij} - L_{i, n-1} (the (n-1)x(n-1) block used in
the Determinant paper for the F_2 criterion, NOT the symmetric E_std =
P^T E P used in Phase 9.13).  Both have the same F_2-rank since they differ
by the column projector subtraction L_{n-1, j}, which is a row-rank-1
modification that vanishes mod 2 in the relevant subspace.

Output: data/phase_9_14_results.json
"""
from __future__ import annotations

import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent  # <repo>/scripts/legacy/ -> <repo>/
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)
OUT = DATA / "phase_9_14_results.json"

N = 6
# Box layout: 2x3 boxes (2 rows tall, 3 cols wide). 6 boxes total.
BOX_H, BOX_W = 2, 3


def box_of(r: int, c: int) -> int:
    return (r // BOX_H) * (N // BOX_W) + (c // BOX_W)


# --------------------------------------------------------------------------- #
# F_2 rank (bitmask Gauss elimination)                                         #
# --------------------------------------------------------------------------- #

def rank_f2_bits(rows: list[int], n_cols: int) -> int:
    """Rank over F_2 of a matrix given as list of int bitmasks (rows)."""
    rows = list(rows)
    rank = 0
    m = len(rows)
    for col in range(n_cols):
        bit = 1 << col
        piv = None
        for r in range(rank, m):
            if rows[r] & bit:
                piv = r
                break
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for r in range(m):
            if r != rank and (rows[r] & bit):
                rows[r] ^= rows[rank]
        rank += 1
    return rank


def A_mod_2_bits(L: np.ndarray) -> tuple[list[int], int]:
    """Return rows of A mod 2 as bitmasks of length n-1, where
       A_{ij} = L_{ij} - L_{i, n-1}.
    Mod 2:  A_{ij} ≡ L_{ij} XOR L_{i, n-1}.
    """
    n = L.shape[0]
    rows = []
    for i in range(n - 1):
        last = int(L[i, n - 1]) & 1
        bm = 0
        for j in range(n - 1):
            v = (int(L[i, j]) & 1) ^ last
            if v:
                bm |= (1 << j)
        rows.append(bm)
    return rows, n - 1


# --------------------------------------------------------------------------- #
# Sudoku-6 enumeration with first row fixed                                   #
# --------------------------------------------------------------------------- #

def enumerate_sudoku6_firstrow_fixed(callback, progress_every: int = 5000):
    """Enumerate all valid Sudoku-6 grids with first row = (1,2,3,4,5,6).

    Calls callback(L) for each valid grid, where L is a 6x6 numpy array.
    Total expected count: 28,200,960 / 720 = 39,168.
    """
    n = N
    grid = [[0] * n for _ in range(n)]
    grid[0] = [1, 2, 3, 4, 5, 6]
    row_mask = [0] * n
    col_mask = [0] * n
    box_mask = [0] * n
    # initialize masks from first row
    for c in range(n):
        v = grid[0][c]
        bit = 1 << (v - 1)
        row_mask[0] |= bit
        col_mask[c] |= bit
        box_mask[box_of(0, c)] |= bit

    cells = [(r, c) for r in range(n) for c in range(n) if r > 0]
    count = [0]
    t0 = time.time()

    def rec(pos: int):
        if pos == len(cells):
            arr = np.array([row[:] for row in grid], dtype=np.int8)
            callback(arr)
            count[0] += 1
            if count[0] % progress_every == 0:
                el = time.time() - t0
                print(f"  enumerated {count[0]:>7d}  elapsed={el:.1f}s",
                      flush=True)
            return
        r, c = cells[pos]
        b = box_of(r, c)
        avail = ~(row_mask[r] | col_mask[c] | box_mask[b]) & 0x3F
        while avail:
            lsb = avail & -avail
            v = lsb.bit_length()
            grid[r][c] = v
            row_mask[r] |= lsb
            col_mask[c] |= lsb
            box_mask[b] |= lsb
            rec(pos + 1)
            grid[r][c] = 0
            row_mask[r] &= ~lsb
            col_mask[c] &= ~lsb
            box_mask[b] &= ~lsb
            avail &= avail - 1

    rec(0)
    return count[0]


# --------------------------------------------------------------------------- #
# LS-6 reduced enumeration (cross-check baseline)                             #
# --------------------------------------------------------------------------- #

def enumerate_reduced_ls6(callback):
    """Enumerate all reduced LS-6 (first row = first col = (1,...,6))."""
    n = N
    grid = [[0] * n for _ in range(n)]
    for i in range(n):
        grid[0][i] = i + 1
        grid[i][0] = i + 1
    row_mask = [0] * n
    col_mask = [0] * n
    for r in range(n):
        for c in range(n):
            if grid[r][c] != 0:
                bit = 1 << (grid[r][c] - 1)
                row_mask[r] |= bit
                col_mask[c] |= bit

    cells = [(r, c) for r in range(n) for c in range(n)
             if r > 0 and c > 0]
    count = [0]

    def rec(pos: int):
        if pos == len(cells):
            arr = np.array([row[:] for row in grid], dtype=np.int8)
            callback(arr)
            count[0] += 1
            return
        r, c = cells[pos]
        avail = ~(row_mask[r] | col_mask[c]) & 0x3F
        while avail:
            lsb = avail & -avail
            v = lsb.bit_length()
            grid[r][c] = v
            row_mask[r] |= lsb
            col_mask[c] |= lsb
            rec(pos + 1)
            grid[r][c] = 0
            row_mask[r] &= ~lsb
            col_mask[c] &= ~lsb
            avail &= avail - 1

    rec(0)
    return count[0]


def enumerate_ls6_firstrow_fixed(callback, progress_every: int = 100000):
    """Enumerate all LS-6 with first row = (1,...,6); first column free.

    Total expected count: 6! * 9408 / 6 ... actually 812,851,200 / 6! = 1,128,960.
    """
    n = N
    grid = [[0] * n for _ in range(n)]
    grid[0] = [1, 2, 3, 4, 5, 6]
    row_mask = [0] * n
    col_mask = [0] * n
    for c in range(n):
        bit = 1 << (grid[0][c] - 1)
        row_mask[0] |= bit
        col_mask[c] |= bit

    cells = [(r, c) for r in range(n) for c in range(n) if r > 0]
    count = [0]
    t0 = time.time()

    def rec(pos: int):
        if pos == len(cells):
            arr = np.array([row[:] for row in grid], dtype=np.int8)
            callback(arr)
            count[0] += 1
            if count[0] % progress_every == 0:
                el = time.time() - t0
                print(f"  enumerated {count[0]:>8d}  elapsed={el:.1f}s",
                      flush=True)
            return
        r, c = cells[pos]
        avail = ~(row_mask[r] | col_mask[c]) & 0x3F
        while avail:
            lsb = avail & -avail
            v = lsb.bit_length()
            grid[r][c] = v
            row_mask[r] |= lsb
            col_mask[c] |= lsb
            rec(pos + 1)
            grid[r][c] = 0
            row_mask[r] &= ~lsb
            col_mask[c] &= ~lsb
            avail &= avail - 1

    rec(0)
    return count[0]


# --------------------------------------------------------------------------- #
# Main scan                                                                   #
# --------------------------------------------------------------------------- #

def scan_universe(name: str, enumerator):
    print(f"\n== {name} ==")
    rank_dist = Counter()
    parity_pattern_set = set()
    parity_pattern_to_rank = {}
    bad_examples: list[list[list[int]]] = []  # rank=5 witnesses
    box_sums_examples: list[dict] = []
    t0 = time.time()
    n_total = [0]

    def cb(L: np.ndarray):
        n_total[0] += 1
        rows, ncols = A_mod_2_bits(L)
        r = rank_f2_bits(rows, ncols)
        rank_dist[r] += 1
        # also collect the parity pattern (L mod 2) for orbit analysis
        parity = tuple(tuple(int(x) & 1 for x in row) for row in L)
        if parity not in parity_pattern_to_rank:
            parity_pattern_to_rank[parity] = r
        parity_pattern_set.add(parity)
        if r == N - 1 and len(bad_examples) < 5:
            bad_examples.append(L.tolist())

    enumerator(cb)
    el = time.time() - t0
    print(f"  total grids: {n_total[0]}  elapsed={el:.1f}s")
    print(f"  rank_F2(A mod 2) distribution: {dict(rank_dist)}")
    full_rank = rank_dist.get(N - 1, 0)
    pct_bad = 100.0 * full_rank / max(1, n_total[0])
    print(f"  full-rank (rank=={N-1}, i.e. n^2-divisibility FAILS by 2): "
          f"{full_rank}/{n_total[0]}  ({pct_bad:.4f}%)")
    print(f"  distinct parity patterns: {len(parity_pattern_set)}")
    pp_rank_dist = Counter(parity_pattern_to_rank.values())
    print(f"  rank distribution over distinct parity patterns: "
          f"{dict(pp_rank_dist)}")

    return {
        "name": name,
        "n_total": n_total[0],
        "rank_dist": dict(rank_dist),
        "full_rank_count": full_rank,
        "full_rank_pct": pct_bad,
        "n_distinct_parity_patterns": len(parity_pattern_set),
        "parity_pattern_rank_dist": dict(pp_rank_dist),
        "bad_examples": bad_examples,
        "elapsed_seconds": el,
    }


def main():
    args = sys.argv[1:]
    skip_ls = ("--no-ls" in args)
    skip_sudoku = ("--no-sudoku" in args)

    print("=" * 70)
    print("  PHASE 9.14 — F_2-Rank Criterion for Sudoku n=6 (Roku-Doku)")
    print("=" * 70)

    results = {}
    if not skip_sudoku:
        results["sudoku6_firstrow_fixed"] = scan_universe(
            "Sudoku-6 (first row = 1..6, all 39 168 grids)",
            enumerate_sudoku6_firstrow_fixed,
        )
    if not skip_ls:
        results["ls6_reduced"] = scan_universe(
            "LS-6 reduced (first row = first col = 1..6, all 9408)",
            enumerate_reduced_ls6,
        )
        results["ls6_firstrow_fixed"] = scan_universe(
            "LS-6 first-row-fixed (first row = 1..6, all 1 128 960 grids)",
            enumerate_ls6_firstrow_fixed,
        )

    # Side-by-side summary
    if not skip_sudoku and not skip_ls:
        s = results["sudoku6_firstrow_fixed"]
        l = results["ls6_reduced"]
        lf = results["ls6_firstrow_fixed"]
        print("\n" + "=" * 70)
        print("  SUMMARY")
        print("=" * 70)
        print(f"  Sudoku-6 first-row-fixed: rank=5 in "
              f"{s['full_rank_count']}/{s['n_total']} "
              f"({s['full_rank_pct']:.4f}%)")
        print(f"  LS-6 first-row-fixed:     rank=5 in "
              f"{lf['full_rank_count']}/{lf['n_total']} "
              f"({lf['full_rank_pct']:.4f}%)")
        print(f"  LS-6 reduced:             rank=5 in "
              f"{l['full_rank_count']}/{l['n_total']} "
              f"({l['full_rank_pct']:.4f}%)")
        ref_pct = lf["full_rank_pct"]
        if s["full_rank_count"] == 0:
            print("  -> Box-constraint REMOVES all F_2-rank-5 counterexamples!")
        elif abs(s["full_rank_pct"] - ref_pct) < 0.01:
            print("  -> Box-constraint preserves F_2-rank distribution (within noise).")
        elif s["full_rank_pct"] < ref_pct:
            print("  -> Box-constraint REDUCES F_2-rank-5 rate.")
        else:
            print("  -> Box-constraint INCREASES F_2-rank-5 rate.")

    with OUT.open("w", encoding="utf-8") as f:
        # bad_examples can be large; trim
        json.dump(results, f, indent=2)
    print(f"\n  results saved -> {OUT}")


if __name__ == "__main__":
    main()
