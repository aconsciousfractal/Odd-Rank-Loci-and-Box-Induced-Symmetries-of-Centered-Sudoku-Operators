"""
Phase 9.13 — Box-Constraint Rank-Parity Forcing (Sudoku n=9).

Question: does the Sudoku box constraint force rank(E_std) to be even
          when L is a valid Sudoku grid (independent of relabeling)?

Protocol:
  9.13a  Sample a diverse corpus of valid Sudoku grids at standard labeling
         (randomized backtracking; independent of the JM/LS sampler).
  9.13b  Compute rank(E_std) and SNF(E_std) for each grid, standard labeling.
  9.13c  Pick one base Sudoku grid and scan ALL 9! = 362880 relabelings γ∈S_9,
         recording rank(E_std) distribution across the orbit.
  9.13e  Test the box-block decomposition of E_std:
         In basis (band-mean ⊕ within-band) for rows and
                  (stack-mean ⊕ within-stack) for columns,
         E_std has structure [[0, A], [B, C]] where the 2x2 top-left block is 0
         (forced by the Sudoku box constraint).  We compute ranks of A, B, C
         on the whole corpus and check whether the global rank is always even.

Output: data/phase_9_13_results.json  and a textual summary to stdout.

E_std convention:
  Let L be 9x9 with entries in {1,...,9}. Define
      S_{i j} = L_{i j} - L_{i,n-1} - L_{n-1,j} + L_{n-1,n-1}    (i,j in 0..n-2)
  Then S = P^T (L - mu J) P with P = [e_i - e_{n-1}] the standard projector,
  so S is an (n-1)x(n-1) integer matrix with the same rank and integer SNF
  as E_std in the Parity/Determinant papers.
"""
from __future__ import annotations

import json
import random
import sys
import time
from collections import Counter
from fractions import Fraction
from itertools import permutations
from pathlib import Path

import numpy as np
from sympy import Matrix as SMatrix, ZZ
from sympy.matrices.normalforms import smith_normal_form as _sympy_snf

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent  # <repo>/scripts/legacy/ -> <repo>/
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)
OUT = DATA / "phase_9_13_results.json"


# --------------------------------------------------------------------------- #
# Sudoku sampler: randomized backtracking                                     #
# --------------------------------------------------------------------------- #

def _box_of(r: int, c: int) -> int:
    return (r // 3) * 3 + (c // 3)


def _solve_random(grid, row_mask, col_mask, box_mask, order, rng):
    """Fill a 9x9 grid via DFS with randomized value choice.

    Args:
      grid:     flat list of 81 ints (values 0..9; 0 = empty)
      row_mask, col_mask, box_mask: bitmasks of digits placed (bit v-1)
      order:    list of 81 cell indices giving the fill order (also randomized)
      rng:      random.Random
    """
    for pos in range(81):
        idx = order[pos]
        if grid[idx] != 0:
            continue
        r, c = divmod(idx, 9)
        b = _box_of(r, c)
        avail = ~(row_mask[r] | col_mask[c] | box_mask[b]) & 0x1FF
        if avail == 0:
            return False
        # randomize digit order
        digits = [v + 1 for v in range(9) if avail & (1 << v)]
        rng.shuffle(digits)
        for v in digits:
            bit = 1 << (v - 1)
            grid[idx] = v
            row_mask[r] |= bit
            col_mask[c] |= bit
            box_mask[b] |= bit
            if _solve_random(grid, row_mask, col_mask, box_mask, order, rng):
                return True
            grid[idx] = 0
            row_mask[r] &= ~bit
            col_mask[c] &= ~bit
            box_mask[b] &= ~bit
        return False
    return True


def sample_sudoku(rng: random.Random):
    """Produce one random-ish valid Sudoku grid at standard labeling.

    Strategy: fill the 3 diagonal 3x3 boxes with random permutations (they
    are independent), then solve the remainder by natural-order backtracking
    with randomized value choice.  This is the standard fast generator.

    Returns a numpy array of shape (9,9) with values 1..9.
    """
    grid = [0] * 81
    row_mask = [0] * 9
    col_mask = [0] * 9
    box_mask = [0] * 9
    # Fill the three diagonal boxes independently with random permutations.
    for bk in range(3):
        perm = list(range(1, 10))
        rng.shuffle(perm)
        t = 0
        for a in range(3):
            for b in range(3):
                r = 3 * bk + a
                c = 3 * bk + b
                v = perm[t]
                t += 1
                bit = 1 << (v - 1)
                grid[9 * r + c] = v
                row_mask[r] |= bit
                col_mask[c] |= bit
                box_mask[bk * 3 + bk] |= bit  # box index
        # NOTE: box index of (r,c) is (r//3)*3 + (c//3); for diagonal it is bk*3+bk
    # Natural row-major order; only empty cells remain.
    order = list(range(81))
    ok = _solve_random(grid, row_mask, col_mask, box_mask, order, rng)
    if not ok:
        raise RuntimeError("sudoku backtracking failed unexpectedly")
    return np.array(grid, dtype=int).reshape(9, 9)


def is_sudoku(L: np.ndarray) -> bool:
    if L.shape != (9, 9):
        return False
    s = set(range(1, 10))
    for i in range(9):
        if set(L[i]) != s or set(L[:, i]) != s:
            return False
    for br in range(3):
        for bc in range(3):
            vals = {int(L[3 * br + a, 3 * bc + b]) for a in range(3) for b in range(3)}
            if vals != s:
                return False
    return True


# --------------------------------------------------------------------------- #
# E_std = S as (n-1) x (n-1) integer matrix                                   #
# --------------------------------------------------------------------------- #

def E_std_matrix(L: np.ndarray):
    """Return E_std as a tuple of tuples (n-1) x (n-1), integer entries."""
    n = L.shape[0]
    last_row = int(L[n - 1, n - 1])
    rows = []
    for i in range(n - 1):
        Li_last = int(L[i, n - 1])
        row = []
        for j in range(n - 1):
            row.append(int(L[i, j]) - Li_last - int(L[n - 1, j]) + last_row)
        rows.append(tuple(row))
    return tuple(rows)


# --------------------------------------------------------------------------- #
# Exact rank (Bareiss-style Gaussian elimination over Q)                      #
# --------------------------------------------------------------------------- #

def exact_rank(M):
    rows = len(M)
    cols = len(M[0]) if rows else 0
    A = [[Fraction(M[i][j]) for j in range(cols)] for i in range(rows)]
    r = 0
    for c in range(cols):
        piv = None
        for k in range(r, rows):
            if A[k][c] != 0:
                piv = k
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv_piv = A[r][c]
        for k in range(rows):
            if k == r or A[k][c] == 0:
                continue
            f = A[k][c] / inv_piv
            for j in range(c, cols):
                A[k][j] -= f * A[r][j]
        r += 1
    return r


# --------------------------------------------------------------------------- #
# Smith Normal Form via SymPy                                                 #
# --------------------------------------------------------------------------- #

def snf_invariants(M):
    """Return the list of diagonal invariants d_1 | d_2 | ... (nonzero entries)."""
    S = SMatrix([list(row) for row in M])
    D = _sympy_snf(S, domain=ZZ)
    invs = []
    r = min(D.rows, D.cols)
    for i in range(r):
        d = int(D[i, i])
        if d != 0:
            invs.append(abs(d))
    return invs


def snf_signature(invs):
    return "[" + ",".join(str(x) for x in invs) + "]"


# --------------------------------------------------------------------------- #
# Box-block decomposition                                                     #
# --------------------------------------------------------------------------- #
#
# Basis for R^9 (rows or cols):
#   band-mean subspace (dim 3): b_0 = (1,1,1,0,0,0,0,0,0)/sqrt(3), ...
#   within-band subspace (dim 6): orthogonal complement (per band)
#
# After quotienting by the overall all-ones vector (row-sum=0 constraint),
# band-mean becomes 2-dim (2 independent band-mean vectors summing to 0),
# within-band stays 6-dim.
#
# E_std lives in (2+6) x (2+6). The Sudoku box constraint P_row E P_col = 0
# (where P is projection onto 3-dim band-mean) forces the 2x2 top-left block
# of E_std in this basis to be zero.

def make_change_of_basis():
    """Return 8x9 matrix Q whose rows form an orthonormal basis of
       {v in R^9 : sum v = 0}, split as first 2 rows = band-mean-with-zero-sum,
       last 6 rows = within-band (zero-mean inside each band).
    """
    # band-mean, zero-sum: 2 basis vectors in the 3-dim band-mean space,
    # orthogonal to (1,...,1).
    bm = np.zeros((2, 9), dtype=float)
    # b1 = band0 - band1
    bm[0, 0:3] = 1.0
    bm[0, 3:6] = -1.0
    # b2 = (band0 + band1 - 2*band2)/sqrt(...) — orthogonal to b1 and to 1.
    bm[1, 0:3] = 1.0
    bm[1, 3:6] = 1.0
    bm[1, 6:9] = -2.0
    # within-band, per band: 2 orthogonal zero-mean vectors on {0,1,2}
    wb_patterns = np.array([
        [1.0, -1.0, 0.0],
        [1.0, 1.0, -2.0],
    ])
    wb = np.zeros((6, 9), dtype=float)
    for k in range(3):
        for t in range(2):
            wb[2 * k + t, 3 * k:3 * k + 3] = wb_patterns[t]
    # normalize rows
    Q = np.vstack([bm, wb])
    for i in range(Q.shape[0]):
        Q[i] /= np.linalg.norm(Q[i])
    return Q


_Q_CACHE = make_change_of_basis()


def exact_box_sum_check(L: np.ndarray):
    """Exact integer check of the Sudoku-Box Band Lemma.

    For every Sudoku grid L, each 3x3 box contains {1,...,9}, hence sums to
    45 = 9*mu where mu = 5.  Equivalently, every 3x3 box of E = L - 5J sums
    to 0.  This is an integer identity, not a floating-point coincidence.

    Equivalent reformulation: P_rowband . E . P_colstack^T = 0 where
    P_rowband averages within rows of a band and P_colstack averages within
    columns of a stack; the composition picks out box-sums.

    Returns dict with:
        'all_45'   : True iff all 9 boxes sum to 45 in L,
        'all_zero' : True iff all 9 boxes sum to 0 in (L - 5J),
        'box_sums' : 3x3 list of integer box-sums of L.
    """
    L = np.asarray(L, dtype=int)
    sums = [[0] * 3 for _ in range(3)]
    for br in range(3):
        for bc in range(3):
            s = 0
            for a in range(3):
                for b in range(3):
                    s += int(L[3 * br + a, 3 * bc + b])
            sums[br][bc] = s
    all_45 = all(sums[br][bc] == 45 for br in range(3) for bc in range(3))
    all_zero = all(sums[br][bc] - 9 * 5 == 0 for br in range(3) for bc in range(3))
    return {"all_45": all_45, "all_zero": all_zero, "box_sums": sums}


def box_block_decomp(L: np.ndarray):
    """Decompose E = L - 5 J in the (band-mean / within-band) basis.

    Returns a dict with keys
        'tl' (2x2), 'tr' (2x6), 'bl' (6x2), 'br' (6x6),
        'rank_full', 'rank_br', 'max_abs_tl'.
    """
    n = L.shape[0]
    assert n == 9
    E = L.astype(float) - 5.0
    Q = _Q_CACHE  # 8 x 9
    Et = Q @ E @ Q.T  # 8 x 8
    tl = Et[:2, :2]
    tr = Et[:2, 2:]
    bl = Et[2:, :2]
    br = Et[2:, 2:]
    return {
        "tl": tl,
        "tr": tr,
        "bl": bl,
        "br": br,
        "max_abs_tl": float(np.max(np.abs(tl))),
        "rank_full": int(np.linalg.matrix_rank(Et, tol=1e-8)),
        "rank_br": int(np.linalg.matrix_rank(br, tol=1e-8)),
    }


# --------------------------------------------------------------------------- #
# Phase 9.13a + 9.13b: corpus + rank/SNF scan                                 #
# --------------------------------------------------------------------------- #

def phase_9_13_ab(n_samples: int, seed: int = 20260424):
    rng = random.Random(seed)
    print(f"\n== 9.13a+b: Sudoku corpus rank/SNF scan ({n_samples} grids) ==")
    rank_dist = Counter()
    snf_dist = Counter()
    grids_info = []
    box_rank_full_dist = Counter()
    box_rank_br_dist = Counter()
    box_tl_max = 0.0
    exact_box_ok_count = 0
    t0 = time.time()
    for k in range(n_samples):
        L = sample_sudoku(rng)
        assert is_sudoku(L)
        S = E_std_matrix(L)
        r = exact_rank(S)
        invs = snf_invariants(S)
        rank_dist[r] += 1
        sig = snf_signature(invs)
        snf_dist[sig] += 1
        # Exact integer Band Lemma check (all 9 boxes sum to 45).
        bs = exact_box_sum_check(L)
        if bs["all_45"] and bs["all_zero"]:
            exact_box_ok_count += 1
        else:
            raise AssertionError(
                f"Box-sum check FAILED on grid {k}: sums={bs['box_sums']}"
            )
        bx = box_block_decomp(L)
        box_rank_full_dist[bx["rank_full"]] += 1
        box_rank_br_dist[bx["rank_br"]] += 1
        box_tl_max = max(box_tl_max, bx["max_abs_tl"])
        grids_info.append({
            "idx": k,
            "rank": r,
            "snf": sig,
            "rank_br": bx["rank_br"],
            "max_abs_tl": bx["max_abs_tl"],
        })
        if (k + 1) % 100 == 0:
            elapsed = time.time() - t0
            print(
                f"  [{k+1}/{n_samples}] elapsed={elapsed:.1f}s "
                f"rank_dist={dict(rank_dist)}"
            )
            sys.stdout.flush()
    elapsed = time.time() - t0
    print(f"  DONE  elapsed={elapsed:.1f}s")
    print(f"  rank(E_std) distribution: {dict(rank_dist)}")
    odd = sum(v for r, v in rank_dist.items() if r % 2 == 1)
    print(f"  ODD rank count: {odd}/{n_samples}")
    print(f"  SNF signatures (top 10):")
    for sig, cnt in snf_dist.most_common(10):
        print(f"    {sig}: {cnt}")
    print(f"  exact box-sum check (Band Lemma, integer identity): "
          f"{exact_box_ok_count}/{n_samples} grids have all 9 boxes summing to 45")
    print(f"  box-decomposition: TL-block max|entry| across corpus = {box_tl_max:.2e}")
    print(f"  rank(br) [6x6 within-within block] distribution: {dict(box_rank_br_dist)}")
    return {
        "n_samples": n_samples,
        "seed": seed,
        "rank_dist": dict(rank_dist),
        "odd_rank_count": odd,
        "snf_dist": dict(snf_dist),
        "box_rank_full_dist": dict(box_rank_full_dist),
        "box_rank_br_dist": dict(box_rank_br_dist),
        "box_tl_max_abs": box_tl_max,
        "exact_box_sum_ok_count": exact_box_ok_count,
        "grids_info": grids_info,
    }


# --------------------------------------------------------------------------- #
# Phase 9.13c: all 9! relabelings on one base Sudoku                          #
# --------------------------------------------------------------------------- #

def phase_9_13_c(L_base: np.ndarray, max_perms: int | None = None):
    print("\n== 9.13c: scan all 9! relabelings on base Sudoku ==")
    total = 362880
    rank_dist = Counter()
    odd_examples = []
    t0 = time.time()
    count = 0
    # Precompute centered base: E_base = L_base - 5 (values -4..4)
    # For relabeling gamma, (gamma o L) has E_{new} = gamma_values[L_ij - 1] - 5.
    # But gamma values are 1..9 permuted; after centering they're a permutation
    # of {-4,-3,...,4}.  So centered entry for value v is gamma_centered[v-1].
    vals = list(range(1, 10))
    for perm in permutations(vals):
        gc = np.array(perm, dtype=int) - 5  # centered image
        # build E_new centered: E_new[i,j] = gc[L_base[i,j] - 1]
        E_new = gc[L_base - 1]  # 9x9 int in {-4..4}
        # compute S via projection on last row/col (= rank-equivalent to 2E_std
        # etc.; here E_new IS the centered version, so S is an (n-1)x(n-1)
        # integer matrix with entries in {-8..8})
        # S_{ij} = E_new[i,j] - E_new[i,8] - E_new[8,j] + E_new[8,8]
        last_row = E_new[8, 8]
        S = (
            E_new[:8, :8]
            - E_new[:8, 8:9]
            - E_new[8:9, :8]
            + last_row
        )
        r = exact_rank([list(row) for row in S.tolist()])
        rank_dist[r] += 1
        if r % 2 == 1 and len(odd_examples) < 5:
            odd_examples.append({"perm": list(perm), "rank": r})
        count += 1
        if max_perms is not None and count >= max_perms:
            break
        if count % 20000 == 0:
            elapsed = time.time() - t0
            print(
                f"  [{count}/{total}] elapsed={elapsed:.1f}s "
                f"rank_dist={dict(rank_dist)}"
            )
            sys.stdout.flush()
    elapsed = time.time() - t0
    print(f"  DONE  count={count}/{total} elapsed={elapsed:.1f}s")
    print(f"  rank distribution over relabelings: {dict(rank_dist)}")
    odd = sum(v for r, v in rank_dist.items() if r % 2 == 1)
    print(f"  ODD rank count: {odd}/{count}")
    if odd_examples:
        print(f"  example odd-rank relabelings:")
        for ex in odd_examples[:3]:
            print(f"    perm={ex['perm']} rank={ex['rank']}")
    return {
        "base_grid": L_base.tolist(),
        "n_perms": count,
        "rank_dist": dict(rank_dist),
        "odd_rank_count": odd,
        "odd_examples": odd_examples,
    }


# --------------------------------------------------------------------------- #
# Main                                                                        #
# --------------------------------------------------------------------------- #

def main():
    args = sys.argv[1:]
    N_CORPUS = int(args[0]) if len(args) >= 1 else 500
    DO_RELABEL = ("--relabel" in args)
    MAX_PERMS_ARG = [a for a in args if a.startswith("--max-perms=")]
    max_perms = int(MAX_PERMS_ARG[0].split("=")[1]) if MAX_PERMS_ARG else None

    print("=" * 70)
    print("  PHASE 9.13 — Box-Constraint Rank-Parity Forcing (Sudoku n=9)")
    print("=" * 70)
    print(f"  corpus size: {N_CORPUS}")
    print(f"  relabel scan: {DO_RELABEL}  (max_perms={max_perms})")

    results = {}
    results["phase_9_13_ab"] = phase_9_13_ab(N_CORPUS)

    if DO_RELABEL:
        # Base grid: classical "first Sudoku" (base1 from twisted_sudoku_lattice)
        L_base = np.array([
            [5, 3, 4, 6, 7, 8, 9, 1, 2],
            [6, 7, 2, 1, 9, 5, 3, 4, 8],
            [1, 9, 8, 3, 4, 2, 5, 6, 7],
            [8, 5, 9, 7, 6, 1, 4, 2, 3],
            [4, 2, 6, 8, 5, 3, 7, 9, 1],
            [7, 1, 3, 9, 2, 4, 8, 5, 6],
            [9, 6, 1, 5, 3, 7, 2, 8, 4],
            [2, 8, 7, 4, 1, 9, 6, 3, 5],
            [3, 4, 5, 2, 8, 6, 1, 7, 9],
        ], dtype=int)
        assert is_sudoku(L_base)
        results["phase_9_13_c"] = phase_9_13_c(L_base, max_perms=max_perms)

    with OUT.open("w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    print(f"\n  results saved -> {OUT}")


if __name__ == "__main__":
    main()
