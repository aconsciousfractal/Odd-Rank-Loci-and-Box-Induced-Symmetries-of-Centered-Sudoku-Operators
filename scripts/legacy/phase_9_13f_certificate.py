"""
Phase 9.13f — 4-way certificate for an odd-rank Sudoku relabeling.

Given a permutation gamma in S_9 and a base Sudoku L, verify by four
independent methods that rank(E_std(gamma o L)) = 7 (odd):

  (1) Bareiss fraction-free determinant of a 7x7 minor != 0.
  (2) Q-RREF (Fraction-based Gaussian elimination) returns rank 7.
  (3) rank over F_p = 7 for p in {101, 1009, 10007, 100003}.
  (4) SymPy's .rank() returns 7.

We also verify det(E_std) = 0 exactly.

Inputs are loaded from data/phase_9_13_results.json (the odd_examples list).
Output: data/phase_9_13f_certificate.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
from sympy import Matrix as SMatrix

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "data"  # <repo>/scripts/legacy/ -> <repo>/data/
RESULTS = DATA / "phase_9_13_results.json"
OUT = DATA / "phase_9_13f_certificate.json"


def centered_S(L: np.ndarray) -> np.ndarray:
    """Return S = P^T (L - mu J) P as (n-1)x(n-1) integer matrix."""
    n = L.shape[0]
    E = L - (n + 1) // 2 if n % 2 == 1 else L - (n + 1) / 2
    E = (L * 2 - (n + 1))  # integer-scaled for odd n (doubles rank-equivalent)
    # But for odd n we actually keep the integer version L-5; it is already int.
    E = L.astype(int) - (n + 1) // 2
    # This is integer for odd n since (n+1)/2 is integer.
    last = int(E[n - 1, n - 1])
    S = np.zeros((n - 1, n - 1), dtype=int)
    for i in range(n - 1):
        for j in range(n - 1):
            S[i, j] = int(E[i, j]) - int(E[i, n - 1]) - int(E[n - 1, j]) + last
    return S


def bareiss_det(A: np.ndarray) -> int:
    """Fraction-free Bareiss determinant for an integer square matrix."""
    n = A.shape[0]
    M = A.astype(object).copy()  # Python ints
    for i in range(n):
        M[i, :] = [int(x) for x in M[i, :]]
    sign = 1
    prev = 1
    for k in range(n - 1):
        # find pivot
        if M[k, k] == 0:
            piv = -1
            for r in range(k + 1, n):
                if M[r, k] != 0:
                    piv = r
                    break
            if piv < 0:
                return 0
            M[[k, piv], :] = M[[piv, k], :]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i, j] = (M[i, j] * M[k, k] - M[i, k] * M[k, j]) // prev
            M[i, k] = 0
        prev = M[k, k]
    return int(sign * M[n - 1, n - 1])


def rank_Q(A: np.ndarray) -> int:
    """Rank via Gaussian elimination over Q (Fraction)."""
    M = [[Fraction(int(x)) for x in row] for row in A]
    rows = len(M)
    cols = len(M[0]) if rows else 0
    r = 0
    for c in range(cols):
        piv = None
        for k in range(r, rows):
            if M[k][c] != 0:
                piv = k
                break
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        for k in range(rows):
            if k == r or M[k][c] == 0:
                continue
            f = M[k][c] / M[r][c]
            for j in range(c, cols):
                M[k][j] -= f * M[r][j]
        r += 1
    return r


def rank_Fp(A: np.ndarray, p: int) -> int:
    """Rank of A mod p via full Gaussian elimination in F_p."""
    rows, cols = A.shape
    M = [[int(x) % p for x in row] for row in A]
    r = 0
    for c in range(cols):
        piv = None
        for k in range(r, rows):
            if M[k][c] != 0:
                piv = k
                break
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], -1, p)
        for k in range(rows):
            if k == r or M[k][c] == 0:
                continue
            f = (M[k][c] * inv) % p
            for j in range(c, cols):
                M[k][j] = (M[k][j] - f * M[r][j]) % p
        r += 1
    return r


def find_nonzero_minor(A: np.ndarray, r: int):
    """Find a r-by-r submatrix with nonzero Bareiss determinant."""
    from itertools import combinations
    rows, cols = A.shape
    for R in combinations(range(rows), r):
        for C in combinations(range(cols), r):
            sub = A[np.ix_(R, C)]
            d = bareiss_det(sub)
            if d != 0:
                return {"rows": list(R), "cols": list(C), "det": d}
    return None


def main():
    with RESULTS.open("r", encoding="utf-8") as f:
        res = json.load(f)
    base_grid = np.array(res["phase_9_13_c"]["base_grid"], dtype=int)
    odd_examples = res["phase_9_13_c"]["odd_examples"]
    if not odd_examples:
        print("No odd-rank examples found — nothing to certify.")
        return
    print("=" * 70)
    print("  PHASE 9.13f — 4-way certificate for odd-rank Sudoku relabelings")
    print("=" * 70)
    certificates = []
    for ex in odd_examples[:3]:
        perm = np.array(ex["perm"], dtype=int)  # 1-indexed image
        expected_rank = int(ex["rank"])
        gamma = np.concatenate([[0], perm])  # index by value 1..9
        L_new = gamma[base_grid]
        S = centered_S(L_new)
        assert S.shape == (8, 8)
        n = 8
        print(f"\n-- perm = {perm.tolist()}  (expected rank={expected_rank}) --")
        # (1) determinant of full 8x8 must be zero
        det_full = bareiss_det(S)
        print(f"  [1] Bareiss det(S_8x8) = {det_full}  (expect 0)")
        # (2) Q-RREF rank
        r_Q = rank_Q(S)
        print(f"  [2] rank over Q (Fraction RREF) = {r_Q}")
        # (3) rank mod p
        mod_ranks = {}
        for p in [101, 1009, 10007, 100003]:
            rp = rank_Fp(S, p)
            mod_ranks[p] = rp
            print(f"  [3] rank mod {p:>7} = {rp}")
        # (4) SymPy rank
        r_sympy = int(SMatrix(S.tolist()).rank())
        print(f"  [4] SymPy rank = {r_sympy}")
        # (5) nonzero 7x7 minor as witness of rank >= 7
        witness = find_nonzero_minor(S, 7)
        if witness is not None:
            print(
                f"  [W] witness 7x7 minor rows={witness['rows']} "
                f"cols={witness['cols']} det={witness['det']}"
            )
        else:
            print("  [W] NO nonzero 7x7 minor found (unexpected)")
        consistent = (
            det_full == 0
            and r_Q == expected_rank
            and r_sympy == expected_rank
            and all(v == expected_rank for v in mod_ranks.values())
            and witness is not None
        )
        print(f"  -> CERTIFIED: {consistent}")
        certificates.append({
            "perm": ex["perm"],
            "expected_rank": expected_rank,
            "S_matrix": S.tolist(),
            "det_Bareiss": det_full,
            "rank_Q": r_Q,
            "rank_sympy": r_sympy,
            "rank_Fp": {str(p): r for p, r in mod_ranks.items()},
            "witness_minor": witness,
            "certified": consistent,
        })
    with OUT.open("w", encoding="utf-8") as f:
        json.dump({"certificates": certificates, "base_grid": base_grid.tolist()},
                  f, indent=2)
    print(f"\n  certificates saved -> {OUT}")


if __name__ == "__main__":
    main()
