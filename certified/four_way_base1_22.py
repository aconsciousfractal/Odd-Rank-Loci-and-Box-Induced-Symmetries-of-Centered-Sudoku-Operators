"""four_way_base1_22.py
==========================
Generate the four-way exact certificate for every one of the 22 odd-rank
relabelings of base1, namely:

  (1) exact integer rank of E_{gamma, base1} = gamma(base1) - 5*J via Sympy,
  (2) a primitive integer kernel vector u in V_std, sign-normalized,
  (3) the integer Smith Normal Form (SNF) elementary-divisor list,
  (4) the F_2 Smith Normal Form elementary-divisor list (rank distribution
      modulo 2).

This is the explicit data for Appendix A of this paper.

Output: ``four_way_base1_22.json`` next to this script.

Runtime: a few seconds.
"""
from __future__ import annotations

import json
from math import gcd
from pathlib import Path
from typing import List

import numpy as np
from sympy import Matrix, Rational
from sympy.matrices.normalforms import smith_normal_form

HERE = Path(__file__).resolve().parent

BASE1 = np.array(
    [
        [5, 3, 4, 6, 7, 8, 9, 1, 2],
        [6, 7, 2, 1, 9, 5, 3, 4, 8],
        [1, 9, 8, 3, 4, 2, 5, 6, 7],
        [8, 5, 9, 7, 6, 1, 4, 2, 3],
        [4, 2, 6, 8, 5, 3, 7, 9, 1],
        [7, 1, 3, 9, 2, 4, 8, 5, 6],
        [9, 6, 1, 5, 3, 7, 2, 8, 4],
        [2, 8, 7, 4, 1, 9, 6, 3, 5],
        [3, 4, 5, 2, 8, 6, 1, 7, 9],
    ],
    dtype=np.int64,
)


def primitive_kernel(M: Matrix) -> List[int]:
    """Return a primitive (gcd 1, sign-normalized) integer kernel vector in
    V_std (sum = 0) of an integer 9x9 matrix M of rank 7.  The full kernel
    is 2-dimensional (the all-ones vector is always in ker E for centered
    Sudoku operators); the intersection ker(E) cap V_std is 1-dimensional,
    obtained as the nullspace of [E; 1^T]."""
    A = M.col_join(Matrix([[1] * M.cols]))
    null = A.nullspace()
    assert len(null) == 1, f"expected 1-dim kernel in V_std, got {len(null)}"
    v = null[0]
    # Clear denominators.
    denom = 1
    for x in v:
        denom = denom * Rational(x).q // gcd(denom, Rational(x).q)
    vint = [int(Rational(x) * denom) for x in v]
    # Reduce by gcd.
    g = 0
    for x in vint:
        g = gcd(g, abs(x))
    if g > 1:
        vint = [x // g for x in vint]
    # Sign-normalize: first nonzero entry positive.
    for x in vint:
        if x != 0:
            if x < 0:
                vint = [-y for y in vint]
            break
    return vint


def smith_diagonal_int(M: Matrix) -> List[int]:
    """Integer SNF elementary divisors as a list of ints."""
    S = smith_normal_form(M)
    n, m = S.shape
    out = []
    for i in range(min(n, m)):
        d = int(S[i, i])
        out.append(d)
    return out


def smith_diagonal_F2(M: np.ndarray) -> List[int]:
    """Smith-style invariants over F_2: every nonzero invariant is 1, the
    rank over F_2 is the count of 1's, and the rest are 0.  Returns a list
    of length n = nrows = ncols (assumed square)."""
    A = M.copy() % 2
    n, m = A.shape
    rank_F2 = 0
    r = 0
    for c in range(m):
        # Find pivot
        piv = None
        for k in range(r, n):
            if A[k, c] == 1:
                piv = k
                break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        for k in range(n):
            if k != r and A[k, c] == 1:
                A[k] = (A[k] + A[r]) % 2
        rank_F2 += 1
        r += 1
        if r == n:
            break
    # Invariants: rank_F2 ones, then zeros.
    out = [1] * rank_F2 + [0] * (n - rank_F2)
    return out


def relabel(L: np.ndarray, gamma: List[int]) -> np.ndarray:
    g = np.asarray(gamma, dtype=np.int64)
    return g[L - 1]


def main() -> int:
    rec_path = HERE / "base1_22_perms_recovered.json"
    with rec_path.open("r", encoding="utf-8") as f:
        rec = json.load(f)
    perms: List[List[int]] = rec["odd_perms"]
    assert len(perms) == 22, f"expected 22, got {len(perms)}"

    certs = []
    for perm in perms:
        Lg = relabel(BASE1, perm)
        E = (Lg - 5).tolist()
        Mtx = Matrix(E)
        rk = int(Mtx.rank())
        assert rk == 7, f"perm {perm}: rank {rk} != 7"
        u = primitive_kernel(Mtx)
        snf_Z = smith_diagonal_int(Mtx)
        snf_F2 = smith_diagonal_F2(np.array(E, dtype=np.int64))
        rank_F2 = sum(1 for d in snf_F2 if d == 1)
        certs.append(
            {
                "perm": list(perm),
                "exact_rank_Z": rk,
                "primitive_kernel_vector_V_std": u,
                "kernel_vector_sum_zero": sum(u) == 0,
                "snf_Z_diagonal": snf_Z,
                "rank_Z_from_snf": sum(1 for d in snf_Z if d != 0),
                "snf_F2_diagonal": snf_F2,
                "rank_F2": rank_F2,
            }
        )

    out = {
        "phase": "App A four-way exact certificate for 22 base1 odd-rank perms",
        "base_grid": BASE1.tolist(),
        "n_certs": len(certs),
        "four_way_certificates": certs,
        "method": (
            "Sympy exact integer rank + Sympy nullspace primitive kernel + "
            "Sympy smith_normal_form (Z) + custom F_2 row-reduction (F_2 SNF)."
        ),
    }
    out_path = HERE / "four_way_base1_22.json"
    with out_path.open("w", encoding="utf-8") as f:
        json.dump(out, f, sort_keys=True, indent=2)
    print(f"[OK] wrote {out_path.name}  ({len(certs)} four-way certs)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
