"""Phase 9.18.S9c — F_p-rank bridge for the 22 odd-rank gammas at n=9.

For each of the 22 gamma producing rk_Q(E_gamma) = 7 on base1, compute:
  rk_F_2(E_gamma)
  rk_F_3(E_gamma)
and correlate with the SNF signature (already computed in 9.13.S9 results).

Identity to verify (CORRECTED — external review 2026-04-25):
  rk_F_p(E_gamma) = #{ d_i in SNF(E_gamma) : d_i != 0 and p does not divide d_i }
                  = rk_Q(E_gamma) - #{ d_i nonzero : p | d_i }
                  = 7 - #{ d_i in SNF : p | d_i }       (these matrices are 8x8, rk_Q = 7)

So:
  - SNF with no even factor => rk_F_2 = 7 (max possible at rk_Q = 7), NOT 8.
  - SNF with one even factor => rk_F_2 = 6.
  - SNF with two even factors => rk_F_2 = 5.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np


def rank_mod_p(M, p):
    """Gaussian elimination over GF(p)."""
    A = [[int(x) % p for x in row] for row in M]
    rows = len(A)
    cols = len(A[0]) if rows else 0
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
        # modular inverse
        a = A[r][c]
        # extended euclid
        def inv(a, p):
            t, newt = 0, 1
            r0, newr = p, a
            while newr != 0:
                q = r0 // newr
                t, newt = newt, t - q * newt
                r0, newr = newr, r0 - q * newr
            return t % p
        ai = inv(a, p)
        # normalize row r
        A[r] = [(x * ai) % p for x in A[r]]
        for k in range(rows):
            if k == r or A[k][c] == 0:
                continue
            f = A[k][c]
            A[k] = [(A[k][j] - f * A[r][j]) % p for j in range(cols)]
        r += 1
    return r


def build_S(L_base, perm):
    L_base = np.asarray(L_base, dtype=int)
    gc = np.array(perm, dtype=int) - 5
    E = gc[L_base - 1]
    last = E[8, 8]
    S = E[:8, :8] - E[:8, 8:9] - E[8:9, :8] + last
    return S


def main():
    here = Path(__file__).resolve().parents[2]  # <repo>/scripts/legacy/ -> <repo>/
    with open(here / "data" / "phase_9_13_S9_results.json", "r", encoding="utf-8") as f:
        d = json.load(f)
    L_base = np.array(d["base_grid"], dtype=int)
    elems = d["elements"]

    rows = []
    for e in elems:
        S = build_S(L_base, e["perm"])
        Sl = S.tolist()
        r2 = rank_mod_p(Sl, 2)
        r3 = rank_mod_p(Sl, 3)
        rows.append({
            "perm": e["perm"],
            "rk_Q": e["rank"],
            "snf": e["snf_signature"],
            "rk_F2": r2,
            "rk_F3": r3,
            "snf_has_factor_2": any(x % 2 == 0 for x in e["snf"]),
            "snf_has_factor_3": any(x % 3 == 0 for x in e["snf"]),
        })

    print("=== Phase 9.18.S9c — F_p-rank bridge for the 22 odd gammas ===\n")
    print(f"{'perm':<35s} {'snf':<25s} {'rkQ':>3s} {'rkF2':>5s} {'rkF3':>5s}")
    for r in rows:
        print(f"{str(r['perm']):<35s} {r['snf']:<25s} "
              f"{r['rk_Q']:>3d} {r['rk_F2']:>5d} {r['rk_F3']:>5d}")

    # Distributions
    print("\n--- Distributions ---")
    print(f"  rk_F2 multiset: {dict(Counter(r['rk_F2'] for r in rows))}")
    print(f"  rk_F3 multiset: {dict(Counter(r['rk_F3'] for r in rows))}")

    print("\n--- Cross-tab: SNF-has-factor-2 vs rk_F2 ---")
    pairs2 = Counter((r["snf_has_factor_2"], r["rk_F2"]) for r in rows)
    for k, v in sorted(pairs2.items()):
        print(f"  has2={k[0]!s:<5s} rk_F2={k[1]} : {v}")

    print("\n--- Cross-tab: SNF-has-factor-3 vs rk_F3 ---")
    pairs3 = Counter((r["snf_has_factor_3"], r["rk_F3"]) for r in rows)
    for k, v in sorted(pairs3.items()):
        print(f"  has3={k[0]!s:<5s} rk_F3={k[1]} : {v}")

    out = {
        "phase": "9.18.S9c",
        "n_elements": len(rows),
        "elements": rows,
        "rk_F2_distribution": {str(k): v for k, v in Counter(r["rk_F2"] for r in rows).items()},
        "rk_F3_distribution": {str(k): v for k, v in Counter(r["rk_F3"] for r in rows).items()},
    }
    dst = here / "data" / "phase_9_18_S9c_results.json"
    with open(dst, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved -> {dst}")


if __name__ == "__main__":
    main()
