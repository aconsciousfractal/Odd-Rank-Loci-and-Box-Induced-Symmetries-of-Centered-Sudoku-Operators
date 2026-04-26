"""
PHASE 9.15 — Algebraic Characterization of Odd-Rank Relabelings at n=7
========================================================================

Open problem (b) of the Parity paper: at n=7 over the non-cyclic base
L_V1, exactly 30/5040 = 0.595% of symbol relabelings gamma in S_7 give
rank_F2(E_std) = 5 (odd). What is the algebraic structure of this set?

Tasks:
  9.15a  Extract the 30 odd-rank gammas + invariants
         (cycle type, parity, fix(m=4), action on complementary pairs)
  9.15b  Test coset/subgroup structure (left/right cosets in S_7,
         orbit under stabiliser of m=4)
  9.15c  Repeat over a second base L_V3 (= L_7 of the paper, std-label CE)
  9.15d  Extend to n=9 base L_Sud (sample, not exhaustive: 9!=362880)
  9.15e  Formulate algebraic necessary condition

Outputs:
  data/phase_9_15_results.json
  data/phase_9_15_log.txt
"""
from __future__ import annotations

import itertools
import json
import time
from fractions import Fraction
from pathlib import Path

# --------------------------------------------------------------------- #
# Bases                                                                 #
# --------------------------------------------------------------------- #

L_V1 = [
    [3, 5, 2, 7, 1, 4, 6],
    [6, 3, 1, 5, 7, 2, 4],
    [1, 6, 7, 4, 2, 3, 5],
    [7, 4, 3, 2, 5, 6, 1],
    [2, 1, 4, 6, 3, 5, 7],
    [5, 2, 6, 1, 4, 7, 3],
    [4, 7, 5, 3, 6, 1, 2],
]

# Standard-labeling counterexample at n=7 (paper's L_7).
L_V3 = [
    [3, 7, 4, 6, 5, 2, 1],
    [5, 6, 7, 2, 1, 3, 4],
    [2, 5, 1, 7, 6, 4, 3],
    [7, 4, 5, 1, 3, 6, 2],
    [6, 1, 2, 3, 4, 7, 5],
    [4, 3, 6, 5, 2, 1, 7],
    [1, 2, 3, 4, 7, 5, 6],
]

L_SUDOKU_9 = [
    [5, 3, 4, 6, 7, 8, 9, 1, 2],
    [6, 7, 2, 1, 9, 5, 3, 4, 8],
    [1, 9, 8, 3, 4, 2, 5, 6, 7],
    [8, 5, 9, 7, 6, 1, 4, 2, 3],
    [4, 2, 6, 8, 5, 3, 7, 9, 1],
    [7, 1, 3, 9, 2, 4, 8, 5, 6],
    [9, 6, 1, 5, 3, 7, 2, 8, 4],
    [2, 8, 7, 4, 1, 9, 6, 3, 5],
    [3, 4, 5, 2, 8, 6, 1, 7, 9],
]


# --------------------------------------------------------------------- #
# Exact rank for E_std                                                  #
# --------------------------------------------------------------------- #

def build_E_std(L, n):
    m = (n + 1) // 2
    d = n - 1
    E = [[L[i][j] - m for j in range(n)] for i in range(n)]
    out = [[0] * d for _ in range(d)]
    for i in range(d):
        for j in range(d):
            out[i][j] = E[i][j] - E[n - 1][j] - E[i][n - 1] + E[n - 1][n - 1]
    return out


def rref_rank_Q(M):
    n = len(M)
    if n == 0:
        return 0
    m = len(M[0])
    A = [[Fraction(x) for x in row] for row in M]
    rank = 0
    for col in range(m):
        pivot_row = None
        for row in range(rank, n):
            if A[row][col] != 0:
                pivot_row = row
                break
        if pivot_row is None:
            continue
        A[rank], A[pivot_row] = A[pivot_row], A[rank]
        scale = A[rank][col]
        for j in range(m):
            A[rank][j] /= scale
        for i in range(n):
            if i == rank:
                continue
            factor = A[i][col]
            if factor != 0:
                for j in range(m):
                    A[i][j] -= factor * A[rank][j]
        rank += 1
    return rank


def rank_relabel(L, gamma, n):
    """gamma is a list/tuple of length n: gamma[k-1] = image of symbol k."""
    Lg = [[gamma[L[i][j] - 1] for j in range(n)] for i in range(n)]
    return rref_rank_Q(build_E_std(Lg, n))


# --------------------------------------------------------------------- #
# Permutation invariants                                                #
# --------------------------------------------------------------------- #

def cycle_type(perm):
    """perm: tuple of length n with values in 1..n. Returns sorted tuple of cycle lengths."""
    n = len(perm)
    visited = [False] * n
    cycles = []
    for i in range(n):
        if visited[i]:
            continue
        j = i
        L = 0
        while not visited[j]:
            visited[j] = True
            j = perm[j] - 1
            L += 1
        cycles.append(L)
    return tuple(sorted(cycles, reverse=True))


def perm_parity(perm):
    """Sign of perm: 0 = even, 1 = odd (number of transpositions mod 2)."""
    n = len(perm)
    visited = [False] * n
    inv_count = 0
    for i in range(n):
        if visited[i]:
            continue
        j = i
        L = 0
        while not visited[j]:
            visited[j] = True
            j = perm[j] - 1
            L += 1
        inv_count += (L - 1)
    return inv_count % 2


def fixes_center(perm, n):
    m = (n + 1) // 2
    return perm[m - 1] == m


def preserves_complementary_pairs(perm, n):
    """True iff for all k, gamma(n+1-k) = n+1-gamma(k)
       (i.e. gamma commutes with the involution k -> n+1-k)."""
    for k in range(1, n + 1):
        if perm[n - k] != n + 1 - perm[k - 1]:
            return False
    return True


def swaps_complementary_pairs(perm, n):
    """How many of the (n-1)/2 complementary pairs {k,n+1-k} (k<m) are
       sent to a complementary pair {gamma(k), gamma(n+1-k)} which is
       also of the form {j, n+1-j}?"""
    m = (n + 1) // 2
    cnt = 0
    for k in range(1, m):
        a = perm[k - 1]
        b = perm[n - k]
        if a + b == n + 1:
            cnt += 1
    return cnt


def perm_compose(a, b):
    """(a o b)(k) = a(b(k))."""
    n = len(a)
    return tuple(a[b[k] - 1] for k in range(n))


def perm_inverse(a):
    n = len(a)
    out = [0] * n
    for k in range(n):
        out[a[k] - 1] = k + 1
    return tuple(out)


# --------------------------------------------------------------------- #
# 9.15a + 9.15b: scan + invariants                                      #
# --------------------------------------------------------------------- #

def scan_relabelings_n7(L, label):
    print(f"\n== Scanning all 7! = 5040 relabelings over base {label} ==")
    t0 = time.time()
    odd = []
    rank_dist = {}
    for perm in itertools.permutations(range(1, 8)):
        r = rank_relabel(L, perm, 7)
        rank_dist[r] = rank_dist.get(r, 0) + 1
        if r % 2 == 1:
            odd.append(perm)
    print(f"  rank_dist = {rank_dist}")
    print(f"  odd-rank count = {len(odd)}")
    print(f"  elapsed = {time.time() - t0:.1f}s")
    return odd, rank_dist


def catalog_perms(perms, n, label):
    print(f"\n--- Invariants of {len(perms)} odd-rank perms over {label} ---")
    cycle_types = {}
    parities = {0: 0, 1: 0}
    fix_center = 0
    pair_preserve = 0
    pair_swap_counts = {}

    for p in perms:
        ct = cycle_type(p)
        cycle_types[ct] = cycle_types.get(ct, 0) + 1
        parities[perm_parity(p)] += 1
        if fixes_center(p, n):
            fix_center += 1
        if preserves_complementary_pairs(p, n):
            pair_preserve += 1
        sc = swaps_complementary_pairs(p, n)
        pair_swap_counts[sc] = pair_swap_counts.get(sc, 0) + 1

    print(f"  cycle types: {cycle_types}")
    print(f"  parity: even={parities[0]}, odd={parities[1]}")
    print(f"  fix center m={(n + 1) // 2}: {fix_center}/{len(perms)}")
    print(f"  preserve complementary pairs: {pair_preserve}/{len(perms)}")
    print(f"  pairs sent to pairs (count distribution): {pair_swap_counts}")
    return {
        "count": len(perms),
        "cycle_types": {str(k): v for k, v in cycle_types.items()},
        "parity_even": parities[0],
        "parity_odd": parities[1],
        "fix_center": fix_center,
        "preserve_complementary_pairs": pair_preserve,
        "pair_swap_distribution": pair_swap_counts,
    }


# --------------------------------------------------------------------- #
# 9.15b: coset / subgroup tests                                         #
# --------------------------------------------------------------------- #

def is_subgroup(S, n):
    """Return True if S (set of perm tuples) is a subgroup of S_n."""
    e = tuple(range(1, n + 1))
    if e not in S:
        return False
    for a in S:
        if perm_inverse(a) not in S:
            return False
    for a in S:
        for b in S:
            if perm_compose(a, b) not in S:
                return False
    return True


def left_coset_structure(S, n):
    """If S is a left coset gH for some subgroup H, then S^{-1} S = H.
       Return (is_coset, H_or_None, repr_g_or_None)."""
    if not S:
        return False, None, None
    g = next(iter(S))
    g_inv = perm_inverse(g)
    H = {perm_compose(g_inv, h) for h in S}
    if is_subgroup(H, n) and len(H) == len(S):
        return True, H, g
    return False, None, None


def union_of_cosets_structure(S, n):
    """Try all subgroups generated by subsets of S of size <= 4 to find
       maximal subgroup H such that S is a union of cosets of H.
       Heuristic; return largest found."""
    e = tuple(range(1, n + 1))
    # The cosets gH partition S iff: for all g in S, gH subset S; and S = union.
    # Equivalently: H subset S^{-1} S, and S = (union of left cosets).
    # We compute: H_max = largest subgroup such that gH subset S for all g in S.
    # Equivalently, h in H iff for all g in S, g*h in S.
    candidates = []
    S_set = set(S)
    for h in itertools.permutations(range(1, n + 1)):
        ok = True
        for g in S:
            if perm_compose(g, h) not in S_set:
                ok = False
                break
        if ok:
            candidates.append(h)
    H = set(candidates)
    is_subg = is_subgroup(H, n)
    return {
        "right_stabiliser_size": len(H),
        "is_subgroup": is_subg,
        "elements_first_few": [list(h) for h in list(H)[:5]],
    }


# --------------------------------------------------------------------- #
# 9.15d: scan over n=9 (sample only, 9! = 362880)                        #
# --------------------------------------------------------------------- #

def scan_n9_full(L, label):
    """Exhaustive 9! = 362880 scan — about 3-5 minutes."""
    print(f"\n== Exhaustive scan all 9! = 362880 relabelings over base {label} ==")
    t0 = time.time()
    odd = []
    rank_dist = {}
    cnt = 0
    for perm in itertools.permutations(range(1, 10)):
        r = rank_relabel(L, perm, 9)
        rank_dist[r] = rank_dist.get(r, 0) + 1
        if r % 2 == 1:
            odd.append(perm)
        cnt += 1
        if cnt % 20000 == 0:
            print(f"  scanned {cnt:,d}/362880  elapsed={time.time() - t0:.1f}s "
                  f"odd-rank-so-far={len(odd)}")
    print(f"  rank_dist = {rank_dist}")
    print(f"  odd-rank count = {len(odd)}")
    print(f"  elapsed = {time.time() - t0:.1f}s")
    return odd, rank_dist


# --------------------------------------------------------------------- #
# Main                                                                  #
# --------------------------------------------------------------------- #

def main():
    import sys
    do_n9 = "--n9" in sys.argv
    print("=" * 72)
    print("  PHASE 9.15 — Algebraic characterization of odd-rank relabelings")
    print(f"  (n=9 exhaustive scan: {'ENABLED' if do_n9 else 'disabled, pass --n9 to enable'})")
    print("=" * 72)

    out = {}

    # -- 9.15a + b on L_V1 ----------------------------------------------
    odd_V1, dist_V1 = scan_relabelings_n7(L_V1, "L_V1 (n=7, non-cyclic)")
    cat_V1 = catalog_perms(odd_V1, 7, "L_V1")
    cos_V1 = union_of_cosets_structure(odd_V1, 7)
    print(f"  union-of-right-cosets analysis: {cos_V1}")
    out["L_V1"] = {
        "rank_dist": dist_V1,
        "odd_perms": [list(p) for p in odd_V1],
        "catalog": cat_V1,
        "coset_structure": cos_V1,
    }

    # -- 9.15c on L_V3 (paper's standard-labeling CE) -------------------
    odd_V3, dist_V3 = scan_relabelings_n7(L_V3, "L_V3 (n=7, std-label CE)")
    cat_V3 = catalog_perms(odd_V3, 7, "L_V3")
    cos_V3 = union_of_cosets_structure(odd_V3, 7)
    print(f"  union-of-right-cosets analysis: {cos_V3}")
    out["L_V3"] = {
        "rank_dist": dist_V3,
        "odd_perms": [list(p) for p in odd_V3],
        "catalog": cat_V3,
        "coset_structure": cos_V3,
    }

    # -- Compare V1 vs V3 odd sets --------------------------------------
    s1 = {tuple(p) for p in odd_V1}
    s3 = {tuple(p) for p in odd_V3}
    inter = s1 & s3
    print("\n--- L_V1 vs L_V3 odd-perm comparison ---")
    print(f"  |odd(V1)| = {len(s1)},  |odd(V3)| = {len(s3)},  |intersection| = {len(inter)}")
    out["V1_vs_V3"] = {
        "odd_V1": len(s1),
        "odd_V3": len(s3),
        "intersection": len(inter),
        "intersection_perms": [list(p) for p in inter],
    }

    # -- 9.15d on n=9 (exhaustive, optional) -----------------------------
    if do_n9:
        odd_V2, dist_V2 = scan_n9_full(L_SUDOKU_9, "L_Sud (n=9, std-label)")
        if odd_V2:
            cat_V2 = catalog_perms(odd_V2, 9, "L_Sud (n=9)")
            cos_V2 = union_of_cosets_structure(odd_V2, 9) if len(odd_V2) <= 5000 else None
            out["L_Sud_n9"] = {
                "rank_dist": dist_V2,
                "odd_perms_count": len(odd_V2),
                "odd_perms_sample": [list(p) for p in odd_V2[:50]],
                "catalog": cat_V2,
                "coset_structure": cos_V2,
            }
        else:
            out["L_Sud_n9"] = {"rank_dist": dist_V2, "odd_perms_count": 0}
    else:
        out["L_Sud_n9"] = "skipped (pass --n9 to enable)"

    # -- 9.15e: enrichment within B_3 = stabiliser of (17)(26)(35) ------
    print("\n--- Enrichment within complementary-pair stabiliser B_3 (|B_3|=48) ---")
    B3 = [p for p in itertools.permutations(range(1, 8))
          if preserves_complementary_pairs(p, 7)]
    print(f"  |B_3| = {len(B3)} (expected 48)")
    for label, L in [("L_V1", L_V1), ("L_V3", L_V3)]:
        b3_odd = []
        b3_dist = {}
        for p in B3:
            r = rank_relabel(L, p, 7)
            b3_dist[r] = b3_dist.get(r, 0) + 1
            if r % 2 == 1:
                b3_odd.append(p)
        baseline_rate = (30 if label == "L_V1" else 56) / 5040
        b3_rate = len(b3_odd) / len(B3)
        print(f"  {label}: rank dist on B_3 = {b3_dist}; "
              f"odd-rank in B_3 = {len(b3_odd)}/{len(B3)} ({b3_rate:.4f}); "
              f"vs baseline {baseline_rate:.4f}; "
              f"enrichment = {b3_rate / baseline_rate:.2f}x")
        out.setdefault("B3_enrichment", {})[label] = {
            "rank_dist": b3_dist,
            "odd_rank_count": len(b3_odd),
            "rate_in_B3": b3_rate,
            "baseline_rate": baseline_rate,
            "enrichment_ratio": b3_rate / baseline_rate,
            "odd_perms_in_B3": [list(p) for p in b3_odd],
        }

    # -- Are V1 cap V3 perms all in B_3? --------------------------------
    inter_in_B3 = sum(1 for p in inter if preserves_complementary_pairs(p, 7))
    print(f"\n  V1 cap V3 (12 perms) of which in B_3: {inter_in_B3}/{len(inter)}")
    out["V1_cap_V3_in_B3"] = {"in_B3": inter_in_B3, "total": len(inter)}

    # -- Save -----------------------------------------------------------
    out_path = Path(__file__).resolve().parent.parent.parent / "data" / "phase_9_15_results.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\n  results saved -> {out_path}")

    # -- Summary --------------------------------------------------------
    print()
    print("=" * 72)
    print("  SUMMARY")
    print("=" * 72)
    print(f"  L_V1: {len(odd_V1)} odd-rank perms; "
          f"cycle types = {cat_V1['cycle_types']}; "
          f"parity even/odd = {cat_V1['parity_even']}/{cat_V1['parity_odd']}; "
          f"complementary-pair preservers = {cat_V1['preserve_complementary_pairs']}")
    print(f"  L_V3: {len(odd_V3)} odd-rank perms; "
          f"cycle types = {cat_V3['cycle_types']}; "
          f"parity even/odd = {cat_V3['parity_even']}/{cat_V3['parity_odd']}; "
          f"complementary-pair preservers = {cat_V3['preserve_complementary_pairs']}")
    print(f"  V1 cap V3: {len(inter)}")
    if do_n9 and isinstance(out.get("L_Sud_n9"), dict):
        print(f"  L_Sud (n=9): {out['L_Sud_n9']['odd_perms_count']} odd-rank perms over 362880")


if __name__ == "__main__":
    main()
