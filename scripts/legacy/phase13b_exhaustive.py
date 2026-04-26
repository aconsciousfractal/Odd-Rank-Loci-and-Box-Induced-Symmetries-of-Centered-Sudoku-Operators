#!/usr/bin/env python3
"""
Phase 13b — Exhaustive Verification of D13.3 + D13.4

Runs on ALL 711 FR and ALL 1473 NM (rank-9, min_ov ≥ 1) from Phase 12 corpus.

D13.3 exhaustive: For every FR pattern, find at least one degrading switch.
                   For every NM pattern, find at least one upgrading switch.
                   Compute exact degradation rate per FR pattern.

D13.4 exhaustive: Compute SNF for all 711 FR + all 1473 NM.
                   Verify complete disjunction.
"""

import numpy as np
from collections import Counter
from fractions import Fraction
import json, time, os, sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, '..', 'results', 'phase12')
N, K = 10, 5

# ═══════════════════════════════════════════════════════════════════
# UTILITIES
# ═══════════════════════════════════════════════════════════════════

def rank_f2(P):
    M = P.copy() % 2
    n_rows, n_cols = M.shape
    r = 0
    for c in range(n_cols):
        piv = None
        for row in range(r, n_rows):
            if M[row, c] % 2:
                piv = row; break
        if piv is None: continue
        M[[r, piv]] = M[[piv, r]]
        pv = M[r].copy()
        for row in range(n_rows):
            if row != r and M[row, c] % 2:
                M[row] = (M[row] + pv) % 2
        r += 1
    return r

def snf(M_int):
    """Smith Normal Form. Returns sorted list of invariant factors (nonzero only)."""
    n = len(M_int)
    m = len(M_int[0])
    A = [[int(M_int[i][j]) for j in range(m)] for i in range(n)]

    size = min(n, m)
    for s in range(size):
        found = True
        while found:
            found = False
            piv_r, piv_c, piv_v = -1, -1, 0
            for i in range(s, n):
                for j in range(s, m):
                    if A[i][j] != 0:
                        if piv_v == 0 or abs(A[i][j]) < abs(piv_v):
                            piv_r, piv_c, piv_v = i, j, A[i][j]
            if piv_v == 0:
                break
            if piv_r != s:
                A[s], A[piv_r] = A[piv_r], A[s]
            if piv_c != s:
                for i in range(n):
                    A[i][s], A[i][piv_c] = A[i][piv_c], A[i][s]
            for i in range(s+1, n):
                if A[i][s] != 0:
                    q = A[i][s] // A[s][s]
                    for j in range(m):
                        A[i][j] -= q * A[s][j]
                    if A[i][s] != 0:
                        found = True
            for j in range(s+1, m):
                if A[s][j] != 0:
                    q = A[s][j] // A[s][s]
                    for i in range(n):
                        A[i][j] -= q * A[i][s]
                    if A[s][j] != 0:
                        found = True
            if not found:
                for i in range(s+1, n):
                    for j in range(s+1, m):
                        if A[i][j] % A[s][s] != 0:
                            for jj in range(m):
                                A[s][jj] += A[i][jj]
                            found = True
                            break
                    if found:
                        break
        if A[s][s] < 0:
            for j in range(m):
                A[s][j] = -A[s][j]

    diag = [abs(A[i][i]) if i < min(n,m) else 0 for i in range(min(n,m))]
    return tuple(d for d in sorted(diag) if d != 0)

# ═══════════════════════════════════════════════════════════════════
# LOAD DATA
# ═══════════════════════════════════════════════════════════════════

def load_corpus():
    with open(os.path.join(DATA_DIR, 'full_rank_corpus.json')) as f:
        fr = json.load(f)
    with open(os.path.join(DATA_DIR, 'deficient_controls.json')) as f:
        de = json.load(f)
    return fr, de

# ═══════════════════════════════════════════════════════════════════
# D13.3 EXHAUSTIVE: 1-switch thickness
# ═══════════════════════════════════════════════════════════════════

def exhaustive_switch_check(records, target_rank_change, label):
    """For each pattern, enumerate ALL valid elementary switches.
    Count how many change the rank in the desired direction.
    target_rank_change: 'upgrade' (rank < N → N) or 'degrade' (rank N → < N)
    Returns list of (n_valid_switches, n_changing_switches) per pattern.
    """
    n = N
    results = []
    n_fail = 0  # patterns where NO switch achieves the change

    for idx, rec in enumerate(records):
        P_orig = np.array(rec['matrix'], dtype=np.int8)
        n_valid = 0
        n_change = 0

        for i in range(n):
            for j in range(i+1, n):
                for c in range(n):
                    for d in range(c+1, n):
                        P = P_orig.copy()
                        a, b, e, f = P[i,c], P[i,d], P[j,c], P[j,d]
                        did = False
                        if a == 1 and f == 1 and b == 0 and e == 0:
                            P[i,c], P[i,d], P[j,c], P[j,d] = 0, 1, 1, 0
                            did = True
                        elif b == 1 and e == 1 and a == 0 and f == 0:
                            P[i,c], P[i,d], P[j,c], P[j,d] = 1, 0, 0, 1
                            did = True
                        if did:
                            n_valid += 1
                            rk = rank_f2(P)
                            if target_rank_change == 'upgrade' and rk == N:
                                n_change += 1
                            elif target_rank_change == 'degrade' and rk < N:
                                n_change += 1

        results.append((n_valid, n_change))
        if n_change == 0:
            n_fail += 1

        done = idx + 1
        if done % 50 == 0 or done == len(records):
            pct = 100 * (len(records) - n_fail) / done
            print(f"\r    {label}: {done}/{len(records)}  "
                  f"has_change={done - n_fail}/{done} ({pct:.1f}%)  ", end='', flush=True)

    print()
    return results

def section_d13_3(full_rank, near_miss):
    print("\n" + "=" * 72)
    print("  D13.3 EXHAUSTIVE — 1-switch thickness")
    print("=" * 72)

    # ── FR → degrade (all 711) ──
    print(f"\n  FR → rank < 10 (all {len(full_rank)} patterns, exhaustive 1-switch):")
    fr_results = exhaustive_switch_check(full_rank, 'degrade', 'FR→degrade')

    n_has_deg = sum(1 for vs, nc in fr_results if nc > 0)
    rates = [nc/vs if vs > 0 else 0 for vs, nc in fr_results]
    n_valid_arr = [vs for vs, nc in fr_results]
    print(f"  Patterns with ≥1 degrading switch: {n_has_deg}/{len(full_rank)} "
          f"({100*n_has_deg/len(full_rank):.2f}%)")
    print(f"  Degradation rate: mean={np.mean(rates):.4f}, "
          f"std={np.std(rates):.4f}, min={min(rates):.4f}, max={max(rates):.4f}")
    print(f"  Valid switches: mean={np.mean(n_valid_arr):.1f}")

    if n_has_deg == len(full_rank):
        print(f"  ✅ UNIVERSAL: every FR pattern has at least one degrading switch")
    else:
        fails = [i for i, (vs, nc) in enumerate(fr_results) if nc == 0]
        print(f"  ⚠️ {len(fails)} FR patterns have NO degrading switch: indices {fails[:10]}")

    # ── NM → upgrade (all 1473) ──
    print(f"\n  NM → rank 10 (all {len(near_miss)} patterns, exhaustive 1-switch):")
    nm_results = exhaustive_switch_check(near_miss, 'upgrade', 'NM→upgrade')

    n_has_upg = sum(1 for vs, nc in nm_results if nc > 0)
    up_rates = [nc/vs if vs > 0 else 0 for vs, nc in nm_results]
    print(f"  Patterns with ≥1 upgrading switch: {n_has_upg}/{len(near_miss)} "
          f"({100*n_has_upg/len(near_miss):.2f}%)")
    print(f"  Upgrade rate: mean={np.mean(up_rates):.4f}, "
          f"std={np.std(up_rates):.4f}, min={min(up_rates):.4f}, max={max(up_rates):.4f}")

    if n_has_upg == len(near_miss):
        print(f"  ✅ UNIVERSAL: every NM pattern has at least one upgrading switch")
    else:
        fails = [i for i, (vs, nc) in enumerate(nm_results) if nc == 0]
        print(f"  ⚠️ {len(fails)} NM patterns have NO upgrading switch: indices {fails[:10]}")

    return fr_results, nm_results

# ═══════════════════════════════════════════════════════════════════
# D13.4 EXHAUSTIVE: SNF complete census
# ═══════════════════════════════════════════════════════════════════

def section_d13_4(full_rank, near_miss):
    print("\n" + "=" * 72)
    print("  D13.4 EXHAUSTIVE — SNF complete census")
    print("=" * 72)

    # ── FR: all 711 ──
    print(f"\n  Computing SNF for all {len(full_rank)} FR patterns...")
    fr_snf = []
    for idx, rec in enumerate(full_rank):
        P = np.array(rec['matrix'], dtype=np.int8)
        s = snf(P.tolist())
        fr_snf.append(s)
        if (idx+1) % 100 == 0 or idx+1 == len(full_rank):
            print(f"\r    FR SNF: {idx+1}/{len(full_rank)}", end='', flush=True)
    print()

    fr_snf_dist = Counter(fr_snf)
    print(f"\n  FR SNF types (all {len(full_rank)}):")
    for s, c in fr_snf_dist.most_common():
        print(f"    {s}: {c}")
    fr_snf_set = set(fr_snf_dist.keys())

    # ── NM: all 1473 ──
    print(f"\n  Computing SNF for all {len(near_miss)} NM patterns...")
    nm_snf = []
    for idx, rec in enumerate(near_miss):
        P = np.array(rec['matrix'], dtype=np.int8)
        s = snf(P.tolist())
        nm_snf.append(s)
        if (idx+1) % 200 == 0 or idx+1 == len(near_miss):
            print(f"\r    NM SNF: {idx+1}/{len(near_miss)}", end='', flush=True)
    print()

    nm_snf_dist = Counter(nm_snf)
    print(f"\n  NM SNF types (all {len(near_miss)}):")
    for s, c in nm_snf_dist.most_common():
        print(f"    {s}: {c}")
    nm_snf_set = set(nm_snf_dist.keys())

    # ── Overlap check ──
    overlap = fr_snf_set & nm_snf_set
    print(f"\n  SNF types in FR: {len(fr_snf_set)}")
    print(f"  SNF types in NM: {len(nm_snf_set)}")
    print(f"  Overlap (shared types): {len(overlap)}")
    if overlap:
        print(f"  ⚠️ SHARED SNF TYPES:")
        for s in overlap:
            print(f"    {s}: FR={fr_snf_dist[s]}, NM={nm_snf_dist[s]}")
    else:
        print(f"  ✅ COMPLETE DISJUNCTION: 0 shared SNF types between FR and NM")
        print(f"     (all {len(full_rank)} FR + all {len(near_miss)} NM verified)")

    return fr_snf_dist, nm_snf_dist

# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════

def main():
    print("=" * 72)
    print("  PHASE 13b — EXHAUSTIVE VERIFICATION (D13.3 + D13.4)")
    print("=" * 72)

    t0 = time.time()

    full_rank, deficient = load_corpus()
    near_miss = [r for r in deficient
                 if r['features']['rank'] == 9 and r['features']['min_ov'] >= 1]
    print(f"  FR={len(full_rank)}, NM={len(near_miss)}")

    # D13.4 first (faster — no switch enumeration)
    fr_snf_dist, nm_snf_dist = section_d13_4(full_rank, near_miss)

    # D13.3 exhaustive (slower — switch enumeration for all patterns)
    fr_sw, nm_sw = section_d13_3(full_rank, near_miss)

    elapsed = time.time() - t0
    print(f"\n{'=' * 72}")
    print(f"  PHASE 13b EXHAUSTIVE COMPLETE — {elapsed:.1f}s total")
    print(f"{'=' * 72}")

if __name__ == '__main__':
    main()
