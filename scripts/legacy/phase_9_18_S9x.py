"""Phase 9.18.S9x — Analytic mechanism of the M19 macro-orbit (mass 72).

Context (S9w / S9w-bis / S9y-six summary):
- M19 is the unique base whose odd-rank locus contains a singleton
  orbit of mass ≥ 4 (mass = 72 perms collapse to ONE certificate-orbit
  in C^min/G_template).
- All other 61 bases either (a) only host recurrent sparse/low orbits
  with small mass, or (b) have purely high-height singletons.

GOAL:
  Characterize **analytically** the 72-element subset
  Γ*_M19 := {γ ∈ Γ_M19 : (simpler_side(γ), u_simpler(γ)) ∈ O_M19} ⊂ S_9.

Plan:
(X1) Recover the canonical orbit representative u9* of O_M19 and its
     stabilizer H_u := Stab_{S_9}(u9*) acting on coordinates.

(X2) For every γ ∈ Γ*_M19, find π ∈ S_9 such that π · u9* = u9(γ).
     This decomposition factors through the coset π H_u.

(X3) Test combinatorial hypotheses on Γ*_M19:
     (H-A) Is Γ*_M19 a union of LEFT cosets of some subgroup K ≤ S_9?
     (H-B) Is Γ*_M19 itself a single LEFT coset π_0 K (then |K|=72)?
     (H-C) Is Γ*_M19 closed under multiplication by Stab(L_M19) (the
           autotopism group of the M19 grid)?
     (H-D) Is the set {π_γ H_u : γ ∈ Γ*_M19} a single H_u-coset?

(X4) Identify the small-orbit residue Γ_M19 \\ Γ*_M19 (mass ≤ 28).

(X5) Side mix: which γ have simpler_side='left' vs 'right' inside Γ*_M19.

Output:
  - data/phase_9_18_S9x_results.json
  - data/phase_9_18_S9x_output.txt
"""
from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from itertools import permutations
from pathlib import Path

import numpy as np
from sympy import Matrix, Rational

HERE = Path(__file__).resolve().parents[2]  # <repo>/scripts/legacy/ -> <repo>/
sys.path.insert(0, str(HERE / "scripts" / "legacy"))

from phase_9_18_S9b import build_S  # noqa: E402
from phase_9_18_S9p import canonicalize_full, BAND_PERMS  # noqa: E402
from phase_9_18_S9t_ter import load_bases  # noqa: E402


# ------------------------------------------------------------
# Kernel helpers
# ------------------------------------------------------------

def primitive_kernel_v8(S8):
    """Return primitive integer left- and right-null vectors of an 8x8
    integer matrix S8 of rank 7, both as length-8 numpy int arrays.
    Returns (u8_left, u8_right) — the unique generators up to ±.
    Sign normalized so first nonzero entry is positive.
    """
    M = Matrix(S8.tolist())
    # right null
    rn = M.nullspace()
    if len(rn) != 1:
        return None, None
    v = rn[0]
    # primitive: clear denominators and gcd
    denoms = [Rational(x).q for x in v]
    from math import lcm, gcd
    L = 1
    for d in denoms:
        L = lcm(L, d)
    vint = [int(Rational(x) * L) for x in v]
    g = 0
    for x in vint:
        g = gcd(g, abs(x))
    if g > 1:
        vint = [x // g for x in vint]
    # sign normalize
    for x in vint:
        if x != 0:
            if x < 0:
                vint = [-y for y in vint]
            break
    u8_right = np.array(vint, dtype=int)

    # left null = right null of S^T
    Mt = Matrix(S8.T.tolist())
    ln = Mt.nullspace()
    if len(ln) != 1:
        return None, None
    v = ln[0]
    denoms = [Rational(x).q for x in v]
    L = 1
    for d in denoms:
        L = lcm(L, d)
    vint = [int(Rational(x) * L) for x in v]
    g = 0
    for x in vint:
        g = gcd(g, abs(x))
    if g > 1:
        vint = [x // g for x in vint]
    for x in vint:
        if x != 0:
            if x < 0:
                vint = [-y for y in vint]
            break
    u8_left = np.array(vint, dtype=int)

    return u8_left, u8_right


def lift_v8_to_v9(u8):
    """V_std isomorphism Z^8 ↔ Z^9 ∩ {sum=0}: (u_0..u_7) → (u_0..u_7, -Σ)."""
    u9 = list(int(x) for x in u8)
    u9.append(-sum(u9))
    return tuple(u9)


def primitive_v9(v9):
    """Make a 9-tuple primitive and sign-normalized (first nonzero > 0)."""
    from math import gcd
    g = 0
    for x in v9:
        g = gcd(g, abs(int(x)))
    if g > 1:
        v9 = tuple(int(x) // g for x in v9)
    for x in v9:
        if x != 0:
            if x < 0:
                v9 = tuple(-int(y) for y in v9)
            break
    return tuple(int(x) for x in v9)


# ------------------------------------------------------------
# Load M19 perms and certificates
# ------------------------------------------------------------

def load_m19_certs():
    """Load M19's 100 odd-rank perms with simpler-side cert vectors."""
    with open(HERE / "data" / "phase_9_18_S9m_results.json", "r", encoding="utf-8") as f:
        s9m = json.load(f)
    info = s9m["per_base"]["M19"]
    rows = info["rows"]
    out = []
    for r in rows:
        if r.get("corank>1"):
            continue
        side = r["simpler_side"]
        vobj = r["left_u9" if side == "left" else "right_v9"]
        out.append({
            "perm": tuple(r["perm"]),
            "side": side,
            "u9_simpler": tuple(vobj["vec"]),
            "left_u9": tuple(r["left_u9"]["vec"]),
            "right_v9": tuple(r["right_v9"]["vec"]),
            "snf": tuple(r["snf"]),
            "simpler_class": tuple(r["simpler_class"]),
        })
    return out


# ------------------------------------------------------------
# S_9 group operations on 9-tuples (treated as functions {1..9}→{1..9})
# Convention: γ is a permutation, γ[i]−1 = image of i (0-based).
# Action on a 9-vector v indexed by coordinates i=0..8:
#   (π · v)[i] = v[π^{-1}(i)]    (permutation of coordinate labels)
# ------------------------------------------------------------

def perm_inv(p):
    n = len(p)
    inv = [0] * n
    for i, v in enumerate(p):
        inv[v - 1] = i + 1
    return tuple(inv)


def perm_compose(p, q):
    """Return p ∘ q: i ↦ p(q(i))."""
    return tuple(p[q[i] - 1] for i in range(len(p)))


def perm_act_on_vec(pi, v):
    """((π · v))[i] = v[π^{-1}(i)−1]; here pi is a 9-tuple in 1-based form."""
    inv = perm_inv(pi)
    return tuple(v[inv[i] - 1] for i in range(9))


def find_pi_mapping(u_target, u_source):
    """Find ALL π ∈ S_9 with π · u_source = u_target (coord permutation)."""
    n = 9
    found = []
    # Group source coordinates by value
    val_to_src = defaultdict(list)
    for i, x in enumerate(u_source):
        val_to_src[x].append(i)
    val_to_tgt = defaultdict(list)
    for i, x in enumerate(u_target):
        val_to_tgt[x].append(i)
    if sorted(val_to_src.keys()) != sorted(val_to_tgt.keys()):
        return []
    if any(len(val_to_src[k]) != len(val_to_tgt[k]) for k in val_to_src):
        return []

    # Build permutations: π(src_index) = tgt_index, for each value class
    from itertools import permutations as iperm
    # π is given as 1-based mapping: π[i] = j means i+1 maps to j+1 (i.e. π(i+1)=j+1)
    # We define π such that (π · u_src)[π(i+1)−1] = u_src[i] = u_tgt[π(i+1)−1]
    # That is: π sends coordinates of u_src to coordinates of u_tgt with same value.
    keys = list(val_to_src.keys())
    src_lists = [val_to_src[k] for k in keys]
    tgt_lists = [val_to_tgt[k] for k in keys]

    def gen(idx, mapping):
        if idx == len(keys):
            # mapping: 0-based src_pos → 0-based tgt_pos
            pi_arr = [0] * n
            # As 1-based: π(src+1) = tgt+1
            for s, t in mapping.items():
                pi_arr[s] = t + 1
            found.append(tuple(pi_arr))
            return
        for tgt_perm in iperm(tgt_lists[idx]):
            new = dict(mapping)
            for s, t in zip(src_lists[idx], tgt_perm):
                new[s] = t
            gen(idx + 1, new)

    gen(0, {})
    # Sanity: verify
    valid = []
    for pi in found:
        if perm_act_on_vec(pi, u_source) == tuple(u_target):
            valid.append(pi)
    return valid


def stabilizer_of_vec(v):
    """All π ∈ S_9 fixing v under coordinate permutation."""
    return find_pi_mapping(list(v), list(v))


# ------------------------------------------------------------
# Subgroup tests
# ------------------------------------------------------------

def is_subgroup(elts):
    """Test if a set of 9-tuples is closed under composition and inverse."""
    eltset = set(elts)
    n = len(elts)
    for p in elts:
        if perm_inv(p) not in eltset:
            return False
    for p in elts:
        for q in elts:
            if perm_compose(p, q) not in eltset:
                return False
    return True


def left_coset_test(S, sample_size=None):
    """Given a set S of perms, test if S = π_0 K for some subgroup K.
    K = π_0^{-1} S; return (K, π_0) if K is subgroup, else None.
    """
    S = list(S)
    if not S:
        return None
    pi0 = S[0]
    pi0_inv = perm_inv(pi0)
    K = [perm_compose(pi0_inv, s) for s in S]
    if is_subgroup(K):
        return K, pi0
    return None


def right_coset_test(S):
    S = list(S)
    if not S:
        return None
    pi0 = S[0]
    pi0_inv = perm_inv(pi0)
    K = [perm_compose(s, pi0_inv) for s in S]
    if is_subgroup(K):
        return K, pi0
    return None


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    print("=" * 78)
    print("Phase 9.18.S9x — Analytic mechanism of the M19 macro-orbit (mass 72)")
    print("=" * 78)

    bases = load_bases()
    L_M19 = bases["M19"]
    print(f"\nM19 grid:\n{L_M19}")

    certs = load_m19_certs()
    print(f"\n|Γ_M19| (odd-rank perms with corank=1): {len(certs)}")

    # ---- (X1) Cluster by canonical_full ----
    print()
    print("=" * 78)
    print("(X1) Clustering M19 certs by canonical full G_template orbit")
    print("=" * 78)

    by_orbit = defaultdict(list)
    for c in certs:
        ks, kf = canonicalize_full(c["side"], c["u9_simpler"])
        by_orbit[kf].append(c)
    sizes = sorted([(len(v), k) for k, v in by_orbit.items()], reverse=True)
    print(f"\nNumber of distinct canonical-full orbits: {len(by_orbit)}")
    print(f"\nTop 5 by mass:")
    for sz, k in sizes[:5]:
        print(f"   mass={sz:>3}  rep_key=({k[0]}, {list(k[1])})")

    # The dominant orbit
    O_key = sizes[0][1]
    O_size = sizes[0][0]
    O_certs = by_orbit[O_key]
    print(f"\n=> O_M19 (dominant): canonical (side, vec) = ({O_key[0]}, {list(O_key[1])})")
    print(f"   |O_M19 ∩ Γ_M19| = {O_size}")
    u9_canonical_simpler = tuple(O_key[1])

    # ---- (X2) Pick representative u9 (in V_std∩Z^9) for O_M19 ----
    # O_M19's canonical "simpler-side" form is u9_canonical_simpler (length 9, sum 0).
    # Verify primitive and sign-normalized.
    print()
    print(f"u9* (canonical orbit rep, simpler side):  {u9_canonical_simpler}")
    print(f"   sum = {sum(u9_canonical_simpler)}  (should be 0)")
    print(f"   support = {sum(1 for x in u9_canonical_simpler if x != 0)}")
    print(f"   ||·||_∞ = {max(abs(x) for x in u9_canonical_simpler)}")
    print(f"   value multiset = {Counter(u9_canonical_simpler).most_common()}")

    # ---- (X3) Stabilizer of u9* in S_9 (coordinate permutation) ----
    print()
    print("=" * 78)
    print("(X3) Stab_{S_9}(u9*) — coordinate-permutation stabilizer")
    print("=" * 78)
    H_u = stabilizer_of_vec(u9_canonical_simpler)
    print(f"|Stab(u9*)| = {len(H_u)}")
    # Decompose by value class for sanity
    from math import factorial
    val_mult = Counter(u9_canonical_simpler)
    expected = 1
    for k, m in val_mult.items():
        expected *= factorial(m)
    print(f"Expected from value multiplicities: ∏ m_v! = {expected}")
    print(f"Match: {len(H_u) == expected}")

    # ---- (X4) Decompose each γ ∈ Γ*_M19 as π_γ · u9* ----
    print()
    print("=" * 78)
    print("(X4) For each γ ∈ Γ*_M19, find π_γ ∈ S_9 with π_γ · u9* = u9_simpler(γ)")
    print("=" * 78)

    pi_list = []
    failures = 0
    for c in O_certs:
        u9 = c["u9_simpler"]
        # Adapt: find_pi_mapping(target=u9, source=u9_canonical_simpler)
        pis = find_pi_mapping(list(u9), list(u9_canonical_simpler))
        if not pis:
            failures += 1
            continue
        pi_list.append({"perm": c["perm"], "side": c["side"], "pi": pis[0],
                        "n_pi_options": len(pis)})

    print(f"\nDecomposition successful for: {len(pi_list)} / {len(O_certs)}")
    print(f"Failures: {failures}")
    if pi_list:
        cnts = Counter(x["n_pi_options"] for x in pi_list)
        print(f"Number-of-π options distribution: {dict(cnts)}")
        print(f"   (each γ admits |Stab(u9*)| = {len(H_u)} options — pick canonical = lex min)")

    # ---- (X5) Test left/right coset structure of Γ*_M19 ⊂ S_9 ----
    print()
    print("=" * 78)
    print("(X5) Coset structure of Γ*_M19 ⊂ S_9")
    print("=" * 78)
    Gamma_star = [c["perm"] for c in O_certs]
    print(f"|Γ*_M19| = {len(Gamma_star)}")

    # Test: is Γ*_M19 closed under left mult by some subgroup K?
    # Compute left "stabilizer-of-set" Stab_L = {k : k Γ ⊆ Γ}
    Gset = set(Gamma_star)
    print("\n(H-B) Left coset test: is Γ*_M19 a single LEFT coset π_0 K?")
    res = left_coset_test(Gamma_star)
    if res:
        K, pi0 = res
        print(f"   YES — K of order {len(K)} with anchor π_0 = {pi0}")
        # Identify K: cycle types
        cyc_types = Counter()
        for k in K:
            cyc_types[cycle_type(k)] += 1
        print(f"   K cycle types: {dict(cyc_types)}")
    else:
        print(f"   NO — Γ*_M19 is not a single left coset.")

    print("\n(H-B') Right coset test:")
    res = right_coset_test(Gamma_star)
    if res:
        K, pi0 = res
        print(f"   YES — K of order {len(K)} with anchor π_0 = {pi0}")
        cyc_types = Counter()
        for k in K:
            cyc_types[cycle_type(k)] += 1
        print(f"   K cycle types: {dict(cyc_types)}")
    else:
        print(f"   NO — Γ*_M19 is not a single right coset.")

    # Test (H-A): largest K with K · Γ ⊆ Γ
    print("\n(H-A) Largest K ≤ S_9 with K · Γ*_M19 ⊆ Γ*_M19 (left-stabilizer of set)")
    # The left-stabilizer of a set: {k : k γ ∈ Γ ∀ γ ∈ Γ}
    # = ∩_γ (Γ γ^{-1})
    # = computed as set intersection
    cand_L = None
    for g in Gamma_star:
        ginv = perm_inv(g)
        coset = set(perm_compose(s, ginv) for s in Gamma_star)
        if cand_L is None:
            cand_L = coset
        else:
            cand_L &= coset
        if len(cand_L) <= 1:
            break
    K_L_size = len(cand_L) if cand_L else 0
    print(f"   |left stabilizer K_L| = {K_L_size}")
    K_L = list(cand_L) if cand_L else []
    K_L_subgroup = False
    K_L_cycle_types = {}
    if 1 < K_L_size <= 200:
        K_L_subgroup = is_subgroup(K_L)
        print(f"   Subgroup test: {K_L_subgroup}")
        K_L_cycle_types = dict(Counter(cycle_type(k) for k in K_L))
        print(f"   Cycle types: {K_L_cycle_types}")

    # And right-stabilizer
    print("\n(H-A') Largest K_R with Γ · K_R ⊆ Γ (right-stabilizer)")
    cand_R = None
    for g in Gamma_star:
        ginv = perm_inv(g)
        coset = set(perm_compose(ginv, s) for s in Gamma_star)
        if cand_R is None:
            cand_R = coset
        else:
            cand_R &= coset
        if len(cand_R) <= 1:
            break
    K_R_size = len(cand_R) if cand_R else 0
    print(f"   |right stabilizer K_R| = {K_R_size}")
    K_R = list(cand_R) if cand_R else []
    K_R_subgroup = False
    K_R_cycle_types = {}
    if 1 < K_R_size <= 200:
        K_R_subgroup = is_subgroup(K_R)
        print(f"   Subgroup test: {K_R_subgroup}")
        K_R_cycle_types = dict(Counter(cycle_type(k) for k in K_R))
        print(f"   Cycle types: {K_R_cycle_types}")

    # Are K_L and K_R the same set?
    if K_L and K_R:
        same = (set(K_L) == set(K_R))
        print(f"\n   K_L == K_R as sets ? {same}")
    else:
        same = False

    # ---- (X6) Side-mix in Γ*_M19 ----
    print()
    print("=" * 78)
    print("(X6) Side-mix and small-orbit residue")
    print("=" * 78)
    side_mix = Counter(c["side"] for c in O_certs)
    print(f"\nSide-mix in Γ*_M19: {dict(side_mix)}")

    Gamma_residue = [c for c in certs if c["perm"] not in Gset]
    print(f"\nResidue Γ_M19 \\ Γ*_M19: {len(Gamma_residue)} certs")
    res_orbits = defaultdict(list)
    for c in Gamma_residue:
        ks, kf = canonicalize_full(c["side"], c["u9_simpler"])
        res_orbits[kf].append(c)
    res_sizes = sorted([(len(v), k) for k, v in res_orbits.items()], reverse=True)
    print(f"Number of residue orbits: {len(res_orbits)}")
    print(f"Mass distribution: {Counter(sz for sz, _ in res_sizes).most_common()}")
    print(f"Top 5 residue orbits:")
    for sz, k in res_sizes[:5]:
        print(f"   mass={sz}  ({k[0]}, {list(k[1])})")

    # Save
    out = {
        "phase": "9.18.S9x",
        "M19_grid": L_M19.tolist(),
        "n_certs_M19": len(certs),
        "u9_star": list(u9_canonical_simpler),
        "stab_size": len(H_u),
        "value_multiplicities": {str(k): v for k, v in val_mult.items()},
        "O_M19_size": O_size,
        "O_M19_side_mix": dict(side_mix),
        "n_orbits_M19": len(by_orbit),
        "left_coset_test": "YES" if left_coset_test(Gamma_star) else "NO",
        "right_coset_test": "YES" if right_coset_test(Gamma_star) else "NO",
        "K_L_size": K_L_size,
        "K_R_size": K_R_size,
        "K_L_subgroup": K_L_subgroup,
        "K_R_subgroup": K_R_subgroup,
        "K_L_cycle_types": {str(k): v for k, v in K_L_cycle_types.items()},
        "K_R_cycle_types": {str(k): v for k, v in K_R_cycle_types.items()},
        "K_L_eq_K_R": same,
        "K_L_elements": [list(k) for k in K_L],
        "K_R_elements": [list(k) for k in K_R],
        "residue_size": len(Gamma_residue),
        "residue_n_orbits": len(res_orbits),
        "residue_mass_dist": dict(Counter(sz for sz, _ in res_sizes)),
        "Gamma_star_perms": [list(p) for p in Gamma_star],
    }
    with open(HERE / "data" / "phase_9_18_S9x_results.json", "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: data/phase_9_18_S9x_results.json")


def cycle_type(p):
    n = len(p)
    seen = [False] * n
    cycles = []
    for i in range(n):
        if seen[i]:
            continue
        j = i
        c = 0
        while not seen[j]:
            seen[j] = True
            j = p[j] - 1
            c += 1
        cycles.append(c)
    return tuple(sorted(cycles, reverse=True))


if __name__ == "__main__":
    main()
