"""Phase 9.18.S9q-γ-1 — Estensione universale a n ∈ {5,6,7,8,9,10,11,12,13}.

Promotes Theorem 8.25.sexdecies.1 from n ∈ {7,9,11} to a wider range, using
the combinatorial Stab + BFS-orbit + lattice + depth pipeline of β-1, now
parametric in n.

Predictions (universal dichotomy, β-3 extrapolation):

    M_n = (-1,-1, 0,...,0, 1, 1)    (2 neg, n-4 zeros, 2 pos)
    O_n = (-2,-1, 0,...,0, 1, 2)    (singleton classes apart from zero block)

  |Stab(M_n)| = 8 · (n-4)!   ≅ D_4 × S_{n-4}
  |Stab(O_n)| = 2 · (n-4)!   ≅ C_2 × S_{n-4}
  |orb(M_n)|  = n(n-1)(n-2)(n-3)/4
  |orb(O_n)|  = n(n-1)(n-2)(n-3)
  Λ_{M_n} = Λ_{O_n} = A_{n-1}
  d_{M_n}(orb(O_n)) = 2

Edge cases:
  n=5: zero-class size = 1, S_1 trivial. M_5 still has supp 4 (the only 0 is
       a single position). Predict |Stab(M_5)| = 8, |orb|=15.
  n=6: zero-class size = 2. |Stab(M_6)| = 16, |orb|=45.
  n>=12,13: BFS orb feasible (|orb(M_13)|=4290, |orb(O_13)|=17160).

Output: data/phase_9_18_S9q_gamma1.json + _output.txt
"""
from __future__ import annotations
import json
import time
from collections import Counter
from itertools import permutations, product
from math import factorial
from pathlib import Path

import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form, hermite_normal_form

ROOT = Path(__file__).resolve().parents[2]  # <repo>/scripts/legacy/ -> <repo>/
DATA = ROOT / "data"

# Range to verify. n=4 excluded (no zero block, edge-case for the dichotomy
# narrative — supp 4 = full vector). n=7,9,11 included as regression check.
N_RANGE = [5, 6, 7, 8, 9, 10, 11, 12, 13]


# --------------------------------------------------------------------
# Combinatorial Stab construction (parametric in n).
# --------------------------------------------------------------------

def value_classes(v):
    classes = {}
    for i, val in enumerate(v):
        classes.setdefault(int(val), []).append(i)
    return {k: tuple(sorted(p)) for k, p in classes.items()}


def stabilizer_combinatorial(v):
    """Return list of (perm_tuple, sign) pairs in Stab_{S_n×{±1}}(v).

    σ encodes the image: σ[i] = j means i ↦ j.
    Action: (σ, ε) · v defined as (ε·v)[σ(i)] = v[i] (i.e. v_new[σ(i)] = ε·v[i]).
    Equivalently, after action, the value at position σ(i) is ε·v[i].
    Stab condition: ε·v = v ∘ σ, i.e. ε·v[σ(i)] = v[i] is wrong direction.

    We use convention: applied = ε * v[σ_arr]  (numpy fancy index).
    This means new[i] = ε · v[σ[i]]. For this to equal v: v[σ[i]] = ε^{-1} v[i] = ε·v[i].
    So ε=+1: σ permutes within value-classes of v.
       ε=-1: σ maps class(c) ↔ class(-c).
    """
    v = tuple(int(x) for x in v)
    n = len(v)
    classes = value_classes(v)
    values = sorted(classes.keys())

    elems = []

    # ε = +1: independent perms within each class.
    pos_lists = [classes[c] for c in values]
    perm_lists = [list(permutations(p)) for p in pos_lists]
    for choice in product(*perm_lists):
        sigma = list(range(n))
        for src_positions, dst_positions in zip(pos_lists, choice):
            for s, d in zip(src_positions, dst_positions):
                sigma[s] = d
        elems.append((tuple(sigma), 1))

    # ε = -1: σ maps class(c) → class(-c).
    pos_classes = [c for c in values if c > 0]
    zero_class = classes.get(0, ())

    paired_ok = all(
        (-c) in classes and len(classes[-c]) == len(classes[c])
        for c in pos_classes
    )
    if paired_ok:
        bijection_choices = []
        for c in pos_classes:
            src = classes[c]
            tgt = list(permutations(classes[-c]))
            bijection_choices.append((src, tgt))
        bijection_choices_neg = []
        for c in pos_classes:
            src = classes[-c]
            tgt = list(permutations(classes[c]))
            bijection_choices_neg.append((src, tgt))
        zero_perm_list = list(permutations(zero_class)) if zero_class else [()]

        for pos_imgs in product(*[t for _, t in bijection_choices]):
            for neg_imgs in product(*[t for _, t in bijection_choices_neg]):
                for zero_img in zero_perm_list:
                    sigma = list(range(n))
                    for (src, _), img in zip(bijection_choices, pos_imgs):
                        for s, d in zip(src, img):
                            sigma[s] = d
                    for (src, _), img in zip(bijection_choices_neg, neg_imgs):
                        for s, d in zip(src, img):
                            sigma[s] = d
                    for s, d in zip(zero_class, zero_img):
                        sigma[s] = d
                    elems.append((tuple(sigma), -1))

    # Sanity sample.
    v_arr = np.asarray(v, dtype=np.int64)
    for sigma, eps in elems[: min(20, len(elems))]:
        sigma_arr = np.asarray(sigma, dtype=np.int64)
        applied = eps * v_arr[sigma_arr]
        assert np.array_equal(applied, v_arr), (
            f"Stab violation: σ={sigma}, ε={eps}, applied={applied.tolist()} "
            f"vs v={v_arr.tolist()}"
        )
    return elems


def project_on_subset(stab, subset):
    """Restrict each (σ, ε) ∈ Stab to subset (which must be σ-invariant when ε=+1).

    If ε=-1 and zero-class is non-empty, σ permutes zero-class within itself
    too (since 0 = -0). If ε=-1 and zero-class is empty, no restriction issue.
    """
    subset_index = {p: i for i, p in enumerate(subset)}
    images = set()
    kernel = []
    for sigma, eps in stab:
        # subset is the zero-class → σ must permute it within itself.
        induced = tuple(subset_index[sigma[p]] for p in subset)
        images.add(induced)
        if induced == tuple(range(len(subset))):
            kernel.append((sigma, eps))
    return images, kernel


def order_of_pair(g, n):
    sigma_id = tuple(range(n))
    cur = g
    cnt = 1
    while cur != (sigma_id, 1):
        sa, ea = g
        sb, eb = cur
        cur = (tuple(sa[sb[i]] for i in range(n)), ea * eb)
        cnt += 1
        if cnt > 10000:
            raise RuntimeError("order > 10000")
    return cnt


# --------------------------------------------------------------------
# Orbit BFS (full G*_n via S_n × {±1} generators).
# --------------------------------------------------------------------

def orbit_via_generators(v, n):
    gens = []
    for i in range(n - 1):
        s = list(range(n))
        s[i], s[i + 1] = s[i + 1], s[i]
        gens.append((tuple(s), 1))
    gens.append((tuple(range(n)), -1))
    v_arr = np.asarray(v, dtype=np.int64)
    seen = {tuple(int(x) for x in v_arr)}
    frontier = [v_arr]
    while frontier:
        new = []
        for u in frontier:
            for sigma, eps in gens:
                applied = eps * u[np.asarray(sigma, dtype=np.int64)]
                key = tuple(int(x) for x in applied)
                if key not in seen:
                    seen.add(key)
                    new.append(applied)
        frontier = new
    return np.array([list(t) for t in seen], dtype=np.int64)


# --------------------------------------------------------------------
# Lattice invariants (HNF/SNF).
# --------------------------------------------------------------------

def lattice_invariants(M_int, label):
    M = sp.Matrix(M_int.tolist())
    H = hermite_normal_form(M.T)
    nonzero_cols = [j for j in range(H.shape[1])
                    if any(H[i, j] != 0 for i in range(H.shape[0]))]
    rank = len(nonzero_cols)
    H_red = H[:, nonzero_cols]
    S = smith_normal_form(H_red)
    elem_div = [int(S[i, i]) for i in range(min(S.shape)) if S[i, i] != 0]
    print(f"      [{label}] shape={M_int.shape}  rank={rank}  "
          f"elem_div={elem_div}", flush=True)
    return {"rank": rank, "elementary_divisors": elem_div}


# --------------------------------------------------------------------
# Per-n analysis.
# --------------------------------------------------------------------

def build_M_n(n):
    return tuple([-1, -1] + [0] * (n - 4) + [1, 1])


def build_O_n(n):
    return tuple([-2, -1] + [0] * (n - 4) + [1, 2])


def predicted(n):
    s_nm4 = factorial(n - 4)
    return {
        "stab_M": 8 * s_nm4,
        "stab_O": 2 * s_nm4,
        "orb_M": n * (n - 1) * (n - 2) * (n - 3) // 4,
        "orb_O": n * (n - 1) * (n - 2) * (n - 3),
    }


def analyze_n(n):
    print()
    print("#" * 78)
    print(f"## n = {n}   |G*_n| = {factorial(n)*2:,}   "
          f"|Stab(M)|=8·(n-4)!={8*factorial(n-4)}, "
          f"|Stab(O)|=2·(n-4)!={2*factorial(n-4)}")
    print("#" * 78, flush=True)

    pred = predicted(n)
    M_n = build_M_n(n)
    O_n = build_O_n(n)

    print(f"  M_{n} = {list(M_n)}", flush=True)
    print(f"  O_{n} = {list(O_n)}", flush=True)

    results = {"n": n, "M_rep": list(M_n), "O_rep": list(O_n),
               "predicted": pred, "checks": {}}

    # --- Stab combinatorial ---
    for label, v, pred_stab, ker_name, ker_order in [
        ("M", M_n, pred["stab_M"], "D_4", 8),
        ("O", O_n, pred["stab_O"], "C_2", 2),
    ]:
        t0 = time.time()
        stab = stabilizer_combinatorial(v)
        dt = time.time() - t0
        ok_size = (len(stab) == pred_stab)
        print(f"    [{label}_{n}] |Stab| = {len(stab)}  pred {pred_stab}  "
              f"{'OK' if ok_size else 'MISMATCH'}  ({dt:.2f}s)", flush=True)
        assert ok_size

        zero_pos = list(value_classes(v).get(0, ()))
        if zero_pos:
            images, kernel = project_on_subset(stab, zero_pos)
            ok_im = len(images) == factorial(n - 4)
            ok_ker = len(kernel) == ker_order
            ok_split = len(stab) == len(images) * len(kernel)
            print(f"        |im(ψ→S_{n-4})|={len(images)} "
                  f"(pred {factorial(n-4)})  "
                  f"|ker|={len(kernel)} (pred {ker_order} for {ker_name})  "
                  f"split: {ok_split}", flush=True)
            assert ok_im and ok_ker and ok_split

            ker_orders = Counter(order_of_pair(g, n) for g in kernel)
            if ker_name == "D_4":
                exp = {1: 1, 2: 5, 4: 2}
            else:
                exp = {1: 1, 2: 1}
            ok_orders = (dict(ker_orders) == exp)
            print(f"        ker order-multiset = {dict(sorted(ker_orders.items()))} "
                  f"vs {ker_name} pred {exp}  "
                  f"{'OK' if ok_orders else 'MISMATCH'}", flush=True)
            assert ok_orders
        else:
            # n=4 edge case: zero-class empty. We don't test n=4.
            ok_im = ok_ker = ok_split = ok_orders = True

        results["checks"][f"stab_{label}"] = {
            "size_ok": ok_size,
            "split_ok": ok_split,
            "kernel_orders_ok": ok_orders,
        }

    # --- Orbit BFS ---
    print(f"    BFS orb(M_{n}) …", flush=True)
    t0 = time.time()
    orb_M = orbit_via_generators(M_n, n)
    print(f"      |orb(M_{n})| = {orb_M.shape[0]}  pred {pred['orb_M']}  "
          f"({time.time()-t0:.2f}s)", flush=True)
    assert orb_M.shape[0] == pred["orb_M"]

    print(f"    BFS orb(O_{n}) …", flush=True)
    t0 = time.time()
    orb_O = orbit_via_generators(O_n, n)
    print(f"      |orb(O_{n})| = {orb_O.shape[0]}  pred {pred['orb_O']}  "
          f"({time.time()-t0:.2f}s)", flush=True)
    assert orb_O.shape[0] == pred["orb_O"]

    assert np.all(orb_M.sum(axis=1) == 0)
    assert np.all(orb_O.sum(axis=1) == 0)

    # --- Lattice ---
    inv_M = lattice_invariants(orb_M, f"Λ_{{M_{n}}}")
    inv_O = lattice_invariants(orb_O, f"Λ_{{O_{n}}}")
    sat_M = (inv_M["rank"] == n - 1
             and inv_M["elementary_divisors"] == [1] * (n - 1))
    sat_O = (inv_O["rank"] == n - 1
             and inv_O["elementary_divisors"] == [1] * (n - 1))
    print(f"      Λ_M = A_{n-1}? {sat_M}    Λ_O = A_{n-1}? {sat_O}",
          flush=True)
    results["Lambda_M"] = {**inv_M, "saturates": sat_M}
    results["Lambda_O"] = {**inv_O, "saturates": sat_O}

    # --- Depth d_M(orb(O)) ---
    M_set = {tuple(int(x) for x in r) for r in orb_M}
    Oa = np.array(O_n, dtype=np.int64)
    d1 = tuple(int(x) for x in Oa) in M_set
    d2_witness = None
    if not d1:
        diffs = Oa[None, :] - orb_M
        for i, row in enumerate(diffs):
            if tuple(int(x) for x in row) in M_set:
                d2_witness = (orb_M[i].tolist(), row.tolist())
                break
    if d1:
        depth = 1
    elif d2_witness is not None:
        depth = 2
    else:
        depth = ">=3"
    print(f"      d_{{M_{n}}}(O_{n}) = {depth}  (pred 2)", flush=True)
    results["depth"] = depth
    if d2_witness is not None:
        results["depth_witness"] = {"a": d2_witness[0], "b": d2_witness[1]}

    # --- Per-n bridges ---
    bridges = {
        "stab_M_size": results["checks"]["stab_M"]["size_ok"],
        "stab_O_size": results["checks"]["stab_O"]["size_ok"],
        "stab_M_D4xS": results["checks"]["stab_M"]["kernel_orders_ok"]
                       and results["checks"]["stab_M"]["split_ok"],
        "stab_O_C2xS": results["checks"]["stab_O"]["kernel_orders_ok"]
                       and results["checks"]["stab_O"]["split_ok"],
        "orb_M_size": True,  # already asserted
        "orb_O_size": True,
        "Lambda_M_saturates": sat_M,
        "Lambda_O_saturates": sat_O,
        "depth_eq_2": (depth == 2),
    }
    results["bridges"] = bridges
    results["all_pass"] = all(bridges.values())
    print(f"    ⇒ n={n}: ALL_PASS = {results['all_pass']}", flush=True)
    return results


def main():
    print("=" * 78)
    print("Phase 9.18.S9q-γ-1 — Estensione universale n ∈", N_RANGE)
    print("=" * 78, flush=True)

    all_results = []
    t_global = time.time()
    for n in N_RANGE:
        try:
            res = analyze_n(n)
            all_results.append(res)
        except AssertionError as e:
            print(f"  *** FAIL at n={n}: {e}", flush=True)
            all_results.append({"n": n, "all_pass": False, "error": str(e)})

    dt = time.time() - t_global
    print()
    print("=" * 78)
    print("  γ-1 SUMMARY")
    print("=" * 78)
    print(f"  total time: {dt:.2f}s", flush=True)
    print()
    print(f"  {'n':>3}  {'|Stab(M)|':>10}  {'|Stab(O)|':>10}  "
          f"{'|orb(M)|':>10}  {'|orb(O)|':>10}  "
          f"{'Λ_M':>5}  {'Λ_O':>5}  {'d':>3}  PASS")
    for r in all_results:
        if "predicted" not in r:
            print(f"  {r['n']:>3}  -- ERROR: {r.get('error', '?')}")
            continue
        p = r["predicted"]
        sm = "A" if r["Lambda_M"]["saturates"] else "x"
        so = "A" if r["Lambda_O"]["saturates"] else "x"
        print(f"  {r['n']:>3}  {p['stab_M']:>10}  {p['stab_O']:>10}  "
              f"{p['orb_M']:>10}  {p['orb_O']:>10}  {sm:>5}  {so:>5}  "
              f"{str(r['depth']):>3}  {'OK' if r['all_pass'] else 'FAIL'}")

    universal = all(r.get("all_pass") for r in all_results)
    print()
    if universal:
        print("  ★ Universal pattern verified for n ∈", N_RANGE, ".")
        print("  ★ Theorem 8.25.sexdecies.1 elevated to: ∀ n ∈ {5,..,13}.")
        print("    Conjecture 8.25.sexdecies.3 (∀ n ≥ 5) reinforced.")
    else:
        print("  ✗ Some n failed. Pattern not universal in this range.")

    out = {
        "phase": "9.18.S9q-γ-1",
        "n_range": N_RANGE,
        "results": all_results,
        "universal_pass": universal,
        "elapsed_sec": dt,
    }
    out_path = DATA / "phase_9_18_S9q_gamma1.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\n  saved → {out_path}", flush=True)


if __name__ == "__main__":
    main()
