"""Phase 9.15.A — Identify the abstract group H_V1 ⊂ B_3.

H_V1 has 12 elements with order distribution {1:1, 2:7, 3:2, 6:2}. The
two candidate groups of order 12 are:
  - Z/12        : orders {1:1, 2:1, 3:2, 4:2, 6:2, 12:4}  -- mismatch
  - Z/6 x Z/2   : orders {1:1, 2:3, 3:2, 6:6}             -- mismatch
  - A_4         : orders {1:1, 2:3, 3:8}                  -- mismatch
  - D_6 (dihedral, order 12) : orders {1:1, 2:7, 3:2, 6:2} -- MATCH
  - Dic_3       : orders {1:1, 2:1, 3:2, 4:6, 6:2}        -- mismatch

So order distribution alone identifies H_V1 ≅ D_6 ≅ S_3 × Z/2.

We verify constructively by:
  (1) Pick generators: r of order 6, s of order 2 with srs = r^{-1}.
  (2) Show every element is r^k or r^k · s.
"""
import json
from pathlib import Path
from collections import Counter

DATA = Path(__file__).resolve().parent.parent.parent / "data" / "phase_9_15_results.json"


def perm_compose(a, b):
    return tuple(a[b[k] - 1] for k in range(len(a)))


def perm_inverse(a):
    out = [0] * len(a)
    for k, v in enumerate(a):
        out[v - 1] = k + 1
    return tuple(out)


def perm_order(p):
    e = tuple(range(1, len(p) + 1))
    q = p
    o = 1
    while q != e:
        q = perm_compose(p, q)
        o += 1
    return o


def cycle_type(p):
    n = len(p)
    seen = [False] * n
    cyc = []
    for i in range(n):
        if seen[i]:
            continue
        j = i
        L = 0
        while not seen[j]:
            seen[j] = True
            j = p[j] - 1
            L += 1
        cyc.append(L)
    return tuple(sorted(cyc, reverse=True))


def main():
    data = json.load(DATA.open())
    perms_v1 = [tuple(p) for p in data["B3_enrichment"]["L_V1"]["odd_perms_in_B3"]]
    H = set(perms_v1)
    n = 7
    e = tuple(range(1, n + 1))

    print(f"|H_V1| = {len(H)}")
    orders = Counter(perm_order(p) for p in H)
    print(f"order distribution: {dict(orders)}")

    # Candidate group identification by order distribution
    expected_d6 = {1: 1, 2: 7, 3: 2, 6: 2}
    print(f"\nExpected for D_6 (order 12): {expected_d6}")
    print(f"Match D_6? {dict(orders) == expected_d6}")

    # Constructive verification: find r of order 6, s of order 2 with srs = r^{-1}
    r = next(p for p in H if perm_order(p) == 6)
    print(f"\nr (order 6): {r}")

    # Generate <r>
    cyclic_r = [e]
    cur = r
    while cur != e:
        cyclic_r.append(cur)
        cur = perm_compose(r, cur)
    print(f"<r> has size {len(cyclic_r)} ✓")
    assert len(cyclic_r) == 6
    assert all(x in H for x in cyclic_r)

    # Find s ∈ H \ <r> of order 2 with s r s = r^{-1}
    r_inv = perm_inverse(r)
    s_found = None
    for s in H:
        if s in cyclic_r:
            continue
        if perm_order(s) != 2:
            continue
        if perm_compose(perm_compose(s, r), s) == r_inv:
            s_found = s
            break

    print(f"s (order 2, srs = r^-1): {s_found}")
    assert s_found is not None, "No dihedral generator found — H_V1 is NOT D_6"

    # Verify generated set <r,s> = H
    generated = set()
    for k in range(6):
        rk = cyclic_r[k]
        generated.add(rk)
        generated.add(perm_compose(rk, s_found))
    print(f"|<r,s>| = {len(generated)}; equals H_V1? {generated == H}")
    assert generated == H

    print("\n✅ THEOREM (Phase 9.15.A): H_V1 ≅ D_6 (dihedral group of order 12).")
    print("   Equivalently, H_V1 ≅ S_3 × Z/2.")

    # Cycle-type breakdown (for SVW comparison)
    print("\nCycle-type distribution within H_V1:")
    for ct, cnt in sorted(Counter(cycle_type(p) for p in H).items()):
        print(f"  {ct}: {cnt}")

    # Save certificate
    out = {
        "abstract_group": "D_6 (dihedral, order 12) ≅ S_3 × Z/2",
        "order_distribution": dict(orders),
        "generator_r_order_6": list(r),
        "generator_s_order_2": list(s_found),
        "dihedral_relation_verified": True,
        "generates_full_H_V1": True,
        "cycle_type_distribution": {str(k): v for k, v in Counter(cycle_type(p) for p in H).items()},
    }
    out_path = Path(__file__).resolve().parent.parent.parent / "data" / "phase_9_15A_H_V1_identification.json"
    json.dump(out, out_path.open("w"), indent=2)
    print(f"\nCertificate saved to {out_path}")


if __name__ == "__main__":
    main()
