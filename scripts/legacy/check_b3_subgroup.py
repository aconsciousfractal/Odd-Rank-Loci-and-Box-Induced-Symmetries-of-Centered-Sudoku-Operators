"""Phase 9.15 sub-analysis: subgroup/coset structure within B_3."""
import json
from itertools import combinations
from pathlib import Path

data_path = Path(__file__).resolve().parent.parent.parent / "data" / "phase_9_15_results.json"
data = json.load(data_path.open())


def perm_compose(a, b):
    n = len(a)
    return tuple(a[b[k] - 1] for k in range(n))


def perm_inverse(a):
    n = len(a)
    out = [0] * n
    for k in range(n):
        out[a[k] - 1] = k + 1
    return tuple(out)


def is_subgroup(S, n):
    e = tuple(range(1, n + 1))
    if e not in S:
        return False, 'no identity'
    for a in S:
        if perm_inverse(a) not in S:
            return False, 'no inv'
    for a in S:
        for b in S:
            if perm_compose(a, b) not in S:
                return False, 'not closed'
    return True, 'OK'


pair_imap = [(1, 7), (2, 6), (3, 5)]


def b3_sign_within(p):
    s = 1
    for a, b in pair_imap:
        if p[a - 1] == a and p[b - 1] == b:
            pass
        elif p[a - 1] == b and p[b - 1] == a:
            s *= -1
    return s


def b3_sign_action(p):
    blocks = {1: 0, 7: 0, 2: 1, 6: 1, 3: 2, 5: 2}
    sigma = [blocks[p[a - 1]] for (a, _) in pair_imap]
    invs = sum(1 for i, j in combinations(range(3), 2) if sigma[i] > sigma[j])
    return -1 if invs % 2 else 1


for label in ('L_V1', 'L_V3'):
    perms = [tuple(p) for p in data['B3_enrichment'][label]['odd_perms_in_B3']]
    print(f"\n== {label} ({len(perms)} perms in B_3 with odd rank) ==")
    is_sg, reason = is_subgroup(set(perms), 7)
    print(f"  Is subgroup of S_7? {is_sg} ({reason})")

    sgn_total = [b3_sign_within(p) * b3_sign_action(p) for p in perms]
    pos = sum(1 for s in sgn_total if s > 0)
    neg = sum(1 for s in sgn_total if s < 0)
    print(f"  total-sign distribution: +1:{pos}, -1:{neg}")
    sgn_w = [b3_sign_within(p) for p in perms]
    sgn_a = [b3_sign_action(p) for p in perms]
    print(f"  within-pair sign: +1:{sum(1 for s in sgn_w if s > 0)}, "
          f"-1:{sum(1 for s in sgn_w if s < 0)}")
    print(f"  S_3 action sign: +1:{sum(1 for s in sgn_a if s > 0)}, "
          f"-1:{sum(1 for s in sgn_a if s < 0)}")

    # If is subgroup of size 24, would equal kernel of some Z/2 character.
    # Test: do they all have total-sign = -1? +1?
    if pos == len(perms):
        print("  ALL have total sign +1 -> coincide with B_3 ∩ ker(sign_total)")
    elif neg == len(perms):
        print("  ALL have total sign -1 -> coincide with B_3 \\ ker(sign_total)")
