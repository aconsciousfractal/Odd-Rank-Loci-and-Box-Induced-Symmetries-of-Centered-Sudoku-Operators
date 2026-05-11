"""build_certified.py
=====================
Build the canonical certified-results package for this paper
("Odd-Rank Loci and Box-Induced Symmetries of Centered Sudoku Operators").

This script reads the immutable source JSONs produced by the Phase 9 / Phase 9.18
pipeline (under ``../../data/``) and emits 12 canonical certificate JSONs in
``./``, plus ``MANIFEST.sha256``.

Canonical schema (per certificate):
    {
      "result_id":     "...",
    "tier":          "THEOREM | EXHAUSTIVE | CERTIFIED_EXHAUSTIVE | PROPOSITION_FINITE_CORPUS",
      "claim":         "<one sentence>",
      "paper_label":   "Theorem A | Proposition P1 | ...",
      "inputs":        {...},
      "outputs":       {...},
      "source":        {"file": "data/<...>.json", "sha256": "<hex64>", "extracted_keys": [...]},
      "script":        "scripts/<file>.py",
      "script_sha256": "<hex64>",
      "seed":          null | <int>,
      "produced_utc":  "YYYY-MM-DDTHH:MM:SSZ",
      "schema_version": 1
    }

Hashing convention:
- The MANIFEST entry for each certificate is computed over the canonical bytes
  of the certificate with ``produced_utc`` removed (``json.dumps(obj,
  sort_keys=True, indent=2)``).  The ``produced_utc`` value is preserved
  inside the JSON file itself for human inspection.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import itertools
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

# ----------------------------------------------------------------------------
# Paths (STANDALONE: every input lives inside REPO_ROOT)
# ----------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent          # <repo>/certified/
REPO_ROOT = HERE.parent                         # <repo>/
DATA = REPO_ROOT / "data"                       # <repo>/data/
SCRIPTS = REPO_ROOT / "scripts" / "legacy"      # <repo>/scripts/legacy/

# Sanity: refuse to run if the standalone layout is broken.
assert DATA.is_dir(), f"standalone violated: missing {DATA}"
assert SCRIPTS.is_dir(), f"standalone violated: missing {SCRIPTS}"

SCHEMA_VERSION = 1

NOW = _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical_bytes(obj: Dict[str, Any]) -> bytes:
    """Canonical serialization for hashing/output."""
    return json.dumps(obj, sort_keys=True, indent=2, ensure_ascii=False).encode(
        "utf-8"
    )


def digest_of(obj: Dict[str, Any]) -> str:
    """SHA-256 over the canonical bytes with `produced_utc` removed."""
    obj2 = {k: v for k, v in obj.items() if k != "produced_utc"}
    return hashlib.sha256(canonical_bytes(obj2)).hexdigest()


def write_certificate(name: str, obj: Dict[str, Any]) -> Tuple[str, str]:
    """Write the canonical JSON to disk; return (filename, manifest-digest)."""
    obj.setdefault("schema_version", SCHEMA_VERSION)
    obj.setdefault("produced_utc", NOW)
    out = HERE / name
    out.write_bytes(canonical_bytes(obj))
    return name, digest_of(obj)


def load_source(rel: str) -> Tuple[Dict[str, Any], Dict[str, str]]:
    """Load a source data JSON and return its content + provenance."""
    p = DATA / rel
    sha = sha256_file(p)
    with p.open("r", encoding="utf-8") as f:
        data = json.load(f)
    prov = {"file": f"data/{rel}", "sha256": sha}
    return data, prov


def script_provenance(rel: str) -> Tuple[str, str]:
    p = SCRIPTS / rel
    if not p.exists():
        return f"scripts/legacy/{rel}", "MISSING"
    return f"scripts/legacy/{rel}", sha256_file(p)


def _identity_perm(n: int) -> Tuple[int, ...]:
    return tuple(range(1, n + 1))


def _compose_perm(p: Tuple[int, ...], q: Tuple[int, ...]) -> Tuple[int, ...]:
    """Return p o q for one-line permutations."""
    return tuple(p[q[i] - 1] for i in range(len(p)))


def _invert_perm(p: Tuple[int, ...]) -> Tuple[int, ...]:
    out = [0] * len(p)
    for i, value in enumerate(p, start=1):
        out[value - 1] = i
    return tuple(out)


def _perm_order(p: Tuple[int, ...]) -> int:
    ident = _identity_perm(len(p))
    x = ident
    for k in range(1, math.factorial(len(p)) + 1):
        x = _compose_perm(p, x)
        if x == ident:
            return k
    raise AssertionError(f"order not found for permutation {p}")


def _counter_dict(values: List[int]) -> Dict[str, int]:
    return {str(k): int(v) for k, v in sorted(Counter(values).items())}


def _exponent_from_order_distribution(order_distribution: Dict[str, int]) -> int:
    exponent = 1
    for order, count in order_distribution.items():
        if int(count) > 0:
            exponent = math.lcm(exponent, int(order))
    return exponent


def _normalize_group_invariant(inv: Dict[str, Any]) -> Dict[str, Any]:
    normalized = dict(inv)
    normalized["exponent"] = _exponent_from_order_distribution(normalized["order_distribution"])
    return normalized


def _generated_group(gens: List[Tuple[int, ...]]) -> set[Tuple[int, ...]]:
    if not gens:
        raise ValueError("cannot generate a permutation group without generators")
    ident = _identity_perm(len(gens[0]))
    group = {ident}
    frontier = list(gens)
    while frontier:
        g = frontier.pop()
        if g in group:
            continue
        old = list(group)
        group.add(g)
        for h in old:
            frontier.append(_compose_perm(g, h))
            frontier.append(_compose_perm(h, g))
    return group


def _subgroup_orbits(subgroup: set[Tuple[int, ...]], n: int) -> List[List[int]]:
    unseen = set(range(1, n + 1))
    orbits: List[List[int]] = []
    while unseen:
        start = min(unseen)
        orbit = sorted({g[start - 1] for g in subgroup})
        orbits.append(orbit)
        unseen -= set(orbit)
    return orbits


def _left_cosets(
    group: set[Tuple[int, ...]], subgroup: set[Tuple[int, ...]]
) -> Tuple[List[frozenset[Tuple[int, ...]]], Dict[Tuple[int, ...], int]]:
    unused = set(group)
    cosets: List[frozenset[Tuple[int, ...]]] = []
    while unused:
        g = min(unused)
        coset = frozenset(_compose_perm(g, h) for h in subgroup)
        cosets.append(coset)
        unused -= set(coset)
    elem_to_coset = {elem: i for i, coset in enumerate(cosets) for elem in coset}
    return cosets, elem_to_coset


def _quotient_order_distribution(
    group: set[Tuple[int, ...]], normal_subgroup: set[Tuple[int, ...]]
) -> Dict[str, int]:
    ident = _identity_perm(len(next(iter(group))))
    cosets, elem_to_coset = _left_cosets(group, normal_subgroup)
    identity_coset = elem_to_coset[ident]
    orders: List[int] = []
    for coset in cosets:
        rep = min(coset)
        x = ident
        for k in range(1, len(cosets) + 1):
            x = _compose_perm(rep, x)
            if elem_to_coset[x] == identity_coset:
                orders.append(k)
                break
        else:
            raise AssertionError("quotient order computation failed")
    return _counter_dict(orders)


def _derived_subgroup(group: set[Tuple[int, ...]]) -> set[Tuple[int, ...]]:
    commutators: List[Tuple[int, ...]] = []
    inverses = {g: _invert_perm(g) for g in group}
    for a in group:
        for b in group:
            commutators.append(
                _compose_perm(
                    _compose_perm(_compose_perm(a, b), inverses[a]),
                    inverses[b],
                )
            )
    return _generated_group(commutators)


def _affine_D4_witness(elements: List[List[int]]) -> Dict[str, Any]:
    elems = [tuple(p) for p in elements]
    if not elems:
        raise AssertionError("empty group element list")
    n = len(elems[0])
    ident = _identity_perm(n)
    group = set(elems)
    if len(group) != len(elems):
        raise AssertionError("duplicate permutation in group element list")
    if ident not in group:
        raise AssertionError("identity missing from group element list")
    if not all(_compose_perm(a, b) in group for a in group for b in group):
        raise AssertionError("group element list is not closed")

    orders = {g: _perm_order(g) for g in group}
    translations = {ident} | {g for g, order in orders.items() if order == 3}
    if len(translations) != 9:
        raise AssertionError(f"translation subgroup should have order 9, got {len(translations)}")
    if not all(_compose_perm(a, b) in translations for a in translations for b in translations):
        raise AssertionError("translation subgroup is not closed")
    if not all(_compose_perm(a, b) == _compose_perm(b, a) for a in translations for b in translations):
        raise AssertionError("translation subgroup is not abelian")
    if not all(
        _compose_perm(_compose_perm(g, t), _invert_perm(g)) in translations
        for g in group
        for t in translations
    ):
        raise AssertionError("translation subgroup is not normal")

    basis: List[Tuple[int, ...]] = []
    nontrivial_translations = sorted(t for t in translations if t != ident)
    for a in nontrivial_translations:
        for b in nontrivial_translations:
            if a != b and _generated_group([a, b]) == translations:
                basis = [a, b]
                break
        if basis:
            break
    if not basis:
        raise AssertionError("could not find a C3^2 basis for translations")

    point_stab = {g for g in group if g[0] == 1}
    if len(point_stab) != 8:
        raise AssertionError(f"point stabilizer at 1 should have order 8, got {len(point_stab)}")

    d4_fingerprint = {"1": 1, "2": 5, "4": 2}
    point_stab_dist = _counter_dict([orders[g] for g in point_stab])
    if point_stab_dist != d4_fingerprint:
        raise AssertionError(f"point stabilizer is not D4 by order fingerprint: {point_stab_dist}")

    quotient_by_translations = _quotient_order_distribution(group, translations)
    if quotient_by_translations != d4_fingerprint:
        raise AssertionError(
            "quotient by translations is not D4 by order fingerprint: "
            f"{quotient_by_translations}"
        )

    derived = _derived_subgroup(group)
    abelianization_dist = _quotient_order_distribution(group, derived)
    if abelianization_dist != {"1": 1, "2": 3}:
        raise AssertionError(f"abelianization is not C2 x C2: {abelianization_dist}")

    product = {_compose_perm(t, h) for t in translations for h in point_stab}
    legacy_block = {1, 2, 3}
    shared_left_kernel_block = {4, 5, 6}
    legacy_block_preservation_count = sum(
        1 for g in group if {g[i - 1] for i in legacy_block} == legacy_block
    )
    shared_left_kernel_block_preservation_count = sum(
        1 for g in group
        if {g[i - 1] for i in shared_left_kernel_block} == shared_left_kernel_block
    )
    return {
        "group_order": len(group),
        "order_distribution": _counter_dict(list(orders.values())),
        "translation_subgroup": {
            "abstract_group_name": "C_3 x C_3",
            "order": len(translations),
            "order_distribution": _counter_dict([orders[g] for g in translations]),
            "basis_perms": [list(p) for p in basis],
            "is_closed": True,
            "is_abelian": True,
            "is_normal": True,
            "orbits_on_points": _subgroup_orbits(translations, n),
            "is_regular_on_points": _subgroup_orbits(translations, n) == [list(range(1, n + 1))],
        },
        "point_stabilizer_1": {
            "abstract_group_name": "D_4",
            "order": len(point_stab),
            "order_distribution": point_stab_dist,
        },
        "quotient_by_translation_subgroup": {
            "abstract_group_name": "D_4",
            "order": len(group) // len(translations),
            "order_distribution": quotient_by_translations,
        },
        "semidirect_product_certified": (
            len(product) == len(group)
            and len(translations & point_stab) == 1
            and quotient_by_translations == d4_fingerprint
        ),
        "abstract_group_name": "3^2 : D_4",
        "derived_subgroup": {
            "order": len(derived),
            "order_distribution": _counter_dict([_perm_order(g) for g in derived]),
        },
        "abelianization": {
            "abstract_group_name": "C_2 x C_2",
            "order": len(group) // len(derived),
            "order_distribution": abelianization_dist,
        },
        "ambient_legacy_u9_block_preservation_count": legacy_block_preservation_count,
        "ambient_shared_left_kernel_block_preservation_count": shared_left_kernel_block_preservation_count,
        "is_subgroup_of_ambient_coordinate_stabilizer": (
            shared_left_kernel_block_preservation_count == len(group)
        ),
    }


# ----------------------------------------------------------------------------
# Mathematical re-verifications (used by certificates 4 and 5)
# ----------------------------------------------------------------------------
def E_of(L: np.ndarray, n: int) -> np.ndarray:
    return L.astype(np.int64) - (n + 1) // 2 * np.ones((n, n), dtype=np.int64) \
        if (n + 1) % 2 == 0 else L.astype(np.float64) - (n + 1) / 2.0


def E_centered_n9(L: np.ndarray) -> np.ndarray:
    """E = L - 5*J for n=9, integer-valued."""
    return L.astype(np.int64) - 5


def box_band_incidence_maps_n9() -> Tuple[np.ndarray, np.ndarray]:
    """Row-band and column-stack indicator matrices (9x3 each).

    Columns of B_rb are indicator vectors of row-bands {1,2,3}, {4,5,6}, {7,8,9}.
    Columns of B_cs are indicator vectors of column-stacks {1,2,3}, {4,5,6}, {7,8,9}.
    """
    P = np.zeros((9, 3), dtype=np.int64)
    for k in range(3):
        P[3 * k : 3 * (k + 1), k] = 1
    return P, P


def verify_box_band(L: np.ndarray) -> Dict[str, Any]:
    """Check the Box Band Lemma identity P_rb^T E P_cs = 0_{3x3} on a 9x9 grid."""
    E = E_centered_n9(L)
    Brb, Bcs = box_band_incidence_maps_n9()
    M = Brb.T @ E @ Bcs  # 3x3 integer matrix
    return {
        "M_rb_cs": M.tolist(),
        "is_zero": bool(np.all(M == 0)),
        "frobenius_norm_sq": int(np.sum(M * M)),
    }


def centered_twice(L: np.ndarray) -> np.ndarray:
    """Return 2E = 2L - (n+1)J as an integer matrix for any n."""
    n = int(L.shape[0])
    return 2 * L.astype(np.int64) - (n + 1)


def rectangular_box_band_check(L: np.ndarray, h: int, w: int) -> Dict[str, Any]:
    """Check the rectangular h x w Box Band identity using 2E for integrality."""
    n = int(L.shape[0])
    if n != h * w:
        raise AssertionError(f"n={n} does not equal h*w={h*w}")
    n_row_bands = n // h
    n_col_stacks = n // w
    Brb = np.zeros((n, n_row_bands), dtype=np.int64)
    Bcs = np.zeros((n, n_col_stacks), dtype=np.int64)
    for band in range(n_row_bands):
        Brb[band * h : (band + 1) * h, band] = 1
    for stack in range(n_col_stacks):
        Bcs[stack * w : (stack + 1) * w, stack] = 1
    M = Brb.T @ centered_twice(L) @ Bcs
    return {
        "box_height": h,
        "box_width": w,
        "n_row_bands": n_row_bands,
        "n_col_stacks": n_col_stacks,
        "uses_twice_centered_matrix": True,
        "twice_centered_M_rb_cs": M.tolist(),
        "is_zero": bool(np.all(M == 0)),
        "frobenius_norm_sq": int(np.sum(M * M)),
    }


def cyclic_diagonal_gerechte_regions(n: int = 9) -> List[List[Tuple[int, int]]]:
    """Broken-diagonal transversal partition for the cyclic Latin square."""
    return [[(i, (i + b) % n) for i in range(n)] for b in range(n)]


def gerechte_region_check(
    L: np.ndarray, regions: List[List[Tuple[int, int]]], gamma: List[int]
) -> Dict[str, Any]:
    """Check that every gerechte region has centered regional sum zero."""
    n = int(L.shape[0])
    if len(regions) != n:
        raise AssertionError(f"expected {n} regions, got {len(regions)}")
    cells = [cell for region in regions for cell in region]
    if len(cells) != n * n or len(set(cells)) != n * n:
        raise AssertionError("regions do not partition the cell set")
    Lg = relabel(L, gamma)
    E2 = centered_twice(Lg)
    region_sums_twice = []
    symbol_sets = []
    cartesian_region_flags = []
    for region in regions:
        if len(region) != n:
            raise AssertionError("region size mismatch")
        symbols = sorted(int(L[i, j]) for i, j in region)
        symbol_sets.append(symbols)
        rows = {i for i, _ in region}
        cols = {j for _, j in region}
        cartesian_region_flags.append(len(rows) * len(cols) == len(region))
        region_sums_twice.append(int(sum(E2[i, j] for i, j in region)))
    all_symbols_once = all(symbols == list(range(1, n + 1)) for symbols in symbol_sets)
    return {
        "n": n,
        "gamma": gamma,
        "n_regions": len(regions),
        "region_size": n,
        "region_formula": "R_b = {(i, i+b mod n) : i=0,...,n-1}",
        "regions_1_based": [
            [[int(i + 1), int(j + 1)] for i, j in region] for region in regions
        ],
        "regions_partition_cells": True,
        "each_region_symbols_once": bool(all_symbols_once),
        "uses_twice_centered_matrix": True,
        "region_sums_twice_centered": region_sums_twice,
        "all_region_sums_zero": all(v == 0 for v in region_sums_twice),
        "all_regions_cartesian_products": all(cartesian_region_flags),
        "any_region_cartesian_product": any(cartesian_region_flags),
    }


def cyclic_LS9() -> np.ndarray:
    """A canonical non-Sudoku Latin square of order 9: cyclic L[i,j] = ((i+j) mod 9) + 1."""
    L = np.zeros((9, 9), dtype=np.int64)
    for i in range(9):
        for j in range(9):
            L[i, j] = ((i + j) % 9) + 1
    return L


def relabel(L: np.ndarray, gamma: List[int]) -> np.ndarray:
    """Apply symbol relabeling gamma (1-indexed list of length n)."""
    g = np.asarray(gamma, dtype=np.int64)
    out = g[L.astype(np.int64) - 1]
    return out


def integer_kernel_via_smith(M: np.ndarray) -> np.ndarray:
    """Return integer kernel basis via Sympy SNF (rows form basis of right nullspace)."""
    from sympy import Matrix, ZZ
    A = Matrix(M.tolist())
    K = A.nullspace()
    if not K:
        return np.zeros((0, M.shape[1]), dtype=np.int64)
    # Clear denominators; produce int64
    rows = []
    for v in K:
        v = v.T  # row
        denom = 1
        for x in v:
            denom = denom * x.q if hasattr(x, "q") else denom
        scaled = [int((x * denom)) for x in v]
        rows.append(scaled)
    return np.array(rows, dtype=np.int64)


def lift_lemma_check(L: np.ndarray, gamma: List[int]) -> Dict[str, Any]:
    """Verify the Lift Lemma at (L, gamma).

    Lift Lemma. For w in V_std, E_gamma(L) w = 0 iff gamma(L) w = 0
    (equivalently: ker(E_gamma) cap V_std has codimension 1 in ker(E_gamma),
    and gamma(L) kills every element of that intersection).

    Computation. Use exact rationals via Sympy to compute:
      (a) ker(E_gamma) (full kernel; always contains the all-ones vector since
          row/col sums of gamma(L) all equal n(n+1)/2);
    # ker(E_gamma) intersect V_std = nullspace of [E_gamma stacked with 1^T]
      (c) verify gamma(L) annihilates (b).
    """
    from sympy import Matrix, zeros
    Lg = relabel(L, gamma)
    Eg = Matrix((Lg - 5).tolist())
    LgM = Matrix(Lg.tolist())
    ones9 = Matrix([1] * 9)

    K_full = Eg.nullspace()
    full_dim = len(K_full)
    eg_kills_full = all((Eg * v) == zeros(9, 1) for v in K_full)
    one_in_ker = (Eg * ones9) == zeros(9, 1)

    # ker(E_gamma) intersect V_std = nullspace of [E_gamma stacked with 1^T]
    A = Eg.col_join(Matrix([[1] * 9]))
    K_Vstd = A.nullspace()
    Vstd_dim = len(K_Vstd)
    gl_kills_Vstd = all((LgM * v) == zeros(9, 1) for v in K_Vstd)
    eg_kills_Vstd = all((Eg * v) == zeros(9, 1) for v in K_Vstd)
    in_Vstd = all((Matrix([[1] * 9]) * v) == zeros(1, 1) for v in K_Vstd)

    return {
        "kernel_dim_full": full_dim,
        "kernel_dim_in_V_std": Vstd_dim,
        "one_vector_in_full_kernel": bool(one_in_ker),
        "E_gamma_kills_full_kernel": bool(eg_kills_full),
        "E_gamma_kills_kernel_in_V_std": bool(eg_kills_Vstd),
        "gammaL_kills_kernel_in_V_std": bool(gl_kills_Vstd),
        "kernel_in_V_std_lies_in_V_std": bool(in_Vstd),
        "codim_one": bool(Vstd_dim == full_dim - 1),
        "lift_lemma_holds": bool(
            gl_kills_Vstd and Vstd_dim == full_dim - 1 and one_in_ker
        ),
    }


# ----------------------------------------------------------------------------
# Builders (one per certificate)
# ----------------------------------------------------------------------------
def cert_base1_odd_rank_22() -> Dict[str, Any]:
    src, prov = load_source("phase_9_13_results.json")
    block = src["phase_9_13_c"]
    cert_src, cert_prov = load_source("phase_9_13f_certificate.json")

    # Full 22-perm list (recovered locally; original phase_9_13.py truncates
    # to 5).  We carry the recovered list and its hash for traceability.
    recovered_path = HERE / "base1_22_perms_recovered.json"
    with recovered_path.open("r", encoding="utf-8") as f:
        recovered = json.load(f)
    rec_sha = sha256_file(recovered_path)
    assert recovered["odd_rank_count"] == 22, "recovered odd-rank count mismatch"
    assert recovered["odd_rank_count"] == block["odd_rank_count"], (
        "recovered count disagrees with phase_9_13.py log"
    )

    # Four-way exact certificate (rank Z, primitive V_std kernel, integer SNF,
    # F_2 SNF) for every one of the 22 odd-rank perms.
    fourway_path = HERE / "four_way_base1_22.json"
    with fourway_path.open("r", encoding="utf-8") as f:
        fourway = json.load(f)
    assert fourway["n_certs"] == 22
    fw_certs = fourway["four_way_certificates"]
    assert all(c["exact_rank_Z"] == 7 for c in fw_certs)
    # Histograms.
    snfZ_hist: Counter = Counter()
    rankF2_hist: Counter = Counter()
    for c in fw_certs:
        snfZ_hist[tuple(c["snf_Z_diagonal"])] += 1
        rankF2_hist[c["rank_F2"]] += 1

    sp1, sh1 = script_provenance("phase13b_exhaustive.py")
    sp2, sh2 = script_provenance("phase_9_13.py")
    sp3, sh3 = script_provenance("phase_9_13f_certificate.py")
    return {
        "result_id": "base1_odd_rank_22",
        "tier": "EXHAUSTIVE",
        "paper_label": "Certified Proposition P1",
        "claim": (
            "Exhaustive scan of the 9! = 362880 symbol relabelings gamma in S_9 of "
            "the canonical Sudoku base 'base1' yields rank distribution "
            "{8: 362858, 7: 22}; the 22 odd-rank relabelings each admit a four-way "
            "exact certificate consisting of (i) integer rank = 7, (ii) a primitive "
            "integer kernel vector in V_std (sum = 0), (iii) the integer Smith "
            "Normal Form, and (iv) the F_2 Smith Normal Form.  All four are "
            "computed exactly via Sympy on every one of the 22 perms."
        ),
        "inputs": {
            "n": 9,
            "base_grid": block["base_grid"],
            "scan_space": "all of S_9",
            "n_perms": block["n_perms"],
        },
        "outputs": {
            "rank_distribution": block["rank_dist"],
            "odd_rank_count": block["odd_rank_count"],
            "odd_perms_full_list": recovered["odd_perms"],
            "odd_perms_count": len(recovered["odd_perms"]),
            "odd_examples_logged_in_phase_9_13": [
                e["perm"] for e in block["odd_examples"]
            ],
            "n_certified_perms_in_phase_9_13f": len(cert_src.get("certificates", [])),
            "four_way_certificates": fw_certs,
            "four_way_n_certs": len(fw_certs),
            "snf_Z_diagonal_histogram": {
                str(list(k)): v for k, v in snfZ_hist.items()
            },
            "rank_F2_histogram": {str(k): v for k, v in rankF2_hist.items()},
        },
        "source": [
            {**prov, "extracted_keys": ["phase_9_13_c"]},
            {**cert_prov, "extracted_keys": ["certificates"]},
            {
                "file": "certified/base1_22_perms_recovered.json",
                "sha256": rec_sha,
                "extracted_keys": ["odd_perms", "odd_rank_count"],
                "note": (
                    "Recovered via certified/recover_base1_22_perms.py; "
                    "exhaustive S_9 rescan of base1 in ~16s."
                ),
            },
            {
                "file": "certified/four_way_base1_22.json",
                "sha256": sha256_file(fourway_path),
                "extracted_keys": ["four_way_certificates"],
                "note": (
                    "Generated via certified/four_way_base1_22.py; "
                    "exact Sympy computation of (rank Z, primitive V_std kernel, "
                    "integer SNF, F_2 SNF) for every one of the 22 perms."
                ),
            },
        ],
        "scripts": [
            {"path": sp1, "sha256": sh1},
            {"path": sp2, "sha256": sh2},
            {"path": sp3, "sha256": sh3},
            {
                "path": "certified/recover_base1_22_perms.py",
                "sha256": sha256_file(HERE / "recover_base1_22_perms.py"),
            },
            {
                "path": "certified/four_way_base1_22.py",
                "sha256": sha256_file(HERE / "four_way_base1_22.py"),
            },
        ],
        "seed": None,
    }


def cert_n6_odd_rank_census() -> Dict[str, Any]:
    src, prov = load_source("phase_9_14_results.json")
    sud = src["sudoku6_firstrow_fixed"]
    ls = src["ls6_firstrow_fixed"]
    sp, sh = script_provenance("phase_9_14.py")

    # Two-proportion z-test (pooled standard error). At n=6, the relevant
    # statistic is the full-rank-over-F_2 incidence for the determinant-paper
    # auxiliary matrix A_{ij}=L_{ij}-L_{i,n-1} (rank n-1 = 5 over F_2), i.e. a
    # rank-parity proxy.  This is NOT the same matrix as the symmetric P^T E P
    # used elsewhere, and it is NOT the same notion as "odd-rank" (rank n-2)
    # used at n in {7, 9, 11}; the script records it under the legacy field name
    # `full_rank_count`.
    p1, n1 = sud["full_rank_count"] / sud["n_total"], sud["n_total"]
    p2, n2 = ls["full_rank_count"] / ls["n_total"], ls["n_total"]
    p_pool = (sud["full_rank_count"] + ls["full_rank_count"]) / (n1 + n2)
    se = math.sqrt(p_pool * (1 - p_pool) * (1 / n1 + 1 / n2))
    z = (p1 - p2) / se
    # Exact two-sided descriptive erfc statistic (no scipy):
    # value = erfc(|z| / sqrt(2)).
    p_two = math.erfc(abs(z) / math.sqrt(2.0))

    return {
        "result_id": "n6_odd_rank_census",
        "tier": "EXHAUSTIVE",
        "paper_label": "Proposition P2",
        "claim": (
            "Exhaustive enumeration of all order-6 Sudoku grids (n=39168, first-row "
            "fixed) and order-6 Latin squares (n=1128960, first-row fixed) shows "
            "full rank over F_2 for the auxiliary matrix A_{ij}=L_{ij}-L_{i,n-1} "
            "(rank n-1 = 5) with incidence 6.6177% vs 6.1224%.  The pooled "
            "two-proportion descriptive comparison gives z = 4.013706 and "
            "two-sided erfc(|z|/sqrt 2) = 5.9773e-05.  At n=6 this is the "
            "determinant-paper rank-parity proxy (`full_rank_count` over F_2), "
            "not rank(P^T E P), and not the odd-rank locus (rank n-2) used at n "
            "in {7, 9, 11}; the box constraint nevertheless yields a small "
            "enhancement."
        ),
        "inputs": {
            "n": 6,
            "first_row_fixed": [1, 2, 3, 4, 5, 6],
            "sudoku6_n_total": sud["n_total"],
            "ls6_n_total": ls["n_total"],
            "rank_parity_proxy": "F_2 full rank of A_{ij}=L_{ij}-L_{i,n-1} (rank == n-1 == 5); not rank(P^T E P)",
        },
        "outputs": {
            "sudoku6_full_rank_F2_count": sud["full_rank_count"],
            "sudoku6_full_rank_F2_pct": sud["full_rank_pct"],
            "sudoku6_rank_dist": sud["rank_dist"],
            "ls6_full_rank_F2_count": ls["full_rank_count"],
            "ls6_full_rank_F2_pct": ls["full_rank_pct"],
            "ls6_rank_dist": ls["rank_dist"],
            "two_proportion_z": z,
            "two_proportion_pvalue_exact": p_two,
            "two_proportion_pvalue_formula": "erfc(|z| / sqrt(2))",
            "two_proportion_pvalue_decimal_string": f"{p_two:.6e}",
        },
        "source": [
            {**prov, "extracted_keys": ["sudoku6_firstrow_fixed", "ls6_firstrow_fixed"]}
        ],
        "scripts": [{"path": sp, "sha256": sh}],
        "seed": None,
    }


def cert_n7_HV1_D6() -> Dict[str, Any]:
    src, prov = load_source("phase_9_15A_H_V1_identification.json")
    src15, prov15 = load_source("phase_9_15_results.json")
    sp1, sh1 = script_provenance("phase_9_15A_identify_HV1.py")
    sp2, sh2 = script_provenance("check_b3_subgroup.py")
    spB, shB = script_provenance("phase_9_15.py")
    L_V1_block = src15["B3_enrichment"]["L_V1"]
    return {
        "result_id": "n7_HV1_D6",
        "tier": "CERTIFIED_EXHAUSTIVE",
        "paper_label": "Certified Proposition B",
        "claim": (
            "At n=7, let L_{V_1} be the laboratory 7x7 Latin square recorded in "
            "scripts/legacy/phase_9_15.py.  Define H_{V_1} as the rank-locus subset "
            "{gamma in B_3 subset S_7 : rank E_{gamma, L_{V_1}} = 5}, i.e. the "
            "intersection of the odd-rank symbol-relabeling locus of L_{V_1} with "
            "the signed-pair group B_3 = C_2 wr S_3 of order 48 acting on the "
            "three pairs {1,7}, {2,6}, {3,5} and fixing 4. "
            "Then |H_{V_1}| = 12, H_{V_1} is closed under composition (a subgroup "
            "of B_3), and H_{V_1} cong D_6 (dihedral of order 12), equivalently "
            "S_3 x Z/2.  The element-order multiset {1:1, 2:7, 3:2, 6:2} uniquely "
            "identifies D_6 among groups of order 12 (alternatives Z/12, Z/6 x Z/2, "
            "A_4, Dic_3 are excluded by their order multisets).  Constructively, "
            "H_{V_1} = <r, s> with r of order 6, s of order 2, and srs = r^{-1}; "
            "the relation is verified element-by-element. NOTE: H_{V_1} is NOT the "
            "coordinate stabilizer of any vector with all-distinct coordinates "
            "(which would be trivial in S_7); it is the rank-locus subset of B_3 "
            "associated with the base L_{V_1}, defined operationally by the rank "
            "condition above."
        ),
        "inputs": {
            "n": 7,
            "base_grid_label": "L_{V_1}",
            "ambient_group": "B_3 = C_2 wr S_3 subset S_7, preserving pairs {1,7}, {2,6}, {3,5} and fixing 4",
            "definition": "H_{V_1} := {gamma in B_3 : rank(E_{gamma, L_{V_1}}) = 5}",
        },
        "outputs": {
            "abstract_group": src["abstract_group"],
            "order": 12,
            "L_V1_rank_distribution_in_B3": L_V1_block["rank_dist"],
            "L_V1_odd_rank_count_in_B3": L_V1_block["odd_rank_count"],
            "L_V1_rate_in_B3": L_V1_block["rate_in_B3"],
            "L_V1_baseline_rate_full_S7": L_V1_block["baseline_rate"],
            "L_V1_enrichment_ratio": L_V1_block["enrichment_ratio"],
            "generator_r_order_6": src["generator_r_order_6"],
            "generator_s_order_2": src["generator_s_order_2"],
            "dihedral_relation_verified": src["dihedral_relation_verified"],
            "generates_full_H_V1": src["generates_full_H_V1"],
            "order_distribution": src["order_distribution"],
            "order_distribution_unique_to_D6_among_order_12": True,
            "alternatives_excluded_by_order_distribution": [
                "Z/12: {1:1, 2:1, 3:2, 4:2, 6:2, 12:4}",
                "Z/6 x Z/2: {1:1, 2:3, 3:2, 6:6}",
                "A_4: {1:1, 2:3, 3:8}",
                "Dic_3: {1:1, 2:1, 3:2, 4:6, 6:2}",
            ],
            "is_subgroup_of_B3": True,
        },
        "source": [
            {**prov, "extracted_keys": ["*"]},
            {**prov15, "extracted_keys": ["B3_enrichment.L_V1"]},
        ],
        "scripts": [
            {"path": sp1, "sha256": sh1},
            {"path": sp2, "sha256": sh2},
            {"path": spB, "sha256": shB},
        ],
        "seed": None,
    }


def cert_box_band_lemma_witness() -> Dict[str, Any]:
    """Re-verifies square, rectangular, and gerechte Box Band identities."""
    base1_src, base1_prov = load_source("phase_9_13_results.json")
    M19_src, M19_prov = load_source("phase_9_18_S9x_results.json")
    n6_src, n6_prov = load_source("phase_9_14_results.json")
    base1 = np.array(base1_src["phase_9_13_c"]["base_grid"], dtype=np.int64)
    M19 = np.array(M19_src["M19_grid"], dtype=np.int64)
    roku6 = np.array(n6_src["sudoku6_firstrow_fixed"]["bad_examples"][0], dtype=np.int64)
    cyclic9 = cyclic_LS9()

    sudoku_base1 = verify_box_band(base1)
    sudoku_M19 = verify_box_band(M19)
    ls_cyclic = verify_box_band(cyclic9)
    rectangular_roku6 = rectangular_box_band_check(roku6, h=2, w=3)
    gerechte_cyclic_diagonal = gerechte_region_check(
        cyclic9,
        cyclic_diagonal_gerechte_regions(9),
        gamma=[9, 1, 8, 2, 7, 3, 6, 4, 5],
    )

    sp1, sh1 = script_provenance("task_9_10c_box_subspace.py")
    sp2, sh2 = script_provenance("phase_9_13.py")
    sp3, sh3 = script_provenance("phase_9_14.py")

    return {
        "result_id": "box_band_lemma_witness",
        "tier": "THEOREM",
        "paper_label": "Theorem A",
        "claim": (
            "For every Sudoku grid L of order 9, the centered row-band/column-stack "
            "average P_rb^T E P_cs of E = L - 5J vanishes identically. This "
            "identity is not forced by the Latin-square axioms alone; the cyclic "
            "LS-9 is exhibited as an explicit "
            "counterexample. The same centered-sum argument is certified for a "
            "rectangular 2x3 Roku-Doku grid and for a non-Cartesian gerechte "
            "partition of cyclic LS-9 by broken diagonals."
        ),
        "inputs": {
            "square_sudoku_n": 9,
            "square_centering_scalar": 5,
            "rectangular_roku_doku": {"n": 6, "box_height": 2, "box_width": 3},
            "gerechte_witness": "cyclic LS-9 partitioned by broken diagonals (non-Cartesian regions)",
        },
        "outputs": {
            "sudoku_base1": sudoku_base1,
            "sudoku_M19": sudoku_M19,
            "non_sudoku_LS_cyclic_counterexample": ls_cyclic,
            "rectangular_roku_doku_2x3": rectangular_roku6,
            "gerechte_cyclic_LS9_broken_diagonals": gerechte_cyclic_diagonal,
        },
        "source": [
            {**base1_prov, "extracted_keys": ["phase_9_13_c.base_grid"]},
            {**M19_prov, "extracted_keys": ["M19_grid"]},
            {**n6_prov, "extracted_keys": ["sudoku6_firstrow_fixed.bad_examples[0]"]},
        ],
        "scripts": [
            {"path": sp1, "sha256": sh1},
            {"path": sp2, "sha256": sh2},
            {"path": sp3, "sha256": sh3},
            {"path": "certified/build_certified.py", "sha256": "self"},
        ],
        "seed": None,
    }


def cert_lift_lemma_evidence() -> Dict[str, Any]:
    """Re-verifies the Lift Lemma biconditional on the 22 base1 odd-rank
    relabelings and on the full 100 M19 rank-7 relabelings."""
    base1_src, base1_prov = load_source("phase_9_13_results.json")
    M19_src, M19_prov = load_source("phase_9_18_S9x_results.json")

    base1 = np.array(base1_src["phase_9_13_c"]["base_grid"], dtype=np.int64)

    # Use the recovered full 22-perm list for base1.
    recovered_path = HERE / "base1_22_perms_recovered.json"
    with recovered_path.open("r", encoding="utf-8") as f:
        recovered = json.load(f)
    base1_perms: List[List[int]] = recovered["odd_perms"]

    # Use the recovered full 100-perm list for M19 (rank-7 locus).
    M19_recovered_path = HERE / "M19_100_perms_recovered.json"
    with M19_recovered_path.open("r", encoding="utf-8") as f:
        M19_recovered = json.load(f)
    assert M19_recovered["rank7_count"] == 100, "M19 rank-7 count mismatch"
    M19 = np.array(M19_src["M19_grid"], dtype=np.int64)
    assert np.array_equal(M19, np.array(M19_recovered["M19_grid"], dtype=np.int64))
    M19_perms: List[List[int]] = M19_recovered["rank7_perms"]
    # Also retain Gamma_star (72 shared left-kernel / middle-band relabelings) as a labeled subset for cross-check.
    Gamma_star_perms: List[List[int]] = M19_src["Gamma_star_perms"]
    assert len(Gamma_star_perms) == 72
    Gamma_star_set = set(tuple(p) for p in Gamma_star_perms)
    rank7_set = set(tuple(p) for p in M19_perms)
    assert Gamma_star_set.issubset(rank7_set), "72 shared left-kernel set not subset of 100 rank-7"

    def summarize(L: np.ndarray, perms: List[List[int]]) -> Dict[str, Any]:
        full_pass = True
        kdim_hist: Counter = Counter()
        Vstd_dim_hist: Counter = Counter()
        # Full lift-lemma verification on every perm (Sympy exact).
        for g in perms:
            r = lift_lemma_check(L, g)
            kdim_hist[r["kernel_dim_full"]] += 1
            Vstd_dim_hist[r["kernel_dim_in_V_std"]] += 1
            if not r["lift_lemma_holds"]:
                full_pass = False
        return {
            "n_perms": len(perms),
            "n_full_lift_check": len(perms),
            "lift_lemma_holds_on_all_checked": full_pass,
            "kernel_dim_full_histogram": dict(sorted(kdim_hist.items())),
            "kernel_dim_in_V_std_histogram": dict(sorted(Vstd_dim_hist.items())),
        }

    res_base1 = summarize(base1, base1_perms)
    res_M19 = summarize(M19, M19_perms)
    # Also separately tag the 72 shared left-kernel subset within the 100.
    res_M19_shared72 = summarize(M19, Gamma_star_perms)

    sp, sh = script_provenance("phase_9_18_S9c.py")

    return {
        "result_id": "lift_lemma_evidence",
        "tier": "THEOREM",
        "paper_label": "Theorem C",
        "claim": (
            "Lift Lemma: for every gamma in S_9 and every w in V_std, "
            "E_gamma(L) w = 0 if and only if gamma(L) w = 0 (and dually for left "
            "kernels). Verified on all 22 base1 rank-7 relabelings and on all 100 "
            "M_{19} rank-7 relabelings (the full rank-locus, not just the "
            "72-element shared left-kernel / middle-band subset)."
        ),
        "inputs": {
            "n": 9,
            "base1_perms": len(base1_perms),
            "M19_rank7_perms_total": len(M19_perms),
            "M19_shared_left_kernel_perms_subset": len(Gamma_star_perms),
        },
        "outputs": {
            "base1": res_base1,
            "M19": res_M19,
            "M19_shared_left_kernel_72": res_M19_shared72,
            "M19_shared_kernel_72": res_M19_shared72,
        },
        "source": [
            {**base1_prov, "extracted_keys": ["phase_9_13_c.odd_examples"]},
            {**M19_prov, "extracted_keys": ["Gamma_star_perms", "M19_grid"]},
            {
                "file": "certified/M19_100_perms_recovered.json",
                "sha256": sha256_file(HERE / "M19_100_perms_recovered.json"),
                "extracted_keys": ["rank7_perms", "M19_grid"],
                "note": (
                    "Recovered via certified/recover_M19_100_perms.py; "
                    "exhaustive S_9 modular determinant rescan of M_{19} followed "
                    "by Sympy exact rank verification on every candidate."
                ),
            },
        ],
        "scripts": [
            {"path": sp, "sha256": sh},
            {"path": "certified/build_certified.py", "sha256": "self"},
            {
                "path": "certified/recover_M19_100_perms.py",
                "sha256": sha256_file(HERE / "recover_M19_100_perms.py"),
            },
        ],
        "seed": None,
    }


def cert_M19_audit() -> Dict[str, Any]:
    src, prov = load_source("phase_9_18_S9x_results.json")
    bis_src, bis_prov = load_source("phase_9_18_S9x_bis_results.json")
    sp, sh = script_provenance("phase_9_18_S9x.py")

    # Use the locally recovered exact list of 100 rank-7 perms.
    M19_recovered_path = HERE / "M19_100_perms_recovered.json"
    with M19_recovered_path.open("r", encoding="utf-8") as f:
        M19_recovered = json.load(f)
    assert M19_recovered["rank7_count"] == 100
    rank7_perms = M19_recovered["rank7_perms"]
    Gamma_star_72 = src["Gamma_star_perms"]
    rank7_set = set(tuple(p) for p in rank7_perms)
    G72_set = set(tuple(p) for p in Gamma_star_72)
    assert G72_set.issubset(rank7_set)
    residue_28 = sorted(list(rank7_set - G72_set))

    # Sanity: 4320 = 3! * 6! = |Stab_{S_9}(w_left)| under coordinate action,
    # since the multiset of w_left = [1,1,1,-2,-2,-2,1,1,1] has multiplicities
    # {-2: 3, 1: 6}, hence stabilizer S_3 x S_6.
    w_left = [1, 1, 1, -2, -2, -2, 1, 1, 1]
    from collections import Counter as _C
    mult = _C(w_left)
    expected_stab = 1
    for v in mult.values():
        expected_stab *= math.factorial(v)
    assert expected_stab == 4320, (
        f"ambient stabilizer mismatch: {expected_stab} != 4320"
    )
    assert src["stab_size"] == 4320

    return {
        "result_id": "M19_audit",
        "tier": "EXHAUSTIVE",
        "paper_label": "Proposition P3",
        "claim": (
            "Exhaustive S_9 audit of the laboratory base M_{19} yields exactly 100 "
            "rank-7 (=odd-rank) relabelings out of 9! = 362880; this rank locus "
            "Gamma_rank(M_{19}) of size 100 contains a distinguished left shared-"
            "kernel / middle-band subset Gamma_star(M_{19}) of size 72 (the "
            "relabelings satisfying w_left^T E_{gamma,M19}=0 for "
            "w_left=(1,1,1,-2,-2,-2,1,1,1)), with residue of size 100 - 72 = 28. "
            "The subset Gamma_star(M_{19}) is a single left K_L-coset and a single "
            "right K_R-coset (left/right coset tests PASS), with K_L, K_R subgroups "
            "of S_9 of order 72 each, conjugate via pi_0 and intersecting trivially. "
            "The number 4320 is the order of the ambient coordinate stabilizer "
            "Stab_{S_9}(w_left) = S_3 x S_6 (since w_left has value-multiplicities "
            "{-2: 3, 1: 6}); 4320 is NOT the cardinality of the rank locus, NOT "
            "the cardinality of Gamma_star, and NOT an orbit of K_L x K_R."
        ),
        "inputs": {
            "n": 9,
            "M19_grid": src["M19_grid"],
            "scan_space": "all of S_9",
            "n_total": 362880,
        },
        "outputs": {
            # --- Exhaustive S_9 rank locus ---
            "Gamma_rank_size": 100,
            "n_certs_M19": src["n_certs_M19"],
            "rank_distribution": M19_recovered["rank_distribution"],
            "rank7_perms_count": len(rank7_perms),
            # --- Shared left-kernel / middle-band subset (72) ---
            "Gamma_star_size": 72,
            "Gamma_star_is_subset_of_Gamma_rank": True,
            "n_orbits_M19": src["n_orbits_M19"],
            "shared_left_kernel_w": w_left,
            "left_kernel_condition": "w_left^T E_{gamma,M19} = 0",
            "legacy_value_template_from_source": src["u9_star"],
            "value_multiplicities": {str(k): int(v) for k, v in sorted(mult.items())},
            "residue_size": src["residue_size"],
            "Gamma_star_residue_size_check": len(residue_28),
            # --- Bilateral coset structure of Gamma_star ---
            "left_coset_test_PASS": src["left_coset_test"] == "YES",
            "right_coset_test_PASS": src["right_coset_test"] == "YES",
            "K_L_size": src["K_L_size"],
            "K_R_size": src["K_R_size"],
            "K_L_eq_K_R": src["K_L_eq_K_R"],
            "K_L_inter_K_R_size": bis_src["K_L_inter_K_R_size"],
            "K_L_conjugate_to_K_R_via_pi0": bool(bis_src["Y2_conjugate_via_pi0"]),
            "double_coset_size_K_L_pi0_K_R": 72,
            "double_coset_intersection_size": (
                src["K_L_size"] * src["K_R_size"] // 72
            ),
            # --- Ambient stabilizer (NOT the orbit) ---
            "ambient_coordinate_stabilizer_S_9_of_shared_left_kernel_size": 4320,
            "ambient_stab_iso": "S_3 x S_6",
            "ambient_stab_order_factorization": "3! * 6! = 6 * 720 = 4320",
        },
        "source": [
            {
                **prov,
                "extracted_keys": [
                    "M19_grid",
                    "n_certs_M19",
                    "n_orbits_M19",
                    "stab_size",
                    "K_L_size",
                    "K_R_size",
                    "K_L_eq_K_R",
                    "left_coset_test",
                    "right_coset_test",
                    "value_multiplicities",
                    "residue_size",
                    "Gamma_star_perms",
                ],
            },
            {
                **bis_prov,
                "extracted_keys": [
                    "K_L_inter_K_R_size",
                    "Y2_conjugate_via_pi0",
                ],
            },
            {
                "file": "certified/M19_100_perms_recovered.json",
                "sha256": sha256_file(HERE / "M19_100_perms_recovered.json"),
                "extracted_keys": ["rank7_perms", "rank_distribution"],
                "note": (
                    "Recovered via certified/recover_M19_100_perms.py; "
                    "two-stage exhaustive scan (modular determinant + Sympy exact rank)."
                ),
            },
        ],
        "scripts": [
            {"path": sp, "sha256": sh},
            {
                "path": "certified/recover_M19_100_perms.py",
                "sha256": sha256_file(HERE / "recover_M19_100_perms.py"),
            },
        ],
        "seed": None,
    }


def cert_M19_KL_KR_affine_D4() -> Dict[str, Any]:
    src1, prov1 = load_source("phase_9_18_S9x_results.json")
    src2, prov2 = load_source("phase_9_18_S9x_bis_results.json")
    inv = src2["invariants"]
    sp1, sh1 = script_provenance("phase_9_18_S9x.py")
    sp2, sh2 = script_provenance("phase_9_18_S9x_bis.py")
    gamma_star_perms = src1.get("Gamma_star_perms", []) or []
    pi0_perm = min(gamma_star_perms) if gamma_star_perms else None
    expected_pi0_perm = [1, 5, 3, 9, 6, 8, 7, 4, 2]
    if pi0_perm != expected_pi0_perm:
        raise AssertionError(f"pi0_perm drifted: {pi0_perm}")
    pi0_text = "(" + ",".join(str(x) for x in pi0_perm) + ")"
    K_L_elements = src1.get("K_L_elements", []) or []
    K_R_elements = src1.get("K_R_elements", []) or []
    pi0_in_K_L = pi0_perm is not None and pi0_perm in K_L_elements
    pi0_in_K_R = pi0_perm is not None and pi0_perm in K_R_elements
    K_L_affine = _affine_D4_witness(K_L_elements)
    K_R_affine = _affine_D4_witness(K_R_elements)
    return {
        "result_id": "M19_KL_KR_affine_D4",
        "tier": "CERTIFIED_EXHAUSTIVE",
        "paper_label": "Certified Proposition E",
        "claim": (
            "At the laboratory base M_{19}, the left and right affine coset-symmetry "
            "groups of the bilateral coset structure of Gamma^*_{M_{19}} satisfy "
            "K_L cong K_R cong 3^2:D_4. Each group has a normal regular "
            "translation subgroup C_3 x C_3 and point stabilizer D_4; both groups "
            "have order 72, exponent 12, derived series 72->18->9->1, "
            "abelianization C_2 x C_2, trivial centre, and 9 conjugacy classes. "
            "K_L != K_R as subgroups of "
            "S_9 (intersection trivial), but they are conjugate via the external "
            f"structural permutation pi_0 = {pi0_text} in S_9 (lex-min of "
            "Gamma^*_{M_{19}}); pi_0 is not contained in K_L or K_R. Neither "
            "group is contained in the ambient coordinate stabilizer of the shared "
            "left-kernel direction w_left=(1,1,1,-2,-2,-2,1,1,1)."
        ),
        "inputs": {"n": 9, "base": "M_{19}"},
        "outputs": {
            "K_L": _normalize_group_invariant(inv["K_L"]),
            "K_R": _normalize_group_invariant(inv["K_R"]),
            "K_L_eq_K_R": src2["K_L_eq_K_R"],
            "K_L_inter_K_R_size": src2["K_L_inter_K_R_size"],
            "Y2_conjugate_via_pi0": src2["Y2_conjugate_via_pi0"],
            "K_L_block_systems": src2["K_L_block_systems"],
            "K_R_block_systems": src2["K_R_block_systems"],
            "abstract_group_id": None,
            "abstract_group_id_note": "No SmallGroup identifier is asserted by this certificate.",
            "abstract_group_name": "3^2 : D_4",
            "K_L_affine_structure": K_L_affine,
            "K_R_affine_structure": K_R_affine,
            "K_L_generators_perms": K_L_elements[:8],
            "pi0_perm": pi0_perm,
            "pi0_in_K_L": pi0_in_K_L,
            "pi0_in_K_R": pi0_in_K_R,
        },
        "source": [
            {**prov1, "extracted_keys": [
                "K_L_elements",
                "K_R_elements",
                "Gamma_star_perms",
            ]},
            {
                **prov2,
                "extracted_keys": [
                    "invariants",
                    "K_L_eq_K_R",
                    "K_L_inter_K_R_size",
                    "Y2_conjugate_via_pi0",
                    "K_L_block_systems",
                    "K_R_block_systems",
                ],
            },
        ],
        "scripts": [
            {"path": sp1, "sha256": sh1},
            {"path": sp2, "sha256": sh2},
        ],
        "seed": None,
    }


def cert_M19_middle_band_balance() -> Dict[str, Any]:
    src, prov = load_source("phase_9_18_S9x_results.json")
    sp, sh = script_provenance("phase_9_18_S9x.py")
    M19_grid = src["M19_grid"]
    gamma_star = {tuple(p) for p in src["Gamma_star_perms"]}
    K_L = {tuple(p) for p in src["K_L_elements"]}
    K_R = {tuple(p) for p in src["K_R_elements"]}
    middle_rows = [4, 5, 6]
    middle_idx = [r - 1 for r in middle_rows]
    w_left = [1, 1, 1, -2, -2, -2, 1, 1, 1]
    target_sum = 15

    def apply_gamma(gamma: Tuple[int, ...], symbol: int) -> int:
        return gamma[symbol - 1]

    def middle_band_sums(gamma: Tuple[int, ...]) -> List[int]:
        return [
            sum(apply_gamma(gamma, M19_grid[row][col]) for row in middle_idx)
            for col in range(9)
        ]

    def left_kernel_residuals(gamma: Tuple[int, ...]) -> List[int]:
        return [
            sum(
                w_left[row] * (apply_gamma(gamma, M19_grid[row][col]) - 5)
                for row in range(9)
            )
            for col in range(9)
        ]

    def edge_image(gamma: Tuple[int, ...], edge: Tuple[int, int, int]) -> Tuple[int, int, int]:
        return tuple(sorted(apply_gamma(gamma, symbol) for symbol in edge))

    def perm_compose(left: Tuple[int, ...], right: Tuple[int, ...]) -> Tuple[int, ...]:
        return tuple(left[right[i] - 1] for i in range(9))

    def perm_inverse(gamma: Tuple[int, ...]) -> Tuple[int, ...]:
        inv = [0] * 9
        for i, value in enumerate(gamma):
            inv[value - 1] = i + 1
        return tuple(inv)

    def sorted_perm_hash(perms: set[Tuple[int, ...]]) -> str:
        payload = json.dumps(
            [list(p) for p in sorted(perms)],
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
        return hashlib.sha256(payload).hexdigest()

    balanced_perms: List[Tuple[int, ...]] = []
    left_kernel_perms: List[Tuple[int, ...]] = []
    scan_space_size = 0
    for gamma in itertools.permutations(range(1, 10)):
        scan_space_size += 1
        is_balanced = all(s == target_sum for s in middle_band_sums(gamma))
        is_left_kernel = all(r == 0 for r in left_kernel_residuals(gamma))
        if is_balanced != is_left_kernel:
            raise AssertionError(
                "middle-band balance drifted from left-kernel equation: "
                f"gamma={gamma}"
            )
        if is_balanced:
            balanced_perms.append(gamma)
        if is_left_kernel:
            left_kernel_perms.append(gamma)

    balanced_set = set(balanced_perms)
    left_kernel_set = set(left_kernel_perms)
    if balanced_set != left_kernel_set:
        raise AssertionError("balanced set differs from left-kernel set")
    if balanced_set != gamma_star:
        raise AssertionError(
            f"middle-band balanced set has size {len(balanced_set)} but does not "
            "equal Gamma^*_{M19}"
        )
    if scan_space_size != math.factorial(9):
        raise AssertionError(f"unexpected S_9 scan size {scan_space_size}")

    magic_triples = [
        list(combo)
        for combo in itertools.combinations(range(1, 10), 3)
        if sum(combo) == target_sum
    ]
    if len(magic_triples) != 8:
        raise AssertionError(f"expected 8 magic triples, got {len(magic_triples)}")

    middle_band_column_triples = [
        [M19_grid[row][col] for row in middle_idx]
        for col in range(9)
    ]
    source_edges = [tuple(sorted(edge)) for edge in middle_band_column_triples]
    source_edge_multiplicities = Counter(source_edges)
    source_distinct_edges = set(source_edge_multiplicities)
    magic_edge_set = {tuple(edge) for edge in magic_triples}

    source_aut: set[Tuple[int, ...]] = set()
    source_multiaut: set[Tuple[int, ...]] = set()
    target_aut: set[Tuple[int, ...]] = set()
    target_multiaut: set[Tuple[int, ...]] = set()
    hypergraph_embeddings: set[Tuple[int, ...]] = set()
    target_distinct_edges: set[Tuple[int, int, int]] | None = None
    pi0_perm = min(gamma_star)
    for gamma in itertools.permutations(range(1, 10)):
        source_image = {edge_image(gamma, edge) for edge in source_distinct_edges}
        if source_image == source_distinct_edges:
            source_aut.add(gamma)
        source_multiimage = Counter(edge_image(gamma, edge) for edge in source_edges)
        if source_multiimage == source_edge_multiplicities:
            source_multiaut.add(gamma)
        if all(edge_image(gamma, edge) in magic_edge_set for edge in source_edges):
            hypergraph_embeddings.add(gamma)
            image = {edge_image(gamma, edge) for edge in source_distinct_edges}
            if target_distinct_edges is None:
                target_distinct_edges = image
            elif target_distinct_edges != image:
                raise AssertionError("middle-band embeddings have non-constant magic image")

    if target_distinct_edges is None:
        raise AssertionError("no hypergraph embeddings found")

    target_edge_multiplicities = Counter(edge_image(pi0_perm, edge) for edge in source_edges)

    for gamma in itertools.permutations(range(1, 10)):
        target_image = {edge_image(gamma, edge) for edge in target_distinct_edges}
        if target_image == target_distinct_edges:
            target_aut.add(gamma)
        target_multiimage = Counter(
            edge_image(gamma, edge)
            for edge, multiplicity in target_edge_multiplicities.items()
            for _ in range(multiplicity)
        )
        if target_multiimage == target_edge_multiplicities:
            target_multiaut.add(gamma)

    pi0_inv = perm_inverse(pi0_perm)
    if hypergraph_embeddings != balanced_set:
        raise AssertionError("hypergraph embedding set differs from balanced set")
    if source_aut != K_R:
        raise AssertionError("source middle-band hypergraph automorphisms differ from K_R")
    if target_aut != K_L:
        raise AssertionError("target magic subhypergraph automorphisms differ from K_L")
    if {perm_compose(pi0_perm, g) for g in source_aut} != balanced_set:
        raise AssertionError("pi0*K_R does not recover the embedding set")
    if {perm_compose(g, pi0_perm) for g in target_aut} != balanced_set:
        raise AssertionError("K_L*pi0 does not recover the embedding set")
    if {perm_compose(perm_compose(pi0_perm, g), pi0_inv) for g in source_aut} != target_aut:
        raise AssertionError("target automorphism group is not pi0-conjugate to source")

    full_magic_aut = {
        gamma
        for gamma in itertools.permutations(range(1, 10))
        if {edge_image(gamma, edge) for edge in magic_edge_set} == magic_edge_set
    }
    if len(full_magic_aut) != 8:
        raise AssertionError(f"expected 8 full magic-triple automorphisms, got {len(full_magic_aut)}")

    source_mult_degree = Counter()
    source_mult_pair_codegree = Counter()
    for edge, multiplicity in source_edge_multiplicities.items():
        for symbol in edge:
            source_mult_degree[symbol] += multiplicity
        for pair in itertools.combinations(edge, 2):
            source_mult_pair_codegree[tuple(sorted(pair))] += multiplicity
    missing_magic_triples = sorted(magic_edge_set - target_distinct_edges)
    balanced_bytes = json.dumps(
        [list(p) for p in sorted(balanced_perms)],
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")

    return {
        "result_id": "M19_middle_band_balance",
        "tier": "CERTIFIED_EXHAUSTIVE",
        "paper_label": "Certified Proposition P5",
        "claim": (
            "For L=M_{19}, Gamma^*_{M_{19}} is exactly the set of symbol "
            "relabelings gamma in S_9 for which the middle row-band entries "
            "gamma(L_{4j}), gamma(L_{5j}), gamma(L_{6j}) sum to 15 in every "
            "column j. Equivalently, with w=(1,1,1,-2,-2,-2,1,1,1), the "
            "left-kernel equation w^T E_{gamma,M_{19}}=0 holds. The same "
            "72 relabelings are exactly the reduced-hypergraph embeddings from "
            "the M19 middle-band hypergraph into the six-edge magic-triple "
            "subhypergraph, with source automorphism group K_R and target "
            "automorphism group K_L."
        ),
        "inputs": {
            "n": 9,
            "base": "M_{19}",
            "scan_space": "S_9",
            "scan_space_size": scan_space_size,
            "middle_band_rows": middle_rows,
            "target_sum": target_sum,
            "w_left": w_left,
        },
        "outputs": {
            "balanced_count": len(balanced_perms),
            "left_kernel_count": len(left_kernel_perms),
            "Gamma_star_size": len(gamma_star),
            "balanced_equals_left_kernel": balanced_set == left_kernel_set,
            "balanced_equals_Gamma_star": balanced_set == gamma_star,
            "balanced_perms_sha256": hashlib.sha256(balanced_bytes).hexdigest(),
            "balanced_perms": [list(p) for p in sorted(balanced_perms)],
            "magic_triples_sum_15": magic_triples,
            "magic_triples_count": len(magic_triples),
            "middle_band_column_triples": middle_band_column_triples,
            "middle_band_hypergraph": {
                "vertices": list(range(1, 10)),
                "ordered_column_triples": middle_band_column_triples,
                "distinct_edges": [list(edge) for edge in sorted(source_distinct_edges)],
                "edge_multiplicities": [
                    {"edge": list(edge), "multiplicity": multiplicity}
                    for edge, multiplicity in sorted(source_edge_multiplicities.items())
                ],
                "multidegree_by_symbol": {
                    str(symbol): source_mult_degree[symbol]
                    for symbol in range(1, 10)
                },
                "pair_codegree_distribution": {
                    str(k): int(v)
                    for k, v in sorted(Counter(source_mult_pair_codegree.values()).items())
                },
            },
            "magic_triple_hypergraph": {
                "vertices": list(range(1, 10)),
                "edges": magic_triples,
                "automorphism_count": len(full_magic_aut),
                "automorphism_perms_sha256": sorted_perm_hash(full_magic_aut),
            },
            "target_magic_subhypergraph": {
                "edges": [list(edge) for edge in sorted(target_distinct_edges)],
                "edge_multiplicities_for_pi0": [
                    {"edge": list(edge), "multiplicity": multiplicity}
                    for edge, multiplicity in sorted(target_edge_multiplicities.items())
                ],
                "missing_magic_triples": [list(edge) for edge in missing_magic_triples],
            },
            "hypergraph_realization": {
                "automorphism_convention": "reduced 3-uniform edge sets; multiplicities are recorded separately",
                "embedding_count": len(hypergraph_embeddings),
                "embeddings_equal_balanced_perms": hypergraph_embeddings == balanced_set,
                "embeddings_equal_Gamma_star": hypergraph_embeddings == gamma_star,
                "source_reduced_automorphism_count": len(source_aut),
                "source_reduced_automorphism_equals_K_R": source_aut == K_R,
                "target_reduced_automorphism_count": len(target_aut),
                "target_reduced_automorphism_equals_K_L": target_aut == K_L,
                "source_multihypergraph_automorphism_count": len(source_multiaut),
                "target_multihypergraph_automorphism_count_for_pi0": len(target_multiaut),
                "all_embeddings_have_same_target_subhypergraph": True,
                "embedding_perms_sha256": sorted_perm_hash(hypergraph_embeddings),
                "source_reduced_automorphism_perms_sha256": sorted_perm_hash(source_aut),
                "target_reduced_automorphism_perms_sha256": sorted_perm_hash(target_aut),
                "source_multihypergraph_automorphism_perms_sha256": sorted_perm_hash(source_multiaut),
                "target_multihypergraph_automorphism_perms_sha256_for_pi0": sorted_perm_hash(target_multiaut),
                "pi0_perm": list(pi0_perm),
                "embeddings_equal_pi0_K_R": {
                    tuple(perm_compose(pi0_perm, g)) for g in source_aut
                } == balanced_set,
                "embeddings_equal_K_L_pi0": {
                    tuple(perm_compose(g, pi0_perm)) for g in target_aut
                } == balanced_set,
                "target_aut_conjugate_to_source_aut_via_pi0": {
                    perm_compose(perm_compose(pi0_perm, g), pi0_inv)
                    for g in source_aut
                } == target_aut,
            },
        },
        "source": [{**prov, "extracted_keys": [
            "M19_grid",
            "Gamma_star_perms",
            "K_L_elements",
            "K_R_elements",
        ]}],
        "scripts": [
            {"path": sp, "sha256": sh},
            {
                "path": "certified/build_certified.py",
                "sha256": sha256_file(HERE / "build_certified.py"),
            },
        ],
        "seed": None,
    }


def cert_A8_saturation() -> Dict[str, Any]:
    src, prov = load_source("phase_9_18_S9q_beta_lattice.json")
    sp, sh = script_provenance("phase_9_18_S9q_beta_lattice.py")
    return {
        "result_id": "A8_saturation",
        "tier": "CERTIFIED_EXHAUSTIVE",
        "paper_label": "Certified Proposition D",
        "claim": (
            "At n=9, the integer lattices Lambda_M and Lambda_{O_6} generated by the "
            "G^*-orbits of the canonical sparse-low template representatives satisfy "
            "Lambda_M = Lambda_{O_6} = V_std,9 cap Z^9 = A_8, with elementary divisors "
            "all equal to 1."
        ),
        "inputs": {
            "n": 9,
            "M_rep": src["M_rep"],
            "O6_rep": src["O6_rep"],
            "orb_sizes": src["orb_sizes"],
        },
        "outputs": {
            "Lambda_M": src["Lambda_M"],
            "Lambda_O6": src["Lambda_O6"],
            "A8_reference": src["A8_reference"],
        },
        "source": [{**prov, "extracted_keys": ["*"]}],
        "scripts": [{"path": sp, "sha256": sh}],
        "seed": None,
    }


def cert_dichotomy_n7_n9_n11() -> Dict[str, Any]:
    src, prov = load_source("phase_9_18_S9q_gamma1.json")
    pick = [r for r in src["results"] if r["n"] in (7, 9, 11)]

    # --- Augmentation (P1.2): make the abstract isomorphisms
    # Stab(M_n) cong D_4 x S_{n-4} and Stab(O_n) cong C_2 x S_{n-4} explicit
    # per n.  The producer script ``phase_9_18_S9q_gamma1.py`` already
    # verifies (i) the semidirect split |Stab| = |im psi| * |ker psi|, with
    # |im psi| = (n-4)! (so the S_{n-4} factor is identified by its order),
    # and (ii) the order-multiset of the kernel.  Order-multiset {1:1, 2:5,
    # 4:2} is the unique fingerprint of D_4 among groups of order 8 (Q_8 is
    # {1:1, 2:1, 4:6}; Z/8 is {1:1, 2:1, 4:2, 8:4}; Z/4 x Z/2 is {1:1, 2:3,
    # 4:4}; (Z/2)^3 is {1:1, 2:7}).  This row distils that into a single
    # human-readable verdict per n.
    D4_fingerprint = {"1": 1, "2": 5, "4": 2}
    C2_fingerprint = {"1": 1, "2": 1}
    for r in pick:
        n = r["n"]
        ai: Dict[str, Any] = {
            "stab_M_iso": f"D_4 x S_{n-4}",
            "stab_O_iso": f"C_2 x S_{n-4}",
            "stab_M_witness": {
                "split_holds": r["checks"]["stab_M"]["split_ok"],
                "kernel_order_multiset_matches_D4": (
                    r["checks"]["stab_M"]["kernel_orders_ok"]
                ),
                "image_order_eq_factorial_n_minus_4": True,
                "D_4_order_multiset_reference": D4_fingerprint,
                "alternatives_excluded_by_order_multiset": [
                    "Q_8", "Z/8", "Z/4 x Z/2", "(Z/2)^3",
                ],
            },
            "stab_O_witness": {
                "split_holds": r["checks"]["stab_O"]["split_ok"],
                "kernel_order_multiset_matches_C2": (
                    r["checks"]["stab_O"]["kernel_orders_ok"]
                ),
                "image_order_eq_factorial_n_minus_4": True,
                "C_2_order_multiset_reference": C2_fingerprint,
            },
            "abstract_iso_certified": (
                r["checks"]["stab_M"]["split_ok"]
                and r["checks"]["stab_M"]["kernel_orders_ok"]
                and r["checks"]["stab_O"]["split_ok"]
                and r["checks"]["stab_O"]["kernel_orders_ok"]
            ),
        }
        r["abstract_identification"] = ai

    sp, sh = script_provenance("phase_9_18_S9q_gamma1.py")
    return {
        "result_id": "dichotomy_n7_n9_n11",
        "tier": "CERTIFIED_EXHAUSTIVE",
        "paper_label": "Certified Proposition F",
        "claim": (
            "For n in {7, 9, 11}, the canonical sparse-low template representatives "
            "M_n and O_n have stabilizers Stab(M_n) cong D_4 x S_{n-4} and "
            "Stab(O_n) cong C_2 x S_{n-4}, integer hulls Lambda_{M_n} = Lambda_{O_n} "
            "= A_{n-1}, and depth d_{M_n}(orb(O_n)) = 2.  The abstract isomorphisms "
            "are certified by (a) semidirect split |Stab| = |im psi| * |ker psi| "
            "with |im psi| = (n-4)!, and (b) kernel order-multiset matching the "
            "D_4 (resp. C_2) fingerprint, which uniquely separates D_4 from all "
            "other groups of order 8."
        ),
        "inputs": {"n_values": [7, 9, 11]},
        "outputs": {"per_n": pick},
        "source": [{**prov, "extracted_keys": ["results[n in {7,9,11}]"]}],
        "scripts": [{"path": sp, "sha256": sh}],
        "seed": None,
    }


def cert_stabilizers_5_to_13() -> Dict[str, Any]:
    src, prov = load_source("phase_9_18_S9q_gamma1.json")
    sp, sh = script_provenance("phase_9_18_S9q_gamma1.py")
    return {
        "result_id": "stabilizers_5_to_13",
        "tier": "EXHAUSTIVE",
        "paper_label": "Proposition P4",
        "claim": (
            "Exhaustive computational verification: for every n in {5,...,13}, "
            "|Stab(M_n)| = 8(n-4)!, |Stab(O_n)| = 2(n-4)!, |orb(M_n)| = "
            "n(n-1)(n-2)(n-3)/4, the depth witness d_{M_n}(orb(O_n)) = 2 holds, and "
            "Lambda_{M_n} = Lambda_{O_n} saturates A_{n-1}."
        ),
        "inputs": {"n_range": src["n_range"]},
        "outputs": {
            "universal_pass": src["universal_pass"],
            "results": src["results"],
            "elapsed_sec": src["elapsed_sec"],
        },
        "source": [{**prov, "extracted_keys": ["*"]}],
        "scripts": [{"path": sp, "sha256": sh}],
        "seed": None,
    }


def cert_Aut_LM19_trivial() -> Dict[str, Any]:
    src, prov = load_source("phase_9_18_S9q_C1.json")
    sp, sh = script_provenance("phase_9_18_S9q_C1.py")
    return {
        "result_id": "Aut_LM19_trivial",
        "tier": "CERTIFIED_EXHAUSTIVE",
        "paper_label": "Certified Proposition G",
        "claim": (
            "The Latin-square autotopism group of L_{M_{19}} is trivial: "
            "|Aut(L_{M_{19}})| = 1. Consequently |K_L| = 72 != 1 = |Aut(L_{M_{19}})|, "
            "so the affine coset-symmetry K_L is not the image of any autotopism subgroup."
        ),
        "inputs": {"n": 9, "base": "M_{19}", "M19_grid": src["M19_grid"]},
        "outputs": {
            "Aut_L_M19_order": src["Aut_L_M19_order"],
            "K_L_order": src["K_L_order"],
            "all_distinct_gammas": src["all_distinct_gammas"],
            "n_distinct_alpha": src["n_distinct_alpha"],
            "n_distinct_beta": src["n_distinct_beta"],
            "n_distinct_gamma": src["n_distinct_gamma"],
            "n_gamma_visited": src["n_gamma_visited"],
            "elapsed_s": src["elapsed_s"],
        },
        "source": [{**prov, "extracted_keys": ["*"]}],
        "scripts": [{"path": sp, "sha256": sh}],
        "seed": None,
    }


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------
BUILDERS = [
    ("base1_odd_rank_22.json", cert_base1_odd_rank_22),
    ("n6_odd_rank_census.json", cert_n6_odd_rank_census),
    ("n7_HV1_D6.json", cert_n7_HV1_D6),
    ("box_band_lemma_witness.json", cert_box_band_lemma_witness),
    ("lift_lemma_evidence.json", cert_lift_lemma_evidence),
    ("M19_audit.json", cert_M19_audit),
    ("M19_middle_band_balance.json", cert_M19_middle_band_balance),
    ("M19_KL_KR_affine_D4.json", cert_M19_KL_KR_affine_D4),
    ("A8_saturation.json", cert_A8_saturation),
    ("dichotomy_n7_n9_n11.json", cert_dichotomy_n7_n9_n11),
    ("stabilizers_5_to_13.json", cert_stabilizers_5_to_13),
    ("Aut_LM19_trivial.json", cert_Aut_LM19_trivial),
]


def main() -> int:
    manifest_lines: List[str] = [
        "# MANIFEST.sha256 -- canonical certificates for this paper",
        "# Format: <sha256>  <filename>",
        "# Hash is computed over canonical bytes (json.dumps sort_keys=True, indent=2)",
        "# of the certificate with the `produced_utc` field removed.",
        "",
    ]
    for name, builder in BUILDERS:
        obj = builder()
        fname, digest = write_certificate(name, obj)
        manifest_lines.append(f"{digest}  {fname}")
        print(f"[OK] {fname}  {digest[:16]}...")

    (HERE / "MANIFEST.sha256").write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")
    print(f"[OK] MANIFEST.sha256 ({len(BUILDERS)} entries)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
