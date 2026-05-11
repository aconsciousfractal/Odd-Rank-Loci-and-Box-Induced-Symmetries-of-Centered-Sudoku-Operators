"""verify_all.py
================
Quick claim-level verifier for this certified package certified package.

What it checks (no source-data re-derivation; for that see
``reproduce_all.py``):

    1. MANIFEST.sha256 parses and lists exactly 12 entries.
  2. Each certificate JSON is valid, contains the required schema fields
     (``result_id``, ``tier``, ``paper_label``, ``claim``, ``inputs``,
     ``outputs``, ``source``, ``scripts``, ``schema_version``,
     ``produced_utc``), and ``schema_version == 1``.
  3. SHA-256 over canonical bytes (with ``produced_utc`` removed) matches
     the digest in MANIFEST.sha256.
    4. Lightweight semantic spot-checks on the 12 ``outputs`` payloads:
       * base1_odd_rank_22: 22 odd perms, all of length 9 over 1..9;
         four_way_certificates count == 22.
       * n6_odd_rank_census: z-stat ~ 4.013706 (tight tol),
         exact two-sided p = erfc(|z|/sqrt(2)) within 1e-9.
       * n7_HV1_D6: |H_V1|=12, dihedral relation verified, subgroup of B_3.
             * box_band_lemma_witness: square Sudoku cases vanish, cyclic LS-9
                 standard boxes fail, rectangular Roku-Doku vanishes, and the
                 non-Cartesian gerechte broken-diagonal witness has zero region sums.
        * lift_lemma_evidence: holds on 22 base1 + 100 M19 + 72 shared left-kernel.
       * M19_audit: |Gamma_rank|=100, |Gamma_star|=72, residue=28,
            shared_left_kernel_w=[1,1,1,-2,-2,-2,1,1,1], |K_L|=|K_R|=72,
         K_L cap K_R = {e}, K_L conjugate to K_R via pi_0,
         left/right coset tests PASS, double_coset_size = 72,
            ambient_stab_S9(w_left) = 4320 = |S_3 x S_6|.
             * M19_middle_band_balance: Gamma^*_{M19} equals the 72 relabelings
                 whose middle row-band column sums are all 15, equivalently
                w^T E=0 for w=(1,1,1,-2,-2,-2,1,1,1); the same 72 are
                the reduced-hypergraph embeddings from the M19 middle-band
                hypergraph to a six-edge magic-triple subhypergraph, with
                Aut(source)=K_R and Aut(target)=K_L.
             * M19_KL_KR_affine_D4: |K_L|=|K_R|=72, both 3^2:D_4
                 affine semidirect-product witnesses, abelianization C_2 x C_2,
                 and not contained in Stab_{S_9}(w_left).
       * A8_saturation: Lambda_M = Lambda_O6 = A_8.
       * dichotomy_n7_n9_n11: n=7,9,11 keys present.
       * stabilizers_5_to_13: keys for n in {5..13}.
       * Aut_LM19_trivial: |Aut(L_M19)| == 1.

Usage::

    python verify_all.py            # full quick check
    python verify_all.py --quick    # alias (default)
    python verify_all.py --quiet    # only emit FAIL lines

Exit code: 0 if every check passes, 1 otherwise.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "MANIFEST.sha256"
REPO_ROOT = HERE.parent  # <repo>/certified -> <repo>/

REQUIRED_FIELDS = (
    "result_id",
    "tier",
    "paper_label",
    "claim",
    "inputs",
    "outputs",
    "source",
    "scripts",
    "schema_version",
    "produced_utc",
)

ALLOWED_TIERS = {
    "THEOREM",
    "EXHAUSTIVE",
    "EMPIRICAL",
    "OBSERVED TEMPLATE",
    "CONJECTURE",
    "CERTIFIED_EXHAUSTIVE",
}


def canonical_bytes(obj: Dict[str, Any]) -> bytes:
    return json.dumps(
        obj, sort_keys=True, indent=2, ensure_ascii=False
    ).encode("utf-8")


def manifest_hash(obj: Dict[str, Any]) -> str:
    obj2 = {k: v for k, v in obj.items() if k != "produced_utc"}
    return hashlib.sha256(canonical_bytes(obj2)).hexdigest()


def parse_manifest() -> List[Tuple[str, str]]:
    out: List[Tuple[str, str]] = []
    for raw in MANIFEST.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) != 2:
            raise ValueError(f"malformed MANIFEST line: {raw!r}")
        out.append((parts[0], parts[1]))
    return out


# ---------------------------------------------------------------------------
# Semantic spot-checks per certificate.
# Each returns (ok: bool, detail: str).
# ---------------------------------------------------------------------------

def chk_base1_odd_rank_22(o: Dict[str, Any]) -> Tuple[bool, str]:
    out = o["outputs"]
    perms = out["odd_perms_full_list"]
    cnt = out["odd_perms_count"]
    if cnt != 22 or len(perms) != 22:
        return False, f"expected 22 perms, got count={cnt}, list={len(perms)}"
    for p in perms:
        if sorted(p) != list(range(1, 10)):
            return False, f"non-permutation entry: {p}"
    if out["odd_rank_count"] != 22:
        return False, f"odd_rank_count={out['odd_rank_count']}"
    if out.get("four_way_n_certs") != 22:
        return False, f"four_way_n_certs={out.get('four_way_n_certs')}"
    fwc = out.get("four_way_certificates")
    if not isinstance(fwc, list) or len(fwc) != 22:
        return False, f"four_way_certificates len={len(fwc) if isinstance(fwc, list) else type(fwc).__name__}"
    # SNF/F_2 distribution sanity (paper §4.3): NOT all torsion-free, NOT
    # uniform F_2 rank. The published paper text references the actual
    # multi-class distribution; if a single SNF class or a single F_2 rank
    # appears for all 22, the certificate has drifted.
    from collections import Counter
    snf_classes = Counter(tuple(c["snf_Z_diagonal"]) for c in fwc)
    f2_ranks = Counter(c["rank_F2"] for c in fwc)
    if len(snf_classes) < 2:
        return False, f"§4.3: SNF over Z is not multi-class as paper claims: {dict(snf_classes)}"
    if len(f2_ranks) < 2:
        return False, f"§4.3: F_2 rank is uniform; paper claims a non-trivial histogram: {dict(f2_ranks)}"
    # Per-cert sanity: all ranks integer 7, kernel sums to 0, SNF length 9.
    for c in fwc:
        if c["exact_rank_Z"] != 7:
            return False, f"perm {c['perm']}: exact_rank_Z={c['exact_rank_Z']}"
        if not c.get("kernel_vector_sum_zero", False):
            return False, f"perm {c['perm']}: kernel not in V_std"
        if len(c["snf_Z_diagonal"]) != 9:
            return False, f"perm {c['perm']}: SNF length != 9"
        if c["rank_Z_from_snf"] != 7:
            return False, f"perm {c['perm']}: SNF rank mismatch"
    return True, (
        f"22 perms rank-7; SNF classes={len(snf_classes)} "
        f"(max={snf_classes.most_common(1)[0][1]}); F_2 ranks={dict(f2_ranks)}"
    )


def chk_n6_odd_rank_census(o: Dict[str, Any]) -> Tuple[bool, str]:
    import math
    out = o["outputs"]
    z = float(out["two_proportion_z"])
    if abs(z - 4.013706) > 1e-4:
        return False, f"z={z} not ~ 4.013706"
    if out["sudoku6_full_rank_F2_count"] != 2592:
        return False, f"sudoku6 count={out['sudoku6_full_rank_F2_count']}"
    if out["ls6_full_rank_F2_count"] != 69120:
        return False, f"ls6 count={out['ls6_full_rank_F2_count']}"
    p_exact = float(out["two_proportion_pvalue_exact"])
    p_recompute = math.erfc(abs(z) / math.sqrt(2.0))
    if abs(p_exact - p_recompute) > 1e-9:
        return False, f"p_exact={p_exact} disagrees with erfc(|z|/sqrt2)={p_recompute}"
    if out.get("two_proportion_pvalue_formula") != "erfc(|z| / sqrt(2))":
        return False, "pvalue_formula label mismatch"
    return True, f"z={z:.6f}, p={p_exact:.3e} = erfc(|z|/sqrt(2))"


def chk_n7_HV1_D6(o: Dict[str, Any]) -> Tuple[bool, str]:
    out = o["outputs"]
    if out.get("order") != 12:
        return False, f"|H_V1|={out.get('order')}"
    if not out.get("dihedral_relation_verified"):
        return False, "dihedral relation not verified"
    if not out.get("generates_full_H_V1"):
        return False, "generators do not generate full H_V1"
    if not out.get("is_subgroup_of_B3"):
        return False, "H_V1 not subgroup of B_3"
    rd = out.get("L_V1_rank_distribution_in_B3", {})
    if rd.get("5") != 12:
        return False, f"rank-5 in B_3 = {rd.get('5')}"
    return True, "|H_V1|=12, D_6, subgroup of B_3, rank-5 count=12"


def chk_box_band(o: Dict[str, Any]) -> Tuple[bool, str]:
    out = o["outputs"]
    if not out["sudoku_base1"]["is_zero"]:
        return False, "base1 M_rb_cs is not zero"
    if not out["sudoku_M19"]["is_zero"]:
        return False, "M19 M_rb_cs is not zero"
    cnt = out["non_sudoku_LS_cyclic_counterexample"]
    if cnt["is_zero"]:
        return False, "non-Sudoku counterexample is zero (should be nonzero)"
    if cnt["frobenius_norm_sq"] <= 0:
        return False, f"non-Sudoku Frob^2={cnt['frobenius_norm_sq']}"
    rect = out.get("rectangular_roku_doku_2x3")
    if not rect or not rect.get("is_zero"):
        return False, f"rectangular Roku-Doku check failed: {rect}"
    if rect.get("box_height") != 2 or rect.get("box_width") != 3:
        return False, f"rectangular dimensions drifted: {rect}"
    ger = out.get("gerechte_cyclic_LS9_broken_diagonals")
    if not ger:
        return False, "missing gerechte cyclic LS9 witness"
    if not ger.get("regions_partition_cells") or not ger.get("each_region_symbols_once"):
        return False, f"invalid gerechte regions: {ger}"
    if not ger.get("all_region_sums_zero"):
        return False, f"gerechte region sums nonzero: {ger.get('region_sums_twice_centered')}"
    if ger.get("all_regions_cartesian_products") or ger.get("any_region_cartesian_product"):
        return False, "gerechte witness should be non-Cartesian"
    return True, (
        f"base1=0, M19=0, rectangular 2x3=0, gerechte diagonal=0, "
        f"standard cyclic boxes Frob^2={cnt['frobenius_norm_sq']}"
    )


def chk_lift_lemma(o: Dict[str, Any]) -> Tuple[bool, str]:
    out = o["outputs"]
    b1 = out["base1"]
    m19 = out["M19"]
    shared = out.get("M19_shared_left_kernel_72") or out.get("M19_shared_kernel_72")
    if not (b1["lift_lemma_holds_on_all_checked"] and b1["n_perms"] == 22):
        return False, f"base1 fail: {b1}"
    if not (m19["lift_lemma_holds_on_all_checked"] and m19["n_perms"] == 100):
        return False, f"M19 fail (need 100): {m19}"
    if shared is None or not (shared["lift_lemma_holds_on_all_checked"] and shared["n_perms"] == 72):
        return False, f"shared left-kernel subset fail: {shared}"
    if b1["n_full_lift_check"] != 22:
        return False, f"base1 not full-checked: {b1['n_full_lift_check']}"
    if m19["n_full_lift_check"] != 100:
        return False, f"M19 not full-checked: {m19['n_full_lift_check']}"
    if shared["n_full_lift_check"] != 72:
        return False, f"shared left-kernel subset not full-checked: {shared['n_full_lift_check']}"
    return True, "Lift Lemma holds on 22 base1 + 100 M19 + 72 shared left-kernel"


def chk_M19_audit(o: Dict[str, Any]) -> Tuple[bool, str]:
    out = o["outputs"]
    if out.get("Gamma_rank_size") != 100:
        return False, f"|Gamma_rank|={out.get('Gamma_rank_size')}"
    if out.get("Gamma_star_size") != 72:
        return False, f"|Gamma_star|={out.get('Gamma_star_size')}"
    if out.get("residue_size") != 28:
        return False, f"residue={out.get('residue_size')}"
    if out.get("K_L_size") != 72 or out.get("K_R_size") != 72:
        return False, f"|K_L|={out.get('K_L_size')}, |K_R|={out.get('K_R_size')}"
    if out.get("K_L_inter_K_R_size") != 1:
        return False, f"K_L cap K_R={out.get('K_L_inter_K_R_size')}"
    if not out.get("K_L_conjugate_to_K_R_via_pi0"):
        return False, "K_L not conjugate to K_R via pi_0"
    if not out.get("left_coset_test_PASS") or not out.get("right_coset_test_PASS"):
        return False, "coset tests not PASS/PASS"
    if out.get("double_coset_size_K_L_pi0_K_R") != 72:
        return False, f"|K_L pi_0 K_R|={out.get('double_coset_size_K_L_pi0_K_R')}"
    ambient_stab = out.get("ambient_coordinate_stabilizer_S_9_of_shared_left_kernel_size")
    if ambient_stab != 4320:
        return False, f"ambient_stab={ambient_stab}"
    if out.get("ambient_stab_iso") != "S_3 x S_6":
        return False, f"ambient_stab_iso={out.get('ambient_stab_iso')}"
    w_left = out.get("shared_left_kernel_w")
    if w_left != [1, 1, 1, -2, -2, -2, 1, 1, 1]:
        return False, f"shared_left_kernel_w={w_left} (expected [1,1,1,-2,-2,-2,1,1,1])"
    legacy = out.get("legacy_value_template_from_source")
    if legacy is not None and legacy != [-2, -2, -2, 1, 1, 1, 1, 1, 1]:
        return False, f"legacy value template drifted: {legacy}"
    return True, (
        "|Gamma_rank|=100, |Gamma_star|=72, residue=28, |K_L|=|K_R|=72, "
        "K_L cap K_R={e}, K_L=pi_0 K_R pi_0^{-1}, shared-left ambient_stab=4320=|S_3xS_6|"
    )


def chk_M19_middle_band_balance(o: Dict[str, Any]) -> Tuple[bool, str]:
    inp = o["inputs"]
    out = o["outputs"]
    if inp.get("scan_space") != "S_9" or inp.get("scan_space_size") != 362880:
        return False, f"scan space mismatch: {inp.get('scan_space')} {inp.get('scan_space_size')}"
    if inp.get("middle_band_rows") != [4, 5, 6]:
        return False, f"middle_band_rows={inp.get('middle_band_rows')}"
    if inp.get("target_sum") != 15:
        return False, f"target_sum={inp.get('target_sum')}"
    if inp.get("w_left") != [1, 1, 1, -2, -2, -2, 1, 1, 1]:
        return False, f"w_left={inp.get('w_left')}"
    if out.get("balanced_count") != 72 or out.get("left_kernel_count") != 72:
        return False, f"balanced={out.get('balanced_count')}, left_kernel={out.get('left_kernel_count')}"
    if out.get("Gamma_star_size") != 72:
        return False, f"Gamma_star_size={out.get('Gamma_star_size')}"
    if not out.get("balanced_equals_left_kernel"):
        return False, "balanced set differs from left-kernel set"
    if not out.get("balanced_equals_Gamma_star"):
        return False, "balanced set differs from Gamma^*_{M19}"
    perms = out.get("balanced_perms", [])
    if len(perms) != 72:
        return False, f"balanced_perms length={len(perms)}"
    for p in perms:
        if not isinstance(p, list) or sorted(p) != list(range(1, 10)):
            return False, f"bad balanced permutation: {p}"
    triples = out.get("magic_triples_sum_15", [])
    if out.get("magic_triples_count") != 8 or len(triples) != 8:
        return False, f"magic triple count={out.get('magic_triples_count')} len={len(triples)}"
    for tri in triples:
        if len(tri) != 3 or len(set(tri)) != 3 or sum(tri) != 15:
            return False, f"bad magic triple: {tri}"
    hsrc = out.get("middle_band_hypergraph", {})
    if len(hsrc.get("distinct_edges", [])) != 6:
        return False, f"middle-band distinct edges={hsrc.get('distinct_edges')}"
    multiplicities = sorted(item.get("multiplicity") for item in hsrc.get("edge_multiplicities", []))
    if multiplicities != [1, 1, 1, 2, 2, 2]:
        return False, f"middle-band edge multiplicities={multiplicities}"
    if set(hsrc.get("multidegree_by_symbol", {}).values()) != {3}:
        return False, f"middle-band multidegrees={hsrc.get('multidegree_by_symbol')}"
    target = out.get("target_magic_subhypergraph", {})
    if len(target.get("edges", [])) != 6 or len(target.get("missing_magic_triples", [])) != 2:
        return False, f"target magic subhypergraph malformed: {target}"
    for tri in target.get("edges", []) + target.get("missing_magic_triples", []):
        if len(tri) != 3 or sum(tri) != 15:
            return False, f"bad target/missing magic triple: {tri}"
    magic_h = out.get("magic_triple_hypergraph", {})
    if magic_h.get("automorphism_count") != 8:
        return False, f"full magic-triple aut count={magic_h.get('automorphism_count')}"
    hg = out.get("hypergraph_realization", {})
    expected_flags = [
        "embeddings_equal_balanced_perms",
        "embeddings_equal_Gamma_star",
        "source_reduced_automorphism_equals_K_R",
        "target_reduced_automorphism_equals_K_L",
        "all_embeddings_have_same_target_subhypergraph",
        "embeddings_equal_pi0_K_R",
        "embeddings_equal_K_L_pi0",
        "target_aut_conjugate_to_source_aut_via_pi0",
    ]
    for flag in expected_flags:
        if not hg.get(flag):
            return False, f"hypergraph realization flag failed: {flag}={hg.get(flag)}"
    if hg.get("embedding_count") != 72:
        return False, f"hypergraph embedding count={hg.get('embedding_count')}"
    if hg.get("source_reduced_automorphism_count") != 72 or hg.get("target_reduced_automorphism_count") != 72:
        return False, (
            f"reduced hypergraph aut counts source={hg.get('source_reduced_automorphism_count')} "
            f"target={hg.get('target_reduced_automorphism_count')}"
        )
    if hg.get("source_multihypergraph_automorphism_count") != 36:
        return False, f"source multihypergraph aut count={hg.get('source_multihypergraph_automorphism_count')}"
    if hg.get("target_multihypergraph_automorphism_count_for_pi0") != 36:
        return False, f"target multihypergraph aut count={hg.get('target_multihypergraph_automorphism_count_for_pi0')}"
    if hg.get("pi0_perm") != [1, 5, 3, 9, 6, 8, 7, 4, 2]:
        return False, f"hypergraph pi0={hg.get('pi0_perm')}"
    return True, (
        "Gamma^*_{M19}=72 middle-band balanced relabelings; reduced hypergraph "
        "embeddings source->magic-subgraph equal Gamma^*, Aut(source)=K_R, Aut(target)=K_L"
    )


def chk_M19_KL_KR(o: Dict[str, Any]) -> Tuple[bool, str]:
    out = o["outputs"]
    KL, KR = out["K_L"], out["K_R"]
    if KL["order"] != 72 or KR["order"] != 72:
        return False, f"|K_L|={KL['order']}, |K_R|={KR['order']}"
    if KL["exponent"] != 12 or KR["exponent"] != 12:
        return False, f"exponent K_L={KL['exponent']}, K_R={KR['exponent']}"
    for side, grp in (("K_L", KL), ("K_R", KR)):
        expected_exp = 1
        for order, count in grp.get("order_distribution", {}).items():
            if int(count) > 0:
                expected_exp = math.lcm(expected_exp, int(order))
        if grp["exponent"] != expected_exp:
            return False, f"{side} exponent inconsistent with order distribution"
    if KL["derived_series_orders"] != [72, 18, 9, 1]:
        return False, f"derived series K_L={KL['derived_series_orders']}"
    if KL["abelianization_order"] != 4 or KL["center_order"] != 1:
        return False, "K_L abelianization/centre mismatch"
    if KL["n_conjugacy_classes"] != 9:
        return False, f"|cc(K_L)|={KL['n_conjugacy_classes']}"
    if out.get("abstract_group_name") != "3^2 : D_4":
        return False, f"abstract_group_name={out.get('abstract_group_name')}"
    if out.get("abstract_group_id") is not None:
        return False, "SmallGroup identifier must not be asserted without a certified GAP witness"
    d4_fingerprint = {"1": 1, "2": 5, "4": 2}
    c3sq_fingerprint = {"1": 1, "3": 8}
    c2sq_fingerprint = {"1": 1, "2": 3}
    for side in ("K_L", "K_R"):
        aff = out.get(f"{side}_affine_structure")
        if not isinstance(aff, dict):
            return False, f"missing {side}_affine_structure"
        if aff.get("abstract_group_name") != "3^2 : D_4":
            return False, f"{side} affine name={aff.get('abstract_group_name')}"
        trans = aff.get("translation_subgroup", {})
        if trans.get("abstract_group_name") != "C_3 x C_3" or trans.get("order") != 9:
            return False, f"{side} translation subgroup mismatch: {trans}"
        if trans.get("order_distribution") != c3sq_fingerprint:
            return False, f"{side} translation order distribution={trans.get('order_distribution')}"
        if not (trans.get("is_closed") and trans.get("is_abelian") and trans.get("is_normal")):
            return False, f"{side} translation subgroup lacks closed/abelian/normal witness"
        if not trans.get("is_regular_on_points"):
            return False, f"{side} translation subgroup is not regular on 9 points"
        point = aff.get("point_stabilizer_1", {})
        if point.get("abstract_group_name") != "D_4" or point.get("order") != 8:
            return False, f"{side} point stabilizer mismatch: {point}"
        if point.get("order_distribution") != d4_fingerprint:
            return False, f"{side} point stabilizer order distribution={point.get('order_distribution')}"
        quotient = aff.get("quotient_by_translation_subgroup", {})
        if quotient.get("abstract_group_name") != "D_4" or quotient.get("order") != 8:
            return False, f"{side} quotient mismatch: {quotient}"
        if quotient.get("order_distribution") != d4_fingerprint:
            return False, f"{side} quotient order distribution={quotient.get('order_distribution')}"
        ab = aff.get("abelianization", {})
        if ab.get("abstract_group_name") != "C_2 x C_2" or ab.get("order") != 4:
            return False, f"{side} abelianization mismatch: {ab}"
        if ab.get("order_distribution") != c2sq_fingerprint:
            return False, f"{side} abelianization order distribution={ab.get('order_distribution')}"
        if not aff.get("semidirect_product_certified"):
            return False, f"{side} semidirect product witness failed"
        if aff.get("is_subgroup_of_ambient_coordinate_stabilizer"):
            return False, f"{side} incorrectly certified inside Stab_S9(w_left)"
        block_count = aff.get("ambient_shared_left_kernel_block_preservation_count")
        if block_count in (None, 72):
            return False, f"{side} ambient block-preservation guard failed"
    if out.get("K_L_inter_K_R_size") != 1:
        return False, f"K_L cap K_R={out.get('K_L_inter_K_R_size')}"
    # pi_0 must be present, must be a permutation of {1..9}, must NOT lie
    # in K_L or K_R, and the conjugacy flag Y2_conjugate_via_pi0 must hold.
    pi0 = out.get("pi0_perm")
    if not isinstance(pi0, list) or sorted(pi0) != list(range(1, 10)):
        return False, f"pi0_perm absent or not a permutation: {pi0}"
    expected_pi0 = [1, 5, 3, 9, 6, 8, 7, 4, 2]
    if pi0 != expected_pi0:
        return False, f"pi0_perm={pi0}, expected {expected_pi0}"
    if out.get("pi0_in_K_L", True):
        return False, "pi_0 must be EXTERNAL to K_L (paper Remark 8.6)"
    if out.get("pi0_in_K_R", True):
        return False, "pi_0 must be EXTERNAL to K_R (paper Remark 8.6)"
    if not out.get("Y2_conjugate_via_pi0", False):
        return False, "pi_0 conjugacy K_L = pi_0 K_R pi_0^{-1} not certified"
    if out.get("K_L_eq_K_R", True):
        return False, "K_L_eq_K_R must be False (paper Certified Proposition 8.5)"
    gens = out.get("K_L_generators_perms", [])
    if pi0 in gens:
        return False, "pi_0 must NOT be one of the K_L generators"
    return True, (
        "K_L,K_R: order 72, exp 12, deriv 72->18->9->1, affine 3^2:D4, "
        "ab C2xC2, not inside Stab_S9(w_left); "
        f"pi_0={pi0} external; Y2 conj via pi_0"
    )

def chk_A8_saturation(o: Dict[str, Any]) -> Tuple[bool, str]:
    out = o["outputs"]
    for key in ("Lambda_M", "Lambda_O6"):
        L = out[key]
        if L["rank"] != 8:
            return False, f"{key} rank={L['rank']}"
        if L["elementary_divisors"] != [1] * 8:
            return False, f"{key} divisors={L['elementary_divisors']}"
        if not L["saturates_V_std"]:
            return False, f"{key} does not saturate V_std"
    return True, "Lambda_M = Lambda_O6 = A_8 (rank 8, all divisors 1)"


def chk_dichotomy(o: Dict[str, Any]) -> Tuple[bool, str]:
    out = o["outputs"]
    seen = {row["n"]: row for row in out["per_n"]}
    for n in (7, 9, 11):
        if n not in seen:
            return False, f"missing n={n}"
        if not seen[n].get("all_pass"):
            return False, f"n={n} all_pass=False"
        ai = seen[n].get("abstract_identification")
        if not ai or not ai.get("abstract_iso_certified"):
            return False, f"n={n} abstract iso not certified"
        if ai.get("stab_M_iso") != f"D_4 x S_{n-4}":
            return False, f"n={n} stab_M_iso label wrong"
        if ai.get("stab_O_iso") != f"C_2 x S_{n-4}":
            return False, f"n={n} stab_O_iso label wrong"
    return True, "n=7,9,11 dichotomy + abstract iso (D_4/C_2 x S_{n-4}) certified"


def chk_stabilizers(o: Dict[str, Any]) -> Tuple[bool, str]:
    out = o["outputs"]
    seen = {row["n"]: row for row in out["results"]}
    for n in range(5, 14):
        if n not in seen:
            return False, f"missing n={n}"
        if not seen[n].get("all_pass"):
            return False, f"n={n} all_pass=False"
    return True, "n=5..13 all pass dichotomy bridges"


def chk_Aut_M19(o: Dict[str, Any]) -> Tuple[bool, str]:
    out = o["outputs"]
    if out.get("Aut_L_M19_order") != 1:
        return False, f"|Aut|={out.get('Aut_L_M19_order')}"
    if out.get("K_L_order") != 72:
        return False, f"|K_L|={out.get('K_L_order')}"
    if out.get("n_gamma_visited") != 362880:
        return False, f"n_gamma_visited={out.get('n_gamma_visited')}"
    return True, "|Aut(L_M19)|=1, |K_L|=72, exhaustive (9!)"


SEMANTIC: Dict[str, Any] = {
    "base1_odd_rank_22.json":     chk_base1_odd_rank_22,
    "n6_odd_rank_census.json":    chk_n6_odd_rank_census,
    "n7_HV1_D6.json":             chk_n7_HV1_D6,
    "box_band_lemma_witness.json": chk_box_band,
    "lift_lemma_evidence.json":   chk_lift_lemma,
    "M19_audit.json":             chk_M19_audit,
    "M19_middle_band_balance.json": chk_M19_middle_band_balance,
    "M19_KL_KR_affine_D4.json":    chk_M19_KL_KR,
    "A8_saturation.json":         chk_A8_saturation,
    "dichotomy_n7_n9_n11.json":   chk_dichotomy,
    "stabilizers_5_to_13.json":   chk_stabilizers,
    "Aut_LM19_trivial.json":      chk_Aut_M19,
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="default mode (claim-level only); accepted for symmetry "
                         "with reproduce_all.py --full")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    failures = 0
    log = (lambda *a, **kw: None) if args.quiet else print

    entries = parse_manifest()
    if len(entries) != 12:
        print(f"FAIL  MANIFEST has {len(entries)} entries, expected 12")
        return 1

    seen = set()
    for digest, name in entries:
        path = HERE / name
        if not path.exists():
            print(f"FAIL  {name}: file missing")
            failures += 1
            continue
        try:
            obj = json.loads(path.read_text(encoding="utf-8"))
        except Exception as e:
            print(f"FAIL  {name}: JSON parse error: {e}")
            failures += 1
            continue
        seen.add(name)

        # 1. Schema fields.
        missing = [f for f in REQUIRED_FIELDS if f not in obj]
        if missing:
            print(f"FAIL  {name}: missing fields {missing}")
            failures += 1
            continue
        if obj["schema_version"] != 1:
            print(f"FAIL  {name}: schema_version={obj['schema_version']}")
            failures += 1
            continue
        if obj["tier"] not in ALLOWED_TIERS:
            print(f"FAIL  {name}: tier={obj['tier']!r} not in {ALLOWED_TIERS}")
            failures += 1
            continue

        # 2. Hash.
        got = manifest_hash(obj)
        if got != digest:
            print(f"FAIL  {name}: sha256 mismatch")
            print(f"        expected {digest}")
            print(f"        got      {got}")
            failures += 1
            continue

        # 3. Semantic spot-check.
        check = SEMANTIC.get(name)
        if check is None:
            log(f"OK    {name}  [no semantic check registered]")
            continue
        try:
            ok, detail = check(obj)
        except Exception as e:
            print(f"FAIL  {name}: semantic check raised: {e}")
            failures += 1
            continue
        if not ok:
            print(f"FAIL  {name}: {detail}")
            failures += 1
            continue
        log(f"OK    {name}  [{obj['tier']}] {detail}")

    # ---- Standalone check: every script/source referenced must live inside REPO_ROOT.
    log("\n[verify_all] standalone-completeness check")
    referenced_scripts: Dict[str, str] = {}
    referenced_sources: Dict[str, str] = {}
    for digest, name in entries:
        path = HERE / name
        if not path.exists():
            continue
        obj = json.loads(path.read_text(encoding="utf-8"))
        for s in obj.get("scripts", []):
            if isinstance(s, dict) and "path" in s:
                referenced_scripts.setdefault(s["path"], s.get("sha256", ""))
        for s in obj.get("source", []):
            if isinstance(s, dict) and "file" in s:
                referenced_sources.setdefault(s["file"], s.get("sha256", ""))
    standalone_failures = 0
    for rel, expected_sha in sorted(referenced_scripts.items()):
        target = (REPO_ROOT / rel)
        if not target.exists():
            print(f"FAIL  standalone: script missing inside repo: {rel}")
            standalone_failures += 1
            continue
        if expected_sha and expected_sha not in ("MISSING", "self"):
            h = hashlib.sha256()
            with target.open("rb") as f:
                for chunk in iter(lambda: f.read(1 << 16), b""):
                    h.update(chunk)
            got = h.hexdigest()
            if got != expected_sha:
                print(f"FAIL  standalone: script sha256 mismatch: {rel}")
                print(f"        expected {expected_sha}")
                print(f"        got      {got}")
                standalone_failures += 1
                continue
        log(f"OK    standalone-script {rel}")
    for rel, expected_sha in sorted(referenced_sources.items()):
        target = (REPO_ROOT / rel)
        if not target.exists():
            print(f"FAIL  standalone: source data missing inside repo: {rel}")
            standalone_failures += 1
            continue
        if expected_sha and expected_sha != "MISSING":
            h = hashlib.sha256()
            with target.open("rb") as f:
                for chunk in iter(lambda: f.read(1 << 16), b""):
                    h.update(chunk)
            got = h.hexdigest()
            if got != expected_sha:
                print(f"FAIL  standalone: source sha256 mismatch: {rel}")
                print(f"        expected {expected_sha}")
                print(f"        got      {got}")
                standalone_failures += 1
                continue
        log(f"OK    standalone-source {rel}")
    failures += standalone_failures

    # Paper-grep guard: forbid drift back to the falsified blanket claim.
    paper_tex = REPO_ROOT / "paper" / "main.tex"
    if paper_tex.exists():
        forbidden = [
            ("cokernel is torsion-free of rank $1$",
             "§4.3 falsified blanket SNF claim"),
            ("agreeing with the integer rank on every",
             "§4.3 falsified blanket F_2 claim"),
            ("$\\pi_0$ is one of the eight generators",
             "Remark 8.6: pi_0 wrongly placed inside K_L"),
            ("odd-rank phenomenon at\n$n = 6$",
             "§5 wording: should be F_2-rank parity proxy"),
        ]
        try:
            text = paper_tex.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            text = paper_tex.read_text(encoding="latin-1")
        paper_failures = 0
        for needle, why in forbidden:
            if needle in text:
                print(f"FAIL  paper-grep: forbidden phrase present — {why}")
                print(f"        needle: {needle!r}")
                paper_failures += 1
        if paper_failures == 0:
            log(f"OK    paper-grep paper/main.tex (no forbidden phrases)")
        failures += paper_failures

    if failures == 0:
        print(f"\n[verify_all] ALL CHECKS PASSED ({len(entries)} certificates, "
              f"{len(referenced_scripts)} scripts, {len(referenced_sources)} sources)")
        return 0
    print(f"\n[verify_all] {failures} FAILURE(S)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
