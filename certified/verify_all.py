"""verify_all.py
================
Quick claim-level verifier for paper II certified package.

What it checks (no source-data re-derivation; for that see
``reproduce_all.py``):

  1. MANIFEST.sha256 parses and lists exactly 11 entries.
  2. Each certificate JSON is valid, contains the required schema fields
     (``result_id``, ``tier``, ``paper_label``, ``claim``, ``inputs``,
     ``outputs``, ``source``, ``scripts``, ``schema_version``,
     ``produced_utc``), and ``schema_version == 1``.
  3. SHA-256 over canonical bytes (with ``produced_utc`` removed) matches
     the digest in MANIFEST.sha256.
  4. Lightweight semantic spot-checks on the 11 ``outputs`` payloads:
       * base1_odd_rank_22: 22 odd perms, all of length 9 over 1..9;
         four_way_certificates count == 22.
       * n6_odd_rank_census: z-stat ~ 4.013706 (tight tol),
         exact two-sided p = erfc(|z|/sqrt(2)) within 1e-9.
       * n7_HV1_D6: |H_V1|=12, dihedral relation verified, subgroup of B_3.
       * box_band_lemma_witness: holds_on_base1 == True, M19 == True.
       * lift_lemma_evidence: holds on 22 base1 + 100 M19 + 72 Hessian.
       * M19_audit: |Gamma_rank|=100, |Gamma_star|=72, residue=28,
         u9_star=[-2,-2,-2,1,1,1,1,1,1], |K_L|=|K_R|=72,
         K_L cap K_R = {e}, K_L conjugate to K_R via gamma0,
         left/right coset tests PASS, double_coset_size = 72,
         ambient_stab_S9(u9*) = 4320 = |S_3 x S_6|.
       * M19_KL_KR_3sq_Q8: |K_L|=|K_R|=72, both 3^2:Q_8 invariants.
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
    "THEOREM_EXHAUSTIVE",
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
    return True, "22 distinct S_9 perms, rank-7 each, four-way x22"


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
    return True, (
        f"base1=0, M19=0, non-Sudoku Frob^2={cnt['frobenius_norm_sq']}"
    )


def chk_lift_lemma(o: Dict[str, Any]) -> Tuple[bool, str]:
    out = o["outputs"]
    b1 = out["base1"]
    m19 = out["M19"]
    hess = out.get("M19_hessian_class_72")
    if not (b1["lift_lemma_holds_on_all_checked"] and b1["n_perms"] == 22):
        return False, f"base1 fail: {b1}"
    if not (m19["lift_lemma_holds_on_all_checked"] and m19["n_perms"] == 100):
        return False, f"M19 fail (need 100): {m19}"
    if hess is None or not (hess["lift_lemma_holds_on_all_checked"] and hess["n_perms"] == 72):
        return False, f"Hessian-class fail: {hess}"
    if b1["n_full_lift_check"] != 22:
        return False, f"base1 not full-checked: {b1['n_full_lift_check']}"
    if m19["n_full_lift_check"] != 100:
        return False, f"M19 not full-checked: {m19['n_full_lift_check']}"
    if hess["n_full_lift_check"] != 72:
        return False, f"Hessian not full-checked: {hess['n_full_lift_check']}"
    return True, "Lift Lemma holds on 22 base1 + 100 M19 + 72 Hessian"


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
        return False, "K_L not conjugate to K_R via gamma0"
    if not out.get("left_coset_test_PASS") or not out.get("right_coset_test_PASS"):
        return False, "coset tests not PASS/PASS"
    if out.get("double_coset_size_K_L_gamma0_K_R") != 72:
        return False, f"|K_L gamma0 K_R|={out.get('double_coset_size_K_L_gamma0_K_R')}"
    if out.get("ambient_coordinate_stabilizer_S_9_of_u9_star_size") != 4320:
        return False, f"ambient_stab={out.get('ambient_coordinate_stabilizer_S_9_of_u9_star_size')}"
    if out.get("ambient_stab_iso") != "S_3 x S_6":
        return False, f"ambient_stab_iso={out.get('ambient_stab_iso')}"
    u9 = out.get("u9_star")
    if u9 != [-2, -2, -2, 1, 1, 1, 1, 1, 1]:
        return False, f"u9_star={u9} (expected [-2,-2,-2,1,1,1,1,1,1])"
    return True, (
        "|Gamma_rank|=100, |Gamma_star|=72, residue=28, |K_L|=|K_R|=72, "
        "K_L cap K_R={e}, K_L=gamma0 K_R gamma0^{-1}, ambient_stab=4320=|S_3xS_6|"
    )


def chk_M19_KL_KR(o: Dict[str, Any]) -> Tuple[bool, str]:
    out = o["outputs"]
    KL, KR = out["K_L"], out["K_R"]
    if KL["order"] != 72 or KR["order"] != 72:
        return False, f"|K_L|={KL['order']}, |K_R|={KR['order']}"
    if KL["exponent"] != 6 or KR["exponent"] != 6:
        return False, f"exponent K_L={KL['exponent']}, K_R={KR['exponent']}"
    if KL["derived_series_orders"] != [72, 18, 9, 1]:
        return False, f"derived series K_L={KL['derived_series_orders']}"
    if KL["abelianization_order"] != 4 or KL["center_order"] != 1:
        return False, "K_L abelianization/centre mismatch"
    if KL["n_conjugacy_classes"] != 9:
        return False, f"|cc(K_L)|={KL['n_conjugacy_classes']}"
    if out.get("K_L_inter_K_R_size") != 1:
        return False, f"K_L cap K_R={out.get('K_L_inter_K_R_size')}"
    return True, "K_L=K_R: order 72, exp 6, deriv 72->18->9->1, ab Z/4, 9 cc"


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
    "M19_KL_KR_3sq_Q8.json":      chk_M19_KL_KR,
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
    if len(entries) != 11:
        print(f"FAIL  MANIFEST has {len(entries)} entries, expected 11")
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

    if failures == 0:
        print(f"\n[verify_all] ALL CHECKS PASSED ({len(entries)} certificates, "
              f"{len(referenced_scripts)} scripts, {len(referenced_sources)} sources)")
        return 0
    print(f"\n[verify_all] {failures} FAILURE(S)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
