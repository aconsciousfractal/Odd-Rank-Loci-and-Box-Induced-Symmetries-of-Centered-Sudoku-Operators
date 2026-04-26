"""reproduce_all.py
====================
Full reproducibility verifier for paper II certified package.

This script re-runs the certificate builders (``build_certified.py``) and
confirms that every regenerated certificate matches the on-disk version
**byte for byte** modulo the ``produced_utc`` timestamp.

It does NOT re-execute the upstream producer scripts under ``scripts/``
(those already wrote the source JSONs under ``data/``); it does re-import
``build_certified`` and call each builder, exercising every assertion
inside (full SymPy lift-lemma check on the 22 base1 + 100 M_{19} + 72
Hessian odd-rank relabelings, Box-Band identity on base1 and M_{19},
kernel-via-SNF on the cyclic LS-9 counterexample, etc.).

Usage::

    python reproduce_all.py            # full check
    python reproduce_all.py --full     # alias

Exit code: 0 if every certificate reproduces byte-equally, 1 otherwise.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, Tuple

HERE = Path(__file__).resolve().parent

# Map result_id -> (certificate filename, builder function name in build_certified).
BUILDERS: Dict[str, Tuple[str, str]] = {
    "base1_odd_rank_22":   ("base1_odd_rank_22.json",      "cert_base1_odd_rank_22"),
    "n6_odd_rank_census":  ("n6_odd_rank_census.json",     "cert_n6_odd_rank_census"),
    "n7_HV1_D6":           ("n7_HV1_D6.json",              "cert_n7_HV1_D6"),
    "box_band_lemma":      ("box_band_lemma_witness.json", "cert_box_band_lemma_witness"),
    "lift_lemma":          ("lift_lemma_evidence.json",    "cert_lift_lemma_evidence"),
    "M19_audit":           ("M19_audit.json",              "cert_M19_audit"),
    "M19_KL_KR_3sq_Q8":    ("M19_KL_KR_3sq_Q8.json",       "cert_M19_KL_KR_3sq_Q8"),
    "A8_saturation":       ("A8_saturation.json",          "cert_A8_saturation"),
    "dichotomy_n7_n9_n11": ("dichotomy_n7_n9_n11.json",    "cert_dichotomy_n7_n9_n11"),
    "stabilizers_5_to_13": ("stabilizers_5_to_13.json",    "cert_stabilizers_5_to_13"),
    "Aut_LM19_trivial":    ("Aut_LM19_trivial.json",       "cert_Aut_LM19_trivial"),
}


def canonical_bytes(obj: Dict[str, Any]) -> bytes:
    return json.dumps(
        obj, sort_keys=True, indent=2, ensure_ascii=False
    ).encode("utf-8")


def strip_utc(obj: Dict[str, Any]) -> Dict[str, Any]:
    return {k: v for k, v in obj.items() if k != "produced_utc"}


def manifest_hash(obj: Dict[str, Any]) -> str:
    return hashlib.sha256(canonical_bytes(strip_utc(obj))).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true",
                    help="full reproducibility check (default)")
    args = ap.parse_args()
    _ = args  # currently the only mode

    # Add the certified folder to sys.path so we can import build_certified
    # without invoking its __main__.
    sys.path.insert(0, str(HERE))
    bc = importlib.import_module("build_certified")

    failures = 0
    print(f"[reproduce_all] re-deriving {len(BUILDERS)} certificates...")
    t0 = time.time()
    for result_id, (fname, builder_name) in BUILDERS.items():
        path = HERE / fname
        if not path.exists():
            print(f"FAIL  {fname}: missing on disk")
            failures += 1
            continue
        on_disk = json.loads(path.read_text(encoding="utf-8"))

        # Re-derive.
        try:
            t1 = time.time()
            builder = getattr(bc, builder_name)
            regen = builder()
            # Builders return a dict without schema_version / produced_utc;
            # match the convention used by ``write_cert`` in build_certified.
            regen.setdefault("schema_version", bc.SCHEMA_VERSION)
            elapsed = time.time() - t1
        except Exception as e:
            print(f"FAIL  {fname}: builder raised: {e}")
            failures += 1
            continue

        # Compare canonical bytes with produced_utc removed.
        regen_hash = manifest_hash(regen)
        disk_hash = manifest_hash(on_disk)
        if regen_hash != disk_hash:
            print(f"FAIL  {fname}: hash mismatch  ({elapsed:.1f}s)")
            print(f"        on-disk : {disk_hash}")
            print(f"        regen   : {regen_hash}")
            # Try to identify the first differing top-level key for a hint.
            for k in sorted(set(regen) | set(on_disk)):
                if k == "produced_utc":
                    continue
                if regen.get(k) != on_disk.get(k):
                    print(f"        first_diff_key: {k!r}")
                    break
            failures += 1
            continue

        print(f"OK    {fname}  [{regen['tier']}]  reproduces byte-equally  "
              f"({elapsed:.1f}s)")

    total = time.time() - t0
    if failures == 0:
        print(f"\n[reproduce_all] ALL {len(BUILDERS)} CERTIFICATES REPRODUCE  "
              f"(total {total:.1f}s)")
        return 0
    print(f"\n[reproduce_all] {failures} FAILURE(S)  (total {total:.1f}s)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
