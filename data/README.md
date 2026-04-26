# Data manifest

This directory contains upstream JSON artifacts imported from the Phase 9 source workspace.
The public paper does not claim that every historical Phase 9 artifact is certified here.
Instead, `certified/verify_all.py` is the source of truth for the paper-facing dependency set.

## Paper-facing certified sources

The following 10 files in this directory are referenced by the certified JSONs and are checked by `certified/verify_all.py` through the certificate source hashes:

- `phase_9_13_results.json`
- `phase_9_13f_certificate.json`
- `phase_9_14_results.json`
- `phase_9_15_results.json`
- `phase_9_15A_H_V1_identification.json`
- `phase_9_18_S9q_beta_lattice.json`
- `phase_9_18_S9q_C1.json`
- `phase_9_18_S9q_gamma1.json`
- `phase_9_18_S9x_bis_results.json`
- `phase_9_18_S9x_results.json`

The complete verified source set has 13 artifacts: these 10 files plus three recovered or cross-check JSONs under `certified/`.

## Supplemental archival sources

The following files are retained as supplemental historical artifacts from the source workspace. They are not direct source dependencies of the 11 certified paper certificates in the current release:

- `phase_9_13_S9_results.json`
- `phase_9_18_S9c_results.json`
- `phase_9_18_S9m_results.json`

They are included for provenance continuity, not as additional paper claims.