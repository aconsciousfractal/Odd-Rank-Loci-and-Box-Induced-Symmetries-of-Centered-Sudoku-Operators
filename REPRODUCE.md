# Reproduce the paper artifacts

This guide records the exact commands and canonical artifact paths for

- `paper/main.tex`
- `paper/Odd-Rank-Loci-and-Box-Induced-Symmetries-of-Centered-Sudoku-Operators.pdf`

inside this repository. The repository is fully self-contained — no
external folder is read by any script in this package.

## Environment

- Python 3.11+
- `numpy`, `sympy`, `scipy`, `matplotlib`
- A working LaTeX distribution (MiKTeX 25.x or TeX Live 2025+) with
  `pdflatex` and `bibtex` on `PATH`

Install Python deps:

```bash
pip install -r requirements.txt
```

## Canonical artifact paths

| Artifact                                           | Canonical path                                       |
|----------------------------------------------------|------------------------------------------------------|
| 12 certificate JSONs                               | `certified/<result_id>.json`                         |
| Manifest of certificate sha256s                    | `certified/MANIFEST.sha256`                          |
| Recovered base1 22-perm list                       | `certified/base1_22_perms_recovered.json`            |
| Four-way base1 audit (rank/SNF Z + F2 / kernel)    | `certified/four_way_base1_22.json`                   |
| Recovered M19 100-perm list                        | `certified/M19_100_perms_recovered.json`             |
| 13 verified source artifacts                       | 10 under `data/`, 3 recovered/audit JSONs under `certified/` |
| 15 producer scripts                                | `scripts/legacy/*.py`                                |
| Paper figures (5 PDFs + producing scripts)         | `paper/figures/`                                     |
| Release paper PDF                                  | `paper/Odd-Rank-Loci-and-Box-Induced-Symmetries-of-Centered-Sudoku-Operators.pdf` |

## Core commands (run from repository root)

### 1. Rebuild the certificate package

```bash
python certified/build_certified.py
```

Reads the SHA-pinned source artifacts under `data/` and `certified/`,
plus `scripts/legacy/*.py` (sha-stamped provenance only — the upstream
data are already on disk), runs every claimed invariant in process
(Lift Lemma over $\mathbb{Q}$ via Sympy on
22 + 100 + 72 perms; Box-Band identity on Sudoku vs non-Sudoku; full SNF
of the cyclic LS-9 counterexample; etc.), and writes 12 certificate
JSONs plus `certified/MANIFEST.sha256`.

### 2. Verify the certificate package

```bash
python certified/verify_all.py
```

Performs four layers of checks:

1. Schema validation (required fields, allowed tiers, schema_version=1).
2. Manifest sha256 match against canonical serialization
   (with `produced_utc` removed).
3. Semantic spot-checks per certificate (rank counts, $H_{V_1}$ order,
   $z$-statistic, lift coverage, $|K_L| = 72$, $|\mathrm{Aut}(L_{M_{19}})| = 1$, …).
4. **Standalone-completeness audit**: every `scripts[].path` and
   `source[].file` referenced inside any certificate must exist inside
   this repository with sha256 matching the recorded provenance.

A passing run prints

```
[verify_all] ALL CHECKS PASSED (12 certificates, 18 scripts, 13 sources)
```

### 3. Re-derive every certificate byte-for-byte

```bash
python certified/reproduce_all.py
```

Re-imports `build_certified.py` and re-runs every builder, comparing
the canonical JSON serialization (with `produced_utc` removed) to the
on-disk certificate. A passing run prints

```
[reproduce_all] ALL 12 CERTIFICATES REPRODUCE
```

Runtime is machine-dependent; on the audit machine this is about 30-35 seconds, dominated by the middle-band balance certificate.

### 4. Recover the perm lists from upstream data (optional)

```bash
python certified/recover_base1_22_perms.py
python certified/recover_M19_100_perms.py
python certified/four_way_base1_22.py
```

These scripts re-derive `certified/base1_22_perms_recovered.json`,
`certified/M19_100_perms_recovered.json`, and
`certified/four_way_base1_22.json` from `data/phase_9_13_results.json`
and `data/phase_9_18_S9x_results.json`. They are exhaustive
$S_9 = 362\,880$ scans and take ~15 s each. After running them, rerun
`build_certified.py` to refresh `MANIFEST.sha256`.

### 5. Rebuild the paper figures (optional)

```bash
cd paper/figures
python fig1_sudoku_projectors.py
python fig2_rank_distribution.py
python fig3_D6_cayley.py
python fig4_A8_saturation.py
python fig5_M19_coset_structure.py
```

Each script writes the corresponding `.pdf` next to itself. Figure
regeneration uses `matplotlib>=3.7` from `requirements.txt`; the
scripts freeze PDF `CreationDate` and `ModDate` metadata so repeated
runs are byte-stable.

### 6. Compile the paper

From the repository root:

```bash
cd paper
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
cp main.pdf Odd-Rank-Loci-and-Box-Induced-Symmetries-of-Centered-Sudoku-Operators.pdf
```

Output: `paper/Odd-Rank-Loci-and-Box-Induced-Symmetries-of-Centered-Sudoku-Operators.pdf` (currently 22 pages). The intermediate `paper/main.pdf` is a local LaTeX build product and is ignored by git.

## Standalone test

To confirm the package does not silently rely on its current parent
directory, copy the entire repository to a fresh location with no
relation to its current location and rerun steps 1–3 + 6:

```bash
cp -r . /tmp/standalone-test/      # or equivalent on Windows
cd /tmp/standalone-test/
python certified/build_certified.py
python certified/verify_all.py
python certified/reproduce_all.py
cd paper && pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex && cp main.pdf Odd-Rank-Loci-and-Box-Induced-Symmetries-of-Centered-Sudoku-Operators.pdf
```

All checks should pass with output identical to the in-place run.

## Notes

- `certified/seeds/seedlist.txt` records the seeds used by the
  exhaustive scans where applicable.
- The producer scripts under `scripts/legacy/` retain their historical
  phase-numbered filenames so that the cross-imports inside the
  research log keep working without renaming. Their meaningful purpose
  is documented in [`scripts/SCRIPT_INDEX.md`](scripts/SCRIPT_INDEX.md).
- The certificate JSONs themselves are the canonical paper-facing
  artifact: every numerical claim in the paper points to the `result_id`
  of one of the 12 entries in `certified/MANIFEST.sha256`.
