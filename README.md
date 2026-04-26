# Odd-Rank Loci and Box-Induced Symmetries of Centered Sudoku Operators

> Companion code, data, and certificates for the paper by O. Babanskyy (2026).
>
> The paper studies the rank-deficient locus of the centered Sudoku
> operator $E_\gamma = \gamma(L) - 5J$ acting on the standard subspace
> $V_{\mathrm{std}} = \{w \in \mathbb{Z}^9 : \sum w_i = 0\}$. It establishes
> four certified results: (i) a Sudoku-specific projector identity
> (Box Band Lemma), (ii) a Lift Lemma into the root lattice $A_{n-1}$,
> (iii) a $\{M_n, O_n\}$ orbit dichotomy with stabilizers
> $D_4 \times S_{n-4}$ vs $C_2 \times S_{n-4}$, and (iv) a Hessian
> symmetry $K_L \cong K_R \cong 3^2{:}Q_8$ at the laboratory base
> $M_{19}$, with $|\mathrm{Aut}(L_{M_{19}})| = 1$.

---

## Repository layout

```
.
├── paper/                       # LaTeX sources for the paper
│   ├── main.tex
│   ├── notation.tex
│   ├── refs.bib
│   ├── main.pdf                 # 21 pp build (regenerable)
│   └── figures/                 # 5 publication PDFs + producing scripts
├── certified/                   # Self-contained certificate package
│   ├── build_certified.py       # Builds 11 certificate JSONs from data/
│   ├── verify_all.py            # Schema + manifest + semantic + standalone audit
│   ├── reproduce_all.py         # Byte-identical re-derivation of every cert
│   ├── recover_base1_22_perms.py
│   ├── recover_M19_100_perms.py
│   ├── four_way_base1_22.py
│   ├── MANIFEST.sha256          # 11 result_id -> sha256 lines
│   ├── 11 *.json                # Certificate files
│   ├── base1_22_perms_recovered.json
│   ├── four_way_base1_22.json
│   ├── M19_100_perms_recovered.json
│   ├── README.md                # Tier conventions, schema, audit policy
│   ├── TIER_LEDGER.md
│   ├── SIBLING_NOTATION.md
│   ├── UNIVERSALITY_LOCK.md
│   └── seeds/seedlist.txt
├── scripts/
│   ├── SCRIPT_INDEX.md          # Cryptic-name → meaningful-purpose map
│   └── legacy/                  # 15 producer scripts (legacy phase_9_* names)
├── data/                        # 13 source JSONs ingested by the certificate builder
├── README.md                    # This file
├── REPRODUCE.md                 # Step-by-step reproduction guide
├── CITATION.cff
├── requirements.txt
├── .gitignore
└── LICENSE                      # MIT
```

The repository is **fully self-contained**: every file referenced by a
certificate JSON (provenance scripts and source data) lives inside the
repository and is checked by `certified/verify_all.py` against the
recorded sha256.

## Quick start

```bash
git clone https://github.com/aconsciousfractal/Odd-Rank-Loci-and-Box-Induced-Symmetries-of-Centered-Sudoku-Operators.git
cd Odd-Rank-Loci-and-Box-Induced-Symmetries-of-Centered-Sudoku-Operators

python -m venv .venv
.venv\Scripts\activate           # Windows
# source .venv/bin/activate      # Linux / macOS

pip install -r requirements.txt

# Verify all 11 certificates (schema + manifest + semantic + standalone audit)
python certified/verify_all.py

# Re-derive every certificate byte-for-byte
python certified/reproduce_all.py

# Rebuild the manifest (writes 11 cert JSONs and MANIFEST.sha256)
python certified/build_certified.py
```

## Dependencies

| Package | Version | Purpose                                            |
|---------|---------|----------------------------------------------------|
| numpy   | ≥ 1.24  | Array algebra, exhaustive S_9 and B_3 scans        |
| sympy   | ≥ 1.12  | Exact rank over ℚ, integer SNF, Gauss–Jordan/ℤ     |
| scipy   | ≥ 1.10  | `erfc` (two-proportion p-value), optional helpers  |

Tested with Python 3.11 / 3.12 / 3.14.

## Reproducing the paper results

See [`REPRODUCE.md`](REPRODUCE.md) for the full artifact map and exact
commands. The minimal three-step sequence is:

1. `python certified/build_certified.py` — regenerates 11 certificate
   JSONs and `MANIFEST.sha256` from `data/` and `scripts/legacy/`.
2. `python certified/verify_all.py` — schema validation, manifest
   sha256 check, semantic spot-checks, and **standalone-completeness
   audit** (every script and source file referenced in any certificate
   must exist inside the repository with matching sha256).
3. `pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex`
   inside `paper/`.

## Certificate tier discipline

Each certificate JSON declares an explicit tier:

- **THEOREM** — proved on paper, audited code witness.
- **THEOREM_EXHAUSTIVE** — proved by exhaustive enumeration over a
  finite group (e.g. $9! = 362\,880$ relabelings), no Monte-Carlo.
- **EXHAUSTIVE** — claim is statistical-style but the population is
  finite and fully scanned.
- **EMPIRICAL** — observed phenomenon used as evidence; the paper
  does not lean on it for a theorem.
- **OBSERVED TEMPLATE** — pattern observed in scope; restricted use.

The full ledger is in [`certified/TIER_LEDGER.md`](certified/TIER_LEDGER.md).

## Provenance audit

`certified/verify_all.py` runs a **standalone-completeness audit**
that, for every certificate, checks:

1. All required schema fields are present.
2. The recorded sha256 in `MANIFEST.sha256` matches the canonical
   serialization of the certificate (with `produced_utc` removed).
3. Each declared semantic invariant (rank counts, group orders,
   z-statistic, lift coverage, …) re-derives correctly.
4. **Every** `scripts[].path` and `source[].file` referenced by the
   certificate is a real file inside this repository, with sha256
   matching the recorded provenance.

A passing run prints
`ALL CHECKS PASSED (11 certificates, N scripts, M sources)`.

## Cryptic legacy script names

The historical producer scripts under `scripts/legacy/` carry phase
labels (`phase_9_13.py`, `phase_9_18_S9x.py`, …) inherited from the
research log. Their meaningful purpose is documented in
[`scripts/SCRIPT_INDEX.md`](scripts/SCRIPT_INDEX.md). The certificate
package never depends on the names directly: every dependency is
recorded as a sha256-stamped provenance entry.

## License

MIT. See [`LICENSE`](LICENSE).
