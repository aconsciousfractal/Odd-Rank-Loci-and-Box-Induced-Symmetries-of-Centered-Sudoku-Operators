"""Phase 9.18.S9q-C1 — Open Question C1: enumerazione autotopismi di L_{M_19}.

Un autotopismo di un Latin square L (9×9) è una tripla (α, β, γ) ∈ S_9³
tale che L[α(i)][β(j)] = γ(L[i][j]) per ogni i,j.

Strategia di enumerazione (riferimento McKay-Meynert-Myrvold 2007):
  - Backtracking sui simboli γ ∈ S_9 (9! = 362880 candidati).
  - Per ogni γ fissato, calcolare se esistono (α,β) tali che la
    permutazione di simboli γ sia "indotta" da una riga/col-permutazione.
  - Per ogni γ, la matrice "γ-rinominata" L_γ[i][j] = γ(L[i][j]) deve
    essere ottenibile da L tramite riga × col permutazione.
  - Test: due Latin squares sono row-col-equivalent sse hanno la stessa
    classe di equivalenza sotto azione di S_n × S_n; per ogni γ
    proviamo il "rectangle matching".

Approccio efficiente: dato che L_M19 è uno Sudoku-grid, il gruppo
autotopismi di Sudoku è già limitato. Per la singola LS, |Aut| è
tipicamente ≤ qualche centinaio.

Implementazione: per ogni permutazione γ di simboli (9!), calcoliamo
L_γ. Poi per ogni α ∈ S_9 (9!) calcoliamo (α applicata a righe di L_γ)
e cerchiamo β tale che L_γ[α(i)] = L[i] modulo permutazione delle
colonne β. Costo nominale 9! × 9! = 1.3·10^11 — troppo.

Smart: **fissato γ, fissato α**, la β è determinata univocamente (se
esiste) dal vincolo che ogni riga deve mappare correttamente. Quindi
per ogni γ fixed, per ogni α fixed (9!), determiniamo β e verifichiamo.
Costo 9!² = 1.3·10^11 — ancora troppo. 

Ulteriore smart: iterando γ × α, β è UNICA per la prima riga (basta
che la prima riga dopo (α,γ) sia una permutazione di una riga di L,
e quella permutazione di colonne deve essere coerente sulle altre 8).
Quindi: per ogni (γ, α) fissati, β è determinata dalla riga 0; poi
verifichiamo le 8 righe restanti. Costo 9!² × 9 = 1.2·10^12. Ancora
troppo.

ULTIMATE smart: il numero di autotopismi di una LS specifica è
calcolabile in O(n!·n²) se si usa: per ogni γ (9!), build L_γ; per
ogni riga r di L (9 scelte), determina la permutazione α_r di righe
e β_r di colonne tale che la riga 0 di L_γ "diventa" riga r di L; poi
verifica.

Ma più drasticamente: per ogni γ ∈ S_9, calcolare L_γ. Il numero di
righe di L_γ uguali a qualche riga di L (come MULTISET di celle (col,
val)) determina rapidamente se γ è ammissibile. Se si, fissato γ,
esiste una BIJEZIONE α : righe(L) → righe(L_γ); le 9 scelte di
α(0)=r_0 determinano poi β.

Implementazione qui (forza bruta semplificata, sufficiente per LS 9×9):

  for γ in S_9 (9!):
      L_γ = γ ∘ L (renumber symbols)
      # cerco α,β tali che L_γ[α(i),β(j)] = L[i,j]
      # equiv: α envia righe di L su righe di L_γ, β colonne similmente
      try permutations α come matchings
      ...

Per L_M19 specifica, ci aspettiamo |Aut| modesto (< 10^4).

ALTERNATIVA: uso sympy + sage-like via gruppo permutazioni. In assenza
di pacchetti specifici, implemento direttamente.
"""
from __future__ import annotations
import json
import time
from itertools import permutations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]  # <repo>/scripts/legacy/ -> <repo>/
DATA = ROOT / "data"


def is_latin_square(L, n=9):
    L = np.asarray(L)
    for i in range(n):
        if set(L[i]) != set(range(1, n + 1)):
            return False
        if set(L[:, i]) != set(range(1, n + 1)):
            return False
    return True


def apply_autotopism(L, alpha, beta, gamma, n=9):
    """L'[i][j] = γ(L[α^{-1}(i), β^{-1}(j)])"""
    L = np.asarray(L)
    Lp = np.zeros_like(L)
    a_inv = [0] * n
    b_inv = [0] * n
    for i, v in enumerate(alpha):
        a_inv[v] = i
    for j, v in enumerate(beta):
        b_inv[v] = j
    for i in range(n):
        for j in range(n):
            Lp[i, j] = gamma[L[a_inv[i], b_inv[j]] - 1] + 1
    return Lp


def find_autotopisms_via_gamma_alpha(L, n=9, max_to_record=10000):
    """Enumera autotopismi (α,β,γ) di L.

    Strategia: per ogni γ ∈ S_n, per ogni mappa α : righe(L) → righe(L_γ),
    determina β unica e verifica.

    Complexity: O(n! · n! · n²). Per n=9 = 1.3·10^11 — troppo.
    Restringo: per ogni γ, applica γ ai simboli (L_γ).
    Calcola insieme di righe di L (come tuple) e di L_γ. Se intersezione
    vuota → skip γ. Altrimenti, per ogni assegnamento α(0) tale che
    riga 0 di L è "colonna-permutabile" alla α(0)-esima riga di L_γ,
    determina β univoca e verifica le 8 righe restanti.
    """
    L = np.asarray(L)
    rows_L = [tuple(L[i]) for i in range(n)]
    rows_L_set = set(rows_L)

    auts = []
    n_gamma_skipped = 0
    n_gamma_visited = 0
    t_start = time.time()

    # γ permutation: γ[v-1] is the new symbol for v
    for gamma in permutations(range(n)):
        n_gamma_visited += 1
        if n_gamma_visited % 36288 == 0:  # every 10%
            elapsed = time.time() - t_start
            print(f"    γ visited: {n_gamma_visited}/362880 "
                  f"({100*n_gamma_visited/362880:.0f}%) "
                  f"auts found: {len(auts)} "
                  f"elapsed: {elapsed:.1f}s", flush=True)
        # Compute L_γ rows
        L_gamma_rows = []
        for r in rows_L:
            new_r = tuple(gamma[v - 1] + 1 for v in r)
            L_gamma_rows.append(new_r)

        # Build per-row "column permutation type" via sorting
        # Two rows are col-permutable iff they are permutations of each other
        # since they're both Latin square rows, they're ALL perms of {1..n}.
        # So any row of L_γ is col-permutable to any row of L. Det. β
        # uniquely once we fix α(0) and α maps to a specific row.

        # For each candidate α(0) (= row r0 of L_γ that we want to be image
        # of row 0 of L), we DETERMINE β uniquely from row 0:
        # β(j) = position in L_γ_row[r0] of L[0,j]
        # Then verify each subsequent row.

        for r0 in range(n):
            target_row = L_gamma_rows[r0]  # this is L_γ[r0]
            source_row = rows_L[0]  # L[0]
            # We want L_γ[α(0), β(j)] = L[0, j], i.e.
            # target_row[β(j)] = source_row[j]
            # So β(j) = position in target_row where target_row equals
            # source_row[j].
            pos_in_target = {v: i for i, v in enumerate(target_row)}
            try:
                beta = tuple(pos_in_target[source_row[j]]
                             for j in range(n))
            except KeyError:
                continue

            # Now we must find α: 0 → r0 fixed; for each i, find α(i) =
            # row of L_γ such that L_γ[α(i), β(j)] = L[i, j] for all j.
            # I.e. for each row i of L, the row r of L_γ that — after
            # column permutation β^{-1} — equals L[i].
            # Equivalent: define γ_inv·L_γ-row = some row of L permuted
            # by β. So for each row r of L_γ, compute permuted-by-β:
            # apply_β(L_γ[r])[j] = L_γ[r, β(j)]
            # then this should equal L[i] for some i = α^{-1}(r)... i.e.
            # we're matching rows.

            valid = True
            alpha = [0] * n
            alpha[0] = r0
            used = {r0}
            for i in range(1, n):
                # Look for row r of L_γ such that L_γ[r, β(j)] = L[i, j]
                target_pattern = rows_L[i]
                found_r = None
                for r in range(n):
                    if r in used:
                        continue
                    perm_row = tuple(L_gamma_rows[r][beta[j]]
                                     for j in range(n))
                    if perm_row == target_pattern:
                        found_r = r
                        break
                if found_r is None:
                    valid = False
                    break
                alpha[i] = found_r
                used.add(found_r)

            if valid:
                aut = (tuple(alpha), tuple(beta), tuple(gamma))
                auts.append(aut)
                if len(auts) >= max_to_record:
                    print(f"    HIT max_to_record = {max_to_record}",
                          flush=True)
                    return auts, n_gamma_visited
    return auts, n_gamma_visited


def main():
    print("=" * 78)
    print("Phase 9.18.S9q-C1 — Aut(L_{M_{19}}) autotopismi espliciti")
    print("=" * 78, flush=True)

    src = DATA / "phase_9_18_S9x_results.json"
    with open(src, "r", encoding="utf-8") as f:
        d = json.load(f)
    L = np.array(d["M19_grid"])
    print(f"  L_M19 grid:")
    for row in L:
        print(f"    {row.tolist()}", flush=True)
    print(f"  is Latin square 9×9: {is_latin_square(L)}", flush=True)

    print()
    print("  Enumerating autotopisms via γ-α decomposition…", flush=True)
    print("  (worst case 9!·9 = 3.3M (α,γ)-tuples to check; "
          "early-exit per row mismatch)", flush=True)
    t0 = time.time()
    auts, n_visited = find_autotopisms_via_gamma_alpha(L)
    elapsed = time.time() - t0
    print()
    print(f"  Found |Aut(L_M19)| = {len(auts)} autotopismi", flush=True)
    print(f"  γ visited: {n_visited}", flush=True)
    print(f"  elapsed: {elapsed:.2f}s", flush=True)

    # Sanity: order divides 9!^3? Should divide |G|=9!·9!·9! but for a
    # specific LS, |Aut| << 9!^3 typically.
    fact9 = 362880
    print(f"  |Aut| / 9! = {len(auts)/fact9:.6f}", flush=True)
    print(f"  9! = {fact9}", flush=True)

    # Project to symbol-component γ
    gammas_used = set(a[2] for a in auts)
    alphas_used = set(a[0] for a in auts)
    betas_used = set(a[1] for a in auts)
    print()
    print(f"  distinct γ (symbol perms) in Aut: {len(gammas_used)}",
          flush=True)
    print(f"  distinct α (row perms) in Aut: {len(alphas_used)}",
          flush=True)
    print(f"  distinct β (col perms) in Aut: {len(betas_used)}",
          flush=True)

    # Verify: identity should be there
    e_id = tuple(range(9))
    assert (e_id, e_id, e_id) in set(auts), "identity missing!"
    print(f"  identity (e,e,e) ∈ Aut ✓", flush=True)

    # Symbol projection — what symbol permutations stabilize L?
    # This is the relevant subgroup for our K_L analysis.
    print()
    print("  Symbol-projection γ-component as sym-stabilizer:",
          flush=True)
    print(f"    distinct γ: {len(gammas_used)} (each may appear with "
          f"multiple (α,β))", flush=True)

    # Compare with K_L = 72 (affine coset symmetry)
    print()
    print(f"  K_L (from S9x) order: 72", flush=True)
    print(f"  |Aut(L_M19)| order:   {len(auts)}", flush=True)
    print(f"  ratio: {len(auts)/72:.3f}", flush=True)

    out = {
        "phase": "9.18.S9q-C1",
        "M19_grid": L.tolist(),
        "Aut_L_M19_order": len(auts),
        "n_distinct_gamma": len(gammas_used),
        "n_distinct_alpha": len(alphas_used),
        "n_distinct_beta": len(betas_used),
        "K_L_order": 72,
        "elapsed_s": elapsed,
        "n_gamma_visited": n_visited,
        # Sample of autotopismi
        "sample_autotopisms_first_20": [
            {"alpha": list(a), "beta": list(b), "gamma": list(g)}
            for a, b, g in auts[:20]
        ],
        "all_distinct_gammas": [list(g) for g in sorted(gammas_used)],
    }
    out_path = DATA / "phase_9_18_S9q_C1.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\n  saved → {out_path}", flush=True)


if __name__ == "__main__":
    main()
