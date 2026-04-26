# SIBLING_NOTATION.md — notational reconciliation with the determinant-divisibility paper

> **Purpose.** One-page consistency check between the notation of paper II
> (this folder) and that of the sibling working draft
> `Determinant Divisibility of Centered Latin Squares`.
> This standalone note supersedes the internal pre-draft task file used
> during manuscript preparation.
>
> **Acceptance.** Every symbol that appears in both papers must denote
> the same mathematical object up to (a) trivially equivalent indexing
> conventions and (b) the relabeling decoration $\gamma$. Every
> divergence is listed below with its rationale.

---

## A. Side-by-side table

| Symbol | Sibling (Det.-divisibility) | Paper II (this paper) | Same object? | Notes |
|---|---|---|---|---|
| $L$ | $n \times n$ Latin square, $L_{ij} \in \{1,\dots,n\}$, indexing $0 \le i, j \le n-1$ | $n \times n$ Latin square; for $n=9$ also a Sudoku grid; indexing $1 \le i, j \le n$ | **Yes** | Pure indexing offset; both papers refer to the same matrix. |
| $m$ | $m = (n+1)/2$ | $m_n = (n+1)/2$ | **Yes** | Subscript added in paper II for clarity at variable $n$. |
| $J$ | $n \times n$ all-ones matrix | $J_n$ | **Yes** | Same. |
| $E$ | $E = L - mJ$ (centered matrix) | $E_{\gamma,L} = \gamma(L) - m_n J_n$ | **Yes** with relabeling | Paper II generalizes $E$ by composing with a symbol relabeling $\gamma \in S_n$; setting $\gamma = \mathrm{id}$ recovers the sibling's $E$. Sign convention is identical. |
| $\one$ | $\one = (1,\dots,1)^\top$, with $E \one = \mathbf{0}$, $\one^\top E = \mathbf{0}$ | identical | **Yes** | The two annihilator identities propagate to $E_{\gamma,L}$ since $\gamma$ permutes entries within each row/column block. |
| $V_{\mathrm{std}}$ | $V_{\mathrm{std}} = \ker(\one^\top) \subset \mathbb{R}^n$ (real vector space) | $V_{\mathrm{std},n} = \{ x \in \mathbb{Q}^n : \one^\top x = 0 \}$ (rational vector space) | **Yes**, with field choice | Paper II works over $\mathbb{Q}$ for SNF and rank purposes. The sibling works over $\mathbb{R}$ for determinant arguments. The hyperplane is the same; rank and kernel are field-stable here because $E_{\gamma,L} \in M_n(\mathbb{Z})$. |
| $A_{n-1}$ | not named explicitly (used as $V_{\mathrm{std}} \cap \mathbb{Z}^n$ implicitly in §6) | $A_{n-1} := V_{\mathrm{std},n} \cap \mathbb{Z}^n$ (root lattice) | **Yes** | Paper II names the integer hull as the root lattice $A_{n-1}$; sibling uses it without the name. |
| basis of $V_{\mathrm{std}}$ | $b_i = e_i - e_{n-1}$, $0 \le i \le n-2$ (last column as anchor, 0-based) | $P_n = (e_i - e_n)_{i=1}^{n-1}$ (last column as anchor, 1-based) | **Yes** | Same anchor (the last standard basis vector); only the index range differs by 1. |
| $P$ | $P = [b_0 \mid \cdots \mid b_{n-2}] \in \mathbb{R}^{n \times (n-1)}$ | $P_n$ (same matrix, 1-based) | **Yes** | Same matrix. |
| $G$ | $G = P^\top P = I_{n-1} + J_{n-1}$, $\det(G) = n$ | (used in spirit, not labelled) | **Same identity available** | Both papers exploit $G = I + J$ with $\det(G) = n$. Paper II does not need $G$ explicitly because it focuses on $\rank$ / SNF rather than $\det$. |
| Gram-projected matrix | $\widehat{E} = P^\top E P$, called $\Estd = GA$ in the sibling | $\widehat{E}_{\gamma,L} = P_n^\top E_{\gamma,L} P_n$ | **Yes** | Identical formula; paper II carries the $\gamma$ decoration. |
| $A$ (difference matrix) | $A_{ij} = L_{ij} - L_{i,n-1}$, $0 \le i, j \le n-2$, with $\widehat{E} = G A$ | not used | — | Paper II does not use the difference-matrix factorization; rank / SNF / kernel statements are made directly on $\widehat{E}_{\gamma,L}$. |
| $\det(\widehat{E}) = n \det(A)$ | Lemma "Gram factorization" (sibling §2) | not invoked | — | Paper II does not reduce to determinants. The two results are orthogonal: sibling computes $\det$, paper II computes $\rank$, $\ker$, and the stabilizer / lattice action. |
| $\rank$, $\ker$ | (not central to sibling) | $\rank(E_{\gamma,L}) = \rank(E_{\gamma,L} \mid_{V_{\mathrm{std},n}})$ (since the all-ones direction is always in the kernel) | — | Paper II's rank statements are independent of basis and equal to the rank of the Gram-projected $\widehat{E}_{\gamma,L}$ minus zero (no further drop), because $G$ is invertible over $\mathbb{Q}$. |
| SNF | not used | $\SNF(\widehat{E}_{\gamma,L})$ | — | Paper II default for SNF is $\widehat{E}_{\gamma,L}$, as locked in [`paper/notation.tex`](../paper/notation.tex). |
| $\Stab$, $\Aut$ | $\Aut$ refers to autotopism in informal language | $\Stab_{G^*}(v)$ for $v \in V_{\mathrm{std},n}$; $\Aut(L)$ = Latin-square autotopism group; $K_L$, $K_R$ = **Hessian symmetry groups** of $E_{\gamma,L}$ (definitely **not** stabilizers, **not** autotopisms) | **Distinct objects** | The sibling does not define $K_L$, $K_R$ — these are purely paper-II constructs. The disambiguation between $K_L / K_R$ vs $\Aut(L_{M_{19}})$ is exactly Theorem G of paper II. |

---

## B. Disambiguation flags (paper II must avoid implicit identification)

1. **Field choice.** Paper II uses $\mathbb{Q}$ throughout; the sibling
   uses $\mathbb{R}$. For all integer-matrix statements (rank, kernel,
   SNF), the choice is immaterial. Paper II must **not** invoke any
   real-analytic or topological argument (orthogonal projection,
   eigenvalue inequality, etc.) without explicitly switching to
   $\mathbb{R}$.

2. **Indexing.** Paper II indexes $1, \dots, n$; sibling indexes
   $0, \dots, n-1$. When formulas are quoted across papers, the index
   shift must be applied. The basis vector $b_i$ in the sibling
   corresponds to $e_{i+1} - e_n$ in paper II (column $i+1$ minus last).
   The matrix $P$ is the same matrix once indices are aligned.

3. **$E$ vs $E_\gamma$.** The sibling's $E$ is paper II's
   $E_{\mathrm{id}, L}$. When citing the sibling's results inside
   paper II, the $\gamma$ subscript must be **explicitly** absent
   ($\gamma = \mathrm{id}$); otherwise the sibling's identity is being
   misapplied.

4. **$\widehat{E}$ vs $\Estd$.** These are the same matrix. Paper II
   uses the hat to be visually compatible with $\widehat{*}$ for "Gram-
   projected" while keeping the symbol $E$ free for the un-projected
   centered operator.

5. **$K_L$, $K_R$ vs $\Aut(L)$.** These are **different** groups in
   paper II:
   - $K_L$, $K_R$ are subgroups of $S_n$ acting on rows / columns of
     the Hessian symmetry of the bilateral coset structure of
     $\Gamma^*_{M_{19}}$ (paper II §8).
   - $\Aut(L)$ is the Latin-square autotopism group of $L$ in the
     classical (Dénes–Keedwell) sense.
   Paper II Theorem G proves $|\Aut(L_{M_{19}})| = 1 \neq 72 = |K_L|$,
   forbidding identification of $K_L$ with the image of any
   autotopism subgroup. The sibling does not define either object.

---

## C. Verbatim translation cheat-sheet

| If the sibling writes … | Paper II writes … |
|---|---|
| "the Latin square $L$" | "the Latin (resp. Sudoku) square $L$" |
| "the centered matrix $E = L - mJ$" | "the centered operator at $\gamma = \mathrm{id}$, i.e.\ $E_{\mathrm{id}, L} = L - m_n J_n$" |
| "the Gram-projected matrix $\Estd$" | "$\widehat{E}_{\gamma=\mathrm{id}, L} = P_n^\top (L - m_n J_n) P_n$" |
| "$\det(\Estd) = n \det(A)$" | (cite as Lemma in sibling; not re-proved) |
| "$V_{\mathrm{std}}$" | "$V_{\mathrm{std},n}$" |
| "the basis $b_i = e_i - e_{n-1}$, $0 \le i \le n-2$" | "the basis $P_n = (e_i - e_n)_{i=1}^{n-1}$" (after $i \mapsto i-1$) |

---

## D. Acceptance checklist

- [x] Same definition of $L$ (modulo indexing).
- [x] Same definition of $E$ (paper II generalizes by $\gamma$).
- [x] Same hyperplane $V_{\mathrm{std}}$ (modulo field).
- [x] Same Gram-projected representative $\widehat{E} = P^\top E P$.
- [x] Sign conventions agree ($m = (n+1)/2$, $E = L - mJ$, $E\one = 0$,
  $\one^\top E = 0$).
- [x] No symbol collision with sibling's $A$ (paper II uses $A_{n-1}$
  for the root lattice — distinct from sibling's difference matrix $A$).
- [x] Disambiguation of $K_L, K_R$ vs $\Aut(L)$ explicit; cross-link
  to Theorem G of paper II.

Last refreshed: 2026-04-26 (P1.3 closure).
