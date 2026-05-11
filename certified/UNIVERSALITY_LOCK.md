# UNIVERSALITY_LOCK.md — §9 post-promotion wording lock

> **Purpose.** Lock the corrected §9 structure after the former
> universal dichotomy conjecture was promoted to a theorem. The paper
> now proves the all-$n$ stabilizer and lattice-saturation statements
> directly, while keeping the orbit-distance statement finite-range.

---

## §9.1 — EXHAUSTIVE block, $5 \le n \le 13$

**Proposition P4.** For every integer $n$ with $5 \le n \le 13$, the
certificate `stabilizers_5_to_13.json` verifies the stabilizer orders,
orbit sizes, lattice saturation, and orbit-distance statement

$$d_{M_n}(\mathrm{orb}_{G^*_n}(O_n))=2.$$

**Tier.** EXHAUSTIVE. This finite certificate remains useful as an
independent audit and as the only source for the orbit-distance claim.

---

## §9.2 — CERTIFIED finite analytic block, $n \in \{7,9,11\}$

**Certified Proposition F.** For $n \in \{7,9,11\}$, the certificate
`dichotomy_n7_n9_n11.json` records the abstract identifications

$$\mathrm{Stab}_{G^*_n}(M_n) \cong D_4 \times S_{n-4},$$

$$\mathrm{Stab}_{G^*_n}(O_n) \cong C_2 \times S_{n-4},$$

and the lattice saturation

$$\Lambda_{M_n}=\Lambda_{O_n}=A_{n-1}.$$

**Tier.** CERTIFIED_EXHAUSTIVE.

---

## §9.3 — THEOREM block, all $n \ge 5$

**Theorem H (General-$n$ dichotomy).** For every integer $n \ge 5$,

$$\mathrm{Stab}_{G^*_n}(M_n) \cong D_4 \times S_{n-4},$$

$$\mathrm{Stab}_{G^*_n}(O_n) \cong C_2 \times S_{n-4},$$

and

$$\Lambda_{M_n}=\Lambda_{O_n}=A_{n-1}.$$

**Proof basis.** Direct coordinate-value stabilizer analysis plus the
root-difference argument producing every $e_i-e_j$ from orbit-vector
differences. No exhaustive computation is needed for the theorem.

---

## Guard rails

- The all-$n$ theorem covers stabilizers and lattice saturation only.
- The all-$n$ theorem does **not** claim
  $d_{M_n}(\mathrm{orb}(O_n))=2$ for every $n \ge 5$.
- The distance statement remains certified only for $5 \le n \le 13$.
- Abstract, introduction, §9, discussion, and this lock must not refer
  to a universal dichotomy conjecture after this promotion.

Last refreshed: 2026-04-26 (affine-D4 / theorem promotion pass).
