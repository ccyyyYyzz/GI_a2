# RANK_LEDGER_TEST — gate CFG-A (Commutant–Fiber Gate, exact Schur ledger)

**Verdict: PASS** (CFG-A PASS, R2 PASS). Runtime ~9 s, CPU. Preregistered; frozen block written
before compute. Extends SYNTHESIS.md rows #4/#1 with the GPT R47 ruling
(`docs/ROUND63_GPT_ROUND47_RULING_RAW.md` Q2.2 / Q4.1 CFG-A).

## Claim

The stable object is the pair `(symmetry rep ρ, self-adjoint commutant A_G)`; the medium
computes exactly the invariant block-Gram algebra on the **base** and is exactly blind to the
**commutant gauge orbit** (hidden fiber). The rank ledger has three *distinct* rows:
`r_alg = dim A_G`, `d_base = rank D Φ_G(x)` (generic), `d_fiber = dim V − d_base`.

## Derivation (independent; agrees with GPT closed forms)

For the even cyclic group `Z_n` regular real rep over `m` channels — 2 real isotypic blocks
(`k=0`, `k=n/2`) + `n/2−1` complex-conjugate blocks, each multiplicity `m`:

```
r_alg (n,m) = 2·m(m+1)/2 + (n/2−1)m²  = (n/2)m² + m
d_base(n,m) = 2m + (n/2−1)(2m−1)      = nm − n/2 + 1
d_fiber(n,m)= dim V − d_base          = n/2 − 1   (one Fourier phase per complex block)
```

The slogan "rank = isotypic count" is exact **only** in the multiplicity-free case (`Z_24`,
m=1); the general invariant is the symmetric-commutant dimension above. My derivation and GPT's
`r_alg=(n/2)m²+m` agree, and both reproduce the synthesis integers.

## Frozen predictions vs measured (all EXACT integers)

| family | `r_alg` | `d_base` | `d_fiber` | measured (r_alg / d_base / d_fiber) |
|---|---:|---:|---:|---|
| `Z_24 ⊗ I_1` | 13 | 13 | 11 | 13 / 13 / 11 ✓ |
| `Z_12 ⊗ I_2` | 26 | 19 | 5  | 26 / 19 / 5  ✓ |
| `Z_8 ⊗ I_3`  | 39 | 21 | 3  | 39 / 21 / 3  ✓ |
| free `R^24`  | 300| 24 | 0  | 300 / 24 / 0 ✓ |

## CFG-A bars

- **item 1 (r_alg):** every measured rank == predicted, across 3 seeds and SVD-tol interval
  `[1e-6 … 1e-10]`. PASS.
- **item 2 (d_base):** generic base-Jacobian `rank[G_l x]` == 13/19/21/24, seed- and tol-stable. PASS.
- **item 3 (d_fiber):** `d − d_base` == 11/5/3/0 AND explicitly **constructed** equal-base states
  `x → exp(θ J_k) x` (per-complex-block Fourier-phase generators `J_k = (c_k s_kᵀ − s_k c_kᵀ)⊗I_m`)
  preserve **every** admissible response to numerical floor: max first-order response derivative
  ≤ 1.1e-14, max finite-angle (θ=0.7) response change ≤ 7.2e-15; constructed fiber-tangent rank
  matches. PASS.
- **item 4 (transverse):** an off-fiber unit perturbation changes at least one admissible response
  by O(1) (1.2–1.6). PASS.
- **item 5 (stability):** all integers invariant across the declared SVD-tol interval and seeds. PASS.
- **R2 (companion):** independent-reseed L-sweep on a `Z_24`-invariant quadratic task — accuracy
  knee at **exactly L = 13** (acc `L=12 → 0.897`, `L=13 → 0.963 = L=20` plateau). PASS.
- **exploratory (no bar):** each masked single-site diagonal-gain defect adds +1 admissible dof
  (13 → 14 → 15 → 16 → 17); one broken grain = one new invariant coordinate.

## Interpretation

The strong Noether complementarity wording is licensed: a symmetric disordered medium computes
its symmetry's invariant algebra on a base of dimension `d_base` and perfectly hides a continuous
`d_fiber`-dimensional commutant gauge orbit (the common Fourier phases). Because CFG-A passes in
full, the flagship may use the base/fiber complementarity claim, not merely the narrower
"commutant-rank selection rule."
