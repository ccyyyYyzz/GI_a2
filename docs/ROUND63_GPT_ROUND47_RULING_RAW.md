# R47 — The flagship: disorder as a quotient machine

**Reference brief:** [`docs/ROUND63_GPT_ROUND47_QUESTION.md`](https://github.com/ccyyyYyzz/GI_a2/blob/ae2331a/docs/ROUND63_GPT_ROUND47_QUESTION.md) at [`ae2331a`](https://github.com/ccyyyYyzz/GI_a2/commit/ae2331abed1e7064c43c3969d4f58c3403bccc4c).  
**Merged divergence:** [`results/round63_next/FIVE_FRONT_DIVERGENCE/SYNTHESIS.md`](https://github.com/ccyyyYyzz/GI_a2/blob/ae2331a/results/round63_next/FIVE_FRONT_DIVERGENCE/SYNTHESIS.md).  
**Two-paper cap:** binding. The frozen PRL Letter remains Paper 1. This issue governs the only remaining flagship.

## Executive ruling

> **Conditional freeze on a Noether-complementarity protagonist, with the Fog Bundle as the narrative spine.**
>
> **Symmetry turns a disordered medium into a quotient machine: its quadratic record exposes the block-Gram invariants permitted by the symmetry representation and is exactly blind to the complementary unitary-orbit coordinates; the same base can carry calibration-free information while motion on the hidden fibers pays the curvature tax.**

This is candidate **(b)**, sharpened. The curvature tax is not demoted; it becomes the dynamics of the hidden fiber. The passed Fog quotient channel becomes the positive information current on the base. The QUAD kill becomes the mechanism-layer discovery: the physically admissible computer did not fail randomly—it collapsed to the exact self-adjoint commutant of the stationary symmetry.

The architecture is conditional on one named gate:

> **COMMUTANT–FIBER GATE (CFG).** The parameter-free algebra ranks, generic base dimensions, and hidden-fiber dimensions in §2.4 must all land exactly, and the joint Fog Bundle run must simultaneously retain quotient capacity and zero conditional fiber information.

If CFG passes, this is the flagship and the preferred title is **Disorder as a Quotient Machine**. If any exact integer or zero-information bar fails, Noether complementarity is not softened into a vague slogan: the paper falls back cleanly to a curvature-tax-centered Optica article with the passed Fog quotient channel as its second pillar.

The operative correction to the merged synthesis is equally important:

> **Say yes to one representation-theoretic ledger; say no to one numerological rank.**

The QUAD cone rank `13`, the generic quotient dimension, the hidden-fiber dimension, the Fog `M_eff=12.15`, and the coherence crossover plateau are different functionals of one factorized experiment. Some may agree in a special multiplicity-free cell, but they are not definitionally equal.

---

# Q1 — Protagonist and spine

## Q1.1 Ranking

1. **Noether complementarity as theorem; Fog Bundle as spine. — FIRST, conditional on CFG.**  
   It absorbs the QUAD negative, the Fog channel, the CBT moonshot, and the coherence/certification ledgers without reading like five campaigns.

2. **Curvature tax as headline; Fog Bundle as spine. — CLEAN FALLBACK.**  
   The moonshot is genuinely new and has survived 22/22 checks plus a second non-optical model. It is the fallback if the exact commutant/fiber conjecture misses. It is not first choice because it leaves the rank-13 QUAD discovery and the saturating quotient channel as attached chapters rather than consequences of one mechanism.

3. **Finite-`p` exponent-flow / 0.454 closure. — NOT A PROTAGONIST.**  
   The closure is excellent science and belongs in the finite-data chapter, but the noncentral-chi-square boundary is classical. Its value is that it removes an anomaly, not that it creates the flagship idea.

## Q1.2 Dinner-table sentence

> **A symmetric fog performs a quotient of reality: it transmits the invariant coordinates and perfectly erases the conjugate phases; on that same record, hidden motion is free only along the fiber and pays a curvature tax when it tries to leave.**

A more technical one-sentence thesis for the manuscript is:

> **For a symmetry-complete quadratic medium, the observable scene coordinates are the block Gram matrices of the self-adjoint commutant, while the complementary right-unitary coordinates carry zero conditional information; complete scrambling then realizes positive calibration-free capacity on the quotient and exact dynamical cloaking on the fiber.**

The qualifier **symmetry-complete quadratic medium** is load bearing. This is not a claim about all polynomial invariants, all moments, or every single fixed medium realization.

## Q1.3 The theorem the paper actually needs

Let a compact group `G` act orthogonally on a real scene space `V` through `rho(g)`. Define the self-adjoint commutant

```text
A_G = { K in Sym(V) : K rho(g) = rho(g) K for every g in G }.    (1.1)
```

A stationary Gaussian or second-order medium law can expose quadratic scene functionals

```text
q_K(x)=x^T K x,   K in A_G.                                    (1.2)
```

Assume the physically admissible kernel family spans `A_G`. Decompose the real representation into isotypic blocks

```text
V ~= direct_sum_alpha [ F_alpha^(m_alpha) tensor_F W_alpha ],   (1.3)
```

where `F_alpha` is `R`, `C`, or `H` according to the real, complex, or quaternionic representation type, and `m_alpha` is the multiplicity. Schur’s lemma gives

```text
A_G ~= direct_sum_alpha Herm_{m_alpha}(F_alpha) tensor I_{W_alpha}.  (1.4)
```

Write the scene block as a matrix `X_alpha` over `F_alpha`. Then every admissible quadratic response has the form

```text
q_B(x)=tr[B X_alpha X_alpha^dagger],                            (1.5)
```

so the complete quadratic record depends on the scene through the block-Gram map

```text
Phi_G(x) = ( X_alpha X_alpha^dagger )_alpha.                    (1.6)
```

### Commutant–Fiber Complementarity Theorem

Under the spanning assumption above:

1. **Computed base.** The vector of all admissible quadratic responses is equivalent to `Phi_G(x)`; no other quadratic scene coordinate is available.
2. **Hidden fiber.** On each regular-rank stratum, two scenes have the same base exactly when their isotypic matrices differ by right orthogonal/unitary/symplectic transformations on the irreducible factors, up to the usual stabilizers and discrete signs:

   ```text
   Y_alpha = X_alpha U_alpha,
   U_alpha in O/U/Sp.                                          (1.7)
   ```

3. **Information complementarity.** If the full record law factors through this base,

   ```text
   X -> Phi_G(X) -> Y,                                         (1.8)
   ```

   then for every prior and every record length,

   ```text
   I(fiber coordinate ; Y | Phi_G(X)) = 0,                     (1.9)
   I(X;Y)=I(Phi_G(X);Y).                                       (1.10)
   ```

4. **Capacity reduction.** The channel capacity over scenes equals the capacity over the realizable quotient base:

   ```text
   sup_{P_X} I(X;Y) = sup_{P_B, B in Phi_G(V)} I(B;Y).          (1.11)
   ```

### Proof sketch

Equation (1.4) is the self-adjoint part of the commutant decomposition. Equation (1.5) shows that evaluating every `B` is equivalent to knowing the Gram matrix `X_alpha X_alpha^dagger`. Equality of Gram matrices gives the right-unitary relation by polar decomposition, with partial isometries on rank-deficient strata. The Markov and mutual-information identities then follow immediately.

### Scope correction that must appear in the paper

For a finite group such as `Z_24`, the hidden fiber is generally **larger than the original finite translation orbit**. The stationary quadratic law hides the phase-orientation group generated by the commutant decomposition. In the multiplicity-free cyclic case this is the Fourier-phase torus (plus discrete sign/reflection choices), not merely the 24 translated copies. The paper must call this the **commutant gauge orbit**, not casually identify it with the original `Z_24` orbit.

This precision is what prevents the flagship from collapsing into “the power spectrum loses phase,” which is classical.

## Q1.4 What is new versus classical

The ingredients—Schur decomposition, commutants, Gram factorization, power-spectrum phase loss, and group-invariant statistics—are classical. The paper may not sell Schur’s lemma as a new theorem.

The candidate new physical package is the closed square:

> **one declared disordered record + parameter-free symmetry rank + positive quotient capacity + exact fiber privacy + curvature-governed hidden dynamics.**

The passed measurements make that square unusually sharp:

- physically admissible quadratic laws collapsed to rank `13`, not the generic `300`;
- the Fog ratio channel retained `2.3219` bits while the equal-DC mean detector remained at AUC `0.50` under 40 dB medium variation;
- the CBT great-circle cloak stayed at AUC `0.5` through `16384` banks while radial and under-actuated controls were detected;
- the three finite-sample ledgers were separately validated.

The theorem is the mechanism; the simultaneous positive/zero-information physics is the contribution.

## Q1.5 Preferred titles

1. **Disorder as a Quotient Machine** — preferred after CFG.
2. **Symmetry Makes Disorder Compute Invariants and Hide Orbits**
3. **The Fog Bundle: Capacity on the Base, Zero Information on the Fiber**
4. **What Symmetric Disorder Computes—and Must Forget**
5. **Noether Complementarity in Disordered Sensing** — technically safe, editorially colder.

Fallback title if CFG fails:

> **The Curvature Tax of Undetectability in Scrambled Measurements**

---

# Q2 — The One-Rank-Object bet

## Q2.1 Binding ruling

> **Architect around one representation object, not one scalar rank.**

The stable object is the pair

```text
(symmetry representation rho, self-adjoint commutant A_G),      (2.1)
```

plus the observation map from its quotient coordinates into the record. The paper should carry a **rank ledger with distinct rows**, not equate every measured `~13`.

The minimum ledger is:

1. `r_alg = dim A_G` — number of linearly independent admissible quadratic kernels;
2. `d_base = rank D Phi_G(x)` on the generic stratum — number of locally independent quotient coordinates;
3. `d_fiber = dim V - d_base` — dimension of the continuously hidden state fiber;
4. `g_rec` or `M_eff` — a real-valued record-space information conductance, such as `sum tau_j^2`, not an integer algebra rank;
5. `q_coh` — the coherence/scene/code bottleneck that caps how many quotient directions can reach the record.

The numerical coincidence `M_eff=12.15≈13` is **not** allowed to carry the architecture. `M_eff` is the Fisher-whitened number of effective covariance looks for estimating a chosen quotient coordinate. `13` is the dimension of a stationary quadratic kernel space in a different representation. They may share an upstream symmetry limitation, but equality is not automatic.

## Q2.2 Exact cyclic predictions

For the regular real representation of an even cyclic group `Z_n`, repeated over `m` channels, the self-adjoint commutant contains two real isotypic blocks and `n/2-1` complex-conjugate blocks. Therefore

```text
r_alg(n,m)
 = 2 * m(m+1)/2 + (n/2-1)m^2
 = (n/2)m^2 + m.                                                (2.2)
```

On the generic nonzero stratum, the base consists of two real rank-one Gram matrices and `n/2-1` complex rank-one Gram matrices, so

```text
d_base(n,m)
 = 2m + (n/2-1)(2m-1)
 = nm - n/2 + 1,                                                (2.3)

d_fiber(n,m)=n/2-1.                                            (2.4)
```

The exact predictions are:

| family | `r_alg` | `d_base` | continuous `d_fiber` |
|---|---:|---:|---:|
| `Z_24 ⊗ I_1` | **13** | **13** | **11** Fourier phases |
| `Z_12 ⊗ I_2` | **26** | **19** | **5** common Fourier phases |
| `Z_8 ⊗ I_3` | **39** | **21** | **3** common Fourier phases |
| free `R^24` quadratic cone | **300** | **24** | **0** continuous dimensions, only the global sign if signed scenes are allowed |

This reveals why the raw “one rank” conjecture is too coarse. The exact-rank tests `13/26/39/300` probe `r_alg`. They do **not** by themselves prove the base dimensions or fiber dimensions in the last two columns.

## Q2.3 The structural bet and its gate

CFG must therefore score both levels:

- exact matrix-span ranks `13/26/39/300`;
- exact generic Jacobian ranks `13/19/21/24`;
- exact phase-fiber dimensions `11/5/3/0`, checked by constructing equal-base states and verifying full-law equality;
- stability of the integer results to numerical threshold and independent seeds.

If all land, the paper may state:

> **The symmetry representation predicts, before simulation, both how many quadratic coordinates the medium can compute and how many continuous scene coordinates it must hide.**

If only `r_alg` lands, use the narrower phrase “commutant-rank selection rule”; do not claim full base/fiber complementarity.

## Q2.4 Disposition of the five measured numbers

- **QUAD rank 13:** primary evidence for `r_alg`.
- **Blind-fiber codimension:** predicted by `d_base`, not by `r_alg` in general; must be measured separately.
- **Fog `M_eff=12.15`:** record conductance; retain as a measured Fisher quantity, not an integer rank.
- **Coherence plateau:** upstream bottleneck `q_coh≈R^2`, already independently recovered to 4.8%; it belongs in the rank-ledger section but is not forced equal to the commutant rank.
- **Capacity prelog/alphabet capacity:** an information-theoretic functional of the base channel. The current result saturates a five-symbol alphabet; it does not measure `r_alg`.

Frozen shorthand:

> **One object, several ranks; no numerology.**

---

# Q3 — Chapter map

Target the PRX form first: approximately **9–11 main-text pages, four main figures, no chronology table**, with a focused but deep supplement. The same scientific unit can be compressed to an Optica article if needed.

## Main-text element 1 — Disorder as a quotient machine

**Single claim:** symmetry partitions a disordered record into a computable base and a perfectly hidden fiber.

Content:

- the one-sentence paradox;
- the declared quadratic/stationary scope;
- one schematic of the physical medium and the quotient map;
- cite the frozen Letter for Rank–Jet Separation rather than reproducing its derivation or figures.

**Figure 1 — hero:** a stratified Fog Bundle. The scene space is drawn as fibers above a base of block-Gram invariants. The record arrow carries `2.3219 bits` on a selected base coordinate and `0 bits` conditionally along the fiber. A curved invisible trajectory lies on a fiber; a chord leaves it and pays the tax.

## Main-text element 2 — Commutant–Fiber Complementarity

**Single claim:** the symmetry representation predicts the computed quadratic algebra and the hidden gauge coordinates.

Content:

- Eqs. (1.1)–(1.7) in compressed form;
- cyclic closed-form ranks;
- exact `13/26/39/300` and `13/19/21/24` tests if CFG passes;
- the QUAD general-computer kill appears here as the positive mechanism result: arbitrary quadratic computing dies because stationary physics restricts the medium to the commutant.

**Figure 2:** representation blocks and exact rank ledger. It should show predicted versus measured integers—not a gallery of four campaigns.

## Main-text element 3 — Capacity on the base, zero information on the fiber

**Single claim:** one complete-scrambling record supports a nontrivial calibration-free quotient channel while revealing no scene coordinate within the selected fiber.

Content:

- exact F-ratio channel;
- 40 dB nuisance cancellation;
- absolute decoder collapse;
- mean AUC at chance for equal-DC symbols;
- alphabet capacity `log2(5)=2.3219` bits;
- the joint Fog Bundle conditional-MI result after CFG.

**Figure 3:** one two-sided figure: left, base capacity/BER under medium swings; right, fiber classifier/conditional MI at chance. The visual sentence is **positive capacity and zero privacy leakage in the same record**.

## Main-text element 4 — Dynamics on a hidden fiber

**Single claim:** static fiber privacy has a dynamical law—the curvature tax—and cross-fiber transfer pays an endpoint information toll.

Content:

- CBT acceleration tax and great-circle equality trajectory;
- 22/22 moonshot outcome;
- second non-optical model;
- endpoint-toll theorem from Q5 below;
- local versus global and anchor conditions.

**Figure 4:** curvature-tax kink and exact cloak, plus a small endpoint-toll inset. Do not reuse the PRL Letter’s jet-slope figure.

## Main-text element 5 — Algebraic rank, record conductance, and finite-data certification

**Single claim:** what exists, what reaches the record, and what can be certified are different bottlenecks.

Content:

- coherence crossover `sqrt(q_fit)=17.99` versus `R=18.9`;
- `M_eff=12.15` as conductance rather than algebra rank;
- TLSG `1:sqrt(p):p` as a classical but experimentally validated finite-data mechanism;
- the `0.4538` exponent is closed by the exact finite-`p` noncentral-chi-square boundary and is not a new exponent.

This section gets no standalone conceptual figure unless typesetting requires one; the relevant curves can be a compact lower panel of Fig. 2 or supplement.

## Main-text element 6 — Physical scope and boundaries

**Single claim:** the quotient principle survives physical optics only inside declared symmetry and nuisance conditions.

Content:

- wave-optics and COMSOL scope in one paragraph;
- medium amplitude must be cancelled or anchored;
- stationary quadratic law and Gaussian/second-order scope;
- D3/D5 specificity dent remains disclosed; IT7 hardware veto remains killed;
- no claim that one fixed medium spans the full commutant;
- no higher-order completeness claim.

## What is absorbed from `paper2_optica`

Main text receives only:

- one concise physical wall/commutant motivation;
- the fresh-calibration transfer point as evidence that a declared fiber can be certified;
- jet transmutation as a corollary explaining how hardware leakage moves a trajectory off the ideal fiber.

Everything else from the completed certification manuscript moves to supplement:

- full Statistical Noether proposition and wall periodic table;
- blindness-capacity SVD curves, DPSS alternative, noncomposability;
- TLSG derivation and bootstrap;
- complete wave twin and COMSOL ledger, including the unfavorable grain discrepancy;
- KT/IT chronology and all killed hardware routes;
- sealed D3/D5 tables and IT7 negative;
- calibrated leak law and PWM apparatus requirement;
- all figure-source custody tables.

Demotion is not deletion. Under the two-paper cap the supplement is the archive; the main text is the quotient argument.

---

# Q4 — Minimal pre-freeze decisive runs

## Q4.1 Named gate: COMMUTANT–FIBER GATE

The architecture must not freeze before the following three items land.

### CFG-A — exact Schur ledger (`RANK_LEDGER_TEST`, already running)

**Budget:** CPU minutes; effectively `0 GPU-day`.

Binding bars:

1. algebra ranks exactly `13/26/39/300`;
2. generic base-Jacobian ranks exactly `13/19/21/24`;
3. constructed fiber motions preserve every admissible quadratic response to numerical floor;
4. a perturbation transverse to each fiber changes at least one response;
5. numerical rank stable over a declared singular-value threshold interval and independent seeds.

Any miss kills the strong complementarity wording. No tolerance around an integer is permitted except floating-point rank certification.

### CFG-B — joint Fog Bundle run

**Budget ceiling:** `0.75 L4-GPU-day`.

Use the same complete-scrambling record for the base symbols and CBT orbit schedules. Frozen bars:

1. five-symbol quotient capacity within `0.02 bit` of the separately validated `2.3219-bit` result;
2. conditional fiber information upper confidence bound `<0.01 bit`;
3. orbit classifier performance statistically indistinguishable from `1/M` for every tested bank count and exact-cloak schedule;
4. quotient decoding invariant to the existing 40 dB medium-amplitude sweep;
5. Fog-derived efficient information predicts the CBT paired-`d'` coefficients within `10%` without refitting.

This is the decisive “same record, two ledgers” experiment. Separate positive and privacy experiments are not enough.

### CFG-C — FOG degree-of-freedom/floor audit

**Budget ceiling:** `0.20 L4-GPU-day`.

Required because `12.15≈13` is currently an attractive coincidence.

- compare the exact weighted-quadratic-form law with the effective-chi-square approximation across `T`;
- report the full visibility spectrum, `tr K`, `tr K^2`, Satterthwaite degree, and exact Imhof/Davies tail;
- require capacity prediction within `0.02 bit` and calibrated tail error within the predeclared probability bar;
- test the bucket-mean flicker and darkness floor named in the synthesis.

The spine survives if the quotient channel remains positive and fiber privacy remains zero. Only the “one-rank corollary” dies if the exact dof is not 13.

### NCX2 finite-boundary closure

**Budget:** CPU / negligible GPU; already running.

This is required only to close the finite-data section:

- reproduce the measured `0.4538` exponent from the exact 50%-power boundary;
- hit the preregistered FPR `0.01` slope and the `p=2080` local slope.

Failure removes the cross-campaign exponent-normal-form paragraph; it does not kill the quotient architecture.

## Q4.2 Runs explicitly not required before freeze

- **No `t^6` jet rung.** Beautiful, but the CBT result is already complete enough and a third rung does not decide the spine.
- **No detector tournament.** CBT tests the full record law, exact divergences, Monte Carlo, CUSUM, and a second model. Another zoo is defensive accumulation.
- **No serial-bottleneck 3×3 grid before freeze.** The current coherence gate already establishes a measurable crossover. The 3×3 grid becomes supplement work only if the rank section needs a cleaner harmonic-composition figure.
- **No TLSG F-ratio rebase.** Ruled separately in Q6.
- **No 13-lane spectral MIMO run as an architecture gate.** Run it only after freeze if CFG-C shows that a multi-lane capacity corollary can be stated without conflating algebra rank and record conductance.

## Q4.3 Total budget

Mandatory hard ceiling:

```text
CFG-A + NCX2          CPU / negligible
CFG-B                 0.75 L4-GPU-day
CFG-C                 0.20 L4-GPU-day
contingency           0.30 L4-GPU-day
TOTAL                 1.25 L4-GPU-days maximum
```

There is no repair round. Once these results land, the architecture freezes.

---

# Q5 — The endpoint toll

## Q5.1 Ruling

> **A precise endpoint toll is provable. The broader “cloak is a prison” recurrence slogan is false and must be retired.**

Exact invisibility confines a path to one quotient fiber. But a path can move to another point on the same connected fiber and stop there with zero information toll; it need not be recurrent. The great-circle cloak is periodic only because it is the constant-speed minimum-acceleration equality trajectory.

What is true is:

> **Within-fiber transport is free; crossing fibers incurs an endpoint divergence that no schedule can erase.**

## Q5.2 Endpoint Divergence Toll Theorem

Let the record law factor through a quotient coordinate `b=Phi(x)`. Let the null state have base `b0`. A controlled schedule produces a full record law `P_gamma` and ends at a state with base `b1`. Assume the endpoint is observed for `H` independent banks. Let `P_0` be the full record under the constant null schedule.

For every Rényi order `alpha in (0,1)` and hence for Chernoff information,

```text
D_alpha(P_gamma || P_0)
 >= H D_alpha(P_{b1} || P_{b0}),                               (5.1)

C(P_gamma,P_0)
 >= H C(P_{b1},P_{b0}).                                       (5.2)
```

### Proof

Project the full record onto its final `H` banks. Rényi divergence obeys data processing, so the divergence of the full record is at least the divergence of that marginal. The endpoint marginal is `P_{b1}^{⊗H}` under the schedule and `P_{b0}^{⊗H}` under the null, so Rényi divergences scale by `H`. Chernoff information is the supremum over the corresponding order-`alpha` exponents, preserving the inequality. `square`

### Equality condition

Equality is achieved when every pre-endpoint observation remains on the initial fiber and the only law change is the final `H`-bank hold at `b1`. Therefore the bound is sharp if instantaneous cross-fiber transfer is permitted. Dynamics, speed, or acceleration limits can add a schedule-dependent path toll, but no universal extra term exists without specifying those constraints.

### Gaussian `Q` specialization

For `Sigma(q)=R+qH`, the exact endpoint toll is

```text
H * max_{0<=s<=1} 1/2 [
  s log det Sigma0 + (1-s) log det Sigma1
  + log det( s Sigma0^{-1} + (1-s) Sigma1^{-1} )
].                                                               (5.3)
```

For a small permanent quotient shift `Delta q`,

```text
C_endpoint
 = H (I_q/8) (Delta q)^2 + o((Delta q)^2).                      (5.4)
```

## Q5.3 Frozen wording

Use:

> **Exact cloaking permits arbitrary transport within a connected blindness fiber, but any transfer to a different quotient value pays a sharp endpoint information toll: if the new state is held for `H` banks, the full-record Chernoff information is at least `H` times the endpoint-law Chernoff information, independent of the preceding schedule.**

Do not use:

- “all invisible motion is recurrent”;
- “every permanent transfer pays” without saying **cross-fiber**;
- “schedule-independent control cost.”

The toll is an information lower bound, not a mechanical work law.

---

# Q6 — Ownership rulings

## Q6(a) The `1.6–1.8` pre-crossover slope

> **Owned by this flagship’s rank/finite-size section, not by a separate memory-effect paper.**

The two-paper cap eliminates a third publication unit. The present evidence says the slope is a finite-crossover chord, not a new universal exponent:

- the harmonic code/coherence bottleneck already predicts downward curvature before the nominal `M*`;
- the exact finite-`p` normal form has closed the analogous TLSG anomaly;
- the named coherence gate is already passed by the primary rank agreement, not by the soft slope bar.

Frozen handling:

> **Report `1.61` as a finite-window effective slope of the measured crossover and fit the full crossover law; do not call it a coherence-spectrum exponent or spectroscopy observable unless a preregistered residual remains after the parameter-free normal form.**

A synthetic beta sweep may be run after freeze for the supplement. It does not control the architecture.

## Q6(b) Re-basing TLSG on F-ratio quotients

> **Do not unfreeze TLSG.**

TLSG has a clean 7/7 frozen result, but its theorem skeleton is classical and belongs in the mechanism layer. Re-basing its arms on F-ratios would:

- change a completed certification experiment;
- blur ownership between the Fog base channel and the finite-data certificate;
- open a separate drift-proof certification program that the two-paper cap cannot support.

Ownership is frozen as follows:

- **Fog Bundle owns quotient cancellation and F-ratio communication.**
- **TLSG owns the original `1:sqrt(p):p` finite-certification mechanism and fresh-calibration transfer.**
- The flagship may state in Discussion that quotient statistics are a future route to drift-resistant certificates, but no new TLSG arm or headline is authorized.

---

# Q7 — Venue

## Q7.1 Conditional first submission

> **PRX first only if CFG passes in full. Optica first if it does not.**

PRX explicitly seeks fundamental discoveries, paradigm shifts, and new connections across physics and adjacent disciplines ([official scope and criteria](https://journals.aps.org/prx/about)). The full CFG-positive package can plausibly meet that bar because it would contain:

- a general compact-group/commutant base–fiber theorem;
- parameter-free exact ranks across multiple representation families;
- simultaneous positive quotient capacity and zero fiber information in one record;
- a validated dynamical curvature tax on the hidden fiber;
- a non-optical replication;
- a physical scattering realization rather than only abstract invariant theory.

Without the exact rank/fiber closure, the package is still strong optics but no longer a clear physics paradigm. Optica is then the honest first target because it is explicitly a high-impact venue for fundamental and applied optics and photonics ([official scope](https://opg.optica.org/content/journal/about/item/optica/)).

## Q7.2 Ladder

### CFG passes

```text
PRX -> Optica -> Physical Review Applied
```

### CFG misses but CBT + Fog Bundle remain intact

```text
Optica -> Physical Review Applied -> Optics Express
```

Do not submit the flagship until the frozen PRL Letter is public as a preprint and is cited explicitly. Paper 1 owns Rank–Jet Separation and the basic scrambled `Q` channel; Paper 2 owns commutant complementarity, quotient capacity, curvature-controlled fiber dynamics, and finite-data/base–fiber closure.

## Q7.3 What most raises the ceiling

Within the present simulation-only constraint, the single most valuable addition is **exact CFG closure across the four representation families plus the same-record Fog Bundle run**. That converts a compelling configuration into a parameter-free physical principle.

Beyond the binding sim-only constraint, the largest ceiling raiser by far would be one physical experiment showing both sides simultaneously:

- positive quotient-channel capacity through a real dynamic diffuser or MMF;
- conditional scene/fiber information at the statistical floor on the same raw record.

No additional simulation would raise the venue ceiling as much as that experiment.

---

# Frozen flagship contract

## Preferred title after CFG

> **Disorder as a Quotient Machine**

## One-sentence thesis

> **A symmetry-complete quadratic medium exposes the block-Gram quotient fixed by its self-adjoint commutant and hides the complementary unitary-orbit coordinates exactly; the same base supports calibration-free communication, while hidden motion on the fibers obeys the curvature and endpoint-toll laws.**

## Forbidden claims

- no claim that Schur’s lemma, power-spectrum phase loss, or maximal invariants are new;
- no claim that one fixed medium computes the full invariant algebra unless its kernel family spans the commutant;
- no claim that the original group orbit equals the hidden fiber—the correct object is the commutant gauge orbit;
- no equality claim `M_eff=13` without the dof audit;
- no universal higher-moment completeness;
- no “all invisible motion is recurrent” or universal mechanical toll;
- no resurrection of the general quadratic computer;
- no new TLSG quotient arms;
- no third-paper memory-effect line.

## Final lock

The architecture freezes when CFG-A/B/C and the NCX2 closure report. Thereafter:

- no additional protagonist search;
- no new main-text campaigns;
- no chronology in the manuscript;
- all remaining work is theorem compression, the four figures, source custody, prior-art fencing, and supplement assembly.

The second paper then has one protagonist and one paradox:

> **The fog can carry everything the symmetry permits—and absolutely nothing about the coordinates the symmetry erases.**
