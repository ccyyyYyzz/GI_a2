# Probe B - does the O4-paired schedule's Fisher-side advantage survive to image level?

*Generated 2026-07-22 09:51:30Z | Python 3.11.5 | numpy 1.24.4 | runtime 23.9s | script `results/round63_next/SCHEDULE_PROBE/probe_b.py`*

ROUND63-NEXT software-method hunt, Probe B. Outcome-blind in design; evaluated only at image level (PSNR of actual reconstructions). READ-ONLY on all inputs.

## Verdict

At the frozen OU correlation time (tc=2 frames), pooled over 4 scenes x 3 drift magnitudes (12 cells, 5 paired seeds each):

| image-level PSNR gain | blind LS | drift-aware joint |
|---|---|---|
| **paired - ordered** | +1.34 dB (range +0.18..+2.01) | +1.35 dB (range +0.24..+2.02) |
| **paired - random**  | +0.10 dB (range -0.02..+0.24) | +0.10 dB (range -0.07..+0.25) |
| (random - ordered)   | +1.24 dB (range +0.14..+1.95) | +1.25 dB (range +0.18..+1.93) |

**The Fisher-side advantage splits in two at image level.** The naive-ordering penalty *survives*: the O4-paired schedule beats naive bank-order by +1.34 dB pooled (up to +2.0 dB at severe drift, +2.6 dB in the slow-drift sweep). Here the real bank order is only mildly adversarial - its Schur-loss elevation over random is ~2.0x on these patterns, not scatter_verify's synthetic 27.8x (which used a deliberately monotonic carrier ordering) - yet it still costs >1 dB. But the paired-over-**random** advantage - the headline 63.6 -> 0.0 three-orders-of-magnitude Fisher win of the exact moment-zeroing construction - **does not survive**: it is +0.10 dB (blind) / +0.10 dB (aware) pooled at the frozen correlation, i.e. within seed noise, and in the slow-drift sweep the aware estimator's paired-vs-random gap goes slightly *negative* (the estimator already removes the low-frequency drift the paired schedule was built to cancel, so the exact-zero moment buys nothing extra, while random keeps the measurement diversity the paired schedule spends on literal pattern duplication). Almost all of the recoverable image-level benefit is captured by *any* interleaving; the specific paired construction adds nothing over a random permutation. This is the RLMI-bridge pattern: a large Fisher-side margin between two already-good designs vanishes at image level.

## Frozen model and parameters

Lifted verbatim from `code/scatter/verify_o1_o4.py` (the committed Fisher-side verification, results/scatter_verify/):

- **Drift**: multiplicative gain `a_n = exp(l_n)`, stationary OU log-gain, exact AR(1) `l_n = mu + phi(l_(n-1)-mu) + sigma_l*sqrt(1-phi^2)*eps_n` (the process the paper's Tauchen grid discretizes). Frozen mean `mu = -0.5 sigma_l^2` (E[a]=1 exactly). Frozen correlation time `tc = 2.0` frames (lag-1 corr `exp(-1/tc) = 0.607`). Drift magnitude swept mild/moderate/severe `sigma_l in [0.05, 0.15, 0.4]` (the O1 verification used 0.30).
- **Measurement**: `Y_n ~ Poisson(a_n * s_n)`, `s_n = m_n . x` (dose-scaled), exactly the paper's conditionally-Poisson bucket model. NO dead time / jitter (chapter separation - drift only).
- **Patterns**: 972 SCAT32 rows of the deployed bank `results/round63_bridge/library/L0.npz` (sha256 `cf04a6cda34b0d30...`), each duplicated into a matched pair -> **N = 1944 frames**. Identical multiset for all arms.
- **Scenes**: 4 frozen 32x32 bridge scenes ['bridge_contour_1', 'bridge_microtex_1', 'bridge_twopop_1', 'bridge_control_2'] (values in [0,1]).
- **Dose**: pattern scale x4 (mean bucket count ~2200/frame), so shot noise (~2%) is subdominant to drift (5-40%) by design - this isolates the pattern-ORDER effect, which is the probe's question. Lower dose would add shot noise equally across arms and cannot change the ordering conclusion.
- **Reconstruction (identical for all arms)**: (i) drift-**blind** ridge LS `x_hat = (A^T A + lam I)^-1 A^T Y`; (ii) drift-**aware** joint estimator - the exact O4.1/O4.2 model `y = M x + D_s H beta`, profiling out the nuisance beta by Schur complement (Woodbury), 3 relinearizations, declared low-frequency modes H = orthonormal Legendre orders [1, 2, 3, 4, 5, 6] (DC excluded = gauge). Ridge `lam = 0.1 * mean(diag(A^T A))`, R = I, identical for every arm/scene/drift. One Cholesky reused throughout.
- **Schedules** (order only; identical patterns, identical design exposure, 5 paired seeds sharing one OU path per seed): **ordered** = bank order at increasing time; **random** = frozen uniform permutation per seed; **paired** = the mirror-chop that zeroed the moment in scatter_verify (pair j's two copies to slots half-1-j and half+j, i.e. symmetric +/-t about centre).

## Moment-condition metric per schedule (confirms the construction matches scatter_verify)

Schur-complement loss `tr(I_xx - I_x|gain)` and cross-block max `|M^T R^-1 D_s H|_max`, modes H=[t, t^3], R=I - exactly the scatter_verify block2(iii) metric, here evaluated on each real scene x. scatter_verify reference (synthetic): ordered 1769.9 >> random 63.6 > paired 0.0.

| scene | ordered loss | random loss | paired loss | ordered mom-max | paired mom-max |
|---|---|---|---|---|---|
| bridge_contour_1 | 1.464e+05 | 7.225e+04 | 0.00e+00 | 4.19e+05 | 2.47e-10 |
| bridge_microtex_1 | 1.440e+05 | 7.286e+04 | 0.00e+00 | 1.49e+05 | 1.02e-10 |
| bridge_twopop_1 | 1.441e+05 | 7.492e+04 | 0.00e+00 | 2.26e+05 | 1.53e-10 |
| bridge_control_2 | 1.440e+05 | 7.222e+04 | 0.00e+00 | 1.88e+05 | 1.02e-10 |

Paired zeroes the moment to machine precision on every real scene (loss = 0.0, mom-max = 0.0), and ordered > random > paired holds throughout - the image-level pattern set reproduces the frozen Fisher-side ordering. Note the real bank-order ordered/random Schur-loss ratio is only ~2.0x (vs scatter_verify's synthetic 27.8x): the deployed SCAT32 bank order is far less adversarial than the monotonic-carrier ordering used in the synthetic demo, so the image-level ordered penalty measured below is a conservative (mild-ordering) estimate. (Absolute magnitudes differ from the synthetic reference because s_n and the dose scale differ; the ordering and the exact paired-zero are what transfer.)

Noiseless reconstruction ceilings (this operator, no drift): contour_1 14.6 dB, microtex_1 22.7 dB, twopop_1 18.2 dB, control_2 16.1 dB.

## Primary table - PSNR (dB), frozen tc=2, 5 paired seeds

Each cell: mean PSNR over 5 seeds. Deltas are seed-paired (per-seed difference, then mean +/- sd over seeds).

### BLIND estimator

| scene | sigma | ordered | random | paired | paired-ordered | paired-random |
|---|---|---|---|---|---|---|
| contour_1 | 0.05 | 12.38 | 12.92 | 13.07 | +0.70+/-0.23 | +0.15+/-0.09 |
| contour_1 | 0.15 | 6.45 | 8.02 | 8.16 | +1.72+/-0.36 | +0.14+/-0.13 |
| contour_1 | 0.4 | -1.63 | 0.33 | 0.39 | +2.01+/-0.45 | +0.06+/-0.31 |
| microtex_1 | 0.05 | 20.12 | 20.72 | 20.70 | +0.58+/-0.21 | -0.02+/-0.15 |
| microtex_1 | 0.15 | 14.16 | 15.77 | 15.78 | +1.62+/-0.35 | +0.01+/-0.36 |
| microtex_1 | 0.4 | 6.09 | 8.05 | 8.06 | +1.97+/-0.30 | +0.02+/-0.37 |
| twopop_1 | 0.05 | 16.35 | 16.75 | 16.96 | +0.62+/-0.18 | +0.21+/-0.18 |
| twopop_1 | 0.15 | 11.19 | 12.64 | 12.88 | +1.69+/-0.36 | +0.24+/-0.26 |
| twopop_1 | 0.4 | 3.40 | 5.28 | 5.40 | +2.00+/-0.47 | +0.12+/-0.32 |
| control_2 | 0.05 | 15.27 | 15.41 | 15.45 | +0.18+/-0.12 | +0.04+/-0.11 |
| control_2 | 0.15 | 11.98 | 12.92 | 13.09 | +1.11+/-0.30 | +0.17+/-0.29 |
| control_2 | 0.4 | 5.07 | 6.88 | 6.93 | +1.86+/-0.45 | +0.05+/-0.28 |

### AWARE estimator

| scene | sigma | ordered | random | paired | paired-ordered | paired-random |
|---|---|---|---|---|---|---|
| contour_1 | 0.05 | 12.32 | 12.92 | 13.05 | +0.74+/-0.24 | +0.14+/-0.08 |
| contour_1 | 0.15 | 6.46 | 8.04 | 8.18 | +1.72+/-0.37 | +0.15+/-0.16 |
| contour_1 | 0.4 | -1.59 | 0.34 | 0.42 | +2.02+/-0.44 | +0.09+/-0.33 |
| microtex_1 | 0.05 | 20.10 | 20.73 | 20.66 | +0.57+/-0.20 | -0.07+/-0.14 |
| microtex_1 | 0.15 | 14.19 | 15.80 | 15.79 | +1.60+/-0.32 | -0.01+/-0.36 |
| microtex_1 | 0.4 | 6.14 | 8.06 | 8.09 | +1.95+/-0.25 | +0.03+/-0.36 |
| twopop_1 | 0.05 | 16.34 | 16.76 | 16.97 | +0.64+/-0.18 | +0.22+/-0.19 |
| twopop_1 | 0.15 | 11.21 | 12.66 | 12.91 | +1.69+/-0.33 | +0.25+/-0.28 |
| twopop_1 | 0.4 | 3.44 | 5.30 | 5.44 | +2.00+/-0.42 | +0.14+/-0.33 |
| control_2 | 0.05 | 15.23 | 15.40 | 15.46 | +0.24+/-0.12 | +0.06+/-0.11 |
| control_2 | 0.15 | 11.97 | 12.92 | 13.11 | +1.13+/-0.31 | +0.19+/-0.30 |
| control_2 | 0.4 | 5.09 | 6.90 | 6.96 | +1.87+/-0.43 | +0.06+/-0.28 |

## Secondary - correlation-time sensitivity (bridge_contour_1, sigma=0.15)

The frozen tc=2 is *fast* relative to N=1944 frames, so the OU drift is dominated by high-frequency residual r (per the paper's model l = H beta + r) that no pattern order can cancel. This sweep varies tc to show the boundary: as the drift slows (more of it lives in the low-order modes H), the aware estimator recovers more and the ordered penalty grows - but the paired-vs-random gap stays sub-dB and sign-unstable throughout.

| tc (frames) | ord blind | rnd blind | pair blind | pair-rnd blind | ord aware | rnd aware | pair aware | pair-rnd aware |
|---|---|---|---|---|---|---|---|---|
| 2 | 6.45 | 8.02 | 8.16 | +0.14+/-0.13 | 6.46 | 8.04 | 8.18 | +0.15+/-0.16 |
| 10 | 6.11 | 8.34 | 8.37 | +0.03+/-0.26 | 6.31 | 8.48 | 8.54 | +0.06+/-0.25 |
| 50 | 6.25 | 8.62 | 8.62 | -0.00+/-0.95 | 7.20 | 9.37 | 9.20 | -0.17+/-0.64 |
| 200 | 7.24 | 9.29 | 9.87 | +0.59+/-0.77 | 9.51 | 11.16 | 11.04 | -0.12+/-0.48 |
| 1000 | 9.74 | 11.12 | 11.29 | +0.17+/-1.13 | 12.17 | 12.76 | 12.68 | -0.09+/-0.18 |

## What materially drives the conclusion (honesty flags)

- **Correlation time vs acquisition length.** The single biggest lever. The O4 schedule cancels only the *low-order* drift moments (Thm O4-A, eq. O4.4; H are 'low-frequency or OU/KL drift modes', main_scatter.tex L416). For fast drift (frozen tc=2 << N) most of the OU power is high-frequency residual that no order cancels, so scheduling is nearly irrelevant; even in the slow-drift limit the paired-vs-random gap never exceeds ~0.6 dB. The ordered-vs-interleaved penalty is the only robust effect and it too shrinks toward the very-slow limit.
- **Estimator.** Under blind LS the paired schedule's odd-moment cancellation shows as a small (~+0.1 dB) edge over random; under the drift-aware estimator that edge disappears or reverses, because the estimator already removes the modeled low-frequency drift, and the paired schedule's literal pattern duplication costs a little measurement diversity that random keeps. Neither estimator turns the Fisher-side paired-over-random margin into an image-level win.
- **Regularization / dose.** Ridge lam and dose set the absolute PSNR and the drift-vs-shot balance but are held identical across arms; they move all three arms together and do not create or erase the (non-)ordering effect (verified across lam_scale in {0.03,0.1} and the tc sweep).
- **Matched-pair duplication** halves distinct patterns (rank <= 972 for 1024 pixels); this is the scatter_verify construction and is identical for all arms, so it caps absolute PSNR but not the comparison.

## One-paragraph verdict

Image-level dB gain of O4-paired scheduling, pooled over 4 scenes x 3 drift magnitudes at the frozen OU correlation (tc=2): **over naive ordered, +1.34 dB blind / +1.35 dB aware** (real, and the image-level survivor of the ordered-over-random Schur-loss elevation - ~2.0x on these real patterns, milder than scatter_verify's synthetic 27.8x); **over random, +0.10 dB blind / +0.10 dB aware** - i.e. essentially zero (within seed noise), and drifting slightly negative for the aware estimator as the drift slows (secondary sweep). The paired schedule's exact moment-zeroing (Fisher-side loss 63.6 -> 0.0 vs random, verified to machine precision here too) therefore **does not survive to image level**: random interleaving captures all the recoverable benefit, and the three-orders-of-magnitude Fisher advantage of the exact O4 construction over a plain shuffle is an information-geometry artifact with no reconstruction payoff in this regime. The honest, deployable takeaway is narrower and still useful: **do not acquire in bank/raster order - interleave (randomly is enough)**; the extra machinery of the exact paired chop is not worth it at image level.

## R29 DEV gate (formal)

R29 ruling (docs/ROUND63_GPT_ROUND29_RULING_RAW.md, sec 5.1) designates this probe the official pre-campaign DEV gate for the top-ranked method lane **DOPS-GI (drift-orthogonal pattern scheduling)**. The gate reads the **drift-BLIND** reconstruction only - no estimator rescue permitted. Evaluated on the frozen 6-scene DEV cohort at the hardest drift cell (sigma=0.4, tc=2.0), 5 paired seeds. Artifact: `DEV_GATE_VERDICT.json` (script `dev_gate.py`).

Per scene: `Q_base = max(Q_ordered, Q_random)` (5-seed mean, blind); `dQ = Q_paired - Q_base`.

| scene | Q_ordered | Q_random | Q_paired | Q_base | dQ (dB) |
|---|---|---|---|---|---|
| contour_1 | -1.63 | 0.33 | 0.39 | 0.33 | +0.060 |
| microtex_1 | 6.09 | 8.05 | 8.06 | 8.05 | +0.018 |
| twopop_1 | 3.40 | 5.28 | 5.40 | 5.28 | +0.122 |
| control_2 | 5.07 | 6.88 | 6.93 | 6.88 | +0.048 |
| contour_3 | -1.65 | 0.30 | 0.33 | 0.30 | +0.030 |
| microtex_3 | 2.51 | 4.54 | 4.48 | 4.54 | -0.052 |

Gain-free control (sigma=0), loss = `Q_base - Q_paired` per scene:

| scene | contour_1 | microtex_1 | twopop_1 | control_2 | contour_3 | microtex_3 |
|---|---|---|---|---|---|---|
| loss (dB) | +0.091 | -0.002 | +0.023 | +0.033 | +0.025 | -0.005 |

**Gate conditions (all three required for PASS):**

- (a) median dQ >= 1.0 dB: median = **0.039 dB** -> **FAIL**
- (b) >= 5/6 scenes dQ > 0: **5/6** positive -> **PASS**
- (c) control (sigma=0) worst loss <= 0.5 dB: worst = **0.091 dB** -> **PASS**

### VERDICT: DEV_GATE_FAIL

The exact O4-paired construction does not beat the better of the two simple baselines (ordered, random) by the required 1.0 dB median - it beats it by **0.04 dB** median (max +0.12 dB, min -0.05 dB). Per R29 sec 5.1 a fail means **no DOPS-GI campaign is launched; the reserve lane (CPL-GI micro-probe) activates instead.** Conditions (b) and (c) pass (paired is marginally >= the best baseline on 5/6 scenes and never hurts under zero drift), but the effect size is two orders of magnitude below the gate threshold. This is not softened: the drift-orthogonal *ordering* delivers no material image-level advantage over a plain random interleave in this regime.

**Two descriptive by-products (honest, not rescues):**

1. **Anti-naive-ordering (prior-art-adjacent).** paired - ordered = **+1.34 dB** pooled (primary grid), **+1.96 dB** at the hard cell (6 scenes). Avoiding bank/raster acquisition order is worth ~1-2 dB, but this is the interleaving effect (prior art), captured equally by a random permutation - it does not license the DOPS-specific paired construction.
2. **Correlation time is the dominant lever.** At the frozen tc=2 (<< N=1944) most OU power is uncancellable high-frequency residual r (paper model l = H beta + r); the low-order moment cancellation the paired schedule performs therefore has little image-level purchase. The paired-vs-best-baseline gap never exceeds ~0.6 dB anywhere in the tc sweep (secondary table).

