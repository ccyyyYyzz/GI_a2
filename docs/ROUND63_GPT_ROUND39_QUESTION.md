# R39 — Beyond-modulator-band imaging via medium-covariance: theorems, sealed probe, and ruling

Housekeeping. R38's CONDITIONAL GO probe was overtaken by events, all consistent
with its own analysis: (i) the pathwise blind route was killed empirically at its
own kill tree — E4/E6 showed refinement seeded at the TRUE medium drifts to
wrong-null solutions (0.67–0.80) with tight multi-start agreement on the SAME
wrong solution (std 0.006), i.e. §1.5 global non-identifiability made concrete;
chronology randomization (the operator's own GI_a1 d=1 remedy) does NOT lift to
d=48 even with a verified-whitened carrier; (ii) the covariance/marginal route
your eq (1.18) pointed to was built and RESURRECTED the direction: injective
(0.054 vs exact target), collusion structurally impossible (multi-start std
0.012–0.017), consistent (bar-crossing at T≈2048 moment-matching), photon-robust
(lags shot-free; 1e4≈clean), mismatch-robust (F5-style full table: 20% basis
rotation +3%, tau 2x +22%, sigma_f ±33% ≤+7%, dim ±8 ≤+7% — all under the 25%
kill line); marginal-likelihood MLE (Kalman, z integrated out exactly) reduces
the temporal price ~4x (bar-crossing T≤512 at 32×32; T≤128 at 16×16);
(iii) the fresh-pattern comparator (your §2.5) was run EARLY and honestly:
in-band, fresh known patterns dominate everything (8x scene, 9/9 cells across
sigma_f∈{0.3,0.6,1.0}×photons{1e3..1e5}, and they read the medium law too) —
the universal software-only flagship claim is dead per kill-tree node 7;
(iv) the direction then pivoted to the configuration that is F4-IMMUNE BY
PHYSICS, which is what R39 is about. Also: third-order cumulants FAIL to buy
patterns (M=12/20 no unlock; info decays with moment order — the MOLT death
law re-measured); structured Fourier-lattice patterns collapse the covariance
estimator (three independent tests) — incoherent random codes are load-bearing.

## The pivot: beyond the modulator's band

The programmable modulator's spatial-frequency band is a hardware wall (DMD
pixel pitch). A scattering medium's speckle grain can be an order of magnitude
finer. In the declared geometry (thin dynamic scattering screen at/near the
object plane, multiplicative: b = <P ⊙ w_t, x>), the medium's spectral band
extends the reachable scene spectrum per the empirically confirmed APERTURE LAW:

    recoverable spectrum = pattern band ⊕ medium spectral band
    (verified: in-coverage rel. err 0.073 vs out-coverage 1.213 — 17× separation;
     band-limited incoherent random patterns, matched law, T=4096, 16×16)

Beyond-band duel at equal budget (16×16, pattern band fx,fy≤3, medium band
i+j≤6, deliberate beyond-band scene content, sigma_f=0.3, 1e5 photons/bucket):

| arm | beyond-band rel. err |
|---|---|
| ANY patterns + mean route | 1.000 exactly (band-limited operators cannot leave their band) |
| FRESH patterns + covariance route (equal budget, first-pass estimator) | 0.955 / 0.957 |
| FIXED bank + covariance route (moment matching) | 0.650 / 0.547 / 0.480 @ T=512/1024/2048 |
| FIXED bank + MLE | 0.505 / 0.439 / 0.481 |
| ORACLE (medium known) | 0.389 / 0.221 / **0.091** |

E8 observation (central new fact): the beyond-band ORACLE improves strongly
with T (0.09 at T=2048 — the information is rich), but BLIND extraction
STALLS at ~0.44–0.48 and stops tracking the oracle (in-band, MLE tracked the
oracle closely; beyond-band it does not). The blind/oracle gap (~5× in energy)
is the open question: a fundamental profiling tax specific to beyond-band
modes (scene fine modes and medium fine modes cohabit the same sideband
covariance entries), or estimator immaturity (nonconvex plateau, gauge)?
Ask 3 below is upgraded to demand this accounting explicitly.

Reading: (a) a CHANNEL-level impossibility — beyond-band content is carried
exclusively by the medium-covariance channel; the first-moment channel is
mathematically blind; (b) a DESIGN-level result — within the covariance
channel, repeated/concentrated codes decisively beat fresh codes at equal
budget (the concentration principle: repeated pairs average the same
covariance entries; fresh spreads budget over single-sample products). Caveat
disclosed: the fresh+cov arm used a first-pass estimator; the sealed probe
must give it best-effort form (GLS weighting) before the design half is frozen.

## Prior-art fence (hostile scout, completed)

Verdict PARTIALLY_OCCUPIED, not killed. Occupied and to be ceded head-on:
unknown-speckle covariance statistics beating a resolution band (blind-SIM
Mudry 2012; PE-SIMS Yeh 2017; Chaigne–Bossy Optica 2017 — the mechanistic
twin, photoacoustic, transducer ARRAY; SOFI 2009); medium injecting
beyond-aperture frequencies (scattering superlens Vellekoop 2010/Choi 2011;
SAI Leonetti 2019); medium-statistics-as-resolver HBT framing (Chandra 2023;
Bertolotti 2012/Katz 2014). OPEN conjunction (the defensible square): true 0-D
bucket (beyond-band content lives ONLY in the scalar temporal covariance) ×
separate band-limited modulator as the exceeded reference band × statistics-only
declared law (no TM, no reference arm, no wavefront sensor, no speckle
calibration) × the provable impossibility separation for the mean channel ×
the aperture-law theorem. The scout's frozen claim sentence is on file; the
impossibility separation is the single strongest novelty lever (no precedent
proves its baseline CANNOT reach the recovered content).

## Interpretation frame (for the manuscript's spine)

The medium is the unmeasured second arm of an intensity interferometer: the
bucket autocorrelation is the fringe; only the arm's coherence STATISTICS are
needed (HBT's own logic — correlation is immune to unknowable realizations);
the medium's spatial band is a synthesized aperture bolted onto the pattern
bank; integration time (independent medium realizations, one per t_c) is the
currency; detector dead-time physics caps the correlator bandwidth (the M1
operating map becomes the ceiling of the design window). Why realizations are
unknowable but the law suffices: expectation over realizations converts
unknown products into known-law quadratic forms — the unknowns never enter the
solved equations; the law is stable, measurable from the same record, and E5d
shows it is only needed roughly.

## R39 asks

1. **Aperture-law theorem.** Exact statement and conditions for: recoverable
   spectrum = pattern band ⊕ medium spectral band, for the bucket-covariance
   channel — genericity/coverage-multiplicity conditions, the gauge (static
   smooth envelope), and the Fisher-weighted refinement (per-frequency
   information density across the coverage set, not bare support).
2. **Channel-impossibility theorem.** Formal version of the 1.000 wall: any
   estimator using only first moments of band-limited known patterns (fresh or
   reused, any schedule) has zero information on beyond-band content; state
   exactly which channels (covariance and higher) carry it, and the
   information ordering (we measured cumulant-order decay — floors ≥3 starved).
3. **Sample-complexity law AND the blind/oracle gap.** T_eff(N, M, band
   geometry, d_w, sigma_f) for the MLE; the beyond-band recoverable fraction
   vs budget; where the concentration principle enters (repeated-pair
   averaging vs fresh single-sample products). CRITICALLY: explain the E8
   stall — oracle reaches 0.091 while blind saturates at ~0.44–0.48. Compute
   the profiled Fisher for beyond-band modes specifically (scene fine modes
   vs medium fine modes cohabiting sideband covariance entries): is there a
   fundamental beyond-band profiling floor, and at what value? If the floor
   is ~0.4, the quantitative claim must be scoped to it; if there is no
   floor, name the estimator upgrade that closes the gap. This number decides
   the strength tier of the entire manuscript.
4. **Geometry validity and mismatch.** The multiplicative thin-screen-near-
   object model is the declared home (shower-glass regime). State the honest
   validity boundary (what fraction of convolutive mixing breaks it) and the
   F5-grade mismatch tests the sealed probe must include (partial convolution,
   wrong band, wrong law class).
5. **Sealed-probe spec (the exam).** Freeze it in your R34/R38 style, now at
   PUBLICATION scale: N=64×64 minimum, natural images + synthetic beyond-band
   witnesses, arms = {fixed+MLE (production), fixed+moment, fresh+mean (the
   wall), fresh+cov BEST-EFFORT (GLS), oracle ceiling, classic averaging},
   sealed banks, frozen bars (coverage of the aperture-law prediction map,
   beyond-band NMSE thresholds, mismatch degradation ≤25%, multi-start
   agreement, equal-budget accounting), kill tree, no repair round. Include
   the band-extension factor to claim (2x linear?) and photon/realization
   budgets stated in physical units for the declared geometry.
6. **Sim-only referee defense.** The manuscript must be publishable and
   credible WITHOUT any experiment (operator's binding constraint; experiments
   and collaboration come only after completion). Architect the defense: what
   the preregistration/disclosure apparatus must show, the Chaigne–Bossy
   head-on paragraph, the forbidden-claims list, and the experimental-transfer
   section written as a precise falsifiable proposal (rotating fine diffuser
   near object; DMD band-limited; bucket) that we do NOT depend on.
7. **Blunt judgment** against the operator's bar (surprising inversion,
   made-possible-at-all effect, robust, universal-within-declared-domain,
   software-only, honest): eye-lighter or boring? If the Fisher accounting in
   (3) predicts the beyond-band information is too thin at publication-scale
   N, say KILL now.
8. **Manuscript consequence.** If the probe passes: is this THE flagship (one
   protagonist: the aperture law + beyond-band method, with M1 as the
   design-window ceiling and the killed universal claim disclosed as the
   in-band F4 map), with the program's prior negative arcs as supplement? Or
   a standalone letter plus the map paper? One-paper discipline binds; state
   the architecture seed (≤7pp, ≤4 figs, DEFINE→LAW→WALL→METHOD?).

Deliver as a GitHub issue titled R39. The machinery is warm; every number
above is reproducible from the repo's FOG_DMD_PROBE / FOG_DMD_E6 directories
and the session scratchpad scripts.
