# Beyond-Modulator-Band Direction — RISK GATES + P1 FISHER PROGNOSIS

**Date:** 2026-07-24 · **Status:** DEV-GRADE gates + prelaunch Fisher kill. No git commit.
**Scope:** the three fast-kill risk gates (G1/G2/G3) that must clear *before* any sealing, plus the
R39-ruling-mandated **P1 profiled-Fisher prognosis** (the prelaunch go/no-go). All ran locally on
the RTX 4060 (free; jobs ≪ 1 h), which is the coordinator's own condition for local execution; no
live Colab calls this session.

## Headline

> **P1 profiled-Fisher prognosis FAILS at the R39 frozen geometry → kill-tree node 2: STOP before
> reconstruction.** The beyond-band information is *not floored* (no exact null; NRMSE_CRB scales as
> `T_eff^(-1/2)`; profiling tax ≈ 0), but the **rate is too thin**: only **5% of the primary
> `1.0–1.8 k_p` annulus reaches Fisher SNR ≥ 3** at `T_eff=4096`, `n=1e4` (bar requires ≥ 70%). The
> usable Fisher-weighted aperture is **≈ 1.1× k_p, not the claimed 1.8×**. This is independently
> corroborated by the empirical G3 washout (viable band-extension depth ≈ 1 shell ≈ 1.1× k_p).

A prelaunch kill here saves the entire sealed campaign — the process working as designed. One
caveat: P1 is a **single analytic implementation** (eq 4.1); the ruling mandates a Monte-Carlo
covariance cross-check before the archive decision is final (§7.4). The signal is strong and
concordant with G3, so the recommendation is **STOP / archive the aperture+wall theorems as a
materials-bank note**, pending that cross-check.

---

## P1 — profiled-Fisher prognosis (`p1_fisher.py`, `P1_results.json`)

Frozen geometry (ruling §6.1): `N=4096`, `k_p=5` (real Fourier box, in-band dof 121, medium
db 120), `k_w=k_p`, claim `≤1.8 k_p`, support `≤2 k_p`, `M=128` signed band-limited codes,
`σ_f=0.30`, `n=1e4`. Beyond-band annulus (1.0–1.8 k_p) dof = 240. Exact profiled covariance Fisher
`J_B = I_ββ − I_βη I_ηη^† I_ηβ` (eq 4.1–4.2), nuisance η = in-band scene modes + mean channel +
law amplitude.

**Primary point (`T_eff=4096`, `n=1e4`, `σ_f=0.30`):**

| scene | exact null? | min-eig(J̄_B) | η_prof p10 / med | f_rec (SNR≥3) | NRMSE_CRB |
|---|---|---|---|---|---|
| witness (uniform annulus) | no | 7.4e-3 | 0.99 / 1.00 | **0.05** | 0.140 |
| cameraman | no | 6.1e-3 | 0.93 / 1.00 | 0.42 | 0.312 |
| coins | no | 1.2e-2 | 0.98 / 0.99 | 0.71 | 0.124 |
| moon | no | 8.3e-3 | 1.00 / 1.00 | 0.10 | 0.772 |

**T_eff sweep (witness):** NRMSE_CRB `0.560 → 0.396 → 0.280 → 0.198 → 0.140` for
`T_eff=256..4096` (clean `T_eff^(-1/2)`, no floor); f_rec `0.01 → 0.05` (grows far too slowly —
reaching 70% would need ~10²× more independent medium states).

**P1 bar (ruling §6.5):**

| check | threshold | result | pass |
|---|---|---|---|
| no exact profiled null | min-eig > 0 | 7.4e-3 | ✅ |
| η_prof p10 | ≥ 0.10 | 0.988 | ✅ |
| η_prof median | ≥ 0.25 | 0.996 | ✅ |
| NRMSE_CRB synthetic | ≤ 0.25 | 0.140 | ✅ |
| NRMSE_CRB natural median | ≤ 0.35 | 0.312 | ✅ |
| **f_rec (≥70% modes SNR≥3)** | **≥ 0.70** | **0.05** | ❌ |

**Verdict: P1_FAIL** — on the mode-coverage (f_rec) bar only. Interpretation: the support aperture
reaches 2·k_p and the *few* strong near-`k_p` modes carry enough energy that aggregate NRMSE_CRB
looks fine (0.14–0.31), but the **per-mode Fisher-weighted aperture collapses toward the edge**
(§2.3 `m(ξ)→0` near the rim), so the *claimed* `1.8×` annulus is mostly sub-SNR at feasible
integration. This is the "information too thin at publication scale" outcome the ruling named as the
immediate pre-launch kill.

---

## G1 — FINE-BAND MISMATCH (`gate_g1.py`, `G1_results.json`)  ·  verdict: SURVIVE-with-caveat

16×16 beyond-band cell (PB=3, medium 1≤i+j≤6), matched fixed+cov moment baseline = **0.516**
(reproduces R39's 0.47–0.79). Estimator fed a perturbed declared band; degradation vs matched:

| perturbation | beyond-band err | degradation |
|---|---|---|
| matched (band 1..6, flat) | 0.516 | — |
| (a) rotated 10% within fine annulus | 0.527 | **+2%** |
| (a) rotated 20% within fine annulus | 0.615 | **+19%** |
| (b) too narrow (declare i+j≤5) | 0.770 | +49% |
| (c) too wide (declare i+j≤7) | 0.687 | +33% |
| (d) wrong radial profile (true 1/f, declared flat) | **1.114** | **+116%** |

**Literal kill rule** (rotation 10–20% > 50%): **NOT fired** (+2% / +19%). But two loud caveats:
(b) band-*extent* knowledge is borderline load-bearing (+49%, at the 50% line); and **(d) radial-
profile mismatch collapses catastrophically (+116%, err > 1.0 = worse than nothing)**. Caveat on
(d): a steep true profile also physically *starves* the fine channel (less high-radial medium
power), so it conflates declared-law mismatch with signal weakening — but the practical consequence
stands: if the real medium's fine-speckle power rolls off radially (physically expected), beyond-
band recovery breaks. **Radial spectral slope is a load-bearing, mandatory mismatch axis** (ruling
P5 lists band/slope ±20%). This directly foreshadows the P1 f_rec thinness.

## G2 — GEOMETRY MIXING (`gate_g2.py`, `G2_results.json`)  ·  verdict: PASS

Multiplicative thin-screen vs partial convolutive propagation (small Gaussian PSF, σ=0.8 px).
Estimator keeps the pure multiplicative model; beyond-band err vs mixing α:

| α | blend `A=(1-α)(P·w)+α·blur(P·w)` | convolutive-medium variant |
|---|---|---|
| 0.00 | 0.516 (+0%) | 0.516 (+0%) |
| 0.10 | 0.514 (−0%) | 0.507 (−2%) |
| 0.25 | 0.545 (+6%) | 0.530 (+3%) |
| 0.50 | 0.583 (+13%) | 0.552 (+7%) |

**Kill rule** (collapse > 0.9 at α=0.1): **NOT fired** (0.514). Graceful degradation to +13% at
α=0.5. The multiplicative model is robust to partial convolutive mixing at the toy scale — the
declared home is not brittle to modest propagation. (Note: near the Fisher-thin aperture edge the
same mismatch can dominate per §5 eq 5.1 — G2 is only reassuring where the Fisher is healthy.)

## G3 — SCALE + NATURAL IMAGES, shell-resolved depth (`gate_g3.py`, `G3_results.json`)  ·  verdict: KILL

64×64 (N=4096), pattern band fx,fy≤8, medium 1≤i+j≤16, `M=128`, `T=4096`. Note: this gate tested a
**more aggressive** geometry (medium band far exceeding pattern) than the ruling's conservative
`k_w=k_p`; read it as the outer-envelope stress.

- **Part A — aperture law (in/out coverage separation): 18.76×** (in-coverage 0.250, out-coverage
  4.697). The Minkowski-support aperture law **holds strongly at scale** — the phenomenon is real.
- **Part B — band-extension DEPTH** (shell Δk = max(i,j)−PB, per-shell NMSE; toy moment depth = 1):

| spectrum / estimator | Δk=1 | Δk=2 | Δk=3 | Δk=4 | depth (≤0.3) | depth (≤0.5) |
|---|---|---|---|---|---|---|
| flat, moment (T=4096) | 0.468 | 0.835 | 0.812 | 1.145 | **0** | 1 |
| decaying (k⁻²), moment | 1.464 | 1.430 | 1.507 | 1.760 | **0** | **0** |
| flat, two-stage MLE (T=1024) | 0.567 | 1.266 | 1.009 | 3.758 | **0** | 0 |

Naturals (superres-zone NMSE): cameraman 1.918, coins 0.768 (Δk1 0.404), moon 15.7. The first
shell degrades vs the toy (0.318 → 0.468); **depth does NOT grow with scale/M_eff**. Under a
**realistic decaying medium spectrum the effect washes out entirely** (all shells > 1.0). The MLE
does not deepen. **Verdict KILL** for the beyond-band *method* at publication scale; the aperture
*law* survives (18.8×).

---

## Synthesis of the four probes

| probe | tests | outcome |
|---|---|---|
| aperture LAW (G3-A, P1 null/η_prof) | does the Minkowski-support aperture exist at scale | **YES** — 18.8× separation, positive profiled Fisher, no floor |
| beyond-band RATE (P1 f_rec, G3-B) | is the usable aperture as wide as claimed (1.8×) | **NO** — usable ≈ 1.1× k_p; f_rec 5%; depth ≈ 1 shell |
| robustness — geometry (G2) | multiplicative vs convolutive | PASS to α=0.5 (where Fisher is healthy) |
| robustness — spectral profile (G1(d), G3 decaying) | flat vs realistic decaying medium | **FRAGILE / washout** under radial rolloff |

The four are **mutually consistent**: the support aperture is real and unfloored, but the *Fisher-
weighted usable* aperture is a thin ≈1.1× rim that (i) does not grow with scale and (ii) collapses
under a realistic decaying medium spectrum. This is the exact `1.8×`-claim failure the ruling's P1
was designed to catch pre-launch.

## Recommendation

Per the ruling kill-tree (node 2), **P1_FAIL → STOP before reconstruction**; do not build/run the
lifted-GLS confirmatory (P1 gates it — the incomplete estimator is correct, not a gap). Archive the
**aperture-law + first-moment-impossibility theorems** as the materials-bank note (the already-
authorized law/boundary fallback). Before finalizing: run the ruling's mandated **MC-covariance
cross-check** of the Fisher (§7.4), and — if the operator wants to salvage a positive claim —
recompute P1 restricting the claim to `≤1.25 k_p` (the near shell the Fisher does support), which
may clear P1 as a *modest* (not 1.8×) extension.
