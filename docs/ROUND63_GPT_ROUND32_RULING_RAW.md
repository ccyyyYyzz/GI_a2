# R32 ruling (GitHub issue #24, raw)

Title: R32
Posted: 2026-07-22T12:57:18Z

---

# R32 — Fast strategy ruling: MOLT survives only as a time-limited, dose-expensive sparse-support method

**Reference checkpoint:** commit [`5ed16d2`](https://github.com/ccyyyYyzz/GI_a2/commit/5ed16d2), especially `results/round63_next/SATURATION_JAILBREAK_PROBE/`.

## One decision

> **CONDITIONAL GO. Adopt equal wall-clock time as the primary resource frame for a narrowly declared, dose-insensitive and switching-limited MOLT campaign. Keep matched-photon performance as a mandatory adverse audit. MOLT is dead for photon-starved or dose-limited applications, and the dense-mask `p2` arm is retired now.**

Equal-photon does **not** universally govern every imaging experiment. If the sample is static and dose-insensitive, the source can supply the required flux within its declared peak/average-power limits, and detector gates fit inside the mask dwell that already dominates acquisition time, then photons are the traded resource and time is the scarce resource. That is a legitimate operating regime, not accounting fraud.

But the trade must be stated exactly: MOLT is **not photon-efficient** and the second operator is not free. The current checkpoint shows that dense 50% masks fail the matched-photon bar by orders of magnitude: median `d' = 0.00563`, about `2.84e5 x` budget for median `d'=3`; 10% `p2` maps require roughly `5e5–6.85e5 x`; and even the favorable dense null pair needs `9,311 x` under optimal allocation. Those results permanently kill dense-mask MOLT economics.

The same checkpoint also establishes that the mathematical object is real: the log-curvature estimator is accurate at about `5e-7` relative error noiselessly; the finite-`C` covariance and saturation ridge are correct; elementary-symmetric subtraction is unusable; the joint Jacobian doubles rank `51 -> 102`; and on known support `s=80>K=51`, the linear reconstruction gives `27.8 dB` while the joint operator gives `296.7 dB` noiselessly. Therefore the right strategic response is not to kill the theorem, but to move to the sparse-mask regime where the dose law is favorable.

## Resource-frame ruling

The primary comparison is **equal end-to-end wall-clock time**, not merely equal detector-gate duration. A valid equal-time experiment must match all of the following:

- total acquisition time, including mask updates, settling, power changes and readout;
- total number of pattern slots, unless an explicitly measured switching slack absorbs the extra sweep;
- the same source peak-power and average-power ceilings;
- the same detector recovery, crosstalk and no-fire guards;
- the same hardware and scene.

The linear comparator must be allowed to use the **same total time and available photon flux**. It may spend all of those photons reducing linear-bucket noise. MOLT earns a win only when the additional nonlinear operator closes a prior/support fiber that no amount of linear averaging can remove.

Every result must also carry a photon ledger: incident photons, detected primary-photoelectron opportunities, peak fluence, average power and photon multiple relative to the original linear campaign. The matched-photon result remains a required adverse panel and must be described as a failure. Equal-time is the primary frame only because the declared application is time-limited and dose-insensitive; it does not erase the dose cost.

**Immediate resource kill:** if the proposed sparse sweep exceeds the declared sample/source dose envelope, requires a higher peak or average power than the equal-time comparator is permitted, or increases measured end-to-end time by more than 10%, the equal-time framing fails and MOLT dies for this campaign.

## Frozen honest headline

Use this form:

> **Sparse MOLT trades photon dose for sensing rank: a microcell-occupancy sweep on support-targeted sparse masks supplies a mask-indexed quadratic operator in addition to ordinary linear buckets. At an expected dose floor of order `30 x` the original linear campaign for 10%-precision `p2` when `n_eff <= 32`, the repeated gates can fit within one mask dwell in a MHz-gated, switching-limited system, allowing known or certified low-dimensional support fibers to close when linear sensing alone remains nonidentifiable.**

Do **not** freeze the shorter statement “jailbreak at `2K>=s`” without its conditions. For a common mask bank `M` measuring both `Mx` and `M(x^2)`, the generic local/known-support threshold is indeed

```text
rank([M_S ; 2 M_S diag(x_S)]) = s,
```

with the dimension-counting threshold `2K >= s`. This is a local or known-support statement, not a global unknown-support recovery theorem.

For the actual proposed hybrid with a dense linear bank `M_D` and a sparse swept bank `M_S`, the correct certificate is

```text
rank([
  M_D|_S ;
  M_S|_S ;
  2 M_S|_S diag(x_S)
]) = s.                                                     (R32.1)
```

The sparse sweep also supplies its own linear `p1` equation; do not discard it. The nominal generic row budget is therefore at most `K_D + 2 K_S`, subject to coverage and conditioning. For a general `d`-dimensional prior tangent `B_x`, replace the support restriction by multiplication with `B_x` and require full rank `d`.

The word **jailbreak** is authorized only when this restricted Jacobian is numerically certified to close the frozen support/prior fiber. It is forbidden as a dense-scene or prior-free claim.

## Minimal campaign skeleton

### 1. Frozen scope

- static, nonnegative, dose-insensitive scenes;
- either a known support `S`, or a predeclared low-dimensional model with a computable tangent/active subspace;
- no claim for biological, photobleaching, eye-safety-limited, astronomy, remote low-light, heating-limited or otherwise photon-starved imaging;
- `p1` and `p2` only; `p3+` remains killed.

### 2. Acquisition

**Dense linear bank.** Acquire ordinary low-occupancy buckets with a dense bank `M_D`. This provides the robust linear image/support channel.

**Sparse MOLT bank.** Use binary masks with row weight at most 32, preferably targeted to the known or linear-estimated support. Each displayed sparse mask receives:

- a low-power anchor for `p1`, power scale and dark response;
- an intermediate audit level near `t=1` during model validation;
- a curvature level in the `t≈2.8–3.5` ridge, with the exact allocation nuisance-profiled rather than copied from the known-`p1` ideal.

Fit log-survival/power-sum coordinates. Do not use `e2` subtraction. Retain both the sparse-bank `p1` and `p2` equations.

**Equal-time comparator.** Give the linear baseline the same complete pattern/time schedule, photon flux and hardware limits. Replace the sparse swept slots by additional linear masks or the strongest predeclared linear design; do not compare against an artificially photon-starved linear arm.

### 3. Reconstruction and claim target

Use the joint count likelihood, or a validated Gaussian approximation with the full cross-level covariance. The primary endpoint is support/prior-restricted recovery and null-pair separation, not arbitrary dense-image PSNR.

The production certificate must report:

- rank and smallest Fisher-weighted singular value of (R32.1), or its prior-tangent analogue;
- support coverage of the sparse bank;
- nuisance-profiled `p2` uncertainty;
- adversarial residual-fiber search;
- reconstruction against the equal-time linear comparator.

### 4. Binding go/no-go bars for the finishing sparse run

Proceed to one confirmatory campaign only if all hold:

1. `n_eff <= 32` on the actual sparse masks and target scenes, with measured `C_eff` and effective mode number.
2. A 10% `p2` map is achieved at no more than `40 x` the original linear detected-photon budget per mask; the expected landing point is about `30 x`.
3. A favorable order-one `Delta p2/p2` null-pair discrimination reaches `d' >= 3` at no more than `4 x` the original linear dose; the expected landing point is about `2.6 x`.
4. The entire sweep fits the measured system schedule with at most 10% wall-clock overhead, including mask changes and settling—not merely microcell gates.
5. The physical P0/P1 conditions from R31 survive: diagonal/calibrated incoherent mixing, stable `C_eff`, no unmodelled refiring/over-saturation, and held-out sweep prediction.
6. The frozen known-support or prior-tangent Jacobian is full rank and adequately conditioned, and the noisy joint reconstruction materially beats the strongest equal-time linear comparator.

Failure of any one bar kills MOLT as the flagship. The conditional Laplace theorem remains a materials-bank result.

## Claim fence

### Authorized

- MOLT creates a second sensing operator by deliberately using microcell occupancy saturation.
- Sparse support-targeted masks improve the economics through the `n_eff^2` dose law.
- The joint operator can double local rank for a common mask bank and can close known or certified low-dimensional prior fibers.
- In switching-limited, dose-insensitive systems, the photon-expensive sweep may add little wall-clock time.
- The method is a **photon-for-rank trade**.

### Forbidden

- universal or dense-scene null-space jailbreak;
- exact recovery at 5% sampling without a support/prior certificate;
- photon efficiency, dose neutrality, zero-cost sensing or zero calibration;
- “negligible time” based only on MHz gate arithmetic without end-to-end timing;
- hiding the matched-photon failure;
- extrapolating the `s=80` noiseless witness to unknown-support noisy recovery.

## Final ruling

**MOLT does not die today. Dense-mask MOLT dies today.** Continue only with the sparse, support/prior-closed, equal-time campaign above. Its publishable thesis is not that saturation beats linear GI for the same photons; it is that, when photons are cheap but pattern time and hardware channels are scarce, a detector defect can be deliberately purchased with dose to create an additional sensing operator and close a certified low-dimensional fiber.