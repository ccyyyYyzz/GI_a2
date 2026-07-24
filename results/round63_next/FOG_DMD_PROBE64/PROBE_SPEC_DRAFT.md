# Beyond-Modulator-Band Sealed Probe — SPEC DRAFT (R39)

> ⚠️ **SUPERSEDED IN PART BY THE R39 RULING** (`docs/ROUND63_GPT_ROUND39_RULING_RAW.md`, landed
> after this draft). The ruling changes: geometry (medium band `k_w=k_p`, pattern `|k|≤5`, primary
> claim `≤1.8 k_p`, bounded not lognormal medium, photon ceiling `1e4`), estimator (retire EM-MLE →
> lifted-covariance GLS + Fisher-scored refinement), arms (seven, §6.3), banks (six, §6.2), and adds
> **P1 profiled-Fisher prognosis as the prelaunch go/no-go**. **P1 FAILED** (see `RISK_GATES_REPORT.md`
> / `P1_results.json`): kill-tree node 2 → STOP before reconstruction. This machinery therefore stands
> as the *structure* to realign IF the operator revives at a narrower claim (`≤1.25 k_p`) after the
> mandated MC cross-check. Sections below reflect the pre-ruling draft geometry; treat the ruling as
> authoritative where they differ.

**Status:** DRAFT machinery, build-only. The confirmatory is **UNRUN**. Every bar threshold is an
`[R39]` placeholder in `bars64.py::THRESHOLDS`; the external referee's ruling sets them, and
**freezing is a one-edit step** (edit those nine numbers and nothing else).

**Purpose.** Adjudicate the claim: *the scattering medium's fine speckle band is a synthesized
aperture bolted onto a band-limited modulator, so beyond-modulator-band scene content is
recoverable from the bucket temporal covariance alone — no reference arm, no transmission matrix,
no wavefront sensor, no speckle calibration.* The mean (first-moment) channel is provably blind to
that content; only the covariance channel carries it.

This spec is the exam. It is designed so a **kill at any dev-grade risk gate (G1/G2/G3) or any
sealed bar (B1–B6) ends the flagship claim with no repair round.**

---

## 1. Frozen geometry (the declared home)

| choice | value | rationale |
|---|---|---|
| scene grid | `n=64`, `N=4096` | publication scale (R39 ask 5 minimum) |
| pattern band | `fx,fy ≤ PB=8` | the DMD/modulator pixel-limit wall (~81 band DOF) |
| medium band | `1 ≤ i+j ≤ 16` (2D-DCT), `db≈152` | fine speckle beyond the modulator band |
| reachable (aperture law) | pattern box ⊕ medium band ≈ radial `24` | Minkowski sum; the super-resolution zone is coverage ∖ pattern-box |
| forward model | `b_{i,t} = ⟨P_i ⊙ w_t, x⟩` + Poisson | multiplicative thin dynamic screen **at/near the object plane** (shower-glass regime) |
| medium law | lognormal OU: `w_t = exp(U z_t − v/2)`, `z_t` OU | positive, `E[w]=1` per pixel; the declared low-dim class |
| OU time | `τ = 8` (`ρ = e^{-1/8}`) | independent-realization currency; lags ℓ≥1 shot-free |
| field RMS | `σ_f = 0.30`, per-coeff sd `= σ_f√(N/db)` | pixelwise RMS fixed independent of `db` (ruling normalization) |
| patterns | `M = 128` band-limited incoherent **nonneg** random; sealed seed `10` | incoherent random codes are load-bearing (Fourier-lattice codes collapse the estimator) |
| budget | `T = 4096` epochs, `1e5` photons/bucket, `3` replicates | `[R39]` |
| metric | beyond-band NMSE on the super-resolution zone (2D-DCT power, scalar-gauge) | isolates the claimed content |

All choices are shared with the 16×16 evidence cell (`bb_common.py`) and the E5a/E7a production
solvers (`../FOG_DMD_PROBE/fog_tracker.py`) so numbers are directly comparable across scales.

The declared geometry is the **binding validity boundary**: the multiplicative thin-screen model is
exact only for a screen at the object plane. Departure toward convolutive (finite propagation)
mixing is a mismatch axis (G2 / bar B3), **not** the home model.

---

## 2. The six arms (`arms64.py`)

| arm | patterns | readout | role |
|---|---|---|---|
| **fixed+MLE** (production) | fixed sealed bank | moment init → Kalman-EM marginal likelihood | the deployable estimator |
| **fixed+moment** | fixed sealed bank | E5a covariance moment matching | robustness / global basin |
| **fresh+mean** (THE WALL) | fresh each epoch | mean route | channel-impossibility control (must stay ≈1.0) |
| **fresh+cov (GLS)** | fresh each epoch | covariance route, GLS-weighted | design-half comparator (concentration principle) |
| **oracle** | fixed bank | medium known | nondeployable ceiling |
| **classic averaging** | fixed bank | time-averaged bucket inversion | ordinary GI baseline (band-limited → wall) |

**Equal budget (binding).** Every arm spends identical `M × T` exposures at identical photons/bucket;
only the pattern strategy (fixed vs fresh) and readout (mean/cov/MLE) differ. The **concentration
principle** claim is exactly `fresh+cov` beyond-band NMSE ≥ `fixed+moment` at equal budget (bar B5).

**Seed hygiene.** Disjoint partition salts `{calibration:810, confirmatory:800, mismatch:820,
oracle:830}` (distinct from the DLGI 700/710/720); independent `SeedSequence` children for
(medium, Poisson); all arms share ONE medium realization per record (fair paired comparison).

---

## 3. Sealed banks (`bank_gen64.py`)

Four disjoint, blind-generated partitions; **no scene, no seed, no natural image crosses banks**
(verified: determinism + all-distinct PASS).

| partition | salt base | role (touched ONLY for) | instances/stratum |
|---|---|---|---|
| calibration | 810000 | aperture-law prediction map + `coeff_sd` gauge check | 2 |
| confirmatory | 800000 | final beyond-band endpoints | 3 |
| mismatch | 820000 | F5 declared-law / geometry mismatch bars | 2 |
| oracle | 830000 | oracle-ceiling anchor (disjoint from confirmatory) | 2 |

**Strata:** `witness` (in-band random + explicit super-res-zone spikes — the load-bearing beyond-band
stratum), `natural` (a fixed skimage image, globally-unique per bank slot, downsampled 64×64,
sha256-recorded), `onef` (1/f^β proxy), `smooth` (in-band control), `texture` (band-pass, super-res
annulus stress). Every image is 64×64 float64 ∈ [0,1], byte-reproducible, sha256 in the manifest.

---

## 4. Bars (`bars64.py`, R34-style binding kill bars)

| bar | statement | `[R39]` threshold (placeholder) |
|---|---|---|
| **B1** aperture-map coverage | measured recoverable spectrum matches the aperture-law map (in-coverage recovered, out-coverage stays blind) | separation ≥ 3×; in-coverage err ≤ 0.70 |
| **B2** beyond-band recovery | fixed+MLE beyond-band NMSE below cap AND beats the wall by a margin | NMSE ≤ 0.70; wall-margin ≥ 0.20 |
| **B3** mismatch degradation (F5) | every mismatch axis degrades ≤ cap | degradation ≤ 25% |
| **B4** multi-start agreement | covariance data-only multi-start std small (injective / collusion-free) | std ≤ 0.05 |
| **B5** equal-budget + concentration | byte-identical exposures/photons across arms; fixed ≥ fresh | exposures byte-equal; concentration ≥ 0 |
| **B6** channel impossibility | fresh+mean AND classic-avg stay at the wall (~1.0) | wall NMSE ≥ 0.90 |

**Kill tree.** Any bar kill ends the flagship claim (demote to a theorem/map note; no repair round).
B2m/B6 are intentionally coupled: a wall breach trips both (defense in depth). The bar arithmetic is
unit-tested on synthetic passing AND failing inputs for every bar (`bars64.py --selftest`, ALL PASS).

### 4.1 F5 mismatch axes (bar B3, `sealed_shard_runner.run_mismatch_cell`)

`fine_rotation` (declared band rotated 10/20% within the fine annulus), `band_width` (declared
±2 in `i+j`), `radial_profile` (true 1/f vs declared flat), `geometry_alpha` (convolutive blend
α∈{0.1,0.25,0.5}), `wrong_tau` (true τ/2, 2τ), `wrong_sigma_f` (true 0.2/0.4). Data use the true
axis value; the estimator keeps the frozen declared law. Degradation = beyond-band NMSE relative to
the matched estimator at the same budget. **Any axis > 25% is an F5 kill (simulator-law artifact).**

---

## 5. Dev-grade risk gates (run BEFORE sealing — `RISK_GATES_REPORT.md`)

The three fast-kill gates are prerequisites; a kill here saves the entire sealed campaign.

- **G1 fine-band mismatch** (`gate_g1.py`) — the sharpest risk.
- **G2 geometry mixing** (`gate_g2.py`) — multiplicative vs convolutive.
- **G3 scale + natural images** (`gate_g3.py`) — does the law/effect survive to 64×64.

See `RISK_GATES_REPORT.md` for verdicts. The sealed bars inherit the gate outcomes: the radial-
profile axis (G1(d)) and the scale washout (G3) are the load-bearing risks the confirmatory must
resolve; their `[R39]` thresholds must be set with those results in view.

---

## 6. Fleet plan (5-shard Colab lane, `sealed_probe_planner.py` + `sealed_shard_runner.py`)

Built in the proven 5-session lane format from the start, reusing the `code/round63/colab/`
templates (`make_bundle` / `live_create_sessions` / `live_launch_all` / `session_driver` +
`remote_lane` / `live_rebind` / `live_watch` / `live_fetch_all` / `live_release`).

- **Shard↔account parity (FROZEN for the official run):** even shards → pro1 (2×L4), odd shards →
  pro2 (3×L4). 5 shards: `{0:pro1, 1:pro2, 2:pro1, 3:pro2, 4:pro2}`.
- **VM-autonomous driver:** task list in the shard manifest on the VM; heartbeat JSON + per-cell
  checkpoint (CSV + `_meta.json`) so token/network loss never kills computation; META-as-truth
  resumability (re-invocation skips done cells); wall-budget guard returns 2 (resumable).
- **Watchdog cycle order:** rebind → heartbeat check → idempotent fetch (per COLAB guide).
- **Target:** once R39 freezes the bars, the confirmatory launches on 5 VMs and completes in 2–4 h.
- **On completion:** release all five VMs (`state.client.unassign`), kill keep-alives, report
  `assignments=NONE`.

Discretionary/dev usage burns pro2 first (units expire sooner). The dev-grade gates ran **locally**
on the RTX 4060 (free, jobs ≪ 1 h) per the coordinator's own local-vs-Colab condition; no live Colab
calls were made this session.

---

## 7. What is frozen vs. `[R39]`-pending

**Frozen now:** geometry (§1), the six arms + equal-budget accounting (§2), the sealed banks +
seed hygiene (§3), the bar *structure* and F5 axes (§4), the fleet parity (§6).

**`[R39]`-pending (one-edit freeze):** the nine bar thresholds in `bars64.py::THRESHOLDS`; the band-
extension factor to *claim* (2× linear reach is the design target; the dev gates measure the
achieved reach); the physical-unit photon/realization budget statement for the declared geometry.
