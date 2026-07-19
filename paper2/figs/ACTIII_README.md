# Paper-2 figures — R20 presentation set (GitHub issue #12)

The R20 reader-side ruling replaced the paper-2 figure hierarchy. The
reader-facing set is now produced by **`code/round63/figs/fig_paper2_r20.py`**
(shared style in `code/round63/figs/style_r20.py`). Regenerate from repo root
`D:\GI_another`:

```
D:/Anacondar/anaconda3/python.exe code/round63/figs/fig_paper2_r20.py
```

Individual figures: `... fig_paper2_r20.py mech jitter speed audit actiiib`.

| file | role (R20 item) | column | status |
|---|---|---|---|
| `fig_mechanism.pdf` | Fig. 1 hero: chain + det-vs-jitter + jitter cap (M3) | double | **new** |
| `fig_jitter_validation.pdf` | Fig. 2: peak-vs-`c_v` + engineering warning (M4) | double | **new** |
| `fig_speed_results.pdf` | positive-results: recon + fixed-dwell + elapsed (M5+M15) | double | **new** |
| `fig_audit_supp.pdf` | SUPPLEMENT: guard path + certificate status (M6) | double | **new** |
| `actiii_b.pdf` | SUPPLEMENT: compact dose-only headroom (M6) | single | restyled |
| `actiii_a.pdf`, `actiii_c.pdf` | merged into `fig_audit_supp` | — | **retired** |
| `actiii_d.pdf` | folded into `fig_speed_results` | — | **retired** |
| `actiii_e.pdf` | resource-corner schematic | — | **CUT** (kept for record) |
| `fig_mechanism` (old, serif) | replaced by the new `fig_mechanism` | — | **retired** |

**Retirement note.** `actiii_a/c/d/e` and the old serif `fig_mechanism` were
produced by `fig_actiii_panels.py`, which is now **superseded** and no-ops on a
bare run (it would otherwise clobber the new `fig_mechanism.pdf`/`actiii_b.pdf`).
The retired PDFs remain on disk for provenance; drop them from the main-figure
list. `actiii_e` is cut entirely — its two-corner message is carried by the
two-resource annotation in `fig_speed_results` and by prose.

**Style (R20 M10 / SF7 / SF8).** Sans-serif; reduced palette — neutral gray =
safe/baseline/deterministic, one saturated blue = ridge/high-flux/benefit, one
secondary vermillion = jitter/negative/cost/censoring. Each panel carries a
1–3-word title or result label and at most one causal callout; provenance,
selection rules, IDs, hashes and paths live in the PROVENANCE sections below,
not in the journal captions. **Final absolute text sizes and single/double
widths are approximate** (built at 8.6 cm / 17.8 cm guesses); a resize +
font-match pass follows the Optica template port (M9).

---

## Ready-to-paste journal captions

Install these in `main_m1.tex` (the `.tex` was not modified). Symbol glossaries
and selection rules are deliberately trimmed to the PROVENANCE sections.

### Fig. 1 — `fig_mechanism.pdf` (double column, ~155 words)

> **Hidden dead-time jitter caps the information-optimal flux.**
> (a) Computational single-pixel imaging model (labelled *simulation model*, not
> a bench experiment): a global power multiplier $\Phi$ scales the source; a DMD
> displays scattered-sparse SCAT32 patterns $a_i$ onto the object $x$; a
> non-paralyzable photon-counting bucket detector returns the counts $N_i$ that
> feed the renewal quasi-likelihood (RQL) reconstruction $\hat{x}$.
> (b) After each registered count the detector is held blind; a fixed hold
> $\tau$ (top) and a random hold $B_i$ (bottom) censor different arrivals
> ($\times$), so hidden hold variance destroys count information.
> (c) Count-information efficiency $J(\bar\rho)$ versus detector load at
> $\nu=2000$ for jitter $c_v\in\{0,0.02,0.05\}$: the deterministic optimum sits
> at $\bar\rho^{*}\approx22$, but $c_v=0.05$ moves it to $\bar\rho^{*}\approx5.7$
> and operating the jittered detector at the deterministic ridge forfeits about
> 55\% of the information. Curves, exact $c=0$ law and long-window jitter law;
> markers, Monte-Carlo.

### Fig. 2 — `fig_jitter_validation.pdf` (double column, ~130 words)

> **Numerical verification of the jitter cap.**
> (a) Information-optimal load $\bar\rho^{*}$ versus hold-time jitter $c_v$ at
> $\nu=2000$ on log–log axes. Two independent Monte-Carlo runs (150k and 400k
> samples) and the preregistered matched-asymptotic prediction fall on the
> pooled fitted slope $-0.658$, consistent with the theorem's $-2/3$ power law
> (reference line); marker and line styles separate theorem, prediction and
> Monte-Carlo.
> (b) Engineering warning: the fraction of count information retained when a
> jittered detector is operated at the zero-jitter ridge $\bar\rho\approx22$
> collapses below one-half by $c_v=0.05$ and to a few percent by $c_v=0.2$, so
> the operating load must track the detector's own jitter rather than the
> deterministic optimum.

### Fig. 3 (results) — `fig_speed_results.pdf` (double column, ~170 words)

> **One global power control trades incident dose for elapsed-time speed and
> scene-dependent fixed-dwell quality.**
> (a) Reconstructions (truth, safe $\bar\rho=0.6$, ridge $\bar\rho^{*}$) for the
> median and worst fixed-dwell scenes, with radiometric PSNR printed on each
> crop.
> (b) Per-scene fixed-dwell quality gain $\Delta Q$ at $\nu=2000$: the campaign
> median is $+1.87$\,dB and 19 of 24 scenes are positive; negatives are shown.
> (c) Reconstruction quality versus elapsed optical time $T_{\mathrm{opt}}$ for
> three disclosed scenes: the ridge policy reaches the Q90 quality target far
> sooner than the safe comparator (median case), while a high-flux quality
> collapse (fast-censored) and a flat, uninformative safe range are shown as
> boundary cases.
> (d) Elapsed-time speedup $S_{\mathrm{gate}}$ across all 24 scenes grouped by
> family: median $19.13$ (family-stratified lower bound $18.33$, shaded), 21 of
> 24 above the gate. The single global power control delivers these benefits at
> a $37.1\times$ incident-dose and $2.6\times$ detected-count cost; this is a
> fixed-dwell power-for-time trade, not an equal-photon-budget comparison.

### Supplement — `fig_audit_supp.pdf` (double column, ~120 words)

> **Development-only geometry probe and descriptive certificate status.**
> (a) Along the descending exact-rounding path for the development glyph scene,
> every adaptive OED-DT fraction $\alpha>0$ breaches the dose band and, by
> $\alpha\approx0.06$, the spectral floor, while the A-risk cap never binds, so
> no feasible scene-adapted point exists (development-only evidence).
> (b) Terminal status of the 480 full-stack certificate cells by $(\nu,b)$
> anchor: 299 counterexample, 181 numerically unresolved, and no certified cell
> (reported descriptively; no near-optimality claim).
> (c) Empirical distribution of the finite counterexample primal gap per task
> dimension, which exceeds the $10^{-2}$ near-optimality bar by about two orders
> of magnitude.

### Supplement — `actiii_b.pdf` (single column, ~70 words)

> **Constructive dose-only geometry headroom (development only).** For each of
> the six development families the dose-only primal gap per task dimension is
> positive at both anchors (right axis: $D$-efficiency lower bound $1.6$–$2.5\times$),
> establishing the existence of scene-adapted headroom under the $\pm5\%$
> per-pixel dose band. The construction is support-preserving; it is not a
> confirmatory endpoint, a population estimate, or a dose-only optimum.

---

## PROVENANCE — data sources, selection rules, and exact values

Every plotted value is transcribed from a frozen artifact; no new Monte-Carlo or
campaign computation was run.

### `fig_mechanism` (M3)
- **(a)** drawing; the SCAT32 mask inset is an illustrative scattered-sparse
  Bernoulli pattern (fixed RNG). No real maze/reconstruction thumbnails, dial
  halo, or decorative ticks (removed per M3).
- **(b)** two illustrative arrival trains with a fixed hold `τ` and lognormal
  jittered holds `B_i` (fixed RNG, display only).
- **(c) data.** Monte-Carlo markers are the frozen `ν=2000, n_mc=100k` scalar
  Fisher grid from `results/round63_study2/jitter_sfi_v2_nu2000.log` for
  `c_v ∈ {0, 0.02, 0.05}` (transcribed as `MC_RHO`/`MC_J`). Smooth curves: the
  `c=0` finite-window **exact** count law `code/round63/physics.fisher_exact`
  (peaks at `(6ν)^{1/3}−2/3 = 22.23`, matching the measured peak `22.297`,
  `J=0.9452`); the `c_v=0.02, 0.05` curves are the analytic long-window law
  `J∞ = ρ̄/[(1+ρ̄)(1+c_v²ρ̄²)]` (matches the measured points to <1% at ν=2000).
  Markers: deterministic optimum `ρ̄*≈22.3`; `c_v=0.05` optimum
  `ρ̄*≈5.7` (long-window optimum `c²ρ²(1+2ρ)=1`, Colab 400k peak `5.700`);
  `~55%` = `1 − J(22.3; c_v=0.05)/J_peak(c_v=0) = 1 − 0.4262/0.9452`. Deployed
  ridge load `ρ̄=22.25` (faint blue line; `main_m1.tex` operating table).

### `fig_jitter_validation` (M4)
- **(a) data.** MC run 1 = local zoom `n_mc=150k` continuous (quad-interpolated)
  peaks `results/round63_study2/jitter_zoom.log`: `(0.02,9.645) (0.05,5.154)
  (0.1,3.752) (0.2,2.106)`. MC run 2 = Colab `400k` peaks (pooled-fit record,
  `docs/R14_PREREGISTERED_PREDICTIONS.md`): `(0.02,9.610) (0.05,5.700)
  (0.1,3.399) (0.2,2.051)`. Preregistered matched-asymptotic prediction (R14
  exact optimum): `(0.02,10.4) (0.05,5.69) (0.1,3.50) (0.2,2.30) (0.3,1.63)`.
  Fitted slope `−0.658` = pooled 8-measurement log-log fit; reference line is
  the theorem exponent `−2/3` anchored at the `2^{-1/3} c^{-2/3}` cap constant.
- **(b) data.** Retained fraction = `J(ρ̄=22.297; c_v)/J(ρ̄=22.297; c_v=0)` from
  the same `nu2000` 100k grid: `c_v = 0/0.02/0.05/0.1/0.2/0.3 →
  100/82.0/45.1/17.0/5.1/2.6 %`.

### `fig_speed_results` (M5 + M15)
- **Fixed-dwell `ΔQ` strip (b).** Per-image `dQ` = 5-seed-mean
  `PSNR_rad(RIDGE-SCAT32) − PSNR_rad(SCAT32-060)` at `ν=2000`, transcribed from
  the corrected table `results/round63_m1/CORRECTION_2026-07-19/`
  `M1_VERDICTS_SPEC_CORRECTED_R19.json` (median `1.86692` dB, 19/24 `> 0`).
  Points coloured by sign (blue `≥0`, vermillion `<0`).
- **Reconstructions (a).** GT | SCAT32-060 | RIDGE-SCAT32 seed-0 crops from
  `results/round63_m1/ACTIII_D_RECON.npz`; the printed PSNR is the seed-0
  realization, the row `dQ` label is the 5-seed mean. Disclosed rows (frozen
  selection): median-`dQ` = `m1_chirp_3`, worst-`dQ` = `m1_contour_1` (taste-call
  3: median + worst; best `m1_maze_1` optional/supplement).
- **Elapsed-time curves (c).** Seed-mean `PSNR_rad` versus `optical_time_s`
  (= elapsed `T_opt = M·ν·τ`) per dwell (`ν ∈ {5..2000}`) from the imaging
  shards `results/round63_m1/shards/*.csv`; arm identity from the shard-id
  prefix (`SCAT32-SAFE`, `RIDGE-SCAT32`), **not** the CSV `arm` column (which
  holds the estimator name `RQL`). Curves are equal-weight block-collapse PAVA
  (unit test `[3,2,1]→[2,2,2]`, R19 §2/D3) in `log T_opt`; the Q90 target is
  `safe_min + 0.9·(safe_max − safe_min)` from the safe PAVA curve (reproduces
  the frozen per-scene targets, e.g. chirp_3 `7.4049`). Disclosed scenes:
  median-speed `m1_chirp_3` (`S=19.26`), `FAST_RIGHT_CENSORED` `m1_contour_1`
  (`S=0`), `SAFE_RANGE_UNINFORMATIVE` `m1_microtexture_1` (`S=1`, fitted safe
  range `0.273 dB < 0.50 dB`).
- **`S_gate` strip (d).** Per-image `S_gate`, censoring status and Q90 target
  from the corrected R19 table (median `19.127`, family-stratified LB `18.328`,
  21/24 `> 1`; gate `S=3`). `FAST_RIGHT_CENSORED` = down-triangle, `S=0`;
  `SAFE_RANGE_UNINFORMATIVE` = open square, `S=1`.
- **Two-resource annotation (M15).** Benefit `+1.87 dB` / `19.1×`; cost
  `37.1×` incident dose and `2.6×` detected counts (`main_m1.tex` operating
  table, RIDGE-SCAT32 row). Mandatory R17/R19 statement — fixed-dwell,
  power-for-time, not equal-photon-budget — printed adjacent.
- **Post-hoc note.** The `ν·ρ` incident-exposure sensitivity (median `0.279`,
  no verdict) is NOT plotted here; it belongs to the provenance/supplement.

### `fig_audit_supp` (M6) and `actiii_b` (M6)
- **(a)** guard trajectory cache `results/round63_m1/ACTIII_A_GUARD_TRAJECTORY.json`
  (frozen `oed_design_v5` mixture scan on the DEV glyph fast anchor; single
  predeclared DEV scene, whole path shown; `PATH_FEASIBLE_ALPHA = None`).
- **(b)/(c)** the 480 certificate cells `results/round63_m1/shards/M1_CERT_*.csv`
  (comma-safe overflow reader): 299 COUNTEREXAMPLE / 181 NUMERICAL_UNRESOLVED /
  0 CERTIFIED; by anchor 77/67/84/71 counterexample. ECDF over the 299 finite
  counterexample `primal_gap_lower_per_r` values (median `1.43`); `-inf` cells
  excluded; bar `G_stack/r ≤ 10^{-2}`.
- **`actiii_b`** the frozen 12-cell replication table
  `results/round63_m1/R18_GAP_PROBE_REPLICATION.md`
  (`DOSE_ONLY_PRIMAL_GAP_PER_R`), all six families × two anchors; right axis
  `exp(gap)`. Support-preserving existence result; development-only, no gate.

### Frozen-language preserved (R20 §9)
`development only` / `descriptive` / `0/480 certified` / `no near-optimality
claim` remain on the supplement figures. Corrected R19 values are used
throughout; the original nonconformant analyzer outputs stay only in the
correction bundle (`CORRECTION_2026-07-19/`), never beside the plotted numbers.

---

## Flags for the Optica template port (M9)

- **Resize pass required.** All widths are 8.6 cm / 17.8 cm guesses; re-run at
  the true single/double-column widths after the `opticajnl` port and repeat a
  100%-scale legibility check. Minimum final in-panel text ≥ 7.5–8 pt.
- **Font match.** Built in DejaVu Sans; swap to the journal sans (or matched
  serif) at port and unify across both papers (SF7).
- **`fig_speed_results` crops.** The reconstruction crops are readable but
  modest; enlarge if the ported column width allows (M14 spirit).
- **Reconstruction cache.** `ACTIII_D_RECON.npz` currently holds median/worst
  (+best) crops; if a third (best `m1_maze_1`) row is wanted in-main, it is
  already cached.
