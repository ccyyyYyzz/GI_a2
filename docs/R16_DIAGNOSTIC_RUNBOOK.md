# R16 estimator diagnostic — runbook

Purpose: decide whether the replicated ~-6% end-point low bias in the
jitter-capped ridge-law peaks (`docs/R14_PREREGISTERED_PREDICTIONS.md`) is an
artifact of the histogram log-pmf / finite-difference Fisher estimator, or a
physical/theoretical residual. Prescribed by
`docs/ROUND63_GPT_ROUND16_RULING_RAW.md` §5.2–5.4.

The diagnostic replaces the histogram + central-difference Fisher estimator with
the R16 single-load conditional-score / missing-information identity, which uses
**no histogram and no finite differencing**:

> `U = N - rho L`, `I_N = Var{E(U|N)} = E N - rho^2 E{Var(L|N)}`  (R16 §5.2)

It is **constant-free**: it fits nothing and proposes no ridge-law constant
(R16 §5.4: "No new constant is to be fitted from the present eight peak values
before this diagnostic.").

## Files

- `code/round63/jitter_score_diag.py` — local CPU smoke driver: `--toy`,
  `--toy-poisson`, `--compare`. The full identity derivation and both toy
  closed forms are in its module docstring.
- `code/round63/jitter_score_diag_colab.py` — self-contained Colab driver
  (upload this ONE file). Runs the §5.3 machinery at Colab scale.

## Local validation (already run, for reference)

```
python code/round63/jitter_score_diag.py --toy          # closed-form check
python code/round63/jitter_score_diag.py --toy-poisson  # J=rho anchor
python code/round63/jitter_score_diag.py --compare --nu 200 --cv 0.05 --rho 5.65
```

Toy (single-count truncated channel, exact `I_N = rho^2 e^{-rho}/(1-e^{-rho})`):
`I_MI` lands within ~1% at `N_MC=20k`, mean relative error over 8 seeds
`+0.02% +/- 0.13%` (straddles zero — unbiased in sign; over 40 seeds MI `+0.01%`,
CS `-0.13%`). These are machinery checks only; **not** science.

## Colab run

1. Upload `jitter_score_diag_colab.py` to the Colab session (no other files
   needed; it is self-contained, imports only `numpy` + optional `scipy`).
2. Choose a **CPU high-RAM** runtime (no GPU). Verify with a smoke by editing the
   header to `N_MC = 20_000` if desired, then restore.
3. Run:
   ```
   !python jitter_score_diag_colab.py
   ```
4. Outputs land in `/content/r16_diag/`:
   `r16_diag_summary.json` (all J-curves, peaks, tail masses, alpha/dlog sweeps)
   and `r16_diag_curves.npz`. Download these before the runtime recycles.

What it computes at every load on the R16 §5.3 grids
(`cv=0.02`: 31 pts on [8.5,11.5]; `cv=0.20`: 29 pts on [1.75,2.45]; plus a `c=0`
calibration grid around the exact `rho_*=22.2543`), from **one** common-random-
number draw set per load:

- `J_MI` (missing-information score), `J_CS` (bias-corrected conditional-score
  cross-check) — no histogram, no finite difference;
- `J_HIST` — the audited histogram estimator over the R16 sweeps
  `dlog in {0.0025,0.005,0.01,0.02}` and floor `alpha/N`,
  `alpha in {0,0.05,0.10,0.50,1.00}`;
- peak of each estimator by three-point quadratic (log-rho) **and** shape-
  preserving cubic (pchip), per §5.3;
- how the histogram peak drifts with `alpha` and with `dlog`.

### Scale note (read before treating output as the frozen result)

This prep run uses `N_MC = 400_000` per load (the value set for this task). R16
§5.3's **frozen** spec is `2,000,000` frames per load with 200 block-bootstrap
peak replicates and full `(N,L)` recording. To promote this to the frozen run,
in the header set `N_MC = 2_000_000` and `SAVE_RAW_NL = True`, and add the
block-bootstrap peak loop. At 400k the per-point bootstrap of the peak location
is **not** included here (only per-point estimator CIs in the local `--compare`);
add it for the frozen run.

### Expected wall time

CPU high-RAM Colab, `N_MC=400_000`: roughly **30–60 min** total (≈20–30 s per
load × 79 loads; dominated by generating the `(N_MC × m_max)` exponential array,
`m_max ≈ 2000–2600`). Peak per-chunk memory ≈ 0.5 GB with `CHUNK=8000`. The
frozen `2e6` run is ≈5× longer (≈2.5–4 h) and may need session splitting per
`cv`; the two end points are independent and can run in separate sessions.

## Decision rule — quoted verbatim from the ruling (R16 §5.4)

The ruling's decision rule is **qualitative / pattern-based**. It names target
load values but gives **no explicit numeric tolerance** for "near" / "reproduces",
and **no single scalar pass/fail threshold**. Quoted verbatim:

> - If `I_MI/I_CS` peak near `10.22` and `2.17`, while the histogram peak moves with `alpha` or sample size, the residual is a tail-floor artifact.
> - If the histogram peak moves quadratically with `dlog^2` but is stable to `alpha`, extrapolate to `dlog=0` and attribute it to finite differencing.
> - If the score-identity estimator independently reproduces `9.62` and `2.08`, the discrepancy is physical/theoretical; then the next task is a genuinely uniform second-order renewal Edgeworth calculation, not another fitted constant.
> - No new constant is to be fitted from the present eight peak values before this diagnostic.

Mapping to the three verdicts:

- **Estimator artifact confirmed (tail-floor)** — score peaks (`J_MI`, `J_CS`)
  sit at the R16 predictions `10.22` (cv=0.02) and `2.17` (cv=0.20), **and** the
  audited histogram peak (`hist_peak_vs_alpha`, `tail_mass`) moves with `alpha`
  or with `N_MC`.
- **Finite-differencing artifact** — histogram peak moves as `dlog^2`
  (`hist_peak_vs_dlog`) but is stable to `alpha`; extrapolate to `dlog=0`.
  (Note R16 §5.1 pre-judges this a *secondary* suspect: the Gaussian benchmark
  predicts only a 0.2–0.5% effect and the *wrong* sign for the observed left
  shift.)
- **Physical / theoretical residual** — the score-identity estimator
  **independently reproduces** the measured `9.62` and `2.08`; then the next step
  is a uniform second-order renewal Edgeworth calculation, **not** a fitted
  constant.

The ruling also fixes the estimator construction and tail handling (R16 §5.2):

> For bins with `k_n>=2`, use [the `I_MI` / `I_CS` forms]. ... Pool only the
> extreme left and right tails whose individual bin counts are below 50, and
> report their total probability. The two estimators must agree within
> Monte-Carlo uncertainty.

The driver reports `tail_mass` (total probability in bins with count < 50) at
every load and prints the `J_MI − J_CS` gap, so the "agree within Monte-Carlo
uncertainty" condition can be checked directly.

## What NOT to do (R16 §6, frozen claim discipline)

Not permitted yet: a new empirical coefficient replacing `2^(-1/3)`; a
lognormal-specific `4/3` multiplier; claiming the residual proves a new
universality subclass; treating histogram-Fisher peak locations as exact without
this score-identity diagnostic.
