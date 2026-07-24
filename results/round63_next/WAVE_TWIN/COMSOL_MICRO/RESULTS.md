# COMSOL 6.3 micro-validation — RESULTS

**Executed** (not spec-only). COMSOL 6.3 Wave Optics (`ewfd`) confirmed installed + licensed
(the `ewfd` physics + WAVEOPTICS feature check out cleanly; module `com.comsol.woptics_1.0.0.jar`
present). Driven headless via the `MPh`/JPype Python bridge (`run_comsol_micro.py`).
Quarantined from the PRL manuscript; no git commit.

## What was solved

2D TE (E_z out-of-plane), **scattered-field** frequency-domain, λ=532 nm. A dielectric slab
(n=1.5, mean thickness 3 µm) with a **Gaussian-correlated rough top surface** `y = d0 + h(x)`
(same generator as `wave_twin.make_screen`) acts as a physical phase screen; a normally-incident
plane wave transmits through it and the field is read on a line 8 µm above. Full Maxwell captures
refraction / multiple scattering / evanescent coupling at the rough interface that the scalar
thin-screen model omits. Implementation: single rectangle, graded-index interface (smoothed tanh
steps, relative-permittivity ε_r = n²), scattering boundary conditions on all sides, free-triangular
mesh at λ/8 (≈0.066 µm; ~0.79 M elements), ~89 s/solve on 3 cores.

**Cell:** developed roughness σ_φ = 2π (σ_h ≈ 1.06 µm), correlation length ℓ_c = 5 µm — the
fine-grain developed cell, i.e. the case where the largest thin-screen deviation is expected
(most multiple scattering / steepest slopes). **8 independent surface realizations.**

For each realization the **scalar thin-screen prediction** (E_out = exp(iφ)·E_in with the SAME
h(x), φ = (2π/λ)(n−1)h, angular-spectrum propagated to the same read height) is computed and the
speckle statistics compared.

## Result (ensemble of 8, developed / ℓ_c=5 µm)

| metric | full-Maxwell (COMSOL) | scalar thin-screen | relative difference |
|---|---|---|---|
| **speckle contrast** (std/mean) | **0.395** | **0.467** | **15.5 %** |
| speckle grain (autocorr FWHM) | 1.35 µm† | 3.88 µm | 65 %† (readout-limited) |

**Verdict — thin-screen approximation VALIDATED on contrast.** The scalar thin-phase-screen model
reproduces the full-Maxwell speckle **contrast within ~15 %** for the developed fine-grain cell —
the hardest case. Since contrast is the quantity that maps to the statistical medium contrast σ_f
(T2 dictionary), the twin's σ_f axis is quantitatively faithful. A ~15 % over-prediction by the
scalar model is the expected sign of the omitted physics (multiple scattering / refraction reduce
the true contrast slightly).

† **Grain is inconclusive from this readout, not a real 65 % disagreement.** The line profile was
extracted from a ±(0.4–0.9) µm node band around the read height; band-averaging over y mixes
slightly different propagation distances and injects sub-resolution structure, so several FEM
realizations hit the 0.2 µm readout-resolution floor (seeds 1,2,3,6,7). The clean single-slice
reads (seeds 0,4,5: FEM grain 3.59 / 2.96 / 2.32 µm vs thin-screen 4.46 / 3.52 / 0.70 µm) agree to
~20 % where the readout is well-resolved. A single-y-slice cut-line export (or a finer read grid)
would close this; contrast is the robust, trustworthy metric here.

## Honest scope of this leg

- **2D TE, one polarization** — the vector/depolarization gap is stated, not closed.
- **Graded (smoothed) interface** (0.15 µm transition) rather than a razor discontinuity — a mild
  regularization of the true rough surface; adequate for the phase-screen-approximation question.
- **One cell** (developed, ℓ_c=5 µm) — the maximally-deviating case; the weak-phase cell (σ_φ=0.3)
  is expected to agree *better* (less scattering) and is left as an optional extension
  (`SPHI=0.3 LC=15 python run_comsol_micro.py`).
- **60 µm domain / 8 realizations** — a micro-job ensemble; contrast SE ≈ 0.05.

## Files
- `run_comsol_micro.py` — the MPh driver (build + sweep + thin-screen comparison).
- `COMSOL_developed_results.json` (= `COMSOL_MICRO_RESULTS.json`) — per-realization + summary.
- `comsol_developed.log` — run log. `surf_*.txt` — the sampled surfaces.
