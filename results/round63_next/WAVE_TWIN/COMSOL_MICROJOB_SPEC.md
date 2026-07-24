# COMSOL 6.3 micro-validation spec — the thin-phase-screen approximation

**Purpose.** The wave twin models the diffuser as a **scalar thin random phase screen**
(a multiplicative `exp(iφ(r))` applied in one plane, then scalar angular-spectrum
propagation). This is the load-bearing approximation behind every twin result. This job
is a **first-principles full-Maxwell check** of that approximation: solve 2D rough-surface
scattering for a Gaussian-correlated dielectric profile and compare the **far-field
speckle statistics** (contrast, grain/angular-correlation) against what the thin-screen
model predicts for the same surface. **Spec only — do not run unless the Wave Optics or RF
module license is confirmed present (see `SPEC_ONLY_LICENSE.md` if absent).**

Scope is deliberately tiny (a "micro-job"): a tens-of-µm domain, one polarization, a
handful of surface realizations. It validates the *physics of the approximation*, not the
full bench. Time cap: ~2 h wall-clock.

---

## 1. Physics module and study

- **Module:** Wave Optics (`ewfd`, Electromagnetic Waves, Frequency Domain) — preferred;
  RF module (`emw`) is an acceptable substitute at these sub-mm scales.
- **Study:** Frequency Domain at f = c/λ, λ = 532 nm (f = 5.635e14 Hz).
- **Dimension:** 2D (TE: out-of-plane E_z; the scalar-diffraction analog). One polarization
  only — the vector/polarization gap is stated as a twin limitation, not closed here.
- **Formulation:** scattered-field. Background field = the incident illumination; solve for
  the scattered field so PML truncation is clean.

## 2. Geometry (tens-of-µm domain)

- **Domain:** 60 µm (x, transverse) × 40 µm (z, propagation), plus PML frames.
- **Rough interface:** a dielectric boundary at mid-domain between vacuum (n=1) and a weak
  dielectric (n=1.5, borosilicate-like), with a **Gaussian-correlated height profile**
  `h(x)`: zero-mean, RMS roughness σ_h, Gaussian autocorrelation length ℓ_c.
  - Realize `h(x)` as an interpolation function (import a 1D sampled profile generated
    exactly as in `wave_twin.make_screen`: filtered white noise, Gaussian kernel width ℓ_c,
    scaled to σ_h). Deform the boundary with a **Deformed Geometry / boundary displacement**
    or build the curve from the interpolation function.
  - **Two roughness regimes** (match the twin): weak `σ_φ = (2π/λ)(n−1)σ_h ≈ 0.3 rad`
    (→ σ_h ≈ 0.3·λ/(2π·0.5) ≈ 51 nm) and developed `σ_φ = 2π` (→ σ_h ≈ 1.06 µm).
  - **Correlation lengths:** ℓ_c ∈ {5, 15} µm (the accessible fine-grain regime; 50 µm
    optional if time permits).

## 3. Sources and boundary conditions

- **Incident field:** a normally-incident Gaussian beam (waist ~20 µm, ≪ domain, ≫ ℓ_c) or a
  plane wave with a super-Gaussian transverse apodization to avoid edge diffraction.
- **PML:** on all four outer frames (≥ 8 layers, λ/2 stretch) to absorb outgoing/scattered
  waves — mandatory to read a clean far field.
- **No periodic BCs** (a finite illuminated patch, not an infinite grating) — the twin's
  diffuser is finite/quasi-static.

## 4. Mesh (the accuracy-critical choice)

- **Maximum element size:** ≤ λ/10 in vacuum = 53 nm; ≤ λ/(10·n) = 35 nm in the dielectric.
- **Boundary-layer refinement** along the rough interface: ≤ λ/20 (27 nm) and enough to
  resolve σ_h and ℓ_c (≥ 8 elements per ℓ_c → ℓ_c=5 µm needs ≤ 0.6 µm there — the λ/10 rule
  already dominates).
- Estimated DOF: ~1–3 M (2D, one field component) — fits a single frequency-domain solve in
  a few minutes to tens of minutes on a workstation; the cost is the **realization sweep**.

## 5. Parametric sweep (the Monte Carlo)

- **Surface realizations:** ≥ 30 (target 50) independent `h(x)` seeds per (σ_φ, ℓ_c) cell —
  the ensemble that produces speckle statistics. Use `Parametric Sweep` over a seed
  parameter driving the interpolation-function import (regenerate `h(x)` per seed).
- **Cells:** {weak, developed} × {ℓ_c 5, 15 µm} = 4 cells × ~40 seeds ≈ 160 solves.
  At ~30–45 s/solve (2D, ~1–3 M DOF) → ~1.5–2 h. **This sets the 2 h cap; drop ℓ_c=15 or
  cut to 30 seeds if over.**

## 6. Outputs and the comparison to the thin-screen model

For each realization, on a far-field arc / a line monitor at z = z_far (e.g. 30 µm past the
surface, inside the physical domain before PML):

1. **Intensity `I(x) = |E_z|²`** along the monitor line.
2. **Speckle contrast** `C = std_x(I)/mean_x(I)` (ensemble-averaged over seeds). Expect
   `C → 1` developed, `C ≪ 1` weak — the twin's T2 result.
3. **Angular / spatial speckle grain** = FWHM of the intensity autocorrelation along the
   monitor (ensemble-averaged) → compare to the twin's grain and to `λ z_far / D_illum`.
4. **Far-field angular spectrum** (near-to-far-field transform, `emw.Efar`/`ewfd.Efar`):
   the scattering angular width `θ_s` vs the thin-screen prediction `θ_s ≈ σ_φ·λ/(2π ℓ_c)`.

**Thin-screen reference (compute alongside, no COMSOL needed):** apply
`E_out(x) = exp(i φ(x))·E_in(x)` with the SAME `h(x)` (φ = (2π/λ)(n−1)h), propagate to
z_far by scalar angular spectrum, and extract the same three statistics. The validation
metric is the **relative difference in (contrast, grain, θ_s)** between full-Maxwell and
thin-screen, ensemble-averaged.

## 7. Pass / interpret

- **PASS (approximation validated):** contrast and grain agree within ~10–15 %, and θ_s
  within ~20 %, for both weak and developed regimes at ℓ_c ≥ λ. Then the twin's thin-screen
  diffuser is quantitatively faithful for the transverse-scale conclusions (T2–T5).
- **EXPECTED DEVIATIONS (report, don't hide):** the thin-screen model omits (i) multiple
  scattering within the roughness, (ii) shadowing at large `σ_h`/steep slopes, (iii)
  polarization/depolarization (absent in 2D TE), (iv) evanescent/near-field coupling within
  ~λ of the surface. These grow with σ_φ and with ℓ_c → λ; the developed + fine-ℓ_c cell is
  where the largest deviation is expected and is the most informative point.
- **Deliverable:** `COMSOL_MICRO/RESULTS.md` — a 4-row table (cell × {contrast, grain, θ_s})
  with full-Maxwell vs thin-screen vs relative error, plus one E_z speckle field image per
  regime.

## 8. Provenance / constraints

- COMSOL install located at `D:\Program Files\COMSOL\COMSOL63\Multiphysics`.
- **License gate:** run only if a Wave Optics (`ewfd`) or RF (`emw`) feature is licensed
  (checked via `comsol.exe` product list / a headless `mphstart` probe). If not licensed,
  this file stands as the runbook and `SPEC_ONLY_LICENSE.md` records the gate.
- **Batch invocation (if licensed):**
  `comsolbatch.exe -inputfile rough_screen.mph -outputfile out.mph -study std1 -batchlog b.log`
  with the model built from this spec (or via a `.java`/Model-Builder script).
- Quarantined from the PRL manuscript (R43 green light); no git commit; results feed the
  bench/collaboration phase only.
