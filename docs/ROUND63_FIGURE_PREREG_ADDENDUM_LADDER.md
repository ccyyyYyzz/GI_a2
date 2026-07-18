# Figure prereg ADDENDUM — round-8 six-panel ladder figure (Study 2)

Frozen 2026-07-18 (Berlin), while the k=16 primary grid was still RUNNING
(~48/90 cells) and BEFORE the frozen analyzer (`study2_runner.py --analyze`)
had been executed: no S_gate value, no bootstrap CI, no primary fixed-budget
median, and no per-image k=16 verdict existed when this file was written.

Honest disclosure of what HAD been seen at freeze time: the five no-gate
modes (robustness/controls/fluxmap/separation/scattering) were complete and
their per-k / per-(ρ,ν) AGGREGATE means had been viewed
(results/round63_study2/MODES_DESCRIPTIVE_SUMMARY.md, commit b0a70ad). No
per-family and no per-image reconstruction or metric comparison had been
viewed for any Study-2 stage; no Study-2 reconstruction had been displayed
as an image. Subject choices below therefore cannot have been steered by
per-subject outcomes.

## The six panels (round-8 ruling §Q3) and their remaining display freedom

| panel | content | display freedom remaining → frozen choice |
|---|---|---|
| (a) | exact FI map + cube-root ridge | none (committed artifact fig_a_ridge_map) |
| (b) | k=512/32/16/1 masks + equal-dose bookkeeping | mask index → **row 0 of B^(k)** (first pattern), deterministic |
| (c) | actual C_u and Γ vs k | none (all fluxmap rows aggregated; Γ shown at (ρ=0.6, ν=2000)) |
| (d) | safe/fast PSNR + CNR vs dwell, resolution target | subject + operating points → frozen below |
| (e) | S_gate + fixed-budget gain by occupancy | none (all 24 images; frozen gate/secondary statistics only) |
| (f) | dynamic scattering in S_det units | none (all scattering cells shown) |

## Panel (d) frozen choices

- **Subject: `detail32_spokes_0`** (family spokes, replicate 0, generator seed
  632000+100·f+0). Rationale, a priori: the round-8 ruling names "a USAF
  target"; the frozen Study-2 scene families (round-8 §2, no post-freeze
  additions allowed) contain no USAF chart, and the spokes family (radial
  Siemens-star-like resolution target) is the closest frozen analog. The
  replicate-0 instance is chosen by index convention, matching the Study-1
  figure prereg convention (detail_glyph_0), not by inspection of any metric.
- **Curves**: seed-mean (5 confirmatory measurement seeds) PSNR_rad vs ν
  (all nine frozen dwell points) and CNR vs ν, at k=16, RQL arm, for
  ρ̄=0.05 (safe) and ρ̄=0.6 (fast). CNR cells recorded NA are omitted from
  the mean without substitution.
- **No peak-location claim**: per round-8 §12, the panel must NOT assert the
  image optimum coincides with the scalar Fisher ridge; the caption states
  the ridge is a count-information result and the image optimum an outcome.

## Panel (b) frozen choices

Masks = row 0 of the actual campaign matrices B^(k) (same generator call as
the runner, kind "sparsek", frozen nonce), k ∈ {512, 32, 16, 1} left→right,
rendered binary without interpolation; annotations limited to occupancy %,
on-pixel multiplier n/k, and the equal-dose bookkeeping line
(A=(n/k)B, Φ=ρ̄/τ).

## Inheritance

The Study-1 figure prereg (docs/ROUND63_FIGURE_PREREG.md) remains in force
unchanged for the Study-1 figure of the two-study manuscript. Nothing in
this addendum modifies any gate, endpoint, cohort, or analysis rule; it
constrains DISPLAY only.
