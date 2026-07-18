# ROUND63 figure preregistration (F1 freeze artifact; spec D2.3 §1/§6)

Frozen BEFORE any confirmatory reconstruction exists. Subjects were chosen by
viewing TRUTH images only (generation-time PNGs); no reconstruction of any
confirmatory instance has been run at freeze time (round-7 §1.2/§7).

## Main figure panel (b) — image array (5 columns per subject row)

Columns (unchanged from D2 §6): safe/full-dwell, safe/short-dwell, high-flux
naive (POISSON-LIN), pre-corrected (PRECORRECT), RQL. Radiometric display
range unified per subject row (min 0, max = 1.2 × truth max, frozen).

Subjects (frozen):
1. **detail_glyph_0** (DETAIL-24 primary cohort, family glyph, seed 631000)
   — full frame 64×64 PLUS zoom crop **x ∈ [2, 42), y ∈ [24, 40)** (40×16,
   the middle glyph row "FBPEW") shown as an inset.
2. **stl_02** (NATURAL-24 secondary cohort) — full frame 64×64, no crop
   (regime-boundary subject; caption carries the round-7 §4 interpretation
   sentence).

Operating points for the columns (frozen): safe/full = (ρ=0.05, ν=2000);
safe/short = (ρ=0.05, ν=50); high-flux columns = (ρ=0.6, ν=50). Seed 0.

## Panel (a): exact count-FI map — generator already frozen
results/round63_theory/fig_a_ridge_map.pdf (committed artifact + generator).

## Panel (c): PSNR_rad vs log T_opt curve families
Subjects: the DETAIL-24 cohort seed-mean curves (RQL, both ρ), with the Q90
crossing markers; per-family thin lines + cohort median band. No subject
selection freedom remains (all 24 shown).

## Panel (d): empirical S_gate vs ρ̄ + exact-FI time-efficiency envelope
All DETAIL-24 images at ρ̄ ∈ {0.05, 0.3, 0.6, 1, 2}; envelope = Fisher-rate
ratio curve; censored images rendered as bounds per the frozen taxonomy
markers. τ±10%/afterpulse/dark bands from S3 anchors.

No figure subject, crop, column layout, or display-range rule may change
after F1 (spec D2.3 §1 blacklist "最终图像/crop").
