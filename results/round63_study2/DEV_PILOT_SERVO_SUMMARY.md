# Ridge-servo dev pilot — summary (v1 + v2, dev-legal, no gate)

Frozen production path verbatim; Lblob translate patterns; dev maze +
dev glyph; rho_t in {0.6, 2, 6, 22}; nu in {200, 2000}; 2 seeds.
v1 servo hit a divide-by-zero on maze (exact-zero background blocks in
the dev instance -> infinite gains); v2 adds the physical gain-range
guard (0.2-5x). Third pre-freeze guard found by a pilot (after p=2
overshoot and C_u-proxy insufficiency).

## Headline (mode "fixed": just raise global load toward the ridge)

fixed-dwell PSNR_rad, maze:

| rho_t | 0.6 | 2 | 6 | 22 |  ridge rho*(nu) |
|---|---|---|---|---|---|
| nu=200 | 10.64 | 16.95 | **17.32** | 17.05 | 9.9 |
| nu=2000 | 10.64 | 21.36 | **22.81** | 22.63 | 22.2 |

+6.7 dB (nu=200) and +12.2 dB (nu=2000) at equal dwell beyond the
campaign's rho<=2 range; peak locations consistent with the cube-root
ridge. The scalar information map is real at the image level: the
campaigns explored only the foothills.

glyph is prior-limited at this contrast (gains ~ +0.2-0.3 dB; the map
predicts little and delivers little — also consistent).

## The refutation (mode "servo": uniformize every pattern's load at rho_t)

maze: servo LOSES by up to ~7 dB (15.79 vs 22.81 at nu=2000). Diagnosis:
uniformization throttles exactly the bright (structure-viewing) patterns
and cannot raise the near-zero-support ones (gain cap x tiny u). On
scenes with sparse support, structural information lives in the high-u
patterns; equalizing loads starves them.

glyph: servo ~ fixed with a crossing (servo wins slightly only at
rho_t=22, +0.08 dB) — the Jensen-curvature prediction (spread helps on
the convex rising side of J, uniformity near the concave peak) appears
exactly where directional values are homogeneous.

## Refined statement (what survives)

1. ZEROTH ORDER (the prize): set GLOBAL flux so informative patterns sit
   near the ridge — worth +6..12 dB on photon-limited scenes at equal
   dwell.
2. SELF-MEDICATION: under fixed flux, per-pattern load auto-tracks
   brightness, and brightness tracks structural information — the
   direction-load coupling partially water-fills itself. Naive "no
   control" is near-optimal on sparse-support scenes.
3. CORRECTED CONTROL LAW: explicit per-pattern gain should implement
   information WATER-FILLING over the ridge kernel (drive
   structure-viewing patterns toward rho*, starve blank ones) — NOT
   load uniformization. Uniform servo is the special case for
   homogeneous directional value and can lose badly outside it.
4. Guards required in any frozen design rule: finite gain range,
   zero-support handling, and (from earlier pilots) weight-dispersion
   cap; C_u alone is not a sufficient design proxy.

All claims here are dev-level descriptive evidence for the follow-up
method line; nothing enters the OE manuscript (R9 scope lock).
