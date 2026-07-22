# R63 DEV bridge scenes — generation notes (Task T-A)

16 fresh DEV-only bridge scenes, 32×32, values in `[0,1]`, generated for the
R24 §4 image-level bridge (`docs/ROUND63_BRIDGE_BUILD_PLAN.md` interface
contract; cohort spec `docs/ROUND63_GPT_ROUND24_RULING_RAW.md` §3.1). Four
groups × four instances: `contour`, `twopop`, `microtex`, `control`.

## Blindness statement

Generation was **blind**. No reconstruction, no arm (SCAT32/RIDGE/FW/RGLI/…),
and no bucket-signal quality metric (PSNR, CNR, bucket contrast `C_u`,
information gap, admission, …) was computed at any point during generation.
The only quantities computed are **image-domain statistics** (mean, std,
min/max, RMS contrast, Michelson contrast, radial spectral centroid), used for
parameter targeting and reporting. These are simple image statistics, not
reconstruction or bucket evaluations. The scenes are frozen (this manifest +
`sha256`) before any arm evaluation; committing the frozen manifest is done by
the owner, not here.

## Files

Per scene `bridge_{group}_{0..3}`:
- `{scene_id}.npz` — field `x`, a 32×32 `float64` array in `[0,1]` (authoritative scene).
- `{scene_id}.png` — 8-bit preview at native 32×32 resolution.
- `BRIDGE_SCENES.json` — manifest: `{scene_id, group, gen_params, sha256, stats}`.
  - `sha256` = SHA-256 of the `.npz` file bytes (freeze fingerprint of the committed file).
  - `gen_params.x_sha256` = SHA-256 of `x.tobytes()` — a reproducible content hash
    that is independent of the `.npz` zip envelope timestamps. Deterministic
    regeneration from the documented seeds reproduces all 16 `x_sha256` exactly
    (verified).

## Seeds (base 650000+, disjoint from 632900 DEV and 633000+ confirmatory)

| group    | scene_ids                | seeds                          | notes |
|----------|--------------------------|--------------------------------|-------|
| contour  | bridge_contour_0..3      | 650000, 650001, 650002, 650003 | frozen contour generator |
| twopop   | bridge_twopop_0..3       | 650100, 650101, 650102, 650103 | R = 4, 8, 16, 32 (one each) |
| microtex | bridge_microtex_0..3     | 650200, 650201, 650202, 650203 | frozen microtexture generator |
| control  | bridge_control_0..3      | 650300, 650301, 650302, 650303 | maze, maze, glyph, glyph |

## Generators

### contour, microtex, control — frozen detail24 family generators (verbatim)

Reused verbatim (import only, no modification) from
`code/round63/detail24.py`, via the pure function
`detail24._generate01(family, seed, side=32)`, which returns a `float64` image
in `[0,1]` as a deterministic function of `(family, seed)` and writes nothing:

- **contour** → `detail24._gen_contour` (`_GEN["contour"]`). Thin
  ellipse/polygon/curve outlines on a nonzero background
  (`detail24.PARAMS["contour"]`, unchanged).
- **microtex** → `detail24._gen_microtexture` (`_GEN["microtexture"]`).
  Rectified band-pass (Difference-of-Gaussians) of a seeded Rademacher ±1
  field, peak-normalized (`detail24.PARAMS["microtexture"]`, unchanged). The
  realized DoG σ-pair per instance (the family's scale/amplitude knob) is
  recorded in `gen_params.realized_sigma_pair`. Realized:
  `bridge_microtex_0/1 = [1.5, 3.0]`, `_2 = [1.2, 2.4]`, `_3 = [0.8, 1.6]` —
  predominantly the larger-σ (lower-frequency, lower-amplitude) end; all four
  have low mean intensity (0.22–0.34), i.e. amplitude at the low end.
- **control** → `detail24._gen_maze` (`_GEN["maze"]`, scenes 0–1) and
  `detail24._gen_glyph` (`_GEN["glyph"]`, scenes 2–3), standard frozen
  parameters (`detail24.PARAMS["maze"]`, `PARAMS["glyph"]`, unchanged).

`gen_params.params_snapshot` embeds the frozen `PARAMS[family]` block for each
scene so the exact ranges are captured in the manifest.

### twopop — analytic two-population witness (R23 Theorem T2-B)

Constructed analytically per `docs/ROUND63_GPT_ROUND23_RULING_RAW.md` §3.2
(Theorem T2-B: bright/uninformative vs dim/informative witness). Each 32×32
scene has two spatially disjoint populations:

- **Population H (bright, low task leverage)** — a large, spatially *smooth*
  bright plateau: a disk of radius `0.34·side` centered near
  `(0.50·side, 0.40·side)` (± a small seeded jitter), Gaussian-smoothed with
  `σ = 2.5` and scaled to `bright_level = 0.80`. High mean intensity, negligible
  fine structure (the "bright but uninformative" population, `Q_H ≈ 0`).
- **Population L (dim, task-carrying)** — thin fine-structure strokes at
  intensity `dim_level`, placed only in the plateau complement (`plateau ≤ 0.5`),
  so the two populations are disjoint. The strokes are a fine cross-hatch:
  vertical `((c+φ_v) mod p)=0`, horizontal `((r+φ_h) mod p)=0`, and diagonal
  `((r+c+φ_d) mod (p+2))=0`, with small period `p ∈ {3,4}` (seeded) and seeded
  phases. This is the low-intensity fine structure that carries the task
  (`Q_L = 1`).
- **Background** = `0.02`.

Construction (per pixel):

```
img = max( bg , bright_level · plateau )        # smooth bright region, plateau∈[0,1]
img = dim_level  where (strokes AND plateau≤0.5) # dim fine structure, disjoint
img = clip(img, 0, 1)
```

**Brightness ratio and amplitude.** The two-population brightness ratio `R`
(Theorem T2-B: `b_H/b_L = R`) is set by construction as
`R = bright_level / dim_level`, i.e. `dim_level = bright_level / R = 0.80 / R`.
The four instances sweep `R ∈ {4, 8, 16, 32}` (one each). The fine-structure
amplitude is therefore `dim_level / bright_level = 1/R ∈ {0.25, 0.125, 0.0625,
0.03125}`, which straddles the nominal `~0.1–0.2 of the bright level` guideline
(the guideline corresponds to `R ≈ 5–10`; the sweep deliberately brackets it on
both sides to probe the `Ω(1)` rank-one gap as `R` grows). The measured
brightness ratio (plateau-core mean / dim-stroke mean) tracks the target closely:

| scene            | target R | plateau-core mean | dim-stroke mean | measured ratio |
|------------------|---------:|------------------:|----------------:|---------------:|
| bridge_twopop_0  |        4 |            0.7796 |          0.2000 |           3.90 |
| bridge_twopop_1  |        8 |            0.7787 |          0.1000 |           7.79 |
| bridge_twopop_2  |       16 |            0.7762 |          0.0500 |          15.52 |
| bridge_twopop_3  |       32 |            0.7774 |          0.0250 |          31.10 |

(The measured ratio is slightly below the nominal because the plateau-core mean
is `≈ 0.78` rather than the `0.80` peak after Gaussian edge-smoothing; the dim
level is exactly on target.) Because the scene power is dominated by the bright
smooth plateau, the twopop radial spectral centroid is low (0.09–0.11) even
though the task-carrying strokes are high-frequency — the intended
"brightness ≠ task-value" structure of the witness.

Full per-instance parameters (`bg_level`, `bright_level`, `dim_level`,
`disk_center`, `stroke_period`, `stroke_phases`, measured means/ratio) are in
`BRIDGE_SCENES.json → scenes[…].gen_params`.

## Contour contrast vs the M1 contour family (verification)

The bridge `contour` scenes use the identical frozen `detail24._gen_contour`
generator with fresh seeds, so they are the same family as the M1 contour
instances by construction. Confirmed by matching RMS-contrast ranges (computed
in-memory on the float generator output; no cache written, no reconstruction):

- bridge contour RMS contrast: **[0.789, 1.141]** (4 scenes).
- M1 contour family RMS contrast: **[0.808, 1.510]** (M1 conf seeds
  633400–633403 + dev seed 632904).

The bridge contours fall within, and at the low end of, the M1 contour-family
contrast distribution — comparable / low range as required.

## Per-scene image statistics

| scene_id          | group    |  mean |   std |   min |   max | RMS contrast | Michelson | spec. centroid |
|-------------------|----------|------:|------:|------:|------:|-------------:|----------:|---------------:|
| bridge_contour_0  | contour  | 0.383 | 0.436 | 0.074 | 1.000 |        1.141 |     0.862 |          0.329 |
| bridge_contour_1  | contour  | 0.534 | 0.436 | 0.127 | 1.000 |        0.816 |     0.775 |          0.390 |
| bridge_contour_2  | contour  | 0.531 | 0.440 | 0.118 | 1.000 |        0.829 |     0.789 |          0.292 |
| bridge_contour_3  | contour  | 0.538 | 0.424 | 0.149 | 1.000 |        0.789 |     0.741 |          0.326 |
| bridge_twopop_0   | twopop   | 0.342 | 0.259 | 0.020 | 0.800 |        0.758 |     0.951 |          0.111 |
| bridge_twopop_1   | twopop   | 0.298 | 0.288 | 0.020 | 0.800 |        0.965 |     0.951 |          0.092 |
| bridge_twopop_2   | twopop   | 0.283 | 0.303 | 0.020 | 0.800 |        1.068 |     0.951 |          0.094 |
| bridge_twopop_3   | twopop   | 0.267 | 0.309 | 0.020 | 0.800 |        1.158 |     0.951 |          0.101 |
| bridge_microtex_0 | microtex | 0.291 | 0.196 | 0.000 | 1.000 |        0.674 |     0.999 |          0.424 |
| bridge_microtex_1 | microtex | 0.221 | 0.174 | 0.000 | 1.000 |        0.784 |     0.999 |          0.383 |
| bridge_microtex_2 | microtex | 0.284 | 0.208 | 0.000 | 1.000 |        0.730 |     1.000 |          0.493 |
| bridge_microtex_3 | microtex | 0.336 | 0.211 | 0.001 | 1.000 |        0.630 |     0.999 |          0.668 |
| bridge_control_0  | control  | 0.238 | 0.426 | 0.000 | 1.000 |        1.788 |     1.000 |          0.376 |
| bridge_control_1  | control  | 0.285 | 0.451 | 0.000 | 1.000 |        1.583 |     1.000 |          0.392 |
| bridge_control_2  | control  | 0.233 | 0.365 | 0.060 | 1.000 |        1.562 |     0.887 |          0.648 |
| bridge_control_3  | control  | 0.155 | 0.284 | 0.060 | 1.000 |        1.826 |     0.887 |          0.354 |

## Reproducibility

Every scene's content is a pure deterministic function of its documented seed
(numpy `default_rng(seed)` only; no wall clock, no global state). Regenerating
`x` from the seeds reproduces all 16 `gen_params.x_sha256` values exactly
(verified in a fresh process). The generation script is kept in the session
scratchpad (outside the repo); the only repo artifacts are the files in
`data/r63_bridge_scenes/`.
