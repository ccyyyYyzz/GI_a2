# BUILD_STATUS — paper_aperture figure layer

**Manuscript:** *The information aperture of a bucket record* (R40 frozen architecture).
**Build:** 2026-07-24. **No git commits** (coordinator reviews).
**Env:** Anaconda base, Python 3.11.5, matplotlib 3.7.2, numpy 1.24.4 (local).

## Pipeline (R40 writing order steps 1–3)

```
repo evidence  ──►  source_tables/extract_source_tables.py  ──►  source_tables/*.{csv,json}
                                                                        │  (RULE ZERO boundary)
                                    figures/*.py  ◄──── read ONLY ──────┘
                                         │
                                         └──►  figures/*.{pdf (vector), png (preview)}
```

Regenerate everything (cwd = repo root `D:\GI_another`):
```
python paper_aperture/source_tables/extract_source_tables.py
python paper_aperture/figures/fig1_atlas.py
python paper_aperture/figures/fig2_certificate.py
python paper_aperture/figures/fig3_boundary.py
```

## Phase 1 — frozen source tables  ✅ COMPLETE (10 tables)

| table | status | notes |
|---|---|---|
| `T_fig2_ridge_curves.csv` | ✅ COPIED | fisher_ridge.csv verbatim (7 ν × 49 ρ) |
| `T_fig2_ridge_points.json` | ✅ COPIED+FROZEN | peaks, ridge law, jitter cap, operating point |
| `T_fig2_cohorts.json` | ✅ FROZEN | +1.87 dB (LB 1.41, 19/24); 19.13× (LB 18.33, 21/24) |
| `T_fig2_closure.csv` | ✅ FROZEN | 5 interventions + materiality bar |
| `T_fig3_aperture.json` | ✅ COPIED+FROZEN | 64×64 18.76× + E8 mean-wall duel; 16×16 pair flagged missing |
| `T_fig3_p1.json` | ✅ COPIED | P1 kill + T_eff^-1/2 sweep |
| `T_fig3_map.csv` | ✅ COPIED | full 81-cell grid |
| `T_fig3_map_meta.json` | ✅ COPIED | grid meta + MC cross-check |
| `T_fig3_pocket.json` | ✅ COPIED | pocket cell; **BRANCH = PENDING** |
| `T_fig1_atlas.json` | ✅ DERIVED | master-atlas summary + derivation notes |

Clean extractions: all ridge/cohort/closure/aperture-64/P1/81-cell/pocket data.
**One missing raw input** (precisely named): the 16×16 aperture pair 0.073/1.213 does
not exist in `D:\GI_another` (see `FIGURE_SOURCES.md` → Missing inputs). Substituted with
the well-sourced G3 64×64 18.76× separation + E8 exact-mean-wall duel.

## Phase 2 — figures  ✅ COMPLETE

| figure | outputs | status |
|---|---|---|
| Fig 1 (gate) | `fig1`, `fig1_pocketA`, `fig1_pocketB` (.pdf+.png) | ✅ 3 design passes |
| Fig 2 (mean certificate) | `fig2` (.pdf+.png) | ✅ ridge + cohorts + closure |
| Fig 3 (covariance boundary) | `fig3`, `fig3_branchA`, `fig3_branchB` (.pdf+.png) | ✅ 4 panels, both branch variants |

Publication settings: vector PDF (`pdf.fonttype 42`, editable text) + 200-dpi PNG
preview; Okabe-Ito colorblind-safe palette; sans-serif; no chartjunk.

### Fig 1 three-pass self-critique (the gate — must read as a MAP without a caption)

- **Pass 1:** legible but cluttered — italic gap text ran under the ✕/◇ markers, the
  "killed 1.8×" label hit the legend, and the mean-vs-covariance contrast was muted
  (dark bar filled ~55% of the pale bar, not an instant "thin rim").
- **Pass 2:** added per-row verdict tags (*support ≈ supply* / *support ≫ supply*),
  hatched the empty "support without supply" gap, moved the legend to a horizontal
  strip below, brightened the usable-edge. Residual: the amber pocket leader crossed
  up into the mean row and the covariance verdict tag was cramped against the wall.
- **Pass 3 (final):** relocated the covariance verdict tag over the empty gap, dropped
  both covariance marker callouts below the row with leaders, pushed the 2× aperture
  arrow to the bottom. **Result:** a referee sees in <10 s — one record, two channels;
  mean fills its support (certified green point), covariance has 2× formal support but
  the extension is empty except one marked pocket and one killed claim. *support is not
  supply.* A map, not a method anthology.

## Known limitations / open items

- **Dense-grid tick labels** in Fig 3c (9×9 atlas) are 5.0–5.5 pt — below the 7 pt
  target, inherent to an 81-cell grid rendered in a main panel. All titles, axis
  labels, annotations and numbers are ≥6 pt (mostly 7–8 pt). Coarsen to macro-blocks
  if the target journal enforces 7 pt on every glyph.
- **Branch-gated variants** (`fig1_pocketA/B`, `fig3_branchA/B`) are prebuilt; when the
  blind lifted-GLS pocket verdict lands, set `T_fig3_pocket.json → BRANCH` and select
  the matching variant as canonical `fig1`/`fig3` (A = POCKET_ACHIEVED filled;
  B = ESTIMATOR_GAP_PERSISTS hollow). Current canonical = PENDING (hollow, neutral).
- **No Optica/IEEE template port** yet (R40 step 9); figures are at 17.8 cm double-column.
- **Captions / one-sentence claims** (R40 step 4) not written — this task delivered the
  figure + source-table layer only.
