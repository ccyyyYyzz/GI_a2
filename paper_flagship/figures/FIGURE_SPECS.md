# FLAGSHIP FIGURE SPECS — transcribed from R35 §6

**Source of record:** `docs/ROUND63_GPT_ROUND35_RULING_RAW.md` §6 (hero) and the
per-section figure specs (§3.2, §3.3, §4.3, §5.5). Figure inventory register:
`docs/FLAGSHIP_MATERIALS_MAP.md` §6.

**HARD RULES (R35 §6.2, §2.1):**
- **Exactly five main figures. No sixth figure.** If a result does not fit these
  five claims it belongs in the supplement.
- Each figure earns exactly one main claim.
- Caption limit 120–160 words; **no provenance narratives in captions** (R35
  cut 13 — captions explain the scientific reading, Methods/Supplement establish
  the audit trail).
- Build in matplotlib per repo conventions (generators live under
  `code/round63/figs/*.py`; style below).

**Repo figure style (from `paper/figs/FIG_MECHANISM_P1_README.md`,
`paper2/figs/ACTIII_README.md`):** sans-serif; neutral **gray = safe**;
**blue = ridge / high-flux**; **vermillion = dead-time / negative**; hero and
closure figures use **≤ 4 colors** (R35 §6.1). Sizes ≈ 8.6 cm (single col) /
17.8 cm (full width); **Optica resize/font pass is gap G6 (pending, M9)**.

**Branch note (R35 §8):** Fig. 1 Stage 3 and Fig. 5 are **PASS-BRANCH only**.
On a DLGI fail, remove Fig. 1 Stage 3 (hero ends at the count-only boundary) and
drop Fig. 5 entirely (four figures, Fig. 4 becomes the culmination).

---

## Fig. 1 — HERO: "the 10-second paper" (limit → reach → close → reuse)

**Assigned section:** §1 Introduction. **Claim:** a dynamic hidden channel
simultaneously limits scene information and can carry medium information.
Full-width, three-stage visual read left → right (R35 §6.1).

### Stage 1 — one hidden-state bucket experiment
Programmable pattern → object → dynamic channel → bucket count `Y_t`. Above the
channel, one hidden-state trace with **two labels only**: *detector live time*
and *medium gain*. Label the optical chain explicitly **"simulation model"** —
**not** a photorealistic bench graphic (R35 §7.1).

### Stage 2 — limit and closure
Split the bucket count into: a **red scene-information arrow** terminating at a
finite ridge curve; a **global source dial** placed on the ridge; a thin **gray
wall** labeled *count-only scene boundary* with **five small unlabeled
intervention ticks** beneath it. **Only two numeric callouts:** `+1.87 dB` and
`19.13×`.

### Stage 3 — reuse  **(PASS-BRANCH ONLY)**
The same raw bucket stream branches into two output cards: **scene**; **medium:
`t_c`, CV**. Between them one compact reciprocity symbol: the two ledgers share
the same normalized confusion spectrum `κ_i²`. If DLGI passes add **0 extra
exposures** and **≤1.5× pilot precision**. *If it fails, this stage is removed,
not relabeled vaguely.*

### Constraints (R35 §6.1)
No equations longer than one line; no algorithm names; no campaign labels; no
image gallery; **no more than four colors**.

### The ten-second referee takeaways
1. Hidden dynamics cap scene information. 2. A global dial reaches the cap and
yields a large speed gain. 3. Count-only adaptations do not materially move the
cap. 4. The same record can instead become a medium instrument.

- **Source materials available now:** `paper/figs/fig_mechanism_p1.pdf` (SPI
  chain + dead-time schematic — Stage 1 seed); `paper2/figs/fig_mechanism.pdf`
  (chain + jitter cap — Stage 2 ridge seed); `results/round63_next/DUAL_LEDGER_PROBE/figs/fig_dual_ledger.png`
  (Stage 3 dual-output seed). The `+1.87 dB`/`19.13×` callouts are frozen.
- **Awaits:** the unified three-stage arc composition itself (**gap G3** — no
  unified hero exists; the three above are separate per-paper heroes). Stage 3
  numbers `0 extra exposures` / `≤1.5×` are gated on the C1–C7 campaign (gap G2).
- **Build plan:** new matplotlib figure, 3 horizontal panels in one `figure`
  (gridspec 1×3), schematic vector art (patches/arrows), ≤4 colors from the repo
  palette; export PDF + PNG at 17.8 cm width. New generator
  `code/round63/figs/fig_flagship_hero.py`.

---

## Fig. 2 — Hidden-state ridge (theorem / ridge figure)

**Assigned section:** §2. **Claim:** posterior hidden-state uncertainty subtracts
scene Fisher information and, under dead-time jitter, produces a finite load
optimum with the frozen scaling law. **Make the reader believe the ridge before
seeing an image — no reconstruction thumbnails here** (R35 §3.2).

- **Panel a:** complete versus observed bucket experiment; show the hidden state
  once and the missing-information subtraction once.
- **Panel b:** the actual jitter-capped information curves `J(ρ;c)` for
  deterministic and two jitter levels, with the deterministic and jittered optima
  marked.
- **Panel c:** blind verification of the two frozen scaling exponents on log–log
  axes, with theory lines and independent simulation points.

- **Source materials available now:** `paper2/figs/fig_mechanism.pdf` panel (c)
  = `J(ρ̄)` det-vs-jitter at ν=2000 (→ panel b almost as-is);
  `paper2/figs/fig_jitter_validation.pdf` = peak-vs-`c_v` slope −0.658 +
  retention warning (→ panel c); `results/round63_theory/fig_a_ridge_map.pdf` =
  exact count-Fisher vs load with ρ*(ν) peak (supporting panel b). Frozen
  numbers: `ρ_c` peaks {22.297, 10.739, 5.173, …}, pooled slope **−0.658** vs
  −2/3 (materials map §2b S4); ρ*(ν) table {4.53…22.16} (materials map §1a).
- **Awaits:** panel (a) complete-vs-observed schematic is **new** (not in either
  paper). Nothing campaign-gated.
- **Build plan:** matplotlib 1×3; panel a schematic (patches), panels b/c line +
  marker plots reusing the `fig_jitter_validation.py` / `fig_a_ridge_map` data
  loaders; log–log axes on c. Generator: extend `code/round63/figs/fig_paper2_r20.py`
  or new `fig_flagship_ridge.py`.

---

## Fig. 3 — Ridge realized in imaging (positive evidence)

**Assigned section:** §3. **Claim:** the theoretically selected global multiplier
improves the frozen image task and substantially reduces elapsed optical time.
The elapsed-time result is the **visual and rhetorical headline**; `+1.87 dB` is
supporting evidence that the speedup is not obtained by lowering the quality
target (R35 §3.3). **No gallery mosaic.**

- **Left third:** two large reconstruction comparisons **only** — one
  representative positive and one boundary case.
- **Upper right:** 24-scene fixed-dwell `ΔQ` strip/forest plot, showing median,
  lower bound, and the zero line.
- **Lower right:** elapsed-time speed distribution with the `3×` gate, median
  `19.13×`, lower bound `18.33×`, and censoring symbols. One small
  quality-versus-time curve for orientation if space permits.

- **Source materials available now:** `paper2/figs/fig_speed_results.pdf` already
  carries (a) recon comparisons, (b) per-scene fixed-dwell ΔQ, (c) quality-vs-time,
  (d) elapsed speedup — this figure is **close to a direct reuse/recrop**. Frozen
  verdict numbers all present (register §9a): +1.87 dB / LB 1.41 / 19 of 24;
  19.13× / LB 18.33 / 21 of 24.
- **Awaits:** relabel internal tokens (`FAST_RIGHT_CENSORED` → reader language,
  R35 §3.3); drop any gallery panel to satisfy "two reconstructions only".
- **Build plan:** adapt `code/round63/figs/fig_paper2_r20.py`
  (`fig_speed_results`) — collapse to the four required elements, strip audit
  labels, keep censoring symbols. Reuse existing frozen-artifact loaders (R19
  corrected bundle).

---

## Fig. 4 — Count-only closure (boundary-discovery figure)

**Assigned section:** §4. **Claim:** five orthogonal preregistered interventions
bracket the remaining design and inference freedoms and find no material
count-only headroom. **Title concept: "Double-sided closure of the count-only
scene channel."**

- **Panel a — intervention map.** Left-to-right pipeline from pattern design to
  count likelihood with the five intervention points placed on it. **Enclose all
  five inside one box labeled *same integrated-count observation*** (teaches that
  the campaigns are orthogonal probes of one class).
- **Panel b — common materiality plot.** Incremental image-quality gain over the
  global-ridge baseline on one dB axis, with a shaded **1 dB materiality region**.
  Two braces: **design side** `+0.68`, `−8.3`, `+0.04`; **inference side**
  `≤0.2` headroom, `+0.41`. Use a **capped arrow** for the estimator-efficiency
  upper bound (not a pretend observed gain). If `−8.3 dB` compresses the positive
  points, use a small **broken-axis inset**, not a second figure.
- **Right-edge visual:** a single arrow leaving the count-only box **sideways**
  (not upward in scene dB), labeled **second scientific output** — the handoff to
  DLGI.

- **Source materials available now:** the five certificate values are frozen
  (materials map §3a): bridge/geometry `+0.680`; allocator `−8.33`; DOPS schedule
  `+0.039`; estimator efficiency `≤0.2` headroom (~98.9–99.5%); CPL `+0.412`.
  No existing figure — the certificates are report/JSON only
  (`SCHEDULE_PROBE`, `CPL_PROBE` have no figs dir).
- **Awaits:** the schematic itself (**gap G4** — no double-sided-bracket figure
  exists). Nothing campaign-gated; all numbers frozen.
- **Build plan:** new matplotlib `fig_flagship_closure.py`. Panel a = schematic
  (patches/arrows + enclosing box). Panel b = horizontal dumbbell/marker plot on
  a dB axis with `axvspan` for the 1 dB band, `annotate` braces, capped-arrow
  marker for the efficiency bound, `brokenaxes`/inset for −8.3. Export PDF+PNG,
  ≤4 colors.

---

## Fig. 5 — Dual output  **(PASS-BRANCH ONLY)**

**Assigned section:** §5. **Claim:** the scene and medium ledgers obey
canonical-confusion reciprocity, and DLGI extracts `t_c` and CV without an
additional acquisition or scene penalty. **Title concept: "One bucket record,
two model-certified products."**

- **Panel a — comparator geometry.** Three compact columns: dedicated
  monitor/DCS-like medium record (medium product, extra channel); pilot-
  interleaved GI (medium product, scene rows sacrificed); DLGI (scene + medium
  from the same raw record). *Pre-empts the DCS/pilot priority question visually.*
- **Panel b — medium precision.** Forest/scatter of pilot-free-to-pilot RMSE
  ratios for `t_c` and CV across the confirmatory grid, with `1.0` and `1.5`
  lines. Oracle ratios as light secondary markers.
- **Panel c — certification.** Coverage versus cell for projected `t_c` and CV
  intervals, with the `[0.92, 0.98]` band and an inset for interval-width ratio.
  **This panel is MANDATORY; without it the phrase "model-certified" is
  forbidden** (R35 §5.5).
- **Panel d — scene and schedule.** Scene `ΔPSNR` against the byte-identical
  linear comparator, with the `−0.2 dB` noninferiority line; marker shape / small
  inset for paired/random/ordered ranking and reciprocity prediction.
- **Edge annotation:** a narrow strip under panels b–d mapping the fast-drift,
  slow-drift, and weak-CV/noise failure regions. **Do not devote a sixth figure
  to the edge map.**

- **Source materials available now:** `results/round63_next/DUAL_LEDGER_PROBE/figs/fig_dual_ledger.png`
  (feasibility 3-panel: t_c precision vs oracle/pilot; no-trade schedule;
  reciprocity) — **feasibility only, not confirmatory**. Probe numbers exist but
  bar 4 (interval coverage) FAILED in the probe (materials map §4a); the
  publishable panels **require the frozen-Neyman C1–C7 campaign**.
- **Awaits:** the entire confirmatory grid (**gap G2** — 27-cell campaign, ≥5000
  calibration + ≥1000 confirmation records/cell, ≥12 fresh scenes; not yet run).
  Panels b/c/d are all `\CAMPAIGN{}`-gated. Whole figure is contingent.
- **Build plan:** new `fig_flagship_dual.py` once the campaign artifact exists;
  gridspec 2×2 + bottom edge strip; forest plot (b), coverage step/scatter with
  `axhspan` band (c), ΔPSNR scatter with noninferiority line (d). Do **not** plot
  gain-path correlation 0.998 (R35 cut 6).

---

## Cross-check against the R35 five-figure inventory (§6.2)

1. Fig. 1 — Paper-in-one-view: limit → reach → close → reuse. ✔
2. Fig. 2 — Hidden-state ridge: missing-info specialization, jitter cap, blind
   scaling. ✔
3. Fig. 3 — Ridge realized in imaging: +1.87 dB and 19.13× with reconstructions. ✔
4. Fig. 4 — Count-only closure: five interventions as one double-sided boundary
   experiment. ✔
5. Fig. 5 — Dual output: reciprocity, pilot/oracle, calibrated intervals, scene
   noninferiority, edge domain. ✔ (PASS-BRANCH)

**No sixth figure.** Discussion (§6) gets **no new visual** (refers back to Figs.
1 and 5). MOLT, operating maps, and all five full campaign tables are
supplement-only.
