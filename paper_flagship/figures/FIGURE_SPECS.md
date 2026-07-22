# FLAGSHIP FIGURE SPECS — R36 four-figure architecture

**Source of record:** `docs/ROUND63_GPT_ROUND36_RULING_RAW.md` §3 (revised figure
architecture) + §1 page-and-claim table. **R36 overrides R35.** The R35
five-figure architecture is preserved at the bottom of this file under
**ARCHIVE** for provenance only — do not build from it.

**HARD RULES (R36 exec + §3):**
- **Exactly four main figures. No fifth.** Section 3 (DECODE) gets **no new
  figure** — Fig. 1 already carries the acquisition architecture.
- Each figure earns one claim.
- **BINDING CUT (R36 §3):** the `+1.87 dB` and `19.13×` callouts and every
  M1/ridge/operating-map number **do not appear in the hero, any main figure, the
  abstract, the conclusion, or the title.** They live in Supplement S4 only.
- Caption 120–160 words; no provenance narrative in captions.
- Build in matplotlib per repo conventions (`code/round63/figs/*.py`).

**Repo figure style** (`paper/figs/FIG_MECHANISM_P1_README.md`,
`paper2/figs/ACTIII_README.md`): sans-serif; neutral **gray = safe**;
**blue = ridge/high-flux**; **vermillion = negative**; hero uses **≤ 4 colors**.
Sizes ≈ 8.6 cm / 17.8 cm; **Optica resize/font pass is gap G6 (pending)**.

**Pre-verdict wording (R36 §3):** use **"two inferred products"** in Fig. 1
before C1–C7 pass; **"model-certified"** only after the campaign clears.

**Fail branch (R36 §7):** a protagonist switch → a **≤3-figure law/boundary
paper** with a different figure set (hidden-state ridge; M1 evidence; double-sided
closure). That set is specified in `../main_failbranch_notes.tex`, not here.

---

## Fig. 1 — HERO: the hidden channel becomes the instrument

**Assigned section:** §1 Introduction. **Claim:** the same hidden medium state
both limits scene inference and leaves a measurable temporal imprint in the
ordinary bucket stream. **One horizontal figure, three visually continuous zones —
not a dashboard** (R36 §3).

### Left — one unchanged acquisition
`patterns → scene x → dynamic medium a_t → bucket stream Y`. Show **one** bucket
detector and **one** time series. Single footer:
> same photons · same exposures · no pilot · no reference arm

### Center — the limit
The hidden path is not observed directly. Show one small **posterior-ambiguity
cloud** around `a_t` and one minimal statement:
> hidden-state uncertainty subtracts scene information

**At most one compact equation.** Do **not** display the dead-time ridge, the
operating map, or any campaign number.

### Right — the reuse
The **same raw `Y`** bifurcates into `scene image x̂` and `medium (t̂_c, ĈV)`.
Connect the two output boxes with a **thin arc labeled "shared canonical
confusion"** — not "conservation" or "information transfer". Visual conclusion:
> one record → two inferred products *(→ "two model-certified products" after C1–C7)*

### Binding cut
**The `+1.87 dB` and `19.13×` callout goes entirely** — not in the hero, any main
figure, abstract, conclusion, or title. The dead-time sibling remains fully in
Supplement S4.

- **Source materials now:** `paper/figs/fig_mechanism_p1.pdf` (acquisition-chain
  seed, left zone); `results/round63_next/DUAL_LEDGER_PROBE/figs/fig_dual_ledger.png`
  (dual-output seed, right zone). Center posterior-ambiguity cloud is **new**.
- **Awaits:** the unified three-zone composition (**gap G3**); the "two
  model-certified products" wording is gated on C1–C7 (**gap G2**).
- **Build plan:** new `code/round63/figs/fig_flagship_hero.py`; gridspec 1×3 with
  continuous flow (shared baseline, connecting arrows across zone boundaries),
  schematic vector art, ≤4 colors, 17.8 cm. **No ridge curve, no numbers.**

---

## Fig. 2 — Reciprocity without an algebra wall

**Assigned section:** §2 (DEFINE). **Claim:** reducing canonical confusion
protects both products. **Three compact elements** (R36 §3); **do not show the
full proof or every matrix identity.**

1. the joint Fisher block `[[A, C], [Cᵀ, B]]`;
2. whitening to `K` and the **shared singular-value loss spectrum** `κ_i` for the
   two Schur complements;
3. **paired / random / ordered** schedules with the preregistered **predicted**
   ordering and the **measured** ordering.

The figure's one sentence:
> Reducing canonical confusion protects both products.

- **Source materials now:** reciprocity structure frozen (R33 Thm 2, materials
  map §4b); schedule evidence `results/round63_next/DUAL_LEDGER_PROBE/t3_reciprocity.json`
  (reciprocity det exact to 1e-15; ‖K‖ schedule-invariant 0.2%; paired best for
  both ledgers). No existing single figure.
- **Awaits:** the schematic composition (**gap G4**-adjacent); the measured
  schedule ordering on the confirmatory grid is **gap G2**-gated (feasibility
  values exist).
- **Build plan:** new `fig_flagship_reciprocity.py`. Element 1 = block-matrix
  schematic (patches). Element 2 = stem/bar plot of `κ_i` with the two Schur
  losses annotated as sharing the spectrum. Element 3 = small paired/random/
  ordered comparison (predicted vs measured markers). Independent of
  reconstruction quality (R36 §5 item 2).

---

## Fig. 3 — Medium product: precision plus calibration

**Assigned section:** §4 (CERTIFY). **Claim (part 1):** the medium product is
competitive and certified. **Two-panel result figure** (R36 §3). **No cell
table.**

- **Left:** pilot-free/pilot **RMSE ratios** for `t_c` and CV across the primary
  grid, with the **`1.5×` line** and an **oracle reference**.
- **Right:** empirical **95% coverage** and **median interval-width ratio**, with
  the **frozen acceptance band**.

Use small multiples or a compact heat map **only if every cell stays readable at
final column width**.

- **Source materials now:** `DUAL_LEDGER_PROBE/figs/fig_dual_ledger.png` panel A
  (t_c precision vs oracle/pilot) is a **feasibility** seed; probe bar-4 coverage
  FAILED (materials map §4a), so the publishable version needs the campaign.
- **Awaits:** the entire confirmatory grid (**gap G2** — 27 cells, ≥5000 calib +
  ≥1000 confirm/cell, ≥12 fresh scenes). All values `\CAMPAIGN{}`.
- **Build plan:** new `fig_flagship_precision.py` once the S8 campaign artifact
  exists; left = forest/scatter of RMSE ratios with `1.5×` line + oracle markers;
  right = coverage points with `axhspan` acceptance band + width-ratio inset. Pilot
  and oracle baselines embedded so the evidence is not self-referential (R36 §5
  item 4).

---

## Fig. 4 — Scene product and validity domain

**Assigned section:** §4 (CERTIFY). **Claim (parts 2–3):** the scene product is
not spent to obtain the medium product, and the theorem predicts the boundary.
**Two-panel result figure** (R36 §3). **No gallery.**

- **Left:** **scene noninferiority** against the byte-identical **plain-linear**
  comparator, with the **`−0.2 dB`** line and the **resource-ledger audit**
  embedded as a small annotation (same photons/exposures/pilots/duration).
- **Right:** the **validated operating domain** in `(Δ/t_c, T/t_c, CV,
  photon-SNR)` with **fast-drift, slow-drift, and weak-fluctuation failure
  regions clearly marked**.

One representative scene pair may be **inset only if** it explains a result the
quantitative panel cannot.

- **Source materials now:** noninferiority + edge data in `DUAL_LEDGER_PROBE/`
  (Δ −0.04 to +1.76 dB, superior 7/9; fast t_c=2 +40%, slow t_c=64 25% floor —
  materials map §4a) are **feasibility**; the certified version needs the campaign.
- **Awaits:** confirmatory grid + edge bank (**gap G2**). Model falsifiability
  (where scalar-gain fails) is part of the claim (R36 §5 item 5).
- **Build plan:** new `fig_flagship_scene.py`; left = ΔPSNR scatter vs linear with
  `−0.2 dB` line + ledger-audit annotation; right = domain map (`pcolormesh`/
  contour) over dimensionless axes with shaded fast/slow/weak failure regions.

---

## Cross-check against the R36 four-figure inventory (§3)

1. Fig. 1 — Hero: the hidden channel becomes the instrument (limit → reuse). ✔
2. Fig. 2 — Reciprocity without an algebra wall. ✔
3. Fig. 3 — Medium product: precision + calibration. ✔
4. Fig. 4 — Scene product + validity domain. ✔

**No fifth figure.** §3 (DECODE) and §5 (Scope) carry **no new visual**. All
ridge/M1 operating maps, five closure tables, and MOLT are supplement-only.

---
---

# ARCHIVE — retired R35 five-figure architecture (do NOT build)

Kept for provenance only. Superseded by the R36 four-figure architecture above.
The R35 spine was LIMIT → REACH → CLOSE → REUSE with five figures; R36 rejected it
as a program summary and cut to one DLGI protagonist.

- **R35 Fig. 1 — Hero "10-second paper":** three-stage limit→reach→close→reuse
  with `+1.87 dB` and `19.13×` callouts and a five-tick count-only boundary wall.
  *Retired:* R36 removes the M1/ridge callouts (binding cut) and re-centres the
  hero on the two-ledger reuse; the boundary is now one frozen text paragraph
  (§2), not a hero stage.
- **R35 Fig. 2 — Hidden-state ridge:** missing-info specialization + jitter cap +
  blind scaling verification. *Retired to Supplement S4* (dead-time sibling); the
  ridge law is no longer a main figure.
- **R35 Fig. 3 — Ridge realized in imaging:** `+1.87 dB` / `19.13×` reconstructions
  + forest plots. *Retired to Supplement S4.*
- **R35 Fig. 4 — Count-only closure:** five interventions as a double-sided
  boundary experiment. *Retired:* the closure is now the frozen four-sentence
  boundary paragraph in the Introduction (main text) + Supplement S6 (full
  evidence); no main figure. **(This figure returns as fail-branch Fig. 3.)**
- **R35 Fig. 5 — Dual output:** reciprocity + pilot/oracle + calibrated intervals
  + scene noninferiority + edge domain (five panels). *Split under R36* into
  **Fig. 3 (precision + calibration)** and **Fig. 4 (scene + validity domain)**,
  with reciprocity/schedule moved to **Fig. 2**.

Mapping R35→R36: R35 Fig 1→R36 Fig 1 (recut); R35 Fig 2/3→S4; R35 Fig 4→text §2 +
S6 (and fail-branch Fig 3); R35 Fig 5→R36 Figs 2+3+4.
