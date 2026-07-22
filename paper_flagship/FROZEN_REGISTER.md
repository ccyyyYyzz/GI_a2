# FROZEN-WORDING REGISTER — flagship manuscript

> **RULE (custody-controlled).** Every string in this register is **frozen**: it
> may be *installed verbatim* but **may NOT be paraphrased, softened, or
> re-numbered** during the writing phase without an explicit unfreeze from the
> operator. These are the load-bearing sentences and numbers the architecture
> ruling (R35) and the prior GPT rounds signed off on. If a frozen string feels
> awkward in draft, the surrounding prose is wrong, not the string. When a
> transition "feels stitched," a frozen beam has been weakened or a
> supplement-only detail has leaked into the main text (R35 §11).

**Sources:** `docs/ROUND63_GPT_ROUND35_RULING_RAW.md` (R35 clauses) and
`docs/FLAGSHIP_MATERIALS_MAP.md` §9 (register carried forward). Provenance
commits are those recorded in the materials map header and §-tables.

---

## A. R35 architectural frozen sentences (new in this register)

### A1. Thesis — PASS branch (R35 §1.1) — *the manuscript contract*
> Hidden dynamics set both the operating limit and the extra scientific output
> of bucket imaging: they impose a finite information ridge, leave no material
> scene-recovery headroom across five preregistered count-only escape tests, and
> imprint a model-certified measurement of the medium on the same unmodified
> record.

*Every main-text paragraph must advance one of its three verbs: **set**,
**close**, **imprint**.* Source: R35 §1.1 (issue #27, `cfdb77b` reference brief).

### A2. Thesis — FAIL branch (R35 §1.2)
> Hidden dynamics set a finite information-optimal ridge in bucket imaging, and
> double-sided preregistered closure tests show that global power control is the
> practical optimum of the declared integrated-count scene-recovery class.

*Do not preserve "second product" language by weakening "certified." Switch
cleanly.* Source: R35 §1.2.

### A3. Introduction opening (R35 §3.1) — installed in `main.tex` §1
> A bucket camera records one number per pattern, but that number is produced by
> a dynamic hidden channel.

Second sentence:
> Detector live time and medium fluctuations are usually treated as distortions
> to be corrected; here they determine both what the scene channel can never
> recover and what the same record can measure instead.

Para-1 landing phrase (must end paragraph 1): **"ghost imaging becomes its own
medium monitor."** Source: R35 §3.1.

### A4. Section transitions — load-bearing beams (R35 §11 + §3.x "transition out")
Each installed verbatim in `main.tex` at the named boundary.

| Beam | String | Location |
|---|---|---|
| Intro → §2 | "We first quantify the information removed by the hidden state before asking how that boundary should be used." | end §1 |
| End of theory | "The hidden-state identity therefore yields an operating point, not an estimator." | end §2 |
| §2 → §3 | "The theorem identifies a finite operating point; the next question is whether one scalar control can realize it in imaging." | end §2 |
| End of M1 | "The global dial works; the remaining question is whether richer control of the same scalar count can materially improve it." | end §3 |
| §3 → §4 | "A global dial reaches the ridge, but that does not yet show that richer count-only design or inference cannot materially surpass it." | end §3 |
| End of closure | "Within the declared count-only class it cannot, which redirects the problem from extracting more scene information to reusing the hidden channel as a second output." | §4 beam (comment) |
| End of DLGI | "The same hidden state that fixes the scene-information boundary can therefore become a model-certified medium measurement without a second acquisition." | end §5 |
| §5 → Discussion | "The hidden channel therefore has a double role: it fixes the scene-information boundary and, when its model is validated, supplies a second measurement without a second acquisition." | Discussion (comment) |

Source: R35 §11, §3.1, §3.2, §3.3, §5.6.

### A5. Count-only closure — three-paragraph scope + pivot (R35 §4.2)
Scope sentence (frozen):
> Together these tests do not prove a universal optimum over every imaginable
> optical measurement; they close the declared integrated-count class from both
> the design and software sides.

Pivot sentence (frozen, immediately after):
> Once the scalar scene channel is closed, the remaining question is not how to
> manufacture another fraction of a decibel from the same count, but whether the
> hidden channel can become a second output.

Editorial reframe (frozen phrase, R35 §4.1): call the five negatives **"a
preregistered closure experiment for the integrated-count scene channel"** —
never "five failed methods." Source: R35 §4.1–§4.2.

### A6. Prior-art placement — DCS priority in three sentences (R35 §5.3)
Drafted in `main.tex` §5 from these frozen beats (order fixed):
1. DCS/DWS already estimate medium dynamics from **dedicated temporal intensity
   fluctuations**.
2. Pilot/interpolated GI monitoring already inserts **reference measurements** to
   estimate gain.
3. The claimed distinction is **simultaneous recovery of an unknown scene and
   quantitative medium hyperparameters from the same ordinary programmable bucket
   record, with no pilot/reference acquisition and with a joint confusion
   certificate**.

Source: R35 §5.3. *Do not postpone this distinction to a defensive Discussion.*

---

## B. Register carried from the materials map §9 (pre-existing frozen strings)

### B1. R19 verdict sentences (paper-2 §5–6, §S5) — install verbatim in §3
Provenance: `results/round63_m1/CORRECTION_2026-07-19/M1_VERDICTS_SPEC_CORRECTED_R19.json`
(R19-corrected).

- **Fixed-dwell:** "ridge control improved median fixed-dwell radiometric quality
  over the matched $0.60$ comparator by **+1.87 dB**, with a family-stratified
  lower bound of **1.41 dB** and **19 of 24** scenes positive."
- **Resource:** "The same operating point used **37.1×** the incident dose and
  **2.6×** the detected counts."
- **Speed:** "the median conservative speed ratio was **19.13×**, the
  family-stratified nested-bootstrap 2.5th-percentile bound was **18.33×**, and
  **21 of 24** scenes had $S_{gate}>1$."
- **Exact frozen values (§S6):** $S_{gate}$ median = **19.127043091646133**,
  LB = **18.328492357080282**.
- **Frame:** "These are fixed-dwell and elapsed-time benefits bought with power,
  not photon-efficiency gains."
- **R19 correction record** (Methods / S4 + Table S-correction): the 3
  discrepancies (axis $\log(\nu\rho)\to\log T_{opt}$; family-label → six-fixed-
  strata; nonstandard isotonic → PAVA+censoring) and "**agreed on all 18 audited
  numerical outputs**." Source-of-record `M1_VERDICTS_SPEC_CORRECTED_R19.json` +
  `ANALYSIS_CORRECTION_DISCLOSURE.md`.

### B2. Paper-1 frozen novelty wording (round-4/5/8, PLACEHOLDER_LEDGER)
- **NO "first" anywhere.** Ridge claim scoped to: ideal nonparalyzable +
  active-start + scalar integrated count + exact finite-window FI + principal
  ridge + $(\rho,\nu)$ asymptotics.
- Intro Para-4 "safe novelty sentence" (verbatim): "…spatial occupancy is varied
  under equal mean detector load and equal incident dose, revealing a
  contrast-controlled boundary between multiplex-limited and photon-limited
  high-flux single-pixel imaging."
- Panel (d) may NOT claim coincidence with $\rho^*$; discussion MUST carry the
  peak-irradiance + bench-validation limitation sentence.
- R12 relation-to-prior-work paragraph (§4.4) frozen verbatim; **Grönberg
  full-text sign-off outstanding** (S5/S10 hard item, gap G5).

### B3. R33 theorem names & forbidden claims
- Theorem name: **Canonical-Confusion Ledger Reciprocity Theorem** (distinct from
  method name **DLGI**).
- MOLT theorem name: **Masked Laplace–Power-Sum Theorem** (NOT
  "elementary-symmetric-polynomial theorem").
- **Forbidden (DLGI):** no "information conservation"; no "free medium
  information"; no "paired trades scene into medium"; no additive conservation
  law. Correct wording: "two projections of the same posterior ambiguity, coupled
  through the posterior cross-covariance $K$."
- **Forbidden novelty "firsts"** (R33 §5 + R31 §6): first info partition /
  nuisance-orthogonal design; first single-pixel/optical medium-dynamics or
  $t_c$ measurement (DCS/DWS own it); first calibration-as-science-product (radio
  selfcal); first telemetry medium ID (AO); first blind joint scene/gain
  recovery.

### B4. R34 authorized DLGI claim (until campaign passes) — install verbatim in §5
> "DLGI extracts a competitive point estimate of medium correlation time and
> fluctuation strength from the same bucket record used for scene reconstruction,
> with no added acquisition and no detected scene penalty; model-calibrated $t_c$
> uncertainty remains under confirmation."

- **Do not use "two certified products" yet.** On pass, restore: "One bucket
  stream, two model-certified products…" ("model-certified" + validated
  scalar-gain domain MANDATORY).

### B5. Act-III master verdict (`docs/SOFTWARE_SATURATION_VERDICT.md`) — verbatim
> "No current software candidate clears the operator's novelty × simplicity ×
> generality × image-positive bar. The software layer for these frozen channels
> is saturated at the requested effect size (≥1 dB deployable image gain)."

Identity explanation (verbatim):
> "Information lost = E[Var(hidden state | record)]; the bucket-count record is
> too thin to identify the hidden state…and five independent mechanisms measured
> the bound as binding."

Bridge closing (verbatim):
> "…the certified global operating law (ridge tracking, +1.87 dB / 19.13×, M1) is
> the practical optimum of this regime, now bracketed from both sides by
> preregistered evidence."

### B6. MOLT authorized/forbidden fence (R32) — if MOLT kept (S9 only)
- **Authorized:** "photon-for-rank trade"; sparse support-targeted masks; second
  sensing operator via microcell occupancy; may add little wall-clock in
  switching-limited systems.
- **Forbidden:** universal/dense null-space jailbreak; exact 5% recovery without
  support/prior certificate; photon efficiency / dose neutrality / zero-cost /
  zero-calibration; "negligible time" on MHz arithmetic alone; hiding
  matched-photon failure.

---

## C. Frozen certificate numbers (closure §4 — drafted, not placeholders)
From materials map §3a; installed in `main.tex` §4 paragraphs 2–3.

| Cert | Lane | Governing number | Bar |
|---|---|---|---|
| 1 | RLMI bridge (geometry) | **+0.680 dB** < 1.0 | ≥1 dB |
| 2 | RLMI allocator | **−8.33 dB** (harm) | no-harm |
| 3 | DOPS schedule | **+0.039 dB** < 1.0 | ≥1 dB |
| 4 | Estimator efficiency | ~98.9–99.5%, **≤0.2 dB** headroom | headroom |
| 5 | CPL-GI likelihood | **+0.412 dB** < 1.0 | ≥1 dB |

Main text rounds to **+0.68 / −8.3 / +0.04 / ≤0.2 / +0.41** (R35 §4.2 wording).

---

## D. Abstract numeric-callout cap (R35 §2.1, §12)
Abstract may cite **at most three** numbers: **+1.87 dB**, **19.13×**, **≤1.5×
pilot**. Coverage, lower bounds, cohort counts, and the five individual negative
values **stay out of the abstract**.
