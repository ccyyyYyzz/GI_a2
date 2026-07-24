# BUILD STATUS — Paper 2 (Optica): *Engineering What an Optical Sensor Cannot See*

**Program:** ROUND63 · Certified-Blindness / Three-Ledger line.
**Binding spec:** `docs/ROUND63_GPT_ROUND45_RULING_RAW.md` **Part C** (architecture FROZEN).
**Gate:** `results/round63_next/TLSG/` → **TLSG_PASS** (all 7 bars) → §C architecture frozen.
**Repo HEAD at build:** `35d858e` (branch `main`).
**Working dir:** `paper2_optica/` — nothing here touches `paper_prl/` (paper 1 frozen) or sealed artifacts.

---

## Phase log

### Phase 1 — claim matrix + hero figure (THIS PHASE) — COMPLETE

| Deliverable | Status | File |
|---|---|---|
| Claim–source matrix (every quoted number → committed artifact → hash) | **DONE** | `CLAIM_SOURCE_MATRIX.md` |
| Frozen-source extraction (verbatim committed blobs) | **DONE** | `figures/_frozen_sources/*.committed.json` (11 files) |
| Figure style module (Optica; reuse of paper-1 `prl_style.py`) | **DONE** | `figures/o2_style.py` |
| Fig 1 hero "What the sensor cannot see" (3 linked panels) | **DONE** | `figures/fig1_hero.{py,pdf,png}` |
| Fig 2 jet transmutation (spec + built) | **DONE** | `figures/fig2_jet_transmutation.{py,pdf,png}` |
| Fig 3 three routes + calibration transfer (spec + built) | **DONE** | `figures/fig3_routes_transfer.{py,pdf,png}` |
| Build status + Optica limits | **DONE** | this file |

**Discipline held:** all figures read ONLY committed frozen JSONs extracted verbatim to
`figures/_frozen_sources/` (same protocol as paper 1). **No data generation.** All cited
artifact directories verified clean in the working tree.

### Phase 2 — prose (NOT STARTED)
Sections per R45 §C3 (below). Theorem statements transcribed from R45 §B (analytical, no data).

---

## Frozen source provenance (see CLAIM_SOURCE_MATRIX.md for the full map)

| Source | Commit | Role |
|---|---|---|
| `TLSG/TLSG_RESULTS.json` | `35d858e` | three-ledger scaling + calibration transfer (gate) |
| `R44_KILL_TESTS/KT1` | `19338a4` | support/field wall 6.6e-16 + breaker |
| `R44_KILL_TESTS/IT1` | `f13a323` | energy/unitary wall 2.18e-16 + NA-clip ε¹ |
| `R44_KILL_TESTS/KT2` | `05f1029` | jet transmutation (slopes 2→4) |
| `R44_KILL_TESTS/KT3` | `2430166` | statistical projection (overlaps→0, 86–99%) |
| `R44_KILL_TESTS/IT4,IT4b` | `33b8ef2`,`5666d1e` | blindness-capacity spectrum (80/35/24 dims) |
| `R44_KILL_TESTS/IT5` | `c8a0487` | DPSS route 2.4e-4, noncomposability 80→24 |
| `R44_KILL_TESTS/KT6` | `635ff05` | third-ledger boundary (oracle z⁴, blind killed) |
| `R44_KILL_TESTS/IT2,IT3,IT8` | `acfcda0`,`6fb5149`,`6fb5149` | leak-law + PWM apparatus constraint |
| `R44_KILL_TESTS/IT6,IT7` | `dade877`,`e1546ac` | supplement negatives |
| `WAVE_TWIN/*` | `0c3b0d9` | COMSOL 15.5%, M^1.8, matched-cell 768∈[513,2996] |
| `JET_TEST/`, `SCRAMBLE_EXT/` | `1bf29f1`,`ed7a1e0` | cross-paper jet/scrambling anchors (paper-1 frozen) |

---

## Section architecture (R45 §C3) — frozen

1. **Symmetry walls of an optical measurement** — Statistical Noether Wall theorem + two exact walls (§1 / matrix A,B).
2. **Measuring blindness capacity** — Blindness Capacity theorem, dimension–leak frontier, robustness (§2 / matrix C).
3. **Restoring the desired jet** — optical / code-space / statistical routes + noncomposability (§3 / matrix D).
4. **What finite data can certify** — Three-Ledger theorem + KT6 boundary (§4 / matrix E).
5. **Wave-optics transfer and apparatus constraints** — leak law, DPSS/SVD choice, single bucket, PWM (§5 / matrix F).

**No chronology of the 17 tests appears in the article** (R45 §C3). Supplement carries IT6/IT7
negatives, DPSS-vs-SVD noncomposability spectra, étendue and coherence-rank future objects (§C2).

### Frozen theorem set (R45 §C2)
1. Statistical Noether Wall Theorem · 2. Blindness Capacity Theorem (`β_d=σ_d`; `+δ` robustness) ·
3. Three-Ledger Certifiability Theorem (`1:√p:p`) · corollary 4. Jet Transmutation.

### Hero figure self-critique (R45 §C3)
The hero delivers the story in ~10 s to a non-specialist optics reader: a left-to-right causal
spine (**symmetry → capacity → certifiability**) over three quantitative panels — two machine-zero
walls with their breakers; a measured dimension–leak capacity curve (80 dims ≤1e-5, does not
compose); and the `1:√p:p` ladder with KT6 below the blind line and the projected jet above its
channel line. One-line takeaways sit under the two data panels. Verdict: **PASS** (see final report).

---

## Budgets & FORMAT LIMITS (Optica — verified against author resources, 2026-07)

**Venue:** Optica (flagship). Fallback: Physical Review Applied (R45 §C4).

| Item | Optica limit | Source |
|---|---|---|
| **Abstract** | **~100 words** (explicit summary: problem, methods, major results/conclusions) | Optica Publishing Group Journal Style Guide |
| **Research Article length** | **6–8 published pages** | Optica submission guidelines |
| Letter (Optica Express-style short form) | 4 pages | Optica submission guidelines |
| Mini-review | 8–12 pages | Optica submission guidelines |
| Length metric | published pages (figures/tables count toward the page total); LaTeX via `opticameeting`/`opticajnl` class on Overleaf | Optica author portal |
| Figures | vector PDF/EPS or ≥300–600 dpi raster; colorblind-safe encouraged | Optica style guide |

**Corrected assumption:** the task asked whether Optica's abstract limit is ≤35 words — it is **not**;
the Optica abstract target is **~100 words**. The ~10–12 pp main-text target in the task **exceeds** the
Optica **Research Article** guideline of **6–8 pp** (10–12 pp is mini-review length). This is a
**Research Article** (original results), so the recommendation is to target **~8 pp main** with the
supplement carrying the negatives and extra spectra; a 10–12 pp submission is possible but invites an
over-length flag. Logged as **GAP-4** in the matrix (length-budget decision for Phase 2, not an
evidence gap).

Sources:
- [Optica Publishing Group Journal Style Guide](https://opg.optica.org/content/author/portal/item/style-optica-styleguide/)
- [Author Resources: Journal Style Guide](https://opg.optica.org/submit/style/style_traditional_journals.cfm)
- [Author Resources (Optica review)](https://opg.optica.org/submit/review/optica.cfm)

---

## Title candidates (R45 §C4)
Preferred: **Engineering What an Optical Sensor Cannot See.** Alternatives: *Certified Blindness in
Computational Optics*; *Symmetry Walls and Certifiability in Optical Sensing*; *The Blindness Capacity
of an Optical Measurement*; *From Exact Optical Walls to Certifiable Change Detection*.

---

## Open items (for Phase 2 / user)
- **[GAP-4]** Length-budget decision: Optica Research Article 6–8 pp vs task's 10–12 pp target (recommend ~8 pp + supplement).
- **[GAP-1]** Hero cites the committed KT1 field wall **6.6e-16** (R45 §C3 rounds to 7e-16) — confirm caption wording.
- **[GAP-2]** Caption precision: joint-35 (spectrum) vs fresh-draw-80 (single plane) provenance kept distinct.
- **[GAP-3]** Two committed "→24" noncomposability instances (IT5 SVD→DPSS; IT4b Tukey) — do not conflate.
- **USER slots** (deferred, per program convention): author block, affiliation, repo URL wording, funding/URIS acknowledgment, ORCID.
- Theorem-proof transcription from R45 §B1/§B4/§B11 into main text + supplement (analytical; no data).
- Abstract draft (≤~100 words) + cover letter (Phase 2).
