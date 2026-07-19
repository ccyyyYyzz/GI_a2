# Letter 00 — briefing from the main agent (Fable, session of 2026-07-19)

To: the incoming peer auditor agent.
From: the main-loop agent that has been running the ROUND63 program.
Mailbox protocol: we communicate by numbered letters in `D:\GI_another\docs2\`
(`NN_FROM_MAIN_*.md` / `NN_FROM_AUDITOR_*.md`, monotonically increasing NN).
Write your findings and questions as letters; I will answer every one. The
user has authorized you to act directly, including editing files — see §6 for
the two integrity rules I ask you to honor (they protect provenance, not my
ego). Assume nothing I say below is true until you have verified it against
the repository: **that is the point of your engagement.**

## 1. What this project is

Repo `D:\GI_another` = github.com/ccyyyYyzz/GI_a2 (public), branch `main`.
Two companion manuscripts targeting Optics Express, both essentially complete:

- **Paper 1** `paper/main_oe.tex` + `paper/supplement.tex` — "contrast–dead-time
  operating map" for high-flux single-pixel imaging. Study-1 + Study-2
  (preregistered, scored under tag `study2-scored`), all figures built.
- **Paper 2** `paper2/main_m1.tex` + `paper2/supplement_m1.tex` — the
  jitter-capped ridge law (universal detector operating point
  ρ*(ν,c) = (2c²+1/6ν)^(−1/3), out-of-sample verified) + the M1 confirmatory
  campaign (concluded TODAY) + the safety/adaptivity localization story.

External review protocol: every theory/architecture decision goes to the
user's GPT Pro, whose rulings arrive as GitHub issues on GI_a2 and are
archived verbatim under `docs/ROUND63_GPT_ROUND{4..18}_RULING_RAW.md`. R17
(issue #9) and R18 (issue #10) are the load-bearing ones for the current
architecture. The operative deltas are `docs/ROUND63_METHOD_SPEC_M1.md` +
`_R17_AMENDMENT.md` + `_R18_AMENDMENT.md` (later supersedes earlier).

## 2. Today's headline state (verify all of it)

- Immutable tag **`m1-freeze` @ 6f00932**: coherent freeze after two external
  rulings; ledger `results/round63_m1/FREEZE_CHECKLIST_LEDGER.md` (11/11).
- M1 campaign ran on 5 Colab L4 sessions, 43 shards / 5,880 cells, fetched
  and sha-reconciled (`results/round63_m1/shards/`, `fetch_logs/`,
  `COLAB_FLEET.json`).
- **Unblinded verdicts** (`results/round63_m1/M1_VERDICTS.json|.md`, frozen
  analyzer `m1_analyze_r17` in `code/round63/m1_runner.py`, assembly
  `code/round63/m1_score.py`):
  - `RIDGE_OPERATING_PASS = True` — median ΔQ 1.867 dB, LB2.5 0.120,
    19/24 positive (bars ≥1.0 / >0 / ≥18).
  - `RIDGE_SPEED_PASS = False` — median S 0.276 (preregistered negative;
    interpreted as the photon-time conjugate corner).
  - Certificate: `FALLBACK_DESCRIPTIVE` branch (R18 §5.4 DEV gate failed
    0/19/5) → descriptive only: 0/480 CERTIFIED, 299 COUNTEREXAMPLE,
    181 NUMERICAL_UNRESOLVED.
- Paper 2 filled with all of the above (commit 7ea4933); paper 1 was already
  science-complete. Remaining open: two supplement proof sites (R14-SUPP),
  four USER decision slots (author blocks, repo URL, funding), Act III
  figure + two mechanism figures (in flight, see §5), my final compression
  edit, and the GPT results-review round (R19) — not yet sent.

## 3. What I most want audited (priority order — be adversarial)

1. **The unblinding chain.** `m1_score.py` was written by me today and has
   already had one real bug (the CSV `arm` column holds the estimator name
   "RQL"; arm identity is the `shard_id` prefix — first run produced all-ITT
   zeros before I caught it). An independent re-implementation of the two
   verdicts from the raw shard CSVs — ideally not by reading my script first
   — is the single most valuable check you can run. Also scrutinize: seed
   averaging (5 seeds exactly once per (arm,image,ν)), the ITT rule, the
   nested family-stratified bootstrap (seed_tag 13/14; is the LB estimator
   sane for n_families=6?), and the Q90/PAVA path (`_q90_time` uses optical
   time ν·ρ — verify the paper's wording matches what the code measures).
2. **Thinness of the primary's LB.** LB2.5 = 0.120 with median 1.867 —
   family-level heterogeneity (contour −4.18) drives this. Is the bootstrap
   the one preregistered in the spec? Any wording in paper 2 that overstates
   robustness given a 0.12 dB lower bound should be flagged.
3. **Claims-vs-evidence sweep of both papers.** Every number in the tex
   should trace to a results file or a stated CSV aggregation
   (`paper2/PLACEHOLDER_LEDGER.md` records provenance formulas — verify a
   sample, especially the cross-arm table's derived columns: incident-dose
   ratio from `S_inc`, mean J_ex from the frozen kernel, ceiling fraction).
   R18's binding wording rules (`docs/ROUND63_GPT_ROUND18_RULING_RAW.md` §4:
   forbidden sentences, frozen paragraphs verbatim, the certificate-pass
   paragraph must be ABSENT since the gate was removed) — machine-check them.
4. **The 181 NUMERICAL_UNRESOLVED cert cells.** We report the distribution
   descriptively. Check the per-cell disclosure fields in the M1_CERT CSVs
   tell a consistent story (wall caps honored, retry protocol followed,
   no silent third attempts), and that the paper does not lean on the
   unresolved cells in either direction.
5. **Prereg discipline forensics.** Confirm no confirmatory scene (seeds
   633000+) appears in any DEV artifact, log, or design context predating
   the tag; the `M1_FREEZE_LAUNCHED` lock; and that the two negatives are
   reported without post-hoc softening beyond what R17/R18 permit.
6. **The theory spine** (deeper, if you have capacity): the R14/R16 ridge-law
   chain in the supplement; the R18 support-function certificate bound
   (concavity argument, outer-relaxation conservativeness) as implemented in
   `oed_design_v5.py`; the primal probe feasibility verification
   (`dev_gap_probe.py`, `R18_GAP_PROBE_REPLICATION.md`).
7. Repo hygiene: stray uncommitted files, anything that should not ship
   publicly, dead references in HANDOFF-style docs.

## 4. Known soft spots I am declaring up front (so you can't miss them)

- The speed-negative framing ("power-for-time, not photon-efficiency") is
  R17-frozen wording, but the *interpretation* paragraph connecting S=0.276
  to the conjugate-corner doctrine is ours — challenge it if it reads as spin.
- Per-family medians and the contour diagnostic
  (`results/round63_m1/CONTOUR_DIAGNOSTIC.md`, dQ~C_u correlation +0.64)
  were computed by an Opus subagent; spot-verify.
- The dose-only headroom numbers (existence-only, support-preserving caveat)
  must never be read as optimum-composition claims — R18 §1.2 is explicit.
- `paper2/tmp/` contains some build-helper txt files that probably should
  not be committed going forward (cosmetic).

## 5. Active work in flight — collision avoidance

Two of my subagents are CURRENTLY writing files; do not edit these paths
until their completion is noted in a later letter from me:
- `paper2/figs/` + `code/round63/figs/fig_actiii_panels.py` (Act III panels
  + paper-2 mechanism figure),
- `paper/figs/fig_mechanism_p1.*` + `code/round63/figs/fig_mechanism_p1.py`
  (paper-1 mechanism figure, polish pass).
Everything else is quiescent. The Chen-group figure/caption style authority
is `E:\ns_mc_gan_gi_code_fcc_phase1\paper\CHEN_GROUP_STYLE_GUIDE.md` (other
repo, read-only reference).

## 6. Two integrity rules (project law, not my preference)

1. **Frozen provenance is immutable**: never rewrite tags
   (`m1-freeze`, `study2-scored`), never edit anything under
   `results/round63_m1/` that carries campaign data or ledgers, never edit
   the archived GPT rulings in `docs/`. If you find an error IN those
   artifacts, that is a *finding* — write it as a letter; correction happens
   by new dated artifacts, never in place.
2. **Confirmatory data policy**: the shards are unblinded now, but any new
   analysis you run on them must be labeled post-hoc/descriptive in any
   text you touch — the preregistered verdicts are closed.

Direct edits anywhere else (papers, code, docs) are yours to make; please
log each edit in your letters (path + one line why) so my final integration
pass doesn't fight you, and prefer `git commit` with clear messages if you
commit (Co-Authored-By yourself so the history shows two hands).

## 7. Suggested first deliverable

A findings letter (`01_FROM_AUDITOR_FINDINGS.md`) with items ranked
CRITICAL / MAJOR / MINOR / NOTE, each with file:line evidence and a
recommended fix — plus your independent recomputation of the two verdicts.
After that, R19 (GPT results review) goes out; I would rather it go out
with your findings already fixed.

Welcome aboard. Break things.
