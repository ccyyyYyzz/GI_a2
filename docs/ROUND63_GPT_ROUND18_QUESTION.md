# ROUND63 — GPT Round 18 consultation (refreeze GO request + a decisive pre-freeze DEV finding on the certificate class)

Status: the R17 re-architecture is fully implemented and self-tested.
REFREEZE_READY at commit 609de53: all nine §7 items PASS on the outcome-blind
ledger (`results/round63_m1/FREEZE_CHECKLIST_LEDGER.md`) — retirements, hashed
SCAT32 deployed matrix + three modes, verdict emission on positive AND negative
mocks, D_cert SHA (151,344 load + 133,120 gain atoms admitted on DEV), expanded
toys ≤1.8e-9, 43 shards / 5,880 expected cells (480 cert), pre-freeze
confirmatory-access guard exercised. Paper 2 was restructured to R17 §6 with all
NOT-PERMITTED sentences verified absent (commit 8f2b234). No confirmatory
scene, pre-scan count, or endpoint has been opened.

**But we are NOT requesting a blind GO.** The mandated DEV-only certificate
feasibility test (§7 item 6) surfaced huge bounds (G_full/r = 26.7 / 79.1), and
the follow-up discriminator we ran produced a decisive structural finding that
must be adjudicated before the tag. Evidence: `results/round63_m1/
R18_GAP_PROBE_LOG.txt`, probe code `code/round63/dev_gap_probe.py` (this
commit).

## 1. The finding: the dose-only certificate gap is REAL, proven by primal construction

On the DEV scene (632900 class, single scene/seed — replication on further DEV
scenes is running and will be pushed when done), we constructed concrete
feasible designs inside C_dose (budget + ±5% dose band over D_cert, exactly the
R17 §4.3 class) by multiplicative ascent with dose/budget projection:

| cell | SCAT32 logdet | best feasible probe | primal gap/r (LOWER bound on true gap) | D-eff ratio |
|---|---|---|---|---|
| ν=2000, b=0.60 | 3646.43 | 3801.58 | **0.776** | ≈ 2.17× |
| ν=200, b=0.05 | 2788.83 | 2969.35 | **0.903** | ≈ 2.47× |

Probe designs verified budget- and dose-feasible at every recorded point
(dose dev 0.038 ≤ 0.05; loads 0.300/0.60 and 0.025/0.05 caps).

**Composition of the winning designs — the scientific crux:** zero D_gain mass;
all mass at the LOWEST palette load; only ~half the incident budget spent; the
advantage comes entirely from **scene-adapted geometry selection** among
admitted D_load atoms. The γ-audit (§3 of the probe report) confirms this is
not a normalization bug: γ=5 atoms carry emergent loads up to 13.11, are
charged in full in the budget row, the converged dual has an active budget
multiplier θ = 6.755e4 > 0, and the winners simply never touch γ atoms.

The tightened dual side (full 1024 μ-pixels, MU_CAP 1e9 → MU_CAP_ACTIVE=False,
25 CG rounds) reached G_full/r = 106.2 at ν=200/b=0.05 (CG residual 3.5e4 not
fully converged — inflates the bound; 1133 s wall) and hit an LP numerical
failure (HiGHS Status 15, coefficients ~1e13) at ν=2000/b=0.60. The dual is
loose and numerically fragile at full scale — but this is moot for the
structural question because the PRIMAL bound alone is ~80–90× the 1e-2 gate.

## 2. What this means (our reading)

R17 §4.3 excluded A-risk and spectral from the dual target class ("larger class
⇒ stronger certificate"). The DEV constructions show that this larger class
genuinely contains ~2.2–2.5× better designs than balanced SCAT32 — via
geometry selection at low load, which the ±5% dose band does not forbid
(informative, scene-adapted geometries can still tile the dose uniformly).
Consequently:

1. **`DOSE_SAFE_CERT_PASS` as frozen is structurally unpassable**, and we know
   it pre-launch, on DEV, constructively.
2. The R15↔R17 reconciliation becomes sharp: the observed collapse of the
   adaptive arm was enforced by the **conditioning anchors (A-risk ≤ 1.05×,
   spectral ≥ 0.5× vs the balanced reference)** — not by dose uniformity
   alone. Uniform dose constrains where energy lands, not which measurement
   directions are interrogated; the conditioning anchors are what pin the
   design to the balanced family.
3. This is exactly the escape-hatch R17 §1.3 anticipated, with a twist: the
   escape is through geometry selection, not gain coupling. The collapse
   attribution in the paper's Act III ("uniform-dose safety removes that
   degree of freedom") is falsified as worded and must move to the
   conditioning stack.

## 3. Questions for R18 (operative ruling requested; freeze wording for the winner)

- **Q1 — certified class.** Choose: (i) restore A-risk + spectral into the
  certified class, so the certificate targets the FULL deployed safety stack
  (where SCAT32 is plausibly near-optimal and the R15-observed mechanism is
  exactly what is being certified); (ii) keep the dose-only class and
  preregister the expected certificate failure with §4.5 distribution
  reporting only; or (iii) **both — our recommendation**: the gated
  `DOSE_SAFE_CERT_PASS` is redefined over the full-stack class
  C_dose+cond, AND a no-gate descriptive companion reports the dose-only
  gap with the constructive primal design, turning the unpassable gate into
  the sharper localization claim: *dose uniformity alone leaves ~2×
  constructively-realizable geometry-selection headroom; the conditioning
  anchors that guarantee reconstruction trust are what remove it; balanced
  SCAT32 + global ridge control is then certified near-optimal within the
  deployed safety stack.*
- **Q2 — dual construction for C_dose+cond.** If (i)/(iii): prescribe the
  mathematical form for incorporating the A-risk and spectral constraints
  into the dual bound (additional multiplier blocks vs restricted-class atom
  admission), the required toy verification, and whether the 1e-2 bar
  stands for the full-stack class.
- **Q3 — paper wording.** Which R17 §6.2/§6.3 frozen sentences survive
  verbatim? Note the second permitted sentence ("per-pattern **load
  adaptation** provided no certified material D-information advantage…") is
  arguably still true — the headroom is geometry selection, not load
  adaptation — but the §6.3 final claim's "strict uniform-dose safety
  collapses the tested adaptive design space" attributes the collapse to the
  wrong constraint. Please freeze the corrected headline and the Act III
  beat-2 wording.
- **Q4 — numerics and compute budget.** The fragile corner (HiGHS 1e13
  coefficient range at ν=2000/b=0.60) and the 1133 s/cell tight-dual wall
  (480 cells ⇒ ~150 CPU-h, infeasible) both concern the DOSE-ONLY class;
  under (i)/(iii) the gate-carrying certificate targets the full-stack class
  where convergence should be benign. Confirm: per-cell handling remains
  "certified OR reported-failed with distribution" (§4.5), and state the
  acceptable per-cell compute budget / LP-failure disclosure rule for the
  confirmatory run.
- **Q5 — the positive finding.** May the constructive geometry-selection
  result appear in THIS paper as a labeled DEV-only descriptive analysis
  (Act III companion panel + future-work sentence), without violating R17's
  "no further scene or endpoint redesign follows"? We propose no new imaging
  arm now; the follow-up (certified geometry selection under relaxed
  conditioning) is future work.

Caveat disclosed: single DEV scene/seed for the probe; the mechanism is
constructive and scene-generic in form, and multi-scene DEV replication is in
flight and will be pushed to the repo on completion.

Ruling delivery: as before, please post the R18 ruling as a GitHub issue on
`ccyyyYyzz/GI_a2`, referencing this document at this commit.
