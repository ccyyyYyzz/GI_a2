# M1 method spec — R18 amendment (OPERATIVE; supersedes conflicting parts of the R17 amendment)

Authority: GPT R18 ruling, GitHub issue #10, archived verbatim at
`docs/ROUND63_GPT_ROUND18_RULING_RAW.md` (posted 2026-07-19T12:00:02Z, audit
target c67ba07). The R18 ruling text is normative; this amendment is the
implementation map. Commit c67ba07 is NOT authorized for m1-freeze; the 9/9
REFREEZE_READY ledger is superseded. The empirical arms/endpoints
(RIDGE_OPERATING_PASS, RIDGE_SPEED_PASS, three SCAT32 modes, bootstrap rules,
power-for-time disclosures) are UNCHANGED.

## 1. Verdict rename and class (ruling §2)

- `DOSE_SAFE_CERT_PASS` deleted everywhere → `FULL_STACK_CERT_PASS`.
- Certified class C_stack (frozen, §2.2): ξ ≥ 0, 1ᵀξ = 1, budget cᵀξ ≤ b,
  ±5% pixel-dose rows U±ξ ≤ 0, A-risk tr{V(ξ)⁻¹} ≤ 1.05·tr(V_fix⁻¹),
  spectral V(ξ) ⪰ 0.5·V_fix. A-risk/spectral are GLOBAL design constraints —
  atom-level admission for them is rejected. Peak/ceiling/support/weight/
  RQL-trust stay atom-admission guards. Class is comparator- and
  subspace-relative; manuscript must say so.
- Gate: 480/480 cells `CERTIFIED` with G_stack/r ≤ 1e-2. Non-rescuing
  structural secondary, as before.
- Dose-only class: no-gate DEV analysis only, fields
  `DOSE_ONLY_PRIMAL_GAP_PER_R`, `DOSE_ONLY_D_EFFICIENCY_LOWER`,
  `DOSE_ONLY_DUAL_UPPER_IF_AVAILABLE`, `DOSE_ONLY_SOLVER_STATUS`. No
  `DOSE_ONLY_*_PASS` verdict may exist.

## 2. Certificate construction (ruling §3 — normative formulas)

- Support-function bound: d_a = M_main·tr(V_d⁻¹H_a);
  d̄_d = tr[V_d⁻¹(V_d − V_pre − ε₀I)]; G_stack = sup_{ξ∈C_stack}(dᵀξ − d̄_d);
  by logdet concavity Eff_D ≥ exp(−G_stack/r).
- Exact conic representation: spectral LMI V(ξ) − 0.5V_fix ⪰ 0; A-risk via
  Schur block [[V(ξ), I],[I, W]] ⪰ 0 with tr W ≤ 1.05·tr(V_fix⁻¹).
- Approved implementation: certified cutting-plane outer approximation —
  spectral eigenvector cuts zᵀ{V(ξ)−0.5V_fix}z ≥ 0; A-risk convex tangent
  cuts h(V_k) − tr{V_k⁻²(V−V_k)} ≤ 1.05·h(V_fix); LP master with
  budget/dose prices + column generation; full-dictionary reduced-
  sensitivity scan each round. Finite cut sets are OUTER relaxations ⇒ a
  small bound is conservative; a large bound is NOT a counterexample.
- Mandatory per-cell primal discriminator (§3.4): a frozen SUPPORT-EXPANDING
  primal search over C_stack (positive init on all admitted columns or
  Frank–Wolfe atom injection — the R18-flagged fix of the support-preserving
  probe). Actual logdet improvement per r > 1e-2 ⇒ `COUNTEREXAMPLE`.

## 3. Toy suite (ruling §3.5 — required pre-freeze)

Four independent finite-dictionary toys: (1) budget+dose active, cond.
inactive; (2) A-risk only active; (3) spectral only active; (4) all four
active, and the all-active toy MUST contain a dose-feasible design rejected
by conditioning. Each: independent primal solve (grid or independent conic
solver); primal–dual agreement ≤ 1e-8; min PSD residual ≥ −1e-9; normalized
complementarity ≤ 1e-8; full-dictionary scan residual ≤ 1e-8.

## 4. Compute + status protocol (ruling §5)

- Per-cell terminal status ∈ {CERTIFIED, COUNTEREXAMPLE,
  NUMERICAL_UNRESOLVED}; no imputation, no deletion, no quality-based retry.
- Budget: 60 s primal screen → 180 s primary dual solve → one deterministic
  rescaled retry 180 s (Cholesky whitening by load-matched SCAT32; scores
  scaled by r; budget by b; dose by mean deployed dose). CERT_CELL_WALL_CAP
  = 420 s. No third attempt, no threshold change.
- Trigger for retry: unknown status, time limit, nonfinite objective, or
  coefficient range > 1e10 after frozen scaling.
- Disclosure fields per cell: solver_status_primary/retry, wall_seconds,
  coefficient_dynamic_range, n_dictionary_scans, n_arisk_cuts,
  n_spectral_cuts, full_dictionary_scan_residual, primal_gap_lower_per_r,
  dual_gap_upper_per_r, min_generalized_eigenvalue, arisk_ratio,
  MU_CAP_ACTIVE. Active dual cap ⇒ certification only after ×100 cap
  stability rerun changes bound ≤ 1e-4·r.

## 5. Pre-freeze DEV gate (ruling §5.4 — FINAL predetermined stop rule)

Run the full-stack protocol on all six DEV family images, seed 0, four
anchors (ν∈{200,2000} × b∈{0.05,0.60}) = 24 DEV cells. The categorical
certificate stays in the campaign iff 24/24 CERTIFIED, no full-stack primal
counterexample > 1e-2, median wall ≤ 120 s, max wall ≤ 420 s. Otherwise:
remove `FULL_STACK_CERT_PASS` before the tag and launch the empirical ridge
campaign with full-stack + dose-only analyses DESCRIPTIVE ONLY. No new
class, no anchor moves, no endpoint edits — either branch launches.

## 6. Dose-only DEV evidence + replication set (ruling §1.2, §6)

- Paper reports the constructive finding as existence only, with the
  support-preserving-initialization caveat, exact DEV image/seed IDs, the
  two (ν,b) cells, feasibility, lower-bound gap, D-eff lower bound, solver
  init/support restriction. Mandatory label: "Development-only constructive
  analysis; not a confirmatory endpoint or population estimate."
- If replication is included: the set is FROZEN BEFORE reading aggregates —
  all six DEV images, declared seeds, same two anchor cells, every result
  reported without selection.
- Approved future-work sentence: ruling §6 verbatim.

## 7. Paper wording (ruling §4 — binding)

- R17 §6.2 sentence 1 amended (full conditioning stack), sentence 2 does
  not survive as written, sentence 3 amended (not an either/or). New
  four-beat Act III (§4.2). New frozen paragraphs (§4.3) with the
  conditional pass/fail branch. Final headline replaced by §4.4 boxed text;
  after a successful certificate only "is tested for" → "is certified
  within 1% local D-efficiency". Figure/caption/handbook edits per §4.5;
  delete every claim that uniform dose itself collapses adaptive
  concentration.

## 8. Refreeze checklist (ruling §7 — selftest must map 1:1, 10 items)

(1) rename wired through spec/analyzer/paper/manifests/ledger; (2) C_stack
exact; (3) four-toy suite passes; (4) full-dictionary scan + residual checks
implemented; (5) 24-cell DEV gate passes OR fallback branch executed;
(6) dose-only finding stored as no-gate DEV evidence with caveat;
(7) paper wording corrected; (8) empirical modes/endpoints unchanged;
(9) 52+972 accounting, dictionaries, SHAs, ridge guards, access locks pass;
(10) no confirmatory data opened.
