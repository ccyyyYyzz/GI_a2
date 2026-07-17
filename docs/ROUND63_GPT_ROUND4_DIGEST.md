# ROUND63 GPT round-4 digest — the F1-frozen λ_TV/GOF rule + novelty verdict

**Source**: GPT Pro reply 2026-07-18 (~81 min deep run), chat 候选方案评审与建议.
Verbatim re-assembled archive: ROUND63_GPT_ROUND4_REPLY_EXTRACT.md (one ~200-char
DLP-blocked hole, non-load-bearing). This digest is the OPERATIVE spec for F1.

---

## Q1 ruling: "outer AUDIT split + coherent refit bootstrap"

Adopts candidate A's statistical idea, REJECTS our implementation (per-replicate
full K-fold cross-fit on concatenated cross-fitted λ̂). Verdicts on the rest:
our E[D] ≈ M + Σ(μ(λᵢ)−μ(λ̂ᵢ))²/V diagnosis CONFIRMED; B and C alone both
insufficient (r̄/corr misses zero-mean variance-heavy predictive failure; one-SE
is regularization selection, not adequacy); M+c√(2M)+Î has no finite-sample
calibration. Flag renamed **MODEL_FAIL_PREDICTIVE** — semantics: "calibrated
renewal detector model + frozen RQL-TV reconstruction procedure cannot explain
independent audit data". NOT a pure likelihood-class test; at M<n count-only
data cannot strictly separate detector-likelihood misfit from null-space
under-recovery (must be written honestly in the paper).

### Frozen algorithm (all constants final)

1. **Outer split per cell** (image–seed–condition): logical groups = one
   physical pattern (bern50/gam4) or one complementary pair (hadpair, atomic).
   Frozen hash permutation from (cell_id, seed, pattern_kind, 63, 4): first 80%
   groups → DEV, rest 20% → AUDIT. All arms share the identical split. AUDIT
   NEVER touches: λ_max, η selection, initializer, solver tuning. If AUDIT <
   128 logical groups → GOF_STATUS=GOF_UNDERPOWERED, MODEL_FAIL_PREDICTIVE=NA
   (hits S4 small cells, not S2-A).
2. **η\* per arm, DEV only**: grouped K=5 folds inside DEV; grid
   H = {1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1}; λ_max,arm computed
   on DEV ONLY (TV_NULL_REL=1e-3, expand cap 40, log-bisection 26, 60 outer
   iters per fit — full-data λ_max is a leakage path, must move). Arms fit
   their OWN fidelity, all scored by the calibrated NP renewal residual.
   Per-fold MEAN d_k(η) = (1/n_k)Σ_{i∈k} r_i², d̄(η) = mean_k d_k,
   SE_min = sd_k[d_k(η_min)]/√K (call it "CV fold-dispersion heuristic",
   never a confidence interval). Tie rule: η_min = SMALLEST argmin of d̄.
   **η\* = max{η ∈ H : d̄(η) ≤ d̄(η_min) + SE_min}. GOF is fully OUT of the
   acceptance set** (D2's "AND GOF pass" clause is rescinded).
3. **Coherent RQL plugin scene — once per cell, RQL only** (not per baseline
   arm): λ_plugin = η_min,RQL · λ_max,RQL,DEV (η_min NOT η\*: MODEL_FAIL must
   not be triggered by one-SE's deliberate smoothing bias). Fit on DEV from the
   frozen initializer, n_audit = 25 outer iterations → x̂_plugin. Predict the
   untouched AUDIT: λ̂ᵢ = Φ aᵢᵀ x̂_plugin + d̂. Compute D_obs = Σ_AUDIT r²,
   r̄_obs, c_obs = corr(r, ρ̂). KEY: the bootstrap generator must be A·x̂ + d
   from ONE coherent image (concatenated cross-fitted λ̂ never replays the
   inverse problem).
4. **Refit-per-replicate bootstrap**: B = 39, n_audit = 25. For b = 1..39:
   simulate N*(b) ~ exact-renewal(Φ aᵢᵀ x̂_plugin + d̂, T, τ̂) on ALL DEV+AUDIT
   patterns; recompute the data-dependent initializer on synthetic DEV only;
   SAME numeric λ_plugin (no η re-selection, no λ_max recompute); 25 outer
   iters; D*(b) on synthetic AUDIT.
5. **Exact MC p-value**: p = (1 + #{b : D*(b) ≥ D_obs}) / (B+1);
   **MODEL_FAIL_PREDICTIVE ⇔ p ≤ 0.025** (B=39 → exact 1/40 level). Early
   stop: generate 19 first; if none exceeds D_obs the decision may already be
   determined; else run the remaining 20 — binary-equivalent to full B=39.
6. **Fixed-λ̂ bootstrap demoted to lower-tail leakage diagnostic only**:
   B₀ = 199 exact-renewal sims at the AUDIT λ̂ (no refits);
   p_low = (1 + #{D₀ᵦ ≤ D_obs})/200; LEAKAGE_SUSPECT ⇔ p_low ≤ 0.01. Never an
   upper gate, never changes η\*.
7. **Residual-structure warnings (not in the MODEL_FAIL union)**:
   MEAN_RESIDUAL_WARN ⇔ |r̄_obs| > max_b |r̄*(b)|;
   LOAD_CORR_WARN ⇔ |c_obs| > max_b |c*(b)| (NA if undefined). Paper still
   shows residual-vs-load plots; warnings explain failure mechanisms.
8. **Final reconstruction**: λ_TV,\* = η\* · λ_max,arm,DEV, refit on ALL frames
   (DEV+AUDIT), n_final = 200. MODEL_FAIL changes NOTHING: no η fallback (the
   current E=∅→η_min fallback MUST BE DELETED — GOF must never alter the
   estimator), no cell deletion, no PSNR rescue, no output suppression; it is a
   preregistered adequacy flag in the summary tables.
9. **Cost**: worst case (1+39)×25 = 1000 auxiliary outer iterations per cell
   (≈ 5 full 200-iter fits) + B₀=199 cheap sims.
10. **Continuous-start S3 cells**: no valid independent AUDIT exists
    (afterpulse has unbounded exponential tails). Freeze:
    η\* inherited from the corresponding active-start cell (same arm, image,
    seed, ρ, ν, A); GOF_STATUS=GOF_NA_DEPENDENT; MODEL_FAIL_PREDICTIVE=NA.
11. **RNG keys must be integers**: current rng_for(..., int(round(ctx.T)), ...)
    collapses to 0 for physical T (µs) → different ν share bootstrap streams
    (REAL BUG). Freeze key = (cell_id_hash, nu_integer, tau_ns_integer, ...).

### select_eta.py embarrassment list (all must be fixed at F1)
- B=30 + empirical quantiles / mean±2.5sd has NO level control.
- Refit generator = concatenated cross-fitted λ̂, not coherent A·x̂+d.
- Fixed-λ̂ mean-residual band used as hard gate = inconsistent null.
- η_min band applied to every η (different smoothing ⇒ different null).
- GOF inside the one-SE acceptance set, and failure changes η\* — decouple.
- λ_max computed on full data (AUDIT leaks through the penalty scale).
- int(round(ctx.T)) RNG key bug (above).

## Q2 verdict: Grönberg 2018 + full prior-art sweep

- **Grönberg–Danielsson–Sjölin 2018** (+ SPIE companion + KTH thesis record +
  official Wiley reference list): NO ρ\*(ν)~Cν^{1/3}, NO (6ν)^{1/3}−2/3, NO
  I_N = E[N]−ρ²E[Var(R_ν|N)], NO (ρ,ν) principal ridge of integrated-count FI.
  Their chain: finite-pulse detector physics → analytical count distribution →
  Gaussian approximation validation → spectral-CT task Fisher metrics.
  Companion stops at large-T renewal mean/variance + pileup corrections.
- **Wang 2011** (most dangerous neighbor): CRLB as function of input count
  rate, pileup degrades task performance, ideal nonparalyzable delta-pulse
  analytic arm — but optimizes energy-discrimination imaging, no finite-window
  count-only FI about log λ, no ν^{1/3}.
- **Alvarez 2014**: multivariate Gaussian count CRLB (neighbor of our CLT
  proxy only, no exact-count fallback layer).
- **Bécares–Blázquez 2012** (nuclear): explicit optimal counting rate EXISTS
  (≈1.05e5 cps at representative conditions) — but for dead-time-corrected
  rate relative uncertainty incl. τ-calibration error; not an exact count-only
  Fisher ridge, no ν^{1/3}.
- **Gupta et al. photon-flooded SPAD LiDAR**: closed-form optimal flux from
  ambient background + pileup distortion — different mechanism, not terminal
  phase. **Rapp high-flux lidar**: full detection-time-sequence Markov model —
  sidesteps our count-only bottleneck entirely. **Jorgensen–Johnson 2026**:
  asymptotics for periodic gated event sequences, not single-window integrated
  count. **2025 dToF CRLB works**: optimal flux/pulse width/quantization under
  dead time — the broad "optimal flux" story is CROWDED. **Flow cytometry /
  FFS**: engineering coincidence rules, PCH narrowing — no isomorphic ridge.

**Final novelty verdict**: "No prior cube-root finite-window count-information
ridge was found" SURVIVES, scoped to: ideal nonparalyzable + active-start +
scalar integrated count + exact finite-window FI about log-rate + principal
ridge + (ρ,ν) asymptotics. CANNOT claim broadly that optimal dead-time flux is
new. Safe paper paragraph (use this framing, no "first"):

> "Previous studies established count-rate-dependent CRLB degradation,
> task-specific optimal operating fluxes, and high-flux recovery using
> time-resolved dead-time models. Here, we study a different information
> bottleneck: a finite exposure from which only the integrated nonparalyzable
> count is retained. For this count-only observation, we derive a principal
> finite-window Fisher-information ridge and show that its optimal load obeys
> a cube-root law, ρ\*(ν) = (6ν)^{1/3} − 2/3 + O(ν^{−1/3}), arising from
> terminal detector-phase information loss."

Residual S5 hard item: one final manual provenance sign-off over the actual
2018 PDF / KTH thesis full text (keyword sweep: Fisher, information, count
rate, maximum, optimal, measurement time, dead time, asymptotic) + Paper I's
Fisher formulas and every figure axis.

## Execution deltas vs current repo (F1 checklist)

1. Rewrite select_eta.py to steps 1–11 (single-file rule; delete gof_mode
   tri-state, delete η_min fallback, DEV-only λ_max, integer RNG keys).
2. campaign.py: drop gof_mode cell key; add audit outputs (MODEL_FAIL_PREDICTIVE,
   GOF_STATUS, p_value, LEAKAGE_SUSPECT, warns) once per cell (RQL);
   per-arm rows keep eta_star.
3. Re-run outcome-blind coverage probe against the NEW rule (null false-alarm
   at p≤0.025 over ≥10 correct-model datasets; detection = paralyzable-truth /
   NP-calibrated audit at ρ=0.6, matching S3 semantics) BEFORE freezing.
4. Spec D2 §4 rewrite (D2.1): new rule text + rescind "GOF pass in E";
   prereg wording for MODEL_FAIL_PREDICTIVE semantics + M<n limitation.
5. S1 pilot: gof machinery per new rule; sizing ≈ 85 fits/arm-image ⇒ full
   pilot ≈ 2.5–3 h with 8 local workers.
