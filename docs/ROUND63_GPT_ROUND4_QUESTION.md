# ROUND63 → GPT Pro round 4: one pre-freeze adjudication + one literature hard-item

**Repo**: https://github.com/ccyyyYyzz/GI_a2 (public). You are the standing theory
auditor of this campaign. Context chain (read in this order if you lack context):
`docs/ROUND63_GPT_ROUND3_DIGEST.md` (your own round-3 operational blueprint),
`docs/ROUND63_SPEC_D2.md` (the frozen-candidate spec implementing it),
`code/round63/select_eta.py` (the λ_TV rule you prescribed, implemented literally),
`code/round63/physics.py`, `code/round63/solvers.py`.

Status since round 3: all six hard-gate fixes are DONE and verified (guard sign +
afterpulse carry; exact-PMF via Poisson-logsf + logdiffexp + analytic Ġ_m;
RQL direct convex objective as main arm; real integration test 4/4, it caught the
hadpair meta-key bug you predicted; λ_TV six-step rule implemented; SHA gates in
shard runner). Spec D2 is committed. We are at: S1 pilot → F1 immutable freeze.
**One blocking design question surfaced during implementation — we will not
freeze F1 until you adjudicate it.**

---

## Q1 (BLOCKING): the §4 GOF acceptance band does not cover estimation noise → correct-model MODEL_FAIL misfire

Your round-3 six-step rule, implemented literally: cross-fitted held-out residuals
r_i = (N_i − μ_NP(λ̂_i))/√V_NP(λ̂_i), discrepancy D(η) = Σ r_i², and a GOF gate
from an **exact-renewal parametric bootstrap at the fixed cross-fitted λ̂**
(B=200 re-simulations of N at λ̂, band = [2.5%, 97.5%] of D and mean(r);
E = {η admissible by one-SE AND GOF pass}; E=∅ → MODEL_FAIL).

**Observed defect (outcome-blind — no truth/PSNR consulted):** on a smoke cell
where the data generator IS the assumed model (non-paralyzable renewal, exact
same detector params, 32², M=1500, ρ̄=0.6, RQL arm), the observed D = 1881 while
the bootstrap band tops out at ≈1604 ≈ M + 2.6√(2M). MODEL_FAIL fires for the
CORRECT likelihood. It will fire even harder in the S2-A main regime (M/n = 0.5,
more underdetermined).

**Mechanism (we believe):** cross-fitting makes the held-out estimation error
independent of the held-out counting noise, so

  E[D] ≈ M + Σ_i (μ(λ_i) − μ(λ̂_i))² / V(λ̂_i)

— the second term is the standardized predictive reconstruction error, which is
Θ(M)-scale on underdetermined cells and is NOT present in a bootstrap that
re-simulates counts at a FIXED λ̂ (that band only sees counting noise ~ M ± c√2M).
So the literal band tests "λ̂ ≡ λ", not "the renewal model class is right".

**Candidate repairs (pick one, design a better one, or amend):**

- **(A) Refit-per-replicate bootstrap at η_min** (our lean): simulate N* ~
  renewal(λ̂(η_min)); rerun the SAME K-fold cross-fit on N* at η_min (reduced
  iteration budget); D* from the replicate's own cross-fitted residuals; band =
  quantiles of D*. This makes the null distribution include estimation
  variance + smoothing bias by construction. Cost: B×K fits — feasible at
  B = 25–50 with select_iter ≈ 60 (comparable to the 9-η × K grid itself). The
  cheap fixed-λ̂ B=200 band would be retained only as a LOWER diagnostic
  (D far below it ⇒ leakage/overfit).
  Sub-questions: is B ≈ 30 defensible for a 2.5% tail decision? Or should the
  gate be "D_obs ≤ max(D*_{(B)}, M + c·√(2M) + Î)" with a plug-in inflation
  term to de-noise the small-B quantile?
- **(B) Drop the D band from the gate; keep mean(r) and corr(r, ρ) only** —
  misspecification shows up as systematic mis-centering; reconstruction noise
  is (arguably) mean-zero in r̄ at rate 1/√M. Risk: TV smoothing bias is signed
  per-frame on structured λ-fields; is r̄ really centered under the null?
- **(C) Morozov-direct**: replace the absolute band by the relative one-SE rule
  alone (selection unaffected), and let MODEL_FAIL be triggered only by gross
  detectors (r̄ outside ±k/√M, corr(r,ρ) beyond a bootstrap band). Weakest gate,
  simplest to freeze honestly.

Constraints for your ruling: rule must stay **truth-free / deployment-legal**
(never sees x_true or PSNR); computable within ~5× a single arm fit per cell;
frozen BEFORE F1 with exact constants (K, B, band levels, iteration budgets);
and the MODEL_FAIL semantics of your round-3 digest must survive in some form
(we need an honest escape hatch that does not silently degrade to PSNR-picking).
Please give the FINAL frozen algorithm as numbered steps with all constants.

## Q2 (novelty hard-item, do it now in parallel): Grönberg–Danielsson–Sjölin 2018 full text

Your round-3 verdict left one open provenance check on the cube-root ridge law
ρ*(ν) = (6ν)^{1/3} − 2/3 + O(ν^{−1/3}) (missing-information identity
I_N = E[N] − ρ²·E[Var(R_ν|N)], marginal balance ρ³ = 6ν). Please now check the
FULL TEXT of Grönberg, Danielsson & Sjölin (2018) — "Count statistics of
nonparalyzable photon-counting detectors with nonzero dead time" (and any
companion papers by that group on dead-time information/CRLB, e.g. their
spectral-CT detector line) — equation by equation, and report:
1. Do they (or anyone they cite) derive a finite-window Fisher-information
   ridge in (ρ, ν), any ν^{1/3} scaling, or any statement that count information
   peaks at finite flux and then DECREASES?
2. If they stop at asymptotic CRLB / large-T variance formulas, say exactly
   where their analysis ends relative to ours (finite-ν exact PMF Fisher, ridge
   asymptotics, the missing-information identity).
3. Any other prior art you can find for "information-optimal flux / count-rate
   for dead-time detectors" (nuclear instrumentation, SPAD LiDAR, flow
   cytometry included). We will not write the novelty paragraph until this
   comes back; "no prior ridge law found" must survive your best attempt to
   kill it.

## Reply format

Two clearly separated sections (Q1 ruling with final frozen algorithm; Q2
literature verdict with citations). Flag anything else you see in
`select_eta.py` that would embarrass us at review. Do NOT re-litigate settled
round-3 items unless you found an actual error.
