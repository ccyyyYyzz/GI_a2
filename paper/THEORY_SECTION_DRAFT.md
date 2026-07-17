# Theory section draft — "A finite-window count-information ridge for dead-time detectors"

Status: DRAFT (pre-F1). Source of record: docs/ROUND63_GPT_ROUND3_DIGEST.md §2
(cube-root law final form) + archived REPLY. Novelty paragraph BLOCKED on the
Grönberg 2018 full-text check (GPT round 4 Q2) — placeholder marked ⛔.
Discipline: no "first"; ridge stated as derived asymptotics "supported by exact
numerical evaluation"; never say variance is O(1) at the ridge.

---

## X. Information-optimal flux for a dead-time-limited counting window

### X.1 Setup

A non-paralyzable single-photon detector with dead time τ observes a frame of
duration T. Arrivals form a Poisson process of rate λ (pattern-dependent);
after each registered count the detector is blind for τ. With ν = T/τ (the
number of dead-time slots per frame) and ρ = λτ (the per-slot load), the
count N in a frame follows the exact renewal law

  P(N ≥ m) = F_{Γ(m,λ)}(T − (m−1)τ),  m = 1, 2, …, ⌈ν⌉,

(active start convention), equivalently P(N ≥ m) = P(Pois(z_m) ≥ m) with
z_m = λ(T − (m−1)τ). All information statements below are exact finite-ν
statements computed from this PMF; asymptotic forms are derived expansions
whose constants are supported by exact numerical evaluation.

### X.2 The missing-information identity

Let I_N(ρ, ν) denote the Fisher information of the single-frame count N with
respect to log λ, and J = I_N/ν its per-slot normalization (J = 1 is the
saturation ceiling: one nat-unit of log-rate information per dead-time slot).
Writing R_ν ∈ [0, 1) for the residual dead-time phase left at the window end
(the unobserved fraction of a dead-time period truncated by the frame
boundary), the information splits as

  I_N = E[N] − ρ² · E[ Var(R_ν | N) ].            (missing-information identity)

The first term is the complete-data information of the underlying renewal
record (each registered event contributes one unit in log-λ scale); the second
is the price of not observing the terminal phase: the frame boundary censors
the last dead-time period, and the count alone cannot tell how much of the
window was "spent" in that final blind interval. The terminal phase has
P(R_ν = 0) = 1/(1+ρ) and density ρ/(1+ρ) on (0,1), so Var(R_ν) →
ρ(ρ+4)/[12(1+ρ)²] → 1/12 at high load — the censoring price grows as
ρ²/(12ν) per slot.

### X.3 The cube-root ridge

Expanding J at large ρ and large ν,

  J(ρ, ν) = 1 − 1/ρ − ρ²/(12ν) + 1/ρ² − ρ/(6ν) + ρ⁴/(144ν²) + o(ν^{−2/3}).

The leading tradeoff is between the saturation gap 1/ρ (a slot only carries
full rate-information when it is actually occupied by a completed
detection–dead-time cycle) and the censoring price ρ²/(12ν). The marginal
balance 1/ρ² = ρ/(6ν) gives ρ³ = 6ν, i.e. the information-optimal load

  ρ*(ν) = (6ν)^{1/3} − 2/3 + O(ν^{−1/3}),

with peak per-slot information

  J(ρ*) = 1 − (c⁻¹ + c²/12) ν^{−1/3} + o(ν^{−1/3}),  c = 6^{1/3}
        = 1 − 0.8255 ν^{−1/3} + o(ν^{−1/3}).

Exact evaluation of the renewal PMF gives ρ* = 4.53, 6.16, 7.87, 9.99, 13.77,
17.45, 22.16 at ν = 20, 50, 100, 200, 500, 1000, 2000 (artifact:
results/round63_theory/fisher_ridge.json); the closed form (6ν)^{1/3} − 2/3
gives 4.27, 6.03, 7.77, 9.96, 13.76, 17.50, 22.23 — a 5.9% gap at ν = 20
shrinking to sub-percent by ν = 500, consistent with the O(ν^{−1/3}) remainder.

Two structural remarks. (i) The ridge is NOT a count-variance pinning effect:
at ρ = ρ* the count variance Var N = ν^{1/3}/6^{2/3} still diverges with ν;
counts only pin to the ceiling at the much larger scale ρ ≍ √ν. The
information peak is an information-geometric tradeoff, not a saturation
artifact. (ii) Beyond the ridge (ρ ≥ ρ*) the count information genuinely
DECREASES — no estimator operating on counts, including the exact-likelihood
reference, recovers it. High-flux operation is therefore beneficial only up to
a computable boundary, and the boundary is the paper's operating-point map.

### X.4 Deployment zoning and the CLT trust boundary

The Gaussian (CLT) surrogate information satisfies J_exact/J_CLT ≃ 1 − ρ²/(12ν),
so the surrogate stays within 10% of exact only for ρ ≤ ρ_0.9 ≍ √(1.2 ν) — a
√ν scaling, distinct from the ν^{1/3} ridge. We partition the (ρ, ν) plane
into: RQL deployment (per-pattern ρ_95 ≤ 1 — the regime the campaign
demonstrates), a transition band, an exact-reference mode (short-window and
extreme-saturation reference computations with the exact PMF), and the
information-decreasing region ρ ≥ ρ*(ν). Zoning uses the per-pattern 95th
percentile of ρ, not its mean: Bernoulli patterns spread per-frame load, and
the deployment guarantee must hold for the loaded tail, not the average.

### X.5 Relation to prior work  ⛔ BLOCKED — do not freeze until GPT round-4 Q2 returns

Dead-time count laws and moment corrections are classical (Müller 1973; Yu &
Fessler 2000); CLT surrogates for photon-counting CT are standard (Alvarez
2014); high-flux timing acceleration exists in the SPAD-LiDAR line (Rapp &
Goyal 2019/2021 — "faster at high flux" is not by itself new); dead-time MLE
and activation statistics are recent (Jorgensen 2026). ⛔ PENDING: the
Grönberg–Danielsson–Sjölin 2018 equation-by-equation check. Claim to defend
(as currently searched): a finite-window, count-only, active-start
information RIDGE with the ν^{1/3} law and the missing-information mechanism
has no located precedent. Wording: "We derive a finite-window
count-information ridge…" — no priority adjectives.

### Supplement pointer

Half-page supplement: "Why a full heteroscedastic Gaussian likelihood fails at
high load" — the log-det incentive reverses sign at ρ = 1/2, and at ceiling
counts (Nτ ≥ T) the per-frame Gaussian NLL is non-coercive
(variance-collapse: NLL → −∞), which is the rigorous form of the brightness-
inflation pathology; RQL's convex objective removes the log-det term and its
pathology while keeping the renewal mean. Include the failure/repair image
pair and the radial objective curve.
