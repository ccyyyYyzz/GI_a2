# FOG_QUOTIENT_TEST — the image-erasing Q-only channel through complete scrambling (R46 §3.2–3.3)

**VERDICT: PASS — all four frozen bars pass at full scale.** A completely scrambling medium
erases every spatial symbol yet carries a calibration-free message encoded in a ratio of
quadratic energies: the unknown 40-dB medium amplitude cancels exactly (F-ratio), and a
constant-DC texture alphabet is invisible to the mean channel while covariance carries the
symbol. This survives the R46 §6.3 kill (it has *both* the mean-invariant alphabet *and* exact
nuisance-cancelled capacity), so it is **not** noncoherent amplitude modulation renamed.

Engine: SCRAMBLE_EXT `scramble_toy.py` conventions. Full run `fog_quotient_test.py FULL=1` on
**Colab L4**, 1203 s. Smoke evidence retained in `FOG_QUOTIENT_TEST_SMOKE.json`.

## M_eff diagnostic (coordinator custody item)

The **raw O∘O Hadamard-square participation ratio = 1.47** — this is the *naive* covariance
projection M_eff (DERIVATION.md predicts ~1 for it), and an earlier build that used an unwhitened
energy statistic reported exactly this (M_eff≈1.46). The **Fisher-efficient** `M_eff = tr(K²)`
with `K = Σ₀^{-1/2}(q₀ O∘O) Σ₀^{-1/2}` is **12.15 ≈ 13** (the DERIVATION value); `tr(K) = 16.26 ≈
(k_p+1)² = 16` is the total informative looks. The gap is *naive-projection vs whitened*, not a
scale artifact — the engine config is the intended one. The corrected whitened, shot-offset-
subtracted energy `S_clean/(aQ) ~ χ²_{M_eff·T}` is used throughout.

## Frozen bars (declared before running) — measured at full scale

| Bar | Frozen | Measured | Verdict |
|---|---|---|---|
| 1 log-ratio noise matches log-F law | pooled KS p > 0.05 | **KS p = 0.636, CvM p = 0.76** (ν=777, 200 pooled) | **PASS** |
| 2 ratio BER & MI invariant to unknown 40-dB medium | dMI < 0.05 AND dBER < 0.02 (CRN-paired) + absolute collapses | **dMI = 0.005, dBER = 0.0004**; ratio BER 0.001→0.000, MI 2.316→2.321; absolute BER 0.000→**0.794**, MI 2.322→**0.041** | **PASS** |
| 3 constant-DC texture invisible to mean channel | mean-detector AUC ∈ [0.45,0.55] all ρ (≥4096 rec/ρ) | **AUC = 0.499 / 0.500 / 0.508 / 0.504** (SE 0.006) | **PASS** |
| 4 rates nontrivial under M_eff≈13 bank cost | BA capacity > 0.5 bit some T AND achieved MI > 0.5 bit | **cap(T=16,64,256,1024)=1.54/2.15/2.32/2.32 bit; achieved 5-symbol MI = 2.321 bit** | **PASS** |

Full-scale requirements met (coordinator custody): N_TX = 20,000; AUC pool = 4096 records/ρ;
common-random-number pairing across amplitude conditions; frozen bars unchanged.

## Reading

- **Exact F-ratio channel (BAR1).** With reference and message blocks sharing an unknown medium
  scale `a`, `R = (S₁/ν₁)/(S₀/ν₀) = ρ·F_{ν₁,ν₀}` — `a` cancels exactly. The engine's empirical
  `Z = log R` at ρ=1 matches the exact `log F` law (KS p=0.64), confirming the effective dof is
  `ν = M_eff·T` with the whitened `M_eff≈13`.
- **Calibration-free (BAR2).** Under a 40-dB (10⁴ energy) log-uniform random medium amplitude,
  the exact-ratio decoder is invariant (dMI 0.005, dBER 0.0004, CRN-paired) while the
  absolute-energy decoder — calibrated to a fixed medium — **collapses to near-chance**
  (BER 0.79, MI 0.04). The ratio, not the energy, carries the message.
- **Image-erasing / mean-invariant (BAR3).** The constant-DC texture symbols have exactly equal
  bucket mean (sum x preserved) but different quadratic energy Q; the best a-blind mean detector
  is indistinguishable from chance (AUC 0.50 ± SE across all ρ at 4096 records/ρ). The message
  lives purely in the covariance quotient.
- **Nontrivial rate (BAR4).** Exact F-law transition matrix + Blahut–Arimoto capacity gives
  1.5–2.3 bits per reference/message pair at the measured `M_eff≈13` bank cost; the realized
  5-symbol mutual information is 2.32 bits.

## Consequence

The Fog Quotient Channel is a genuine image-erasing, Q-only covariance communication channel
through complete dynamic scrambling, with exact differential cancellation of the unknown medium
scale and a mean-invariant texture alphabet — the narrow novelty R46 §3.4/§6.3 requires. Not
killed. Bench-phase only; nothing here touches the Letter or sealed artifacts.
