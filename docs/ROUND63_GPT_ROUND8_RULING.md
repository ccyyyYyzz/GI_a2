# GPT round-8 ruling — frozen Study-2 spec (contrast–dead-time phase diagram)

Reconstructed from screenshots (DLP blocked tool extraction); ChatGPT thread
governs on discrepancy. Operative alongside rounds 4–7. Date 2026-07-18.

## Framing (verbatim intent)
Study 2 is a NEW preregistered follow-up, not a reinterpretation or
replacement of Study 1; Study 1 remains the confirmatory dense-pattern
result (24/24 uninformative, S_gate=1, fixed-budget +0.28 dB / 24/24 / LB
+0.23). **New primary hypothesis (boxed):** "At fixed mean detector load and
fixed incident-photon budget, increasing illumination-pattern contrast by
reducing spatial occupancy opens a photon-limited regime in which
dead-time-aware high-flux operation reduces acquisition time."
"Sparse illumination itself is not the claimed invention."

## Frozen constants and protocol
1. **Scene geometry**: 32×32, n=1024, **M=n=1024**; same nine dwell points
   ν∈{5,10,20,50,100,200,500,1000,2000}; ρ₀=0.05, ρ₁=0.60; τ=50 ns, σ_b=0,
   d=0.
2. **Scenes**: the six detail families re-parameterized at 32² (glyph,
   chirp, spokes, maze, contour, band-pass microtexture on frozen nonzero
   background). Confirmatory seeds s_{f,r}=632000+100f+r (f=0..5, r=0..3);
   six dev-only instances 631900+f. **Each generator output must include:
   normalized truth (x≥0, Σx=1); a SIGNAL ROI; a BACKGROUND ROI; all
   generator parameters; image and ROI SHA256.** No image may be rejected
   or regenerated based on TV, Fourier content, estimated contrast, PSNR,
   CNR, speedup, or visual quality.
3. **Pattern ladder**: four physical occupancy levels **k∈{512,32,16,1}**
   (occupancy 50%/3.1%/1.6%/0.1%). Construction (binary matrix B^(k)):
   initialize from k distinct cyclic permutations; apply 20n=20480 accepted
   degree-preserving 2-switches; verify binary entries, every row/column sum
   k, and rank n; if rank-deficient advance a preregistered nonce and repeat
   solely until full rank; record nonce and condition number but impose NO
   condition-number acceptance threshold. k=1 = seeded permutation matrix
   (raster).
4. **CRITICAL photon-budget normalization (the main amendment)** — sparse
   patterns must NOT win merely because they reduce detector count rate:
   on-pixel incident rate **Φ_on(k)=ρ̄n/(kτ)**; equivalently the
   reconstruction code may use **A^(k)=(n/k)·B^(k), Φ=ρ̄/τ**; then for a
   sum-normalized image E_i[λ_iτ]=ρ̄ and ALL occupancies share identical
   mean detector load, incident photons per pattern, and cumulative
   incident dose per object pixel.
5. **Primary confirmatory arm: ONLY k=16 carries the positive gate.** Run
   24 images × 5 measurement seeds × both ρ × all nine ν × M=1024. Primary
   reconstruction = RQL with ALL Study-1 machinery unchanged: c=0.50,
   analytic score-noise TV scale, actual κ_A (≈n/k here), pooled-flux
   initialization, 200 solver iterations, same radiometric metrics. On the
   full k=16 grid also run POISSON-LIN, SAT-POISSON, PRECORRECT, GI
   (display/reference arm).
6. **Primary gate transfers unchanged**: median S_gate ≥ 3; LB2.5>1;
   #{S_gate,j>1} ≥ 18/24; 10,000 nested replicates (seed-resample within
   image → rebuild isotonic curves + censoring → family-resample 4×6 →
   24-image median). Nominal homogeneous-load Fisher-rate ratio
   0.6/(1.6)/(0.05/1.05)=7.875 → 3× median stays a reasonable lower
   engineering threshold.
7. **Fixed-budget secondary unchanged** (ΔQ_j at ν=2000; median ≥1 dB,
   LB>0, ≥18/24; cannot rescue a failed primary).
8. **k=32 robustness arm**: same RQL grid, all 24 images, 5 seeds;
   preregistered robustness/replication arm with NO confirmatory gate.
9. **Dense (k=512) and raster (k=1) controls**: demonstration/controls (no
   gate). Note: "positive Hadamard/cosine multiplexing can be inferior to
   point scanning under photon noise" is an ESTABLISHED effect (cited), not
   a Study-2 discovery.
10. **Mechanism variables (descriptive, never gates)**: do NOT claim
    contrast=1/√k automatically (scene-dependent). Record actual
    bucket-energy contrast C_u=sd_i(u_i)/mean_i(u_i); ideal fixed-k
    Var(u)=((n−k)/(k(n−1)))(n‖x‖²−1); first-order contrast-to-count-noise
    index **Γ(ρ,ν)=C_u·√(νρ/(1+ρ))** (Γ≪1 buried; Γ≳1 photon-limited
    structural regime observable); also the pattern-aware ideal ρ^pattern.
11. **CNR co-secondary (no gate)**: their exact definition with
    generator-frozen signal/background ROIs; record NA if denominator zero
    (no tuned floor). Report two photon axes: S_det=(1/Mn)ΣN_i (Hao's
    measured-photons convention, 0.096–0.44 in their experiments) and
    S_inc=(1/Mn)Σλ_iT (required: dead time makes S_det nonlinear).
12. **Descriptive flux map (no gate)**: one frozen representative per
    family; ρ∈{0.05,0.1,0.3,0.6,1,2}, ν∈{20,200,2000}, 3 seeds, all four k;
    report PSNR, CNR, C_u, Γ, ceiling-count fraction, detected/incident
    photon rates. Do NOT preregister that image CNR peaks exactly at the
    scalar Fisher ridge (ridge is a scalar count-information result; the
    spatial image optimum is an outcome).
13. **Mechanism-separation run (MANDATORY, closest-prior killer)**: on the
    six mechanism representatives, ν∈{20,200,2000}, 3 seeds — separate
    (1) known sparse mitigation: turn off pixels and thereby REDUCE detector
    rate (unnormalized flux) vs (2) the Study-2 mechanism: retain the same
    mean load and photon budget while using sparsity to raise bucket
    contrast. Prevents "this is just the known consequence of lowering the
    count rate".
14. **Dynamic-scattering secondary (G5 approved NO-GATE only)**: primary
    k=16 geometry on six frozen family representatives; ρ∈{0.05,0.6},
    ν∈{20,200,2000}, 3 seeds; multiplicative scaling factors applied BEFORE
    dead time: α_t=exp(z_t−σ_z²/2), z_t=0.95·z_{t−1}+ε_t, stationary
    CV(α)=0.20. Collaboration-ready (their scaling-factor problem).
15. **Freeze discipline**: dev instances may not alter k, scene families,
    ρ, ν, c, endpoints or gates. After the Study-2 freeze: **if k=16 fails,
    there is no further geometry, image-class or endpoint redesign**; a
    k=32 or raster success remains secondary.

## Q2 — prior-art ruling (risk REAL and closer than "sparse SPI exists")
- Zhang et al.: pile-up in single-photon NIR SPI (random+Hadamard),
  image quality vs count rate, sparsities 50%→3.1% — "sufficiently sparse
  patterns mitigated pile-up by REDUCING the detector photon rate".
- Liu et al.: sparse measurement matrix + pile-up correction in
  single-pixel compressive LiDAR. Han et al.: tunable sparse Hadamard
  speckle patterns (sparse-Hadamard SPI is an explicit method family).
- NOT defensible: "we introduce sparse illumination to mitigate pile-up";
  "we first study sparse-pattern dead-time SPI"; "we first show count rate
  has an optimum in SPI".
- Defensible-difference table: prior work lowers detector photon count at
  effectively unchanged irradiance to AVOID pile-up, with empirical curves
  and selected acquisitions and no boundary theorem; frozen Study 2 holds
  mean detector load FIXED (on-pixel intensity rises as n/k), tests whether
  contrast makes HIGH load usable, exact renewal model + (ρ,ν) map,
  preregistered time-to-quality endpoint, explicit C_u/Γ/FI ridge/ladder.

## Q3 — paper structure ruling
Two-study structure is STRONGER than an isolated positive. Retitle:
**"When high flux helps single-pixel imaging: a contrast–dead-time phase
diagram"** (or "Pattern contrast governs high-flux single-pixel imaging
with dead-time-limited detectors"). Study 1 appears in the MAIN TEXT
immediately before Study 2 — not supplement, not "pilot that led to the
answer"; its frozen negative and the +0.28 dB uniform positive are part of
the scientific result. Limitations must state: sparse operation requires
much higher peak irradiance and still needs bench validation.
**Main six-panel ladder figure**: (a) exact FI map + cube-root ridge;
(b) representative k=512/32/16/1 masks with equal-dose bookkeeping;
(c) actual C_u and Γ versus k; (d) safe/fast PSNR and CNR versus dwell for
a USAF target; (e) S_gate and fixed-budget gain by occupancy; (f)
dynamic-scattering result in Hao-compatible S_det units. Do not promise
panel (d) peaks precisely at ρ*.
**Final paper-level novelty sentence (safe version, verbatim)**: "Prior
work has used sparse patterns to reduce detector count rates and alleviate
pile-up. Here, instead, spatial occupancy is varied under equal mean
detector load and equal incident dose, revealing a contrast-controlled
boundary between multiplex-limited and photon-limited high-flux
single-pixel imaging."
