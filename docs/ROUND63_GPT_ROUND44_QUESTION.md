# R44 — Divergence round over the wave-optics twin: the wall leaks, the channels cross, the family beckons

The manuscript is BUILT: R43's green-light protocol executed end-to-end. Phase 2
compiled (main 6pp: 4pp core + End Matter A/B/C; supplement 5pp S1–S5), both
§7.3 reads completed (general-physics: PASS_WITH_FIXES 0/4/6; hostile
statistical-optics: PASS_WITH_FIXES 1/5/7 — sim-labeling, G≻0 qualifier,
"specificity" word, fence gaps, 2-point CUSUM disclosure; 21/22 custody
spot-checks exact, refs 16/16 verified), all adjudicated repairs applied and
committed @ e6213b1. The Letter is submit-ready minus user items (authors/
affiliation/URL) and ONE pending paragraph: the bench-transfer wording — which
is what this round is partly for.

Per the standing quarantine (your R43 §7): the wave twin generates NO manuscript
data. Its findings feed (a) the transfer paragraph (the single allowed
write-back), (b) the bench protocol, (c) the next paper.

## What the twin measured (results/round63_next/WAVE_TWIN/, uncommitted-yet)

Scalar band-limited angular-spectrum twin of the planned bench: DMD (10.8 µm
micromirrors, fill 0.92, 64×64 code on 32 µm macro-pixels, band k_p=5 enforced)
→ thin Gaussian-correlated phase screen (l_c ∈ {5,15,50} µm; developed σφ=2π
and weak 0.3 rad; rotation = lateral shift, quasi-static per bank) → static
64×64 object with a toggleable beyond-band feature (k∈[6,9] annulus, 2%/5%
energy) → single bucket, Poisson 1e4 pe. Energy conservation 1e-6; complex64
validated against float64 to 5 sig figs on the wall metric.

**T1 — the mean wall LEAKS in real optics.** Statistical model: exact zero
(5.7e-16). Physical: signed bucket = |E⁺|²−|E⁻|² with a± = (1±s)/2; after
propagation E is complex and the cancellation breaks; micromirror pixelation
leaks even at contact. Measured relative in-band leak at z1 = 0, 1, 2, 5, 10,
20 mm: 6.8e-3, 1.2e-2, 1.6e-2, 3.2e-2, 5.4e-2, 8.9e-2. Pixelation floor ~0.7%
at conjugate imaging; monotone growth with z1.

**T2 — the (k_w, σ_f) grid axes are physically COUPLED.** Object-plane grain is
illumination-diffraction-limited (grain ≈ λ·z2/D_illum): with the 2 mm
illumination aperture only fine grain (kwf ≳ 4) is reachable — the sealed
grid's coarse-medium cells are unrealizable in this geometry. CONTACT
DEGENERACY: a pure phase screen at z2=0 leaves intensity unchanged — zero
covariance signal at contact; contrast (σ_f) develops only with z2.

**T3 — aperture law CONFIRMED.** Usable covariance scene-frequency edge tracks
grain (27 → 20 → 13 as grain coarsens 8 → 48 → 96 µm): B_p ⊕ B_w reproduced in
real diffraction.

**T4 — the multiplicative regime is SUB-MILLIMETRE.** Developed speckle:
pointwise-model blend α=0.10/0.25/0.50 at z2 = 0.065/0.163/0.327 mm; weak
screen holds to ~1.4–2 mm. Beyond = the scrambling regime (your rank-one Q
world) — where the sentinel survives by construction.

**T5 — end-to-end sentinel (CAVEAT: M=32 codes vs sealed M=128; an M-scaling
rerun is queued — treat absolute gaps as provisional).** Near-contact
(z2=1 mm): the MEAN channel beats covariance ~3× (mean_d′/cov_d′ ≈ 2.9–3.1) —
the T1 leak is operationally significant exactly where contact degeneracy
starves covariance; T_det(d′=5, 2%) ~184k banks. Mid (z2=5 mm): covariance
beats mean ~1.6×; 5% change AUC 0.96 at the sealed T=453 budget; 2% ~12.4k
banks.

An internal 5-lens divergence workflow (higher/finer/stranger/next-paper/bench)
is running in parallel; your R44 will be merged against it. You have repeatedly
out-derived it — derive first, and give your own innovations (standing
operator rule).

## R44 asks

1. **The exact-leak-support question (derive!).** The intensity spectrum of a
   band-limited field is supported in twice the field band. Is the T1 leak
   spectrally CONFINED to [k_p, 2k_p] (plus a micromirror-sinc pedestal with
   its own known comb structure)? If yes: a witness placed beyond 2k_p
   restores an EXACT physical wall — our current witness (k∈[6,9], k_p=5) sits
   exactly inside the suspect band, so this single theorem would rewrite the
   bench design. State assumptions (pixelation, fill factor, z1 Fresnel
   phase), derive the support statement or its failure, give the decisive twin
   test (leak vs witness-annulus radius sweep, ≤60 GPU-min).
2. **The two-channel phase diagram.** Mean-leak channel vs covariance channel
   dominance as a function of (z1, z2): derive the crossover law (contact
   degeneracy starves covariance as z2→0 while leak feeds mean; conversion
   completes and covariance wins as z2 grows). Is there a "physical jet order"
   — the contact order of the physically-realized channel — interpolating the
   statistical theorem's exact classes? Any exact sum rule across the two
   moment layers (leak gained = wall lost, quantitatively)?
3. **Transfer-paragraph frozen wording (the write-back).** Rule the sentence
   set for the Letter's bench-transfer paragraph given T1/T2/T4/T5: must state
   (i) mean-channel blindness is exact at code-algebra level, physically
   near-exact with leak ≤ X(z1) suppressed by near-conjugate relay +
   apodization (and, if ask 1 lands, by witness placement beyond 2k_p);
   (ii) the thin-screen declared law maps to z2 ≲ 0.3 mm (developed) / ~2 mm
   (weak), with the scrambling regime beyond — sentinel surviving throughout;
   (iii) no implication that the sealed grid's latencies or coarse-grain cells
   map 1:1 to this bench geometry. Freeze the exact 3–4 sentences.
4. **Next-paper ruling.** The twin machinery makes the memory-effect
   interpolation family (scene-Fisher rank + sentinel d′ vs z2, from thin
   endpoint through the scrambling endpoint) demonstrable end-to-end in real
   diffraction physics. Is THIS the next paper (the connecting family between
   the Letter's two anchored endpoints)? Verdict + campaign sketch (sweeps,
   cells, GPU budget, the one hero figure, what theorem the data would
   demand). Rank against: physical specificity recalibration (closing D3/D5
   with physically-modeled medium events), depth localization by channel
   crossover, leak-as-asset two-channel sentinel.
5. **Your own ranked innovations (≥5, standing rule).** Anything the five
   internal lenses and asks 1–4 miss — surprise × rigor × reach scored, each
   with mechanism, exact-statement sketch, decisive ≤60-GPU-min kill test on
   the existing engines, nearest prior art.

Constraints: manuscript frozen except the transfer paragraph; sim-only until
publication; dead lines stay dead (beyond-band imaging, third-order cumulants,
code design, pathwise estimation, Rank-1 duality, DLGI, RLMI). Deliver as a
GitHub issue titled R44. The twin artifacts commit next push; the fix-pass
Letter is at e6213b1.
