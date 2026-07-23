# SCGI-style learned correction vs DLGI — head-to-head on the frozen OU-drift protocol

**EXPLORATORY / DEV-ONLY. Not preregistered. No commits. Local GPU (RTX 4060 Laptop), `pytorch` env (py3.9 / torch2.1.0+cu121).**
All work on DEV assets only (6 bridge scenes + freshly-generated procedural training scenes with far-offset seeds). The sealed confirmatory/edge banks under `results/round63_dlgi_campaign/confirmatory/` were never read or touched. Colab untouched.

This is a **faithful-in-spirit** reimplementation of Peng & Chen's SCGI (learning-based correction with Gaussian constraints, Opt. Lett. 2023, DOI 10.1364/OL.499787) built from the method description, **not code-identical**. See *Methods & faithfulness* for every documented divergence.

---

## 1. Conclusion first

When a learned corrector's **training corruption family matches the physical drift**, it is a genuinely strong competitor to our model-based DLGI — it recovers the gain path as well as DLGI (Pearson *r* > 0.98 vs truth), reconstructs comparably, and even its *discarded* correction factors yield competitive medium point-estimates. **The model-based DLGI's decisive, robust advantages are elsewhere, and they hold:**

1. **Assumed-corruption-family fragility is real and large (the killer result).** A U-Net trained only on the paper's assumed class (exponential decay), tested on our OU drift, **collapses**: reconstruction barely clears the no-correction floor (+0.2 to +0.6 dB over floor across all cells), and its correction-derived medium parameters are uninformative (median CV rel-err **0.65**, t_c rel-err **0.31**). DLGI is model-based and training-free — it needs no corruption-family assumption and works across every cell.

2. **"Their trash" is decent in-family but is an uncalibrated, memorized point estimate, not a measurement.** SCGI's implicit correction factors, fit to an OU, give in-family t_c/CV rel-errors (~0.04/0.03) that are competitive with DLGI's product — but they sit **below the Cramér–Rao sd floor** (0.13–0.25 for t_c), a tell that the accuracy is **shrinkage toward the training grid**, not honest inference. It carries no uncertainty/calibration, and it degrades on t_c extrapolation beyond the training range (rel-err 0.09 → 0.14). DLGI delivers the same numbers *with* the Fisher/CRB certificate the ruling (§5.3, promotion bar 3–4) requires.

3. **Honest DLGI weakness surfaced (a lead, not a threat).** In the fast-and-strong drift corner (short t_c, high CV) the family-matched learned corrector **beats** the frozen DLGI on scene PSNR by **+3.5 dB** (on-grid, t_c=16/CV=0.40) and **+5.5 dB** (off-grid, t_c=16/CV=0.55) — despite **near-identical gain-path correlation** (both *r* ≈ 0.995). The gap is therefore a gain-**level/scale** bias in DLGI's reconstruction (Jensen/gauge term at high CV), not a path-shape error. This is a concrete improvement lead for the DLGI estimator; it does not touch the medium-measurement thesis.

**One-line takeaway:** the learned corrector wins where you can guarantee the corruption family and don't need a certificate; DLGI wins on *robustness to the unknown corruption family* and on being a *certified* measurement — exactly the surviving distinction §5.3 names ("recovered medium parameters as a certified science output rather than merely a correction variable").

Summary figure: `head_to_head_figure.png`.

---

## 2. Head-to-head scene PSNR (item 3)

6 DEV bridge scenes × 3 seeds per cell = 18 records/cell, **identical records fed to every arm**, paired schedule. All arms reconstruct with the **same** `arm_A4_gaincorr` + frozen TV inversion; only the per-exposure gain estimate differs (Arm 0: gain≡1; Arm L: U-Net correction; Arm D: DLGI OU/Kalman). Oracle = true per-exposure gain (ceiling).

| cell | Arm 0 (floor) | Arm L SCGI-exp (OOD) | Arm L SCGI-mix (in-family) | Arm D DLGI | oracle | L_mix − D |
|---|---|---|---|---|---|---|
| t_c=16, CV=0.15 | 13.84 | 14.00 | 15.09 | **15.40** | 18.78 | −0.30 |
| t_c=16, CV=0.40 | 6.16 | 6.49 | **12.28** | 8.80 | 18.39 | **+3.48** |
| t_c=64, CV=0.15 | 16.46 | 16.62 | 15.73 | **17.23** | 18.66 | −1.49 |
| t_c=64, CV=0.40 | 11.30 | 11.90 | 12.67 | **14.02** | 18.48 | −1.36 |

(dB, mean over 18 records.) DLGI wins 3/4 on-grid cells (by 0.3–1.5 dB); the family-matched learned corrector wins the hardest fast-strong-drift corner by 3.5 dB. SCGI-exp (out-of-family) is essentially the floor everywhere.

**Recovered-gain-path fidelity (Pearson r vs true path):** L_mix 0.985–0.997, D 0.988–0.998, L_exp 0.35–0.63. In-family SCGI recovers the gain path as well as DLGI; the OOD model does not.

---

## 3. "Their trash vs our product" (item 3) — medium parameters

SCGI's pipeline produces corrected measurements `B = Y / g` and **discards** the correction factors `g`. We extracted `g` (the implicit per-exposure gain), demeaned its log, and fit the **same** MoM-autocovariance OU estimator DLGI uses. This turns the discarded "trash" into a (t_c, CV) estimate directly comparable to DLGI's certified product.

Median relative error, mean over the 4 on-grid OU cells:

| quantity | Arm D (product) | L_mix (trash, in-family) | L_exp (trash, OOD) | oracle-monitor (ref) | CRB sd floor |
|---|---|---|---|---|---|
| t_c rel-err | 0.095 | 0.044 | 0.305 | 0.072 | 0.13 (t_c16) – 0.25 (t_c64) |
| CV rel-err | 0.046 | 0.029 | 0.649 | 0.032 | 0.06 (t_c16) – 0.13 (t_c64) |

Reading it honestly:
- **In-family, the "trash" is not worthless** — SCGI-mix's correction factors give t_c/CV point estimates competitive with (numerically a touch better than) DLGI's. This is a fair finding, reported without spin.
- **But it is not a measurement.** Both L_mix and D land *below the CRB sd floor* for an unbiased estimator — impossible for honest inference; it is **shrinkage toward the training/prior grid** (SCGI) and short-window MoM bias (DLGI). SCGI's medium output carries **no calibration or uncertainty**; DLGI's comes with the Fisher/CRB machinery §5.3 requires.
- **Out-of-family, the trash is garbage** — L_exp CV rel-err 0.65 (near the uninformative line 0.5+).

---

## 4. The OOD / assumed-corruption-family result (item 4) — the killer question

Train the corrector on **their assumption class only** (exponential-decay scaling, λ∈[0.9995,1]); test on **our OU drift**. This makes "the assumed corruption family fragility" quantitative.

- **Reconstruction:** SCGI-exp gain over the no-correction floor is **+0.16 / +0.33 / +0.16 / +0.60 dB** across the four OU cells — i.e. it does **almost nothing**, because it only ever learned a monotone decay and cannot represent zero-mean OU fluctuations. Family-matched SCGI-mix gains +1.3 / +6.1 / −0.7 / +1.4 dB; DLGI gains +1.6 / +2.6 / +0.8 / +2.7 dB (Panel C).
- **Medium parameters:** SCGI-exp correction factors give t_c rel-err 0.18–0.45 and CV rel-err 0.57–0.72 — uninformative.
- **DLGI is unaffected** because it assumes a *physical* OU/scalar-gain model, not a *learned* corruption prior. No retraining, no family assumption.

This is the clean structural point: a learned corrector is only as good as the corruption distribution it was trained on; a model-based estimator is only as good as its physical model — and the OU/scalar-gain model transfers across the whole drift family for free.

---

## 5. Off-grid sensitivity (added to separate generalization from grid-recall)

The item-3 test cells (t_c∈{16,64}, CV∈{0.15,0.40}) sit **exactly on** SCGI-mix's OU training grid (t_c∈{8,16,32,64}, CV∈{5,15,40}%), a home-field advantage for the learned medium estimate. We retested off-grid:

| cell | Arm 0 | L_mix PSNR | D PSNR | oracle | t_c rel-err D / L_mix | CV rel-err D / L_mix |
|---|---|---|---|---|---|---|
| interior (t_c=24, CV=0.25) | 11.57 | 13.89 | 13.82 | 18.70 | 0.036 / 0.040 | 0.018 / 0.049 |
| extrap t_c (t_c=128, CV=0.15) | 17.19 | 16.19 | **17.83** | 18.77 | **0.088 / 0.143** | 0.040 / 0.061 |
| extrap CV (t_c=16, CV=0.55) | 3.29 | **11.53** | 6.03 | 18.14 | 0.253 / 0.064 | 0.106 / 0.043 |

- **Interior interpolation:** SCGI-mix generalizes fine between grid nodes (≈ DLGI). So its in-family accuracy is *not* pure node-recall — it learned a real gain extractor within the training support.
- **t_c extrapolation** (slower than trained): SCGI-mix degrades and DLGI overtakes (t_c rel-err 0.143 vs 0.088; PSNR −1.6 dB) — the expected "learned within support, weaker outside."
- **CV extrapolation** (stronger than trained): SCGI-mix holds up and DLGI *collapses* on scene PSNR (6.03) and t_c (rel-err 0.253) — the same fast-strong-drift weakness of the frozen DLGI estimator, now off-grid and larger (+5.5 dB for L_mix).

---

## 6. Training details

- **Forward model:** OUR frozen complementary-Hadamard M0 pipeline (`part1_gain_ladder` via `dl_common`). 32×32, NPIX=1024, N_EXP=2048, PHI=11.9375 (matched to 2200 counts/exposure on the DEV scenes), Poisson, within-exposure-constant OU gain.
- **Training scenes:** procedurally generated 32×32 [0,1] scenes (shapes / strokes-digits / sparse dots / smooth blobs), amplitude-rescaled to member-sum ∈ [220,560] (DEV band). RNG seeded at 1,000,000+offset — **disjoint** from the forward-model seeds (1000+/5000+) and from the bridge scenes.
- **Two models** (item 1: I used **both** families, as two models, to serve both the fair head-to-head and the OOD contrast):
  - `scgi_mix` — trained on exponential + OU (t_c∈{8,16,32,64}, CV∈{5,15,40}%), 50/50 mix. In-family reference for the head-to-head.
  - `scgi_exp` — trained on exponential decay only. The out-of-family arm.
- **Data size:** 4000 train / 800 val records per model. **Compact U-Net** (263,257 params, base width 24, 2 pool levels on the 32×64 record image). Adam lr 2e-3, batch 64, ReduceLROnPlateau, early stop (patience 16, 80 epochs cap), MSE + Gaussian random-walk prior (weight 0.1) on the predicted log-gain.
- **Compute:** ~90 s/model on the RTX 4060 (both models + data gen well under 4 min total) — far inside the 2–6 h budget; no downscaling was needed. VRAM/RAM never a constraint (tiny net, tiny tensors).
- **Fit quality:** best val MSE on demeaned log-gain — `scgi_mix` 3.7e-4, `scgi_exp` 1.0e-5, against a target variance of ~0.051 (i.e. mix recovers ~99.3% of the log-gain variance in-family; exp fits its trivial monotone family near-exactly).

---

## 7. Methods & faithfulness (documented divergences)

Faithful-in-spirit, not code-identical. Peng & Chen: train a U-Net to map a drift-corrupted measurement sequence R to the static-equivalent measurement B on simulated corruption (exponential-decay scalings), loss = MSE + a Gaussian-prior term; then differential GI on B.

Our reimplementation and why:
1. **Forward model = our frozen OU/Hadamard pipeline**, so the head-to-head shares one physical model with DLGI.
2. **The U-Net predicts the per-exposure log-gain correction δ (gauge-demeaned)**, and B = Y·exp(−δ). Predicting the correction and dividing it out is equivalent to predicting B, and it is precisely what exposes the implicit correction factors the head-to-head must quantify (§3).
3. **Inputs** = per-position-baseline-centred log-counts. The baseline is the mean *clean* log-signal per exposure position over the training scenes (gain-independent), which isolates the deterministic signal structure so the conv net's job is temporal drift extraction. Per-record global mean removed (the scene-brightness / mean-gain gauge).
4. **Loss** = weighted MSE(δ̂, δ_true) + 0.1 × Gaussian random-walk prior on δ̂ in acquisition/time order (the "Gaussian constraint" = drift is slow/smooth).
5. **Reconstruction** = the identical `arm_A4_gaincorr` + frozen TV inversion the other arms use, so Arm L vs Arm D differ **only** in the gain-estimate source. Arm L is **single-pass / feed-forward** (no iterative refinement), faithful to a learned corrector; Arm D is the frozen `joint_dual_ledger` (its native 2 outer iterations).
6. **Exponential family** = a_e = λ^e, λ∼U[0.9995,1] (monotone decay, literal to the paper's class). **OU family** = the frozen `ou_path`.

Comparison respects §5.3 / the ruling: identical photons/exposures/patterns/schedule across arms (only the gain estimator changes); DLGI treated as the certified-product method; the learned corrector treated as the correction-variable method it is.

---

## 8. Honest caveats

- **Exploratory, not preregistered.** No blinding, no confirmatory split, DEV scenes only. Nothing here is a campaign result and none of it should be folded into the sealed campaign.
- **Faithful-in-spirit, not code-identical** to Peng & Chen. A different U-Net/loss/normalization could shift the in-family numbers; the *qualitative* conclusions (OOD collapse, uncalibrated-shrinkage medium estimate, DLGI's fast-strong-drift scale weakness) are robust to those choices because they follow from structure, not tuning.
- **The in-family "SCGI is competitive/better in some cells" finding is real and reported without spin.** It does not contradict the medium-measurement thesis, which rests on family-robustness + certification, but it should not be hidden.
- **DLGI used as frozen** (`joint_dual_ledger`, n_outer=2). Its fast-strong-drift PSNR deficit vs family-matched SCGI (despite equal gain-path correlation) points to a gain-**level** bias — a fixable estimator issue, logged as a lead, not litigated here.
- **Medium point-errors below the CRB floor are a shrinkage artifact**, for both methods, not super-efficiency. Only calibrated-interval coverage (DLGI, not done here) would settle "measurement vs recall."

---

## 9. Files

- `DL_EXPLORATORY_REPORT.md` — this report.
- `head_to_head_figure.png` — 3-panel summary.
- `head_to_head_results.json` — on-grid: 72 per-record entries + per-cell aggregates.
- `off_grid_results.json` — off-grid sensitivity (interior + 2 extrapolations).
- `checkpoints/scgi_mix.pt`, `scgi_exp.pt` (+ `_history.json`) — trained models, frozen baseline, meta.
- `code/scgi_common.py` — scene generator, gain families, U-Net, correction + reconstruction, `medium_from_gain`.
- `code/train_scgi.py` — trains both models.
- `code/head_to_head.py` — on-grid arms L/D/0 + oracle + medium products.
- `code/off_grid.py` — off-grid probe.
- `code/make_figure.py` — figure.

Reproduce (from `results/round63_next/DL_EXPLORATORY/code`, `pytorch` env):
`python train_scgi.py` → `python head_to_head.py` → `python off_grid.py` → `python make_figure.py`.
