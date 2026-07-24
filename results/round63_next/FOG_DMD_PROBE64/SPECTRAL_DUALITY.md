# Detection-Estimation Spectral Duality (R41 Rank-1, the PRL lever)
**Verdict: COLLAPSE_LOOSE** — estimation collapse $R^2$=0.09, detection collapse $R^2$=1.000.
Corrected physical-|P|x-shot profiled Fisher engine, all 81 cells, witness scene, $T_{\mathrm{eff}}$=4096, $n$=10$^4$ e$^-$/bucket. No reconstruction.
## The claim
Estimation is governed by the **inverse/lower tail** of the efficient Fisher spectrum $J_B$; isotropic (energy-spread) detection by its **trace**:
```
E[d'^2]/(T_eff eps^2 ||x||^2) = tr(J_B)/d        (detection reads the trace)
R_est >= tr(J_B^+)/T_eff                          (estimation reads the tail)
P := [tr(J_B)/d][tr(J_B^+)/d] >= 1,  = iff isotropic
```
## 1. Engine self-check
- recomputed $f_{rec}$ vs frozen map: max abs err 0.0005, median 0.0003.
- recomputed tr$(J)/d$ vs detection-map `lam_mean_bar`: max rel err 0.0e+00. The spectrum is the same object that produced both frozen maps.

## 2. The two collapses
| collapse | predictor | $R^2$ | Spearman |
|---|---|---:|---:|
| **Estimation** $f_{rec}$ | lower-tail mass (unweighted/energy-spread) | 0.085 | -0.360 |
| Estimation NRMSE$_{CRB}$ | tr$(J^+)/d$ (log-log) | 1.000 | +1.000 |
| **Detection** $T_{det}(2\%)$ energy-spread | tr$(J)/d$ (log-log) | 1.000 | -1.000 |
| Detection $T_{det}(2\%)$ worst-mode | $\lambda_{min}$ (log-log, dual) | 1.000 | — |

**Discriminant controls** (the axes do not cross-talk):

| test | $R^2$ | reads |
|---|---:|---|
| $f_{rec}$ vs tr$(J)/d$ | 0.036 | estimation does NOT read the trace |
| $T_{det}$avg vs $\lambda_{min}$ | 0.580 | isotropic detection does NOT read the tail |
| $T_{det}$worst vs tr$(J)/d$ | 0.580 | worst-mode detection does NOT read the trace |

Scene-weighted lower-tail mass is $\equiv 1-f_{rec}$ by construction ($R^2$=1.000); reported as an identity check, not an independent collapse.

## 3. Duality inequality  P = [tr(J)/d][tr(J^+)/d] >= 1
- Holds in **81/81** cells. Range $P\in[1.05,\,2.6]$; all spectra full-rank (no exact null), so the AM-HM gap is exact.
- **Most isotropic** ($P$=1.05): k^-2 kw1x c1.25 sf1.0 — eff-rank/d=0.97, $f_{rec}$=0.23, $T_{det}(2\%)$=467 banks.
- **Spikiest** ($P$=2.6): k^-2 kw4x c1.8 sf0.3 — eff-rank/d=0.55, $f_{rec}$=0.32, $T_{det}(2\%)$=2340 banks.

## 4. Qualitative high-P prediction (Spearman across 81 cells)
| relation | Spearman | prediction |
|---|---:|---|
| $P$ vs $f_{rec}$ | -0.488 | spiky -> worse imaging (expect <0) |
| $P$ vs NRMSE$_{CRB}$ | +0.716 | spiky -> worse NRMSE (expect >0) |
| $P$ vs $T_{det}(2\%)$ | +0.505 | ruling expected <0 (spiky=strong detector) — **CONTRADICTED**, sign is positive |
| eff-rank/d vs $f_{rec}$ | +0.476 | isotropy -> better imaging (expect >0) |

## Honest ledger — does (9.1) earn the central-theorem slot?
**No. It stays a lemma.** The ruling's promotion test is whether the maps *collapse onto one spectral-anisotropy curve*. They do not, for three concrete reasons:

1. **The two exact collapses are definitional, not predictive.** $T_{det}(2\%)$ vs tr$(J)/d$ ($R^2$=1.000) and NRMSE$_{CRB}$ vs tr$(J^+)/d$ ($R^2$=1.000) are both perfect because the detection deflection and the CRB are *built from these very spectral functionals* (eq 3.3 and the profiled CRB), with $\|x\|^2$ fixed across cells. They confirm the algebra of the duality but are not independent evidence. The independent detection check is the sealed Monte-Carlo probe (R41 section 4), not this analytic map.
2. **The one genuinely non-circular estimation collapse is loose.** $f_{rec}$ — a scene-weighted *count* of illuminated modes — does not collapse onto any scene-independent lower-tail mass ($R^2$=0.085, Spearman=-0.36). $f_{rec}$ is not a spectral invariant; it depends on which modes the witness happens to occupy. The spectral estimation quantity that *is* invariant is tr$(J^+)$, and that collapse is the definitional one in point 1.
3. **The anisotropy product $P$ does not deliver the promised task trade-off.** $P$ correctly anti-predicts imaging ($P$ vs $f_{rec}$ Spearman -0.49; $P$ vs NRMSE +0.72), but has the **wrong sign for detection** ($P$ vs $T_{det}(2\%)$ = +0.50; the ruling predicted spiky $\to$ *stronger* detector, i.e. negative). Concretely the most isotropic cell detects in 467 banks while the spikiest needs 2340 — spiky cells are worse at **both** tasks. The reason: across this physical grid, anisotropy and total Fisher supply co-vary — spiky cells also have a smaller trace — so 'spiky = strong detector' (true only at *fixed* trace) never materialises.

**Consequence for the flagship.** The dead-estimation / live-detection dissociation is real, but it is driven by the **demand threshold** (imaging requires every mode above SNR 3; detection pools the trace), not by spectral anisotropy $P$. Eq (9.1) is an exact per-cell AM-HM identity worth one lemma/remark; it is not the paper's organising principle and is not a PRL-grade central object on this evidence. Promote the demand-threshold mechanism, keep (9.1) as support.

## 5. Per-cell table (sorted by duality product P, spiky -> isotropic)
| cell | $d$ | tr$(J)/d$ | $\lambda_{min}$ | $P$ | eff-rank/d | low-tail | $f_{rec}$ | NRMSE | $T_{det}2\%$avg |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| k^-2 kw4x c1.8 sf0.3 | 240 | 0.0223 | 0.0010 | 2.6 | 0.55 | 0.19 | 0.32 | 0.28 | 2340 |
| k^-2 kw2x c1.8 sf0.3 | 240 | 0.0249 | 0.0012 | 2.5 | 0.56 | 0.14 | 0.31 | 0.26 | 2096 |
| k^-2 kw1x c1.8 sf0.3 | 240 | 0.0285 | 0.0015 | 2.4 | 0.57 | 0.11 | 0.34 | 0.23 | 1834 |
| k^-2 kw4x c1.8 sf0.6 | 240 | 0.0348 | 0.0022 | 2.2 | 0.59 | 0.07 | 0.27 | 0.20 | 1499 |
| k^-2 kw2x c1.8 sf0.6 | 240 | 0.0375 | 0.0025 | 2.2 | 0.60 | 0.05 | 0.23 | 0.20 | 1392 |
| k^-2 kw1x c1.8 sf0.6 | 240 | 0.0412 | 0.0030 | 2.1 | 0.61 | 0.05 | 0.24 | 0.18 | 1265 |
| k^-2 kw4x c1.8 sf1.0 | 240 | 0.0420 | 0.0031 | 2.1 | 0.61 | 0.05 | 0.20 | 0.18 | 1241 |
| k^-2 kw2x c1.8 sf1.0 | 240 | 0.0446 | 0.0035 | 2.1 | 0.62 | 0.04 | 0.20 | 0.17 | 1171 |
| k^-2 kw1x c1.8 sf1.0 | 240 | 0.0482 | 0.0039 | 2.0 | 0.62 | 0.02 | 0.16 | 0.17 | 1083 |
| k^-1 kw4x c1.8 sf0.3 | 240 | 0.0130 | 0.0012 | 1.8 | 0.69 | 0.23 | 0.14 | 0.30 | 4006 |
| k^-1 kw2x c1.8 sf0.3 | 240 | 0.0192 | 0.0020 | 1.7 | 0.72 | 0.11 | 0.13 | 0.24 | 2715 |
| k^-1 kw4x c1.8 sf0.6 | 240 | 0.0220 | 0.0024 | 1.7 | 0.72 | 0.08 | 0.15 | 0.22 | 2372 |
| k^-1 kw1x c1.8 sf0.3 | 240 | 0.0264 | 0.0030 | 1.6 | 0.73 | 0.05 | 0.11 | 0.20 | 1977 |
| k^-1 kw4x c1.8 sf1.0 | 240 | 0.0274 | 0.0033 | 1.6 | 0.74 | 0.05 | 0.11 | 0.20 | 1904 |
| k^-1 kw2x c1.8 sf0.6 | 240 | 0.0283 | 0.0035 | 1.6 | 0.74 | 0.05 | 0.10 | 0.19 | 1840 |
| k^-1 kw1x c1.8 sf0.6 | 240 | 0.0352 | 0.0045 | 1.6 | 0.75 | 0.02 | 0.06 | 0.17 | 1484 |
| k^-1 kw2x c1.8 sf1.0 | 240 | 0.0333 | 0.0043 | 1.5 | 0.75 | 0.02 | 0.08 | 0.17 | 1567 |
| k^-1 kw1x c1.8 sf1.0 | 240 | 0.0395 | 0.0052 | 1.5 | 0.75 | 0.02 | 0.04 | 0.16 | 1321 |
| flat kw4x c1.8 sf0.3 | 240 | 0.0057 | 0.0009 | 1.5 | 0.75 | 0.59 | 0.03 | 0.41 | 9125 |
| flat kw4x c1.8 sf0.6 | 240 | 0.0125 | 0.0020 | 1.4 | 0.77 | 0.17 | 0.04 | 0.27 | 4184 |
| flat kw1x c1.8 sf0.3 | 240 | 0.0270 | 0.0041 | 1.4 | 0.79 | 0.02 | 0.03 | 0.19 | 1935 |
| flat kw1x c1.8 sf1.0 | 240 | 0.0381 | 0.0056 | 1.4 | 0.79 | 0.00 | 0.04 | 0.16 | 1371 |
| flat kw1x c1.8 sf0.6 | 240 | 0.0345 | 0.0053 | 1.4 | 0.79 | 0.02 | 0.03 | 0.16 | 1512 |
| flat kw2x c1.8 sf0.3 | 240 | 0.0151 | 0.0025 | 1.4 | 0.78 | 0.10 | 0.06 | 0.25 | 3449 |
| flat kw4x c1.8 sf1.0 | 240 | 0.0173 | 0.0028 | 1.4 | 0.78 | 0.08 | 0.06 | 0.23 | 3009 |
| flat kw2x c1.8 sf0.6 | 240 | 0.0235 | 0.0040 | 1.4 | 0.79 | 0.03 | 0.07 | 0.20 | 2222 |
| flat kw2x c1.8 sf1.0 | 240 | 0.0281 | 0.0048 | 1.4 | 0.80 | 0.02 | 0.06 | 0.18 | 1856 |
| k^-2 kw4x c1.5 sf0.3 | 104 | 0.0397 | 0.0094 | 1.3 | 0.82 | 0.00 | 0.24 | 0.15 | 1314 |
| k^-2 kw2x c1.5 sf0.3 | 104 | 0.0440 | 0.0107 | 1.3 | 0.82 | 0.00 | 0.25 | 0.14 | 1186 |
| k^-2 kw1x c1.5 sf0.3 | 104 | 0.0498 | 0.0126 | 1.3 | 0.83 | 0.00 | 0.26 | 0.13 | 1047 |
| k^-2 kw4x c1.5 sf0.6 | 104 | 0.0601 | 0.0161 | 1.3 | 0.84 | 0.00 | 0.21 | 0.12 | 868 |
| k^-2 kw2x c1.5 sf0.6 | 104 | 0.0643 | 0.0176 | 1.3 | 0.84 | 0.00 | 0.24 | 0.11 | 811 |
| k^-2 kw1x c1.5 sf0.6 | 104 | 0.0703 | 0.0197 | 1.3 | 0.85 | 0.00 | 0.21 | 0.11 | 742 |
| k^-2 kw4x c1.5 sf1.0 | 104 | 0.0715 | 0.0202 | 1.2 | 0.85 | 0.00 | 0.19 | 0.11 | 729 |
| k^-2 kw2x c1.5 sf1.0 | 104 | 0.0755 | 0.0216 | 1.2 | 0.85 | 0.00 | 0.19 | 0.10 | 691 |
| k^-2 kw1x c1.5 sf1.0 | 104 | 0.0810 | 0.0236 | 1.2 | 0.85 | 0.00 | 0.21 | 0.10 | 644 |
| k^-1 kw4x c1.5 sf0.3 | 104 | 0.0204 | 0.0069 | 1.2 | 0.89 | 0.00 | 0.18 | 0.19 | 2561 |
| k^-1 kw2x c1.5 sf0.3 | 104 | 0.0295 | 0.0106 | 1.2 | 0.90 | 0.00 | 0.18 | 0.16 | 1766 |
| k^-1 kw4x c1.5 sf0.6 | 104 | 0.0336 | 0.0123 | 1.1 | 0.90 | 0.00 | 0.22 | 0.15 | 1551 |
| k^-1 kw1x c1.5 sf0.3 | 104 | 0.0401 | 0.0150 | 1.1 | 0.90 | 0.00 | 0.19 | 0.14 | 1302 |
| k^-1 kw4x c1.5 sf1.0 | 104 | 0.0414 | 0.0156 | 1.1 | 0.91 | 0.00 | 0.17 | 0.13 | 1259 |
| k^-1 kw2x c1.5 sf0.6 | 104 | 0.0427 | 0.0162 | 1.1 | 0.91 | 0.00 | 0.19 | 0.13 | 1221 |
| flat kw4x c1.5 sf0.3 | 104 | 0.0080 | 0.0033 | 1.1 | 0.90 | 0.27 | 0.11 | 0.30 | 6558 |
| k^-1 kw1x c1.5 sf0.6 | 104 | 0.0525 | 0.0203 | 1.1 | 0.91 | 0.00 | 0.18 | 0.12 | 993 |
| k^-1 kw1x c1.5 sf1.0 | 104 | 0.0586 | 0.0227 | 1.1 | 0.91 | 0.00 | 0.18 | 0.11 | 891 |
| k^-1 kw2x c1.5 sf1.0 | 104 | 0.0497 | 0.0192 | 1.1 | 0.91 | 0.00 | 0.18 | 0.12 | 1050 |
| flat kw4x c1.5 sf0.6 | 104 | 0.0174 | 0.0074 | 1.1 | 0.91 | 0.00 | 0.17 | 0.20 | 3004 |
| flat kw4x c1.5 sf1.0 | 104 | 0.0242 | 0.0104 | 1.1 | 0.91 | 0.00 | 0.21 | 0.17 | 2159 |
| flat kw2x c1.5 sf0.3 | 104 | 0.0211 | 0.0091 | 1.1 | 0.91 | 0.00 | 0.23 | 0.18 | 2475 |
| flat kw1x c1.5 sf1.0 | 104 | 0.0537 | 0.0226 | 1.1 | 0.92 | 0.00 | 0.20 | 0.12 | 972 |
| flat kw1x c1.5 sf0.6 | 104 | 0.0487 | 0.0207 | 1.1 | 0.92 | 0.00 | 0.20 | 0.12 | 1071 |
| flat kw1x c1.5 sf0.3 | 104 | 0.0381 | 0.0164 | 1.1 | 0.92 | 0.00 | 0.19 | 0.14 | 1369 |
| flat kw2x c1.5 sf0.6 | 104 | 0.0327 | 0.0143 | 1.1 | 0.92 | 0.00 | 0.23 | 0.15 | 1598 |
| flat kw2x c1.5 sf1.0 | 104 | 0.0390 | 0.0171 | 1.1 | 0.92 | 0.00 | 0.21 | 0.14 | 1336 |
| flat kw4x c1.25 sf0.3 | 48 | 0.0091 | 0.0055 | 1.1 | 0.94 | 0.08 | 0.29 | 0.28 | 5715 |
| k^-2 kw4x c1.25 sf0.3 | 48 | 0.0569 | 0.0266 | 1.1 | 0.96 | 0.00 | 0.25 | 0.11 | 917 |
| flat kw4x c1.25 sf0.6 | 48 | 0.0199 | 0.0122 | 1.1 | 0.94 | 0.00 | 0.46 | 0.19 | 2615 |
| k^-1 kw4x c1.25 sf0.3 | 48 | 0.0255 | 0.0139 | 1.1 | 0.95 | 0.00 | 0.31 | 0.16 | 2043 |
| k^-2 kw2x c1.25 sf0.3 | 48 | 0.0627 | 0.0297 | 1.1 | 0.96 | 0.00 | 0.27 | 0.10 | 832 |
| flat kw2x c1.25 sf0.3 | 48 | 0.0242 | 0.0148 | 1.1 | 0.95 | 0.00 | 0.48 | 0.17 | 2157 |
| flat kw4x c1.25 sf1.0 | 48 | 0.0278 | 0.0170 | 1.1 | 0.95 | 0.00 | 0.48 | 0.16 | 1879 |
| flat kw1x c1.25 sf1.0 | 48 | 0.0620 | 0.0369 | 1.1 | 0.95 | 0.00 | 0.29 | 0.10 | 842 |
| flat kw1x c1.25 sf0.6 | 48 | 0.0563 | 0.0337 | 1.1 | 0.95 | 0.00 | 0.27 | 0.11 | 927 |
| flat kw1x c1.25 sf0.3 | 48 | 0.0440 | 0.0265 | 1.1 | 0.95 | 0.00 | 0.25 | 0.12 | 1185 |
| k^-2 kw1x c1.25 sf0.3 | 48 | 0.0706 | 0.0340 | 1.1 | 0.96 | 0.00 | 0.25 | 0.10 | 739 |
| flat kw2x c1.25 sf0.6 | 48 | 0.0374 | 0.0230 | 1.1 | 0.95 | 0.00 | 0.50 | 0.14 | 1394 |
| flat kw2x c1.25 sf1.0 | 48 | 0.0447 | 0.0274 | 1.1 | 0.95 | 0.00 | 0.54 | 0.12 | 1166 |
| k^-1 kw2x c1.25 sf0.3 | 48 | 0.0366 | 0.0204 | 1.1 | 0.96 | 0.00 | 0.27 | 0.14 | 1427 |
| k^-1 kw4x c1.25 sf0.6 | 48 | 0.0415 | 0.0233 | 1.1 | 0.96 | 0.00 | 0.25 | 0.13 | 1258 |
| k^-1 kw1x c1.25 sf0.3 | 48 | 0.0491 | 0.0278 | 1.1 | 0.96 | 0.00 | 0.27 | 0.12 | 1062 |
| k^-1 kw1x c1.25 sf0.6 | 48 | 0.0637 | 0.0365 | 1.1 | 0.96 | 0.00 | 0.31 | 0.10 | 819 |
| k^-1 kw4x c1.25 sf1.0 | 48 | 0.0507 | 0.0288 | 1.1 | 0.96 | 0.00 | 0.29 | 0.12 | 1029 |
| k^-1 kw1x c1.25 sf1.0 | 48 | 0.0707 | 0.0406 | 1.1 | 0.96 | 0.00 | 0.33 | 0.10 | 738 |
| k^-2 kw4x c1.25 sf0.6 | 48 | 0.0842 | 0.0417 | 1.1 | 0.97 | 0.00 | 0.21 | 0.09 | 619 |
| k^-1 kw2x c1.25 sf0.6 | 48 | 0.0522 | 0.0298 | 1.1 | 0.96 | 0.00 | 0.33 | 0.12 | 1000 |
| k^-1 kw2x c1.25 sf1.0 | 48 | 0.0603 | 0.0347 | 1.1 | 0.96 | 0.00 | 0.27 | 0.11 | 865 |
| k^-2 kw2x c1.25 sf0.6 | 48 | 0.0898 | 0.0449 | 1.1 | 0.97 | 0.00 | 0.19 | 0.09 | 581 |
| k^-2 kw1x c1.25 sf0.6 | 48 | 0.0976 | 0.0492 | 1.0 | 0.97 | 0.00 | 0.27 | 0.08 | 534 |
| k^-2 kw4x c1.25 sf1.0 | 48 | 0.0993 | 0.0503 | 1.0 | 0.97 | 0.00 | 0.21 | 0.08 | 525 |
| k^-2 kw2x c1.25 sf1.0 | 48 | 0.1044 | 0.0532 | 1.0 | 0.97 | 0.00 | 0.23 | 0.08 | 500 |
| k^-2 kw1x c1.25 sf1.0 | 48 | 0.1117 | 0.0572 | 1.0 | 0.97 | 0.00 | 0.23 | 0.08 | 467 |

## Figures
- `duality_collapse_frec.pdf/png` — estimation $f_{rec}$ vs lower-tail mass.
- `duality_collapse_tdet.pdf/png` — detection $T_{det}(2\%)$ vs tr$(J)/d$ (log-log).
