# Attribution Simplex + Calibration-side D1 (dry run) · **PASS**

**R41 §4.6 D1 + §9 Rank 2** · engine `run_d1_cal.py` (sha `b1a93d24e0089b42`) · data
`SIMPLEX_CALIBRATION.json` · **calibration bank only** (confirmatory sealed). Wall 53 s.

---

## Part A — Calibration-side D1 (analytic vs empirical d')

FIXED-COV empirical `d'` vs the corrected-shot analytic `d'` on **calibration scenes** (a witness with
beyond-band structure + the `camera` natural), best/mid/floor cells, `ε∈{1,2,5}%`, 300 records/point,
`T_eff = min(T_det, 4096)`.

**18 points · median ratio 1.022 · range [0.931, 1.140] · median |ARE| 4.0%.**

| check (R41 §4.6 D1) | requirement | result |
|---|---|---|
| `d'_emp/d'_an ∈ [0.80,1.20]` (all points) | all in band | **PASS** (min 0.931, max 1.140) |
| median absolute relative error | ≤ 10% | **PASS** (4.0%) |
| none outside ±30% | all in [0.70,1.30] | **PASS** |

Representative points:

| scene | cell | ε | T_eff | d'_analytic | d'_empirical | ratio | AUC |
|---|---|---|---|---|---|---|---|
| witness | best | 1% | 1893 | 4.99 | 5.02 | 1.00 | — |
| witness | best | 2% | 473 | 4.93 | 5.23 | 1.06 | — |
| witness | mid | 1% | 4096 | 4.15 | 4.73 | 1.14 | — |
| witness | floor | 2% | 1955 | 4.90 | 4.79 | 0.98 | — |
| camera | best | 2% | 305 | 5.17 | 5.65 | 1.09 | — |
| camera | floor | 5% | 259 | 4.77 | 4.78 | 1.00 | — |

The corrected `|P|x`-shot analytic Fisher predicts the actual profiled-covariance detector across
scenes, cells, and anomaly sizes with no systematic optimism (median ratio 1.02, slightly conservative
at the `T_eff=4096`-capped 1% mid/floor points where the analytic under-reads by ≤14%). **D1 (calibration
side) = PASS.** The 27-cell confirmatory D1 is sealed.

---

## Part B — Attribution simplex Gram (R41 §9 Rank 2)

Four tangent spaces in the whitened joint observation space **[mean ⊕ lag-0 cov ⊕ lag-1 cov]**
(joint dimension **16 640**), on the best-cell geometry with a calibration witness scene:

1. `mean/in-band scene` — lights the mean channel `∂m=P·U_in` + a lag-0 covariance signature;
2. `cov/beyond-band scene` — the headline detector (lag-0 cov; mean = 0, the wall);
3. `cov-amplitude` — medium `σ_f` direction (lag-0 cov ∝ C);
4. `cov-lag` — medium correlation-time `τ` direction (lag-1 cross-bank cov ∝ `dρ·C`).

### Raw pairwise canonical correlations (max canonical cosine)

| pair | max cc | reading |
|---|---|---|
| **beyond ↔ amplitude** | **0.039** | **KEY** — beyond-band ⟂ medium amplitude ✓ (`< 0.10`) |
| **beyond ↔ lag** | **0.026** | **KEY** — beyond-band ⟂ medium timescale ✓ (`< 0.10`) |
| beyond ↔ in-band | 0.187 | moderate cov overlap; profiled out + separated by the mean channel |
| in-band ↔ amplitude | 0.988 | strong cov overlap — **separated by the mean channel** (in-band lights `∂m`, amplitude does not) |
| in-band ↔ lag | 0.654 | separated by mean channel + lag block |
| amplitude ↔ lag | 0.662 | separated by lag block (lag-0 vs lag-1) |

**The specificity claim that matters is met intrinsically:** the beyond-band sentinel is orthogonal to
both medium-amplitude (0.039) and medium-timescale (0.026) changes, well under the 0.10 target — a
beyond-band alarm cannot be produced by, or aliased as, a change in the fog's contrast or correlation
time. The large in-band↔amplitude / amplitude↔lag correlations live entirely in the covariance channel
and are resolved operationally by (i) the exact mean channel (in-band lights `∂m`, medium does not) and
(ii) the lag-0 / lag-1 block split (amplitude vs timescale).

### Efficient scores

Projecting each tangent off the union of the other three yields the efficient (nuisance-profiled)
scores. By construction they are **mutually orthogonal** (max off-diagonal efficient canonical
correlation = **0.0**) — the formal "asymptotically independent Gaussian score coordinates." Efficiency
retained after profiling: in-band 0.997, beyond 0.990, lag 0.750, amplitude 0.132 (the amplitude
covariance direction is 99% spanned by the in-band tangent, but the two are separated by the mean
channel — see the classifier).

### Five-class balanced accuracy (bar D3 machinery, LDA)

Nearest-centroid in the pooled-within-class-covariance (LDA) metric over the four score statistics
`{t_mean, t_beyond, t_amp, t_lag}`, on `{H0, in-band, beyond-band, medium-amplitude, medium-timescale}`,
300 records/class:

| T_eff (banks) | balanced accuracy | H0 | in-band | beyond | amplitude | lag |
|---|---|---|---|---|---|---|
| 1200 (15.4 s) | 0.863 | 0.687 | 0.690 | 0.937 | 1.00 | 1.00 |
| **2048 (26.2 s)** | **0.925** | 0.863 | 0.797 | 0.967 | 1.00 | 1.00 |

**Balanced accuracy = 0.925 ≥ 0.90 at T=2048 banks** (within the 4096 cap). Medium-amplitude and
medium-timescale are perfectly separated (recall 1.00); the beyond-band target recall is 0.967. The
residual H0↔in-band confusion (both "no beyond-band anomaly") does not touch the sentinel's
specificity and continues to shrink with `T_eff`.

---

## Bottom line

The corrected-shot analytic detector is calibrated (D1 median ratio 1.02, |ARE| 4%). The attribution
simplex delivers the Rank-2 target: the beyond-band sentinel is intrinsically orthogonal to the two
medium confounders (0.039 / 0.026 < 0.10), the efficient scores are independent coordinates, and the
five-class balanced accuracy reaches 0.925 within the bank cap. **Calibration D1 + simplex = PASS.**
The confirmatory D1/D3 endpoints remain sealed.
