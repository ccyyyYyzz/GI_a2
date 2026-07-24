# BOLD PUSH — three attacks on the detection sentinel's limits

**Date:** 2026-07-24 · **Engine:** `bold_push.py` (corrected-shot Fisher, same as validated ROC) ·
**Data:** `BOLD_PUSH.json`, `BOLD_PUSH_frontier.png`, `BOLD_PUSH_latency.png`. Local RTX 4060 (~30 s).
Every headline **MC-validated (≥100–150 trials)**. Bank model: `T_eff` = independent banks,
12.8 ms/bank; times in banks and seconds.

---

## WHAT WAS UNLOCKED

### 🔓 UNLOCK 1 — the 0.5 % floor is **broken** (Attack 2). It was integration time, not physics.
The prior map said "0.5 % below floor at `T≤4096`." Push `T` to 16384 banks (**3.5 min**) and the
best cell's frontier crosses it cleanly:

| T (banks) | 256 | 512 | 1024 | 2048 | 4096 | 8192 | **16384** |
|---|---|---|---|---|---|---|---|
| ε_min best, energy-spread | 3.29 % | 2.33 % | 1.65 % | 1.16 % | 0.82 % | 0.58 % | **0.41 %** |
| ε_min best, matched-target | 2.97 % | 2.10 % | 1.49 % | 1.05 % | 0.74 % | 0.53 % | **0.37 %** |
| ε_min median, energy-spread | 4.62 % | 3.26 % | 2.31 % | 1.63 % | 1.15 % | 0.82 % | 0.58 % |

**0.5 % is reached at** best cell **≈9 000–11 000 banks (116–142 s)**, median cell **≈13 000–22 000
banks (169–279 s)** — depending on matched-target vs energy-spread. **MC-validated:** matched-target
detector, ε = 0.75 %, T = 4000 → empirical `d' = 4.61` (analytic 5.0), AUC = 1.000. The `d'∝√T`
frontier is exactly what the corrected ROC already confirmed, now pushed to the floor. **There is no
0.5 % wall — there is a ~2-minute integration cost.**

### 🔓 UNLOCK 2 — a streaming CUSUM **instrument**, not a batch test (Attack 3).
Turning the batch LR into a sequential change-point detector (CUSUM on per-bank covariance
innovations) gives a **continuous sentinel with an explicit latency–false-alarm curve** — the
top-venue figure. Detection latency at a fixed false-alarm rate (best cell, matched target):

| anomaly ε | latency @ 1 false-alarm/hour | latency @ 1 false-alarm/day |
|:--:|:--:|:--:|
| **5 %** | 83 banks = **1.1 s** | 111 banks = 1.4 s |
| **2 %** | 418 banks = **5.4 s** | 595 banks = 7.6 s |
| **1 %** | 1366 banks = **17.5 s** | 2071 banks = 26.5 s |
| **0.5 %** | 4248 banks = **54.4 s** | 7052 banks = 90.3 s |

**MC-validated (Siegmund ARL formula):** at ARL0 target 500 the empirical false-alarm ARL0 = **590**
and empirical detection latency = **102 banks** vs Siegmund **111** (ratio **0.92**). The formula is
trustworthy; the latency table above is real. **A 2 %-energy beyond-band intrusion is flagged in
5 s at one false alarm per hour; even 0.5 % in under a minute.** This is the "instrument" figure.

---

## WHAT DID NOT MOVE (honest)

### ▪ Attack 1 — detection-optimal code design: **~3 %, not a lever.**
Optimizing the band-limited code bank's spectral weight profile (edge-concentrated, greedy per-shell,
vs flat) to maximize aggregate beyond-band detection Fisher yields only:
- best cell: best parametric profile `k¹` → `d'×1.017`; **greedy per-shell** (weights
  `[1.23, 0.82, 0.55, 0.55, 1.85]` — edge + low concentrated) → `d'×1.031`;
- median cell: flat is already optimal (`×1.00`).

**Why so small:** `M = 128` band-limited codes already **near-completely span the ~121-dim in-band
space**, so reweighting the code spectrum barely changes the code×medium sideband coverage of the
beyond-band shells. The coordinator's edge-concentration intuition is *directionally right* (greedy
does load the band edge, shell 5 weight 1.85) but the ceiling is a ~3 % `d'` gain (≈6 % `T_det`).
Code design is **not** the boundary lever; **integration time (T) is.** Honest negative.

---

## BOLD LEDGER — the sentinel's real operating envelope now

Combining the corrected map with these pushes (best cell, matched-target where a change class is
declared; energy-spread otherwise):

| ε | detect within 4096 banks? | frontier T for d'=5 | streaming latency @1/hr |
|:--:|:--:|:--:|:--:|
| 5 % | yes (worst-mode too) | ~120 banks (1.5 s) | 1.1 s |
| 2 % | yes (95 % of cells) | ~470–740 banks (6–9 s) | 5.4 s |
| 1 % | best cells only | ~1900 banks (24 s) | 17.5 s |
| **0.5 %** | **yes, at T≈9–11 k banks (2 min)** | 9044–11106 banks (116–142 s) | 54 s |

## Boundaries (stated plainly, at the end)
- **0.5 % is now reachable but costs ~2 minutes** of integration (11 k banks) — real, not
  instantaneous; below ~0.3 % would need >16 k banks / matched priors and was not pushed further.
- **Code design is inert (~3 %)** — the near-complete band-limited bank leaves little to optimize;
  do not budget a code-design win.
- **Matched-target beats energy-spread by only ~10–20 %** (λ_max/λ_mean ≈ 1.2–1.3 at the best cell);
  declaring a change-signature class helps modestly, not dramatically.
- These are **detection (aggregate)** results — **imaging (mode coverage) remains empty**; the
  demand-tier lives, the supply-tier does not.
- The **bank-acquisition model** (independent banks) makes multi-lag GLS moot (each bank is already
  one independent sample); a continuous-OU deployment would recover the same frontier via the lags.
- Mean-channel specificity (beyond-band anomaly invisible to the mean, medium-change rejected)
  carries over unchanged from `DETECTION_ROC`.

## Files
`bold_push.py`, `BOLD_PUSH.json`, `BOLD_PUSH_frontier.png`, `BOLD_PUSH_latency.png`. No git commit.
