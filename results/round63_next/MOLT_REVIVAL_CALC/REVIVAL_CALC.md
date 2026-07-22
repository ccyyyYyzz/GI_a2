# MOLT revival desk calculation

**Working dir:** `D:\GI_another` · **Outputs:** `results/round63_next/MOLT_REVIVAL_CALC/`
(`REVIVAL_CALC.md`, `revival_calc.json`, `code/revival_calc.py`). No commits.
Physics imported verbatim from the committed probe
(`SATURATION_JAILBREAK_PROBE/code/saturation_core.py`, `gi_operator.py`); nothing rewritten.

---

## Conclusion (first): **STAY DEAD**

> **MOLT stays dead.** Under the frozen rule — *revive only if the median margin
> `= aggregated d'/3 >= 10x` on at least one honest scene class at an honest
> dose* — no honest class clears the bar. Within the R32 dose envelope
> (`<= x40`/mask, bar 2) the best median margin is **4.13x** (high-contrast
> sparse) and **3.81x** (natural control): squarely in the **2–5x = STAY DEAD**
> band we have watched evaporate three times today. Both revival levers fail on
> their own terms: (H1) fiber ambiguities are **not** order-one on the favorable
> classes (`f_med ≈ 0.13–0.17`, same as natural, not `~1`), and (H2) refresh +
> larger-C lift the margin by only **~4%** at the authorized dose — the quenched
> floor is *subdominant to shot noise at `x40`*, so there is no `10–30x` precision
> to harvest. The only way to reach `~13x` is budget **`x400`** — 10x past the
> R32 dose bar — and even there the "favorable" high-contrast class does **not
> beat the natural control** (12.86x vs 11.86x), so there is no scene-class
> revival, only "spend 10x the dose for `sqrt(10)` more `d'`." The two most
> favorable-looking classes (binary sparse, binary structured) are killed
> outright by the `x^2 = x` identity, and the DL-tangent ("learned manifold")
> route dies with the rest: for every manifold dimension the linear bank can
> afford (`d <= K_D+K_S = 102`) the linear coordinates are already identifiable,
> so MOLT certifies nothing new, and where it could help conditioning (`d=100`)
> the discrimination margin at achievable precision is `<= 2` (`<= 6.4` even at
> `x400`).

**Decision: do NOT open a GPT round. MOLT remains retired.**

---

## 1. The margin table

`margin = (aggregated Fisher-weighted d', R31 §5.5) / 3`, computed on the most
**favorable** box-feasible fiber pair `(x, x+h)` with `h ∈ ker([M_D; M_S])`
(identical buckets on the full linear bank the equal-time comparator uses, so
MOLT's `p2` is the *only* channel that can separate them). Median over the class
instances. Banks per spec: `M_D` = 51 dense 50% masks, `M_S` = 51 32-support
masks (combined linear rank = 102).

| scene class | `f_med` (best-Q) | (i) dead `x40` fixed-W | (ii) +refresh `x40` | (iii) +C=14400 refresh `x40` | (iii) +C=57600 fixed-W `x40` | (iv) +`x400` refresh | (v) best: C=57600 `x400` refresh |
|---|---:|---:|---:|---:|---:|---:|---:|
| **(a) binary sparse** `s∈{40,80,150}` | 0.000 | **0.00** | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| **(b) binary structured** | 0.000 | **0.00** | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| **(c) high-contrast sparse** (9:1) | 0.165 | 3.93 | 4.06 | **4.13** | 3.57 | 12.86 | 12.74 |
| **(d) control (natural bridge)** | 0.130 | 3.55 | 3.75 | **3.81** | 3.30 | 11.86 | 11.76 |
| | | ⟵ **honest dose `<= x40`** ⟶ | | | | ⟵ **`x400`: 10x past R32 dose bar** ⟶ | |

Reading the table:

- **Honest-dose columns (i–iii, `<= x40`):** best is **4.13x** (class c), **3.81x**
  (control). Both in the 2–5x STAY-DEAD band. **No class reaches 10x.**
- **H1 falsified:** the favorable class (c) carries `f_med = 0.165` vs the control's
  `0.130` — *not* order-one, and *barely* above natural. The hoped-for
  "sparse/binary/high-contrast ⇒ order-one `Δp2`" does not happen for the one
  class that is even eligible (binary is killed separately, below).
- **H2 falsified:** dead → refresh at `x40` moves the margin **3.93 → 4.06** (c) and
  **3.55 → 3.75** (d) — a **~3–6%** lift, not `10–30x`. At the authorized dose the
  achievable `p2` SE is shot-noise-limited; the fixed-W *quenched floor is
  subdominant*, so refreshing the speckle recovers almost nothing. Larger `C` at
  fixed budget is **C-independent when refreshed** (iii-refresh ≈ ii) and only
  lowers the fixed-W floor (iii-fixed-W is *below* ii-refresh, 3.57 < 4.06).
- **The only lever that moves the margin is raw dose.** dead → `x400` is
  3.93 → 12.86 = **3.27x ≈ sqrt(10)**: pure photon budget, exactly the naive
  `d' ∝ sqrt(B)` scaling. Reaching 10x *requires* violating the R32 `x40` dose
  bar by an order of magnitude, at which point MOLT is a *worse* photon-for-rank
  trade, and class (c) still does not beat the natural control.

---

## 2. Binary-trick check — kills classes (a) and (b)

For a binary scene `x ∈ {0,1}^P`, `x^2 = x` pixel-wise, so the MOLT curvature
functional collapses onto the ordinary sparse linear bucket:

- **`M_S(x^2) = M_S(x)` exactly** (measured max error **`0.00e+00`**). MOLT's second
  channel is a *bit-for-bit copy* of the sparse linear bucket you already record.
- Hence for any second binary scene `x'`, the MOLT-channel difference
  `Δp2 = M_S(x'^2 - x^2)` equals the **sparse linear bucket difference**
  `M_S(x' - x)` exactly (max error **`0.00e+00`**) — zero information beyond
  ordinary linear GI with the combined bank `[M_D; M_S]`.
- **Known-support R32.1 rank (`s=120`):** `rank[M_D|_S; M_S|_S] = 101`,
  `rank[M_D|_S; M_S|_S; 2 M_S|_S diag(x_S)] = 106` (quad rows add 5). That +5 is a
  *local-Jacobian artifact for continuous perturbations off the binary manifold*;
  the moment you impose the `x^2=x` prior the quadratic **equation** becomes
  `M_S(x)=z` — a duplicate of a linear row — and adds **zero** rank.

**Verdict:** binary scenes don't need MOLT; they need the `x^2 = x` identity (a
free linear/combinatorial prior) plus ordinary linear GI. The saturation sweep
and its extra dose buy nothing. The two most favorable-*looking* classes are
exactly the ones the binary trick kills → **margin `≡ 0` at every rung.**

---

## 3. Competitor check — class (a), `s = 80`

To close a **known** support of size `s = 80` (generic masks, exact):

| route | patterns | extra dose | note |
|---|---:|---|---|
| pure linear (add masks) | **80** (29 beyond `K_D=51`) | ordinary (`~29x B_LIN`) | `s` independent linear rows |
| MOLT, **continuous** scene | **66** (`K_D` + 15 sparse sweeps) | `~15 × x40 = 600x B_LIN` | quad rows give `2` rows/sparse pattern |
| MOLT, **binary** scene | **80** | 600x + wasted | quad rows redundant → **saves 0 patterns** |
| blind CS (unknown support) | 160 (`2s`) – 295 (`s·log2(N/s)`) | ordinary | for reference only |

Pattern-count ratio MOLT/linear = **0.825** (continuous) / **1.000** (binary).
MOLT saves 14 of 80 patterns on a continuous support, paying **~20x** more total
dose to do it; on a binary support it saves nothing. **MOLT wins only when
pattern slots are the strictly binding resource *and* dose is free** — the exact
narrow corner R32 already authorized, and not a revival.

---

## 4. DL-tangent lever — "DL proposes the manifold, MOLT certifies the coordinates"

Proxy for a learned manifold tangent `B_x` (honest stand-in): the leading `d`
PCA directions of a 300-scene class ensemble (high-contrast: its own generator;
natural: smoothed `1/f` random fields, since the frozen bridge generators live
outside the repo — a **linear/global** proxy, not a true nonlinear learned
manifold). Restrict the R32.1 Jacobian to the tangent:
`J_lin = [M_D B_x; M_S B_x]`, `J_joint = [J_lin; 2 M_S diag(x) B_x]`.

| class | `d` | `rank(J_lin)=d`? | `rank(J_joint)=d`? | `σ_min` lin → joint | margin dead `x40` | margin best `x400` |
|---|---:|:--:|:--:|---|---:|---:|
| (c) high-contrast | 20 | **yes** | yes | 2.01 → 2.15 | 3.73 | 12.15 |
| (c) high-contrast | 50 | **yes** | yes | 0.97 → 1.05 | 2.08 | 6.76 |
| (c) high-contrast | 100 | **yes** | yes | 0.077 → 0.237 | 1.95 | 6.35 |
| (d) natural | 20 | **yes** | yes | 2.02 → 2.29 | 0.83 | 2.77 |
| (d) natural | 50 | **yes** | yes | 0.96 → 1.35 | 1.00 | 3.35 |
| (d) natural | 100 | **yes** | yes | 0.049 → 0.424 | 0.81 | 2.73 |

- **Linear already closes the manifold at every `d ∈ {20,50,100}`** (`d <= K_D+K_S
  = 102`). Where DL affords the coordinates, ordinary linear GI already
  identifies them — **MOLT's quadratic rows are redundant for closure.** The
  MOLT-only regime is `d > 102`, i.e. the same photon-for-rank trade at a *larger*
  row deficit.
- MOLT does improve the **conditioning** of the worst tangent direction
  (`σ_min` up ~3–9x at `d=100`), but that never becomes a `≥10x` discrimination
  margin: at the honest dose the DL-tangent margins are `<= 3.73` (and `<= 2.03`
  at `d=100`, where MOLT actually carries the ill-conditioned direction). Even at
  `x400` the `d=100` margins are `6.35` (c) / `2.73` (d) — **below 10x**.

**Verdict:** the DL-tangent revival route dies with the rest. Where MOLT is
needed (large `d`, ill-conditioned coordinate) the margin at achievable precision
is `< 10x`; where the margin is large (`d=20`) MOLT is redundant because linear
already closes it. Closing even a `d=50` tangent to a `10x` margin needs
precision beyond the refresh + large-C ladder — which does not exist at honest
dose.

---

## 5. Method notes (honesty ledger)

- **Physics:** `saturation_core` verbatim. `Q(a;v)=∏_p 1/(1+a v_p)`,
  `L=-log Q = a p1 − a²p2/2 + …`, MOLT operator `p2 = M(x∘x)`, `n_eff=p1²/p2`.
- **Achievable-precision ladder (`rel-se(p2)` per mask):**
  *refreshed-W* = committed annealed FIM at the 3-level ridge design (`three_level_design`
  + `annealed_fim`, order-2, `p1` profiled), `se ∝ 1/sqrt(B)`, C-cancels at fixed
  budget; *fixed-W* = `sqrt(annealed² + quenched_floor²)` with the analytic
  closed-form of the exact S2.1 quench covariance `Cov[r_i,r_j]=(r_{i+j}−r_i r_j)/C`
  propagated through the fit (matches `quenched_floor.json`; floor `∝ 1/sqrt(C_eff)`,
  does not average away).
- **Margin metric** = aggregated Fisher-weighted `d'` over **all** sparse masks
  (`d'² = Σ_k Δp2_k²·I_{p2,k}`), divided by 3, per the task. `f`-summary = median
  `|Δp2|/p2` over the **best-quartile** masks.
- **Fiber pairs** verified to share linear buckets on `[M_D; M_S]` to `< 1e-8`
  relative before use. Continuous classes clamped to `[0.08,0.92]` (the committed
  `bars_dose`/`t3` convention) so box-feasible steps exist — a MOLT-favorable
  choice, and MOLT still fails.
- **Not softened:** the frozen decision rule governs; `x400` rungs are reported
  but flagged as `>10x` past the R32 dose bar and excluded from the revival
  decision. Every honest-dose margin is `<= 4.13x`.
