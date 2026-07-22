"""
T3 -- NULL-SPACE BREAK (the money demonstration / jailbreak certificate).

Construct x and x' = x + delta with delta != 0 EXACTLY in null(M) of the 51
binary masks -> all 51 IDEAL linear buckets identical to machine precision.
delta is a real texture projected onto the null space and rescaled so both
scenes stay in [0,1] (NO clipping -> null-space property preserved exactly).

Then at matched total photon budget, compute the discriminability d' between x
and x' from the saturated sweeps (51 patterns x L levels).  Because the SAME
fixed detector measures both scenes, the quenched-W systematic is common-mode
and cancels in the difference -> d' is photon-noise-limited (favorable framing).

Per-window mean diff  Dmu = C*(Q_il(x) - Q_il(x'));  var ~ C*(Qbar - Q2bar).
d'^2 = sum_il M_il * Dmu^2 / var = C * sum_il M_il * DQ^2/(Qbar-Q2bar).
Budget = sum_il M_il * a_il * C * p1_i (incident photons).
Flat allocation: d' = sqrt(kappa * budget).  Also report optimal allocation
(concentrate on most informative (i,l)) as a best-case upper bound.

BAR: d' >= 3 at a budget <= the full linear campaign's (51 * B_LIN).
"""
import json
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import saturation_core as sc
import gi_operator as op

OUT = os.path.dirname(HERE)
DATA = "data/r63_bridge_scenes"
rng = np.random.default_rng(651_707)

C = 3600
L = 12
B_LIN = 1.0e4                       # linear-arm incident photons / pattern
LIN_CAMPAIGN = op.K * B_LIN        # full linear campaign budget (51 * 1e4)

# --------------------------------------------------------------------------- #
# operator + base scene + null-space delta
# --------------------------------------------------------------------------- #
M = op.build_operator()                                   # 51 x 1024, 50% dens
base = "bridge_microtex_0"                                 # texture-rich
# clamp to interior [0.08,0.92] so both scenes have >=0.08 headroom on each side
# (still an image-relevant [0,1] scene; guarantees the NO-CLIP null-space add).
LO, HI = 0.08, 0.92
x = np.clip(op.load_scene(os.path.join(DATA, base + ".npz")), LO, HI)

# texture to hide in the null space: a high-frequency checkerboard + a second
# microtexture scene (rich high-freq content -> large null-space projection)
side = op.SIDE
ii, jj = np.meshgrid(np.arange(side), np.arange(side), indexing="ij")
checker = (((ii // 2 + jj // 2) % 2) * 2.0 - 1.0).ravel()  # +-1 blocks
tex2 = op.load_scene(os.path.join(DATA, "bridge_microtex_3.npz"))
texture = checker + 1.5 * (tex2 - tex2.mean())

# project onto null(M): delta_raw = (I - M^+ M) texture
Mpinv = np.linalg.pinv(M)
delta_raw = texture - Mpinv @ (M @ texture)
delta_raw = delta_raw - Mpinv @ (M @ delta_raw)            # re-project (clean up)
delta_raw = delta_raw / np.abs(delta_raw).max()           # normalise max|.|=1

# rescale so x + alpha*delta_raw stays in [0,1] with safety margin (NO clip)
pos = delta_raw > 1e-12
neg = delta_raw < -1e-12
alpha_max = np.inf
if pos.any():
    alpha_max = min(alpha_max, np.min((1.0 - x[pos]) / delta_raw[pos]))
if neg.any():
    alpha_max = min(alpha_max, np.min((0.0 - x[neg]) / delta_raw[neg]))
alpha = 0.9 * alpha_max
delta = alpha * delta_raw
xp = x + delta

assert xp.min() >= -1e-9 and xp.max() <= 1 + 1e-9, "x' out of [0,1]"
lin_x = M @ x
lin_xp = M @ xp
max_bucket_diff = float(np.max(np.abs(lin_x - lin_xp)))
rel_bucket_diff = float(max_bucket_diff / np.max(np.abs(lin_x)))
print("null-space delta: ||delta||_2=%.4f  max|delta|=%.4f  alpha=%.4f" %
      (np.linalg.norm(delta), np.abs(delta).max(), alpha))
print("LINEAR BUCKETS x vs x':  max abs diff = %.3e  (rel %.3e)  over %d masks"
      % (max_bucket_diff, rel_bucket_diff, op.K))
print("  -> all %d ideal linear buckets identical to machine precision: %s"
      % (op.K, "YES" if max_bucket_diff < 1e-8 else "NO"))

# --------------------------------------------------------------------------- #
# per-(pattern, level) discrimination efficiency
# --------------------------------------------------------------------------- #
eff_sum = 0.0            # sum ratio  DQ^2/(Qbar-Q2bar)        (for d'^2)
occ_sum = 0.0           # sum a*p1                              (for budget)
eff_per_photon = []     # per-(i,l) d'^2 per incident photon (for optimal alloc)
p2diff = []
for i in range(op.K):
    m = M[i]
    vx = (m * x); vxp = (m * xp)
    vxa = vx[vx > 0]; vxpa = vxp[vxp > 0]
    p1 = float(vx.sum())
    a_lv = sc.default_a_levels(vxa, L=L)                    # sweep for this pattern
    Qx = sc.Q_exp(a_lv, vx); Qxp = sc.Q_exp(a_lv, vxp)
    Q2x = sc.Q_exp(2 * a_lv, vx); Q2xp = sc.Q_exp(2 * a_lv, vxp)
    DQ = Qx - Qxp
    Qbar = 0.5 * (Qx + Qxp); Q2bar = 0.5 * (Q2x + Q2xp)
    var_ppw = (Qbar - Q2bar)                               # per-window var of hatQ * C
    ratio = C * DQ ** 2 / var_ppw                          # d'^2 per window at (i,l)
    incid = a_lv * C * p1                                  # incident/window at (i,l)
    eff_sum += ratio.sum()
    occ_sum += incid.sum()
    eff_per_photon.extend((ratio / incid).tolist())
    p2diff.append(float(sc.power_sums(vxpa, 2)[2] - sc.power_sums(vxa, 2)[2]))

eff_per_photon = np.array(eff_per_photon)
kappa_flat = eff_sum / occ_sum                # d'^2 = kappa_flat * budget (flat M)
kappa_opt = float(eff_per_photon.max())       # best single-(i,l) efficiency
budget_d3_flat = 9.0 / kappa_flat
budget_d3_opt = 9.0 / kappa_opt
print("\nkappa (d'^2 per incident photon):  flat=%.3e   optimal(best i,l)=%.3e"
      % (kappa_flat, kappa_opt))
print("budget for d'=3:  flat=%.3e   optimal=%.3e   (linear campaign=%.3e)"
      % (budget_d3_flat, budget_d3_opt, LIN_CAMPAIGN))
print("d' at linear-campaign budget (%.1e):  flat=%.3f   optimal=%.3f"
      % (LIN_CAMPAIGN, np.sqrt(kappa_flat * LIN_CAMPAIGN),
         np.sqrt(kappa_opt * LIN_CAMPAIGN)))

# d'-vs-budget curve
budgets = np.logspace(np.log10(LIN_CAMPAIGN / 100), np.log10(budget_d3_opt * 30), 40)
dprime_flat = np.sqrt(kappa_flat * budgets)
dprime_opt = np.sqrt(kappa_opt * budgets)

report = dict(
    C=C, L=L, K=op.K, base_scene=base, B_LIN=B_LIN, linear_campaign=LIN_CAMPAIGN,
    delta_l2=float(np.linalg.norm(delta)), delta_max=float(np.abs(delta).max()),
    alpha=float(alpha), max_bucket_diff=max_bucket_diff,
    rel_bucket_diff=rel_bucket_diff,
    linear_buckets_identical=bool(max_bucket_diff < 1e-8),
    mean_p2_diff_per_pattern=float(np.mean(p2diff)),
    kappa_flat=float(kappa_flat), kappa_opt=float(kappa_opt),
    budget_d3_flat=float(budget_d3_flat), budget_d3_opt=float(budget_d3_opt),
    dprime_at_linear_campaign_flat=float(np.sqrt(kappa_flat * LIN_CAMPAIGN)),
    dprime_at_linear_campaign_opt=float(np.sqrt(kappa_opt * LIN_CAMPAIGN)),
    budget_d3_flat_over_campaign=float(budget_d3_flat / LIN_CAMPAIGN),
    budget_d3_opt_over_campaign=float(budget_d3_opt / LIN_CAMPAIGN),
    dprime_curve=dict(budget=budgets.tolist(), flat=dprime_flat.tolist(),
                      optimal=dprime_opt.tolist()),
    BAR="d' >= 3 at budget <= linear campaign (51*B_LIN)",
    PASS=bool(budget_d3_flat <= LIN_CAMPAIGN))

# --------------------------------------------------------------------------- #
# MC validation of analytic d' (LLR test) at a chosen budget, one fixed W
# --------------------------------------------------------------------------- #
# validate at a budget where analytic d' ~ 3 (resolvable; a tiny d' is dominated
# by the LLR estimator's own noise floor and cannot be validated).
print("\nMC validation of d' (fixed W, LLR statistic) at budget -> analytic d'~3:")
budget_mc = budget_d3_flat                                 # analytic d' ~ 3
M_per = max(1, int(round(budget_mc / occ_sum)))            # flat windows/(i,l)
W = sc.draw_W(C, op.P, rng)                                # ONE fixed detector
# precompute per-(i,l) fire probs under x and x'
levels = []
for i in range(op.K):
    m = M[i]
    a_lv = sc.default_a_levels((m * x)[(m * x) > 0], L=L)
    fx = sc.cell_fire_prob(a_lv, m * x, W)                 # (L,C)
    fxp = sc.cell_fire_prob(a_lv, m * xp, W)
    levels.append((fx, fxp))

def draw_and_llr(true_f_list):
    """One campaign: draw S ~ PoissonBinomial(f) per (i,l) over M_per windows
    (Gaussian fast), return LLR = sum log p(S|x')/p(S|x) using Gaussian model."""
    llr = 0.0
    for i in range(op.K):
        fx, fxp = levels[i]
        ftrue = true_f_list[i]
        mu_t = ftrue.sum(1); var_t = (ftrue * (1 - ftrue)).sum(1)
        S = rng.normal(mu_t[:, None], np.sqrt(var_t)[:, None], size=(len(mu_t), M_per))
        mux = fx.sum(1); varx = (fx * (1 - fx)).sum(1)
        muxp = fxp.sum(1); varxp = (fxp * (1 - fxp)).sum(1)
        # Gaussian log-density ratio, summed over windows and levels
        lx = -0.5 * ((S - mux[:, None]) ** 2 / varx[:, None] + np.log(varx[:, None]))
        lxp = -0.5 * ((S - muxp[:, None]) ** 2 / varxp[:, None] + np.log(varxp[:, None]))
        llr += float((lxp - lx).sum())
    return llr

NT = 200
llr_x = np.array([draw_and_llr([lv[0] for lv in levels]) for _ in range(NT)])
llr_xp = np.array([draw_and_llr([lv[1] for lv in levels]) for _ in range(NT)])
dprime_mc = (llr_xp.mean() - llr_x.mean()) / np.sqrt(0.5 * (llr_x.var() + llr_xp.var()))
dprime_analytic_mc = np.sqrt(kappa_flat * (M_per * occ_sum))
print("  M/(i,l)=%d  budget=%.2e  analytic d'=%.3f  MC d'=%.3f"
      % (M_per, M_per * occ_sum, dprime_analytic_mc, dprime_mc))
report["mc_validation"] = dict(budget=float(M_per * occ_sum),
    dprime_analytic=float(dprime_analytic_mc), dprime_mc=float(dprime_mc),
    M_per_level=int(M_per), NT=NT)

with open(os.path.join(OUT, "t3_nullspace_break.json"), "w") as fh:
    json.dump(report, fh, indent=2, default=float)

# save scene arrays for plotting
np.savez(os.path.join(OUT, "t3_scenes.npz"),
         x=x.reshape(side, side), xp=xp.reshape(side, side),
         delta=delta.reshape(side, side))
print("\nVERDICT T3: d'>=3 at budget <= linear campaign? %s"
      "  (flat needs %.2e = %.1fx campaign)"
      % ("PASS" if report["PASS"] else "FAIL",
         budget_d3_flat, budget_d3_flat / LIN_CAMPAIGN))
print("wrote t3_nullspace_break.json, t3_scenes.npz")
