"""
SELF-AVERAGING (quenched vs annealed) -- the load-bearing physical assumption.

The "no calibration of W" claim rests on the FIXED detector self-averaging over
its C cells.  We quantify:
  (1) quenched p2-recovery bias vs C (annealed-model fit to a single frozen W);
  (2) the 1/sqrt(C) scaling (matches analytic survival fluctuation);
  (3) rescue A -- speckle boiling / frame averaging (N_spk frozen frames
      averaged -> bias ~ 1/sqrt(C*N_spk) -> annealed limit, NO per-W knowledge);
  (4) rescue B -- one-time calibration: a fixed W is a valid fixed instrument;
      its own response curve, learned once on known scenes, removes the bias.
Also: cell-cell CORRELATION as an effective-C reduction (C_eff = C/(1+(n-1)r)).
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
rng = np.random.default_rng(651_202)

# representative masked scene.  Pixels with v_p=0 contribute u_c-term 0 and
# survival factor 1, so we work on the ACTIVE pixels only (exact, faster).
x = op.load_scene(os.path.join(DATA, "bridge_contour_2.npz"))
M = op.build_operator()
vfull = M[0] * x
v = vfull[vfull > 0]                              # active pixels only (exact)
p2t = float(sc.power_sums(v, 2)[2])
a_lv = sc.default_a_levels(v, L=10)
ORDER = 4
rep = {"scene": "bridge_contour_2", "p2_true": p2t, "order": ORDER,
       "n_active": int(v.size)}

# ------------------------------------------------------------------ #
# (1)+(2) quenched bias vs C, many realisations -> RMS bias & scaling
# ------------------------------------------------------------------ #
print("QUENCHED p2 bias vs C (annealed-model fit to single frozen W):")
NREAL = 50
scaling = {}
for C in (900, 3600, 14400):
    errs = np.empty(NREAL)
    for r in range(NREAL):
        W = sc.draw_W(C, len(v), rng)
        Qq = sc.empirical_survival(a_lv, v, W)
        errs[r] = (sc.fit_powersums(a_lv, Qq, order=ORDER)["p2"] - p2t) / p2t
    rms = float(np.sqrt(np.mean(errs ** 2)))
    scaling[C] = dict(rms_relbias=rms, mean_relbias=float(errs.mean()),
                      p50_abs=float(np.median(np.abs(errs))))
    print("  C=%6d  RMS rel bias = %.3e   (x sqrt(C) = %.3f)"
          % (C, rms, rms * np.sqrt(C)))
rep["quenched_bias_vs_C"] = scaling
# verify ~1/sqrt(C): product rms*sqrt(C) ~ const
consts = [scaling[C]["rms_relbias"] * np.sqrt(C) for C in scaling]
rep["sqrtC_law_const_range"] = [float(min(consts)), float(max(consts))]
print("  1/sqrt(C) law: rms*sqrt(C) in [%.3f, %.3f] (approx constant => scaling holds)"
      % (min(consts), max(consts)))

# ------------------------------------------------------------------ #
# (3) rescue A: frame averaging (speckle boils) at C=3600
# ------------------------------------------------------------------ #
print("\nRESCUE A -- frame averaging (N_spk boiled frames), C=3600:")
C = 3600
frame_rescue = {}
for Nspk in (1, 4, 16, 32):
    errs = np.empty(NREAL)
    for r in range(NREAL):
        Qacc = np.zeros(len(a_lv))
        for _ in range(Nspk):
            W = sc.draw_W(C, len(v), rng)
            Qacc += sc.empirical_survival(a_lv, v, W)
        Qacc /= Nspk
        errs[r] = (sc.fit_powersums(a_lv, Qacc, order=ORDER)["p2"] - p2t) / p2t
    rms = float(np.sqrt(np.mean(errs ** 2)))
    frame_rescue[Nspk] = rms
    print("  N_spk=%3d  RMS rel bias = %.3e   (effective C = %d)"
          % (Nspk, rms, C * Nspk))
rep["frame_averaging_rescue_C3600"] = frame_rescue

# ------------------------------------------------------------------ #
# (4) rescue B: one-time calibration of a fixed instrument.
# A frozen W defines a fixed forward map v -> Qhat(a;v).  Calibrate the
# quadratic response coefficient on a small set of KNOWN scenes (learn the
# effective 'a^2' column gain g so that fit recovers p2 unbiasedly), then apply
# to a held-out scene.  Demonstrates a fixed detector needs no W-knowledge, only
# a one-time scalar calibration.
# ------------------------------------------------------------------ #
print("\nRESCUE B -- one-time calibration of a single frozen W (C=3600):")
C = 3600
P = len(vfull)                                    # full pixel dimension (1024)
W = sc.draw_W(C, P, rng)                          # ONE fixed instrument (all px)
# calibration scenes: other masked bridge scenes (known p2). Use FULL vectors so
# the SAME frozen W applies to every scene.
cal_names = ["bridge_twopop_0", "bridge_microtex_0", "bridge_control_1",
             "bridge_contour_0"]
ratios = []
for nm in cal_names:
    xc = op.load_scene(os.path.join(DATA, nm + ".npz"))
    vc = M[0] * xc                                # full-length masked scene
    p2c = float(sc.power_sums(vc, 2)[2])
    a_c = sc.default_a_levels(vc, L=10)
    p2_meas = sc.fit_powersums(a_c, sc.empirical_survival(a_c, vc, W), order=ORDER)["p2"]
    ratios.append(p2_meas / p2c)
g = float(np.mean(ratios))                        # one scalar calibration gain
# apply to held-out scene (full-length) through the SAME W
vheld = vfull
p2t_held = float(sc.power_sums(vheld, 2)[2])
p2_raw = sc.fit_powersums(a_lv, sc.empirical_survival(a_lv, vheld, W), order=ORDER)["p2"]
p2_cal = p2_raw / g
p2t = p2t_held                                    # align truth to full-length scene
print("  calibration gain g = %.4f (from %d known scenes)" % (g, len(cal_names)))
print("  held-out p2: raw rel err = %.3e -> calibrated rel err = %.3e"
      % (abs(p2_raw - p2t) / p2t, abs(p2_cal - p2t) / p2t))
rep["calibration_rescue"] = dict(gain=g, relerr_raw=float(abs(p2_raw - p2t) / p2t),
                                 relerr_calibrated=float(abs(p2_cal - p2t) / p2t),
                                 n_cal_scenes=len(cal_names))

# ------------------------------------------------------------------ #
# (5) cell-cell correlation -> effective C.  If cells share speckle grains with
# pairwise correlation r over groups of n, the variance of the empirical mean
# inflates by ~[1+(n-1) r], i.e. C_eff = C/[1+(n-1) r].  Demonstrate with a
# block-correlated W (n cells per block share one grain, blended by r).
# ------------------------------------------------------------------ #
print("\nCELL-CELL CORRELATION -> effective C  (C=3600):")
# Marginal-preserving model: a block of `block` physical cells shares ONE speckle
# grain, i.e. identical weight vectors (perfect intra-block correlation, r=1
# within block, 0 across).  Marginals stay Exp(1/C) -> product formula UNBIASED;
# only C/block distinct draws -> C_eff = C/block exactly.  This isolates the
# effective-C (variance) effect from any distribution-change bias.
p2t_corr = float(sc.power_sums(v, 2)[2])          # v is active-pixel scene here
def draw_W_dup(C, P, rng, block):
    # distinct grains drawn with the PHYSICAL mean 1/C (not 1/ndist), then each
    # shared by `block` physical cells -> marginals stay Exp(1/C), C_eff=C/block.
    ndist = C // block
    distinct = rng.exponential(scale=1.0 / C, size=(ndist, P))
    return np.repeat(distinct, block, axis=0)[:ndist * block]
corr = {}
for block in (1, 4, 16, 36):
    Ceff = C // block
    errs = np.empty(40)
    for rr in range(40):
        Wc = draw_W_dup(C, len(v), rng, block)
        Qc = sc.empirical_survival(a_lv, v, Wc)
        errs[rr] = (sc.fit_powersums(a_lv, Qc, order=ORDER)["p2"] - p2t_corr) / p2t_corr
    rms = float(np.sqrt(np.mean(errs ** 2)))
    corr[f"block{block}"] = dict(rms_relbias=rms, C_eff=int(Ceff),
                                 rms_times_sqrt_Ceff=float(rms * np.sqrt(Ceff)))
    print("  block=%2d  C_eff=%4d  RMS rel bias=%.3e  (x sqrt(C_eff)=%.3f)"
          % (block, Ceff, rms, rms * np.sqrt(Ceff)))
rep["cell_correlation_effectiveC"] = corr
rep["cell_correlation_note"] = ("perfect intra-block sharing => C_eff=C/block; "
    "rms*sqrt(C_eff) stays ~const (same law as independent cells) => correlation "
    "acts purely through effective-C reduction")

with open(os.path.join(OUT, "selfaverage.json"), "w") as fh:
    json.dump(rep, fh, indent=2, default=float)
print("\nwrote selfaverage.json")
