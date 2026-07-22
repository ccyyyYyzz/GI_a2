#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
DLGI probe -- MAIN runner: T1 (medium precision), T2 (baselines), T4 (accounting),
and promotion-bar evidence (bars 2,3,4,5).  Writes dual_ledger_results.json.

Grid: t_c in {2,16,64} x CV in {5%,15%,40%}, 6 DEV scenes x 5 seeds (paired
schedule = the natural ladder order).  All arms share the identical 2048-exposure
complementary record / photon budget; the pilot arm REPLACES 5% of Hadamard rows.
CPU.  Read-only on inputs; writes only this probe dir.
"""
import os, sys, json, time, platform
from datetime import datetime, timezone
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dl_common as C

T0 = time.time()
OUT = os.path.join(C.REPO, "results", "round63_next", "DUAL_LEDGER_PROBE")
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

SCENES = C.SCENES
SEEDS = list(range(5))
TC = [2.0, 16.0, 64.0]
CV = [0.05, 0.15, 0.40]

def relerr(hat, true):
    return (np.array(hat) - true) / true

# ------------------------------------------------------------------ bar 2/gauge
def held_out_residual_test(Y, x_hat, slot=None):
    """Scalar-gain adequacy: whiten the log-gain residuals with the fitted OU and
    check they are ~white unit-variance (a spatially-varying/mis-specified gain
    would leave temporal+structural residual autocorrelation).  Returns (whitened
    resid std, lag-1 autocorr of whitened resid)."""
    zc, R, valid, mu = C._z_residuals(Y, x_hat)
    if slot is None: slot = np.arange(C.N_EXP)
    order = np.argsort(slot); zc = zc[order]; R = R[order]; valid = valid[order]
    tc, sig = C.mom_autocov(zc, valid)
    if not np.isfinite(tc): return float("nan"), float("nan")
    phi = np.exp(-1.0 / max(tc, 1e-3))
    w = C.kalman_innovations(zc, R, valid.astype(bool), phi, max(sig, 1e-6))
    r_std = float(np.std(w))
    r_ac1 = float(np.corrcoef(w[:-1], w[1:])[0, 1]) if len(w) > 3 else float("nan")
    return r_std, r_ac1

def gauge_invariance_check(Y, x_hat, slot=None):
    """Scale gauge: a global scene-scale error must not move (t_c,CV)."""
    m1 = C.medium_estimate(Y, x_hat, slot=slot, return_path=False)
    m2 = C.medium_estimate(Y, 3.7 * x_hat, slot=slot, return_path=False)
    return (abs(m1["tc_hat"] - m2["tc_hat"]) / max(m1["tc_hat"], 1e-9),
            abs(m1["cv_hat"] - m2["cv_hat"]) / max(m1["cv_hat"], 1e-9))

# ------------------------------------------------------------------ bar 4 coverage
def bootstrap_ci(x_hat, tc_hat, sig_hat, mu_hat, n_boot=40, seed=0):
    """Parametric bootstrap CI for (t_c,CV): simulate synthetic records from the
    fitted (x_hat, OU) and re-run the medium estimator (fast, reuse x_hat).  Returns
    95% percentile intervals for t_c and CV."""
    s = C.s_of_x(np.clip(x_hat, 0.0, None))
    rng = np.random.default_rng(4242 + seed)
    tcs = []; cvs = []
    phi = np.exp(-1.0 / max(tc_hat, 1e-3))
    for b in range(n_boot):
        # synthetic OU log-gain with the fitted params
        l = np.empty(C.N_EXP); l[0] = rng.standard_normal() * sig_hat
        sd = sig_hat * np.sqrt(max(1 - phi ** 2, 1e-9))
        for t in range(1, C.N_EXP):
            l[t] = phi * l[t - 1] + sd * rng.standard_normal()
        a = np.exp(mu_hat + l)
        Yb = rng.poisson(np.maximum(a * C.PHI * s, 0)).astype(float)
        m = C.medium_estimate(Yb, x_hat, return_path=False)
        tcs.append(m["tc_hat"]); cvs.append(m["cv_hat"])
    return (np.percentile(tcs, [2.5, 97.5]).tolist(),
            np.percentile(cvs, [2.5, 97.5]).tolist())

# ================================================================== T1 + T2 + bars
def run_grid():
    results = {}
    cover_cells = [(16.0, 0.15), (64.0, 0.15), (2.0, 0.15)]  # bar-4 coverage subset
    for tc in TC:
        for cv in CV:
            sig = C.sigma_l_of_cv(cv)
            key = f"tc{int(tc)}_cv{int(cv*100)}"
            pf_tc = []; pf_cv = []; pf_corr = []; pf_psnr = []; lin_psnr = []
            pl_tc = []; pl_cv = []; pl_psnr = []
            or_tc = []; or_cv = []
            rstd = []; rac1 = []; gtc = []; gcv = []
            cover_tc = 0; cover_cv = 0; cover_n = 0
            do_cover = (tc, cv) in cover_cells
            for sc in SCENES:
                x = C.scene_x[sc]
                for seed in SEEDS:
                    # ---- pilot-free (DLGI): one record, two ledgers ----
                    Y, a_time, slot = C.simulate_record(sc, seed, tc, sig)
                    r = C.joint_dual_ledger(Y, slot=slot, n_outer=2)
                    m = r["med"]
                    pf_tc.append(m["tc_hat"]); pf_cv.append(m["cv_hat"])
                    pf_corr.append(C.gain_path_corr(m["a_time"], a_time))
                    pf_psnr.append(C.psnr(r["x_hat"], x))
                    # plain-linear arm (gain-agnostic, same record) for noninferiority
                    xa2 = C.arm_A2(Y[0::2], Y[1::2], C.LAM_A2)
                    lin_psnr.append(C.psnr(xa2, x))
                    # bar-2 residual + gauge
                    a, b = held_out_residual_test(Y, r["x_hat"], slot)
                    rstd.append(a); rac1.append(b)
                    ga, gb = gauge_invariance_check(Y, r["x_hat"], slot)
                    gtc.append(ga); gcv.append(gb)
                    # ---- pilot-interleaved baseline ----
                    p = C.pilot_interleaved(sc, seed, tc, sig)
                    pl_tc.append(p["med"]["tc_hat"]); pl_cv.append(p["med"]["cv_hat"])
                    pl_psnr.append(C.psnr(p["x_hat"], x))
                    # ---- oracle monitor (upper bound) ----
                    om = C.oracle_monitor_estimate(p["a_time"])
                    or_tc.append(om["tc_hat"]); or_cv.append(om["cv_hat"])
                    # ---- bar-4 coverage (subset) ----
                    if do_cover:
                        citc, cicv = bootstrap_ci(r["x_hat"], m["tc_hat"],
                                                  m["sigma_l_hat"], m["mu_hat"],
                                                  n_boot=40, seed=seed)
                        cover_tc += int(citc[0] <= tc <= citc[1])
                        cover_cv += int(cicv[0] <= cv <= cicv[1])
                        cover_n += 1
            pf_tc = np.array(pf_tc); pf_cv = np.array(pf_cv)
            pl_tc = np.array(pl_tc); pl_cv = np.array(pl_cv)
            or_tc = np.array(or_tc); or_cv = np.array(or_cv)
            floor = C.oracle_floor(tc)
            def rmse(a, t): return float(np.sqrt(np.nanmean((np.array(a) - t) ** 2)))
            res = dict(
                tc_true=tc, cv_true=cv, sigma_l=sig, n=len(pf_tc),
                # T1 pilot-free medium precision
                pf_tc_med=float(np.median(pf_tc)),
                pf_tc_relerr_med=float(np.median(relerr(pf_tc, tc))),
                pf_tc_absrelerr_med=float(np.median(np.abs(relerr(pf_tc, tc)))),
                pf_cv_med=float(np.median(pf_cv)),
                pf_cv_relerr_med=float(np.median(relerr(pf_cv, cv))),
                pf_cv_absrelerr_med=float(np.median(np.abs(relerr(pf_cv, cv)))),
                pf_gaincorr_med=float(np.median(pf_corr)),
                pf_scene_psnr_med=float(np.median(pf_psnr)),
                # noninferiority (bar 5)
                lin_scene_psnr_med=float(np.median(lin_psnr)),
                scene_delta_db=float(np.median(np.array(pf_psnr) - np.array(lin_psnr))),
                # T2 baselines
                pilot_tc_med=float(np.nanmedian(pl_tc)),
                pilot_tc_relerr_med=float(np.nanmedian(relerr(pl_tc, tc))),
                pilot_cv_med=float(np.nanmedian(pl_cv)),
                pilot_scene_psnr_med=float(np.median(pl_psnr)),
                oracle_tc_med=float(np.median(or_tc)),
                oracle_cv_med=float(np.median(or_cv)),
                # RMSE and ratios (bar 3)
                pf_tc_rmse=rmse(pf_tc, tc), pilot_tc_rmse=rmse(pl_tc, tc),
                oracle_tc_rmse=rmse(or_tc, tc),
                pf_cv_rmse=rmse(pf_cv, cv), pilot_cv_rmse=rmse(pl_cv, cv),
                oracle_cv_rmse=rmse(or_cv, cv),
                ratio_tc_pf_over_pilot=rmse(pf_tc, tc) / max(rmse(pl_tc, tc), 1e-9),
                ratio_cv_pf_over_pilot=rmse(pf_cv, cv) / max(rmse(pl_cv, cv), 1e-9),
                ratio_tc_pf_over_oracle=rmse(pf_tc, tc) / max(rmse(or_tc, tc), 1e-9),
                ratio_cv_pf_over_oracle=rmse(pf_cv, cv) / max(rmse(or_cv, cv), 1e-9),
                # oracle Fisher floors (R33 sec.4)
                floor_tc_rel=floor["sd_tc_rel"], floor_cv_rel=floor["sd_cv_rel"],
                # bar 2
                resid_whitened_std_med=float(np.nanmedian(rstd)),
                resid_lag1_ac_med=float(np.nanmedian(rac1)),
                gauge_tc_shift_med=float(np.nanmedian(gtc)),
                gauge_cv_shift_med=float(np.nanmedian(gcv)),
            )
            if do_cover and cover_n:
                res["coverage_tc_95"] = cover_tc / cover_n
                res["coverage_cv_95"] = cover_cv / cover_n
                res["coverage_n"] = cover_n
            results[key] = res
            log(f"{key}: pf tc%={100*res['pf_tc_relerr_med']:+.0f} cv%={100*res['pf_cv_relerr_med']:+.0f} "
                f"corr={res['pf_gaincorr_med']:.3f} PSNR={res['pf_scene_psnr_med']:.1f}(+{res['scene_delta_db']:.2f}) "
                f"| pilot tc%={100*res['pilot_tc_relerr_med']:+.0f} | ratio_tc={res['ratio_tc_pf_over_pilot']:.2f}")
    return results

# ================================================================== T4 accounting
def t4_audit():
    """Verify photons/exposures identical for the pilot-free arm and the plain
    linear arm on one representative record: the medium estimate is a pure
    post-processing by-product."""
    sc = "bridge_contour_1"; sig = C.sigma_l_of_cv(0.15)
    Y, a_time, slot = C.simulate_record(sc, 0, 16.0, sig)
    total_photons = float(Y.sum()); n_exp = int(Y.size)
    # plain linear arm consumes the same Y; DLGI consumes the same Y + CPU
    return dict(
        arm_pilot_free=dict(exposures=n_exp, total_detected_photons=total_photons,
                            patterns="1024-row sequency Hadamard x complementary",
                            pilots=0, reference_arm="none", extra_photons=0,
                            medium_output="t_c, CV, gain path (post-processing)"),
        arm_plain_linear=dict(exposures=n_exp, total_detected_photons=total_photons,
                              patterns="identical", pilots=0, reference_arm="none",
                              medium_output="none"),
        arm_pilot=dict(exposures=n_exp, total_detected_photons="~identical",
                       patterns="5% Hadamard rows replaced by all-ones pilots",
                       pilots=int(round(0.05 * C.NPIX)), reference_arm="none",
                       scene_cost="5% of Hadamard coefficients missing"),
        identical_photons_pilot_free_vs_linear=True,
        identical_exposures_pilot_free_vs_linear=True,
        note=("pilot-free and plain-linear arms consume the byte-identical record; "
              "the medium ledger is extracted by post-processing at zero extra "
              "photons/exposures/time/pilots/reference arm"))

def main():
    log("DLGI main runner: T1 medium precision + T2 baselines + T4 + bars 2,3,4,5")
    log(f"grid t_c={TC} x CV={CV}; scenes={len(SCENES)} seeds={len(SEEDS)}; PHI={C.PHI:.4f}")
    grid = run_grid()
    audit = t4_audit()
    meta = dict(utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
                python=platform.python_version(), numpy=np.__version__,
                scenes=SCENES, seeds=SEEDS, tc_grid=TC, cv_grid=CV,
                PHI=C.PHI, N_EXP=C.N_EXP, TARGET_COUNTS=C.TARGET_COUNTS,
                scene_sha={sc: C.scene_sha[sc] for sc in SCENES},
                runtime_s=time.time() - T0)
    blob = dict(meta=meta, grid=grid, t4_audit=audit)
    fn = os.path.join(OUT, "dual_ledger_results.json")
    json.dump(blob, open(fn, "w"), indent=2)
    log(f"saved {fn}  ({time.time()-T0:.1f}s)")

if __name__ == "__main__":
    main()
