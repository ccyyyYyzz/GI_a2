#!/usr/bin/env python
# cbt_law_plane.py -- CALIBRATION PLANE for the CURVED_BLINDNESS_TEST theorem (JOB 2).
#
# Hero-map data for the manuscript: sweep operating points of the curvature-tax law and verify
# (i) the detectability kink sits at the critical curvature-resource budget A* = s^2/r0, and
# (ii) the quartic-in-time leakage coefficient follows
#          C_pred = (||K||_F^2 / 16) * max(s^2 - A r0, 0)^2
# on the leaky (supercritical) branch, collapsing to 0 on the cloaked (subcritical) branch.
#
# This DRIVER imports the FROZEN exam's world + generator functions VERBATIM (no edits to
# curved_blindness_test.py):  build_world, build_geometry, dQ_along, _paired_scores,
# _paired_stats, plus jet_test's kappa_of / kl_exact / chernoff_star and scramble_toy's
# efficient_weight.  Operating point (r0 scale, omega) maps to build_geometry(amp, s_frac):
#   amp   = BASE_AMP * r0_scale   (r0 = ||xi0||_G is linear in amp)
#   s_frac= omega                 (omega = s/r0 = s_frac exactly)
# The grid ratio chi = s^2 / (A_traj r0) sets the trajectory's acceleration budget
#   A_traj = A* / chi   ->   chi>1: A_traj<A* (under-paid, LEAKY = supercritical);
#                            chi<=1: A_traj>=A* (over-paid, CLOAKED = subcritical, C_pred=0).
#
# Divergence normalisations (verified against jet_test):
#   kl_exact(z)      ~ (||K||^2/4)  z^2   (small z)
#   chernoff_star(z) ~ (||K||^2/16) z^2   (small z; = Bhattacharyya)   <-- matches C_pred prefactor
#   with z = DeltaQ(t) ~ t^2 (s^2 - A r0).  So C_fit := median(|chernoff_star(DeltaQ)| / t^4)
#   is the observable that must equal C_pred on the supercritical branch.
#
# The clamped-envelope arm dQ_along(..., "cap_opt", A_cap) realises the max(.,0) behaviour
# (a = -min(A,A*) xi0/r0), so its Chernoff coefficient is 0 (to the parabolic-cloak residual)
# on the subcritical branch and (s^2 - A r0)^2 on the supercritical branch.
#
# MODES:  analytic (default, all operating points, float64, venue-independent);
#         mc --corner K (one of 4 corners, GPU generator, CRN-paired on-orbit AUC).
# WRITE-SCOPE: results/round63_next/CBT_LAW_PLANE/ ONLY.
import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
CBT = os.path.join(ROOT, "results", "round63_next", "CURVED_BLINDNESS_TEST")
SCRAMBLE = os.path.join(ROOT, "results", "round63_next", "SCRAMBLE_EXT")
JET = os.path.join(ROOT, "results", "round63_next", "JET_TEST")
for p in (CBT, SCRAMBLE, JET):
    sys.path.insert(0, p)

# frozen exam functions, imported VERBATIM (importing the module runs no main())
import curved_blindness_test as cbt  # noqa: E402
from curved_blindness_test import (  # noqa: E402
    build_world, build_geometry, dQ_along, _paired_scores, _paired_stats,
)
from jet_test import kappa_of, kl_exact, chernoff_star  # noqa: E402
from scramble_toy import efficient_weight  # noqa: E402

# ----------------------------------------------------------------- grid + frozen bars
CHI_GRID = [0.25, 0.5, 0.8, 1.0, 1.25, 2.0, 4.0]
R0_SCALES = [0.5, 1.0, 2.0]
OMEGAS = [0.7, 1.4]
BASE_AMP = 0.08
GEOM_SEED = 460
T_JET = np.logspace(-2.3, -1.0, 16)          # deep-jet grid, identical to the frozen exam stage [1]

# 4 MC corners = extreme (r0_scale, omega)
MC_CORNERS = [(0.5, 0.7), (0.5, 1.4), (2.0, 0.7), (2.0, 1.4)]

FROZEN = {
    "description": "Calibration plane of the curvature-tax law; bars preregistered before compute.",
    "grid": {"chi": CHI_GRID, "r0_scale": R0_SCALES, "omega": OMEGAS,
             "base_amp": BASE_AMP, "geom_seed": GEOM_SEED,
             "n_operating_points": len(R0_SCALES) * len(OMEGAS),
             "n_grid_points": len(CHI_GRID) * len(R0_SCALES) * len(OMEGAS)},
    "definitions": {
        "A_star": "s^2 / r0 (critical curvature-resource budget)",
        "A_traj": "A_star / chi (trajectory budget at grid ratio chi)",
        "supercritical": "chi > 1  (A_traj < A_star, leaky branch, C_pred > 0)",
        "subcritical": "chi <= 1 (A_traj >= A_star, cloaked branch, C_pred = 0)",
        "C_pred": "(||K||_F^2 / 16) * max(s^2 - A_traj r0, 0)^2",
        "C_fit": "median(|chernoff_star(DeltaQ_capopt)| / t^4) over the deep-jet t grid",
        "C0": "C_fit of the straight arm (A_traj=0) at the operating point (normaliser)",
        "A_kink": "vertex of a quadratic fit to the raw-cap coefficient(A) V-shape",
    },
    "bars": {
        "B1_median_kink_rel_err_lt": 0.02,
        "B2_median_supercritical_C_ratio_in": [0.9, 1.1],
        "B3_all_supercritical_quartic_R2_gt": 0.99,
        "B4_all_subcritical_C_fit_over_C0_lt": 1e-2,
        "MC_paired_auc_on_orbit_in": [0.49, 0.51],
        "MC_B": 4096, "MC_R_records": 3000, "MC_theta": 0.35,
    },
    "note": "A failed bar is a reported result, not a retry license (preregistration).",
}


# ----------------------------------------------------------------- helpers
def loglog_fit_r2(x, y):
    """LS slope and R^2 of log y vs log x over strictly-positive finite points."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = (x > 0) & (y > 0) & np.isfinite(y)
    if int(m.sum()) < 3:
        return float("nan"), float("nan"), int(m.sum())
    lx, ly = np.log(x[m]), np.log(y[m])
    A = np.vstack([lx, np.ones_like(lx)]).T
    coef, _, _, _ = np.linalg.lstsq(A, ly, rcond=None)
    yhat = A @ coef
    ss_res = float(np.sum((ly - yhat) ** 2))
    ss_tot = float(np.sum((ly - ly.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(coef[0]), float(r2), int(m.sum())


def make_geom(world, r0_scale, omega):
    return build_geometry(world, amp=BASE_AMP * r0_scale, s_frac=omega, seed=GEOM_SEED)


# ----------------------------------------------------------------- analytic arm (one operating point)
def analytic_point(world, r0_scale, omega):
    geom = make_geom(world, r0_scale, omega)
    r0, s, A_star = geom["r0"], geom["s"], geom["A_star"]
    Sig0 = geom["Q0"] * world["H"] + geom["R"]
    kap = kappa_of(Sig0, world["H"])
    K_F2 = float(np.sum(kap ** 2))                       # ||K||_F^2
    t = T_JET

    # --- A_kink: fine raw-cap coefficient V-shape, vertex of a quadratic fit (exact for a parabola)
    Agrid = np.linspace(0.3 * A_star, 1.7 * A_star, 61)
    cap_coef = np.empty(len(Agrid))
    for i, Ac in enumerate(Agrid):
        dq = np.array([dQ_along(world, geom, "cap", tt, A_cap=Ac) for tt in t])
        kl = np.array([kl_exact(z, kap) for z in dq])
        cap_coef[i] = float(np.median(np.abs(kl) / t ** 4))
    p = np.polyfit(Agrid, cap_coef, 2)                   # coef(A) ~ quadratic in A
    A_kink = float(-p[1] / (2.0 * p[0]))
    kink_rel_err = float(abs(A_kink - A_star) / A_star)

    # --- straight-arm normaliser C0 (Chernoff coefficient at A_traj = 0)
    dq0 = np.array([dQ_along(world, geom, "straight", tt) for tt in t])
    ch0 = np.array([chernoff_star(z, kap)[0] for z in dq0])
    C0 = float(np.median(np.abs(ch0) / t ** 4))

    # --- per-chi C_fit on the clamped-envelope (cap_opt) arm
    chi_rows = []
    for chi in CHI_GRID:
        A_traj = A_star / chi
        dq = np.array([dQ_along(world, geom, "cap_opt", tt, A_cap=A_traj) for tt in t])
        ch = np.array([chernoff_star(z, kap)[0] for z in dq])
        slope, r2, npts = loglog_fit_r2(t, np.abs(ch))
        C_fit = float(np.median(np.abs(ch) / t ** 4))
        margin = s ** 2 - A_traj * r0
        C_pred = float((K_F2 / 16.0) * max(margin, 0.0) ** 2)
        supercritical = bool(chi > 1.0)
        chi_rows.append(dict(
            chi=float(chi), A_traj=float(A_traj), A_over_Astar=float(A_traj / A_star),
            margin_s2_minus_Ar0=float(margin), supercritical=supercritical,
            C_fit=C_fit, C_pred=C_pred,
            C_ratio=(float(C_fit / C_pred) if C_pred > 0 else None),
            C_fit_over_C0=(float(C_fit / C0) if C0 > 0 else None),
            quartic_slope=slope, quartic_R2=r2, n_fit_pts=npts,
            dQ=dq.tolist(), chernoff=ch.tolist()))

    return dict(
        r0_scale=float(r0_scale), omega=float(omega), r0=float(r0), s=float(s),
        A_star=float(A_star), omega_check=float(geom["omega"]), Q0=float(geom["Q0"]),
        K_F2=K_F2, C0=C0, A_kink=A_kink, kink_rel_err=kink_rel_err,
        A_grid=Agrid.tolist(), cap_coef=cap_coef.tolist(),
        t=t.tolist(), chi_rows=chi_rows)


def run_analytic(world):
    pts = []
    for r0s in R0_SCALES:
        for om in OMEGAS:
            t0 = time.time()
            pt = analytic_point(world, r0s, om)
            pt["runtime_sec"] = round(time.time() - t0, 2)
            pts.append(pt)
            print("  [analytic] r0_scale=%.2f omega=%.2f | A*=%.4f A_kink=%.4f (rel %.2e) "
                  "| C0=%.3e | %.1fs" % (r0s, om, pt["A_star"], pt["A_kink"],
                                         pt["kink_rel_err"], pt["C0"], pt["runtime_sec"]),
                  flush=True)
    return pts


def evaluate_bars(pts):
    kink_errs = [p["kink_rel_err"] for p in pts]
    B1_val = float(np.median(kink_errs))
    B1 = bool(B1_val < FROZEN["bars"]["B1_median_kink_rel_err_lt"])

    super_ratios = [r["C_ratio"] for p in pts for r in p["chi_rows"]
                    if r["supercritical"] and r["C_ratio"] is not None]
    B2_val = float(np.median(super_ratios)) if super_ratios else float("nan")
    lo, hi = FROZEN["bars"]["B2_median_supercritical_C_ratio_in"]
    B2 = bool(lo <= B2_val <= hi)

    super_r2 = [r["quartic_R2"] for p in pts for r in p["chi_rows"] if r["supercritical"]]
    B3_min = float(np.min(super_r2)) if super_r2 else float("nan")
    B3 = bool(all(r2 > FROZEN["bars"]["B3_all_supercritical_quartic_R2_gt"] for r2 in super_r2))

    sub_ratios = [r["C_fit_over_C0"] for p in pts for r in p["chi_rows"]
                  if (not r["supercritical"]) and r["C_fit_over_C0"] is not None]
    B4_max = float(np.max(sub_ratios)) if sub_ratios else float("nan")
    B4 = bool(all(v < FROZEN["bars"]["B4_all_subcritical_C_fit_over_C0_lt"] for v in sub_ratios))

    return {
        "B1_median_kink_rel_err": {"value": B1_val, "pass": B1},
        "B2_median_supercritical_C_ratio": {"value": B2_val, "n": len(super_ratios), "pass": B2},
        "B3_min_supercritical_quartic_R2": {"value": B3_min, "n": len(super_r2), "pass": B3},
        "B4_max_subcritical_C_fit_over_C0": {"value": B4_max, "n": len(sub_ratios), "pass": B4},
    }


# ----------------------------------------------------------------- MC corner (on-orbit AUC)
def run_mc_corner(world, corner_idx, B=None, R=None, theta=None):
    r0s, om = MC_CORNERS[corner_idx]
    B = B or FROZEN["bars"]["MC_B"]
    R = R or FROZEN["bars"]["MC_R_records"]
    theta = theta or FROZEN["bars"]["MC_theta"]
    geom = make_geom(world, r0s, om)
    W = efficient_weight(world["H"], geom["Q0"], geom["shot_diag"])
    w = geom["omega"]
    tt = theta / w

    def gc_scene(t_):
        xi_t = geom["xi0"] * np.cos(w * t_) + (geom["v"] / w) * np.sin(w * t_)
        return geom["B0"] + xi_t

    seed0 = 63000 + 1000 * corner_idx
    t0 = time.time()
    sc = _paired_scores(world, geom, W, geom["scene0"], gc_scene(tt), B, [B], R, seed0=seed0)
    s0, s1 = sc[B]
    stt = _paired_stats(s0, s1)
    band = FROZEN["bars"]["MC_paired_auc_on_orbit_in"]
    return dict(
        corner_idx=corner_idx, r0_scale=r0s, omega=om, B=B, R_records=R, theta=theta,
        t_on_orbit=float(tt), dQ_exact=float(dQ_along(world, geom, "exact", tt)),
        paired_auc=stt["paired_auc"], paired_dprime=stt["dprime_paired"],
        se_dprime=stt["se_dprime"], unpaired_auc=stt["unpaired_auc"],
        unpaired_auc_ci95=stt["unpaired_auc_ci95"],
        paired_auc_in_band=bool(band[0] <= stt["paired_auc"] <= band[1]),
        dprime_within_2sigma_of_0=stt["dprime_within_2sigma_of_0"],
        seed0=seed0, runtime_sec=round(time.time() - t0, 1))


# ----------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["analytic", "mc"], default="analytic")
    ap.add_argument("--corner", type=int, default=-1, help="MC corner index 0..3 (mc mode)")
    ap.add_argument("--mc-R", type=int, default=0, help="override MC records (smoke)")
    ap.add_argument("--mc-B", type=int, default=0, help="override MC banks (smoke)")
    ap.add_argument("--out", default="", help="output json path (default auto)")
    args = ap.parse_args()

    import scramble_toy as st
    print("[cbt_law_plane] mode=%s device=%s" % (args.mode, st.DEV), flush=True)
    world = build_world()

    if args.mode == "analytic":
        t0 = time.time()
        pts = run_analytic(world)
        bars = evaluate_bars(pts)
        payload = dict(frozen=FROZEN, kind="analytic",
                       config=dict(N_SIDE=cbt.N_SIDE, N=cbt.N, M=cbt.M, PHOT=cbt.PHOT),
                       operating_points=pts, bars=bars,
                       runtime_sec=round(time.time() - t0, 1))
        out = args.out or os.path.join(HERE, "CBT_LAW_PLANE_analytic.json")
        with open(out, "w") as f:
            json.dump(payload, f, indent=2, default=float)
        print("[cbt_law_plane] BARS: " + " ".join(
            "%s=%s" % (k, v["pass"]) for k, v in bars.items()), flush=True)
        print("[cbt_law_plane] wrote %s (%.1fs)" % (out, payload["runtime_sec"]), flush=True)
    else:
        assert 0 <= args.corner < len(MC_CORNERS), "need --corner 0..3"
        res = run_mc_corner(world, args.corner,
                            B=(args.mc_B or None), R=(args.mc_R or None))
        payload = dict(frozen=FROZEN, kind="mc_corner", corner=res)
        out = args.out or os.path.join(HERE, "CBT_LAW_PLANE_mc_corner%d.json" % args.corner)
        with open(out, "w") as f:
            json.dump(payload, f, indent=2, default=float)
        print("[cbt_law_plane] corner %d (r0=%.1f om=%.1f): paired_auc=%.4f in_band=%s "
              "dQ_exact=%.2e (%.1fs) -> %s" % (
                  args.corner, res["r0_scale"], res["omega"], res["paired_auc"],
                  res["paired_auc_in_band"], res["dQ_exact"], res["runtime_sec"], out),
              flush=True)


if __name__ == "__main__":
    main()
