#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
EXPLORATORY / DEV-ONLY.  NOT preregistered.  No commits.  Local CPU/GPU only.

Configurable COPY of the frozen DLGI estimator (dl_common.medium_estimate /
joint_dual_ledger) for the DL_TOOL_EXPLORATORY experiments.  dl_common is imported
READ-ONLY; every variant here is a copy with switches, so the frozen baseline is
never mutated.

The one thing under test in E1 is the gain-LEVEL / scale gauge (a Jensen/gauge
term at high CV that the head-to-head flagged: DLGI loses the fast-strong-drift
corner on scene PSNR despite equal gain-path CORRELATION -> a multiplicative level
error, not a path-shape error).  Switches:

  gauge   : 'counts'  (frozen: weighted mean of z with w=1/R=Y+0.5)
            'uniform' (unweighted mean of z over valid exposures)
  jensen  : True  (frozen: a = exp(mu + m_sm + 0.5 V_sm), posterior MEAN of a)
            False (plain plug-in a = exp(mu + m_sm), posterior MEDIAN of a)
  renorm  : 'none'  (frozen)
            'arith' (impose the KNOWN physics E[a]=1: a <- a / mean(a))
            'geom'  (impose geometric-mean(a)=1)

Everything else (mom_autocov OU fit, RTS smoother, reconstruction arm_A4+TV) is
reused verbatim from dl_common / part1_gain_ladder.
"""
import os, sys
import numpy as np

REPO = "D:/GI_another"
sys.path.insert(0, os.path.join(REPO, "results", "round63_next", "DUAL_LEDGER_PROBE", "code"))
import dl_common as DL

# re-export frozen objects
NPIX, SIDE, N_EXP, PHI = DL.NPIX, DL.SIDE, DL.N_EXP, DL.PHI
Mp, Mm = DL.Mp, DL.Mm
s_of_x = DL.s_of_x
mom_autocov = DL.mom_autocov
kalman_rts = DL.kalman_rts
cv_of_sigma_l = DL.cv_of_sigma_l
sigma_l_of_cv = DL.sigma_l_of_cv
arm_A2, arm_A4_gaincorr = DL.arm_A2, DL.arm_A4_gaincorr
LAM_A2, LAM_TV = DL.LAM_A2, DL.LAM_TV
ou_path = DL.ou_path
scene_x, SCENES = DL.scene_x, DL.SCENES
simulate_record = DL.simulate_record
gain_path_corr = DL.gain_path_corr
psnr = DL.psnr


def _z_residuals_cfg(Y, x_hat, gauge="counts", s_floor_frac=0.02):
    """Copy of dl_common._z_residuals with a switchable scene-scale gauge.
    gauge='counts' reproduces the frozen weighted mean (w=1/R=Y+0.5); 'uniform'
    uses the unweighted mean over valid exposures."""
    s = s_of_x(np.clip(x_hat, 0.0, None))
    s_floor = s_floor_frac * np.median(s[s > 0])
    valid = (s > s_floor).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.log((Y + 0.5) / np.maximum(PHI * s, 1e-12))
    R = 1.0 / (Y + 0.5)
    z = np.where(valid > 0, z, 0.0)
    if gauge == "counts":
        w = valid / (R + 1e-12)               # w ~ (Y+0.5): counts-weighted (frozen)
    elif gauge == "uniform":
        w = valid.copy()
    else:
        raise ValueError(gauge)
    mu_hat = float((w * z).sum() / max(w.sum(), 1e-30))
    return z - mu_hat, R, valid, mu_hat


def medium_estimate_cfg(Y, x_hat, slot=None, return_path=True,
                        gauge="counts", jensen=True, renorm="none"):
    """Configurable copy of dl_common.medium_estimate.  Returns the same dict; the
    a_time / a_hat gain path is built with the chosen gauge / Jensen / renorm."""
    zc_e, R_e, valid_e, mu_hat = _z_residuals_cfg(Y, x_hat, gauge=gauge)
    if slot is None:
        slot = np.arange(N_EXP)
    order = np.argsort(slot)
    zc = zc_e[order]; R = R_e[order]; valid = valid_e[order]
    tc_hat, sig_hat = mom_autocov(zc, valid)
    if not np.isfinite(tc_hat):
        tc_hat, sig_hat = 1e3, 1e-3
    out = dict(tc_hat=float(tc_hat), sigma_l_hat=float(sig_hat),
               cv_hat=cv_of_sigma_l(sig_hat), mu_hat=mu_hat,
               n_valid=int(valid.sum()))
    if return_path:
        phi = np.exp(-1.0 / max(tc_hat, 1e-3))
        ms, Vs = kalman_rts(zc, R, valid.astype(bool), phi, max(sig_hat, 1e-6))
        jterm = 0.5 * Vs if jensen else 0.0
        a_time = np.exp(mu_hat + ms + jterm)
        if renorm == "arith":
            a_time = a_time / max(float(np.mean(a_time)), 1e-30)   # impose E[a]=1
        elif renorm == "geom":
            a_time = a_time / float(np.exp(np.mean(np.log(np.maximum(a_time, 1e-30)))))
        elif renorm != "none":
            raise ValueError(renorm)
        out["a_time"] = a_time
        a_exp = np.empty(N_EXP); a_exp[order] = a_time
        out["a_hat"] = a_exp
    return out


def joint_dual_ledger_cfg(Y, slot=None, n_outer=2,
                          gauge="counts", jensen=True, renorm="none"):
    """Configurable copy of dl_common.joint_dual_ledger."""
    Yp = Y[0::2]; Ym = Y[1::2]
    x_hat = np.clip(arm_A2(Yp, Ym, LAM_A2), 1e-9, None)
    med = medium_estimate_cfg(Y, x_hat, slot=slot, return_path=True,
                              gauge=gauge, jensen=jensen, renorm=renorm)
    for _ in range(max(n_outer - 1, 0)):
        a_hat = med["a_hat"]
        x_hat = np.clip(arm_A4_gaincorr(Yp, Ym, a_hat[0::2], a_hat[1::2], LAM_TV), 0.0, None)
        med = medium_estimate_cfg(Y, x_hat, slot=slot, return_path=True,
                                  gauge=gauge, jensen=jensen, renorm=renorm)
    return dict(x_hat=x_hat, med=med)


def gain_pieces(Y, x_hat, slot=None, gauge="counts"):
    """Expose the smoother internals in TIME order for the E1(c) residual corrector:
    zc (demeaned log-gain obs), R (obs var), valid, mu_hat, ms (RTS-smoothed demeaned
    log-gain = l_hat), Vs (posterior var), tc_hat, sig_hat, order (exposure->time)."""
    zc_e, R_e, valid_e, mu_hat = _z_residuals_cfg(Y, x_hat, gauge=gauge)
    if slot is None:
        slot = np.arange(N_EXP)
    order = np.argsort(slot)
    zc = zc_e[order]; R = R_e[order]; valid = valid_e[order]
    tc_hat, sig_hat = mom_autocov(zc, valid)
    if not np.isfinite(tc_hat):
        tc_hat, sig_hat = 1e3, 1e-3
    phi = np.exp(-1.0 / max(tc_hat, 1e-3))
    ms, Vs = kalman_rts(zc, R, valid.astype(bool), phi, max(sig_hat, 1e-6))
    return dict(zc=zc, R=R, valid=valid, mu_hat=mu_hat, ms=ms, Vs=Vs,
                tc_hat=float(tc_hat), sig_hat=float(sig_hat), order=order)


def a_hat_from_l(ms_time, Vs_time, mu_hat, order, jensen=True, renorm="arith"):
    """Rebuild the per-exposure gain a_hat from a (possibly corrected) smoothed
    demeaned log-gain path ms_time in TIME order, applying Jensen + E[a]=1 renorm the
    same way the estimator does.  Returns a_exp (per exposure)."""
    jterm = 0.5 * Vs_time if jensen else 0.0
    a_time = np.exp(mu_hat + ms_time + jterm)
    if renorm == "arith":
        a_time = a_time / max(float(np.mean(a_time)), 1e-30)
    elif renorm == "geom":
        a_time = a_time / float(np.exp(np.mean(np.log(np.maximum(a_time, 1e-30)))))
    a_exp = np.empty(N_EXP); a_exp[order] = a_time
    return a_exp


def best_scale_psnr(x_hat, x):
    """PSNR after the optimal single scalar rescale kappa* = <x_hat,x>/<x_hat,x_hat>.
    Isolates the SHAPE-limited PSNR; the gap to the raw PSNR is the pure LEVEL/scale
    loss."""
    denom = float(x_hat @ x_hat)
    if denom <= 0:
        return DL.psnr(x_hat, x), 1.0
    k = float(x_hat @ x) / denom
    return DL.psnr(k * x_hat, x), k
