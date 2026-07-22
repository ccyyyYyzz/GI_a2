#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
DUAL-LEDGER IMAGING probe -- shared machinery.
ROUND63-NEXT / GI_a2.  CPU.  Read-only on inputs; writers live in sibling files.

THE CLAIM UNDER TEST.  One standard GI bucket record through a dynamic medium is
simultaneously (1) a camera for the scene x and (2) an instrument measuring the
MEDIUM's dynamics -- correlation time t_c, fluctuation strength CV(a), and the gain
path a_t itself -- at ZERO extra photons / time / pilots / reference arm.  The HSGI
identity  I_scene(unknown a) = I_scene(known a) - M^T E[Cov(a|Y)] M  makes the term
that limits the scene channel identical to the medium channel's signal.

FORWARD MODEL (lifted VERBATIM from
results/round63_next/LADDER_PROBE/part1_gain_ladder.py, itself verbatim from the
frozen CPL probe): 32x32 NPIX=1024; 6 DEV bridge scenes; full 1024-row sequency
Hadamard as 2048 complementary (m+,m-) exposures; PHI matched to 2200 counts/exp;
frozen within-exposure-constant OU gain a_e (model M0, one gain per exposure).  We
use M0 because the MEDIUM channel is about exposure-to-exposure drift.

MEDIUM ESTIMATOR (the new machinery; the ladder smoother ASSUMED (t_c,CV) known --
here they are UNKNOWN and estimated).  Given the record Y and a scene estimate
x_hat:
  z_e = log((Y_e+0.5)/(PHI * s_e)),   s_e = pattern_e . x_hat,   e with s_e>0
is a noisy observation of the log-gain l_e = log a_e, with obs variance
R_e ~ 1/(Y_e+0.5) (delta method).  z_e is thus a length-2048 noisy sample of the
OU log-gain path.  The scene-scale gauge (a constant multiplicative gain is
degenerate with object amplitude) enters z only as a CONSTANT OFFSET, so we PROFILE
IT OUT by demeaning z; the medium parameters (t_c, CV) live entirely in the
FLUCTUATIONS and are gauge-free.  We then fit a stationary zero-mean OU
(phi=exp(-1/t_c), stationary sd sigma_l) with per-sample obs noise R_e by MAXIMISING
THE KALMAN MARGINAL LIKELIHOOD of z (the gain path integrated out) over (t_c,
sigma_l).  CV(a)=sqrt(exp(sigma_l^2)-1).  A method-of-moments autocovariance fit
(which deconvolves the white Poisson noise: gamma_z(h>=1)=sigma_l^2 phi^h,
gamma_z(0)=sigma_l^2+E[R]) is carried as an independent cross-check.
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = "D:/GI_another"
sys.path.insert(0, os.path.join(REPO, "results", "round63_next", "LADDER_PROBE"))
sys.path.insert(0, os.path.join(REPO, "code", "scatter"))
import part1_gain_ladder as L        # frozen forward model (verbatim)

# ---- re-export the frozen objects we reuse -------------------------------
NPIX, SIDE, N_EXP = L.NPIX, L.SIDE, L.N_EXP
PHI = L.PHI
Mp, Mm, MpT, MmT, HT, H = L.Mp, L.Mm, L.MpT, L.MmT, L.HT, L.H
scene_x, scene_sha, SCENES = L.scene_x, L.scene_sha, L.SCENES
ou_path = L.ou_path
arm_A2, arm_A3, arm_A4_gaincorr = L.arm_A2, L.arm_A3, L.arm_A4_gaincorr
rof_denoise, hadamard_inverse, psnr = L.rof_denoise, L.hadamard_inverse, L.psnr
LAM_TV, LAM_A2, LAM_A3 = L.LAM_TV, L.LAM_A2, L.LAM_A3
TARGET_COUNTS = L.TARGET_COUNTS

# grid: cells declared by TARGET CV(a); sigma_l derived so E[a]=1 exactly.
TC_GRID = [2.0, 16.0, 64.0]
CV_GRID = [0.05, 0.15, 0.40]
def sigma_l_of_cv(cv):
    """log-gain sd giving lognormal CV(a)=cv with E[a]=1 (mu=-sigma^2/2)."""
    return float(np.sqrt(np.log(1.0 + cv * cv)))
def cv_of_sigma_l(sig):
    return float(np.sqrt(np.expm1(sig * sig)))

# ---- exposure -> pattern projection s_e = pattern_e . x ------------------
def s_of_x(x):
    """s_e for all 2048 exposures given image x.  even e -> Mp[e//2], odd -> Mm."""
    s = np.empty(N_EXP)
    s[0::2] = Mp @ x
    s[1::2] = Mm @ x
    return s

# ==========================================================================
# Kalman filter / RTS smoother for a stationary zero-mean scalar OU state
# observed with per-sample Gaussian noise; missing obs allowed.
#   l_t = phi l_{t-1} + w,  w~N(0,Q),  Q=sig^2(1-phi^2),  l_1~N(0,sig^2)
#   z_t = l_t + v_t,  v~N(0,R_t)   (valid_t marks a usable observation)
# ==========================================================================
def kalman_loglik(z, R, valid, phi, sig):
    """Exact Gaussian marginal log-likelihood of the observed z (gain path
    integrated out).  O(N).  Used as the MLE objective for (t_c=-1/log phi, sig)."""
    N = z.shape[0]
    Q = sig * sig * (1.0 - phi * phi)
    a = 0.0; P = sig * sig                 # stationary prior N(0, sig^2)
    ll = 0.0
    for t in range(N):
        if t > 0:
            a = phi * a
            P = phi * phi * P + Q
        if valid[t]:
            S = P + R[t]
            innov = z[t] - a
            ll += -0.5 * (np.log(2.0 * np.pi * S) + innov * innov / S)
            K = P / S
            a = a + K * innov
            P = (1.0 - K) * P
    return ll

def kalman_innovations(z, R, valid, phi, sig):
    """Standardized one-step-ahead (held-out) filter innovations nu_t/sqrt(S_t).
    Under the correct scalar-OU-gain model these are iid N(0,1); departures (std!=1
    or nonzero autocorrelation) flag a mis-specified gain model.  Returns the array
    of standardized innovations at valid observations."""
    N = z.shape[0]
    Q = sig * sig * (1.0 - phi * phi)
    a = 0.0; P = sig * sig
    out = []
    for t in range(N):
        if t > 0:
            a = phi * a; P = phi * phi * P + Q
        if valid[t]:
            S = P + R[t]; innov = z[t] - a
            out.append(innov / np.sqrt(S))
            K = P / S; a = a + K * innov; P = (1.0 - K) * P
    return np.array(out)

def kalman_rts(z, R, valid, phi, sig):
    """RTS smoother; returns (m_sm, V_sm) smoothed mean/var of l_t."""
    N = z.shape[0]
    Q = sig * sig * (1.0 - phi * phi)
    ap = np.empty(N); Pp = np.empty(N); af = np.empty(N); Pf = np.empty(N)
    ap[0] = 0.0; Pp[0] = sig * sig
    for t in range(N):
        if t > 0:
            ap[t] = phi * af[t - 1]
            Pp[t] = phi * phi * Pf[t - 1] + Q
        if valid[t]:
            S = Pp[t] + R[t]; K = Pp[t] / S
            af[t] = ap[t] + K * (z[t] - ap[t])
            Pf[t] = (1.0 - K) * Pp[t]
        else:
            af[t] = ap[t]; Pf[t] = Pp[t]
    ms = np.empty(N); Vs = np.empty(N)
    ms[-1] = af[-1]; Vs[-1] = Pf[-1]
    for t in range(N - 2, -1, -1):
        C = phi * Pf[t] / max(Pp[t + 1], 1e-30)
        ms[t] = af[t] + C * (ms[t + 1] - phi * af[t])
        Vs[t] = Pf[t] + C * C * (Vs[t + 1] - Pp[t + 1])
    return ms, np.maximum(Vs, 0.0)

# ---- medium-parameter MLE over (t_c, sigma_l) ----------------------------
from scipy.optimize import minimize_scalar, minimize

def _profile_mean(z, R, valid):
    """weighted mean of valid z (removes the scene-scale gauge = constant offset)."""
    w = valid / (R + 1e-12)
    return float((w * z).sum() / max(w.sum(), 1e-30))

def mom_autocov(zc, valid, Hmax=45):
    """Method-of-moments OU fit from the autocovariance of the demeaned log-gain.
    For an OU/AR(1) log-gain the SIGNAL autocovariance is gamma(h)=sigma_l^2 phi^h,
    so  log gamma(h) = log sigma_l^2 + h log phi  is linear in the lag h.  White
    Poisson observation noise inflates ONLY gamma(0), so a fit over lags h>=1
    analytically DECONVOLVES it (this is why the moment estimator resists the
    scene-residual and Poisson contamination that biases the full-path likelihood).
    The exponential-decay footprint scales with t_c, so the lag window is ADAPTIVE:
    H_eff ~ 2.5 t_c from a noise-free two-lag pilot phi=gamma(2)/gamma(1).
    Returns (t_c, sigma_l) or (nan,nan) if ill-posed."""
    n = len(zc)
    v = zc * valid
    vb = valid.astype(bool)
    g = np.full(Hmax + 1, np.nan)
    for h in range(Hmax + 1):
        m = vb[:n - h] & vb[h:]
        if m.sum() > 10:
            g[h] = float(np.mean(v[:n - h][m] * v[h:][m]))
    if not (np.isfinite(g[1]) and np.isfinite(g[2]) and g[1] > 0 and g[2] > 0):
        return float("nan"), float("nan")
    phi0 = min(max(g[2] / g[1], 1e-3), 0.9995)        # noise-free pilot (lags 1,2)
    tc0 = -1.0 / np.log(phi0)
    # SHORT adaptive window: long lags at high CV carry a nonlinear-log contamination
    # that biases t_c high (verified: identical bias with the true scene, so it is
    # intrinsic, not scene-coupling); a short window near t_c stays in the clean
    # short-lag decay regime.  Bounds keep >=4 lags for noise averaging.
    Heff = int(np.clip(round(1.0 * tc0), 4, 10))
    lags = np.arange(1, Heff + 1); gl = g[1:Heff + 1]
    ok = np.isfinite(gl) & (gl > 0)
    if ok.sum() < 3:                                   # fall back to short window
        lags = np.arange(1, 6); gl = g[1:6]; ok = np.isfinite(gl) & (gl > 0)
        if ok.sum() < 2:
            return float("nan"), float("nan")
    w = 1.0 / lags[ok]                                 # downweight noisy long lags
    slope, intercept = np.polyfit(lags[ok], np.log(gl[ok]), 1, w=w)
    phi = min(max(float(np.exp(slope)), 1e-4), 0.9999)
    sig2 = float(np.exp(intercept))                    # = sigma_l^2 (line value at h=0)
    tc = -1.0 / np.log(phi)
    return float(tc), float(np.sqrt(max(sig2, 1e-12)))

def _z_residuals(Y, x_hat, s_floor_frac=0.02):
    """Form the per-exposure log-gain residual series z_e (gauge removed), obs
    variance R_e and validity mask from record Y and scene estimate x_hat."""
    s = s_of_x(np.clip(x_hat, 0.0, None))
    s_floor = s_floor_frac * np.median(s[s > 0])
    valid = (s > s_floor).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.log((Y + 0.5) / np.maximum(PHI * s, 1e-12))
    R = 1.0 / (Y + 0.5)
    z = np.where(valid > 0, z, 0.0)
    mu_hat = _profile_mean(z, R, valid)          # scene-scale gauge = constant offset
    return z - mu_hat, R, valid, mu_hat

def medium_mle_kalman(Y, x_hat, tc_bounds=(0.3, 400.0), sig_bounds=(3e-3, 1.4)):
    """CROSS-CHECK estimator: maximise the exact Kalman marginal likelihood of the
    log-gain residuals over (t_c, sigma_l) (the gain path integrated out).  Slower
    and, under scene-reconstruction error, positively biased in t_c (it integrates
    low-frequency scene residuals into the correlation estimate) -- documented, not
    the headline."""
    zc, R, valid, mu_hat = _z_residuals(Y, x_hat)
    vb = valid.astype(bool)
    def nll(theta):
        tc = np.exp(theta[0]); sig = np.exp(theta[1]); phi = np.exp(-1.0 / tc)
        return -kalman_loglik(zc, R, vb, phi, sig)
    tc0, sig0 = mom_autocov(zc, valid)
    starts = []
    if np.isfinite(tc0):
        starts.append((np.clip(tc0, *tc_bounds), np.clip(sig0, *sig_bounds)))
    starts += [(2.0, 0.12), (16.0, 0.12), (64.0, 0.12)]
    lo = np.array([np.log(tc_bounds[0]), np.log(sig_bounds[0])])
    hi = np.array([np.log(tc_bounds[1]), np.log(sig_bounds[1])])
    best = None
    for (t0, g0) in starts:
        res = minimize(nll, np.array([np.log(t0), np.log(g0)]), method="L-BFGS-B",
                       bounds=list(zip(lo, hi)), options=dict(maxiter=120, ftol=1e-10))
        if (best is None) or (res.fun < best.fun):
            best = res
    return float(np.exp(best.x[0])), cv_of_sigma_l(float(np.exp(best.x[1])))

def medium_estimate(Y, x_hat, slot=None, return_path=True):
    """PRIMARY medium estimator (the deliverable): (t_c, sigma_l, CV) by the
    adaptive-window autocovariance method of moments on the per-exposure log-gain
    residuals, with the scene-scale gauge profiled out.  `slot[e]` = time slot of
    exposure e (schedule); residuals are re-ordered into TIME order before the OU
    fit (the medium lives in wall-clock time, not acquisition order).  Optionally
    returns the RTS-smoothed gain path in TIME order and per-exposure order."""
    zc_e, R_e, valid_e, mu_hat = _z_residuals(Y, x_hat)
    if slot is None:
        slot = np.arange(N_EXP)
    order = np.argsort(slot)                       # exposures sorted by time slot
    zc = zc_e[order]; R = R_e[order]; valid = valid_e[order]
    tc_hat, sig_hat = mom_autocov(zc, valid)
    if not np.isfinite(tc_hat):
        tc_hat, sig_hat = 1e3, 1e-3               # degenerate: no drift detected
    out = dict(tc_hat=float(tc_hat), sigma_l_hat=float(sig_hat),
               cv_hat=cv_of_sigma_l(sig_hat), mu_hat=mu_hat,
               n_valid=int(valid.sum()))
    if return_path:
        phi = np.exp(-1.0 / max(tc_hat, 1e-3))
        ms, Vs = kalman_rts(zc, R, valid.astype(bool), phi, max(sig_hat, 1e-6))
        a_time = np.exp(mu_hat + ms + 0.5 * Vs)   # keep gauge for fair corr
        out["a_time"] = a_time
        a_exp = np.empty(N_EXP); a_exp[order] = a_time
        out["a_hat"] = a_exp                       # per-exposure (for gain correction)
    return out

# ---- joint dual-ledger estimator: one record -> (x_hat, medium params) ---
def joint_dual_ledger(Y, slot=None, n_outer=2):
    """From ONE bucket record Y (2048 totals): reconstruct the scene AND estimate
    the medium, with NO knowledge of (t_c, CV) and NO pilots.  Alternates:
      image step  : gain-robust recon (A2) then gain-corrected recon (A4) + TV
      medium step : adaptive autocovariance OU fit on the log-gain residuals
    The medium estimate is a pure post-processing by-product of the SAME record."""
    Yp = Y[0::2]; Ym = Y[1::2]
    x_hat = np.clip(arm_A2(Yp, Ym, LAM_A2), 1e-9, None)   # gain-robust init
    med = medium_estimate(Y, x_hat, slot=slot, return_path=True)
    for _ in range(max(n_outer - 1, 0)):
        a_hat = med["a_hat"]                        # per-exposure gain (camera ledger improves)
        x_hat = np.clip(arm_A4_gaincorr(Yp, Ym, a_hat[0::2], a_hat[1::2], LAM_TV),
                        0.0, None)
        med = medium_estimate(Y, x_hat, slot=slot, return_path=True)
    return dict(x_hat=x_hat, med=med)

# ==========================================================================
# schedules (same 2048 exposures / same photons; only the time ORDER changes)
# ==========================================================================
def make_schedule(kind, seed=0):
    """slot[e] = time slot at which exposure e is acquired.  Exposure e even ->
    (row e//2, m+); e odd -> (row e//2, m-).  All schedules share the identical
    exposure/pattern multiset and photon budget -- only order and inference change.
      paired : natural ladder order (complementary pair adjacent -> differential
               cancels the common-mode gain; the O4-A-protecting schedule).
      split  : all m+ first, then all m- (a pair's members 1024 slots apart ->
               differential maximally gain-corrupted; large cross-confusion).
      random : uniform random permutation (intermediate)."""
    if kind == "paired":
        return np.arange(N_EXP)
    if kind == "split":
        slot = np.empty(N_EXP, dtype=int)
        slot[0::2] = np.arange(NPIX)               # m+ -> first half of time
        slot[1::2] = NPIX + np.arange(NPIX)        # m- -> second half
        return slot
    if kind == "random":
        return np.random.default_rng(9000 + seed).permutation(N_EXP)
    raise ValueError(kind)

# ---- simulate one record (frozen M0 gain), schedule-aware ----------------
def simulate_record(sc, seed, tc, sigma_l, schedule="paired", sched_seed=0):
    """One 2048-exposure complementary record with frozen within-exposure-constant
    OU gain sampled in TIME order, then mapped to exposures by the schedule.
    Returns (Y[exposure], a_time[time slot], slot[exposure])."""
    x = scene_x[sc]; s = s_of_x(x)
    grng = np.random.default_rng(1000 + seed)
    a_time = np.ones(N_EXP) if sigma_l <= 0 else ou_path(grng, sigma_l, tc, N_EXP)
    slot = make_schedule(schedule, sched_seed)
    a_exp = a_time[slot]                            # gain multiplying each exposure
    nrng = np.random.default_rng(5000 + seed)
    Y = nrng.poisson(a_exp * PHI * s).astype(np.float64)
    return Y, a_time, slot

def gain_path_corr(a_hat_time, a_time):
    """Pearson correlation of recovered vs true gain path (time order)."""
    if a_hat_time.std() < 1e-12 or a_time.std() < 1e-12:
        return 1.0
    return float(np.corrcoef(a_hat_time, a_time)[0, 1])

# ==========================================================================
# oracle AR(1)/OU Fisher precision floors (R33 sec.4.2, eqs 4.2-4.7)
#   phi = exp(-Delta/t_c), Delta=1 (exposure units), T = N Delta = N_valid
# ==========================================================================
def oracle_floor(tc, N=N_EXP):
    """Cramer-Rao sd floors for a DIRECTLY-OBSERVED high-SNR gain path (the best any
    method could do, incl. the oracle monitor).  Returns relative sd floors."""
    phi = np.exp(-1.0 / tc)
    var_phi = (1.0 - phi ** 2) / N                       # (4.3)
    # t_c = -1/log phi  ->  dt_c/dphi = 1/(phi (log phi)^2) ;  sd_tc/tc via (4.5)
    r = 1.0 / tc
    var_tc_over_tc2 = (np.exp(2 * r) - 1.0) / (N * r * r)  # (4.5)
    var_logCV = (1.0 + phi ** 2) / (2.0 * N * (1.0 - phi ** 2))  # (4.4)
    return dict(sd_tc_rel=float(np.sqrt(var_tc_over_tc2)),
                sd_cv_rel=float(np.sqrt(var_logCV)),          # sd(log CV) ~ sd(CV)/CV
                sd_phi=float(np.sqrt(var_phi)))

# ==========================================================================
# baseline arms (T2): pilot-interleaved and dedicated oracle monitor
# ==========================================================================
def oracle_monitor_estimate(a_time):
    """Dedicated noiseless monitor beam samples a_t at every exposure (UPPER BOUND).
    Fit the OU directly to the clean log-gain path (gauge = subtract sample mean)."""
    l = np.log(np.maximum(a_time, 1e-12))
    l = l - l.mean()
    tc, sig = mom_autocov(l, np.ones_like(l))
    if not np.isfinite(tc):
        tc, sig = 1e3, 1e-3
    return dict(tc_hat=float(tc), cv_hat=cv_of_sigma_l(sig), sigma_l_hat=float(sig))

def ou_fit_irregular(times, l):
    """OU fit from log-gain samples at IRREGULAR times (pilots).  Pools all sample
    pairs, bins by integer lag |t_i-t_j|, and fits gamma(h)=sigma_l^2 phi^h.  Uses
    BOTH the dense within-pair lag-1 and the sparse between-pilot lag, giving the
    pilot baseline its strongest shot.  A noise floor (estimated from the far-lag
    gammas, which are ~0) selects the informative lags so the many near-zero long
    lags do not flatten the decay.  Returns (t_c, sigma_l)."""
    t = np.asarray(times, float); v = np.asarray(l) - np.asarray(l).mean()
    n = len(t)
    dt = np.abs(t[:, None] - t[None, :]).astype(int)
    pr = v[:, None] * v[None, :]
    iu = np.triu_indices(n, k=1)
    lags = dt[iu]; prods = pr[iu]
    if len(lags) == 0:
        return float("nan"), float("nan")
    ul = np.unique(lags)
    g = np.array([prods[lags == h].mean() for h in ul])
    npair = np.array([np.sum(lags == h) for h in ul])
    # the pilot record has two natural lag scales: lag-1 (within a pilot pair) and
    # the between-row spacing dt.  One decorrelation step (lag-1 -> lag~dt) pins the
    # decay; longer lags are noisy (gamma near the floor) and bias t_c.  Cap at
    # ~1.6 dt to keep the clean two-anchor decay.
    dt_row = ul[ul >= 5].min() if np.any(ul >= 5) else ul.max()
    keep = (g > 0) & (npair >= 3) & (ul <= 1.6 * dt_row)
    if keep.sum() < 2:
        keep = (g > 0) & (npair >= 1)
    sel = np.where(keep)[0]
    if len(sel) < 2:
        return float("nan"), float("nan")
    w = npair[sel]                                   # weight by gamma-estimate precision
    slope, intercept = np.polyfit(ul[sel], np.log(g[sel]), 1, w=w)
    phi = min(max(float(np.exp(slope)), 1e-5), 0.99999)
    sig2 = float(np.exp(intercept))
    return float(-1.0 / np.log(phi)), float(np.sqrt(max(sig2, 1e-12)))

def pilot_interleaved(sc, seed, tc, sigma_l, pilot_frac=0.05):
    """PILOT-INTERLEAVED baseline (prior art: interpolated monitoring, Yang OE2018;
    the mandatory strongest baseline per R33 sec.5.3).  Replace pilot_frac of the
    Hadamard rows by all-ones PILOT pairs that directly read the gain (s=sum x).
    Pilots cost scene rows (5% of the Hadamard coefficients go missing).  Total
    exposures/photons identical to the pilot-free arm.  Gain path from the pilot
    reads (dense lag-1 within a pair + sparse lag between pilot rows), OU fit."""
    x = scene_x[sc]; sumx = float(x.sum())
    n_pair = NPIX
    n_pilot = max(int(round(pilot_frac * n_pair)), 2)
    pilot_rows = np.unique(np.linspace(1, n_pair - 1, n_pilot).round().astype(int))
    is_pilot = np.zeros(n_pair, dtype=bool); is_pilot[pilot_rows] = True
    grng = np.random.default_rng(1000 + seed)
    a_time = np.ones(N_EXP) if sigma_l <= 0 else ou_path(grng, sigma_l, tc, N_EXP)
    s_meas = s_of_x(x)
    s_meas[2 * pilot_rows] = sumx                    # m+ pilot all-ones
    s_meas[2 * pilot_rows + 1] = sumx                # m- pilot all-ones
    nrng = np.random.default_rng(5000 + seed)
    Y = nrng.poisson(a_time * PHI * s_meas).astype(np.float64)
    # gain reads at pilot exposures (both members of each pilot pair)
    pe = np.concatenate([2 * pilot_rows, 2 * pilot_rows + 1])
    pe.sort()
    a_reads = Y[pe] / (PHI * sumx)
    l_reads = np.log(np.maximum(a_reads, 1e-9))
    tc_p, sig_p = ou_fit_irregular(pe, l_reads)
    # scene: interpolate gain across all slots, gain-correct, invert Hadamard with
    # the pilot rows' coefficients set to zero (missing = honest 5% scene cost)
    a_interp = np.interp(np.arange(N_EXP), pe, a_reads)
    Yp = Y[0::2]; Ym = Y[1::2]
    c = Yp / np.maximum(a_interp[0::2], 1e-9) - Ym / np.maximum(a_interp[1::2], 1e-9)
    c[pilot_rows] = 0.0                              # missing Hadamard coefficients
    x_hat = np.clip(rof_denoise(hadamard_inverse(c) / PHI, LAM_TV), 0.0, None)
    med = dict(tc_hat=float(tc_p), cv_hat=cv_of_sigma_l(sig_p), sigma_l_hat=float(sig_p))
    return dict(med=med, x_hat=x_hat, n_pilot=int(len(pilot_rows)),
                a_interp=a_interp, a_time=a_time, is_pilot_row=is_pilot)
