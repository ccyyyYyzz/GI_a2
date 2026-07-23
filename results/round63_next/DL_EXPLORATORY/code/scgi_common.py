#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
EXPLORATORY / DEV-ONLY.  NOT preregistered.  No commits.  Local GPU only.

SCGI-style learned-correction machinery for the DL_EXPLORATORY head-to-head
(Arm L) against our model-based DLGI (Arm D) on the frozen OU-drift forward model.

FAITHFUL-IN-SPIRIT reimplementation of Peng & Chen, "Learning-based correction
with Gaussian constraints for ghost imaging through dynamic scattering media"
(Opt. Lett. 2023, DOI 10.1364/OL.499787), from the method description only -- NOT
code-identical.  Their SCGI: train a U-Net to map a drift-corrupted measurement
sequence R to the "static-equivalent" measurement B, on SIMULATED corruption
(exponential-decay scaling factors), loss = MSE + a Gaussian-prior term; then do
ordinary (differential) GI on the corrected measurements B.

OUR reimplementation (documented divergences, all faithful-in-spirit):
  * forward model = OUR frozen complementary-Hadamard M0 pipeline (part1_gain_ladder
    / dl_common), so the head-to-head shares one physical model.
  * the U-Net predicts a per-exposure log-gain CORRECTION delta_e (the smooth
    multiplicative drift), gauge-demeaned; corrected measurement B_e = Y_e*exp(-delta_e).
    Predicting the correction and dividing it out is equivalent to predicting B and
    is what exposes the "implicit correction factors" the head-to-head must quantify.
  * inputs are per-position-baseline-centred log-counts (the baseline = mean clean
    log-signal per exposure position over the training scenes; isolates the pure
    signal structure so the conv net's job is temporal drift extraction).
  * loss = MSE(delta_hat, delta_true) + a Gaussian random-walk prior on delta_hat in
    acquisition/time order (the "Gaussian constraint" that the drift is slow/smooth).
  * reconstruction from corrected measurements uses the IDENTICAL arm_A4_gaincorr+TV
    inversion that Arm D and the oracle use, so Arm L vs Arm D differ ONLY in where
    the per-exposure gain estimate comes from (learned correction vs OU/Kalman).
"""
import os, sys, math
import numpy as np
import torch
import torch.nn as nn

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = "D:/GI_another"
# dl_common re-exports the frozen forward model (which imports part1_gain_ladder)
sys.path.insert(0, os.path.join(REPO, "results", "round63_next", "DUAL_LEDGER_PROBE", "code"))
import dl_common as DL

NPIX, SIDE, N_EXP, PHI = DL.NPIX, DL.SIDE, DL.N_EXP, DL.PHI
Mp, Mm = DL.Mp, DL.Mm
ou_path = DL.ou_path
mom_autocov = DL.mom_autocov
cv_of_sigma_l = DL.cv_of_sigma_l
sigma_l_of_cv = DL.sigma_l_of_cv

# record reshaped as a 2D image for the conv net: 2048 = 32 x 64 (row-major = time order)
IMG_H, IMG_W = 32, 64
assert IMG_H * IMG_W == N_EXP

# ---------------------------------------------------------------- procedural training scenes
def _disk(x, cy, cx, r, val):
    yy, xx = np.ogrid[:SIDE, :SIDE]
    x[(yy - cy) ** 2 + (xx - cx) ** 2 <= r * r] = val

def _rect(x, y0, x0, h, w, val):
    x[y0:y0 + h, x0:x0 + w] = val

def _stroke(x, y0, x0, y1, x1, val, width=1):
    n = int(max(abs(y1 - y0), abs(x1 - x0)) + 1)
    ys = np.linspace(y0, y1, n).round().astype(int)
    xs = np.linspace(x0, x1, n).round().astype(int)
    for yy, xx in zip(ys, xs):
        y_lo, y_hi = max(0, yy - width), min(SIDE, yy + width + 1)
        x_lo, x_hi = max(0, xx - width), min(SIDE, xx + width + 1)
        x[y_lo:y_hi, x_lo:x_hi] = val

def make_scene(rng):
    """One 32x32 [0,1] 'simple' scene: random mix of shapes / strokes / sparse dots /
    smooth blobs on a random background, rescaled toward the DEV brightness band
    (member sum ~ [200,600], matching the 6 bridge scenes so counts stay in-regime)."""
    x = np.zeros((SIDE, SIDE), dtype=np.float64)
    bg = float(rng.uniform(0.0, 0.35))
    x[:] = bg
    kind = rng.integers(0, 4)
    if kind == 0:                                   # shapes
        for _ in range(rng.integers(1, 5)):
            v = float(rng.uniform(0.3, 1.0))
            if rng.random() < 0.5:
                _disk(x, rng.integers(4, 28), rng.integers(4, 28), rng.integers(2, 8), v)
            else:
                _rect(x, rng.integers(0, 24), rng.integers(0, 24),
                      rng.integers(3, 12), rng.integers(3, 12), v)
    elif kind == 1:                                 # strokes / digit-like
        for _ in range(rng.integers(2, 6)):
            _stroke(x, rng.integers(2, 30), rng.integers(2, 30),
                    rng.integers(2, 30), rng.integers(2, 30),
                    float(rng.uniform(0.4, 1.0)), width=int(rng.integers(0, 2)))
    elif kind == 2:                                 # sparse dots
        x[:] = float(rng.uniform(0.0, 0.15))
        n = int(rng.integers(6, 40))
        ys = rng.integers(0, SIDE, n); xs = rng.integers(0, SIDE, n)
        x[ys, xs] = rng.uniform(0.5, 1.0, n)
    else:                                           # smooth blobs
        yy, xx = np.mgrid[:SIDE, :SIDE]
        for _ in range(rng.integers(1, 4)):
            cy, cx = rng.uniform(4, 28), rng.uniform(4, 28)
            s = rng.uniform(3, 9); a = rng.uniform(0.4, 1.0)
            x += a * np.exp(-((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * s * s))
    x = np.clip(x, 0.0, 1.0)
    # rescale amplitude toward a random target sum in the DEV band
    tgt = float(rng.uniform(220.0, 560.0))
    s = x.sum()
    if s > 1e-6:
        x = np.clip(x * (tgt / s), 0.0, 1.0)
    return x.ravel()

# ---------------------------------------------------------------- gain families
def gain_exp(rng):
    """THEIR assumption class: exponential-decay scaling.  a_e = lam^e, lam~U[0.9995,1].
    Random sign/phase kept off (monotone decay) to stay literal to the paper's family."""
    lam = float(rng.uniform(0.9995, 1.0))
    e = np.arange(N_EXP, dtype=np.float64)
    a = lam ** e
    return a / a.mean() * 1.0   # arbitrary level; gauge is demeaned downstream anyway

def gain_ou(rng, tc, cv):
    sig = sigma_l_of_cv(cv)
    return ou_path(rng, sig, tc, N_EXP)

# grids for the two OU/exp training families (item 1)
OU_TC_TRAIN = [8.0, 16.0, 32.0, 64.0]
OU_CV_TRAIN = [0.05, 0.15, 0.40]

def sample_gain(rng, family):
    """family in {'exp','ou'}.  Returns a per-exposure multiplicative gain a_e (>0)."""
    if family == "exp":
        return gain_exp(rng)
    tc = float(rng.choice(OU_TC_TRAIN)); cv = float(rng.choice(OU_CV_TRAIN))
    return gain_ou(rng, tc, cv)

# ---------------------------------------------------------------- one training record
def s_of_x(x):
    s = np.empty(N_EXP)
    s[0::2] = Mp @ x
    s[1::2] = Mm @ x
    return s

def simulate_train_record(x, a_e, rng):
    """Poisson complementary-Hadamard record with per-exposure gain a_e (frozen M0)."""
    s = s_of_x(x)
    Y = rng.poisson(np.maximum(a_e * PHI * s, 0.0)).astype(np.float64)
    return Y, s

# ---------------------------------------------------------------- dataset builder
def build_dataset(n, family, seed, baseline=None):
    """Generate n training records for the given gain `family` ('exp','ou','mix').
    Returns dict with float32 tensors: X (n,1,32,64) net input, D (n,1,32,64) target
    demeaned log-gain, W (n,1,32,64) validity weights, plus the frozen per-position
    signal baseline (n-independent) used to centre inputs.
    Seeds are far-offset (disjoint from the forward-model and test seeds)."""
    rng = np.random.default_rng(1_000_000 + seed)
    Xs = np.empty((n, N_EXP), dtype=np.float32)
    Ds = np.empty((n, N_EXP), dtype=np.float32)
    Ws = np.empty((n, N_EXP), dtype=np.float32)
    logsig = np.empty((n, N_EXP), dtype=np.float64)   # clean log-signal per position
    raw = np.empty((n, N_EXP), dtype=np.float64)      # log(Y+0.5)
    la = np.empty((n, N_EXP), dtype=np.float64)       # demeaned true log-gain
    valid = np.empty((n, N_EXP), dtype=np.float64)
    fams = (["exp", "ou"] if family == "mix" else [family])
    for i in range(n):
        x = make_scene(rng)
        fam = fams[i % len(fams)] if family == "mix" else family
        a_e = sample_gain(rng, fam)
        Y, s = simulate_train_record(x, a_e, rng)
        s_floor = 0.02 * np.median(s[s > 0])
        v = (s > s_floor).astype(np.float64)
        raw[i] = np.log(Y + 0.5)
        logsig[i] = np.log(np.maximum(PHI * s, 1e-9) + 0.5)
        l = np.log(np.maximum(a_e, 1e-12))
        # gauge: demean the true log-gain over valid exposures (weighted uniform)
        lm = (l * v).sum() / max(v.sum(), 1.0)
        la[i] = l - lm
        valid[i] = v
    # frozen per-position signal baseline (mean clean log-signal over records) -- isolates
    # the deterministic signal structure so the net sees ~ log-gain + scene deviation + noise
    if baseline is None:
        baseline = logsig.mean(axis=0)
    for i in range(n):
        xin = raw[i] - baseline
        # remove the per-record global gauge (mean over valid exposures)
        m = (xin * valid[i]).sum() / max(valid[i].sum(), 1.0)
        xin = xin - m
        Xs[i] = xin.astype(np.float32)
        Ds[i] = la[i].astype(np.float32)
        Ws[i] = valid[i].astype(np.float32)
    def to2d(a):
        return torch.from_numpy(a.reshape(n, 1, IMG_H, IMG_W))
    return dict(X=to2d(Xs), D=to2d(Ds), W=to2d(Ws),
                baseline=baseline.astype(np.float64))

# ---------------------------------------------------------------- compact U-Net
class DoubleConv(nn.Module):
    def __init__(self, ci, co):
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv2d(ci, co, 3, padding=1), nn.BatchNorm2d(co), nn.ReLU(inplace=True),
            nn.Conv2d(co, co, 3, padding=1), nn.BatchNorm2d(co), nn.ReLU(inplace=True))
    def forward(self, x): return self.net(x)

class UNet(nn.Module):
    """Compact U-Net for (1,32,64) -> (1,32,64) log-gain correction map."""
    def __init__(self, c=24):
        super().__init__()
        self.e1 = DoubleConv(1, c)
        self.e2 = DoubleConv(c, 2 * c)
        self.pool = nn.MaxPool2d(2)
        self.b = DoubleConv(2 * c, 4 * c)
        self.u2 = nn.ConvTranspose2d(4 * c, 2 * c, 2, stride=2)
        self.d2 = DoubleConv(4 * c, 2 * c)
        self.u1 = nn.ConvTranspose2d(2 * c, c, 2, stride=2)
        self.d1 = DoubleConv(2 * c, c)
        self.out = nn.Conv2d(c, 1, 1)
    def forward(self, x):
        e1 = self.e1(x)
        e2 = self.e2(self.pool(e1))
        b = self.b(self.pool(e2))
        d2 = self.d2(torch.cat([self.u2(b), e2], 1))
        d1 = self.d1(torch.cat([self.u1(d2), e1], 1))
        return self.out(d1)

# ---------------------------------------------------------------- correction + reconstruction (Arm L)
def apply_correction(model, baseline, Y, device):
    """Feed-forward SCGI correction: record Y (2048,) -> per-exposure gain estimate
    g_hat (2048,), gauge-normalised (geometric mean 1).  ONE-shot (no iterative refine),
    faithful to a learned feed-forward corrector."""
    model.eval()
    raw = np.log(Y + 0.5)
    xin = raw - baseline
    # per-record gauge removal on valid exposures (approx via all positions; DC handled by net)
    m = xin.mean()
    xin = (xin - m).astype(np.float32).reshape(1, 1, IMG_H, IMG_W)
    with torch.no_grad():
        d = model(torch.from_numpy(xin).to(device)).cpu().numpy().reshape(N_EXP)
    d = d - d.mean()                      # gauge: geometric-mean-1 gain
    g_hat = np.exp(d)
    return g_hat

def reconstruct_from_gain(Y, g_hat):
    """Differential GI on corrected measurements = IDENTICAL arm_A4 inversion the other
    arms use, with the learned per-exposure gain.  Returns x_hat (1024,)."""
    Yp = Y[0::2]; Ym = Y[1::2]
    return np.clip(DL.arm_A4_gaincorr(Yp, Ym, g_hat[0::2], g_hat[1::2], DL.LAM_TV), 0.0, None)

def medium_from_gain(g_hat):
    """Extract (t_c, CV) from SCGI's implicit correction factors -- 'their trash':
    same MoM autocovariance OU fit used everywhere, on the demeaned log-correction."""
    l = np.log(np.maximum(g_hat, 1e-12))
    l = l - l.mean()
    tc, sig = mom_autocov(l, np.ones_like(l))
    if not np.isfinite(tc):
        tc, sig = 1e3, 1e-3
    return dict(tc_hat=float(tc), sigma_l_hat=float(sig), cv_hat=cv_of_sigma_l(sig))
