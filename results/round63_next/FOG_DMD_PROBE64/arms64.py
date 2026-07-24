#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
THE SIX ARMS of the beyond-modulator-band sealed probe (R39 ask 5).

  1. fixed+MLE (PRODUCTION)  : fixed sealed bank + two-stage moment->Kalman-EM MLE (ft.kalman_em).
  2. fixed+moment            : fixed sealed bank + E5a moment matching (ft.cov_estimate).
  3. fresh+mean (THE WALL)   : fresh band-limited patterns each epoch + mean route. Band-limited
                               first-moment operators CANNOT leave their band -> beyond-band = 1.0.
                               The channel-impossibility control (R39 ask 2).
  4. fresh+cov BEST-EFFORT   : fresh patterns each epoch + covariance route with GLS weighting.
                               The design-half comparator (concentration principle: does the
                               FIXED bank beat FRESH at equal budget?).  R39 mandates best-effort
                               (GLS) form before the design half is frozen.
  5. oracle ceiling          : medium known -> ft.solve_oracle (nondeployable ceiling).
  6. classic averaging       : fixed bank, plain time-averaged bucket inversion (ignores the
                               fluctuation entirely) -- the ordinary ghost-imaging baseline
                               (band-limited -> beyond-band wall).

ALL arms share ONE medium realization W per record (fair paired comparison, like the DLGI shared
gain path).  EQUAL BUDGET: every arm uses M patterns x T epochs exposures at the same photons/
bucket; only the pattern strategy (fixed vs fresh) and readout (mean vs covariance vs MLE) differ.

Seed hygiene mirrors dlgi_common: disjoint partition salts, independent SeedSequence children for
(medium, poisson) so streams never alias and never cross partitions.  Read-only on sealed banks.
CPU/GPU (heavy solves via ft on GPU when available).  [R39] marks referee-frozen choices.
"""
import hashlib

import numpy as np
import torch
from scipy.fftpack import dct, idct

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "FOG_DMD_PROBE"))
import fog_tracker as ft  # noqa: E402

DEV = "cuda" if torch.cuda.is_available() else "cpu"

# disjoint partition salts (no seed crosses partitions) -- distinct from DLGI 700/710/720
_PART_SALT = {"calibration": 810, "confirmatory": 800, "mismatch": 820, "oracle": 830}


def _u32(s):
    return int(hashlib.sha256(str(s).encode()).hexdigest()[:8], 16)


def record_rngs(partition, scene_id, replicate):
    """Independent (medium_rng, poisson_rng) for one record; pure function of
    (partition, scene, replicate); streams never alias, never cross partitions."""
    ss = np.random.SeedSequence([_PART_SALT[partition], _u32(scene_id), int(replicate)])
    m, p = ss.spawn(2)
    return np.random.default_rng(m), np.random.default_rng(p)


class ProbeConfig:
    """Frozen geometry for one probe scale (n=64 production; n=16 for the selftest)."""
    def __init__(self, n=64, PB=8, med_lo=1, med_hi=16, sig_f=0.30, tau=8.0, M=128):
        self.n = n; self.N = n * n; self.PB = PB
        self.med_lo = med_lo; self.med_hi = med_hi; self.sig_f = sig_f
        self.tau = tau; self.rho = float(np.exp(-1.0 / tau)); self.M = M
        band = [(i, j) for i in range(n) for j in range(n) if med_lo <= i + j <= med_hi]
        self.Ub = np.linalg.qr(np.array([self._spike(i, j) for (i, j) in band]).T)[0]
        self.db = self.Ub.shape[1]
        self.coeff_sd = sig_f * np.sqrt(self.N / self.db)
        self.superres, self.pat_box, self.cov = self._masks(band)

    def _spike(self, i, j):
        g = np.zeros((self.n, self.n)); g[i, j] = 1.0
        return idct(idct(g, axis=0, norm="ortho"), axis=1, norm="ortho").ravel()

    def D2(self, a):
        return dct(dct(a.reshape(self.n, self.n), axis=0, norm="ortho"), axis=1, norm="ortho")

    def I2(self, A):
        return idct(idct(A, axis=0, norm="ortho"), axis=1, norm="ortho").ravel()

    def _masks(self, band):
        n, PB = self.n, self.PB
        pat = np.zeros((n, n), bool); pat[:PB + 1, :PB + 1] = True
        cov = np.zeros((n, n), bool)
        for pi in range(PB + 1):
            for pj in range(PB + 1):
                for (ki, kj) in band:
                    for si in (pi + ki, abs(pi - ki)):
                        for sj in (pj + kj, abs(pj - kj)):
                            if si < n and sj < n:
                                cov[si, sj] = True
        return cov & (~pat), pat, cov

    def bl_patterns(self, M_, seed):
        rp = np.random.default_rng(seed); ps = []
        for _ in range(M_):
            Cp = np.zeros((self.n, self.n))
            Cp[:self.PB + 1, :self.PB + 1] = rp.standard_normal((self.PB + 1, self.PB + 1))
            Cp[0, 0] = 0
            f = self.I2(Cp); f /= np.abs(f).max(); ps.append(0.5 + 0.45 * f)
        return np.array(ps)

    def medium_field(self, T, medium_rng):
        return self.medium_field_custom(T, medium_rng)

    def medium_field_custom(self, T, medium_rng, coeff_sd=None, rho=None, gen_sd=None):
        """Lognormal OU medium W (T x N), E[w]=1 per pixel. Overrides let the F5 mismatch axes
        generate a TRUE medium with a different tau (rho), sigma_f (coeff_sd), or a non-flat
        per-mode radial profile (gen_sd), while the estimator still declares the frozen law."""
        rho = self.rho if rho is None else rho
        sd = self.coeff_sd if coeff_sd is None else coeff_sd
        sdv = sd if gen_sd is None else np.asarray(gen_sd, float)
        z = (sdv if gen_sd is not None else sd) * medium_rng.standard_normal(self.db)
        var_n = ((sdv ** 2) if gen_sd is not None else (sd ** 2)) * (self.Ub ** 2)
        var_n = var_n.sum(axis=1) if gen_sd is not None else (sd ** 2) * (self.Ub ** 2).sum(axis=1)
        W = np.zeros((T, self.N))
        for t in range(T):
            W[t] = np.exp(self.Ub @ z - 0.5 * var_n)
            z = rho * z + np.sqrt(1 - rho ** 2) * (sdv if gen_sd is not None else sd) * \
                medium_rng.standard_normal(self.db)
        return W

    def rotate_declared(self, eps, shell_hi=None, n_rot=None):
        """Return an Ub rotated by angle arcsin(eps) toward the just-beyond fine shell
        (medium band (med_hi, shell_hi]) -- the declared-band rotation for the F5 fine_rotation axis."""
        shell_hi = (self.med_hi + max(2, self.med_hi // 4)) if shell_hi is None else shell_hi
        shell = [(i, j) for i in range(self.n) for j in range(self.n)
                 if self.med_hi < i + j <= shell_hi]
        Up = np.linalg.qr(np.array([self._spike(i, j) for (i, j) in shell]).T)[0]
        k = min(self.db, Up.shape[1] if n_rot is None else n_rot)
        Urot = self.Ub.copy().astype(float)
        Urot[:, :k] = np.sqrt(1 - eps ** 2) * self.Ub[:, :k] + eps * Up[:, :k]
        return np.linalg.qr(Urot)[0]

    def beyond_band_nmse(self, xh, x, mask=None):
        m = self.superres if mask is None else mask
        s = np.dot(xh, x) / max(np.dot(xh, xh), 1e-12)
        E = self.D2(s * xh - x) ** 2; X0 = self.D2(x) ** 2
        d = X0[m].sum()
        return float(E[m].sum() / d) if d > 1e-12 else float("nan")

    def shell_mask(self, dk):
        """Delta-k shell: DCT freqs at Chebyshev distance dk beyond the pattern box, in coverage."""
        m = self.PB + dk
        out = np.zeros((self.n, self.n), bool)
        for i in range(self.n):
            for j in range(self.n):
                if max(i, j) == m and self.cov[i, j]:
                    out[i, j] = True
        return out

    def shell_errs(self, xh, x, max_dk=6):
        """Per-shell beyond-band NMSE, Dk=1..max_dk (the band-extension-depth profile)."""
        return {dk: self.beyond_band_nmse(xh, x, self.shell_mask(dk)) for dk in range(1, max_dk + 1)}

    def viable_depth(self, shell_errs, thr=0.30):
        """Largest contiguous Dk (from 1) with NMSE <= thr -- the claimed band-extension depth."""
        d = 0
        for dk in sorted(shell_errs):
            e = shell_errs[dk]
            if np.isfinite(e) and e <= thr:
                d = dk
            else:
                break
        return d

    def depth_mask(self, depth):
        """Union of shells Dk=1..depth (the PRODUCTION clip: content beyond depth is masked out)."""
        m = np.zeros((self.n, self.n), bool)
        for dk in range(1, depth + 1):
            m |= self.shell_mask(dk)
        return m


def _poisson_buckets(mu, poisson_rng, phot):
    scp = phot / mu.mean()
    B = poisson_rng.poisson(np.clip(mu, 0, None) * scp) / scp
    R_diag = np.clip(mu, 1e-9, None) / scp
    return B, R_diag


# ----------------------------------------------------------------------- arms
def arm_oracle(cfg, x, W, poisson_rng, Pfix, phot, lam_x=1e-6):
    mu = np.einsum("mn,tn->tm", Pfix, W * x)
    B, _ = _poisson_buckets(mu, poisson_rng, phot)
    return ft.solve_oracle(Pfix, W, B, np.ones_like(B), lam_x, DEV)


def arm_fixed_moment(cfg, x, W, poisson_rng, Pfix, phot, n_starts=3, steps=4000):
    mu = np.einsum("mn,tn->tm", Pfix, W * x)
    B, _ = _poisson_buckets(mu, poisson_rng, phot)
    rc = ft.cov_estimate(Pfix, cfg.Ub, B, cfg.coeff_sd, cfg.rho,
                         n_starts=n_starts, steps=steps, dev=DEV)
    return rc["per_start_x"][int(np.argmin(rc["objs"]))], rc


def arm_fixed_mle(cfg, x, W, poisson_rng, Pfix, phot, n_starts=3, steps=4000, n_em=25):
    mu = np.einsum("mn,tn->tm", Pfix, W * x)
    B, R_diag = _poisson_buckets(mu, poisson_rng, phot)
    rc = ft.cov_estimate(Pfix, cfg.Ub, B, cfg.coeff_sd, cfg.rho,
                         n_starts=n_starts, steps=steps, dev=DEV)
    x_init = rc["per_start_x"][int(np.argmin(rc["objs"]))]
    xL = ft.kalman_em(Pfix, cfg.Ub, B, R_diag, cfg.rho, cfg.coeff_sd, 1e-4,
                      n_em=n_em, dev=DEV, x_init=x_init)
    return xL


def arm_fresh_mean(cfg, x, W, poisson_rng, phot, pat_seed_base):
    """THE WALL: fresh band-limited patterns + mean route (band-limited -> beyond-band blind)."""
    T, N = W.shape
    AtA = np.zeros((N, N)); Atb = np.zeros(N)
    scp = phot / np.einsum("mn,tn->tm", cfg.bl_patterns(cfg.M, pat_seed_base), W * x).mean()
    for t in range(T):
        Pf = cfg.bl_patterns(cfg.M, pat_seed_base + 1 + t)
        mu = Pf @ (W[t] * x)
        B = poisson_rng.poisson(np.clip(mu, 0, None) * scp) / scp
        AtA += Pf.T @ Pf; Atb += Pf.T @ B
    lam = 1e-3 * np.trace(AtA) / N
    return np.clip(np.linalg.solve(AtA + lam * np.eye(N), Atb), 0, None)


def arm_classic_avg(cfg, x, W, poisson_rng, Pfix, phot):
    """Ordinary ghost-imaging averaging: time-average the buckets, single mean-route inversion
    with the fixed bank (ignores the fluctuation channel -> band-limited -> beyond-band wall)."""
    mu = np.einsum("mn,tn->tm", Pfix, W * x)
    B, _ = _poisson_buckets(mu, poisson_rng, phot)
    b_bar = B.mean(0)
    x0, _, _, _ = np.linalg.lstsq(Pfix, b_bar, rcond=None)
    return np.clip(x0, 0, None)


def arm_fresh_cov_gls(cfg, x, W, poisson_rng, phot, pat_seed_base, iters=2500,
                      batch=64, lr=2e-2, gls=True):
    """BEST-EFFORT fresh+covariance with GLS weighting (R39-mandated design comparator).
    Fresh patterns P_t each epoch; model within-epoch centered outer products
      E[r_i r_j] = coeff_sd^2 (P_t,i (.) x)^T UU^T (P_t,j (.) x),  r=B-P_t x
    fit x by Adam on the off-diagonal residual, weighted (GLS) by 1/(model_var+eps)."""
    T, N = W.shape; S = cfg.M
    Kg = (cfg.coeff_sd ** 2) * (cfg.Ub @ cfg.Ub.T)
    scp = phot / np.einsum("mn,tn->tm", cfg.bl_patterns(S, pat_seed_base), W * x).mean()
    Plist, Blist = [], []
    for t in range(T):
        Pf = cfg.bl_patterns(S, pat_seed_base + 1 + t)
        mu = Pf @ (W[t] * x)
        Blist.append(poisson_rng.poisson(np.clip(mu, 0, None) * scp) / scp)
        Plist.append(Pf)
    # mean-route init (in-band)
    AtA = np.zeros((N, N)); Atb = np.zeros(N)
    for t in range(T):
        AtA += Plist[t].T @ Plist[t]; Atb += Plist[t].T @ Blist[t]
    x0 = np.clip(np.linalg.solve(AtA + 1e-3 * np.trace(AtA) / N * np.eye(N), Atb), 1e-3, None)
    tK = torch.tensor(Kg, device=DEV, dtype=torch.float64)
    tP = torch.stack([torch.tensor(p, device=DEV, dtype=torch.float64) for p in Plist])   # T,S,N
    tB = torch.stack([torch.tensor(b, device=DEV, dtype=torch.float64) for b in Blist])   # T,S
    xp = torch.tensor(x0, device=DEV, dtype=torch.float64, requires_grad=True)
    opt = torch.optim.Adam([xp], lr=lr)
    mask = ~torch.eye(S, dtype=torch.bool, device=DEV)
    for it in range(iters):
        opt.zero_grad(); xr = torch.relu(xp)
        sub = torch.tensor(np.random.default_rng(it).choice(T, min(batch, T), replace=False),
                           device=DEV)
        Ps = tP[sub]; Bs = tB[sub]
        resid = Bs - torch.einsum("tsn,n->ts", Ps, xr)
        A = Ps * xr[None, None, :]
        G_pred = torch.einsum("tsn,nm,trm->tsr", A, tK, A)
        obs = resid[:, :, None] * resid[:, None, :]
        diff = (G_pred - obs)
        if gls:                                   # GLS: downweight high-model-variance entries
            w = 1.0 / (G_pred.detach() ** 2 + 1e-6)
            loss = ((diff ** 2) * w)[:, mask].sum()
        else:
            loss = (diff[:, mask] ** 2).sum()
        loss.backward(); opt.step()
        if it == int(0.72 * iters):
            for g_ in opt.param_groups:
                g_["lr"] = lr * 0.25
    return np.clip(xp.detach().cpu().numpy(), 0, None)


ARMS = ["fixed_mle", "fixed_moment", "fresh_mean", "fresh_cov_gls", "oracle", "classic_avg"]


def run_all_arms(cfg, x, partition, scene_id, replicate, T, phot, steps=4000):
    """Run all six arms on ONE shared medium realization; return {arm: beyond_band_nmse}."""
    med_rng, poi_rng = record_rngs(partition, scene_id, replicate)
    W = cfg.medium_field(T, med_rng)
    Pfix = cfg.bl_patterns(cfg.M, 10)          # the sealed fixed bank (frozen seed [R39])
    pat_base = 7_000_000 + _u32(scene_id) % 1_000_000
    out = {}
    xo = arm_oracle(cfg, x, W, np.random.default_rng(poi_rng.integers(1 << 31)), Pfix, phot)
    out["oracle"] = cfg.beyond_band_nmse(xo, x)
    xm, _ = arm_fixed_moment(cfg, x, W, np.random.default_rng(poi_rng.integers(1 << 31)),
                             Pfix, phot, steps=steps)
    out["fixed_moment"] = cfg.beyond_band_nmse(xm, x)
    xl = arm_fixed_mle(cfg, x, W, np.random.default_rng(poi_rng.integers(1 << 31)),
                       Pfix, phot, steps=steps)
    out["fixed_mle"] = cfg.beyond_band_nmse(xl, x)
    xw = arm_fresh_mean(cfg, x, W, np.random.default_rng(poi_rng.integers(1 << 31)), phot, pat_base)
    out["fresh_mean"] = cfg.beyond_band_nmse(xw, x)
    xc = arm_classic_avg(cfg, x, W, np.random.default_rng(poi_rng.integers(1 << 31)), Pfix, phot)
    out["classic_avg"] = cfg.beyond_band_nmse(xc, x)
    xg = arm_fresh_cov_gls(cfg, x, W, np.random.default_rng(poi_rng.integers(1 << 31)),
                           phot, pat_base)
    out["fresh_cov_gls"] = cfg.beyond_band_nmse(xg, x)
    return out


if __name__ == "__main__":
    # tiny-scale selftest (n=16): confirm all six arms EXECUTE and the qualitative ordering holds
    # (fresh_mean & classic_avg ~ 1.0 wall; oracle near 0; fixed arms recover beyond-band).
    import time
    t0 = time.time()
    cfg = ProbeConfig(n=16, PB=3, med_lo=1, med_hi=6, M=24)
    print(f"selftest cfg n={cfg.n} N={cfg.N} PB={cfg.PB} db={cfg.db} "
          f"superres-frac {cfg.superres.mean():.3f}", flush=True)
    rs = np.random.default_rng(4)
    C = np.zeros((cfg.n, cfg.n)); C[:cfg.PB + 1, :cfg.PB + 1] = rs.standard_normal((cfg.PB + 1, cfg.PB + 1))
    for (i, j) in [(5, 2), (2, 5), (4, 4), (6, 1), (1, 6), (5, 5), (7, 0), (0, 7)]:
        C[i, j] = 2.0 * rs.choice([-1, 1])
    x = cfg.I2(C); x = x - x.min(); x = x / x.max()
    res = run_all_arms(cfg, x, "confirmatory", "confirmatory_witness_0", 0, T=1024, phot=1e5, steps=2500)
    print("beyond-band NMSE by arm on FULL superres zone (n=16, T=1024):", flush=True)
    for a in ARMS:
        print(f"  {a:16s} {res[a]:.3f}", flush=True)
    # The aggregate full-zone metric is dominated by deep corrupting shells (G3 shell-depth finding).
    # The PRODUCTION metric is shell-resolved: check the near shell (Dk=1) where recovery is claimed.
    med_rng, poi_rng = record_rngs("confirmatory", "confirmatory_witness_0", 0)
    W = cfg.medium_field(1024, med_rng); Pfix = cfg.bl_patterns(cfg.M, 10)
    xm, _ = arm_fixed_moment(cfg, x, W, np.random.default_rng(poi_rng.integers(1 << 31)),
                             Pfix, 1e5, steps=2500)
    xw2 = arm_fresh_mean(cfg, x, W, np.random.default_rng(poi_rng.integers(1 << 31)), 1e5,
                         7_000_000 + _u32("confirmatory_witness_0") % 1_000_000)
    sh1 = cfg.shell_mask(1)
    fixed_dk1 = cfg.beyond_band_nmse(xm, x, sh1); wall_dk1 = cfg.beyond_band_nmse(xw2, x, sh1)
    print(f"  [Dk=1 shell] fixed_moment {fixed_dk1:.3f}  vs  fresh_mean(wall) {wall_dk1:.3f}", flush=True)
    ok = (res["fresh_mean"] > 0.9 and res["classic_avg"] > 0.9 and res["oracle"] < 0.5
          and fixed_dk1 < wall_dk1)
    print(f"\nARMS SELFTEST: {'PASS' if ok else 'CHECK'} "
          f"(wall arms>0.9, oracle<0.5, fixed beats wall on the Dk=1 shell)  [{time.time()-t0:.0f}s]", flush=True)
