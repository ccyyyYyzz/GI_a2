#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
THE SEVEN ARMS of the sealed detection probe (R41 sec 4.5), corrected-shot engine.

  1. FIXED-COV       : one repeated incoherent sealed code bank; running M x M bucket covariance and
                       the profiled efficient beyond-band score. The production detector.
  2. FRESH-COV-OPT   : a NEW band-limited code bank each medium state, using the EXACT
                       code-conditioned covariance score for that bank (per-bank matched filter on the
                       bank's own residual outer product) -- NOT a coordinate-changing sample
                       covariance. Full-strength mandatory comparator (F4 immunity is mean-only).
  3. FIXED-MEAN      : mean-only detector on the repeated bank. Exactly blind beyond-band (the wall).
  4. FRESH-MEAN      : strongest equal-budget fresh-pattern mean detector -- exact zero-information
                       control beyond B_p (P_t U_beta = 0 by band-limiting).
  5. TRUE-LAW ORACLE : covariance score built with the GENERATING law; nondeployable ceiling.
  6. CROSSFIT-LAW    : deployable score with the law/tangent estimated from held-out baseline banks.
  7. AMPLITUDE / LAG : dedicated medium diagnostics (attribution): lag-0 amplitude matched filter and
                       lag-1 (cross-bank) matched filter for the medium-correlation-time direction.

Analytic side: each arm's per-effective-bank profiled Fisher (or mean deflection) -> d'(T_eff).
MC side: measured d' from the corrected Gaussian bank generator (FIXED arms) or the fresh-code
per-bank score accumulation (FRESH arms). Read-only on sealed banks; writes nowhere.
"""
import numpy as np
import torch

import sealed_common as sc

DEV = sc.DEV
DT = sc.DT
M = sc.M


# ------------------------------------------------------------------ analytic d' helpers
def dprime_analytic(cell, c, T_eff):
    """Energy/direction d' for a beyond-band delta with coeffs c: d'^2 = T_eff c^T JB c."""
    JB = cell["JB"].cpu().numpy()
    return float(np.sqrt(max(T_eff * float(c @ (JB @ c)), 0.0)))


def T_det_strong(cell, lam, eps):
    """Banks to reach d'>=5 for an eps-energy delta at aggregate curvature lam (lam_mean or lam_max)."""
    return sc.DP_STRONG2 / (lam * eps * eps * cell["xnorm2"])


# ================================================================== ARM 1: FIXED-COV
def fixed_cov_score(cell, Sc, c):
    """Profiled beyond-band statistic t = <Sc, W(c)> - <V0, W(c)> for a batch of sample covariances Sc."""
    W = cell["make_W"](c)
    V0W = float((cell["V0d"] * W).sum().item())
    return (torch.einsum("aij,ij->a", Sc, W)).cpu().numpy() - V0W


def fixed_cov_mc(cell, c, delta_px, T_eff, n_rec, rng0):
    """MC d'/AUC for FIXED-COV on the corrected Gaussian bank generator."""
    Sc0, _ = sc.gen_records(cell, n_rec, T_eff, "H0", None, rng=rng0)
    Sc1, _ = sc.gen_records(cell, n_rec, T_eff, "beyond", delta_px, rng=rng0 + 1)
    t0 = fixed_cov_score(cell, Sc0, c)
    t1 = fixed_cov_score(cell, Sc1, c)
    return _dp_auc(t0, t1)


# ================================================================== ARM 2: FRESH-COV-OPT
def fresh_cov_fisher(cell_geo, c_unit, n_code_draws=8, x_np=None):
    """Analytic per-effective-bank profiled beyond-band Fisher AVERAGED over fresh incoherent code
    banks: E_P[ c^T J_B(P) c ]. Each fresh bank uses its OWN code-conditioned covariance score, so the
    aggregate d'^2(T_eff) = T_eff * E_P[c^T J_B c]. Returns (lam_eff, per_draw)."""
    vals = []
    for s in range(n_code_draws):
        P = sc.signed_codes(seed=5000 + s)
        cell = sc.setup_cell(**cell_geo, x_np=x_np, code_P=P, want_mc=False)
        JB = cell["JB"].cpu().numpy()
        vals.append(float(c_unit @ (JB @ c_unit)))
    return float(np.mean(vals)), vals


_UC_CODE = torch.tensor(sc.band_modes(1, sc.KP), device=DEV, dtype=DT)   # fixed signed-code Fourier basis
_D_UC = _UC_CODE.shape[1]


def _fresh_codes(n_banks, gen):
    """n_banks FRESH signed band-limited code banks as a stacked tensor (n_banks, M, N), row-max
    normalized (physical complementary pairs). Reuses the FIXED code Fourier basis (no per-bank QR)."""
    Ccoef = torch.randn(n_banks, M, _D_UC, device=DEV, dtype=DT, generator=gen)
    P = Ccoef @ _UC_CODE.t()                          # (n_banks, M, N)
    P = P / P.abs().amax(dim=2, keepdim=True)
    return P


def fresh_cov_mc(cell_geo, c_unit, eps, T_eff, n_rec, rng0, x_np=None, chunk=48):
    """MC d' for FRESH-COV-OPT: each bank draws a FRESH code P_t and a fresh medium state; the per-bank
    code-conditioned matched score s_t = <r_t r_t^T - V_t, W_t> is accumulated over banks (NOT a
    cross-bank sample covariance in shared coordinates). Fully on GPU, chunked over banks to bound the
    (banks, M, N) code memory. This is the exact code-conditioned covariance score at full strength."""
    gen = torch.Generator(device=DEV); gen.manual_seed(int(rng0))
    x = torch.tensor(sc.XW_np if x_np is None else np.asarray(x_np, float), device=DEV, dtype=DT)
    Phi_b = torch.tensor(sc.BETA_np[cell_geo["claim"]], device=DEV, dtype=DT)
    xn = float(torch.linalg.norm(x).item())
    c = torch.tensor(c_unit, device=DEV, dtype=DT) * eps * xn
    delta_px = Phi_b @ c
    Um_np, kabs = sc.medium_modes(cell_geo["kwf"] * sc.KP)
    Um = torch.tensor(Um_np, device=DEV, dtype=DT)
    S = torch.tensor(sc.power_spectrum(kabs, cell_geo["shape"], sc.SIG_EFF[cell_geo["sf"]]),
                     device=DEV, dtype=DT)
    sdS = torch.sqrt(S)
    dm = Um.shape[1]
    var_n = (S[None, :] * (Um ** 2)).sum(dim=1)
    eyeM = torch.eye(M, device=DEV, dtype=DT)

    def one_stream(xh, seed):
        g = torch.Generator(device=DEV); g.manual_seed(int(seed))
        total = torch.zeros((), device=DEV, dtype=DT)
        done = 0
        while done < T_eff:
            nb = min(chunk, T_eff - done)
            P = _fresh_codes(nb, g)                                    # (nb, M, N)
            Pabs = P.abs()
            # independent medium states per bank (frozen ledger: independent between banks) -> vectorized
            Z = sdS[None, :] * torch.randn(nb, dm, device=DEV, dtype=DT, generator=g)
            w = torch.exp(Z @ Um.t() - 0.5 * var_n[None, :])          # (nb, N), E[w]=1
            wxh = w * xh[None, :]
            mu = torch.einsum("bmn,bn->bm", P, wxh)                   # signed effective bucket mean
            flux = torch.einsum("bmn,n->bm", Pabs, x)                 # |P| x  (declared-scene throughput)
            Rvar = flux * (flux.mean(dim=1, keepdim=True) / sc.PHOT)  # (nb, M)
            shot = torch.randn(nb, M, device=DEV, dtype=DT, generator=g) * torch.sqrt(Rvar)
            b = mu + shot
            r = b - torch.einsum("bmn,n->bm", P, x)                   # code-conditioned residual
            # per-bank model cov V_t and matched filter W_t (declared law, target direction c)
            A = P * x[None, None, :]                                  # (nb, M, N)
            # KA[b] = Um diag(S) Um^T A[b]^T  -> (nb, N, M), batched
            KA = torch.einsum("nk,bkm->bnm", Um, S[:, None] * torch.einsum("kn,bmn->bkm", Um.t(), A))
            C = torch.einsum("bmn,bnl->bml", A, KA)                   # (nb, M, M)
            V = C + torch.diag_embed(Rvar)
            Vinv = torch.linalg.inv(V + 1e-9 * eyeM[None])
            Bk = torch.einsum("bmn,bnl->bml", (P * delta_px[None, None, :]), KA)
            dV = Bk + Bk.transpose(1, 2)                              # beyond-band cov change (per bank)
            Wb = torch.einsum("bml,blk,bkr->bmr", Vinv, dV, Vinv)     # raw beyond-band matched filter
            # per-bank amplitude profiling: remove the medium-amplitude nuisance (dominant) direction,
            # W_eff = Wb - (I_ba/I_aa) W_amp, so the fresh score matches the PROFILED analytic J_B.
            Wamp = torch.einsum("bml,blk,bkr->bmr", Vinv, C, Vinv)
            I_ba = 0.5 * torch.einsum("bml,blm->b", torch.einsum("bml,blk->bmk", Vinv, dV),
                                      torch.einsum("bkl,blm->bkm", Vinv, C))
            I_aa = 0.5 * torch.einsum("bml,blm->b", torch.einsum("bml,blk->bmk", Vinv, C),
                                      torch.einsum("bkl,blm->bkm", Vinv, C))
            Weff = Wb - (I_ba / (I_aa + 1e-30))[:, None, None] * Wamp
            obs = torch.einsum("bm,bn->bmn", r, r) - V
            total = total + torch.einsum("bmn,bmn->", obs, Weff)
            done += nb
        return float(total.item())

    t0 = np.array([one_stream(x, int(1e6) + 2 * s) for s in range(n_rec)])
    t1 = np.array([one_stream(x + delta_px, int(3e6) + 2 * s) for s in range(n_rec)])
    return _dp_auc(t0, t1)


# ================================================================== ARM 3: FIXED-MEAN
def fixed_mean_mc(cell, delta_px, T_eff, n_rec, rng0):
    """Mean-channel Mahalanobis detector on the repeated bank. Beyond-band delta -> mean shift = P delta,
    which is EXACTLY zero (band-limited codes), so this is blind (d'~0) to beyond-band anomalies."""
    Vinv = cell["Vinvd"]
    m0 = cell["Pt"] @ cell["x"]
    Sc0, Mb0 = sc.gen_records(cell, n_rec, T_eff, "H0", None, rng=rng0)
    Sc1, Mb1 = sc.gen_records(cell, n_rec, T_eff, "beyond", delta_px, rng=rng0 + 1)

    def stat(Mb):
        dmv = Mb - m0[None, :]
        return (T_eff * torch.einsum("ai,ij,aj->a", dmv, Vinv, dmv)).cpu().numpy()
    return _dp_auc(stat(Mb0), stat(Mb1))


# ================================================================== ARM 4: FRESH-MEAN
def fresh_mean_wall(cell_geo, c, T_eff):
    """The exact zero-information control: fresh band-limited codes + mean route. P_t U_beta = 0 for
    all fresh band-limited P_t (Fourier orthogonality), so the mean deflection is identically 0."""
    Phi_b = sc.BETA_np[cell_geo["claim"]]
    leak = []
    for s in range(4):
        P = sc.signed_codes(seed=6000 + s)
        leak.append(float(np.linalg.norm(P @ (Phi_b @ c)) / (np.linalg.norm(P) + 1e-30)))
    return dict(mean_deflection_rel=float(np.mean(leak)), dprime_analytic=0.0, T_eff=T_eff)


# ================================================================== ARM 5: TRUE-LAW ORACLE
def oracle_mc(cell, c, delta_px, T_eff, n_rec, rng0):
    """Covariance score built with the GENERATING law (declared == true). Nondeployable ceiling."""
    return fixed_cov_mc(cell, c, delta_px, T_eff, n_rec, rng0)   # cell built with true law -> oracle


# ================================================================== ARM 6: CROSSFIT-LAW
def crossfit_cell(cell_geo, baseline_banks_np, x_np=None):
    """Deployable detector whose declared medium law is ESTIMATED from held-out baseline banks
    (nuisance/tangent cross-fit) rather than assumed. Returns a cell whose FILTER side uses the
    estimated spectrum while the DATA side keeps the true law -> honest deployable comparator.

    `baseline_banks_np`: (n_banks, M) H0 bucket residuals from a held-out stream; we estimate the
    medium power via the projected bucket covariance and fit a flat/slope spectrum by moment match.
    """
    # estimate realized sigma_f and slope from the baseline bucket covariance eigen-energy
    Sc = np.cov(baseline_banks_np.T)                             # (M,M) empirical bucket covariance
    # project onto the medium-generated subspace to estimate contrast; declare a matched flat law
    tr_sig = max(np.trace(Sc) - M * (Sc.diagonal().min()), 1e-9)  # crude signal/shot split
    # map estimated aggregate power back to an effective sigma_f level (nearest declared level)
    sf_est = min(sc.SIGMA_F_LEVELS, key=lambda s: abs(
        np.log(sc.SIG_EFF[s]) - 0.5 * np.log(max(tr_sig, 1e-9) / M)))
    declared = dict(sf=sf_est, shape="flat", rho=sc.RHO,
                    Um=sc.medium_modes(cell_geo["kwf"] * sc.KP))
    cell = sc.setup_cell(**cell_geo, x_np=x_np, declared=declared)
    cell["sf_est"] = sf_est
    return cell


# ================================================================== ARM 7: AMPLITUDE / LAG
def amplitude_score(cell, Sc):
    """Medium-amplitude statistic t_amp = <Sc, W_amp> - <V0, W_amp> (lights on sigma_f change)."""
    W = cell["W_amp"]
    V0W = float((cell["V0d"] * W).sum().item())
    return (torch.einsum("aij,ij->a", Sc, W)).cpu().numpy() - V0W


def lag_score(cell, Lag1):
    """Medium-correlation (tau) statistic on the lag-1 cross-bank covariance: t_lag = <Lag1, W_lag>.
    Under H0 E[Lag1] = rho*C; a tau change moves rho -> shows up here but is PROFILED OUT of the
    beyond-band (lag-0) score, giving orthogonal attribution."""
    W = cell["W_lag"]
    L1 = float((cell["rho"] * cell["Cd"] * W).sum().item())
    return (torch.einsum("aij,ij->a", Lag1, W)).cpu().numpy() - L1


# ------------------------------------------------------------------ shared stats
def _dp_auc(t0, t1):
    t0 = np.asarray(t0)
    t1 = np.asarray(t1)
    dp = float((t1.mean() - t0.mean()) / (t0.std() + 1e-12))
    auc = float(np.mean(t1[:, None] > t0[None, :]))
    return dict(dprime=round(dp, 3), auc=round(auc, 4))


def _sep(a, b):
    a = np.asarray(a)
    b = np.asarray(b)
    return float((b.mean() - a.mean()) / (0.5 * (a.std() + b.std()) + 1e-12))


ARMS = ["FIXED_COV", "FRESH_COV_OPT", "FIXED_MEAN", "FRESH_MEAN", "TRUE_LAW_ORACLE",
        "CROSSFIT_LAW", "AMPLITUDE", "LAG"]


if __name__ == "__main__":
    import time
    t0 = time.time()
    geo = sc.BEST_CELL
    cell = sc.setup_cell(**geo)
    xn = np.sqrt(cell["xnorm2"])
    eps = 0.02
    c_unit = sc.energy_spread_delta(cell["db"], 1.0, 1.0, seed=100)   # unit energy-spread direction
    c_unit = c_unit / np.linalg.norm(c_unit)
    c = c_unit * eps * xn
    delta_px = cell["Phi_b"].cpu().numpy() @ c
    T_eff = int(round(T_det_strong(cell, cell["lam_mean"], eps)))
    print(f"[selftest] best cell eps={eps*100:.0f}% T_eff={T_eff} banks", flush=True)

    dpa = dprime_analytic(cell, c, T_eff)
    fx = fixed_cov_mc(cell, c, delta_px, T_eff, 200, rng0=11)
    print(f"  FIXED-COV        d'_an={dpa:.2f}  d'_mc={fx['dprime']:.2f}  AUC={fx['auc']:.3f}", flush=True)

    lam_eff, _ = fresh_cov_fisher(geo, c_unit, n_code_draws=6)
    dpa_fresh = float(np.sqrt(T_eff * lam_eff * (eps * xn) ** 2))
    print(f"  FRESH-COV-OPT    d'_an={dpa_fresh:.2f}  (lam_eff={lam_eff:.4e} vs fixed lam_mean={cell['lam_mean']:.4e})", flush=True)

    fm = fixed_mean_mc(cell, delta_px, T_eff, 200, rng0=21)
    print(f"  FIXED-MEAN       d'_mc={fm['dprime']:.2f}  AUC={fm['auc']:.3f}  (expect ~0: the wall)", flush=True)

    wall = fresh_mean_wall(geo, c, T_eff)
    print(f"  FRESH-MEAN       mean-deflection rel={wall['mean_deflection_rel']:.2e}  (expect ~0)", flush=True)

    orc = oracle_mc(cell, c, delta_px, T_eff, 200, rng0=31)
    print(f"  TRUE-LAW ORACLE  d'_mc={orc['dprime']:.2f}  AUC={orc['auc']:.3f}", flush=True)

    # amplitude / lag diagnostics
    Sc0, _, Lag1_0 = sc.gen_records(cell, 200, 800, "H0", None, tau_scale=1.0, rng=41, want_lag=True)
    Sc_amp, _ = sc.gen_records(cell, 200, 800, "H0", None, sf_scale=1.2, rng=42, want_lag=False)
    _, _, Lag1_l = sc.gen_records(cell, 200, 800, "H0", None, tau_scale=2.0, rng=43, want_lag=True)
    ta0 = amplitude_score(cell, Sc0); taA = amplitude_score(cell, Sc_amp)
    tl0 = lag_score(cell, Lag1_0); tlL = lag_score(cell, Lag1_l)
    print(f"  AMPLITUDE  sep(H0,sf+20%)={_sep(ta0,taA):.2f}", flush=True)
    print(f"  LAG        sep(H0,tau*2) ={_sep(tl0,tlL):.2f}", flush=True)
    print(f"[selftest done {time.time()-t0:.0f}s]", flush=True)
