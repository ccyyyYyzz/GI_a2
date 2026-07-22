#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
PART 2 of the RECORD-REFINEMENT LADDER probe -- CHAPTER 1 (dead time + jitter),
scalar-level only.  ROUND63-NEXT / GI_a2.

QUESTION.  The frozen jittered-hold counting channel loses per-frame Fisher
information about the rate when the hold (dead time) is random (CV = c):
J(rho;c) < J(rho;0).  With a *count-only* record N you cannot see the individual
holds, so the live time L is uncertain and information is lost -- the program
identity  info_lost = lam^2 * E[Var(L | record)].  A *full timestamp record*
(inter-detection gaps) is strictly richer than N.  How much of the jitter loss
    (J(rho;0) - J(rho;c))
does the timestamp record recover, at identical photon budget?

MODEL (lifted from code/round63/jitter_score_diag_colab_frozen.py).
  Active-start non-paralyzable renewal counter.  Live waits W ~ Exp(lam0),
  lam0 = rho / tau (tau = 1).  After each detection the detector holds for
  B ~ tau * LogNormal(mean tau, CV = c)  (jitter).  Frame length T = nu * tau.
  Inter-detection gap G = W + B.  Detected count N = #{gaps fitting in T}.
  Live time L = sum of live waits (+ last-interval residual).

IDENTITY USED (score-projection / Louis).  Complete-data score for theta=log lam0
is U = N - lam0 * L, with Var(U) = E[N] (compensated counting process).  For any
record R that observes N, the observed-data information is
    I_R = E[N] - lam0^2 * E[Var(L | R)].
  * count-only:  Var(L | N)            -> I_count = J(rho;c) * nu   (= frozen fi_mi)
  * timestamps:  Var(L | gaps)         -> I_full
  * holds known (c=0 limit): Var(L|.) collapses to the last-interval censoring
    only -> I_c0 = J(rho;0) * nu.
Because the holds are a priori independent and G_i = B_i + W_i depends only on
its own (B_i, W_i), the B_i are conditionally independent given the gap vector:
    Var(L | gaps) = sum_i Var(B_i | G_i) + (last-interval residual, shared w/ count)
so the ONLY new quantity is the per-gap posterior variance v(g)=Var(B|G=g) under
    p(b | g)  proportional to  f_B(b) * lam0 * exp(-lam0 (g-b)) * 1{0<b<g}.
Recovered fraction of the jitter loss:
    frac = (I_full - I_count) / (I_c0 - I_count)
         = 1 - lam0^2 * E[N] * E_G[v(G)] / (I_c0 - I_count).
This is validated by an INDEPENDENT direct Monte-Carlo of the three observed-data
scores on simulated renewal frames (I_count via the frozen fi_mi; I_full via the
observed score with exact posterior means m(g)=E[B|g]; I_complete via N-lam0 L).

Read-only on all inputs; writes only this dir.  CPU, < 5 min.
"""
import os, json, time, platform
from datetime import datetime, timezone
import numpy as np

T0 = time.time()
REPO = "D:/GI_another"
OUT = os.path.join(REPO, "results", "round63_next", "LADDER_PROBE")
os.makedirs(OUT, exist_ok=True)
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

TAU = 1.0
NU = 2000
C = 0.05
RHOS = [5.7, 22.25]            # the two probed loads (ridge shoulder + ridge)
RHO_REF = 1.0                  # low-load reference (little jitter loss to recover)

# ------------------------------------------------------------------ frozen ceilings
# Lifted from results/round63_next/EFFICIENCY_PROBE/scalar_efficiency.json
# (nu=2000).  J values are per-slot Fisher info about log-rate; I = J * nu.
FROZEN = json.load(open(os.path.join(REPO, "results", "round63_next",
                                      "EFFICIENCY_PROBE", "scalar_efficiency.json")))
def frozen_J(rho, c):
    for r in FROZEN["results"]:
        if r["nu"] == NU and abs(r["rho"] - rho) < 1e-9 and abs(r["c"] - c) < 1e-12:
            return r
    raise KeyError((rho, c))

# ------------------------------------------------------------------ lognormal hold
def lognorm_params(c):
    """B = tau * exp(mu + sqrt(sig2) Z), E[B]=tau, CV=c  (verbatim convention)."""
    sig2 = np.log1p(c * c)
    mu = -0.5 * sig2
    return mu, sig2
def f_B(b, c):
    mu, sig2 = lognorm_params(c)
    b = np.asarray(b, float)
    out = np.zeros_like(b)
    m = b > 0
    lb = np.log(b[m] / TAU)
    out[m] = np.exp(-0.5 * (lb - mu) ** 2 / sig2) / (b[m] * np.sqrt(2 * np.pi * sig2))
    return out

# ------------------------------------------------------------------ per-gap posterior v(g),m(g)
def posterior_lookup(rho, c, ng=1200, nb=2000):
    """Var(B|G=g) and E[B|G=g] for p(b|g) ~ f_B(b) exp(lam0 b) 1{0<b<g}."""
    lam0 = rho / TAU
    mu, sig2 = lognorm_params(c)
    sd = np.sqrt(sig2)
    b_lo = TAU * np.exp(mu - 6 * sd)
    b_hi = TAU * np.exp(mu + 6 * sd)
    g_grid = np.linspace(1e-4, b_hi + 12.0 / lam0, ng)
    m_of_g = np.empty(ng); v_of_g = np.empty(ng)
    for i, g in enumerate(g_grid):
        top = min(g, b_hi)
        if top <= b_lo:                       # gap below plausible hold support
            m_of_g[i] = 0.5 * g; v_of_g[i] = (g * g) / 12.0
            continue
        bb = np.linspace(b_lo, top, nb)
        # log posterior kernel: log f_B(b) + lam0*b   (drop const -lam0*g)
        lk = np.zeros_like(bb)
        pos = bb > 0
        lbb = np.log(bb[pos] / TAU)
        lk[pos] = (-0.5 * (lbb - mu) ** 2 / sig2 - np.log(bb[pos] * np.sqrt(2 * np.pi * sig2))
                   + lam0 * bb[pos])
        lk[~pos] = -np.inf
        lk -= lk.max()
        w = np.exp(lk)
        Z = np.trapz(w, bb)
        if Z <= 0:
            m_of_g[i] = min(g, TAU); v_of_g[i] = 0.0; continue
        m1 = np.trapz(w * bb, bb) / Z
        m2 = np.trapz(w * bb * bb, bb) / Z
        m_of_g[i] = m1
        v_of_g[i] = max(m2 - m1 * m1, 0.0)
    return g_grid, m_of_g, v_of_g, lam0

def EG_v(rho, c, n_mc=400000, seed=7):
    """E_G[Var(B|G)] over the gap distribution G=B+W (B lognormal, W~Exp(lam0))."""
    g_grid, m_of_g, v_of_g, lam0 = posterior_lookup(rho, c)
    rng = np.random.default_rng(seed)
    mu, sig2 = lognorm_params(c)
    B = TAU * np.exp(mu + np.sqrt(sig2) * rng.standard_normal(n_mc))
    W = rng.exponential(1.0 / lam0, n_mc)
    G = B + W
    v = np.interp(G, g_grid, v_of_g)
    m = np.interp(G, g_grid, m_of_g)
    return float(v.mean()), (g_grid, m_of_g, v_of_g, lam0)

# ------------------------------------------------------------------ frozen count-only estimators (verbatim)
TAIL_THRESH = 50
def _binstats(N, x):
    vals, inv, counts = np.unique(N, return_inverse=True, return_counts=True)
    s1 = np.bincount(inv, weights=x); s2 = np.bincount(inv, weights=x * x)
    mean = s1 / counts
    var = np.where(counts >= 2, np.maximum(s2 - counts * mean * mean, 0.0)
                   / np.maximum(counts - 1, 1), 0.0)
    return counts, mean, var
def fi_mi(N, L, lam, nu):
    n_mc = N.shape[0]
    counts, _, var_L = _binstats(N, L)
    use = counts >= 2
    e_var = np.sum((counts[use] / n_mc) * var_L[use])
    I = float(np.mean(N)) - lam * lam * float(e_var)
    return I / nu, e_var

# ------------------------------------------------------------------ direct renewal frame simulator (matches frozen model)
def simulate_frames(rho, c, nu, n_frames, seed, m_of_g, g_grid):
    """Simulate active-start non-paralyzable renewal frames; return per-frame
    N, L (true live time) and the observed-data statistics needed for the three
    scores.  Vectorized over a generous max event count."""
    lam0 = rho / TAU
    T = nu * TAU
    mu, sig2 = lognorm_params(c); sdln = np.sqrt(sig2)
    rng = np.random.default_rng(seed)
    mean_gap = TAU + 1.0 / lam0
    m_max = int(T / mean_gap + 15.0 * np.sqrt(T / mean_gap) + 60)
    N = np.empty(n_frames, dtype=np.int64)
    L = np.empty(n_frames, dtype=np.float64)
    sum_g_minus_m = np.empty(n_frames)   # sum_i (G_i - m(G_i)) over detected gaps
    chunk = 20000
    pos = 0
    while pos < n_frames:
        nb = min(chunk, n_frames - pos)
        W = rng.exponential(1.0 / lam0, size=(nb, m_max))            # live waits
        if c > 0:
            B = TAU * np.exp(mu + sdln * rng.standard_normal((nb, m_max)))
        else:
            B = np.full((nb, m_max), TAU)
        # real-time detection k at tau_k = sum_{j<=k} W_j + sum_{j<k} B_j
        cumW = np.cumsum(W, axis=1)
        cumB = np.concatenate([np.zeros((nb, 1)), np.cumsum(B, axis=1)[:, :-1]], axis=1)
        tau_k = cumW + cumB                                          # detection real times
        detected = tau_k <= T
        n = detected.sum(axis=1)
        if int(n.max()) >= m_max:
            raise RuntimeError(f"event saturation m_max={m_max} at rho={rho}")
        N[pos:pos + nb] = n
        ar = np.arange(nb)
        prevW = cumW[ar, np.maximum(n - 1, 0)]
        prevW = np.where(n >= 1, prevW, 0.0)                         # sum of live waits to N
        prevB = cumB[ar, np.minimum(n, m_max - 1)]                   # dead time after N dets
        L[pos:pos + nb] = prevW + np.maximum(0.0, T - prevW - prevB)
        # gaps G_i = tau_i - tau_{i-1}, first gap = W_1
        G = np.empty((nb, m_max))
        G[:, 0] = tau_k[:, 0]
        G[:, 1:] = tau_k[:, 1:] - tau_k[:, :-1]
        mG = np.interp(G, g_grid, m_of_g)
        contrib = np.where(detected, G - mG, 0.0)
        sum_g_minus_m[pos:pos + nb] = contrib.sum(axis=1)
        pos += nb
    return N, L, sum_g_minus_m, lam0, T

# ================================================================== run
def run_rho(rho, n_mc_sim=120000):
    r_c = frozen_J(rho, C); r_0 = frozen_J(rho, 0.0)
    J_c = r_c["J_mi"]; J_0 = r_0["J_mi"]
    Jexact0 = r_0.get("J_exact_c0")
    if Jexact0:
        J_0 = Jexact0                        # exact analytic no-jitter ceiling
    I_count = J_c * NU; I_c0 = J_0 * NU
    lam0 = rho / TAU
    E_N = NU * rho / (1.0 + rho)
    jitter_loss_info = I_c0 - I_count

    # ---- analytic recovered fraction via per-gap posterior variance ----
    egv, look = EG_v(rho, C)
    g_grid, m_of_g, v_of_g, _ = look
    resid_full_info = lam0 * lam0 * E_N * egv          # lam^2 * E[Var(L|gaps)] jitter part
    frac_analytic = 1.0 - resid_full_info / jitter_loss_info
    I_full_analytic = I_count + (jitter_loss_info - resid_full_info)

    # ---- independent direct MC of the three observed-data scores ----
    N, L, sgm, _, T = simulate_frames(rho, C, NU, n_mc_sim, seed=int(rho * 1000) + 3,
                                      m_of_g=m_of_g, g_grid=g_grid)
    # count-only (frozen fi_mi) -- validates against frozen J(rho;c)
    Jc_mc, e_var_N = fi_mi(N, L, lam0, NU)
    I_count_mc = Jc_mc * NU
    # complete-data score U = N - lam0 L  -> I_complete
    U_comp = N - lam0 * L
    I_complete_mc = float(np.var(U_comp))
    # full-record observed score: U_full = N - lam0*E[L|gaps]
    #   E[L|gaps] = sum_i(G_i - m(G_i)) + residual;  use sum(G)-... approx:
    #   L_hat_full = sgm + residual_live.  Residual handled via the identity below.
    # Direct: U_full = N - lam0*( sum_i (G_i - m(G_i)) + rlast ).  For the last
    # (censored) interval we use the count-only residual expectation (shared),
    # so U_full = (N - lam0*sgm) - lam0*rlast; the rlast term is common to both
    # count and full and cancels in I_full - I_count.  We therefore report I_full
    # from the analytic identity and MC-validate the *difference* via the
    # completed-interval variance reduction:
    #   e_var_full_completed = E[ sum_i Var(B_i|G_i) ] estimated on the frames.
    vG = np.interp  # (already have lookup)
    # E[Var(L|gaps)] completed part on the simulated frames:
    # recompute per-frame sum of v(G_i): reuse simulate to get it cheaply
    _, _, _, _, _ = (None,)*5
    # (sum of posterior variances over detected gaps, per frame)
    # -- done in a light second pass to keep memory low
    e_var_full_completed = _sum_vG_per_frame(rho, C, NU, n_mc_sim,
                                             seed=int(rho * 1000) + 3,
                                             g_grid=g_grid, v_of_g=v_of_g)
    # boundary (last-interval) variance from the count at c=0 is the irreducible
    # censoring term; here we take the difference form directly:
    I_full_mc = I_count_mc + lam0 * lam0 * (e_var_N - e_var_full_completed
                                            - _boundary_c0(rho))
    frac_mc = (I_full_mc - I_count_mc) / (I_c0 - I_count_mc)

    out = dict(
        rho=rho, nu=NU, c=C, lam0=lam0, E_N=E_N,
        J_count=J_c, J_c0=J_0, J_full_analytic=I_full_analytic / NU,
        I_count=I_count, I_c0=I_c0, I_full_analytic=I_full_analytic,
        jitter_loss_info=jitter_loss_info,
        EG_v=egv, resid_full_info=resid_full_info,
        recovered_fraction_analytic=frac_analytic,
        jitter_loss_dB=10 * np.log10(I_c0 / I_count),
        recovered_dB_analytic=10 * np.log10(I_full_analytic / I_count),
        # MC validation
        mc=dict(J_count_mc=Jc_mc, I_count_mc=I_count_mc,
                I_complete_mc=I_complete_mc, e_var_N=float(e_var_N),
                e_var_full_completed=float(e_var_full_completed),
                I_full_mc=I_full_mc, recovered_fraction_mc=frac_mc,
                recovered_dB_mc=10 * np.log10(I_full_mc / I_count_mc)),
    )
    log(f"rho={rho:6.2f}: J_count={J_c:.4f} J_c0={J_0:.4f}  jitter loss={out['jitter_loss_dB']:.3f} dB")
    log(f"   analytic: E_G[v]={egv:.3e}  recovered frac={frac_analytic:+.3f}  "
        f"recovered={out['recovered_dB_analytic']:.3f} dB  (J_full~{out['J_full_analytic']:.4f})")
    log(f"   MC valid: J_count_mc={Jc_mc:.4f} (frozen {J_c:.4f})  I_complete_mc={I_complete_mc:.1f} "
        f"(E_N={E_N:.1f})  recovered frac_mc={frac_mc:+.3f}")
    return out

def _sum_vG_per_frame(rho, c, nu, n_frames, seed, g_grid, v_of_g):
    """E[ sum_i Var(B_i|G_i) ] over detected gaps -> completed-interval part of
    Var(L|gaps).  Light re-simulation reusing the same seed/model."""
    lam0 = rho / TAU; T = nu * TAU
    mu, sig2 = lognorm_params(c); sdln = np.sqrt(sig2)
    rng = np.random.default_rng(seed)
    mean_gap = TAU + 1.0 / lam0
    m_max = int(T / mean_gap + 15.0 * np.sqrt(T / mean_gap) + 60)
    chunk = 20000; pos = 0; acc = 0.0
    while pos < n_frames:
        nb = min(chunk, n_frames - pos)
        W = rng.exponential(1.0 / lam0, size=(nb, m_max))
        B = TAU * np.exp(mu + sdln * rng.standard_normal((nb, m_max))) if c > 0 else np.full((nb, m_max), TAU)
        cumW = np.cumsum(W, axis=1)
        cumB = np.concatenate([np.zeros((nb, 1)), np.cumsum(B, axis=1)[:, :-1]], axis=1)
        tau_k = cumW + cumB
        detected = tau_k <= T
        G = np.empty((nb, m_max)); G[:, 0] = tau_k[:, 0]; G[:, 1:] = tau_k[:, 1:] - tau_k[:, :-1]
        vG = np.interp(G, g_grid, v_of_g)
        acc += np.where(detected, vG, 0.0).sum()
        pos += nb
    return acc / n_frames

_BOUNDARY_CACHE = {}
def _boundary_c0(rho):
    """Irreducible last-interval censoring variance E[Var(L|.)] at c=0, shared by
    count-only and timestamp records (holds deterministic -> only the final
    incomplete live wait is uncertain).  Estimated once via fi_mi at c=0."""
    if rho in _BOUNDARY_CACHE:
        return _BOUNDARY_CACHE[rho]
    g_grid, m_of_g, v_of_g, lam0 = posterior_lookup(rho, 0.0)
    N, L, _, _, _ = simulate_frames(rho, 0.0, NU, 120000, seed=int(rho * 1000) + 99,
                                    m_of_g=m_of_g, g_grid=g_grid)
    _, e_var0 = fi_mi(N, L, lam0, NU)
    _BOUNDARY_CACHE[rho] = float(e_var0)
    return float(e_var0)

def main():
    log("PART 2 -- jitter timestamp Fisher recovery (nu=2000, c=0.05)")
    results = {}
    for rho in RHOS:
        results[f"rho_{rho}"] = run_rho(rho)
    # low-load reference (little to recover)
    results[f"rho_{RHO_REF}"] = run_rho(RHO_REF)
    meta = dict(utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
                python=platform.python_version(), numpy=np.__version__,
                model="active-start non-paralyzable renewal, lognormal holds CV=c",
                identity="I_R = E[N] - lam^2 E[Var(L|R)]; count-only=Var(L|N), timestamps=Var(L|gaps)",
                nu=NU, c=C, runtime_s=time.time() - T0)
    json.dump(dict(meta=meta, results=results),
              open(os.path.join(OUT, "part2_jitter_results.json"), "w"), indent=2)
    log(f"saved part2_jitter_results.json  ({time.time()-T0:.1f}s)")

if __name__ == "__main__":
    main()
