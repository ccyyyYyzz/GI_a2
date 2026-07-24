#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
COHERENCE_CROSSOVER -- the R45-NAMED GATE for the memory-effect line (label: the GATE, not
the campaign).  R45 sec B5 (Ore vein 5) + Rank-5 kill test.  PREREGISTERED, frozen bar.

Coherence-rank spectroscopy.  Code-pair features fill a latent R^2 quadratic coherence space:
    p = M(M+1)/2          (symmetric code-pair coordinates)
    q ~ R^2               (latent quadratic coherence products)
    r_eff(M;R) = p q/(p+q+1)                                             (5.2)
    crossover  M* ~ sqrt(2) R                                            (5.4)
    local slope alpha(M) = 2q/(p+q)                                      (5.5)
Consequences: p<<q -> r_eff~p~M^2/2 (unsaturated, slope 2); p>>q -> r_eff~q~R^2 (saturated).

Engine: WAVE_TWIN/wave_twin.py + twin_pool.py (shared bank pool at a fixed geometry).
r_eff(M) = participation ratio of the whitened covariance-Fisher Jacobian (KT4 method).
Independent R = participation ratio of the object-plane field mutual-coherence tensor (KT4).

TWO CONFIGURATIONS:
  (A) NAMED-GATE REFERENCE geometry (z2=5mm, l_c=50um, developed sigma_phi=2pi): R~110 per
      KT4, so M* ~ 155.  R^2 ~ 12000 exceeds the resolvable rank (bank/scene budget), so this
      arm can only show the pre-crossover slope~2; the crossover is NOT reachable here -- exactly
      KT4's finding.  Reported as the reference, not the decisive test.
  (B) REDUCED-R geometry (coarse screen / small illumination) auto-tuned to R~20-40 so the M
      sweep passes THROUGH M* ~ sqrt(2)R and r_eff saturates at ~R^2 -- the decisive crossover.

FROZEN PASS BAR (R45 B5 kill test / task): the crossover-fitted coherence rank sqrt(q_fit)
agrees with the independent KT4 R within 25%.  Corroborating (reported): pre-crossover slope
within 0.2 of 2; post-crossover slope < 0.8 where p/q > 3.  FAIL -> the memory-effect line
stays closed.

WRITE-SCOPE: results/round63_next/COHERENCE_CROSSOVER/ only.  torch/numpy.  No git here.
"""
import os, json, time
import numpy as np
import torch
import wave_twin as wt
import twin_pool as tp
from scipy.optimize import curve_fit

DEV = wt.DEV
SMOKE = os.environ.get("SMOKE", "0") == "1"
t0 = time.time()

# ------------------------------------------------------------------ FROZEN config
Z1 = 10e-3
RIDGE = 1e-4
SK = 256                                        # field sketch points for KT4 R
PASS_TOL = 0.25                                 # frozen: |sqrt(q_fit)-R_KT4|/R_KT4 < 0.25

if SMOKE:
    M_SWEEP_HI = [8, 16, 32, 48]
    M_SWEEP_LO = [4, 8, 12, 16, 24, 32]
    N_BANK_HI = 120
    N_BANK_LO = 200
    SCENE_KHI = 12
else:
    M_SWEEP_HI = [8, 16, 32, 64, 96, 128, 192]
    M_SWEEP_LO = [4, 8, 12, 16, 24, 32, 48, 64, 96]
    N_BANK_HI = 300
    N_BANK_LO = 800                             # > R^2 for reduced-R so saturation is resolvable
    SCENE_KHI = 20                              # scene basis band 0..20 -> ns ~ (41)^2 (> R^2)


# ------------------------------------------------------------------ helpers (KT4 conventions)
def matrix_isqrt(V):
    w, U = torch.linalg.eigh(0.5 * (V + V.T))
    w = torch.clamp(w, min=w.max() * 1e-12)
    return U @ torch.diag(w.rsqrt()) @ U.T


def eff_rank(sv):
    sv = np.asarray(sv, float); sv = sv[sv > 0]
    pr = float(sv.sum() ** 2 / (np.square(sv).sum() + 1e-300))
    cnt = int((sv > 1e-3 * sv.max()).sum())
    return pr, cnt


def cov_fisher_reff(G, f0, V0, U_scene, M, n_bank):
    """r_eff(M) = participation ratio of the whitened covariance-Fisher Jacobian over the scene
    basis (KT4 cov_fisher_rank).  dCov_s = Cov_bank(f0, g_s)+transpose, whitened by V0^{-1/2}."""
    ns = U_scene.shape[1]
    Rh = matrix_isqrt(V0)
    f0c = f0 - f0.mean(0, keepdim=True)
    gmat = (G.reshape(n_bank * M, -1) @ U_scene).reshape(n_bank, M, ns)
    gc = gmat - gmat.mean(0, keepdim=True)
    cross = torch.einsum('bi,bjs->ijs', f0c, gc) / (n_bank - 1)
    dC = cross + cross.transpose(0, 1)
    Jw = torch.einsum('ai,ijs,jb->abs', Rh, dC, Rh).reshape(M * M, ns)
    sv = torch.linalg.svdvals(Jw).cpu().numpy()
    return eff_rank(sv)


def field_coherence_R(Ep, z2, scr, offs, n_bank, pts):
    """Independent coherence rank R = participation ratio of the object-plane field mutual-
    coherence tensor eigenvalues (KT4 method), one code, spatial sketch across banks."""
    Esk = torch.empty(n_bank, SK, device=DEV, dtype=torch.complex64)
    e0 = Ep[0]
    for tk in range(n_bank):
        ph = wt.screen_crop(scr, int(offs[tk, 0]), int(offs[tk, 1]))
        u = wt.propagate(e0 * ph, z2) if z2 > 0 else e0 * ph
        Esk[tk] = u[pts[:, 0], pts[:, 1]].to(torch.complex64)
    Ec = Esk - Esk.mean(0, keepdim=True)
    C = (Ec.conj().T @ Ec) / (n_bank - 1)
    ev = torch.linalg.eigvalsh(0.5 * (C + C.conj().T)).abs().cpu().numpy()
    return eff_rank(ev)


def r_eff_model(M, q):
    p = M * (M + 1.0) / 2.0
    return p * q / (p + q + 1.0)


def fit_q(Ms, reffs):
    Ms = np.asarray(Ms, float); reffs = np.asarray(reffs, float)
    q0 = max(np.median(reffs), 4.0)
    try:
        popt, _ = curve_fit(r_eff_model, Ms, reffs, p0=[q0],
                            bounds=(1.0, 1e7), maxfev=20000)
        q = float(popt[0])
    except Exception:
        # grid fallback
        qs = np.logspace(0, 5, 400)
        err = [np.sum((reffs - r_eff_model(Ms, qq)) ** 2) for qq in qs]
        q = float(qs[int(np.argmin(err))])
    return q


def run_config(name, z2mm, l_c_um, sigma_phi, M_sweep, n_bank, illum_frac=1.0, seed=300):
    """Build shared pool at M_max (nested codes), sweep M, measure r_eff(M) + KT4 R, fit q."""
    z2 = z2mm * 1e-3
    Mmax = max(M_sweep)
    # nested codes: signed_codes(Mmax)[:M] == signed_codes(M) (row-major, fixed seed)
    Ep, Em, S = tp.code_fields_at_diffuser(Mmax, Z1, seed=11)
    # optional illumination aperture: keep only the central illum_frac of the DMD footprint
    if illum_frac < 1.0:
        Ep, Em = _apodize_illum(Ep, Em, illum_frac)
    # scene basis (band 0..SCENE_KHI) large enough that ns > R^2 for the reduced-R config
    U_scene = torch.tensor(wt.band_modes(0, SCENE_KHI), device=DEV, dtype=wt.RDT)
    ns = U_scene.shape[1]
    x64 = wt.witness_scene(4)
    x0t = torch.as_tensor(x64, device=DEV, dtype=wt.RDT)
    # build the full pool ONCE at Mmax
    Gp, Gm = tp.speckle_pool(Ep, Em, z2, sigma_phi, l_c_um, n_bank, seed=seed)
    # KT4 independent R (uses one code, its own screen realization -- same geometry)
    scr = wt.make_screen(l_c_um, sigma_phi, seed=seed)
    offs = wt.bank_offsets(n_bank, scr.shape[0], 300, seed=seed + 1)
    rng = np.random.default_rng(0); HALF = wt.DMD_PIX // 2
    pts = rng.integers(wt._cx - HALF, wt._cx + HALF, size=(SK, 2))
    R_PR, R_cnt = field_coherence_R(Ep, z2, scr, offs, n_bank, pts)

    cells = []
    for M in M_sweep:
        Gp_m = Gp[:, :M]; Gm_m = Gm[:, :M]
        G = (Gp_m - Gm_m).reshape(n_bank, M, -1)
        fp0 = torch.einsum('bmij,ij->bm', Gp_m, x0t)
        fm0 = torch.einsum('bmij,ij->bm', Gm_m, x0t)
        f0 = fp0 - fm0
        scp = float(wt.PHOT / (0.5 * (fp0 + fm0)).mean().clamp(min=1e-30))
        shot0 = ((fp0 + fm0).mean(0)) / scp
        C0 = torch.cov(f0.T)
        if C0.ndim == 0:
            C0 = C0.reshape(1, 1)
        V0 = C0 + torch.diag(shot0) + RIDGE * float(torch.diag(C0).mean()) * torch.eye(M, device=DEV, dtype=wt.RDT)
        (pr, cnt) = cov_fisher_reff(G, f0, V0, U_scene, M, n_bank)
        p = M * (M + 1) / 2
        cells.append(dict(M=M, p_codepairs=p, r_eff_PR=round(pr, 2), r_eff_count=cnt))
        print("   [%s M=%3d] p=%d  r_eff=%.1f (cnt %d)  [%.0fs]" % (name, M, p, pr, cnt, time.time() - t0))
        del G
        if DEV == "cuda":
            torch.cuda.empty_cache()

    Ms = [c["M"] for c in cells]
    reffs = [c["r_eff_PR"] for c in cells]
    q_fit = fit_q(Ms, reffs)
    R_fit = float(np.sqrt(q_fit))
    M_star = float(np.sqrt(2.0) * R_fit)
    # local slopes (finite difference of log r_eff vs log M)
    lM = np.log(Ms); lr = np.log(np.clip(reffs, 1e-9, None))
    slopes = {}
    for i in range(1, len(Ms)):
        slopes[str(Ms[i])] = round(float((lr[i] - lr[i - 1]) / (lM[i] - lM[i - 1])), 3)
    # pre-crossover slope: mean slope where p < q_fit; post-crossover where p > 3 q_fit
    pre = [(lr[i] - lr[i - 1]) / (lM[i] - lM[i - 1]) for i in range(1, len(Ms))
           if Ms[i] * (Ms[i] + 1) / 2 < q_fit]
    post = [(lr[i] - lr[i - 1]) / (lM[i] - lM[i - 1]) for i in range(1, len(Ms))
            if Ms[i - 1] * (Ms[i - 1] + 1) / 2 > 3 * q_fit]
    return dict(
        config=name, z2_mm=z2mm, l_c_um=l_c_um, sigma_phi=round(float(sigma_phi), 3),
        illum_frac=illum_frac, n_bank=n_bank, scene_dim=ns, M_sweep=Ms,
        cells=cells, R_KT4_PR=round(R_PR, 2), R_KT4_count=R_cnt,
        q_fit=round(q_fit, 1), R_fit_sqrt_q=round(R_fit, 2), M_star_pred=round(M_star, 1),
        local_slopes=slopes,
        pre_crossover_slope=round(float(np.mean(pre)), 3) if pre else None,
        post_crossover_slope=round(float(np.mean(post)), 3) if post else None,
        rel_err_Rfit_vs_KT4=round(abs(R_fit - R_PR) / (R_PR + 1e-9), 3))


def _apodize_illum(Ep, Em, frac):
    """Multiply the code fields by a soft central aperture (shrinks the illuminated area ->
    coarser object-plane speckle -> lower coherence rank R)."""
    n = wt.N
    yy, xx = torch.meshgrid(torch.arange(n, device=DEV), torch.arange(n, device=DEV), indexing='ij')
    r = torch.sqrt(((xx - wt._cx).float()) ** 2 + ((yy - wt._cx).float()) ** 2)
    R0 = frac * (wt.DMD_PIX / 2)
    ap = torch.exp(-(r / (R0 + 1e-9)) ** 2).to(Ep.dtype)
    return Ep * ap[None], Em * ap[None]


def calibrate_reduced_R(candidates, n_bank):
    """Measure KT4 R for a few reduced geometries; pick the one closest to R~30 (target 20-40)."""
    Ep, Em, _ = tp.code_fields_at_diffuser(4, Z1, seed=11)
    rng = np.random.default_rng(0); HALF = wt.DMD_PIX // 2
    pts = rng.integers(wt._cx - HALF, wt._cx + HALF, size=(SK, 2))
    out = []
    for (z2mm, l_c, sig, illum) in candidates:
        Epc, Emc = (_apodize_illum(Ep, Em, illum) if illum < 1.0 else (Ep, Em))
        scr = wt.make_screen(l_c, sig, seed=300)
        offs = wt.bank_offsets(n_bank, scr.shape[0], 300, seed=301)
        R_PR, _ = field_coherence_R(Epc, z2mm * 1e-3, scr, offs, n_bank, pts)
        out.append(dict(z2mm=z2mm, l_c=l_c, sigma_phi=round(float(sig), 3), illum=illum,
                        R=round(R_PR, 2)))
        print("   [calib] z2=%.1f l_c=%.0f sig=%.2f illum=%.2f -> R=%.1f [%.0fs]"
              % (z2mm, l_c, sig, illum, R_PR, time.time() - t0))
    # pick closest to 30
    best = min(out, key=lambda d: abs(d["R"] - 30.0))
    return best, out


def main():
    rep = {"test": "COHERENCE_CROSSOVER",
           "label": "R45-NAMED GATE for the memory-effect line (the GATE, not the campaign)",
           "ref": "R45 sec B5 (Ore vein 5) + Rank-5 kill test",
           "frozen": dict(pass_tol=PASS_TOL, smoke=SMOKE,
                          bar="crossover-fitted sqrt(q) within 25% of independent KT4 R")}
    print("[COHERENCE_CROSSOVER] device=%s" % DEV)

    # ---- (A) NAMED-GATE REFERENCE geometry: z2=5mm, l_c=50, developed
    print(" config A (named-gate reference, high-R): z2=5mm l_c=50 developed")
    A = run_config("A_ref_highR", 5.0, 50.0, 2 * np.pi, M_SWEEP_HI, N_BANK_HI, seed=300)
    rep["config_A_named_gate_reference"] = A

    # ---- (B) REDUCED-R geometry: auto-calibrate to R~20-40, then sweep through M*
    print(" config B: calibrating reduced-R geometry (target R~30)")
    cand = [(2.0, 250.0, 0.6, 1.0), (2.0, 400.0, 0.5, 0.5),
            (1.0, 300.0, 0.4, 0.6), (2.0, 200.0, 0.8, 0.7)]
    if SMOKE:
        cand = [(2.0, 250.0, 0.6, 1.0), (2.0, 400.0, 0.5, 0.5)]
    best, calib = calibrate_reduced_R(cand, N_BANK_LO)
    print(" config B chosen: z2=%.1f l_c=%.0f sig=%.2f illum=%.2f (R~%.1f)"
          % (best["z2mm"], best["l_c"], best["sigma_phi"], best["illum"], best["R"]))
    B = run_config("B_reduced_R", best["z2mm"], best["l_c"], best["sigma_phi"],
                   M_SWEEP_LO, N_BANK_LO, illum_frac=best["illum"], seed=300)
    rep["config_B_reduced_R"] = B
    rep["config_B_calibration"] = calib

    # ---- FROZEN PASS BAR on the decisive (reduced-R) config
    rel = B["rel_err_Rfit_vs_KT4"]
    pass_bar = rel < PASS_TOL
    pre_ok = (B["pre_crossover_slope"] is not None) and (abs(B["pre_crossover_slope"] - 2.0) < 0.2)
    post_ok = (B["post_crossover_slope"] is not None) and (B["post_crossover_slope"] < 0.8)
    rep["verdict"] = dict(
        decisive_config="B_reduced_R",
        R_KT4=B["R_KT4_PR"], R_fit_sqrt_q=B["R_fit_sqrt_q"], rel_err=rel,
        pre_crossover_slope=B["pre_crossover_slope"], post_crossover_slope=B["post_crossover_slope"],
        primary_bar="|sqrt(q_fit)-R_KT4|/R_KT4 < %.2f" % PASS_TOL,
        primary_pass=bool(pass_bar),
        corroborating=dict(pre_slope_near_2=bool(pre_ok), post_slope_below_0p8=bool(post_ok)),
        overall="PASS" if pass_bar else "FAIL",
        interpretation=(
            "Coherence-rank spectroscopy VALIDATED: the code-count crossover recovers the "
            "independent field coherence rank within 25%% -- the memory-effect line's named gate "
            "opens." if pass_bar else
            "The crossover-fitted rank disagrees with the independent KT4 R by >25%% -> per R45 "
            "the coherence-rank spectroscopy claim fails and the memory-effect line stays CLOSED "
            "(the algebraic rank<=min(N,R^2) bound is untouched)."))
    rep["runtime_sec"] = round(time.time() - t0, 1)

    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "COHERENCE_CROSSOVER.json"), "w") as f:
        json.dump(rep, f, indent=2)
    print("[COHERENCE_CROSSOVER done %.1fs] overall=%s" % (rep["runtime_sec"], rep["verdict"]["overall"]))
    print("  config B: R_KT4=%.1f  R_fit=sqrt(q)=%.1f  rel_err=%.3f (bar<%.2f)  pre/post slope=%s/%s"
          % (B["R_KT4_PR"], B["R_fit_sqrt_q"], rel, PASS_TOL,
             B["pre_crossover_slope"], B["post_crossover_slope"]))
    try:
        _make_fig(rep)
    except Exception as e:
        print("  (figure skipped: %s)" % e)
    return rep


def _make_fig(rep):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for key, col, lbl in [("config_A_named_gate_reference", "#762a83", "A: named-gate ref (high-R)"),
                          ("config_B_reduced_R", "#1b7837", "B: reduced-R (decisive)")]:
        c = rep[key]
        Ms = np.array([x["M"] for x in c["cells"]], float)
        rf = np.array([x["r_eff_PR"] for x in c["cells"]], float)
        ax[0].loglog(Ms, rf, "o-", color=col, label=lbl)
        mm = np.linspace(Ms.min(), Ms.max(), 100)
        ax[0].loglog(mm, r_eff_model(mm, c["q_fit"]), "--", color=col, alpha=0.6)
        ax[0].axhline(c["R_KT4_PR"] ** 2, ls=":", color=col, alpha=0.5)
        ax[0].axvline(np.sqrt(2) * c["R_KT4_PR"], ls=":", color=col, alpha=0.4)
    ax[0].loglog(Ms, Ms * (Ms + 1) / 2, "k:", alpha=0.4, label="p=M(M+1)/2 (slope 2)")
    ax[0].set_xlabel("M (codes)"); ax[0].set_ylabel("r_eff (cov-Fisher PR)")
    ax[0].set_title("Coherence-rank spectroscopy: r_eff(M)=pq/(p+q+1)")
    ax[0].legend(fontsize=7.5)
    # bar comparison of fitted vs KT4 R
    labels = ["A ref\nR_KT4", "A ref\nsqrt(q_fit)", "B\nR_KT4", "B\nsqrt(q_fit)"]
    A = rep["config_A_named_gate_reference"]; B = rep["config_B_reduced_R"]
    vals = [A["R_KT4_PR"], A["R_fit_sqrt_q"], B["R_KT4_PR"], B["R_fit_sqrt_q"]]
    ax[1].bar(labels, vals, color=["#762a83", "#c2a5cf", "#1b7837", "#a6dba0"])
    ax[1].set_ylabel("coherence rank R")
    ax[1].set_title("B: |sqrt(q_fit)-R_KT4|/R_KT4=%.2f (bar<0.25) -> %s"
                    % (rep["verdict"]["rel_err"], rep["verdict"]["overall"]))
    fig.suptitle("COHERENCE_CROSSOVER -- R45 named gate for the memory-effect line", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "COHERENCE_CROSSOVER.png"), dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
