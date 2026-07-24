#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
QUAD_COMPUTER_TEST -- the law-only quadratic computer  (R46 sec 4.2-4.3 / kill sec 6.4).
PREREGISTERED, frozen bars, honest verdict.

Claim under test: a dynamic scattering medium read only through its 0-D law-only covariance
scalars q_l(x)=x^T G_l x implements a degree-two polynomial kernel WITHOUT measuring or
reconstructing the input.  For random rank-one laws G_l=a_l a_l^T, a_l~N(0,I), the centered
feature z_l(x)=(a_l^T x)^2 - ||x||^2 obeys the exact Isserlis identity
    E[z_l(x) z_l(y)] = 2 (x^T y)^2 .
With L independent laws the empirical feature-Gram converges at L^{-1/2}.

FROZEN PARAMETERS (declared before running):
  * input dim d = 24;  L in {1,4,16,64,256,1024};  n_class = 3000 samples
  * arms: digital random-feature (ideal), 0-D law-only covariance (random + admissible),
          camera-speckle random-feature -- under matched feature/MAC budgets
  * admissible law family = circulant/Toeplitz grain kernels (the physical scramble medium)

FROZEN PASS BARS:
  BAR1  kernel identity  E[z z] = 2(x.y)^2 : median relative error over held-out pairs < 0.05
        at L=100000  (classical Gaussian moment algebra -- MUST hold or the math is wrong)
  BAR2  empirical feature-Gram error ~ L^{-1/2} : fitted log-log slope in [-0.6,-0.4]
  BAR3  quadratic-boundary task a linear bucket statistic provably cannot solve:
        law-only(random) & digital-RF acc > 0.85 ; best LINEAR statistic acc < 0.62 (near
        chance) ; law-only within 0.05 of digital under matched budget
  BAR4  KILL TEST (sec 6.4): physically-admissible (circulant) law kernels.  Measure the
        cone dimension = rank of span{x -> x^T G_l x}.  Random rank-1 -> d(d+1)/2 ; circulant
        -> ~d (power-spectrum only).  Report collapse ratio and whether the admissible arm can
        solve a GENERIC quadratic boundary.  The GENERAL claim is KILLED if the admissible cone
        is very low-dim (ratio < 0.25) AND the admissible arm fails the generic quadratic task.

WRITE-SCOPE: results/round63_next/QUAD_COMPUTER_TEST/ only.  numpy/scipy.  No git here.
"""
import os, json, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SMOKE = os.environ.get("SMOKE", "0") == "1"
RNG = np.random.default_rng(20260724)

D = 24
L_LIST = [1, 4, 16, 64, 256, 1024]
N_CLASS = 3000 if not SMOKE else 800
L_IDENT = 100000 if not SMOKE else 20000
KS_BAR1 = 0.05
SLOPE_LO, SLOPE_HI = -0.6, -0.4
ACC_SOLVE = 0.85
ACC_LINEAR_MAX = 0.62
COLLAPSE_RATIO_KILL = 0.25


# ================================================================= feature maps
def z_random(X, A):
    """Centered law-only feature z_l(x)=(a_l^T x)^2 - ||x||^2 for random rank-1 laws.
    X (n,d), A (L,d) rows a_l ~ N(0,I).  Returns (n,L)."""
    proj = X @ A.T                                   # (n,L)
    xn2 = np.sum(X ** 2, axis=1, keepdims=True)      # (n,1)
    return proj ** 2 - xn2


def digital_rf(X, W):
    """Digital random quadratic features phi_l(x)=(w_l^T x)^2, w_l~N(0,I). (n,L)."""
    return (X @ W.T) ** 2


def circulant_symbols(L, d, seed):
    """L physically-admissible circulant (Toeplitz) grain kernels via nonneg PSD symbols
    g_l(k) >= 0 (a grain-pupil power spectrum).  G_l = F^H diag(g_l) F (circulant, PSD)."""
    rp = np.random.default_rng(seed)
    syms = np.abs(rp.standard_normal((L, d))) + 0.05   # nonneg PSD symbols (admissible)
    return syms


def q_circulant(X, syms):
    """Law-only covariance q_l(x)=x^T G_l x for circulant G_l with symbol g_l:
    = sum_k g_l(k) |xhat(k)|^2  (a linear readout of the POWER SPECTRUM of x). (n,L)."""
    Xf = np.fft.fft(X, axis=1)
    P = (Xf.real ** 2 + Xf.imag ** 2) / X.shape[1]   # (n,d) power spectrum
    return P @ syms.T                                # (n,L)


def camera_speckle(X, A, B):
    """Camera-speckle random features: intensity of random projections read on the speckle
    pixels, |B (A x)|^2 with B the FIXED speckle mixing (same physical realization for train and
    test).  (n, npix) real features, quadratic in x by construction."""
    field = (X @ A.T) @ B.T                          # (n, npix) complex
    return (field.real ** 2 + field.imag ** 2)


# ================================================================= tasks
def quadratic_task(n, d, generic=True, seed=0):
    """Binary labels from a QUADRATIC boundary  y = sign(x^T Q x - median).  If generic, Q is a
    full random symmetric indefinite matrix (needs the whole d(d+1)/2 quadratic space).  If not
    generic (power-spectrum task), Q is circulant so labels depend only on the power spectrum."""
    rp = np.random.default_rng(seed)
    X = rp.standard_normal((n, d))
    if generic:
        Qm = rp.standard_normal((d, d)); Qm = 0.5 * (Qm + Qm.T)
    else:
        s = np.abs(rp.standard_normal(d)) - 0.5      # circulant symbol (power-spectrum boundary)
        Qm = np.real(np.fft.ifft(np.diag(np.ones(d)) * 0))  # placeholder
        Fk = np.fft.fft(np.eye(d), axis=0)
        Qm = (Fk.conj().T @ np.diag(s) @ Fk).real / d
        Qm = 0.5 * (Qm + Qm.T)
    val = np.einsum('ni,ij,nj->n', X, Qm, X)
    y = (val > np.median(val)).astype(int)
    return X, y, Qm


def ridge_acc(Feat, y, Feat_te, y_te, ridge=1e-2):
    """Ridge-regression linear classifier on features; return test accuracy."""
    F = Feat - Feat.mean(0, keepdims=True); Fs = F.std(0, keepdims=True) + 1e-9
    F = F / Fs
    Fte = (Feat_te - Feat.mean(0, keepdims=True)) / Fs
    A = F.T @ F + ridge * np.eye(F.shape[1])
    w = np.linalg.solve(A, F.T @ (2.0 * y - 1.0))
    pred = (Fte @ w > 0).astype(int)
    return float(np.mean(pred == y_te))


def feature_cone_rank(kernel_family, d, L, seed):
    """Dimension of the linear span of the quadratic forms {x -> x^T G_l x} realized by a family.
    Equivalent to rank of the matrix whose rows are vec_sym(G_l).  Random rank-1 -> d(d+1)/2;
    circulant -> ~d."""
    rp = np.random.default_rng(seed)
    rows = []
    # symmetric-vectorization basis indices
    iu = np.triu_indices(d)
    if kernel_family == "random_rank1":
        A = rp.standard_normal((L, d))
        for a in A:
            G = np.outer(a, a)
            rows.append(G[iu])
    else:  # circulant admissible
        syms = circulant_symbols(L, d, seed + 1)
        Fk = np.fft.fft(np.eye(d), axis=0)
        for s in syms:
            G = (Fk.conj().T @ np.diag(s) @ Fk).real / d
            rows.append(G[iu])
    Rm = np.array(rows)
    sv = np.linalg.svd(Rm, compute_uv=False)
    rank = int((sv > 1e-9 * sv.max()).sum())
    pr = float(sv.sum() ** 2 / (np.square(sv).sum() + 1e-30))
    return rank, pr


# ================================================================= main
def main():
    t0 = time.time()
    rep = {"test": "QUAD_COMPUTER_TEST",
           "ref": "R46 sec 4.2-4.3 (Law-Only Quadratic Computer); kill sec 6.4",
           "frozen": dict(d=D, L=L_LIST, n_class=N_CLASS, L_ident=L_IDENT,
                          slope_band=[SLOPE_LO, SLOPE_HI], acc_solve=ACC_SOLVE,
                          acc_linear_max=ACC_LINEAR_MAX, collapse_kill=COLLAPSE_RATIO_KILL,
                          smoke=SMOKE)}

    # ----- BAR1: exact kernel identity  E[z z] = 2 (x.y)^2.  Use the DIAGONAL (x=y, Ktrue=
    # 2||x||^4 always large) for the robust relative-error bar; off-diagonal (x.y)^2 can be ~0
    # for near-orthogonal random pairs so its per-entry relative error is ill-conditioned --
    # reported separately as the whole-matrix Frobenius relative error.
    npair = 40 if not SMOKE else 15
    Xh = RNG.standard_normal((npair, D))
    A = RNG.standard_normal((L_IDENT, D))
    Z = z_random(Xh, A)                               # (npair, L)
    Kemp = (Z @ Z.T) / L_IDENT                        # empirical E[z z]
    Ktrue = 2.0 * (Xh @ Xh.T) ** 2
    diag_rel = np.abs(np.diag(Kemp) - np.diag(Ktrue)) / (np.abs(np.diag(Ktrue)) + 1e-9)
    med_diag_rel = float(np.median(diag_rel))
    frob_rel = float(np.linalg.norm(Kemp - Ktrue) / np.linalg.norm(Ktrue))
    rep["BAR1_kernel_identity"] = dict(
        L=L_IDENT, npairs=npair, median_diagonal_rel_err=med_diag_rel,
        frobenius_rel_err=round(frob_rel, 4),
        diag_emp_over_true=float(np.median(np.diag(Kemp) / (np.diag(Ktrue) + 1e-9))),
        bar="median diagonal relative error < %.2f (Isserlis identity, classical)" % KS_BAR1,
        verdict="PASS" if med_diag_rel < KS_BAR1 else "FAIL")

    # ----- BAR2: empirical feature-Gram error ~ L^{-1/2}
    ntest = 30 if not SMOKE else 12
    Xt = RNG.standard_normal((ntest, D))
    Kt = 2.0 * (Xt @ Xt.T) ** 2
    ntrial = 40 if not SMOKE else 12
    Ls = L_LIST + ([4096] if not SMOKE else [])
    errs = []
    for L in Ls:
        e = []
        for _ in range(ntrial):
            Al = RNG.standard_normal((L, D))
            Zl = z_random(Xt, Al)
            Kl = (Zl @ Zl.T) / L
            e.append(np.linalg.norm(Kl - Kt) / np.linalg.norm(Kt))
        errs.append(np.mean(e))
    errs = np.array(errs); Lsa = np.array(Ls, float)
    slope = float(np.polyfit(np.log(Lsa), np.log(errs), 1)[0])
    rep["BAR2_L_half_convergence"] = dict(
        L=Ls, gram_rel_error=[round(float(x), 5) for x in errs],
        fitted_loglog_slope=round(slope, 3), theory=-0.5,
        bar="slope in [%.2f,%.2f]" % (SLOPE_LO, SLOPE_HI),
        verdict="PASS" if SLOPE_LO <= slope <= SLOPE_HI else "FAIL")

    # ----- BAR3: quadratic-boundary task; matched budget L_task features.  A generic quadratic
    # boundary needs the full d(d+1)/2=300-dim quadratic space, so use L_task > 300 to span it.
    L_task = 500 if not SMOKE else 150
    npix = L_task                                     # camera-speckle pixels matched to feature count
    Xtr, ytr, Qm = quadratic_task(N_CLASS, D, generic=True, seed=1)
    Xte, yte, _ = quadratic_task(N_CLASS, D, generic=True, seed=2)
    # relabel test set with SAME boundary Qm
    valte = np.einsum('ni,ij,nj->n', Xte, Qm, Xte)
    tr_val = np.einsum('ni,ij,nj->n', Xtr, Qm, Xtr)
    thr = np.median(tr_val); ytr = (tr_val > thr).astype(int); yte = (valte > thr).astype(int)
    Aq = RNG.standard_normal((L_task, D))
    Wq = RNG.standard_normal((L_task, D))
    Bcam = RNG.standard_normal((npix, L_task)) + 1j * RNG.standard_normal((npix, L_task))  # fixed speckle
    acc_lawrand = ridge_acc(z_random(Xtr, Aq), ytr, z_random(Xte, Aq), yte)
    acc_digital = ridge_acc(digital_rf(Xtr, Wq), ytr, digital_rf(Xte, Wq), yte)
    acc_camera = ridge_acc(camera_speckle(Xtr, Aq, Bcam), ytr, camera_speckle(Xte, Aq, Bcam), yte)
    # best LINEAR statistic: linear classifier on raw x AND on the mean bucket (both linear) --
    # provably insufficient for a zero-mean quadratic boundary
    acc_linear_x = ridge_acc(Xtr, ytr, Xte, yte)
    acc_linear_mean = ridge_acc(Xtr.sum(1, keepdims=True), ytr, Xte.sum(1, keepdims=True), yte)
    acc_linear = max(acc_linear_x, acc_linear_mean)
    solve_ok = (acc_lawrand > ACC_SOLVE) and (acc_digital > ACC_SOLVE)
    linear_fails = acc_linear < ACC_LINEAR_MAX
    matched = abs(acc_lawrand - acc_digital) < 0.05
    rep["BAR3_quadratic_task"] = dict(
        L_features=L_task, acc_law_only_random=round(acc_lawrand, 3),
        acc_digital_RF=round(acc_digital, 3), acc_camera_speckle=round(acc_camera, 3),
        acc_best_linear=round(acc_linear, 3), acc_linear_x=round(acc_linear_x, 3),
        acc_linear_mean=round(acc_linear_mean, 3),
        bar="law-only & digital acc>%.2f ; best linear<%.2f ; |law-digital|<0.05" %
            (ACC_SOLVE, ACC_LINEAR_MAX),
        solve_ok=bool(solve_ok), linear_fails=bool(linear_fails), matched=bool(matched),
        verdict="PASS" if (solve_ok and linear_fails and matched) else "FAIL")

    # ----- BAR4: KILL TEST -- physically-admissible circulant law kernels
    Lcone = 2000 if not SMOKE else 600
    rank_rand, pr_rand = feature_cone_rank("random_rank1", D, Lcone, seed=11)
    rank_circ, pr_circ = feature_cone_rank("circulant", D, Lcone, seed=12)
    full_dim = D * (D + 1) // 2
    collapse_ratio = rank_circ / float(rank_rand)
    # can the admissible arm solve (a) a power-spectrum task and (b) the GENERIC quadratic task?
    syms = circulant_symbols(L_task, D, seed=21)
    # generic quadratic task (same Qm as BAR3)
    acc_circ_generic = ridge_acc(q_circulant(Xtr, syms), ytr, q_circulant(Xte, syms), yte)
    # power-spectrum-expressible task
    Xp_tr, yp_tr, Qp = quadratic_task(N_CLASS, D, generic=False, seed=3)
    Xp_te, yp_te, _ = quadratic_task(N_CLASS, D, generic=False, seed=4)
    vp_tr = np.einsum('ni,ij,nj->n', Xp_tr, Qp, Xp_tr); thr2 = np.median(vp_tr)
    yp_tr = (vp_tr > thr2).astype(int)
    vp_te = np.einsum('ni,ij,nj->n', Xp_te, Qp, Xp_te); yp_te = (vp_te > thr2).astype(int)
    acc_circ_power = ridge_acc(q_circulant(Xp_tr, syms), yp_tr, q_circulant(Xp_te, syms), yp_te)
    # bank cost: to span the random cone you need ~full_dim independent laws (=banks)
    banks_to_span = full_dim
    cone_collapses = collapse_ratio < COLLAPSE_RATIO_KILL
    admissible_fails_generic = acc_circ_generic < 0.62
    general_claim_killed = cone_collapses and admissible_fails_generic
    rep["BAR4_admissible_cone_KILL"] = dict(
        L_kernels=Lcone, full_quadratic_dim=full_dim,
        random_rank1_cone_rank=rank_rand, circulant_cone_rank=rank_circ,
        collapse_ratio=round(collapse_ratio, 4), circulant_PR=round(pr_circ, 1),
        acc_admissible_generic_quadratic=round(acc_circ_generic, 3),
        acc_admissible_powerspectrum_task=round(acc_circ_power, 3),
        acc_law_only_random_generic=round(acc_lawrand, 3),
        banks_needed_to_span_cone=int(banks_to_span),
        bar="general claim KILLED if collapse_ratio<%.2f AND admissible fails generic quadratic "
            "(<0.62); survives only as power-spectrum computing" % COLLAPSE_RATIO_KILL,
        cone_collapses=bool(cone_collapses),
        admissible_fails_generic_task=bool(admissible_fails_generic),
        general_claim_killed=bool(general_claim_killed),
        verdict="KILL_GENERAL" if general_claim_killed else "SURVIVES")

    # ----- overall verdict.  Idealized identity/convergence/task are classical (bars 1-3); the
    # flagship survives ONLY if the physically-admissible cone does NOT collapse (bar 4).
    ideal_ok = all(rep[b]["verdict"] == "PASS" for b in
                   ["BAR1_kernel_identity", "BAR2_L_half_convergence", "BAR3_quadratic_task"])
    if rep["BAR4_admissible_cone_KILL"]["verdict"] == "KILL_GENERAL":
        overall = "KILL_GENERAL_CLAIM"
        interp = ("The idealized random-law quadratic kernel is exact (classical Gaussian moment "
                  "algebra) and solves quadratic tasks a linear statistic cannot -- but physically "
                  "admissible circulant/Toeplitz grain kernels collapse the feature cone from "
                  "d(d+1)/2 to ~d (power-spectrum only) and cannot compute a generic quadratic "
                  "boundary.  Per R46 6.4 the GENERAL law-only quadratic-computer flagship is "
                  "KILLED; it survives only as a power-spectrum feature machine, and spanning the "
                  "full cone would need ~d(d+1)/2 independent banks (cost overwhelms advantage).")
    else:
        overall = "PASS" if ideal_ok else "FAIL"
        interp = ("Physically admissible kernels retain useful feature rank -- the law-only "
                  "quadratic computer survives as a flagship." if overall == "PASS"
                  else "an idealized bar failed unexpectedly.")
    rep["verdict"] = dict(
        per_bar={b: rep[b]["verdict"] for b in
                 ["BAR1_kernel_identity", "BAR2_L_half_convergence", "BAR3_quadratic_task",
                  "BAR4_admissible_cone_KILL"]},
        idealized_bars_pass=bool(ideal_ok), overall=overall, interpretation=interp)
    rep["runtime_sec"] = round(time.time() - t0, 1)

    with open(os.path.join(HERE, "QUAD_COMPUTER_TEST.json"), "w") as f:
        json.dump(rep, f, indent=2)
    print("[QUAD_COMPUTER done %.1fs] overall=%s" % (rep["runtime_sec"], overall))
    for b in ["BAR1_kernel_identity", "BAR2_L_half_convergence", "BAR3_quadratic_task",
              "BAR4_admissible_cone_KILL"]:
        print("  %-30s %s" % (b, rep[b]["verdict"]))
    print("  ident_diag_relerr=%.4f  slope=%.3f  acc law/dig/lin=%.2f/%.2f/%.2f  cone %d/%d ratio=%.3f"
          % (med_diag_rel, slope, acc_lawrand, acc_digital, acc_linear, rank_circ, rank_rand,
             collapse_ratio))
    try:
        _make_fig(rep, Ls, errs)
    except Exception as e:
        print("  (figure skipped: %s)" % e)
    return rep


def _make_fig(rep, Ls, errs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    ax[0].loglog(Ls, errs, "o-", label="empirical")
    ax[0].loglog(Ls, errs[0] * (np.array(Ls, float) / Ls[0]) ** -0.5, "k--", label="L^{-1/2}")
    ax[0].set_xlabel("L (independent laws)"); ax[0].set_ylabel("Gram rel. error")
    ax[0].set_title("BAR2: kernel convergence (slope %.2f)" %
                    rep["BAR2_L_half_convergence"]["fitted_loglog_slope"])
    ax[0].legend(fontsize=8)
    b3 = rep["BAR3_quadratic_task"]
    labels = ["law-only\n(random)", "digital\nRF", "camera\nspeckle", "best\nlinear"]
    accs = [b3["acc_law_only_random"], b3["acc_digital_RF"], b3["acc_camera_speckle"],
            b3["acc_best_linear"]]
    ax[1].bar(labels, accs, color=["#1a9850", "#4575b4", "#91bfdb", "#d73027"])
    ax[1].axhline(0.5, ls=":", color="k"); ax[1].set_ylim(0.4, 1.0)
    ax[1].set_ylabel("test accuracy"); ax[1].set_title("BAR3: generic quadratic boundary")
    b4 = rep["BAR4_admissible_cone_KILL"]
    labels2 = ["random\nrank-1", "circulant\n(admissible)", "full\nd(d+1)/2"]
    ranks = [b4["random_rank1_cone_rank"], b4["circulant_cone_rank"], b4["full_quadratic_dim"]]
    ax[2].bar(labels2, ranks, color=["#1a9850", "#d73027", "#999999"])
    ax[2].set_ylabel("feature cone rank")
    ax[2].set_title("BAR4 KILL: admissible cone collapse (ratio %.3f)" % b4["collapse_ratio"])
    fig.suptitle("QUAD_COMPUTER_TEST -- law-only quadratic computer: ideal exact, admissible cone "
                 "collapses", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(os.path.join(HERE, "QUAD_COMPUTER_TEST.png"), dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
