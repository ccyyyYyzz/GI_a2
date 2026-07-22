"""
MOLT REVIVAL DESK CALCULATION  (round63_next / MOLT_REVIVAL_CALC)
================================================================================
ONE question, frozen decision rule.  Does MOLT's dead flagship (bar 6: noisy
equal-time joint recon only TIED linear on natural dense scenes, because
box-constrained fiber perturbations carry sub-order-one Delta_p2 vs an achievable
~10% p2 precision) REVIVE on sparse / binary / high-contrast scene classes, where
the revival hypothesis is: (H1) fiber ambiguities carry ORDER-ONE Delta_p2, and
(H2) refresh + larger-C lift the precision by 10-30x -> a >=10x d' margin.

FROZEN DECISION RULE (verbatim from the task):
    REVIVE  iff  median margin >= 10x on at least one HONEST scene class,
    where margin = (aggregated Fisher-weighted d' at achievable precision) / 3.
    margin 2-5x = STAY DEAD.

This module imports the COMMITTED probe physics verbatim (does NOT rewrite it):
    saturation_core : Q_exp, power_sums, design_matrix, three_level_design,
                      annealed_fim, default_a_levels, draw_W, ...
    gi_operator     : P=1024, build_operator, load_scene
from results/round63_next/SATURATION_JAILBREAK_PROBE/code/.

Physics recap (saturation_core docstring):
  masked scene v = m (x) x in [0,1]^P;  survival Q(a;v)=prod_p 1/(1+a v_p);
  L(a)=-log Q = a p1 - a^2 p2/2 + ...;  MOLT's 2nd operator = p2 = M(x (x) x).
  n_eff = p1^2/p2 = beta2^{-1}.  Load t=a p1; ridge t*=2.8214.  B_LIN=1e4.

Everything below is CPU, deterministic (fixed seeds), <~a few minutes.
"""
import json
import os
import sys
import time
import numpy as np

# ---- import the COMMITTED probe physics (verbatim, no rewrite) -------------- #
PROBE_CODE = os.path.join(
    "results", "round63_next", "SATURATION_JAILBREAK_PROBE", "code")
sys.path.insert(0, os.path.abspath(PROBE_CODE))
import saturation_core as sc          # noqa: E402
import gi_operator as op              # noqa: E402

# ---- constants -------------------------------------------------------------- #
SIDE = op.SIDE                        # 32
P = op.P                              # 1024
B_LIN = 1.0e4                         # original linear detected-photon budget/mask
C_BASE = 3600                         # frozen probe microcell count
T_STAR = 2.821439372122079
K_D = 51                             # dense linear bank (task spec)
K_S = 51                             # sparse MOLT bank (task spec)
SUPPORT = 32                         # sparse mask row weight
OUT = os.path.abspath(os.path.join(
    "results", "round63_next", "MOLT_REVIVAL_CALC"))
rng_global = np.random.default_rng(63_0000)

report = {"meta": dict(K_D=K_D, K_S=K_S, C_base=C_BASE, B_LIN=B_LIN,
                       t_star=T_STAR, support=SUPPORT, P=P,
                       decision_rule="REVIVE iff median margin>=10x on an honest "
                       "class; margin=aggregated d'/3; 2-5x=STAY DEAD")}


# =========================================================================== #
# 0.  BANKS  (task spec: M_D = 51 dense 50% masks; M_S = 51 32-support masks)
# =========================================================================== #
M_D = op.build_operator(K=K_D, P=P, density=0.5, seed=651_000)      # 51 x 1024


def build_sparse_bank(K=K_S, support=SUPPORT, seed=651_032):
    r = np.random.default_rng(seed)
    M = np.zeros((K, P))
    for i in range(K):
        M[i, r.choice(P, size=support, replace=False)] = 1.0
    return M


M_S = build_sparse_bank()
L_BANK = np.vstack([M_D, M_S])                                       # 102 x 1024
L_PINV = np.linalg.pinv(L_BANK)


# =========================================================================== #
# 1.  SCENE CLASSES  (32x32, values in [0,1])
# =========================================================================== #
def scene_binary_sparse(s, seed):
    """(a) s random binary supports: mix of single spots and short bars, {0,1}."""
    r = np.random.default_rng(seed)
    img = np.zeros((SIDE, SIDE))
    placed = 0
    while placed < s:
        if r.random() < 0.55:                              # single spot
            img[r.integers(SIDE), r.integers(SIDE)] = 1.0
            placed = int(img.sum())
        else:                                              # short bar (len 2-4)
            ln = int(r.integers(2, 5)); horiz = r.random() < 0.5
            rr, cc = r.integers(SIDE), r.integers(SIDE)
            for k in range(ln):
                a, b = (rr, min(cc + k, SIDE - 1)) if horiz else (min(rr + k, SIDE - 1), cc)
                img[a, b] = 1.0
            placed = int(img.sum())
    # trim to exactly s active pixels
    on = np.argwhere(img > 0)
    if len(on) > s:
        drop = on[r.choice(len(on), size=len(on) - s, replace=False)]
        for a, b in drop:
            img[a, b] = 0.0
    return img.ravel().astype(float)


def scene_binary_structured(kind, seed):
    """(b) barcode / grid / text-like binary patterns, {0,1}."""
    r = np.random.default_rng(seed)
    img = np.zeros((SIDE, SIDE))
    if kind == "barcode":
        cols = r.random(SIDE) < 0.45
        img[:, cols] = 1.0
        img[: r.integers(4, 10), :] = 0.0                  # margin
    elif kind == "grid":
        pr, pc = int(r.integers(3, 6)), int(r.integers(3, 6))
        img[::pr, :] = 1.0; img[:, ::pc] = 1.0
    else:  # text-like: random 3x5 glyph blocks on a lattice
        for gr in range(2, SIDE - 6, 8):
            for gc in range(2, SIDE - 4, 6):
                glyph = (r.random((5, 3)) < 0.5).astype(float)
                img[gr:gr + 5, gc:gc + 3] = glyph
    return img.ravel().astype(float)


def scene_highcontrast_sparse(n_blobs, seed, contrast=10.0):
    """(c) few bright blobs on dark bg, NOT binary, contrast ~10:1, in [0,1]."""
    r = np.random.default_rng(seed)
    bg = 1.0 / contrast                                    # dark background 0.1
    img = np.full((SIDE, SIDE), bg)
    yy, xx = np.mgrid[0:SIDE, 0:SIDE]
    for _ in range(n_blobs):
        cy, cx = r.uniform(3, SIDE - 3), r.uniform(3, SIDE - 3)
        sig = r.uniform(0.9, 2.2)
        amp = r.uniform(0.7, 1.0) * (1.0 - bg)
        img += amp * np.exp(-((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * sig ** 2))
    return np.clip(img, 0.0, 1.0).ravel()


BRIDGE = "data/r63_bridge_scenes"
CONTROL_NAMES = ["bridge_contour_0", "bridge_contour_2", "bridge_twopop_1",
                 "bridge_twopop_3", "bridge_microtex_0", "bridge_microtex_2"]

# Interior clamp for CONTINUOUS classes (c, d): the committed probe clamps scenes
# to leave box-feasible headroom for null perturbations (certificate_r321.py uses
# [0.12,0.72]/[0.15,1.0]; bars_dose/t3 use [0.08,0.92]).  Natural scenes with
# pixels pinned at exactly 1.0 have ZERO upward headroom, which would give a 0-size
# fiber step.  We adopt the [0.08,0.92] convention verbatim.  Binary classes are
# left exactly {0,1} (handled by the identity, no box step needed).
CLAMP_LO, CLAMP_HI = 0.08, 0.92


def clamp_cont(x):
    return np.clip(x, CLAMP_LO, CLAMP_HI)


def load_control(name):
    return clamp_cont(op.load_scene(os.path.join(BRIDGE, name + ".npz")))


# =========================================================================== #
# 2.  ACHIEVABLE-PRECISION LADDER  (rel-se(p2) per sparse mask)
#     refreshed-W  : annealed FIM at the 3-level ridge design, se ~ 1/sqrt(B),
#                    C cancels at fixed budget (independent gates).
#     fixed-W      : sqrt(annealed^2 + quenched_floor^2); the floor ~ 1/sqrt(C_eff)
#                    does NOT average away (Pro appendix C, quenched_floor.json).
# =========================================================================== #
def annealed_relse_p2(v, B, C):
    """rel-se(p2) from the committed annealed FIM (order-2, p1 profiled)."""
    v = v[v > 0]
    if v.size < 2:
        return np.inf
    p1 = float(v.sum()); p2 = float((v ** 2).sum())
    if p2 <= 0:
        return np.inf
    a3, G3 = sc.three_level_design(p1, budget_incident=B, C=C)
    I = sc.annealed_fim(a3, G3, v, C, order=2)
    try:
        se = float(np.sqrt(np.linalg.inv(I)[1, 1]))
    except np.linalg.LinAlgError:
        return np.inf
    return se / p2


def quenched_floor_relse_p2(v, C):
    """fixed-W G->inf quenched floor of rel-se(p2): analytic propagation of the
    exact S2.1 quench covariance Cov[r_i,r_j]=(r_{i+j}-r_i r_j)/C through the
    order-2 ridge fit p = A^+ (-log r).  (= quenched_floor.py, closed form.)"""
    v = v[v > 0]
    if v.size < 2:
        return np.inf
    p1 = float(v.sum()); p2 = float((v ** 2).sum())
    a3, _ = sc.three_level_design(p1, budget_incident=None, C=C)     # ridge levels
    r = sc.Q_exp(a3, v)
    Smat = a3[:, None] + a3[None, :]
    rij = np.array([[sc.Q_exp(np.array([s]), v)[0] for s in row] for row in Smat])
    Sig_r = (rij - np.outer(r, r)) / C
    A = sc.design_matrix(a3, 2)                                      # 3 x 2
    Apinv = np.linalg.pinv(A)                                        # 2 x 3
    Jac = Apinv @ np.diag(-1.0 / r)                                  # d(p1,p2)/dr
    Cov_p = Jac @ Sig_r @ Jac.T
    se_p2 = float(np.sqrt(max(Cov_p[1, 1], 0.0)))
    return se_p2 / p2


def relse_p2(v, B, C, refresh):
    ann = annealed_relse_p2(v, B, C)
    if refresh:
        return ann
    fl = quenched_floor_relse_p2(v, C)
    return float(np.sqrt(ann ** 2 + fl ** 2))


# the ladder rungs used across the whole margin table -------------------------
RUNGS = [
    dict(key="i_dead_C3600_x40_fixedW",   C=3600,  mult=40,  refresh=False,
         label="(i) frozen probe: fixed-W, C=3600, x40  [THE DEAD POINT]"),
    dict(key="ii_refresh_C3600_x40",      C=3600,  mult=40,  refresh=True,
         label="(ii) +speckle refresh, C=3600, x40"),
    dict(key="iii_refresh_C14400_x40",    C=14400, mult=40,  refresh=True,
         label="(iii) +C=14400 (refresh), x40"),
    dict(key="iii_fixedW_C57600_x40",     C=57600, mult=40,  refresh=False,
         label="(iii) +C=57600 (fixed-W), x40"),
    dict(key="iv_refresh_C3600_x400",     C=3600,  mult=400, refresh=True,
         label="(iv) +budget x400 (refresh), C=3600  [>10x over R32 dose bar]"),
    dict(key="v_best_C57600_x400_refresh", C=57600, mult=400, refresh=True,
         label="(v) BEST ladder rung: C=57600, x400, refresh  [>10x over dose bar]"),
]


# =========================================================================== #
# 3.  FIBER AMBIGUITY  +  AGGREGATED FISHER-WEIGHTED d'   (R31 SS5.5)
# =========================================================================== #
def per_mask_stats(x, xp, M):
    """Per sparse mask k: dp2_k = M_k((xp^2 - x^2)), p2_k, n_eff_k, f_k=|dp2|/p2."""
    dp2, p2s, neff, fk, vlist = [], [], [], [], []
    x2, xp2 = x ** 2, xp ** 2
    for k in range(M.shape[0]):
        m = M[k]
        va = (m * x); va = va[va > 0]
        if va.size < 2:
            continue
        p1 = float(va.sum()); p2 = float((va ** 2).sum())
        if p2 <= 0:
            continue
        d = float(np.sum(m * (xp2 - x2)))
        dp2.append(d); p2s.append(p2); neff.append(p1 ** 2 / p2)
        fk.append(abs(d) / p2); vlist.append(va)
    return (np.array(dp2), np.array(p2s), np.array(neff),
            np.array(fk), vlist)


def aggregated_dprime(x, xp, M, B, C, refresh):
    """d'^2 = sum_k dp2_k^2 * I_p2,k  with I_p2,k = 1/se_p2,k(rung)^2 (all masks).
    Also returns median f over the best-quartile masks (by |dp2|/p2)."""
    dp2, p2s, neff, fk, vlist = per_mask_stats(x, xp, M)
    if dp2.size == 0:
        return 0.0, 0.0, 0
    d2 = 0.0
    used = 0
    for d, p2, v in zip(dp2, p2s, vlist):
        rel = relse_p2(v, B * B_LIN, C, refresh)
        se = rel * p2
        if not np.isfinite(se) or se <= 0:
            continue
        d2 += (d / se) ** 2
        used += 1
    # best-quartile f summary (reconstruction uses all masks; f-summary uses top 25%)
    q = max(1, int(np.ceil(0.25 * fk.size)))
    fq = np.sort(fk)[::-1][:q]
    return float(np.sqrt(d2)), float(np.median(fq)), used


def box_scale(x, hd, lo=0.0, hi=1.0):
    """Largest box-feasible step (no clip): x + alpha*h stays in [lo,hi]."""
    h = hd / (np.abs(hd).max() + 1e-30)
    pos = h > 1e-12; neg = h < -1e-12
    amax = np.inf
    if pos.any():
        amax = min(amax, float(np.min((hi - x[pos]) / h[pos])))
    if neg.any():
        amax = min(amax, float(np.min((lo - x[neg]) / h[neg])))
    if not np.isfinite(amax):
        return np.zeros_like(h)
    return 0.9 * amax * h


def project_ker(c):
    """Project c onto ker(L_BANK) = ker([M_D; M_S]) (re-projected for cleanliness)."""
    h = c - L_PINV @ (L_BANK @ c)
    h = h - L_PINV @ (L_BANK @ h)
    return h


def best_fiber_pair(x, seed, ncand=60):
    """Most FAVORABLE box-feasible fiber h in ker([M_D;M_S]): the pair (x, x+h) has
    IDENTICAL buckets on the full linear bank the equal-time comparator uses, so
    MOLT's p2 is the ONLY channel that can separate them.  Rank candidates by the
    aggregated d' at a fixed reference rung (x40, refreshed, C_base)."""
    r = np.random.default_rng(seed)
    lo, hi = float(x.min()), float(x.max())
    lo = max(0.0, lo); hi = min(1.0, hi)
    cands = [x - x.mean()]                                          # scene-aligned
    ii, jj = np.meshgrid(np.arange(SIDE), np.arange(SIDE), indexing="ij")
    cands.append((((ii + jj) % 2) * 2.0 - 1.0).ravel() * x)        # sign-textured x
    for _ in range(ncand):
        cands.append(r.standard_normal(P) * (x > x.mean()))        # bright-support noise
    best = (-1.0, None, 0.0)
    for c in cands:
        hd = project_ker(c)
        if np.abs(hd).max() < 1e-9:
            continue
        h = box_scale(x, hd, 0.0, 1.0)
        if np.abs(h).max() < 1e-12:
            continue
        xp = x + h
        # verify linear buckets identical (the fiber property)
        if float(np.max(np.abs(L_BANK @ xp - L_BANK @ x))) > 1e-8 * (np.abs(L_BANK @ x).max() + 1e-12):
            continue
        d, fmed, _ = aggregated_dprime(x, xp, M_S, 40, C_BASE, True)
        if d > best[0]:
            best = (d, xp, fmed)
    return best[1]


# =========================================================================== #
# 4.  BINARY-TRICK CHECK  (x^2 = x makes MOLT's p2 channel redundant)
# =========================================================================== #
def binary_trick_check(x_binary, seed=7):
    """For binary x: (a) MOLT's measured p2 = M_S(x^2) equals the sparse LINEAR
    bucket M_S(x) bit-for-bit; (b) hence for ANY second binary scene x', the MOLT
    channel difference Delta_p2 = M_S(x'^2 - x^2) equals the sparse LINEAR bucket
    difference M_S(x'-x) exactly -> the quadratic channel carries NO information
    beyond ordinary linear GI with bank [M_D;M_S]; (c) on a known support the R32.1
    quadratic rows 2 M_S|S diag(x_S) add ZERO rank once x^2=x is imposed."""
    assert set(np.unique(x_binary)).issubset({0.0, 1.0})
    p2_meas = M_S @ (x_binary ** 2)
    p1_sparse = M_S @ x_binary
    id_err = float(np.max(np.abs(p2_meas - p1_sparse)))              # ~0 exact

    r = np.random.default_rng(seed)
    xp = x_binary.copy()                                            # a 2nd binary scene
    on = np.flatnonzero(x_binary > 0); off = np.flatnonzero(x_binary == 0)
    if on.size and off.size:                                        # random 0<->1 swaps
        nsw = min(on.size, off.size, 30)
        xp[r.choice(on, nsw, replace=False)] = 0.0
        xp[r.choice(off, nsw, replace=False)] = 1.0
    dp2 = M_S @ (xp ** 2 - x_binary ** 2)
    dp1_sparse = M_S @ (xp - x_binary)
    molt_specific = float(np.max(np.abs(dp2 - dp1_sparse)))          # ~0 exact

    # known-support R32.1 rank: linear joint vs MOLT joint on a support S
    s = 120
    S = np.sort(r.choice(P, size=s, replace=False))
    xS = x_binary[S]
    MD_S, MS_S = M_D[:, S], M_S[:, S]
    Jq = 2 * MS_S * xS[None, :]
    rk_lin = int(np.linalg.matrix_rank(np.vstack([MD_S, MS_S])))
    rk_molt = int(np.linalg.matrix_rank(np.vstack([MD_S, MS_S, Jq])))
    return dict(p2_equals_sparse_p1_maxerr=id_err,
                molt_specific_dp2_maxerr=molt_specific,
                support_s=s, rank_linear_joint=rk_lin, rank_molt_joint=rk_molt,
                quad_rows_add_rank=int(rk_molt - rk_lin),
                verdict="MOLT p2 channel is a bit-copy of the sparse linear bucket for "
                        "binary scenes (x^2=x); nonlinear channel REDUNDANT -> binary "
                        "scenes need the x^2=x identity + ordinary linear GI, not MOLT.")


# =========================================================================== #
# 5.  COMPETITOR CHECK  (class a, s=80): extra linear masks vs MOLT patterns
# =========================================================================== #
def competitor_check(s=80):
    """To close a KNOWN support of size s (generic random masks, exact):
       - pure linear needs s independent rows = s patterns  (K_extra = s - K_D).
       - MOLT (continuous scene): K_D dense + K_S sparse patterns give rank up to
         K_D + 2 K_S; minimal patterns = K_D + ceil((s-K_D)/2).
       - MOLT (BINARY scene): quad rows redundant (x^2=x) -> rank only K_D+K_S ->
         needs s patterns, SAME as linear (no saving).
       - blind CS (unknown support): ~2s (exact-ish) to ~s log2(N/s) (random)."""
    K_extra_linear = s - K_D                                        # additional linear masks
    lin_patterns = s
    molt_patterns_cont = K_D + int(np.ceil((s - K_D) / 2))          # uses 2 rows/sparse pattern
    molt_patterns_bin = s                                           # quad redundant
    cs_2s = 2 * s
    cs_slog = int(np.ceil(s * np.log2(P / s)))
    # per-sparse-mask MOLT dose multiple (n_eff^2 law, favorable sparse n_eff)
    dose_mult_per_sparse = 40.0                                     # R32 bar-2 cap (x40)
    extra_dose_total = (molt_patterns_cont - K_D) * dose_mult_per_sparse
    return dict(support_s=s, K_D=K_D, K_S=K_S,
                pure_linear_patterns=lin_patterns,
                extra_linear_masks_needed=K_extra_linear,
                molt_patterns_continuous=molt_patterns_cont,
                molt_patterns_binary=molt_patterns_bin,
                pattern_ratio_molt_over_linear_cont=float(molt_patterns_cont / lin_patterns),
                pattern_ratio_molt_over_linear_bin=float(molt_patterns_bin / lin_patterns),
                patterns_saved_cont=lin_patterns - molt_patterns_cont,
                molt_extra_dose_xLIN=float(extra_dose_total),
                blind_cs_2s=cs_2s, blind_cs_slog2=cs_slog,
                verdict="MOLT saves %d of %d patterns on a continuous s=%d support, "
                        "at ~%.0fx extra linear dose; on a BINARY support it saves 0 "
                        "patterns (quad rows redundant). MOLT wins only if pattern "
                        "slots are the binding resource AND dose is free."
                        % (lin_patterns - molt_patterns_cont, lin_patterns, s,
                           extra_dose_total))


# =========================================================================== #
# 6.  DL-TANGENT LEVER  (learned-manifold proxy; coordinator addendum)
# =========================================================================== #
def class_tangent_basis(gen_fn, d, n_samples=300, seed=99):
    """Proxy for a LEARNED manifold tangent: leading d PCA directions (pixel-space
    right singular vectors) of an n_samples ensemble of the scene class.  HONEST
    PROXY: a linear/global PCA of a class ensemble, NOT a true nonlinear learned
    manifold; it is a smoothness+class-statistics stand-in."""
    r = np.random.default_rng(seed)
    X = np.stack([gen_fn(r.integers(1 << 30)) for _ in range(n_samples)])   # n x P
    Xc = X - X.mean(0, keepdims=True)
    # right singular vectors of centered ensemble = principal tangent directions
    _, _, Vt = np.linalg.svd(Xc, full_matrices=False)
    return Vt[:d].T                                                 # P x d, orthonormal


def dl_tangent_block(x, B_x, rungs, seed=11):
    """R32.1 closure on a d-dim prior tangent B_x (P x d) at scene x.
      Jlin  = [M_D B_x ; M_S B_x]                       (linear-only tangent Jac)
      Jjnt  = [M_D B_x ; M_S B_x ; 2 M_S diag(x) B_x]   (R32.1 tangent Jac)
    (i) closure: does rank(Jjnt)=d?  does rank(Jlin)=d (=> MOLT redundant)?
    (ii) margin: aggregated MOLT d' for a box-feasible step along the LEAST
         linearly-visible tangent direction (smallest Jlin singular vector),
         at each precision rung, /3."""
    d = B_x.shape[1]
    MDB = M_D @ B_x; MSB = M_S @ B_x
    Jlin = np.vstack([MDB, MSB])                                    # 102 x d
    Jq = (2.0 * M_S * x[None, :]) @ B_x                             # K_S x d
    Jjnt = np.vstack([MDB, MSB, Jq])                               # 153 x d
    rk_lin = int(np.linalg.matrix_rank(Jlin))
    rk_jnt = int(np.linalg.matrix_rank(Jjnt))
    sv_lin = np.linalg.svd(Jlin, compute_uv=False)
    sv_jnt = np.linalg.svd(Jjnt, compute_uv=False)
    smin_lin = float(sv_lin[min(rk_lin, d) - 1]) if d <= Jlin.shape[0] else 0.0
    smin_jnt = float(sv_jnt[min(rk_jnt, d) - 1]) if d <= Jjnt.shape[0] else 0.0

    # least linearly-visible tangent direction (smallest right-sing vector of Jlin)
    _, _, VtL = np.linalg.svd(Jlin, full_matrices=False)
    theta = VtL[-1]                                                 # d-vector
    u = B_x @ theta                                                 # P-vector, unit-ish
    h = box_scale(x, u, 0.0, 1.0)                                   # box-feasible step
    xp = x + h
    margins = {}
    for rg in rungs:
        dpr, _, _ = aggregated_dprime(x, xp, M_S, rg["mult"], rg["C"], rg["refresh"])
        margins[rg["key"]] = float(dpr / 3.0)
    return dict(d=d, rank_linear_tangent=rk_lin, rank_joint_tangent=rk_jnt,
                linear_closes=bool(rk_lin >= d), joint_closes=bool(rk_jnt >= d),
                sigma_min_linear=smin_lin, sigma_min_joint=smin_jnt,
                least_visible_dp2_scale=float(np.max(np.abs(M_S @ (xp ** 2 - x ** 2)))),
                margin_by_rung=margins)


# =========================================================================== #
# DRIVER
# =========================================================================== #
def class_margin_table(name, scenes, is_binary):
    """Median margin per rung over the class's scene instances."""
    print("\n" + "=" * 74 + "\nCLASS %s  (%d instances, binary=%s)"
          % (name, len(scenes), is_binary))
    rows = {rg["key"]: [] for rg in RUNGS}
    fmed_ref = []
    per_scene = []
    for i, x in enumerate(scenes):
        if is_binary:
            # binary: MOLT-specific Delta_p2 == sparse linear bucket diff -> 0 gain
            per_scene.append(dict(instance=i, margin_by_rung={rg["key"]: 0.0 for rg in RUNGS},
                                  f_median_bestquartile=0.0,
                                  note="binary: MOLT p2 channel == sparse linear bucket "
                                       "(x^2=x); MOLT-specific margin identically 0"))
            for rg in RUNGS:
                rows[rg["key"]].append(0.0)
            fmed_ref.append(0.0)
            continue
        xp = best_fiber_pair(x, seed=63_1000 + i)
        if xp is None:
            continue
        _, fmed, _ = aggregated_dprime(x, xp, M_S, 40, C_BASE, True)
        fmed_ref.append(fmed)
        mb = {}
        for rg in RUNGS:
            dpr, _, used = aggregated_dprime(x, xp, M_S, rg["mult"], rg["C"], rg["refresh"])
            mb[rg["key"]] = float(dpr / 3.0)
            rows[rg["key"]].append(dpr / 3.0)
        per_scene.append(dict(instance=i, margin_by_rung=mb,
                              f_median_bestquartile=float(fmed)))
        print("  inst %d: f_med(bestQ)=%.3f  margin[dead]=%.2f  margin[refresh x40]=%.2f"
              "  margin[best x400]=%.2f"
              % (i, fmed, mb["i_dead_C3600_x40_fixedW"],
                 mb["ii_refresh_C3600_x40"], mb["v_best_C57600_x400_refresh"]))
    med = {rg["key"]: (float(np.median(rows[rg["key"]])) if rows[rg["key"]] else 0.0)
           for rg in RUNGS}
    return dict(median_margin_by_rung=med,
                f_median_bestquartile=float(np.median(fmed_ref)) if fmed_ref else 0.0,
                per_scene=per_scene)


def main():
    t0 = time.time()
    print("MOLT REVIVAL DESK CALCULATION")
    print("banks: M_D=%dx%d dense50%%  M_S=%dx%d %d-support  combined-linear rank=%d"
          % (M_D.shape[0], P, M_S.shape[0], P, SUPPORT,
             np.linalg.matrix_rank(L_BANK)))

    # ---- build class instances --------------------------------------------- #
    a_scenes = ([scene_binary_sparse(40, 63_1100 + i) for i in range(3)]
                + [scene_binary_sparse(80, 63_1200 + i) for i in range(3)]
                + [scene_binary_sparse(150, 63_1300 + i) for i in range(3)])
    b_scenes = ([scene_binary_structured("barcode", 63_1400 + i) for i in range(2)]
                + [scene_binary_structured("grid", 63_1420 + i) for i in range(2)]
                + [scene_binary_structured("text", 63_1440 + i) for i in range(2)])
    c_scenes = [clamp_cont(scene_highcontrast_sparse(int(nb), 63_1500 + i))
                for i, nb in enumerate([3, 4, 5, 6, 8, 10])]
    d_scenes = [load_control(n) for n in CONTROL_NAMES]

    # sanity: contrast / binarity
    report["class_descriptions"] = dict(
        a_binary_sparse=dict(n=len(a_scenes), supports=[40, 80, 150],
            binary=bool(all(set(np.unique(s)).issubset({0.0, 1.0}) for s in a_scenes))),
        b_binary_structured=dict(n=len(b_scenes),
            binary=bool(all(set(np.unique(s)).issubset({0.0, 1.0}) for s in b_scenes))),
        c_highcontrast_sparse=dict(n=len(c_scenes),
            binary=False, mean_max_over_min=float(np.median(
                [s[s > 0].max() / s[s > 0].min() for s in c_scenes]))),
        d_control_natural=dict(n=len(d_scenes), names=CONTROL_NAMES, binary=False))

    # ---- 4 margin tables --------------------------------------------------- #
    report["margin_table"] = {
        "a_binary_sparse": class_margin_table("a_binary_sparse", a_scenes, True),
        "b_binary_structured": class_margin_table("b_binary_structured", b_scenes, True),
        "c_highcontrast_sparse": class_margin_table("c_highcontrast_sparse", c_scenes, False),
        "d_control_natural": class_margin_table("d_control_natural", d_scenes, False),
    }

    # ---- binary-trick check ------------------------------------------------ #
    print("\n" + "=" * 74 + "\nBINARY-TRICK CHECK (x^2=x)")
    bt = binary_trick_check(a_scenes[3])                            # an s=80 binary scene
    print("  p2_meas == sparse-p1  max err = %.2e" % bt["p2_equals_sparse_p1_maxerr"])
    print("  MOLT-specific Delta_p2 (== sparse linear bucket diff) max err = %.2e"
          % bt["molt_specific_dp2_maxerr"])
    print("  known-support s=%d: rank(linear joint)=%d  rank(MOLT joint)=%d  quad adds %d"
          % (bt["support_s"], bt["rank_linear_joint"], bt["rank_molt_joint"],
             bt["quad_rows_add_rank"]))
    report["binary_trick_check"] = bt

    # ---- competitor check -------------------------------------------------- #
    print("\n" + "=" * 74 + "\nCOMPETITOR CHECK (class a, s=80)")
    comp = competitor_check(80)
    print("  pure-linear needs %d patterns (%d extra); MOLT(cont) %d; MOLT(binary) %d"
          % (comp["pure_linear_patterns"], comp["extra_linear_masks_needed"],
             comp["molt_patterns_continuous"], comp["molt_patterns_binary"]))
    print("  pattern ratio MOLT/linear: cont=%.3f  binary=%.3f  | blind CS 2s=%d, slog2=%d"
          % (comp["pattern_ratio_molt_over_linear_cont"],
             comp["pattern_ratio_molt_over_linear_bin"],
             comp["blind_cs_2s"], comp["blind_cs_slog2"]))
    report["competitor_check"] = comp

    # ---- DL-tangent lever (coordinator addendum) --------------------------- #
    print("\n" + "=" * 74 + "\nDL-TANGENT LEVER (learned-manifold proxy: class-ensemble PCA)")
    dl = {}
    gen_c = lambda seed: scene_highcontrast_sparse(int(np.random.default_rng(seed).integers(3, 11)), seed)
    tangent_gens = {
        "c_highcontrast_sparse": (gen_c, c_scenes[2]),
        "d_control_natural": (
            # natural proxy ensemble: smoothed 1/f random fields (honest smoothness proxy,
            # since the frozen bridge generators live outside the repo)
            lambda seed: _smoothed_field(seed), d_scenes[0]),
    }
    for cname, (gen_fn, xrep) in tangent_gens.items():
        dl[cname] = {}
        for d in (20, 50, 100):
            B_x = class_tangent_basis(gen_fn, d, n_samples=300, seed=770 + d)
            blk = dl_tangent_block(xrep, B_x, RUNGS)
            dl[cname][f"d{d}"] = blk
            print("  %-22s d=%3d  lin_closes=%s joint_closes=%s  smin_lin=%.2e smin_jnt=%.2e"
                  "  margin[dead]=%.2f margin[best]=%.2f"
                  % (cname, d, blk["linear_closes"], blk["joint_closes"],
                     blk["sigma_min_linear"], blk["sigma_min_joint"],
                     blk["margin_by_rung"]["i_dead_C3600_x40_fixedW"],
                     blk["margin_by_rung"]["v_best_C57600_x400_refresh"]))
    report["dl_tangent_lever"] = dict(
        proxy_note="tangent B_x = leading-d PCA of a 300-scene class ensemble "
                   "(high-contrast: own generator; natural: smoothed 1/f random-field "
                   "proxy, since the frozen bridge generators are outside the repo). "
                   "HONEST: a linear/global manifold stand-in, not a true learned "
                   "nonlinear manifold.",
        blocks=dl)

    # ---- verdict synthesis ------------------------------------------------- #
    honest = ["c_highcontrast_sparse", "d_control_natural"]
    dose_bar_rungs = ["i_dead_C3600_x40_fixedW", "ii_refresh_C3600_x40",
                      "iii_refresh_C14400_x40", "iii_fixedW_C57600_x40"]  # <= x40 honest dose
    best_honest_dose = {}
    best_any = {}
    for c in honest:
        mb = report["margin_table"][c]["median_margin_by_rung"]
        best_honest_dose[c] = max(mb[k] for k in dose_bar_rungs)
        best_any[c] = max(mb.values())
    dl_best = 0.0
    for c in honest:
        for d in ("d20", "d50", "d100"):
            dl_best = max(dl_best, max(report["dl_tangent_lever"]["blocks"][c][d]
                                       ["margin_by_rung"].values()))
    revive_honest_dose = any(v >= 10 for v in best_honest_dose.values())
    revive_any_rung = any(v >= 10 for v in best_any.values()) or dl_best >= 10
    report["verdict"] = dict(
        best_median_margin_within_R32_dose_bar_x40=best_honest_dose,
        best_median_margin_any_ladder_rung_incl_x400=best_any,
        dl_tangent_best_margin=float(dl_best),
        revive_at_honest_dose_x40=bool(revive_honest_dose),
        revive_only_beyond_dose_bar_x400=bool(revive_any_rung and not revive_honest_dose),
        DECISION="REVIVE" if revive_honest_dose else "STAY_DEAD")
    print("\n" + "=" * 74)
    print("VERDICT: %s" % report["verdict"]["DECISION"])
    print("  best median margin within R32 x40 dose bar:", {k: round(v, 2) for k, v in best_honest_dose.items()})
    print("  best median margin any rung (incl x400):   ", {k: round(v, 2) for k, v in best_any.items()})
    print("  DL-tangent best margin: %.2f" % dl_best)
    print("  elapsed %.1fs" % (time.time() - t0))

    with open(os.path.join(OUT, "revival_calc.json"), "w") as fh:
        json.dump(report, fh, indent=2, default=float)
    print("wrote", os.path.join(OUT, "revival_calc.json"))


def _smoothed_field(seed):
    """Smoothed 1/f-ish random field in [0,1] (natural-scene smoothness proxy)."""
    r = np.random.default_rng(int(seed) & 0x7FFFFFFF)
    f = r.standard_normal((SIDE, SIDE))
    F = np.fft.fft2(f)
    ky = np.fft.fftfreq(SIDE)[:, None]; kx = np.fft.fftfreq(SIDE)[None, :]
    k = np.sqrt(ky ** 2 + kx ** 2) + 1e-3
    F = F / k                                                       # 1/f envelope
    img = np.real(np.fft.ifft2(F))
    img = (img - img.min()) / (img.ptp() + 1e-12)
    return img.ravel()


if __name__ == "__main__":
    main()
