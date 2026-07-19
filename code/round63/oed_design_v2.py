"""Dead-time-aware OED v2 -- chain-rule renewal information atoms (R10 ruling).

Upgrades oed_design.py (v1, geometry-only a a^T atoms) per the frozen R10
ruling  docs/ROUND63_GPT_ROUND10_RULING_RAW.md :

Q1  chain-rule atom      H(a) = nu * J(rho(a), nu) * g g^T,
                         g = a / max(a . xhat, floor)
    (the log-rate chain rule at the frozen coarse estimate xhat; J(rho,nu) is
    the exact per-slot renewal count information about log-lambda).
G2  variable load        atom = (shape, translation, p, rho_level);
                         row scaled so predicted load  a . xhat * Phi * tau
                         = rho_level  (convention: Phi*tau = 1, dark = 0, and
                         xhat normalized so a full-field row has load 1 --
                         the Study-2 relative-irradiance units).
Q4.1 V0 term             objective  log det( V0 + M_rows * sum_i xi_i H_i )
                         with V0 = exact information of a frozen deterministic
                         52-row multiscale pre-scan (built here, documented).
Q4.3 guards              G3 peak a_j <= 64 (drop), G4 weights clipped to
                         [1/4, 4] + p in {0,1} only, G5 post-rounding
                         per-pixel dose +/-5% (report only), G6 generalized-
                         eigenvalue spectral floor hook (report + alpha mix).

Because g is invariant to the row's overall scale, the SAME geometry at two
load levels differs only by the scalar nu*J(rho,nu).  J is increasing in rho
below the ridge rho*(nu), so an UNCONSTRAINED optimizer prefers the highest
surviving load; load mixing can only come from G3 drops (high-load atoms on
dim supports) or from the ruling's design-average-load equality constraint
(G2, "design-average load exactly 0.60"), which is NOT implemented here --
the selftest prints an explicit WARNING with this diagnosis if the optimizer
collapses to a single load (anticipated behavior, not hidden).

Solver: multiplicative (Titterington) weight updates  xi_i <- xi_i * d_i/dbar
-- monotone for D-optimality and far faster than vanilla Fedorov-Wynn on a
30k-atom dictionary (v1's FW stalled at gap ~0.96 after 400 iters) -- with a
log-det backtracking guard, support pruning, and Frank-Wolfe atom injection
on stall.  KW certificate with the V0 offset:  at the optimum
max_a tr[V^-1 H(a)] = sum_i xi_i tr[V^-1 H_i]  (= tr[V^-1 (V - V0)]/M_rows),
so   gap = max_a d(a) / dbar - 1.

J source: the committed exact table results/round63_theory/fig_a_sweep.npz
(convention verified against code/round63/fig_a_ridge_map.py:
 J = exact_fisher_analytic(rho, nu, tau=1)/nu, per-slot info about log-lambda),
bilinearly interpolated in (log rho, log nu); a numpy-only exact-PMF fallback
is included and cross-checked in the selftest.

ADD-only file; imports numpy + stdlib + sibling module oed_design (v1).
"""

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import oed_design as v1                      # noqa: E402  (ADD-only reuse)

ROOT = os.path.dirname(os.path.dirname(HERE))
J_NPZ = os.path.join(ROOT, "results", "round63_theory", "fig_a_sweep.npz")
J_CSV = os.path.join(ROOT, "results", "round63_theory", "fisher_ridge.csv")


# ========================================================================== #
#  J(rho, nu): exact per-slot renewal count information about log-lambda      #
# ========================================================================== #
class JTable:
    """Bilinear interpolation of the committed exact sweep in (log rho, log nu)."""

    source = "table:results/round63_theory/fig_a_sweep.npz"

    def __init__(self, npz_path=J_NPZ):
        z = np.load(npz_path)
        self.lr = np.log(z["rho_grid"])          # (481,)
        self.ln = np.log(z["nu_grid"])           # (49,)
        self.J = np.asarray(z["J"])              # (481, 49)

    def __call__(self, rho, nu):
        lr = np.log(np.atleast_1d(np.asarray(rho, dtype=float)))
        ln = math.log(float(nu))
        lr = np.clip(lr, self.lr[0], self.lr[-1])
        ln = min(max(ln, self.ln[0]), self.ln[-1])
        i = np.clip(np.searchsorted(self.lr, lr) - 1, 0, self.lr.size - 2)
        j = min(max(int(np.searchsorted(self.ln, ln)) - 1, 0), self.ln.size - 2)
        tr_ = (lr - self.lr[i]) / (self.lr[i + 1] - self.lr[i])
        tn = (ln - self.ln[j]) / (self.ln[j + 1] - self.ln[j])
        out = ((1 - tr_) * (1 - tn) * self.J[i, j]
               + tr_ * (1 - tn) * self.J[i + 1, j]
               + (1 - tr_) * tn * self.J[i, j + 1]
               + tr_ * tn * self.J[i + 1, j + 1])
        return out if np.ndim(rho) else float(out[0])


class JExactNumpy:
    """Numpy-only exact fallback: J from the renewal-count PMF via the score.

    Active-start non-paralyzable renewal counting in a window of nu slots
    (tau = 1, lam = rho):  G_m = P(Pois(z_m) >= m), z_m = rho*(nu-(m-1));
    p_m = G_m - G_{m+1};  dG_m/dtheta = m*PoisPMF(m; z_m)  (theta = log lam);
    I_theta = sum_m (Gdot_m - Gdot_{m+1})^2 / p_m ;  J = I_theta / nu.
    Same identities as physics.exact_fisher_analytic, without scipy.
    """

    source = "computed:numpy-exact-pmf-score"

    def _log_sf_pois(self, m, z):
        """log P(Pois(z) >= m) for scalar integer m >= 1, scalar z > 0."""
        hi = int(max(m + 60, z + 12.0 * math.sqrt(z) + 60))
        j = np.arange(m, hi + 1, dtype=np.float64)
        lt = j * math.log(z) - z - self._lf[m:hi + 1]
        mx = lt.max()
        return mx + math.log(np.exp(lt - mx).sum())

    def __call__(self, rho, nu):
        scal = np.ndim(rho) == 0
        rhos = np.atleast_1d(np.asarray(rho, dtype=float))
        nu = float(nu)
        m_max = int(nu) + 2
        zmax_all = float(rhos.max()) * nu
        n_lf = int(max(m_max, zmax_all + 12.0 * math.sqrt(max(zmax_all, 1.0)))) + 200
        self._lf = np.concatenate([[0.0],
                                   np.cumsum(np.log(np.arange(1, n_lf + 1)))])
        out = np.empty(rhos.shape)
        for a, r in enumerate(rhos):
            mm = np.arange(0, m_max + 2)
            t = nu - (mm - 1.0)
            z = r * np.maximum(t, 0.0)
            lG = np.full(mm.size, -np.inf)
            lG[0] = 0.0                                   # G_0 = 1
            for m in range(1, mm.size):
                if z[m] > 0:
                    lG[m] = self._log_sf_pois(m, z[m])
            with np.errstate(invalid="ignore"):
                lgd = np.where((mm >= 1) & (z > 0),
                               np.log(np.maximum(mm, 1))
                               + mm * np.log(np.maximum(z, 1e-300)) - z
                               - self._lf[np.minimum(mm, n_lf)],
                               -np.inf)
            Gd = np.exp(lgd)
            num = (Gd - np.roll(Gd, -1))[:-1]
            # p_m in log domain: lG_m + log(1 - exp(lG_{m+1} - lG_m))
            with np.errstate(invalid="ignore", divide="ignore"):
                lp = lG[:-1] + np.log(-np.expm1(lG[1:] - lG[:-1]))
            ok = np.isfinite(lp) & (np.abs(num) > 0)
            with np.errstate(divide="ignore"):
                lt = 2.0 * np.log(np.abs(num[ok])) - lp[ok]
            keep = lt > -745.0
            if not keep.any():
                out[a] = 0.0
                continue
            mx = float(lt[keep].max())
            out[a] = math.exp(mx) * float(np.exp(lt[keep] - mx).sum()) / nu
        return float(out[0]) if scal else out


def get_J_source():
    """Committed exact table if present (preferred), else numpy-only exact."""
    if os.path.exists(J_NPZ):
        return JTable()
    return JExactNumpy()


# ========================================================================== #
#  Frozen 52-row multiscale pre-scan (V0)                                     #
# ========================================================================== #
def prescan_52_supports(side=32):
    """Deterministic balanced multiscale 52-row pre-scan (documented choice).

    * 32 fine rows  : 4x4 block indicators on the even-parity half
                      (checkerboard) of the 8x8 grid of 4x4 blocks;
    * 16 mid rows   : ALL 8x8 block indicators (4x4 grid of them);
    *  4 coarse rows: the four 16x16 full-field quadrants.
    Total 32 + 16 + 4 = 52 rows.  Three dyadic scales; the mid + coarse rows
    tile the frame completely, the fine rows sample it in a checkerboard.
    """
    rows = []
    for bi in range(8):
        for bj in range(8):
            if (bi + bj) % 2 == 0:
                ys, xs = 4 * bi, 4 * bj
                rows.append([(ys + dy) * side + (xs + dx)
                             for dy in range(4) for dx in range(4)])
    for bi in range(4):
        for bj in range(4):
            ys, xs = 8 * bi, 8 * bj
            rows.append([(ys + dy) * side + (xs + dx)
                         for dy in range(8) for dx in range(8)])
    for bi in range(2):
        for bj in range(2):
            ys, xs = 16 * bi, 16 * bj
            rows.append([(ys + dy) * side + (xs + dx)
                         for dy in range(16) for dx in range(16)])
    assert len(rows) == 52
    return [np.asarray(r, dtype=np.int64) for r in rows]


def build_V0(xhat, nu, arm_rho, Jfun, side):
    """Exact chain-rule information of the 52-row pre-scan at the arm's rho.

    Each pre-scan row is a binary indicator scaled to load arm_rho; its
    chain-rule direction is g = a/(a.xhat) (scale-invariant), and each row
    contributes one exposure:  V0 = sum_rows nu*J(arm_rho,nu) * g g^T.
    """
    n = side * side
    V0 = np.zeros((n, n))
    Jv = float(Jfun(arm_rho, nu))
    c = nu * Jv
    for idx in prescan_52_supports(side):
        load = xhat[idx].sum()
        g = np.zeros(n)
        g[idx] = 1.0 / load                     # g . xhat = 1 exactly
        V0[np.ix_(idx, idx)] += c * np.outer(g[idx], g[idx])
    return V0, Jv


# ========================================================================== #
#  Guard G6 hook: relative spectral floor on the top-r subspace               #
# ========================================================================== #
def g6_spectral_floor(V_oed, V_fix, r=200, floor=0.5):
    """Generalized-eigenvalue floor  V_OED >= floor * V_fix  on V_fix's
    top-r eigenspace (subspace surrogate of ruling G6).  Returns dict with
    mu_min, PASS flag, and (if FAIL) the largest alpha such that
    alpha*V_OED + (1-alpha)*V_fix still satisfies the floor there.
    """
    evals, evecs = np.linalg.eigh(V_fix)
    U = evecs[:, -r:]                                   # top-r eigenvectors
    lam = evals[-r:]
    A = U.T @ V_oed @ U
    # B = U^T V_fix U = diag(lam) exactly; whiten by 1/sqrt(lam)
    s = 1.0 / np.sqrt(lam)
    C = A * s[None, :] * s[:, None]
    mu = np.linalg.eigvalsh(C)
    mu_min = float(mu[0])
    ok = mu_min >= floor - 1e-12
    out = {"r": r, "floor": floor, "mu_min": mu_min, "pass": bool(ok)}
    if not ok:
        # projected eigenvalues of the mixture are alpha*mu_i + (1-alpha)
        out["alpha_max"] = float((1.0 - floor) / (1.0 - mu_min))
    return out


# ========================================================================== #
#  Geometry dictionary (G3/G4 guards applied)                                 #
# ========================================================================== #
def _build_geometries(xhat, shapes, powers, side, w_clip, peak_cap,
                      rho_levels, cv_guard=1.0):
    """All (shape, translation, p) chain-rule directions g (load-invariant).

    Per geometry: support IDX (k,), direction values GVAL (k,) with
    g.xhat = 1 exactly, and per-load admissibility ALLOW[geo, load]
    (G3: peak a_j = rho*max(g) <= peak_cap).
    G4: w = clip((xhat_on/mean(xhat_on))^p, 1/4, 4) renormalized to mean 1
    on the support; p is restricted by the caller to {0, 1}.
    Belt guard (v1 frozen lesson): drop geometries with sd(w)/mean(w) > 1.
    """
    IDXp, GVALp, namep, tp, pp = [], [], [], [], []
    t_all = np.arange(side * side)
    n_cv_drop = 0
    clip_frac_num = 0
    clip_frac_den = 0
    for name, offs in shapes.items():
        sidx = v1._shape_support(offs, side)            # (n, k)
        xh_on = xhat[sidx]
        for p in powers:
            if p not in (0, 1):                          # G4: p in {0,1} only
                raise ValueError("G4: powers restricted to {0,1}, got %r" % (p,))
            w0 = (xh_on / xh_on.mean(axis=1, keepdims=True)) ** float(p)
            w1 = np.clip(w0, w_clip[0], w_clip[1])
            clip_frac_num += int((w1 != w0).sum())
            clip_frac_den += w0.size
            w = w1 / w1.mean(axis=1, keepdims=True)      # mean 1 on support
            cv = w.std(axis=1) / w.mean(axis=1)
            keep = cv <= cv_guard
            n_cv_drop += int((~keep).sum())
            load = (w * xh_on).sum(axis=1)               # raw row . xhat
            gval = w / load[:, None]                     # g values: g.xhat = 1
            IDXp.append(sidx[keep])
            GVALp.append(gval[keep])
            namep.extend([name] * int(keep.sum()))
            tp.append(t_all[keep])
            pp.append(np.full(int(keep.sum()), float(p)))
    IDX = np.concatenate(IDXp, axis=0)
    GVAL = np.concatenate(GVALp, axis=0)
    t_arr = np.concatenate(tp, axis=0)
    p_arr = np.concatenate(pp, axis=0)
    gmax = GVAL.max(axis=1)
    rho = np.asarray(rho_levels, dtype=float)
    ALLOW = gmax[:, None] * rho[None, :] <= peak_cap + 1e-12   # G3
    g3 = {"peak_cap": peak_cap,
          "dropped_per_load": {float(r): int((~ALLOW[:, l]).sum())
                               for l, r in enumerate(rho)},
          "n_geometries": int(IDX.shape[0])}
    g4 = {"w_clip": tuple(w_clip),
          "clipped_weight_fraction": clip_frac_num / max(clip_frac_den, 1),
          "cv_guard_dropped": n_cv_drop}
    return IDX, GVAL, namep, t_arr, p_arr, ALLOW, g3, g4


def _build_M(IDX, GVAL, u_geo, n):
    """M = sum_g u_g * g_g g_g^T for geometry weights u (scatter via bincount)."""
    nz = np.nonzero(u_geo)[0]
    if nz.size == 0:
        return np.zeros((n, n))
    I = IDX[nz]
    V = GVAL[nz]
    u = u_geo[nz]
    lin = (I[:, :, None] * n + I[:, None, :]).ravel()
    contrib = (u[:, None, None] * V[:, :, None] * V[:, None, :]).ravel()
    return np.bincount(lin, weights=contrib, minlength=n * n).reshape(n, n)


# ========================================================================== #
#  Main v2 solver                                                             #
# ========================================================================== #
def design_v2(xhat, shapes=None, powers=(0, 1), rho_levels=(0.3, 0.6, 1.0),
              nu=2000.0, M_rows=1024, max_iter=600, tol=0.1, arm_rho=0.6,
              peak_cap=64.0, w_clip=(0.25, 4.0), J_source=None,
              cert_every=10, verbose=False):
    """Chain-rule dead-time-aware D-optimal SPI design (R10 v2).

    Maximizes  log det( V0 + M_rows * sum_i xi_i H_i )  over the design
    simplex, H_i = nu*J(rho_i, nu) * g_i g_i^T  (Q1 chain-rule atoms),
    V0 = frozen 52-row pre-scan information at arm_rho.

    Returns dict: A (M_rows x n), atoms [(shape, t, p, rho, count)],
    gap, gap_traj [(iter, gap)], cert, guards, V (final), V0, xi arrays.
    """
    xhat = np.asarray(xhat, dtype=float).ravel()
    n = xhat.size
    side = int(round(np.sqrt(n)))
    if side * side != n:
        raise ValueError("xhat length must be a perfect square.")
    if xhat.min() < 0:
        raise ValueError("xhat must be non-negative.")
    if shapes is None:
        shapes = v1._default_shapes(side)
    Jfun = J_source if J_source is not None else get_J_source()

    rho = np.asarray(rho_levels, dtype=float)
    nuJ = np.array([nu * float(Jfun(r, nu)) for r in rho])   # atom info scale

    IDX, GVAL, names, t_arr, p_arr, ALLOW, g3, g4 = _build_geometries(
        xhat, shapes, powers, side, w_clip, peak_cap, rho_levels)
    G = IDX.shape[0]
    L = rho.size
    n_atoms = int(ALLOW.sum())

    V0, J_arm = build_V0(xhat, nu, arm_rho, Jfun, side)

    # ---- initial design: uniform over admissible atoms -------------------- #
    xi = np.where(ALLOW, 1.0, 0.0)
    xi /= xi.sum()

    def assemble(xi_):
        u_geo = M_rows * (xi_ * nuJ[None, :]).sum(axis=1)
        M = _build_M(IDX, GVAL, u_geo, n)
        Vr = V0 + M
        eps = 1e-9 * np.trace(Vr) / n
        Vr[np.diag_indices(n)] += eps
        return Vr, eps

    def full_sensitivity(Vinv):
        sub = Vinv[IDX[:, :, None], IDX[:, None, :]]         # (G, k, k)
        dgeo = np.einsum('ai,aij,aj->a', GVAL, sub, GVAL)     # g^T Vinv g
        d = dgeo[:, None] * nuJ[None, :]                      # tr[Vinv H]
        return np.where(ALLOW, d, -np.inf), dgeo

    V, eps = assemble(xi)
    sign, phi = np.linalg.slogdet(V)
    gap_traj = []
    gap = np.inf
    stall = 0
    last_gap = np.inf
    it = 0
    for it in range(1, max_iter + 1):
        Vinv = np.linalg.inv(V)
        d, dgeo = full_sensitivity(Vinv)
        dbar = float((xi * np.where(ALLOW, d, 0.0)).sum())
        dmax = float(d.max())
        gap = dmax / dbar - 1.0
        gap_traj.append((it, gap))
        if verbose and (it % 25 == 0 or it == 1):
            print("    [v2] it=%4d  gap=%.4f  logdet=%.3f" % (it, gap, phi),
                  flush=True)
        if gap < tol:
            break
        # ---- multiplicative update with monotone backtracking ------------- #
        beta = 1.0
        accepted = False
        for _ in range(6):
            xi_new = xi * np.where(ALLOW, np.maximum(d, 0.0) / dbar, 0.0) ** beta
            ssum = xi_new.sum()
            if ssum > 0:
                xi_new /= ssum
                V_new, eps = assemble(xi_new)
                s_new, phi_new = np.linalg.slogdet(V_new)
                if s_new > 0 and phi_new >= phi - 1e-9:
                    accepted = True
                    break
            beta *= 0.5
        if accepted:
            xi, V, phi = xi_new, V_new, phi_new
        # ---- FW atom injection on stall ----------------------------------- #
        if gap > last_gap * 0.999:
            stall += 1
        else:
            stall = 0
        last_gap = gap
        if stall >= 5 and gap >= tol:
            astar = np.unravel_index(int(np.argmax(d)), d.shape)
            gamma = min(0.05, gap / (gap + 1.0))
            for _ in range(6):
                xi_try = (1.0 - gamma) * xi
                xi_try[astar] += gamma
                V_try, eps = assemble(xi_try)
                s_t, phi_t = np.linalg.slogdet(V_try)
                if s_t > 0 and phi_t >= phi - 1e-9:
                    xi, V, phi = xi_try, V_try, phi_t
                    stall = 0
                    break
                gamma *= 0.5
        # ---- prune decayed atoms ------------------------------------------ #
        thr = xi.max() * 1e-16
        xi[xi < thr] = 0.0
        xi /= xi.sum()

    final_gap = float(gap)

    # ---- rounding: largest-remainder apportionment (deterministic) -------- #
    flat = xi.ravel()
    thr = flat.max() * 1e-6
    supp = np.where(flat > thr)[0]
    w_s = flat[supp] / flat[supp].sum()
    quota = w_s * M_rows
    base = np.floor(quota).astype(np.int64)
    rem = int(M_rows - base.sum())
    frac = quota - base
    order = np.lexsort((supp, -frac))
    for i in range(rem):
        base[order[i]] += 1

    Aout = np.zeros((M_rows, n))
    atoms = []
    row = 0
    for j, fa in enumerate(supp):
        c = int(base[j])
        if c == 0:
            continue
        geo, l = divmod(int(fa), L)
        a_row = np.zeros(n)
        a_row[IDX[geo]] = rho[l] * GVAL[geo]        # load = rho[l] exactly
        Aout[row:row + c] = a_row
        atoms.append((names[geo], int(t_arr[geo]), float(p_arr[geo]),
                      float(rho[l]), c))
        row += c
    assert row == M_rows, (row, M_rows)

    # ---- G5: per-pixel cumulative dose (post-rounding report, no fix) ----- #
    dose = Aout.mean(axis=0)
    rel = dose / dose.mean()
    g5 = {"max_rel_dev": float(np.abs(rel - 1.0).max()),
          "band": 0.05,
          "pass": bool(np.abs(rel - 1.0).max() <= 0.05)}

    cert = {"n": n, "final_gap": final_gap, "dict_size": n_atoms,
            "n_geometries": G, "iters": it, "eps": float(eps),
            "logdet": float(phi), "converged": bool(final_gap < tol),
            "J_source": Jfun.source, "nu": float(nu),
            "rho_levels": [float(r) for r in rho],
            "nuJ_levels": [float(x) for x in nuJ],
            "J_arm_prescan": float(J_arm), "arm_rho": float(arm_rho)}
    guards = {"G3": g3, "G4": g4, "G5": g5}

    return {"A": Aout, "atoms": atoms, "gap": final_gap,
            "gap_traj": gap_traj, "cert": cert, "guards": guards,
            "V": V, "V0": V0, "xi": xi,
            "_geom": (IDX, GVAL, names, t_arr, p_arr, ALLOW, nuJ)}


# ========================================================================== #
#  Baseline helper: uniform single-shape design through the SAME v2 atoms     #
# ========================================================================== #
def baseline_V(xhat, shape_name, out, M_rows=1024, p_sel=0.0, loads=None):
    """V = V0 + M_rows * mean H over one shape's admissible p = p_sel atoms
    (uniform over translations x admissible loads), same V0/atoms as `out`."""
    IDX, GVAL, names, t_arr, p_arr, ALLOW, nuJ = out["_geom"]
    names = np.asarray(names)
    sel = (names == shape_name) & (p_arr == p_sel)
    xi = np.zeros(ALLOW.shape)
    mask = ALLOW & sel[:, None]
    if loads is not None:
        lm = np.isin(np.arange(ALLOW.shape[1]),
                     [list(out["cert"]["rho_levels"]).index(x) for x in loads])
        mask = mask & lm[None, :]
    xi[mask] = 1.0
    xi /= xi.sum()
    n = xhat.size
    u_geo = M_rows * (xi * nuJ[None, :]).sum(axis=1)
    M = _build_M(IDX, GVAL, u_geo, n)
    Vb = out["V0"] + M
    epsb = 1e-9 * np.trace(Vb) / n
    Vb[np.diag_indices(n)] += epsb
    return Vb, int(mask.sum())


# ========================================================================== #
#  Selftest                                                                   #
# ========================================================================== #
if __name__ == "__main__":
    import csv
    import time
    from collections import Counter

    t_all0 = time.time()
    SIDE, N = 32, 1024
    NU = 2000.0
    RHOS = (0.3, 0.6, 1.0)

    # ----- scene (matches dev_pilot_L3v2 / v1 selftest) -------------------- #
    yy, xx = np.mgrid[0:SIDE, 0:SIDE]
    scene = np.ones((SIDE, SIDE))
    for (cy, cx) in [(8, 9), (10, 23), (22, 7), (24, 21)]:
        scene += 0.6 * np.exp(-(((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * 3.5 ** 2)))
    x = (scene / scene.sum()).ravel()
    coarse = x.reshape(SIDE, SIDE).reshape(SIDE // 4, 4, SIDE // 4, 4).mean(axis=(1, 3))
    xhat = np.repeat(np.repeat(coarse, 4, axis=0), 4, axis=1).ravel()
    print("[v2 selftest] xhat: n=%d sum=%.6f cv=%.3f" %
          (xhat.size, xhat.sum(), xhat.std() / xhat.mean()), flush=True)

    # ----- J source + cross-checks ----------------------------------------- #
    Jfun = get_J_source()
    print("\n[v2 selftest] === J source: %s ===" % Jfun.source, flush=True)
    ok_j = True
    if os.path.exists(J_CSV):
        rows = []
        with open(J_CSV, newline="") as f:
            rd = csv.reader(f)
            next(rd)
            for r in rd:
                if float(r[0]) == NU:
                    rows.append((float(r[1]), float(r[2])))
        worst = 0.0
        for target in RHOS:
            rr, jj = min(rows, key=lambda t: abs(math.log(t[0] / target)))
            jt = float(Jfun(rr, NU))
            rel = abs(jt - jj) / jj
            worst = max(worst, rel)
            print("  interp vs fisher_ridge.csv @ nu=%d rho=%.4f: "
                  "csv=%.6e interp=%.6e rel=%.2e" % (NU, rr, jj, jt, rel))
        ok_j = worst < 1e-3
        print("  table-interp cross-check (rel < 1e-3)      -> %s"
              % ("PASS" if ok_j else "FAIL"))
    Jnum = JExactNumpy()
    worst_fb = 0.0
    for r in RHOS:
        jt, jn = float(Jfun(r, NU)), float(Jnum(r, NU))
        rel = abs(jt - jn) / jt
        worst_fb = max(worst_fb, rel)
        print("  fallback exact vs table @ rho=%.1f: table=%.6e "
              "numpy=%.6e rel=%.2e" % (r, jt, jn, rel))
    ok_fb = worst_fb < 1e-2
    print("  numpy-fallback cross-check (rel < 1e-2)    -> %s"
          % ("PASS" if ok_fb else "FAIL"))
    print("  J monotone on [0.3,1.0] at nu=2000: J=%.4f < %.4f < %.4f"
          % (Jfun(0.3, NU), Jfun(0.6, NU), Jfun(1.0, NU)))

    # ----- run the v2 design ----------------------------------------------- #
    print("\n[v2 selftest] === design_v2(xhat, nu=%d, rho_levels=%s) ==="
          % (NU, list(RHOS)), flush=True)
    t0 = time.time()
    out = design_v2(xhat, rho_levels=RHOS, nu=NU, J_source=Jfun, verbose=True)
    t_design = time.time() - t0
    cert, guards = out["cert"], out["guards"]
    print("  final KW gap (V0-offset) = %.6f  (tol 0.1, %s)"
          % (out["gap"], "REACHED" if cert["converged"] else "iter-capped"))
    print("  iterations = %d   dictionary = %d admissible atoms "
          "(%d geometries x %d loads - G3 drops)"
          % (cert["iters"], cert["dict_size"], cert["n_geometries"],
             len(RHOS)))
    print("  nu*J per load: " + ", ".join(
        "rho=%.1f -> %.2f" % (r, c)
        for r, c in zip(cert["rho_levels"], cert["nuJ_levels"])))
    ok_gap = out["gap"] < 0.1

    # ----- (i) guard reports ----------------------------------------------- #
    print("\n[v2 selftest] === (i) guard reports ===", flush=True)
    g3, g4, g5 = guards["G3"], guards["G4"], guards["G5"]
    print("  G3 peak<=%.0f : dropped per load %s  (of %d geometries)"
          % (g3["peak_cap"],
             {("%.1f" % k): v for k, v in g3["dropped_per_load"].items()},
             g3["n_geometries"]))
    peak_final = out["A"].max()
    ok_g3 = peak_final <= 64.0 + 1e-9
    print("  G3 final design peak a_j = %.3f (<= 64)      -> %s"
          % (peak_final, "PASS" if ok_g3 else "FAIL"))
    print("  G4 clip [1/4,4]: clipped weight fraction = %.2e; "
          "cv-guard drops = %d; powers = {0,1}"
          % (g4["clipped_weight_fraction"], g4["cv_guard_dropped"]))
    print("  G5 per-pixel dose +/-5%%: max rel dev = %.3f    -> %s"
          % (g5["max_rel_dev"], "PASS" if g5["pass"] else "FAIL (reported, not fixed)"))

    # row-load check: each row's a.xhat must equal its atom's rho_level
    loads = out["A"] @ xhat
    expect = np.concatenate([[a[3]] * a[4] for a in out["atoms"]])
    rel_load = np.abs(loads - expect) / expect
    ok_load = rel_load.max() <= 1e-9
    print("  row load a.xhat == rho_level (rel<=1e-9): max rel = %.2e -> %s"
          % (rel_load.max(), "PASS" if ok_load else "FAIL"))

    # ----- (ii) certificate ------------------------------------------------ #
    print("\n[v2 selftest] === (ii) KW certificate gap < 0.1 ===")
    print("  gap = %.6f                                  -> %s"
          % (out["gap"], "PASS" if ok_gap else "FAIL"))

    # ----- (iii) log-det vs single-shape baselines (same v2 atoms) --------- #
    print("\n[v2 selftest] === (iii) log det(V0+M) vs baselines ===",
          flush=True)
    ld_opt = cert["logdet"]
    Vb_s, ns_ = baseline_V(xhat, "scatter16", out)
    Vb_l, nl_ = baseline_V(xhat, "Lblob6x6", out)
    _, ld_s = np.linalg.slogdet(Vb_s)
    _, ld_l = np.linalg.slogdet(Vb_l)
    print("  optimized v2 design            : %.4f" % ld_opt)
    print("  pure scatter16 (p=0, %4d atoms): %.4f" % (ns_, ld_s))
    print("  pure Lblob6x6  (p=0, %4d atoms): %.4f" % (nl_, ld_l))
    ok_ld = (ld_opt >= ld_s - 1e-6) and (ld_opt >= ld_l - 1e-6)
    print("  optimized >= both baselines                 -> %s"
          % ("PASS" if ok_ld else "FAIL"))

    # ----- G6 hook vs the scattered fixed baseline ------------------------- #
    g6 = g6_spectral_floor(out["V"], Vb_s, r=200, floor=0.5)
    print("\n[v2 selftest] G6 hook vs scatter16 baseline (top-r=%d): "
          "mu_min = %.4f -> %s%s"
          % (g6["r"], g6["mu_min"], "PASS" if g6["pass"] else "FAIL",
             "" if g6["pass"] else
             ("; alpha-mixture needed: alpha <= %.4f" % g6["alpha_max"])))

    # ----- (iv) load histogram --------------------------------------------- #
    print("\n[v2 selftest] === (iv) selected-atom composition ===")
    n_rows_by = Counter()
    n_rows_load = Counter()
    for (nm, t, p, r, c) in out["atoms"]:
        n_rows_by[(nm, p, r)] += c
        n_rows_load[r] += c
    for key in sorted(n_rows_by):
        print("  %-12s p=%.1f rho=%.1f : %4d rows"
              % (key[0], key[1], key[2], n_rows_by[key]))
    print("  load histogram: " + ", ".join(
        "rho=%.1f -> %d rows" % (r, n_rows_load.get(r, 0)) for r in RHOS))
    shares = np.array([n_rows_load.get(r, 0) for r in RHOS]) / 1024.0
    if shares.max() > 0.98:
        print("  WARNING: optimizer did not mix load levels "
              "(%.0f%% of rows at rho=%.1f)."
              % (100 * shares.max(), RHOS[int(shares.argmax())]))
        print("  DIAGNOSIS: the chain-rule direction g = a/(a.xhat) is "
              "invariant to the row's load scaling, so the same geometry at "
              "two loads differs ONLY by the scalar nu*J(rho,nu); J is "
              "strictly increasing on [0.3,1.0] at nu=2000 (ridge rho* ~ "
              "%.1f), so the unconstrained optimizer always prefers the "
              "highest load that survives G3. Load mixing requires either "
              "G3 to drop rho=1.0 atoms on the needed dim-support geometries "
              "or the ruling's G2 design-average-load equality constraint "
              "(design-average exactly 0.6), which is NOT implemented in v2 "
              "(amendment hook)." % ((6 * NU) ** (1 / 3) - 2 / 3))
    else:
        print("  load levels mixed (max share %.0f%%)." % (100 * shares.max()))

    # ----- (v) wall time --------------------------------------------------- #
    wall = time.time() - t_all0
    ok_time = wall < 600.0
    print("\n[v2 selftest] === (v) wall time ===")
    print("  design_v2 = %.1f s; selftest total = %.1f s (< 600 s) -> %s"
          % (t_design, wall, "PASS" if ok_time else "FAIL"))

    hard = {"J-cross-check": ok_j, "J-fallback": ok_fb, "gap<0.1": ok_gap,
            "G3-final-peak": ok_g3, "row-loads": ok_load,
            "logdet-beats-baselines": ok_ld, "wall<10min": ok_time}
    all_ok = all(hard.values())
    print("\n[v2 selftest] hard checks: %s" % hard)
    print("[v2 selftest] G5 dose report: %s (report-only, per spec)"
          % ("PASS" if g5["pass"] else "FAIL max dev %.3f" % g5["max_rel_dev"]))
    print("[v2 selftest] RESULT: %s"
          % ("ALL HARD CHECKS PASS" if all_ok else "*** FAILURE ***"))
