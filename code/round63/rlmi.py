"""RLMI -- Robust Library Mixture Illumination allocator + materializer.

T-C of the ROUND63 DEV bridge.  NORMATIVE spec:
`docs/ROUND63_GPT_ROUND25_RULING_RAW.md` (issue #17); summary in
`docs/ROUND63_BRIDGE_BUILD_PLAN.md` T-C.  This module implements the R25 Step-3
eight-step pipeline (R25 SS13):

  1. produce S=16 declared parametric-bootstrap scenario matrices {x_s};
  2. compute each scenario's bank-simplex A-risk oracle  R*_s;
  3. solve the standardized maximin relative-regret program (R25.1);
  4. lexicographic closest-to-knob (e_0) tie-break;
  5. exact 972-row materialization from the union design measure by constrained
     rounding + deterministic exchange (R25 SS9);
  6. recompute guards, realized regret, and the simplex KKT certificate (R25.2);
  7. L0 fallback ONLY on solver / materialization / guard failure, w/ a flag;
  8. (acquire once + RQL -- outside this module; T-D).

The method is a robust library MIXTURE, not an abstaining router.  No runtime
0.01/1.02/2%/hysteresis thresholds exist (R25 SS13).

--------------------------------------------------------------------------------
DESIGN NOTES / interpretation register (quote + interpretation of every R25 /
falsifier ambiguity; the build plan's I-register discipline).
--------------------------------------------------------------------------------

KERNEL (on-branch evaluation; falsifier addendum docs2/14).
  docs2/14 requires every diagnostic to be evaluated ON-BRANCH and offers a
  choice: "either clip the evaluation to the declared branch OR use the exact
  finite-nu kernel with its true (non-concave) shape in the F_k prediction while
  keeping the CERTIFICATE claims branch-fenced."  We adopt the R23 closed-form
  finite-c dead-time information kernel (R23 SS1, eq. for J):

        J(rho; c) = rho / [ (1+rho) (1 + c^2 rho^2) ] ,

  evaluated per frame with an explicit nu-frame multiplier (V += nu*J(rho)*Q).
  Interpretation & justification for choosing this over `physics.fisher_exact`:
    (i)  it is the exact kernel the R23/R24 branch + residual algebra is written
         in -- it supplies the *closed-form* J', J'' and hence a machine-exact
         alpha*(rho), the strong-concavity modulus, and the argmax rho_c that the
         branch-fenced certificate needs;
    (ii) it carries the hidden hold-time CV c, required by the hard decision
         cell (c=0.05); `physics.fisher_exact` is the c=0 renewal kernel only;
    (iii) it reproduces the c->0, moderate-rho renewal information to O(rho^2),
         so the F_k *prediction* uses the true (for rho>rho_c, DECREASING) shape
         and never fabricates information beyond the ridge.
  BRANCH: on 0<=rho<=rho_c, J'>=0 and J''<0 (monotone-concave); rho_c is the
  positive root of c^2 rho^2 (1+2 rho) = 1 (rho_c=+inf for c=0).  A bank whose
  operating loads cross rho_c violates the R23 branch hypothesis: we emit
  BRANCH_VIOLATION and report the eps_1^2/(2m) bound UNAVAILABLE (never
  vacuous), exactly as docs2/14 SS2 requires.  The F_k / A-risk assembly still
  uses J's true value at those loads.

H0 (R24 SS2.3 step 2):  H0(x) = F_pre(x) + lambda_TV R_delta(x) + eps I.
  F_pre  = balanced-52 pre-scan information at operating loads (reproduces
           oed_design_v4.balanced_prescan_52 / V0_prescan; documented below).
  R_delta= frozen differentiable TV/IRLS (lagged-diffusivity) curvature of
           isotropic 2-D TV:  R_delta = Dx^T diag(w) Dx + Dy^T diag(w) Dy with
           per-pixel isotropic weight w_p = 1/sqrt(gx_p^2+gy_p^2+delta^2)
           (reproduces solvers._grad/_div; delta = tv_smooth).
  All projected onto the declared reconstruction subspace B (r-dim); W=I_r.

SCENARIOS (R25 SS13.1 / R24 step 3): S frozen parametric-bootstrap draws
  x_s = clip_[0,1]( xhat + H0(xhat)^{-1/2} z_s ),  z_s ~ N(0, I) with a frozen
  per-scenario seed.  H0(xhat)^{-1} is the RQL local (Laplace) covariance -- the
  inverse of the pre-scan+prior curvature.  Box projection keeps draws physical
  (documented approximation; the Laplace law is the local Gaussian, the [0,1]
  clip is a projection onto the physical scene box).

OBJECTIVE (R25.1): min_{w in simplex} max_s R_s(w)/R*_s, R_s(w)=tr[W X_s(w)^-1],
  X_s(w)=H0(x_s)+sum_k w_k F_{k,s}, F_{k,s}=972 sum_r p_{k,r} H_s(a_{k,r}).
  t_cont = the achieved worst relative regret; t_cont>=1.

CERTIFICATE (R25.2): least-favorable scenario multipliers alpha_s>=0, sum=1,
  supported on active scenarios (r_s=t); with the positive marginal value
  g_{s,k}=tr[W X_s^-1 F_{k,s} X_s^-1]/R*_s,  sum_s alpha_s g_{s,k} = lambda on
  supp(w), <= lambda off support; complementarity alpha_s (r_s - t)=0.

MATERIALIZATION (R25 SS9): merge duplicate atoms (weight addition), target
  p_a = 972 xi(w)(a), largest-remainder rounding to total=972, then a
  DETERMINISTIC EXCHANGE enforcing incident budget, +-5% dose band, and
  per-row peak/admission.  No online full-dictionary OED.  Final block order
  randomized with a frozen seed.

Only numpy + scipy are used.  This module is self-contained by design (the
frozen campaign modules pull in heavy deps via patterns.gi_core); the pre-scan
and exchange constructions are re-implemented here and documented as
reproductions of oed_design_v4's interfaces, per the T-C "else implement a clean
deterministic exchange with documented move rules" allowance.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field, asdict
from typing import Optional

import numpy as np
from scipy.linalg import cho_factor, cho_solve, cholesky, solve_triangular, eigh
from scipy.optimize import brentq, nnls

# ============================================================================ #
#  1.  finite-c dead-time information kernel (R23 SS1) -- exact closed form     #
# ============================================================================ #

def kernel_J(rho, c):
    """J(rho;c) = rho / [(1+rho)(1+c^2 rho^2)]  (per-frame information rate)."""
    rho = np.asarray(rho, dtype=np.float64)
    return rho / ((1.0 + rho) * (1.0 + c * c * rho * rho))


def kernel_Jp(rho, c):
    """J'(rho;c) = [1 - c^2 rho^2 (1+2 rho)] / [(1+rho)^2 (1+c^2 rho^2)^2]."""
    rho = np.asarray(rho, dtype=np.float64)
    c2 = c * c
    num = 1.0 - c2 * rho * rho * (1.0 + 2.0 * rho)
    den = (1.0 + rho) ** 2 * (1.0 + c2 * rho * rho) ** 2
    return num / den


def kernel_Jpp(rho, c):
    """J''(rho;c), exact quotient derivative of J' (see module header).
    N = 1 - c^2 rho^2 - 2 c^2 rho^3 ;  D = (1+rho)^2 (1+c^2 rho^2)^2 ;
    J'' = (N' D - N D') / D^2 with
    N'  = -2 c^2 rho - 6 c^2 rho^2 ,
    D'  = 2(1+rho)(1+c^2 rho^2)(1 + 2 c^2 rho + 3 c^2 rho^2)."""
    rho = np.asarray(rho, dtype=np.float64)
    c2 = c * c
    N = 1.0 - c2 * rho * rho - 2.0 * c2 * rho ** 3
    Np = -2.0 * c2 * rho - 6.0 * c2 * rho * rho
    one_c = 1.0 + c2 * rho * rho
    D = (1.0 + rho) ** 2 * one_c ** 2
    Dp = 2.0 * (1.0 + rho) * one_c * (1.0 + 2.0 * c2 * rho + 3.0 * c2 * rho * rho)
    return (Np * D - N * Dp) / (D * D)


def kernel_rho_c(c):
    """Argmax of J (root of c^2 rho^2 (1+2 rho) = 1); +inf for c<=0."""
    if c <= 0:
        return np.inf
    f = lambda r: c * c * r * r * (1.0 + 2.0 * r) - 1.0
    hi = 1.0
    while f(hi) < 0:
        hi *= 2.0
        if hi > 1e12:
            return np.inf
    return float(brentq(f, 0.0, hi, xtol=1e-12, rtol=1e-14))


# ============================================================================ #
#  2.  pre-scan + TV curvature (reproduce oed_design_v4 / solvers interfaces)   #
# ============================================================================ #

def balanced_prescan_52(side=32):
    """52-row balanced pre-scan (reproduces oed_design_v4.balanced_prescan_52):
    32 paired-4x4-block rows (amp n/32), 16 8x8-block rows (amp n/64), 4 quadrant
    rows (amp n/256).  Mean row == all-ones and per-pixel dose exactly equal."""
    n = side * side
    P = np.zeros((52, n))
    r = 0
    for q in range(32):
        for blk in (q, q + 32):
            by, bx = divmod(blk, 8)
            ys, xs = 4 * by, 4 * bx
            for dy in range(4):
                for dx in range(4):
                    P[r, (ys + dy) * side + (xs + dx)] = n / 32.0
        r += 1
    for blk in range(16):
        by, bx = divmod(blk, 4)
        ys, xs = 8 * by, 8 * bx
        for dy in range(8):
            for dx in range(8):
                P[r, (ys + dy) * side + (xs + dx)] = n / 64.0
        r += 1
    for blk in range(4):
        by, bx = divmod(blk, 2)
        ys, xs = 16 * by, 16 * bx
        for dy in range(16):
            for dx in range(16):
                P[r, (ys + dy) * side + (xs + dx)] = n / 256.0
        r += 1
    assert r == 52
    return P


def default_prescan(side):
    """Balanced multiscale-block pre-scan for arbitrary side (power of 2).  For
    side==32 this returns the exact 52-row campaign construction
    (balanced_prescan_52); otherwise a generic multiscale block set (positive
    loads, full coverage) used by the unit tests."""
    if side == 32:
        return balanced_prescan_52(32)
    n = side * side
    rows = []
    scales = sorted({max(1, side // 4), max(1, side // 2), side})
    for bs in scales:
        nb = side // bs
        for by in range(nb):
            for bx in range(nb):
                row = np.zeros(n)
                for dy in range(bs):
                    for dx in range(bs):
                        row[(by * bs + dy) * side + (bx * bs + dx)] = n / (bs * bs)
                rows.append(row)
    return np.array(rows)


def _grad_ops(side):
    """Dense forward-difference operators (Dx, Dy) matching solvers._grad
    (Neumann/zero boundary).  Returns (Dx, Dy) each (n, n)."""
    n = side * side
    Dx = np.zeros((n, n))
    Dy = np.zeros((n, n))
    for i in range(side):
        for j in range(side):
            p = i * side + j
            if i + 1 < side:                       # gx[i,j] = x[i+1,j]-x[i,j]
                Dx[p, (i + 1) * side + j] += 1.0
                Dx[p, p] -= 1.0
            if j + 1 < side:                       # gy[i,j] = x[i,j+1]-x[i,j]
                Dy[p, i * side + (j + 1)] += 1.0
                Dy[p, p] -= 1.0
    return Dx, Dy


def tv_curvature(x, side, Dx, Dy, tv_smooth=1e-3):
    """Frozen differentiable TV/IRLS (lagged-diffusivity) curvature R_delta(x):
    R_delta = Dx^T diag(w) Dx + Dy^T diag(w) Dy, isotropic per-pixel reweight
    w_p = 1/sqrt(gx_p^2 + gy_p^2 + delta^2).  (Hessian of the Charbonnier-smoothed
    isotropic TV, dropping the always-PSD rank-1 term -- the standard IRLS design
    curvature used only for local scoring.)"""
    gx = Dx @ x
    gy = Dy @ x
    w = 1.0 / np.sqrt(gx * gx + gy * gy + tv_smooth * tv_smooth)
    return (Dx.T * w) @ Dx + (Dy.T * w) @ Dy


# ============================================================================ #
#  3.  bank schema (mixable approximate-design measure xi_k; R25 SS8)           #
# ============================================================================ #

@dataclass
class Bank:
    """Mixable approximate-design measure xi_k = {(a_{k,r}, p_{k,r})}.

    rows    : (R, n) nonneg row atom pool (relative-irradiance load-profile rows).
    weights : (R,)   normalized p_{k,r} (sum == 1).
    load    : (R,)   predicted operating load a.xhat of each atom (meta).
    incident: (R,)   incident cost of each atom (meta; default row-sum).
    peak    : (R,)   per-atom peak (meta; default row max).
    admission:(R,)   bool admission mask (meta; default True).
    name    : bank id (L0..L7).  is_knob: True for L0 (the fallback corner e_0)."""
    name: str
    rows: np.ndarray
    weights: np.ndarray
    load: np.ndarray
    incident: np.ndarray
    peak: np.ndarray
    admission: np.ndarray
    is_knob: bool = False
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        self.rows = np.asarray(self.rows, dtype=np.float64)
        w = np.asarray(self.weights, dtype=np.float64)
        s = w.sum()
        self.weights = w / s if s > 0 else w
        R = self.rows.shape[0]
        if self.load is None:
            self.load = self.rows.sum(axis=1)
        if self.incident is None:
            self.incident = self.rows.sum(axis=1)
        if self.peak is None:
            self.peak = self.rows.max(axis=1)
        if self.admission is None:
            self.admission = np.ones(R, dtype=bool)
        self.load = np.asarray(self.load, dtype=np.float64)
        self.incident = np.asarray(self.incident, dtype=np.float64)
        self.peak = np.asarray(self.peak, dtype=np.float64)
        self.admission = np.asarray(self.admission, dtype=bool)


def load_bank(path, name=None, is_knob=False):
    """Load a T-B bank npz.  Per-row metadata is read from the npz arrays
    (atom_load/atom_incident/atom_peak/atom_admission); the R27 §2 exact-972
    admission witness (admission_witness_mult, aligned to rows) is stored in
    meta['witness_mult'].  meta['admission'] is the R27 admission OBJECT (dict),
    NOT a per-row mask -- do not coerce it to bool."""
    d = np.load(path, allow_pickle=True)
    files = set(d.files)
    rows = np.asarray(d["rows"], dtype=np.float64)
    weights = np.asarray(d["weights"], dtype=np.float64)
    meta = {}
    if "meta" in files:
        m = d["meta"]
        meta = json.loads(m.item() if hasattr(m, "item") else str(m))

    def npz(k, dt=np.float64):
        return np.asarray(d[k], dtype=dt) if k in files else None

    admission = None
    if "atom_admission" in files:
        admission = np.asarray(d["atom_admission"], dtype=bool)
    wit = npz("admission_witness_mult", np.int64)
    if wit is not None:
        meta = dict(meta)
        meta["witness_mult"] = wit                     # (R,) int, aligned to rows
    return Bank(name=name or meta.get("name", "L?"), rows=rows, weights=weights,
                load=npz("atom_load"), incident=npz("atom_incident"),
                peak=npz("atom_peak"), admission=admission,
                is_knob=is_knob, meta=meta)


# ============================================================================ #
#  4.  scenario generation (frozen Laplace / parametric-bootstrap draws)        #
# ============================================================================ #

@dataclass
class ScenarioSet:
    xhat: np.ndarray
    draws: np.ndarray            # (S, n)
    seeds: list
    cov_source: str


def build_H0_full(xhat, side, nu, c, lambda_TV, eps, rho_prescan=0.60,
                  tv_smooth=1e-3, prescan_rows=None, Dx=None, Dy=None):
    """Full-space RQL-TV local curvature H0(xhat) = F_pre(xhat) + lambda_TV
    R_delta(xhat) + eps I  (n x n), the Laplace precision used by make_scenarios
    (R24 SS2.3 step 2, unprojected).  Convenience for T-D so the scenario
    covariance is built with the SAME pre-scan + TV curvature as the allocator."""
    xhat = np.asarray(xhat, dtype=np.float64)
    n = xhat.size
    if prescan_rows is None:
        prescan_rows = default_prescan(side)
    if Dx is None or Dy is None:
        Dx, Dy = _grad_ops(side)
    pl = prescan_rows @ xhat
    rho_bar = rho_prescan / max(float(pl.mean()), 1e-12)
    Hpre = np.zeros((n, n))
    for m in range(prescan_rows.shape[0]):
        lm = float(pl[m])
        if lm <= 1e-12:
            continue
        q = prescan_rows[m] / lm
        Hpre += (nu * float(kernel_J(rho_bar * lm, c))) * np.outer(q, q)
    Rd = tv_curvature(xhat, side, Dx, Dy, tv_smooth)
    H = Hpre + lambda_TV * Rd + eps * np.eye(n)
    return 0.5 * (H + H.T)


def make_scenarios(xhat, H0_full, n_scenarios=16, seed_base=650000,
                   clip=(0.0, 1.0), jitter_scale=1.0):
    """S frozen parametric-bootstrap draws x_s = clip(xhat + L z_s), L L^T =
    (jitter_scale^2) H0_full^{-1} (RQL local / Laplace covariance from the
    pre-scan+prior curvature).  Per-scenario frozen seed = seed_base + s.

    Uses the eigendecomposition of H0_full (SPD): H0^{-1/2} = U diag(ev^-1/2) U^T.
    The seed base 650000 is the bridge cohort base (disjoint from 632900 DEV /
    633000 confirmatory), documented in the build plan interface contract."""
    xhat = np.asarray(xhat, dtype=np.float64)
    ev, U = np.linalg.eigh(0.5 * (H0_full + H0_full.T))
    ev = np.maximum(ev, 1e-12)
    half_inv = U * (1.0 / np.sqrt(ev))            # columns scaled -> H0^{-1/2}=half_inv@U.T
    draws = np.empty((n_scenarios, xhat.size))
    seeds = []
    for s in range(n_scenarios):
        seed = int(seed_base + s)
        seeds.append(seed)
        z = np.random.default_rng(seed).standard_normal(xhat.size)
        pert = jitter_scale * (half_inv @ (U.T @ z))
        draws[s] = np.clip(xhat + pert, clip[0], clip[1])
    return ScenarioSet(xhat=xhat, draws=draws, seeds=seeds,
                       cov_source="laplace_H0inv")


# ============================================================================ #
#  5.  per-scenario information matrices (projected onto subspace B; R25 SS8)    #
# ============================================================================ #

@dataclass
class ScenarioMatrices:
    H0: list                     # [S] r x r  local RQL-TV curvature (projected)
    F: list                      # [S][K] r x r per-bank 972-row information
    Rstar: np.ndarray            # [S] bank-simplex A-risk oracle values
    vstar: list                  # [S] oracle simplex weights
    W: np.ndarray                # r x r task weight (I_r)
    branch: dict                 # per-bank branch-violation flags/loads
    kernel_c: float
    nu: float
    r: int


def _proj_row_info(row, x, B, nu, c, ridge_floor=1e-12):
    """nu*J(rho)*(B^T q)(B^T q)^T for one physical row atom at scene x, on-branch:
    rho = a.x (operating load), q = a/rho.  Returns (r x r, rho, in_branch?)."""
    load = float(row @ x)
    if load <= ridge_floor:
        r = B.shape[1]
        return np.zeros((r, r)), 0.0
    q = row / load
    g = B.T @ q
    return (nu * float(kernel_J(load, c))) * np.outer(g, g), load


def build_scenario_matrices(scenarios, banks, B, nu, c, lambda_TV, eps,
                            side, prescan_rows=None, rho_prescan=0.60,
                            tv_smooth=1e-3, W=None, M_rows=972,
                            Dx=None, Dy=None, freeze_directions=True):
    """Assemble {H0_s, F_{k,s}} projected onto B, plus per-bank branch flags.

    F_{k,s} = M_rows * sum_r p_{k,r} H_s(a_{k,r})  (R25 SS8).
    H0_s    = B^T[F_pre(x_s) + lambda_TV R_delta(x_s) + eps I]B  (R24 SS2.3 step2).
    On-branch: every kernel is evaluated at the atom's own operating load a.x_s;
    if a bank's loads cross rho_c a BRANCH_VIOLATION is recorded (bound
    UNAVAILABLE downstream) but J keeps its true value in F.

    freeze_directions: local-linearization (freeze q = a/(a.xhat) across the S
    small Laplace draws so B^T q is precomputable; only the scalar kernel weight
    J(a.x_s) varies per scenario).  This is the SS13 latency precompute; set
    False for the exact per-scene direction (used in the certificate recompute)."""
    xhat = scenarios.xhat
    draws = scenarios.draws
    S = draws.shape[0]
    K = len(banks)
    r = B.shape[1]
    if W is None:
        W = np.eye(r)
    if prescan_rows is None:
        prescan_rows = default_prescan(side)
    if Dx is None or Dy is None:
        Dx, Dy = _grad_ops(side)
    rho_c = kernel_rho_c(c)
    I_r = np.eye(r)

    # pre-scan flux scale so mean pre-scan operating load == rho_prescan (at xhat)
    pl = prescan_rows @ xhat
    rho_bar = rho_prescan / max(float(pl.mean()), 1e-12)

    # frozen projected atom directions g = B^T (a/(a.xhat)) per bank
    frozen_g, frozen_load0 = [], []
    for bk in banks:
        loads0 = bk.rows @ xhat
        loads0 = np.maximum(loads0, 1e-12)
        G = (bk.rows / loads0[:, None]) @ B          # (R, r) projected directions
        frozen_g.append(G)
        frozen_load0.append(loads0)

    H0_list, F_list = [], []
    branch = {bk.name: {"violation": False, "n_beyond": 0, "frac_beyond": 0.0,
                        "rho_c": float(rho_c), "max_load": 0.0} for bk in banks}
    for s in range(S):
        x = draws[s]
        # F_pre(x_s): balanced-52 pre-scan at operating loads (V0_prescan form)
        Hpre = np.zeros((r, r))
        p_loads = prescan_rows @ x
        for m in range(prescan_rows.shape[0]):
            lm = float(p_loads[m])
            if lm <= 1e-12:
                continue
            rho_m = rho_bar * lm
            g = B.T @ (prescan_rows[m] / lm)
            Hpre += (nu * float(kernel_J(rho_m, c))) * np.outer(g, g)
        # R_delta(x_s) projected
        Rd = B.T @ tv_curvature(x, side, Dx, Dy, tv_smooth) @ B
        H0 = Hpre + lambda_TV * Rd + eps * I_r
        H0_list.append(0.5 * (H0 + H0.T))

        Fk = []
        for kk, bk in enumerate(banks):
            if freeze_directions:
                loads = bk.rows @ x
                loads = np.maximum(loads, 1e-12)
                Jw = nu * kernel_J(loads, c)                     # (R,)
                wgt = M_rows * bk.weights * Jw                   # per-atom weight
                G = frozen_g[kk]
                Fks = (G.T * wgt) @ G                            # r x r
            else:
                Fks = np.zeros((r, r))
                loads = bk.rows @ x
                for rr in range(bk.rows.shape[0]):
                    Hr, _ = _proj_row_info(bk.rows[rr], x, B, nu, c)
                    Fks += M_rows * bk.weights[rr] * Hr
            Fk.append(0.5 * (Fks + Fks.T))
            # branch bookkeeping (max over scenarios)
            nbeyond = int(np.sum(loads > rho_c))
            if nbeyond > 0:
                branch[bk.name]["violation"] = True
            branch[bk.name]["n_beyond"] = max(branch[bk.name]["n_beyond"], nbeyond)
            branch[bk.name]["frac_beyond"] = max(
                branch[bk.name]["frac_beyond"], nbeyond / bk.rows.shape[0])
            branch[bk.name]["max_load"] = max(branch[bk.name]["max_load"],
                                              float(loads.max()))
        F_list.append(Fk)

    sm = ScenarioMatrices(H0=H0_list, F=F_list, Rstar=None, vstar=None, W=W,
                          branch=branch, kernel_c=c, nu=nu, r=r)
    return sm


# ============================================================================ #
#  6.  linear-algebra primitives + bank-simplex oracle (R25.1, R25 SS4)         #
# ============================================================================ #

def _spd_inv(X):
    """Symmetric-PD inverse via Cholesky (raises on non-PD)."""
    cf = cho_factor(0.5 * (X + X.T), lower=True, check_finite=False)
    return cho_solve(cf, np.eye(X.shape[0]), check_finite=False)


def arisk(X, W):
    """R = tr[W X^{-1}]  (A-risk / estimator risk on the declared subspace)."""
    Xi = _spd_inv(X)
    return float(np.einsum("ij,ji->", W, Xi)), Xi


def _arisk_and_grads(X, Fk, W):
    """R = tr[W X^-1] and the K marginal gradients dR/dv_k = -tr[W X^-1 F_k X^-1].

    Fast: precompute XiWXi = X^-1 W X^-1 once (2 r^3), then each gradient is the
    Frobenius inner product <F_k, XiWXi> (O(r^2)), using
    tr[W X^-1 F_k X^-1] = <F_k, X^-1 W X^-1>."""
    Xi = _spd_inv(X)
    R = float(np.einsum("ij,ji->", W, Xi))
    XiWXi = Xi @ W @ Xi
    g = np.array([-float(np.einsum("ij,ij->", Fk[k], XiWXi)) for k in range(len(Fk))])
    return R, Xi, XiWXi, g


def _assemble_X(H0, Fk, w):
    X = H0.copy()
    for k in range(len(Fk)):
        if w[k] != 0.0:
            X += w[k] * Fk[k]
    return X


def _arisk_curve(X, dX, W):
    """Coefficients (lam, d) s.t. tr[W (X + g dX)^-1] = sum_i d_i / (1 + g lam_i),
    via the generalized eig of dX wrt X (X SPD, X+g dX SPD for g in [0,1]).  Makes
    the FW line search O(r) per evaluation (one Cholesky + one eig, no per-g
    inversion)."""
    L = cholesky(0.5 * (X + X.T), lower=True, check_finite=False)
    Y = solve_triangular(L, dX, lower=True, check_finite=False)
    M = solve_triangular(L, Y.T, lower=True, check_finite=False).T
    M = 0.5 * (M + M.T)
    lam, U = eigh(M)
    YW = solve_triangular(L, W, lower=True, check_finite=False)
    Wt = solve_triangular(L, YW.T, lower=True, check_finite=False).T
    Wt = 0.5 * (Wt + Wt.T)
    d = np.einsum("ir,ij,jr->r", U, Wt, U)      # diag(U^T (L^-1 W L^-T) U) >= 0
    return lam, d


def bank_simplex_oracle(H0, Fk, W, tol=1e-9, max_iter=200):
    """R*_s = min over the K-simplex of R_s(v)=tr[W X_s(v)^-1] (R25 SS4).

    Convex (matrix-fractional).  Primary solver: SLSQP on the K-simplex with the
    analytic gradient dR/dv_k = -tr[W X^-1 F_k X^-1] (few evals for K=8); Frank-
    Wolfe with the O(r) eig line search is the robust fallback.  Returns
    (Rstar, v*)."""
    K = len(Fk)
    v0 = np.full(K, 1.0 / K)
    try:
        from scipy.optimize import minimize

        def fun(v):
            R, _ = arisk(_assemble_X(H0, Fk, v), W)
            return R

        def jac(v):
            _, _, _, ng = _arisk_and_grads(_assemble_X(H0, Fk, v), Fk, W)
            return ng                              # dR/dv_k

        res = minimize(fun, v0, jac=jac, method="SLSQP",
                       bounds=[(0.0, 1.0)] * K,
                       constraints=[{"type": "eq",
                                     "fun": lambda v: v.sum() - 1.0,
                                     "jac": lambda v: np.ones(K)}],
                       options={"maxiter": 100, "ftol": 1e-12})
        v = _simplex_proj(res.x)
        R, _ = arisk(_assemble_X(H0, Fk, v), W)
        if np.isfinite(R):
            return R, v
    except Exception:
        pass
    return _bank_simplex_oracle_fw(H0, Fk, W, tol=tol, max_iter=max_iter)


def _bank_simplex_oracle_fw(H0, Fk, W, tol=1e-9, max_iter=200):
    """Frank-Wolfe fallback for the bank-simplex oracle (O(r) eig line search)."""
    K = len(Fk)
    v = np.full(K, 1.0 / K)
    X = _assemble_X(H0, Fk, v)
    R, Xi, XiWXi, grad = _arisk_and_grads(X, Fk, W)
    for _ in range(max_iter):
        j = int(np.argmin(grad))                 # LMO: vertex e_j
        d = -v.copy(); d[j] += 1.0
        gap = -float(grad @ d)
        if gap <= tol * max(1.0, R):
            break
        Xj = _assemble_X(H0, Fk, np.eye(K)[j])
        dX = Xj - X
        lam, dc = _arisk_curve(X, dX, W)
        curve = lambda g: float(np.sum(dc / (1.0 + g * lam)))   # convex on [0,1]
        lo, hi = 0.0, 1.0
        gr = 0.5 * (np.sqrt(5.0) - 1.0)
        a, b = hi - gr * (hi - lo), lo + gr * (hi - lo)
        fa, fb = curve(a), curve(b)
        for _ in range(50):
            if fa < fb:
                hi, b, fb = b, a, fa
                a = hi - gr * (hi - lo); fa = curve(a)
            else:
                lo, a, fa = a, b, fb
                b = lo + gr * (hi - lo); fb = curve(b)
            if hi - lo < 1e-9:
                break
        g_star = 0.5 * (lo + hi)
        v = (1 - g_star) * v + g_star * np.eye(K)[j]
        X = _assemble_X(H0, Fk, v)
        R, Xi, XiWXi, grad = _arisk_and_grads(X, Fk, W)
    return R, v


def _normalize_W(sm):
    """Numerical conditioning (T-D, engineering; criterion-invariant).

    The standardized-maximin decision depends on W only through the RATIOS
    r_s(w)=R_s(w)/R*_s, which are invariant to a global positive rescaling of W
    (R_s and R*_s both scale).  Real dead-time information matrices are ~1e10, so
    R(uniform) ~ 1e-7, where the oracle SLSQP (absolute ftol) and the maximin
    both stall (observed on the T-B banks: R*_s over-estimated -> t_cont<1, and
    the maximin runs its full iteration budget, ~200 s).  We rescale W so
    R(uniform) ~ O(1); t_cont, w*, alpha, the KKT residual and every A_j are
    numerically unchanged (well-scaled synthetic fixtures: kappa ~ 1, no-op)."""
    K = len(sm.F[0])
    X0 = _assemble_X(sm.H0[0], sm.F[0], np.full(K, 1.0 / K))
    R0, _ = arisk(X0, sm.W)
    if np.isfinite(R0) and R0 > 0:
        sm.W = sm.W / R0
    return sm


def compute_oracles(sm, tol=1e-10):
    """Fill sm.Rstar, sm.vstar with the per-scenario bank-simplex oracles.

    W is first normalized (criterion-invariant conditioning; see _normalize_W)."""
    sm = _normalize_W(sm)
    S = len(sm.H0)
    Rstar = np.empty(S)
    vstar = []
    for s in range(S):
        R, v = bank_simplex_oracle(sm.H0[s], sm.F[s], sm.W, tol=tol)
        Rstar[s] = R
        vstar.append(v)
    sm.Rstar = Rstar
    sm.vstar = vstar
    return sm


# ============================================================================ #
#  7.  standardized maximin (R25.1) + KKT certificate (R25.2)                   #
# ============================================================================ #

def _rs_and_grad(H0, Fk, W, Rstar, w):
    """r_s(w)=R_s(w)/R*_s and its gradient dr_s/dw_k = -g_{s,k},
    g_{s,k}=tr[W X^-1 F_k X^-1]/R*_s.  Returns (r_s, g_vec[K], R_s).
    Fast marginal via <F_k, X^-1 W X^-1> (see _arisk_and_grads)."""
    X = _assemble_X(H0, Fk, w)
    R, Xi, XiWXi, ng = _arisk_and_grads(X, Fk, W)
    g = (-ng) / Rstar               # g_{s,k} = tr[W X^-1 F_k X^-1]/R*_s (positive)
    return R / Rstar, g, R


@dataclass
class MaximinResult:
    w: np.ndarray
    t_cont: float
    alpha: np.ndarray
    lam: float
    kkt_residual: float
    active_scenarios: np.ndarray
    support: np.ndarray
    r_s: np.ndarray
    g: np.ndarray                # (S, K) marginal values at w*
    converged: bool
    n_iter: int


def standardized_maximin(sm, tol=1e-9, tie_tol=1e-7, outer_iter=120,
                         inner_iter=120, seed=0, refine_maxiter=None,
                         fast_linalg=False):
    """Solve min_{w in simplex} max_s r_s(w), r_s=R_s(w)/R*_s (R25.1).

    Primal-dual saddle formulation (exact, gives the least-favorable alpha):
        min_{w in Delta_K} max_{alpha in Delta_S}  sum_s alpha_s r_s(w).
    The dual function g(alpha) = min_{w in Delta_K} sum_s alpha_s r_s(w) is
    CONCAVE in alpha with supergradient r(w(alpha)); we maximize it over the
    S-simplex by Frank-Wolfe (LMO = e_{argmax_s r_s}, step 2/(k+2)).  The inner
    convex problem min_w sum_s alpha_s r_s(w) is solved to tight accuracy by
    projected-gradient with backtracking (grad = -sum_s alpha_s g_s(w)).

    At the saddle, w* = w(alpha*) minimizes sum_s alpha*_s r_s over the simplex,
    so the KKT stationarity (R25.2)  sum_s alpha*_s g_{s,k} = lambda on supp(w),
    <= lambda off support, holds by construction -- alpha* IS the least-favorable
    scenario distribution.  K=8, S=16 -> small + robust.

    fast_linalg (R27 §3.2, output-equivalent): evaluate rs_g with BATCHED linear
    algebra over the S scenarios (one stacked (S,r,r) inverse instead of S scipy
    Cholesky-solves).  Same objective, same math; verified against the reference
    (fast_linalg=False) to the R27 equivalence bar (|t|<=1e-8, ||w||_1<=1e-6, KKT
    within tol, IDENTICAL 972-row realization) in test_rlmi.  Default off."""
    H0, F, W, Rstar = sm.H0, sm.F, sm.W, sm.Rstar
    S, K = len(H0), len(F[0])

    if fast_linalg:
        H0s = np.stack([np.ascontiguousarray(h) for h in H0])       # (S,r,r)
        Fs = np.stack([np.stack(F[s]) for s in range(S)])           # (S,K,r,r)
        Rst = np.asarray(Rstar, dtype=np.float64)

        def rs_g(w):
            X = H0s + np.einsum("k,skij->sij", w, Fs)               # (S,r,r)
            X = 0.5 * (X + np.transpose(X, (0, 2, 1)))
            Xi = np.linalg.inv(X)                                   # batched
            R = np.einsum("ij,sji->s", W, Xi)                       # tr[W Xi]
            XiWXi = Xi @ W @ Xi                                     # (S,r,r)
            gpos = np.einsum("skij,sij->sk", Fs, XiWXi)             # tr[F XiWXi]
            return R / Rst, gpos / Rst[:, None]
    else:
        def rs_g(w):
            rs = np.empty(S); gg = np.empty((S, K))
            for s in range(S):
                rs[s], gg[s], _ = _rs_and_grad(H0[s], F[s], W, Rstar[s], w)
            return rs, gg

    def inner_min(alpha, w0):
        """argmin_{w in Delta_K} sum_s alpha_s r_s(w) (convex), warm-started."""
        w = w0.copy()
        rs, gg = rs_g(w)
        phi = float(alpha @ rs)
        step = 1.0
        for _ in range(inner_iter):
            grad = -(alpha[:, None] * gg).sum(axis=0)      # d/dw sum alpha_s r_s
            gn = np.linalg.norm(grad)
            if gn < 1e-14 or not np.all(np.isfinite(grad)):
                break
            improved = False
            for _bt in range(40):
                w_new = _simplex_proj(w - step * grad)
                rs_n, gg_n = rs_g(w_new)
                phi_n = float(alpha @ rs_n)
                if phi_n <= phi - 1e-4 * step * (grad @ (w - w_new)):
                    if phi - phi_n < 1e-14 * max(1.0, abs(phi)):
                        w, rs, gg, phi = w_new, rs_n, gg_n, phi_n
                        improved = False       # converged: stop outer inner-loop
                        break
                    w, rs, gg, phi = w_new, rs_n, gg_n, phi_n
                    step = min(step * 1.5, 1e6)   # cap growth (avoid overflow)
                    improved = True
                    break
                step *= 0.5
            if not improved or step < 1e-12:
                break
        return w, rs, gg, phi

    # PRIMARY primal solver: SLSQP on the epigraph  min_{w,t} t  s.t.
    # r_s(w) <= t, w in simplex (analytic Jacobians), from the uniform start.
    # The problem is convex (r_s convex), so SLSQP converges to the global min in
    # few evals.  Refine budget scales down with r (each eval costs S r^3): full
    # accuracy at small r (the certified tests), bounded at large r (benchmark).
    if refine_maxiter is None:
        r = sm.r
        refine_maxiter = 200 if r <= 64 else (60 if r <= 128 else 30)
    w = np.full(K, 1.0 / K)
    w = _slsqp_epigraph(rs_g, w, S, K, maxiter=refine_maxiter)
    rs, gg = rs_g(w)
    t_cont = float(rs.max())

    # robustness fallback: if SLSQP stalled (t below the best oracle-consistent
    # bound is impossible; t must be >= 1), run the primal-dual FW and keep the
    # better primal.
    n_iter = 1
    if not np.isfinite(t_cont) or t_cont < 1.0 - 1e-6:
        alpha = np.full(S, 1.0 / S)
        wf = np.full(K, 1.0 / K)
        for k in range(outer_iter):
            n_iter += 1
            wf, rsf, ggf, phi = inner_min(alpha, wf)
            j = int(np.argmax(rsf))
            gap = float(rsf[j] - alpha @ rsf)
            if gap <= tol:
                break
            gamma = 2.0 / (k + 3.0)
            alpha = (1.0 - gamma) * alpha
            alpha[j] += gamma
        wf = _slsqp_epigraph(rs_g, wf, S, K, maxiter=refine_maxiter)
        rsf, ggf = rs_g(wf)
        if rsf.max() <= t_cont or not np.isfinite(t_cont):
            w, rs, gg, t_cont = wf, rsf, ggf, float(rsf.max())

    support = np.where(w > tie_tol)[0]
    # active-set detection: scenarios within an adaptive band of t
    band = max(1e-7, 1e-6 * t_cont)
    active = np.where(rs >= t_cont - band)[0]
    if active.size == 0:
        active = np.array([int(rs.argmax())])
    alpha_out, lam_out, res_out = _polish_alpha(gg, active, support, rs,
                                                t_cont, K, S)
    # widen the active band once if the certificate is loose (a barely-inactive
    # scenario is carrying the least-favorable mass)
    if res_out > 1e-6:
        band2 = max(band, 1e-3 * t_cont)
        active2 = np.where(rs >= t_cont - band2)[0]
        a2, l2, r2 = _polish_alpha(gg, active2, support, rs, t_cont, K, S)
        if r2 < res_out:
            alpha_out, lam_out, res_out, active = a2, l2, r2, active2
    conv = res_out <= 1e-6
    return MaximinResult(w=w, t_cont=t_cont, alpha=alpha_out, lam=lam_out,
                         kkt_residual=res_out, active_scenarios=active,
                         support=support, r_s=rs, g=gg, converged=bool(conv),
                         n_iter=n_iter)


def _slsqp_epigraph(rs_g, w0, S, K, maxiter=150):
    """Refine w by SLSQP on  min_{w,t} t  s.t. t - r_s(w) >= 0, sum w = 1, w>=0.
    Analytic Jacobians (dr_s/dw = -g_s).  Falls back to w0 on any failure."""
    try:
        from scipy.optimize import minimize
    except Exception:
        return w0
    # 1-entry memo: SLSQP evaluates the constraint value and Jacobian separately
    # at the same w, so cache the last rs_g(w) (each call costs S r^3).
    _cache = {"key": None, "val": None}

    def rs_g_c(w):
        key = w.tobytes()
        if _cache["key"] != key:
            _cache["key"] = key
            _cache["val"] = rs_g(w)
        return _cache["val"]

    rs0, _ = rs_g_c(w0)
    z0 = np.concatenate([w0, [float(rs0.max())]])

    def unpack(z):
        return z[:K], z[K]

    def obj(z):
        return z[K]

    def obj_jac(z):
        g = np.zeros(K + 1); g[K] = 1.0
        return g

    def con_ineq(z):
        w, t = unpack(z)
        rs, _ = rs_g_c(w)
        return t - rs                       # >= 0

    def con_ineq_jac(z):
        w, t = unpack(z)
        rs, gg = rs_g_c(w)
        J = np.zeros((S, K + 1))
        J[:, :K] = gg                        # d(t - r_s)/dw_k = -dr_s/dw_k = g_{s,k}
        J[:, K] = 1.0
        return J

    def con_eq(z):
        w, t = unpack(z)
        return np.array([w.sum() - 1.0])

    def con_eq_jac(z):
        J = np.zeros((1, K + 1)); J[0, :K] = 1.0
        return J

    bounds = [(0.0, 1.0)] * K + [(None, None)]
    try:
        res = minimize(obj, z0, jac=obj_jac, bounds=bounds,
                       constraints=[{"type": "ineq", "fun": con_ineq,
                                     "jac": con_ineq_jac},
                                    {"type": "eq", "fun": con_eq,
                                     "jac": con_eq_jac}],
                       method="SLSQP", options={"maxiter": maxiter, "ftol": 1e-12})
        w_new = _simplex_proj(res.x[:K])
        rs_new, _ = rs_g(w_new)
        rs_old, _ = rs_g(w0)
        if rs_new.max() <= rs_old.max() + 1e-10:
            return w_new
    except Exception:
        pass
    return w0


def _kkt_from_alpha(alpha, gg, support, rs, t, K):
    """R25.2 residual for a given alpha: stationarity on/off support +
    complementarity + simplex."""
    agk = (alpha[:, None] * gg).sum(axis=0)          # sum_s alpha_s g_{s,k}
    if support.size:
        lam = float(agk[support].mean())
        stat_supp = float(np.max(np.abs(agk[support] - lam)))
    else:
        lam = float(agk.max()); stat_supp = 0.0
    off = np.setdiff1d(np.arange(K), support)
    stat_off = float(np.max(np.maximum(agk[off] - lam, 0.0))) if off.size else 0.0
    compl = float(np.max(alpha * np.abs(rs - t))) if alpha.size else 0.0
    simplex = abs(float(alpha.sum()) - 1.0)
    return alpha, lam, stat_supp + stat_off + compl + simplex


def _polish_alpha(gg, active, support, rs, t, K, S, alpha_fw=None):
    """Least-favorable alpha on active scenarios equalizing sum_s alpha_s g_{s,k}
    over supp(w) (R25.2), alpha>=0, sum=1; NNLS on the equalization system."""
    A = active if active.size else np.arange(S)
    P = support if support.size else np.arange(K)
    G = gg[np.ix_(A, P)]                              # |A| x |P|
    nA = A.size
    rows_ls, rhs = [], []
    for a_i in range(len(P) - 1):
        rows_ls.append(G[:, a_i] - G[:, a_i + 1]); rhs.append(0.0)
    big = 1e3
    rows_ls.append(big * np.ones(nA)); rhs.append(big)
    M = np.array(rows_ls); b = np.array(rhs)
    try:
        aA, _ = nnls(M, b)
    except Exception:
        aA = np.ones(nA)
    ss = aA.sum()
    aA = aA / ss if ss > 0 else np.full(nA, 1.0 / nA)
    alpha = np.zeros(S); alpha[A] = aA
    return _kkt_from_alpha(alpha, gg, support, rs, t, K)


def _simplex_proj(v):
    """Euclidean projection onto the probability simplex (Duchi et al. 2008).
    Non-finite inputs are sanitized to the uniform point (defensive)."""
    v = np.asarray(v, dtype=np.float64)
    if not np.all(np.isfinite(v)):
        return np.full(v.size, 1.0 / v.size)
    u = np.sort(v)[::-1]
    css = np.cumsum(u) - 1.0
    ind = np.arange(1, v.size + 1)
    cond = u - css / ind > 0
    rho = ind[cond][-1]
    theta = css[cond][-1] / rho
    return np.maximum(v - theta, 0.0)


# ============================================================================ #
#  8.  lexicographic closest-to-knob tie-break (R25 P-B, eq. w-dagger)          #
# ============================================================================ #

def tie_break_to_knob(sm, mm, knob_index=0, tie_tol=1e-6):
    """Among (near-)minimizers within tie_tol of t_cont, pick w closest to e_0
    in l2 (R25.1 step 4; strictly convex secondary problem).

    Projected-gradient descent of ||w - e0||^2 restricted to the near-optimal
    set { w in simplex : max_s r_s(w) <= t_cont + tie_tol }; the feasibility of
    the near-optimal set is maintained by a projected step that only accepts
    moves keeping max_s r_s within tolerance."""
    H0, F, W, Rstar = sm.H0, sm.F, sm.W, sm.Rstar
    S, K = len(H0), len(F[0])
    e0 = np.zeros(K); e0[knob_index] = 1.0
    t_bar = mm.t_cont + tie_tol

    def maxr(w):                                   # value-only (Cholesky)
        m = 0.0
        for s in range(S):
            R, _ = arisk(_assemble_X(H0[s], F[s], w), W)
            m = max(m, R / Rstar[s])
        return m

    max_it = 120 if sm.r <= 64 else 40
    w = mm.w.copy()
    step = 0.4
    rejects = 0
    for _ in range(max_it):
        # projection of e0 onto the near-optimal band: descend ||w-e0||^2
        # (direction -2(w-e0)), accept only moves that stay within the band.
        cand = _simplex_proj(w - step * 2.0 * (w - e0))
        if (maxr(cand) <= t_bar
                and np.linalg.norm(cand - e0) < np.linalg.norm(w - e0) - 1e-12):
            w = cand
            rejects = 0
        else:
            step *= 0.5
            rejects += 1
            # unique / separated optimum: consecutive rejections mean e0 is
            # outside the near-optimal band -> w* stands (no flat set to exploit).
            if rejects >= 5 or step < 1e-7:
                break
    return w


# ============================================================================ #
#  9.  R23 rank-one residual disclosure (eps_1, m, eps_1^2/2m; branch-fenced)    #
# ============================================================================ #

def r23_residual(bank, x, B, nu, c, W=None, eps=1e-9, knob=True):
    """Per-scenario R23 disclosure for the KNOB bank at scene x (R23 SS2.5):
      leverage  ell_i = q_i^T V^{-1} q_i  (D-opt leverage on the projected V),
      h_i       = ell_i * nu * J'(rho_i),
      eps_1^2   = sum_i (h_i - hbar)^2                (T1.10 budget-only interior),
      m         = min_i (nu*(-J''(rho_i))) * ell_i    (R23 kernel-curvature modulus),
      U_1       = eps_1^2 / (2 m)                     (T1.6) IF on-branch (m>0 &
                  all rho_i <= rho_c); else BRANCH_VIOLATION -> UNAVAILABLE.
    V = eps I + sum_i nu J(rho_i) q_i q_i^T over the bank rows (projected)."""
    r = B.shape[1]
    rho_c = kernel_rho_c(c)
    loads = np.maximum(bank.rows @ x, 1e-12)
    Q = (bank.rows / loads[:, None]) @ B                  # (R, r) projected q_i
    V = eps * np.eye(r) + (Q.T * (nu * kernel_J(loads, c))) @ Q
    Vi = _spd_inv(V)
    ell = np.einsum("ir,rs,is->i", Q, Vi, Q)              # leverage
    Jp = kernel_Jp(loads, c)
    h = ell * (nu * Jp)
    hbar = float(h.mean())
    eps1_sq = float(np.sum((h - hbar) ** 2))
    Jpp = kernel_Jpp(loads, c)
    m_vec = (nu * (-Jpp)) * ell
    m = float(m_vec.min())
    on_branch = bool(np.all(loads <= rho_c) and m > 0)
    out = {"eps_1": float(np.sqrt(max(eps1_sq, 0.0))), "eps_1_sq": eps1_sq,
           "m": m, "max_load": float(loads.max()), "rho_c": float(rho_c),
           "n_beyond": int(np.sum(loads > rho_c)), "on_branch": on_branch}
    if on_branch:
        out["U1_bound"] = eps1_sq / (2.0 * m)
        out["branch_status"] = "ON_BRANCH"
    else:
        out["U1_bound"] = None
        out["branch_status"] = "BRANCH_VIOLATION"
    return out


# ============================================================================ #
#  10.  materialization: union measure -> exact 972 rows (R25 SS9-SS10)          #
# ============================================================================ #

def _atom_key(row, decimals=6):
    return hashlib.sha1(np.round(row, decimals).tobytes()).hexdigest()


@dataclass
class Materialization:
    rows: np.ndarray                 # (972, n) exact acquired rows (block-ordered)
    counts: np.ndarray               # (A,) per-union-atom counts
    atom_rows: np.ndarray            # (A, n) union atom pool
    atom_bankw: np.ndarray           # (A, K) continuous bank contribution weights
    what: np.ndarray                 # (K,) realized bank fractions n_{.,k}/972
    guards: dict
    flag: str                        # "OK" | "MIXTURE_MATERIALIZATION_FAIL"
    n_exchange: int
    gamma: np.ndarray = None         # (A,) per-atom source-power setting (R25 SS8)


def _union_measure(banks, w):
    """Merge duplicate atoms across banks (weight addition; R25 SS10).  Returns
    (atom_rows (A,n), xi (A,), bankw (A,K), incident (A,), peak (A,), adm (A,))."""
    K = len(banks)
    acc = {}
    for k, bk in enumerate(banks):
        for r in range(bk.rows.shape[0]):
            key = _atom_key(bk.rows[r])
            if key not in acc:
                acc[key] = {"row": bk.rows[r], "xi": 0.0, "bw": np.zeros(K),
                            "inc": float(bk.incident[r]), "peak": float(bk.peak[r]),
                            "adm": bool(bk.admission[r])}
            contrib = w[k] * float(bk.weights[r])
            acc[key]["xi"] += contrib
            acc[key]["bw"][k] += contrib
    keys = sorted(acc.keys())                     # deterministic atom order
    A = len(keys)
    n = banks[0].rows.shape[1]
    rows = np.zeros((A, n)); xi = np.zeros(A); bw = np.zeros((A, K))
    inc = np.zeros(A); peak = np.zeros(A); adm = np.zeros(A, dtype=bool)
    for i, key in enumerate(keys):
        rows[i] = acc[key]["row"]; xi[i] = acc[key]["xi"]; bw[i] = acc[key]["bw"]
        inc[i] = acc[key]["inc"]; peak[i] = acc[key]["peak"]; adm[i] = acc[key]["adm"]
    s = xi.sum()
    if s > 0:
        xi = xi / s
        bw = bw / s
    return rows, xi, bw, inc, peak, adm


def _largest_remainder(target, total):
    """Largest-remainder (Hamilton) apportionment to integer counts summing to
    total; deterministic tie-break by index."""
    base = np.floor(target).astype(np.int64)
    rem = int(total - base.sum())
    if rem > 0:
        frac = target - base
        order = np.lexsort((np.arange(target.size), -frac))
        base[order[:rem]] += 1
    elif rem < 0:
        frac = target - base
        order = np.lexsort((np.arange(target.size), frac))
        for i in range(-rem):
            j = order[i]
            if base[j] > 0:
                base[j] -= 1
    return base


def _int_dose_descent(counts, rows, supp, M_rows, target, max_moves=1500):
    """Strengthen the integer stage past the local-min of the pairwise exchange:
    single-count gradient descent on the L2 dose penalty (drop from the atom most
    aligned with the over-dose, add to the admissible atom most aligned with the
    under-dose).  Bounded; helps compact/bar/multiscale banks reach the ~1/k
    structured-drop floor.  (Scattered/ring/weighted measures can stall at a
    combinatorial local minimum — those honestly fall through to the guard
    check and the L0 fallback; see the T-D phase-1.5 report.)"""
    # bounded; skip for very large atom pools (weighted FW measure) where the
    # per-move matvec is expensive and the greedy stalls anyway.
    if rows.shape[0] > 5000:
        return counts
    dose = (counts[:, None] * rows).sum(0).astype(np.float64)
    best_dev = float(np.abs(dose / dose.mean() - 1.0).max()) if dose.mean() > 0 else np.inf
    best = counts.copy()
    for _ in range(max_moves):
        dbar = dose.mean()
        if dbar <= 0:
            break
        rel = dose / dbar - 1.0
        dev = float(np.abs(rel).max())
        if dev < best_dev:                              # track best (never worse)
            best_dev = dev; best = counts.copy()
        if dev <= target:
            break
        score = rows @ rel
        src = int(np.argmax(np.where(counts > 0, score, -np.inf)))
        dst = int(np.argmin(np.where(supp, score, np.inf)))
        if src == dst or score[src] - score[dst] <= 1e-9:
            break
        counts[src] -= 1; counts[dst] += 1
        dose += rows[dst] - rows[src]
    return best


def _power_repair(counts, rows, M_rows, dose_band, lo=0.94, hi=1.06, iters=800):
    """R25 SS8 per-row SOURCE-SETTING repair: bounded per-atom power trim
    gamma_a in [lo, hi] that pulls per-pixel dose into the +-dose_band band,
    minimal-norm (starts at 1).  Fast projected-gradient (O(nnz)/iter) on the
    band-violation penalty; returns gamma over ALL atoms (1 off support) and the
    achieved dose_dev.  A least-squares power repair per R25 SS8 per-row settings;
    only closes the residual granularity gap the integer stage leaves ~1/k."""
    use = np.where(counts > 0)[0]
    if use.size == 0:
        return np.ones(counts.shape[0]), np.inf
    Ru = counts[use][:, None] * rows[use]              # (A_used, n)
    g = np.ones(use.size)

    def dev_of(gv):
        d = Ru.T @ gv
        dm = d.mean()
        return np.inf if dm <= 0 else float(np.abs(d / dm - 1.0).max())

    cur = dev_of(g)
    best = (cur, g.copy())
    step = None
    for _ in range(iters):
        if cur <= dose_band:
            break
        dose = Ru.T @ g
        rel = dose / dose.mean() - 1.0
        viol = np.sign(rel) * np.maximum(np.abs(rel) - dose_band + 0.004, 0.0)
        grad = Ru @ viol
        gn = float(np.abs(grad).max())
        if gn < 1e-14:
            break
        if step is None:
            step = 0.7 / gn
        # backtracking line search: accept only steps that reduce dev
        accepted = False
        for _bt in range(20):
            gtry = np.clip(g - step * grad, lo, hi)
            dtry = dev_of(gtry)
            if dtry < cur - 1e-9:
                g, cur = gtry, dtry
                if cur < best[0]:
                    best = (cur, g.copy())
                step *= 1.3
                accepted = True
                break
            step *= 0.5
        if not accepted:
            break
    gfull = np.ones(counts.shape[0])
    gfull[use] = best[1]
    return gfull, best[0]


def materialize(banks, w, xhat, M_rows=972, dose_band=0.05,
                incident_budget=None, peak_cap=1536.0, side=32,
                block_seed=650777, max_exchange=4000, power_repair=True,
                gamma_lo=0.94, gamma_hi=1.06):
    """Exact 972-row materialization from the union design measure (R25 SS9).

    Steps: merge duplicate atoms (weight addition) -> target p_a = M_rows*xi(w) ->
    largest-remainder rounding to total M_rows -> deterministic exchange enforcing
    admission/peak, incident budget, and +-band per-pixel dose -> frozen-seed
    block-order randomization.  Returns a Materialization; sets flag
    MIXTURE_MATERIALIZATION_FAIL if no compliant exact design is found."""
    K = len(banks)
    # R27 §2: a PURE bank realization IS the bank's frozen exact-972 admission
    # witness (dose-band-compliant under the bank's nominal powers).  Use it
    # verbatim rather than re-deriving it -- guarantees the manifest-certified
    # design and gamma=1 (nominal powers).
    wa = np.asarray(w, dtype=np.float64)
    nz = np.where(wa > 1e-12)[0]
    if nz.size == 1 and abs(float(wa[nz[0]]) - 1.0) < 1e-9:
        wm = (banks[int(nz[0])].meta or {}).get("witness_mult")
        if wm is not None:
            mv = _materialize_from_witness(banks, int(nz[0]), np.asarray(wm),
                                           M_rows, dose_band, incident_budget,
                                           peak_cap, block_seed)
            if mv is not None:
                return mv
    rows, xi, bw, inc, peak, adm = _union_measure(banks, w)
    A = rows.shape[0]
    n = rows.shape[1]

    # hard admission: zero any inadmissible / over-peak atom, renormalize target
    ok = adm & (peak <= peak_cap + 1e-9)
    if not ok.any():
        return _mat_fail(banks, K, "MIXTURE_MATERIALIZATION_FAIL")
    xi_ok = np.where(ok, xi, 0.0)
    if xi_ok.sum() <= 0:
        return _mat_fail(banks, K, "MIXTURE_MATERIALIZATION_FAIL")
    xi_ok = xi_ok / xi_ok.sum()
    # exchange destinations = admissible atoms IN the mixture support (xi>0); the
    # realized design stays a materialization of the union measure (R25 SS9), so
    # attribution to banks is well-defined and t_real reflects the same mixture.
    supp = ok & (xi > 1e-15)
    counts = _largest_remainder(M_rows * xi_ok, M_rows)
    counts[~ok] = 0
    if counts.sum() != M_rows:                    # repair rounding on admissible set
        deficit = M_rows - counts.sum()
        adm_idx = np.where(ok)[0]
        frac = M_rows * xi_ok - counts
        order = adm_idx[np.argsort(-frac[adm_idx])]
        i = 0
        while deficit != 0 and order.size:
            j = order[i % order.size]
            if deficit > 0:
                counts[j] += 1; deficit -= 1
            elif counts[j] > 0:
                counts[j] -= 1; deficit += 1
            i += 1

    def dose_dev(cnts):
        dose = (cnts[:, None] * rows).sum(axis=0) / M_rows
        dmean = dose.mean()
        if dmean <= 0:
            return np.inf, dose
        return float(np.abs(dose / dmean - 1.0).max()), dose

    def incident_total(cnts):
        return float((cnts * inc).sum())

    n_exchange = 0
    B_inc = incident_budget

    # -- budget exchange: move a count from the highest-incident atom to the
    #    lowest-incident admissible atom (deterministic; R25 SS9 constrained rounding)
    if B_inc is not None:
        while incident_total(counts) > B_inc + 1e-9 and n_exchange < max_exchange:
            hi = np.where(counts > 0, inc, -np.inf)
            src = int(np.argmax(hi))
            lo = np.where(supp, inc, np.inf)
            dst = int(np.argmin(lo))
            if src == dst or inc[src] <= inc[dst] + 1e-12:
                break
            counts[src] -= 1; counts[dst] += 1
            n_exchange += 1

    # -- dose-band exchange: while over band, move one count from an atom
    #    covering the most over-dosed pixel to an admissible atom covering the
    #    most under-dosed pixel, choosing the move that most reduces the L2 dose
    #    penalty (deterministic argmin; bounded).
    dev, dose = dose_dev(counts)
    cover = rows > 0
    guard_it = 0
    while dev > dose_band and guard_it < max_exchange:
        guard_it += 1; n_exchange += 1
        dmean = dose.mean()
        rel = dose / dmean - 1.0
        j_over = int(np.argmax(rel)); j_under = int(np.argmin(rel))
        srcs = np.where((counts > 0) & cover[:, j_over])[0]
        dsts = np.where(supp & cover[:, j_under])[0]
        if srcs.size == 0 or dsts.size == 0:
            break
        # pick the (src,dst) most reducing sum of squared band-excess
        band = dose_band
        def pen(dvec):
            rl = np.abs(dvec / dvec.mean() - 1.0)
            return float((np.maximum(rl - band + 0.002, 0.0) ** 2).sum())
        cur = pen(dose)
        best = (0.0, None)
        for s_ in srcs[:32]:
            for d_ in dsts[:32]:
                if s_ == d_:
                    continue
                nd = dose + (rows[d_] - rows[s_]) / M_rows
                dp = pen(nd) - cur
                if dp < best[0] - 1e-15:
                    best = (dp, (s_, d_))
        if best[1] is None:
            break
        s_, d_ = best[1]
        counts[s_] -= 1; counts[d_] += 1
        dev, dose = dose_dev(counts)

    # -- materializer v2 (R25 SS8/SS9): if the deterministic exchange left the
    #    dose above the band, (i) strengthen the integer stage by gradient
    #    descent toward the ~1/k structured-drop floor, then (ii) apply the
    #    bounded per-row source-power repair gamma in [gamma_lo, gamma_hi].
    dev, dose = dose_dev(counts)
    dev_pre = dev
    gamma = np.ones(A)
    if power_repair and dev > dose_band:
        counts = _int_dose_descent(counts, rows, supp, M_rows,
                                   target=max(dose_band * 0.9, 0.045))
        dev, dose = dose_dev(counts)
        gamma, dev = _power_repair(counts, rows, M_rows, dose_band,
                                   lo=gamma_lo, hi=gamma_hi)

    # rows carrying the per-row source setting (gamma * atom row)
    rows_eff = gamma[:, None] * rows

    # final guard recompute ON THE TRIMMED (source-set) rows
    def dose_dev_eff(cnts):
        d = (cnts[:, None] * rows_eff).sum(0) / M_rows
        dm = d.mean()
        return (np.inf if dm <= 0 else float(np.abs(d / dm - 1.0).max()))
    dev = dose_dev_eff(counts)
    total_ok = int(counts.sum()) == M_rows
    peak_real = float((rows_eff[counts > 0].max()) if (counts > 0).any() else 0.0)
    inc_real = float((counts * inc * gamma).sum())
    inc_ok = (B_inc is None) or (inc_real <= B_inc + 1e-6)
    dose_ok = dev <= dose_band + 1e-9
    adm_ok = bool(np.all(counts[~ok] == 0))
    peak_ok = peak_real <= peak_cap + 1e-9
    guards = {"exact_count": int(counts.sum()), "count_ok": bool(total_ok),
              "dose_dev": dev, "dose_dev_pretrim": dev_pre, "dose_ok": bool(dose_ok),
              "dose_band": dose_band, "incident_total": inc_real,
              "incident_budget": B_inc, "incident_ok": bool(inc_ok),
              "peak": peak_real, "peak_ok": bool(peak_ok), "peak_cap": peak_cap,
              "admission_ok": adm_ok, "n_exchange": n_exchange,
              "gamma_min": float(gamma[counts > 0].min()) if (counts > 0).any() else 1.0,
              "gamma_max": float(gamma[counts > 0].max()) if (counts > 0).any() else 1.0,
              "gamma_std": float(gamma[counts > 0].std()) if (counts > 0).any() else 0.0}
    if not (total_ok and dose_ok and inc_ok and adm_ok and peak_ok):
        fail = _mat_fail(banks, K, "MIXTURE_MATERIALIZATION_FAIL")
        fail.guards.update(guards)
        return fail

    # expand to 972 rows + frozen-seed block-order randomization (R25 SS10);
    # each acquired row carries its atom's source setting gamma (rows_eff).
    expanded = np.repeat(np.arange(A), counts)
    perm = np.random.default_rng(block_seed).permutation(expanded.size)
    block = expanded[perm]
    acquired = rows_eff[block]
    what = np.zeros(K)
    for i in range(A):
        if counts[i] > 0 and bw[i].sum() > 0:
            what += counts[i] * bw[i] / bw[i].sum()
    what = what / M_rows
    return Materialization(rows=acquired, counts=counts, atom_rows=rows,
                           atom_bankw=bw, what=what, guards=guards, flag="OK",
                           n_exchange=n_exchange, gamma=gamma)


def materialize_fw_relaxed(rows, weights, M_rows=972, peak_cap=1536.0,
                           incident_budget=None, block_seed=650777):
    """R28 §3 dose-RELAXED FW diagnostic materialization
    (NONDEPLOYABLE_DOSE_RELAXED_FW_DIAGNOSTIC).  Largest-remainder integer
    multiplicities from the FW design measure at its FROZEN nominal amplitudes
    (preserving the FW information direction -- NO coverage/Sinkhorn dose
    uniformization), with the SAME exact-972 count, row-wise peak cap, incident-
    budget cap, nonnegativity and admission guards, but NO +-5% per-pixel dose
    band.  Records the FULL realized per-pixel dose profile + max deviation +
    source-power distribution.  Returns a Materialization (flag OK) or None
    (peak / incident / count guard violated -> caller falls back to L0)."""
    rows = np.asarray(rows, dtype=np.float64)
    w = np.asarray(weights, dtype=np.float64)
    if rows.shape[0] == 0 or w.sum() <= 0:
        return None
    counts = _largest_remainder(M_rows * (w / w.sum()), M_rows)
    use = np.where(counts > 0)[0]
    if use.size == 0 or int(counts.sum()) != M_rows or (rows[use] < -1e-12).any():
        return None
    peak = float(rows[use].max())
    inc = float((counts * rows.sum(axis=1)).sum())
    if peak > peak_cap + 1e-9 or (incident_budget is not None
                                  and inc > incident_budget + 1e-6):
        return None
    dose = (counts[:, None] * rows).sum(0) / M_rows     # per-pixel dose (relaxed)
    dm = dose.mean()
    dev = float(np.abs(dose / dm - 1.0).max()) if dm > 0 else float("inf")
    power = rows[use].sum(axis=1)                        # per-row source power
    guards = {"exact_count": int(counts.sum()), "count_ok": True,
              "dose_dev": dev, "dose_band": None, "dose_relaxed": True,
              "dose_max": float(dose.max()), "dose_min": float(dose.min()),
              "dose_mean": float(dm),
              "power_min": float(power.min()), "power_max": float(power.max()),
              "power_mean": float(power.mean()), "power_std": float(power.std()),
              "incident_total": inc, "incident_budget": incident_budget,
              "incident_ok": bool(incident_budget is None or inc <= incident_budget + 1e-6),
              "peak": peak, "peak_ok": bool(peak <= peak_cap + 1e-9), "peak_cap": peak_cap,
              "admission_ok": True, "n_exchange": 0,
              "label": "NONDEPLOYABLE_DOSE_RELAXED_FW_DIAGNOSTIC"}
    expanded = np.repeat(np.arange(rows.shape[0]), counts.astype(np.int64))
    perm = np.random.default_rng(block_seed).permutation(expanded.size)
    acquired = rows[expanded[perm]]
    return Materialization(rows=acquired, counts=counts, atom_rows=rows,
                           atom_bankw=None, what=None, guards=guards, flag="OK",
                           n_exchange=0, gamma=np.ones(rows.shape[0]))


def _coverage_select(supports, n_pix, M_rows):
    """Coverage-aware greedy SELECTION of M_rows atoms (lifted from the drop rule
    of oed_design_v5.fixed_dose_scat32): repeatedly add the atom whose support
    covers the least-covered pixels (min sum of running coverage over its
    support), keeping per-pixel coverage as uniform as possible.  Deterministic
    (index tie-break).  supports = binary (A, n_pix)."""
    A = supports.shape[0]
    cov = np.zeros(n_pix)
    chosen = np.zeros(A, dtype=bool)
    order = np.arange(A)
    picks = []
    for _ in range(min(M_rows, A)):
        score = supports @ cov                          # sum coverage over support
        score[chosen] = np.inf
        j = int(order[np.argmin(score)])                # lowest -> index tie-break
        chosen[j] = True
        picks.append(j)
        cov += supports[j]
    return np.array(picks, dtype=np.int64)


def witness_materialize(rows, weights, M_rows=972, dose_band=0.05, w_lo=0.25,
                        w_hi=4.0, peak_cap=1536.0, incident_budget=None,
                        x_design=None, budget=0.60, block_seed=650777,
                        n_iter=6000, band_margin=0.003):
    """T-B admission-witness construction for an online k_eff>=32 design (R27 §2;
    lifted from oed_design_v5.fixed_dose_scat32 -- coverage-aware greedy selection
    + deterministic Sinkhorn nominal-amplitude balancing, the same machinery that
    produced the frozen library witnesses).  `rows` = unit-indicator geometry
    patterns (1.0 on each k_eff>=32 support), `weights` the design measure.

    Steps: (1) coverage-aware select M_rows atoms (uniform per-pixel coverage);
    (2) multiplicatively balance per-atom NOMINAL amplitudes c in [w_lo, w_hi]
    (default [1/4,4], the frozen bank power range) until per-pixel dose is within
    +-dose_band; (3) global scale so the mean operating load at x_design equals
    `budget` (dose_dev is scale-invariant); (4) verify peak<=peak_cap, incident
    budget, exact M_rows, c in range.  Seed-free deterministic (R27 §2.3).
    Returns a Materialization (gamma = c) or None (caller -> L0 fallback)."""
    rows = np.asarray(rows, dtype=np.float64)
    w = np.asarray(weights, dtype=np.float64)
    A, n = rows.shape
    if A == 0 or w.sum() <= 0:
        return None
    supports = (rows > 0).astype(np.float64)
    if A > M_rows:                                      # coverage-aware selection
        sel = _coverage_select(supports, n, M_rows)
    elif A == M_rows:
        sel = np.arange(A)
    else:                                               # too few atoms to admit
        return None
    R = supports[sel]                                   # unit-amplitude base rows
    counts = np.ones(sel.size, dtype=np.int64)          # one exposure each
    Rb = R
    deg = np.maximum(Rb.sum(axis=1), 1.0)
    c = np.ones(sel.size)
    for _ in range(n_iter):
        dose = c @ R
        dm = dose.mean()
        if dm <= 0:
            return None
        rel = dose / dm
        if np.abs(rel - 1.0).max() <= dose_band - band_margin:
            break
        c = np.clip(c / np.power((Rb @ rel) / deg, 0.7), w_lo, w_hi)
    rows_eff = c[:, None] * R
    # global load scaling: mean operating load at x_design == budget (scale-free
    # for dose_dev; sets absolute loads/incident/peak).
    if x_design is not None:
        loads = rows_eff @ np.asarray(x_design, dtype=np.float64)
        ml = float(loads.mean())
        if ml > 0:
            rows_eff = rows_eff * (budget / ml)
    dose = rows_eff.sum(0)
    dev = float(np.abs(dose / dose.mean() - 1.0).max()) if dose.mean() > 0 else np.inf
    peak = float(rows_eff.max()) if rows_eff.size else 0.0
    inc = float(rows_eff.sum())
    ok = (dev <= dose_band + 1e-9 and peak <= peak_cap + 1e-9
          and (incident_budget is None or inc <= incident_budget + 1e-6)
          and w_lo - 1e-9 <= c.min() and c.max() <= w_hi + 1e-9)
    guards = {"exact_count": int(counts.sum()), "count_ok": True,
              "dose_dev": dev, "dose_dev_pretrim": dev,
              "dose_ok": bool(dev <= dose_band + 1e-9), "dose_band": dose_band,
              "incident_total": inc, "incident_budget": incident_budget,
              "incident_ok": bool(incident_budget is None or inc <= incident_budget + 1e-6),
              "peak": peak, "peak_ok": bool(peak <= peak_cap + 1e-9), "peak_cap": peak_cap,
              "admission_ok": True, "n_exchange": 0, "gamma_min": float(c.min()),
              "gamma_max": float(c.max()), "gamma_std": float(c.std()),
              "from_witness_construct": True}
    if not ok:
        return None
    perm = np.random.default_rng(block_seed).permutation(sel.size)
    acquired = rows_eff[perm]
    gamma_full = np.ones(A); gamma_full[sel] = c
    return Materialization(rows=acquired, counts=counts, atom_rows=R,
                           atom_bankw=None, what=None, guards=guards, flag="OK",
                           n_exchange=0, gamma=gamma_full)


def _materialize_from_witness(banks, k, wm, M_rows, dose_band, incident_budget,
                              peak_cap, block_seed):
    """R27 §2 pure-bank realization from the frozen exact-972 admission witness
    (multiplicities aligned to the bank rows, nominal powers, gamma=1).  Guards
    are recomputed and re-verified post-hoc; returns None (fall through to the
    general materializer) if the witness is malformed or fails a guard."""
    K = len(banks)
    bk = banks[k]
    counts = np.asarray(wm, dtype=np.int64)
    rows = bk.rows
    if counts.shape[0] != rows.shape[0] or int(counts.sum()) != M_rows:
        return None
    dose = (counts[:, None] * rows).sum(0) / M_rows
    dmean = dose.mean()
    dev = float(np.abs(dose / dmean - 1.0).max()) if dmean > 0 else np.inf
    inc_row = bk.incident if bk.incident is not None else rows.sum(axis=1)
    inc_real = float((counts * inc_row).sum())
    peak_real = float(rows[counts > 0].max()) if (counts > 0).any() else 0.0
    dose_ok = dev <= dose_band + 1e-9
    inc_ok = (incident_budget is None) or (inc_real <= incident_budget + 1e-6)
    peak_ok = peak_real <= peak_cap + 1e-9
    guards = {"exact_count": int(counts.sum()), "count_ok": True,
              "dose_dev": dev, "dose_dev_pretrim": dev, "dose_ok": bool(dose_ok),
              "dose_band": dose_band, "incident_total": inc_real,
              "incident_budget": incident_budget, "incident_ok": bool(inc_ok),
              "peak": peak_real, "peak_ok": bool(peak_ok), "peak_cap": peak_cap,
              "admission_ok": True, "n_exchange": 0, "gamma_min": 1.0,
              "gamma_max": 1.0, "gamma_std": 0.0, "from_witness": True}
    if not (dose_ok and inc_ok and peak_ok):
        return None                                     # let the general path try
    expanded = np.repeat(np.arange(rows.shape[0]), counts)
    perm = np.random.default_rng(block_seed).permutation(expanded.size)
    acquired = rows[expanded[perm]]
    what = np.zeros(K); what[k] = 1.0
    return Materialization(rows=acquired, counts=counts, atom_rows=rows,
                           atom_bankw=None, what=what, guards=guards, flag="OK",
                           n_exchange=0, gamma=np.ones(rows.shape[0]))


def _mat_fail(banks, K, flag):
    """L0 realization fallback (R25 SS9): return the frozen knob bank's exact
    rows.  The knob bank L0 is materialized by taking its atoms at their integer
    expected counts (it is already a compliant exact design)."""
    knob = next((b for b in banks if b.is_knob), banks[0])
    M = 972
    tgt = M * knob.weights
    cnts = _largest_remainder(tgt, M)
    expanded = np.repeat(np.arange(knob.rows.shape[0]), cnts)
    acquired = knob.rows[expanded][:M]
    what = np.zeros(K); what[0] = 1.0
    return Materialization(rows=acquired, counts=cnts, atom_rows=knob.rows,
                           atom_bankw=None, what=what,
                           guards={"fallback": "L0"}, flag=flag, n_exchange=0)


# ============================================================================ #
#  11.  post-materialization: realized regret + KKT recompute (R25.2, SS6, SS9)  #
# ============================================================================ #

def realized_regret(mat, scenarios, B, nu, c, sm, side, lambda_TV, eps,
                    xhat=None, freeze_directions=True, prescan_rows=None,
                    rho_prescan=0.60, tv_smooth=1e-3, Dx=None, Dy=None):
    """t_real = max_s R_s(X_s^real)/R*_s where X_s^real = H0_s + sum over the 972
    materialized rows of H_s(row) (R25 SS9 post-materialization).

    Direction convention matches the allocation (freeze_directions): with frozen
    directions q = a/(a.xhat) the realized F uses the SAME projected atom
    geometry as R*_s and F_{k,s}, so t_real - t_cont is a pure rounding/exchange
    measure (and the L0 fallback certifies t_real >= 1).  With exact per-scene
    directions q = a/(a.x_s), t_real is the fully honest realized risk."""
    W = sm.W; Rstar = sm.Rstar
    S = scenarios.draws.shape[0]
    rows = mat.rows
    if freeze_directions and xhat is not None:
        l0 = np.maximum(rows @ xhat, 1e-12)
        Gfroz = (rows / l0[:, None]) @ B                 # frozen projected dirs
    r_real = np.empty(S)
    for s in range(S):
        x = scenarios.draws[s]
        loads = np.maximum(rows @ x, 1e-12)
        if freeze_directions and xhat is not None:
            Q = Gfroz
        else:
            Q = (rows / loads[:, None]) @ B
        Freal = (Q.T * (nu * kernel_J(loads, c))) @ Q
        X = sm.H0[s] + Freal
        R, _ = arisk(X, W)
        r_real[s] = R / Rstar[s]
    return float(r_real.max()), r_real


# ============================================================================ #
#  12.  top-level pipeline + disclosures (R25 SS13)                              #
# ============================================================================ #

@dataclass
class RLMIConfig:
    n_scenarios: int = 16
    nu: float = 2000.0
    c: float = 0.05                 # hidden hold-time CV (hard decision cell)
    lambda_TV: float = 1.0
    eps: float = 1e-6               # H0 ridge
    tv_smooth: float = 1e-3
    rho_prescan: float = 0.60
    M_rows: int = 972
    dose_band: float = 0.05
    incident_budget: Optional[float] = None
    peak_cap: float = 1536.0
    side: int = 32
    seed_scenario: int = 650000
    seed_block_order: int = 650777
    jitter_scale: float = 1.0
    oracle_tol: float = 1e-10
    maximin_tol: float = 1e-9
    tie_tol: float = 1e-6
    knob_index: int = 0
    freeze_directions: bool = True
    compute_loo: bool = False        # R25 SS3 Delta_LOO (S extra solves; off=fast)
    fast_linalg: bool = False        # R27 SS3.2 batched maximin (output-equivalent)
    use_cache: bool = False          # R27 SS3.2 scenario/oracle cache (per state)


# ---- R27 SS3.2 scenario-matrix + oracle cache (output-equivalent) ---------- #
#  Keyed by the identical (xhat, banks, B, W, physics) state; stores the
#  ScenarioSet + solved ScenarioMatrices so repeated allocator calls (equivalence
#  checks, delta_LOO, arms sharing the pre-scan state) skip make_scenarios +
#  build_scenario_matrices + compute_oracles.  On a HIT the returned objects are
#  bit-identical -> the maximin trajectory and 972-row realization are identical.
#  For the bridge GRID each cell is a distinct state (cache miss) -> no per-cell
#  route-latency change; disabled by default for that reason.
_SM_CACHE = {}
_SM_CACHE_MAX = 8


def _state_key(xhat, banks, B, W, config):
    h = hashlib.sha256()
    h.update(np.ascontiguousarray(xhat, dtype=np.float64).tobytes())
    h.update(np.ascontiguousarray(B, dtype=np.float64).tobytes())
    h.update(np.ascontiguousarray(W, dtype=np.float64).tobytes())
    for bk in banks:
        h.update(bk.name.encode())
        h.update(np.ascontiguousarray(bk.rows, dtype=np.float64).tobytes())
        h.update(np.ascontiguousarray(bk.weights, dtype=np.float64).tobytes())
    for fld in ("n_scenarios", "nu", "c", "lambda_TV", "eps", "tv_smooth",
                "rho_prescan", "M_rows", "side", "seed_scenario", "jitter_scale",
                "oracle_tol", "freeze_directions"):
        h.update(repr(getattr(config, fld)).encode())
    return h.hexdigest()


def run_rlmi(xhat, banks, B, H0_full, config: RLMIConfig, W=None,
             prescan_rows=None, Dx=None, Dy=None):
    """Full R25 Step-3 pipeline.  Returns a disclosures dict (R25 SS13).

    Inputs
      xhat     : (n,) pre-scan estimate in [0,1].
      banks    : list[Bank] (banks[knob_index].is_knob == True; the L0 corner).
      B        : (n, r) declared reconstruction subspace (orthonormal columns).
      H0_full  : (n, n) full-space pre-scan+prior curvature for the Laplace cov.
      config   : RLMIConfig.
    """
    import time
    t0 = time.perf_counter()
    flag = "OK"
    r = B.shape[1]
    if W is None:
        W = np.eye(r)
    side = config.side

    if Dx is None or Dy is None:
        Dx, Dy = _grad_ops(side)

    # -- 1-2. scenarios + scenario matrices + oracles (R27 SS3.2 cacheable) -- #
    ck = _state_key(xhat, banks, B, W, config) if config.use_cache else None
    if ck is not None and ck in _SM_CACHE:
        scen, sm = _SM_CACHE[ck]                       # bit-identical reuse
    else:
        scen = make_scenarios(xhat, H0_full, n_scenarios=config.n_scenarios,
                              seed_base=config.seed_scenario,
                              jitter_scale=config.jitter_scale)
        sm = build_scenario_matrices(scen, banks, B, config.nu, config.c,
                                     config.lambda_TV, config.eps, side,
                                     prescan_rows=prescan_rows,
                                     rho_prescan=config.rho_prescan,
                                     tv_smooth=config.tv_smooth, W=W,
                                     M_rows=config.M_rows, Dx=Dx, Dy=Dy,
                                     freeze_directions=config.freeze_directions)
        try:
            sm = compute_oracles(sm, tol=config.oracle_tol)
        except Exception as e:               # noqa
            return _solver_fail(banks, config, "SOLVER_FAIL", str(e), scen, sm)
        if ck is not None:
            if len(_SM_CACHE) >= _SM_CACHE_MAX:
                _SM_CACHE.pop(next(iter(_SM_CACHE)))
            _SM_CACHE[ck] = (scen, sm)

    # -- 3-4. standardized maximin + tie-break ------------------------------ #
    try:
        mm = standardized_maximin(sm, tol=config.maximin_tol,
                                  tie_tol=config.tie_tol,
                                  fast_linalg=config.fast_linalg)
        w_star = tie_break_to_knob(sm, mm, knob_index=config.knob_index,
                                   tie_tol=config.tie_tol)
    except Exception as e:               # noqa
        return _solver_fail(banks, config, "SOLVER_FAIL", str(e), scen, sm)

    t_alloc = time.perf_counter() - t0

    # -- R23 residual disclosures (branch-fenced) --------------------------- #
    knob = banks[config.knob_index]
    r23 = [r23_residual(knob, scen.draws[s], B, config.nu, config.c,
                        eps=config.eps) for s in range(config.n_scenarios)]
    on_branch_all = all(d["on_branch"] for d in r23)
    if on_branch_all:
        U1 = max(d["U1_bound"] for d in r23)
        U1_status = "ON_BRANCH"
    else:
        U1 = None
        U1_status = "BRANCH_VIOLATION"

    # -- 5. materialization ------------------------------------------------- #
    mat = materialize(banks, w_star, xhat, M_rows=config.M_rows,
                      dose_band=config.dose_band,
                      incident_budget=config.incident_budget,
                      peak_cap=config.peak_cap, side=side,
                      block_seed=config.seed_block_order)
    if mat.flag != "OK":
        flag = mat.flag

    # -- 6. realized regret + KKT recompute --------------------------------- #
    t_real, r_real = realized_regret(mat, scen, B, config.nu, config.c, sm, side,
                                     config.lambda_TV, config.eps, xhat=xhat,
                                     freeze_directions=config.freeze_directions,
                                     prescan_rows=prescan_rows,
                                     rho_prescan=config.rho_prescan,
                                     tv_smooth=config.tv_smooth, Dx=Dx, Dy=Dy)

    dLOO = None
    if config.compute_loo:
        try:
            dLOO, _ = delta_loo(sm, w_star, config)
        except Exception:               # noqa
            dLOO = None

    wall = time.perf_counter() - t0
    what = mat.what
    A_cont = 1.0 - float(w_star[config.knob_index])
    A_real = 1.0 - float(what[config.knob_index])
    ent = _entropy(what)

    disc = {
        "flag": flag,
        "w_star": w_star.tolist(),
        "w_realized": what.tolist(),
        "alpha_s": mm.alpha.tolist(),
        "t_cont": mm.t_cont,
        "t_real": t_real,
        "materialization_regret": t_real - mm.t_cont,
        "A_cont": A_cont,
        "A_realized": A_real,                    # R25.4  A_j = 1 - what_0
        "entropy_realized": ent,                 # reported, NOT gated (R25 SS11)
        "active_banks_cont": mm.support.tolist(),
        "active_banks_realized": np.where(what > 1e-9)[0].tolist(),
        "active_scenarios": mm.active_scenarios.tolist(),
        "lambda_kkt": mm.lam,
        "kkt_residual_cont": mm.kkt_residual,
        "r_s": mm.r_s.tolist(),
        "r_s_realized": r_real.tolist(),
        "R_star": sm.Rstar.tolist(),
        "branch": sm.branch,
        "eps1_m_U1": {"per_scenario": r23, "U1": U1, "status": U1_status},
        "guards": mat.guards,
        "n_exchange": mat.n_exchange,
        "wall_time_s": wall,
        "alloc_time_s": t_alloc,                 # route latency (Gate D; ex-recon)
        "delta_LOO": dLOO,                       # R25 SS3 (None unless compute_loo)
        "config": asdict(config),
    }
    return disc, mat, sm, scen, mm


def delta_loo(sm, w_star, config):
    """R25 SS3 leave-one-draw allocation spread  Delta_LOO = max_s ||w* - w*_{-s}||_1.
    Re-solves the maximin S times, each dropping one scenario (oracles R*_s reused).
    Off by default in run_rlmi (S extra solves)."""
    S = len(sm.H0)
    spreads = []
    for s in range(S):
        keep = [i for i in range(S) if i != s]
        sub = ScenarioMatrices(H0=[sm.H0[i] for i in keep],
                               F=[sm.F[i] for i in keep], Rstar=sm.Rstar[keep],
                               vstar=None, W=sm.W, branch={}, kernel_c=sm.kernel_c,
                               nu=sm.nu, r=sm.r)
        mm_s = standardized_maximin(sub, tol=config.maximin_tol,
                                    tie_tol=config.tie_tol)
        w_s = tie_break_to_knob(sub, mm_s, knob_index=config.knob_index,
                                tie_tol=config.tie_tol)
        spreads.append(float(np.abs(w_star - w_s).sum()))
    return float(max(spreads)), spreads


def _entropy(w):
    w = np.asarray(w, dtype=np.float64)
    p = w[w > 0]
    if p.size <= 1:
        return 0.0
    return float(-(p * np.log(p)).sum() / np.log(len(w)))


def _solver_fail(banks, config, flag, msg, scen=None, sm=None):
    mat = _mat_fail(banks, len(banks), flag)
    disc = {"flag": flag, "error": msg,
            "w_star": (np.eye(len(banks))[config.knob_index]).tolist(),
            "w_realized": mat.what.tolist(), "guards": mat.guards}
    return disc, mat, sm, scen, None


# ============================================================================ #
#  synthetic bank builder (unit-test fixtures; NOT the T-B banks)               #
# ============================================================================ #

def synth_bank(name, rows, weights=None, is_knob=False, admission=None,
               peak=None, incident=None, xhat=None):
    rows = np.asarray(rows, dtype=np.float64)
    R = rows.shape[0]
    if weights is None:
        weights = np.full(R, 1.0 / R)
    load = rows @ xhat if xhat is not None else rows.sum(axis=1)
    return Bank(name=name, rows=rows, weights=np.asarray(weights, float),
                load=load, incident=incident, peak=peak, admission=admission,
                is_knob=is_knob)


if __name__ == "__main__":
    print("rlmi.py -- RLMI allocator + constrained materializer (R25).")
    for cc in (0.0, 0.05, 0.1):
        print("  rho_c(c=%.2f) = %s" % (cc, kernel_rho_c(cc)))
