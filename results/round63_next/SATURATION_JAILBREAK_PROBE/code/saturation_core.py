"""
SATURATION-CODED BUCKET DETECTION -- core physics.

Working name: "saturation jailbreak". A SiPM bucket detector is C binary
microcells; the per-cell (Poisson >=1 photon) saturation acts BEFORE the
spatial sum, so the swept-power response encodes NONLINEAR functionals of the
masked scene v = m (x) x, i.e. power sums p_k = sum_p v_p^k, not just the
linear bucket p_1.

============================================================================
NORMALIZATION (stated once, used everywhere -- this is where silent errors live)
============================================================================
Scene x in [0,1]^P (P = 1024 for 32x32).  Binary mask m in {0,1}^P.
Masked scene v = m (elementwise) x  in [0,1]^P.

Per-cell speckle weight vector w_c in R_+^P.  BASE MODEL (fully developed
speckle, one grain/cell): w_cp i.i.d. Exponential with MEAN 1/C, across cells
c=1..C AND pixels p.  (Gamma(k) generalisation below.)

Per-cell mean photon count (occupancy)         u_c = rho * w_c^T v .
Per-cell fire probability (Poisson, >=1 photon) f_c = 1 - exp(-u_c) .
Bucket output S = sum_c Bernoulli(f_c)  in {0,..,C}   (Poisson-binomial).

EFFECTIVE COUPLING  a := rho * E[w_cp] = rho / C .   (physical rho = a * C.)
We parametrise every sweep by `a`.  Mean per-cell occupancy ~ a * p_1.

ENSEMBLE (annealed / self-averaged) SURVIVAL:
    Q(a; v) = E_w[ exp(-rho w^T v) ] = prod_p 1/(1 + a v_p) .
Measured non-fire fraction 1 - E[S]/C -> Q(a; v) as C->inf (LLN over cells).

KEY IDENTITY (log / power-sum form -- clean, well-conditioned):
    L(a) := -log Q(a; v) = sum_p log(1 + a v_p)
          = sum_{k>=1} (-1)^{k+1} a^k p_k / k
          = a p_1 - a^2 p_2 / 2 + a^3 p_3 / 3 - ...
    with power sums p_k = sum_p v_p^k .
The ordinary linear bucket = p_1.  The saturation sweep also gives
p_2 = sum v_p^2 (and higher) -> a genuine null-space break.

GAMMA(k) GENERALISATION (grain/cell mismatch): w ~ Gamma(shape k, mean 1/C):
    Q_k(a; v) = prod_p (1 + a v_p / k)^{-k},
    L_k(a)    = a p_1 - a^2 p_2 / (2k) + ...     (p_2 term scales as 1/k;
    k=1 => Exponential; k->inf => linear detector, signal vanishes.)

VARIANCE (Poisson-binomial, self-averaged):
    Var(S)/C -> E_w[f(1-f)] = Q(a) - Q(2a) =: Q - Q2 .
"""
import numpy as np


# --------------------------------------------------------------------------- #
# Ensemble (annealed) closed forms
# --------------------------------------------------------------------------- #
def Q_exp(a, v):
    """Ensemble survival Q(a;v) = prod_p 1/(1+a v_p), exponential speckle (k=1).

    a may be a scalar or array (broadcast over leading axis); v is 1-D."""
    a = np.asarray(a, float)
    v = np.asarray(v, float)
    logQ = -np.log1p(a[..., None] * v).sum(axis=-1)
    return np.exp(logQ)


def Q_gamma(a, v, k):
    """Ensemble survival for Gamma(shape k, mean 1/C) speckle weights."""
    a = np.asarray(a, float)
    v = np.asarray(v, float)
    logQ = (-k * np.log1p(a[..., None] * v / k)).sum(axis=-1)
    return np.exp(logQ)


def power_sums(v, kmax):
    """p_k = sum_p v_p^k for k=1..kmax  (index 0 is NaN placeholder)."""
    v = np.asarray(v, float)
    return np.array([np.nan] + [np.sum(v ** k) for k in range(1, kmax + 1)])


def L_series_pk(a, pks, k_gamma=1):
    """L(a) via truncated power-sum series (cross-check only).
    pks = [nan, p1, p2, ...].  L = sum_k (-1)^{k+1} a^k p_k /(k*kg^{k-1})."""
    a = np.asarray(a, float)
    out = np.zeros_like(a, dtype=float)
    for k in range(1, len(pks)):
        out += ((-1.0) ** (k + 1)) * a ** k * pks[k] / (k * k_gamma ** (k - 1))
    return out


# --------------------------------------------------------------------------- #
# Monte-Carlo simulators
# --------------------------------------------------------------------------- #
def draw_W(C, P, rng, k_gamma=1):
    """Draw a speckle weight matrix W (C x P). MEAN of each entry = 1/C.
    Exponential (k_gamma=1) or Gamma(shape=k_gamma, mean 1/C)."""
    if k_gamma == 1:
        return rng.exponential(scale=1.0 / C, size=(C, P))
    return rng.gamma(shape=k_gamma, scale=1.0 / (C * k_gamma), size=(C, P))


def empirical_survival(a, v, W):
    """QUENCHED empirical non-fire fraction (1/C) sum_c exp(-rho w_c^T v) for a
    FIXED realised W.  rho = a*C.  Returns array shaped like a."""
    a = np.atleast_1d(np.asarray(a, float))
    C = W.shape[0]
    rho = a * C
    Uc = W @ v                                        # (C,)
    surv = np.exp(-rho[:, None] * Uc[None, :]).mean(axis=1)
    return surv


def cell_fire_prob(a, v, W):
    """f_c(a) fire prob per cell for FIXED W. Returns (len(a), C)."""
    a = np.atleast_1d(np.asarray(a, float))
    C = W.shape[0]
    rho = a * C
    Uc = W @ v
    return 1.0 - np.exp(-rho[:, None] * Uc[None, :])


def simulate_bucket(a, v, W, M, rng, fast=False):
    """Full noisy measurement. For FIXED W, each of M windows: cell c fires
    ~Bernoulli(f_c).  Returns S (len(a), M) bucket counts and f (len(a), C).

    fast=True: draw S from the Gaussian approximation to the Poisson-binomial,
    S ~ Normal(sum f_c, sum f_c(1-f_c)) (exact mean & variance; valid for large
    C, and consistent with the CRB's Gaussian model).  Much faster for large M."""
    a = np.atleast_1d(np.asarray(a, float))
    C = W.shape[0]
    f = cell_fire_prob(a, v, W)
    if fast:
        mu = f.sum(axis=1)                         # (L,)
        var = (f * (1 - f)).sum(axis=1)            # (L,)
        S = rng.normal(mu[:, None], np.sqrt(var)[:, None], size=(len(a), M))
        S = np.clip(np.round(S), 0, C).astype(np.int64)
        return S, f
    S = np.empty((len(a), M), dtype=np.int64)
    for li in range(len(a)):
        S[li] = (rng.random((M, C)) < f[li][None, :]).sum(axis=1)
    return S, f


# --------------------------------------------------------------------------- #
# Estimation: recover power sums from a survival sweep
# --------------------------------------------------------------------------- #
def design_matrix(a_levels, order, k_gamma=1):
    """A[:,k-1] = (-1)^{k+1} a^k /(k * kg^{k-1}) so that L = A @ [p1..p_order]."""
    a_levels = np.asarray(a_levels, float)
    cols = [((-1.0) ** (k + 1)) * a_levels ** k / (k * k_gamma ** (k - 1))
            for k in range(1, order + 1)]
    return np.stack(cols, axis=1)


def fit_powersums(a_levels, Q_meas, order=3, k_gamma=1, weights=None):
    """Recover (p1,p2,...) by least squares of L=-log Q vs truncated series."""
    Q_meas = np.asarray(Q_meas, float)
    L = -np.log(Q_meas)
    A = design_matrix(a_levels, order, k_gamma)
    if weights is None:
        coef, *_ = np.linalg.lstsq(A, L, rcond=None)
    else:
        Wd = np.sqrt(weights)
        coef, *_ = np.linalg.lstsq(A * Wd[:, None], L * Wd, rcond=None)
    p = {f"p{k}": float(coef[k - 1]) for k in range(1, order + 1)}
    p["_resid"] = float(np.sqrt(np.mean((A @ coef - L) ** 2)))
    p["_coef"] = coef
    return p


def default_a_levels(v, L=10, occ_lo=0.05, occ_hi=3.0):
    """Log-spaced sweep of `a` so mean per-cell occupancy a*p1 in [occ_lo,occ_hi]."""
    p1 = float(np.sum(v))
    return np.logspace(np.log10(occ_lo / p1), np.log10(occ_hi / p1), L)


def fit_es_poly(a_levels, Q_meas, order, weights=None):
    """EXACT elementary-symmetric estimator.  Fit
        1/Q(a) = 1 + a e_1 + a^2 e_2 + ... + a^order e_order
    (e_0 = 1 fixed).  For a scene with <= `order` active pixels this is EXACT
    (1/Q is a degree-P polynomial in a).  Returns dict e1..e_order and derived
    p1,p2 via Newton (p1=e1, p2=e1^2-2 e2)."""
    a_levels = np.asarray(a_levels, float)
    y = 1.0 / np.asarray(Q_meas, float) - 1.0
    A = np.stack([a_levels ** k for k in range(1, order + 1)], axis=1)
    if weights is None:
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    else:
        Wd = np.sqrt(weights)
        coef, *_ = np.linalg.lstsq(A * Wd[:, None], y * Wd, rcond=None)
    e = {f"e{k}": float(coef[k - 1]) for k in range(1, order + 1)}
    e["p1"] = e["e1"]
    e["p2"] = e["e1"] ** 2 - 2.0 * e["e2"] if order >= 2 else np.nan
    e["_coef"] = coef
    return e


# --------------------------------------------------------------------------- #
# Noise / variance / resource closed forms
# --------------------------------------------------------------------------- #
def var_S_over_C(a, v):
    """Self-averaged Var(S)/C = E_w[f(1-f)] = Q(a) - Q(2a)."""
    return Q_exp(a, v) - Q_exp(2.0 * a, v)


def incident_photons(a, v, C):
    """Expected incident photons per window on whole detector = rho*p1 = a*C*p1.
    (RESOURCE, includes photons wasted on already-fired cells.)"""
    return a * C * float(np.sum(v))


# --------------------------------------------------------------------------- #
# R31-ruling additions (MOLT): coherent counterexample, exact covariance,
# annealed FIM in power-sum coords, three-level design.
# NOTE: throughout, the sweep variable `a` here == the ruling's rho (W mean 1);
# dimensionless load t = a * p1; incident photons/gate = a*C*p1 = t*C.
# --------------------------------------------------------------------------- #
def coherent_survival_mc(a, v, C, rng):
    """MC no-fire fraction under the COHERENT-collapse model (ruling S1.3):
    E_c = sum_p h_cp sqrt(v_p), h~CN(0,1) -> T_c=|E_c|^2 ~ Exp(p1).  The curve
    depends only on p1: r(a)=1/(1+a p1).  Returns array over `a`."""
    a = np.atleast_1d(np.asarray(a, float))
    v = np.asarray(v, float)
    Nn = v.size
    h = (rng.standard_normal((C, Nn)) + 1j * rng.standard_normal((C, Nn))) / np.sqrt(2)
    E = h @ np.sqrt(v)
    T = np.abs(E) ** 2                                  # (C,)  ~Exp(p1)
    return np.array([np.exp(-ai * T).mean() for ai in a])


def coherent_survival_exact(a, v):
    """r(a)=1/(1+a p1) -- coherent-collapse closed form (p2-independent)."""
    a = np.atleast_1d(np.asarray(a, float))
    return 1.0 / (1.0 + a * float(np.sum(v)))


def r_cov_exact(ai, aj, v, C, kind="quenched"):
    """Exact finite-C covariance of r_C (ruling S2.1):
       Cov[r_C(ai), r_C(aj)] = (r(ai+aj) - r(ai) r(aj)) / C.
    With kind='pergate_var' returns the ADDITIONAL per-gate binomial variance
    term diag piece (r_j - r(2 a_j))/C (per single gate G=1)."""
    r_i = Q_exp(np.array([ai]), v)[0]
    r_j = Q_exp(np.array([aj]), v)[0]
    r_ij = Q_exp(np.array([ai + aj]), v)[0]
    return (r_ij - r_i * r_j) / C


def annealed_fim(a_levels, G_levels, v, C, order):
    """Annealed binomial FIM in power-sum coords (ruling S3.2):
        I_{lm} = C sum_j G_j (r_j/(1-r_j)) (-1)^{l+m} a_j^{l+m}/(l m).
    Returns (order x order) matrix for params (p1,...,p_order)."""
    v = np.asarray(v, float); v = v[v > 0]
    a_levels = np.asarray(a_levels, float)
    G_levels = np.asarray(G_levels, float)
    r = Q_exp(a_levels, v)
    I = np.zeros((order, order))
    for j in range(len(a_levels)):
        dlog = np.array([((-1.0) ** l) * a_levels[j] ** l / l
                         for l in range(1, order + 1)])
        I += C * G_levels[j] * (r[j] / (1 - r[j])) * np.outer(dlog, dlog)
    return I


def three_level_design(p1, t_anchor=0.1, t_audit=1.0, t_ridge=3.54,
                       alloc=(0.245, 0.0, 0.755), budget_incident=None, C=3600):
    """R31 robust three-level sweep in load t (ruling S4.3).  Returns
    (a_levels, G_levels).  If budget_incident given, allocate incident photons
    by `alloc` fractions and convert to integer gate counts G_j = frac*B/(t_j*C).
    The audit level gets a small floor even at alloc 0."""
    ts = np.array([t_anchor, t_audit, t_ridge])
    a_levels = ts / p1
    if budget_incident is None:
        return a_levels, None
    alloc = np.array(alloc, float)
    if alloc[1] == 0:
        alloc = np.array([alloc[0], 0.03, alloc[2]])       # small audit floor
        alloc = alloc / alloc.sum()
    G = np.maximum(1, np.round(alloc * budget_incident / (ts * C))).astype(int)
    return a_levels, G
