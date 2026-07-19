"""R16 estimator diagnostic: single-load conditional-score / missing-information
Fisher information for the jittered nonparalyzable renewal channel.

This implements the diagnostic prescribed verbatim by GPT ruling R16
(`docs/ROUND63_GPT_ROUND16_RULING_RAW.md`, section 5.2). It exists to decide
whether the replicated ~-6% end-point low bias in the jitter-capped ridge-law
peaks (`docs/R14_PREREGISTERED_PREDICTIONS.md`) is an artifact of the histogram
log-pmf / finite-difference Fisher estimator, or a genuine physical residual.

It is CONSTANT-FREE by design: it fits nothing and proposes no new ridge-law
constant. Per R16 section 5.4, no new constant may be fitted before this
diagnostic rules.

================================================================================
Channel (identical generation to jitter_scalar_fi_v2.py, so draws are
comparable draw-for-draw)
================================================================================
Active start; wait W_i ~ Exp(lam); detect (count += 1); dead-hold B_i
(deterministic tau for cv=0, else mean-tau lognormal with the given cv); repeat.
With A_m = sum_{i<=m} unit-exponentials (so W_i = A_i/lam) and D_m = sum_{i<m}B_i,
the m-th detection lands at S_m = A_m/lam + D_m, and it fits the window iff
A_m <= lam*(T - D_m). Counts for several lam therefore come from ONE (A,D) draw
(common random numbers / exact pathwise coupling in lam). Here T = nu*tau, and
rho = lam*tau is the dimensionless load; with tau=1 (used throughout) lam = rho.

================================================================================
The R16 identity (section 5.2) -- no finite difference, no histogram, one load
================================================================================
The simulator knows each frame's total ACTIVE (live) time L. The Poisson
intensity is lam during live time, so lam*L is the compensator of N, and the
complete-path score with respect to theta = log rho is

    U = N - lam*L                                       (tau=1: U = N - rho*L).

Because M = N - lam*L is the compensated counting martingale, its predictable
quadratic variation is the compensator, giving the exact complete-data
information

    I_complete = Var(U) = E[N]  = lam*E[L].

The observed data is the count N only. By the missing-information / total-variance
decomposition (E[U] = 0 at the truth),

    I_complete = Var(E[U|N]) + E[Var(U|N)],
    I_N        = Var(E[U|N]) = I_complete - E[Var(U|N)].

Since U is affine in L within a fixed-N bin, Var(U|N) = lam^2 Var(L|N), hence the
R16 count-only Fisher information

    I_N = E[N] - lam^2 * E[ Var(L|N) ].                 (*)

Monte-Carlo estimators (single load; bins n = distinct count values):
  Missing-information form (R16 5.2):
    I_MI = mean(N) - lam^2 * sum_{n: k_n>=2} (k_n/N_MC) * s^2_{L|n}
  Bias-corrected conditional-score cross-check (R16 5.2):
    I_CS = sum_{n: k_n>=2} (k_n/N_MC) * ( Ubar_n^2 - s^2_{U|n}/k_n )
  s^2 is the unbiased within-bin sample variance. Per R16, the two estimators
  "must agree within Monte-Carlo uncertainty". Per-slot J = I / nu.

L in the simulation variables (nonparalyzable, boundary handled exactly):
  For N>=1, active time = sum of the N pre-detection waits (= A_N/lam) plus, if
  the window ends AFTER the final dead-hold, the residual live tail:
    L = A_N/lam + max(0, T - A_N/lam - D_{N+1}),  D_{N+1} = B_1+...+B_N.
  For N=0 the detector is live the whole window (L = T); the same expression
  reproduces this with A_0/lam = 0 and D_1 = 0.

================================================================================
--toy : exactly-known Fisher information (closed form, derived here)
================================================================================
Single-count regime. Take the window no longer than one dead-hold (tau >= T) so
that after any detection the detector is dead through the end of the window,
forcing N in {0,1}. Concretely nu=1, tau=1, cv=0  =>  T=1, lam=rho, and

    p1 = P(N=1) = P(W_1 <= T) = 1 - e^{-rho T}.

This is a Bernoulli(p1) family in theta = log rho. With d p1/d theta = rho*T*e^{-rho T},

    I_N(theta) = (dp1/dtheta)^2 (1/p1 + 1/p0)
               = (rho T e^{-rho T})^2 * [1/(1-e^{-rho T}) + 1/e^{-rho T}]
               = T^2 * rho^2 * e^{-rho T} / (1 - e^{-rho T}).           (closed form)

This is exact for every rho and reduces (T=1) to  I_toy = rho^2 e^{-rho}/(1-e^{-rho}).
It is a genuine test of the full identity (*): L is deterministic (=T) when N=0
but random when N=1 (L|N=1 = W_1 truncated to [0,T], with strictly positive
variance), so the lam^2*E[Var(L|N)] subtraction is exercised, not bypassed.
Consistency check of (*): E[N] = p1; Var(L|N=1) = Var(W_1|W_1<=T) =
(1/rho^2)[1 - (rho T)^2 e^{-rho T}/(1-e^{-rho T})^2]; substituting reproduces the
boxed closed form. Both routes are hard-coded and asserted equal below.

(A second, weaker exact anchor -- zero dead-hold, B_i=0 => L=T deterministically,
N ~ Poisson(rho*nu), I_N = rho*nu, J = rho -- is available via --toy-poisson but
does NOT exercise the Var(L|N) subtraction, so the single-count regime is the
primary toy.)

Usage:
  python jitter_score_diag.py --toy [--nmc 20000]
  python jitter_score_diag.py --toy-poisson [--nmc 20000]
  python jitter_score_diag.py --compare --nu 200 --cv 0.05 --rho 5.65 [--nmc 20000]

All local runs are small CPU smoke tests (< 2 min, N_MC <= 20k by contract);
the full 2e6-frame diagnostic runs on Colab via jitter_score_diag_colab.py.
"""
import argparse
import sys
import time

import numpy as np

TAU = 1.0
DLOG = 0.01          # finite-difference step of the histogram estimator being audited
HIST_ALPHA = 0.5     # pseudocount numerator of the audited floor (alpha/N), v2 default
TAIL_THRESH = 50     # R16: report total probability of bins with count < 50


# ----------------------------------------------------------------------------
# Simulation core
# ----------------------------------------------------------------------------
def _m_max(nu, rho, dlog):
    """Column budget: safely above the largest attainable count (Gaussian-tight)."""
    mean_ct = nu * rho / (1.0 + rho) * (1.0 + 4.0 * dlog)
    return int(mean_ct + 12.0 * np.sqrt(max(nu, 1.0)) + 50.0)


def _active_time(A, D, lam, N, T, m_max):
    """Exact per-draw live time L at rate lam (see module docstring)."""
    nb = A.shape[0]
    ar = np.arange(nb)
    prev = np.maximum(N - 1, 0)
    A_N = np.where(N >= 1, A[ar, prev], 0.0)          # A_N/lam = sum of N waits
    D_Np1 = D[ar, np.minimum(N, m_max - 1)]           # D_{N+1} = B_1+...+B_N (0 for N=0)
    waits = A_N / lam
    tail = np.maximum(0.0, T - waits - D_Np1)
    return waits + tail                               # N=0 -> 0 + T = T


def simulate(rho, cv, nu, n_mc, seed, want_hist, dlog=DLOG, chunk=4000):
    """Return per-draw arrays from ONE set of (A,B) draws (common random numbers).

    Always returns N0 (count at lam0) and L0 (live time at lam0). When
    want_hist, also returns cm, cp (counts at lam0*exp(-+dlog)) so the audited
    histogram estimator can be run on the *identical* draws draw-for-draw.
    """
    lam0 = rho / TAU
    T = nu * TAU
    m_max = _m_max(nu, rho, dlog)
    rng = np.random.default_rng(seed)
    N0 = np.empty(n_mc, dtype=np.int64)
    L0 = np.empty(n_mc, dtype=np.float64)
    cm = np.empty(n_mc, dtype=np.int64) if want_hist else None
    cp = np.empty(n_mc, dtype=np.int64) if want_hist else None
    if cv > 0:
        sig2 = np.log1p(cv * cv)
        mu = -0.5 * sig2
    lam_m = lam0 * np.exp(-dlog)
    lam_p = lam0 * np.exp(+dlog)
    pos = 0
    while pos < n_mc:
        nb = min(chunk, n_mc - pos)
        A = np.cumsum(rng.exponential(1.0, size=(nb, m_max)), axis=1)
        if cv > 0:
            B = TAU * np.exp(mu + np.sqrt(sig2) * rng.standard_normal((nb, m_max)))
        else:
            B = np.full((nb, m_max), TAU)
        D = np.concatenate([np.zeros((nb, 1)), np.cumsum(B, axis=1)[:, :-1]], axis=1)
        slack = T - D                                  # lam-independent
        n0 = (A <= lam0 * slack).sum(axis=1)
        if int(n0.max()) >= m_max:
            raise RuntimeError(f"count saturated m_max={m_max} at rho={rho}, cv={cv}; "
                               "widen _m_max margin")
        N0[pos:pos + nb] = n0
        L0[pos:pos + nb] = _active_time(A, D, lam0, n0, T, m_max)
        if want_hist:
            cm[pos:pos + nb] = (A <= lam_m * slack).sum(axis=1)
            cp[pos:pos + nb] = (A <= lam_p * slack).sum(axis=1)
        pos += nb
    return N0, L0, cm, cp


# ----------------------------------------------------------------------------
# Estimators
# ----------------------------------------------------------------------------
def _binstats(N, x):
    """Unbiased within-bin means/variances of x grouped by integer count N."""
    vals, inv, counts = np.unique(N, return_inverse=True, return_counts=True)
    s1 = np.bincount(inv, weights=x)
    s2 = np.bincount(inv, weights=x * x)
    mean = s1 / counts
    var = np.where(counts >= 2,
                   np.maximum(s2 - counts * mean * mean, 0.0) / np.maximum(counts - 1, 1),
                   0.0)
    return vals, counts, mean, var


def fi_mi(N, L, lam, nu):
    """R16 missing-information estimator: I_MI = mean(N) - lam^2 E[Var(L|N)]."""
    n_mc = N.shape[0]
    _, counts, _, var_L = _binstats(N, L)
    use = counts >= 2
    e_var = np.sum((counts[use] / n_mc) * var_L[use])
    I = float(np.mean(N)) - lam * lam * float(e_var)
    tail_mass = float(np.sum(counts[counts < TAIL_THRESH]) / n_mc)
    return I, I / nu, tail_mass


def fi_cs(N, L, lam, nu):
    """R16 bias-corrected conditional-score cross-check estimator."""
    n_mc = N.shape[0]
    U = N - lam * L
    _, counts, mean_U, var_U = _binstats(N, U)
    use = counts >= 2
    contrib = (counts[use] / n_mc) * (mean_U[use] ** 2 - var_U[use] / counts[use])
    I = float(np.sum(contrib))
    return I, I / nu


def fi_hist(cm, c0, cp, nu, dlog=DLOG, alpha=HIST_ALPHA):
    """The audited estimator: log-pmf central finite difference with alpha/N floor.

    Byte-for-byte the fi_at() logic of jitter_scalar_fi_v2.py (default alpha=0.5,
    dlog=0.01), exposed with alpha/dlog knobs for the R16 sweeps.
    """
    n_mc = c0.shape[0]
    lo = int(min(cm.min(), c0.min(), cp.min()))
    hi = int(max(cm.max(), c0.max(), cp.max())) + 1
    bins = np.arange(lo, hi + 1)
    eps = alpha / n_mc
    p0 = np.histogram(c0, bins=bins)[0] / n_mc
    pp = np.maximum(np.histogram(cp, bins=bins)[0] / n_mc, eps)
    pm = np.maximum(np.histogram(cm, bins=bins)[0] / n_mc, eps)
    score = (np.log(pp) - np.log(pm)) / (2 * dlog)
    return float(np.sum(p0 * score ** 2)) / nu


# ----------------------------------------------------------------------------
# Bootstrap
# ----------------------------------------------------------------------------
def bootstrap_ci(fn, n_mc, n_boot, seed, alpha=0.05):
    """Percentile CI for a scalar estimator fn(idx) over resampled draw indices."""
    rng = np.random.default_rng(seed)
    reps = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n_mc, size=n_mc)
        reps[b] = fn(idx)
    lo = float(np.quantile(reps, alpha / 2))
    hi = float(np.quantile(reps, 1 - alpha / 2))
    return lo, hi, float(reps.std(ddof=1))


# ----------------------------------------------------------------------------
# Closed forms for the toy
# ----------------------------------------------------------------------------
def toy_I_singlecount(rho, T=1.0):
    """Exact count-only FI (theta=log rho) of the single-count truncated channel."""
    q = rho * T
    return (T * T) * (rho * rho) * np.exp(-q) / (1.0 - np.exp(-q))


def _toy_I_via_identity(rho, T=1.0):
    """Same quantity via the missing-information identity (cross-check of algebra)."""
    q = rho * T
    p1 = 1.0 - np.exp(-q)
    var_trunc = (1.0 / rho ** 2) * (1.0 - (q ** 2) * np.exp(-q) / (1.0 - np.exp(-q)) ** 2)
    return p1 - rho ** 2 * (p1 * var_trunc)             # E[N] - rho^2 E[Var(L|N)]


# ----------------------------------------------------------------------------
# Modes
# ----------------------------------------------------------------------------
def run_toy(n_mc, n_seeds=8, base_seed=20250716):
    rho, nu, cv = 1.5, 1, 0.0                 # T = nu*tau = 1, single-count regime
    closed = toy_I_singlecount(rho, T=nu * TAU)
    ident = _toy_I_via_identity(rho, T=nu * TAU)
    assert abs(closed - ident) < 1e-12, (closed, ident)
    print(f"[toy] single-count regime  rho={rho} nu={nu} tau={TAU} cv={cv}", flush=True)
    print(f"[toy] closed-form I_N (theta=log rho) = {closed:.6f}  "
          f"(direct Bernoulli == missing-info identity: {ident:.6f})", flush=True)
    print(f"[toy] N_MC={n_mc}  seeds={n_seeds}", flush=True)

    # single-seed report with bootstrap CI
    N, L, _, _ = simulate(rho, cv, nu, n_mc, base_seed, want_hist=False)
    I_mi, J_mi, tail = fi_mi(N, L, rho / TAU, nu)
    I_cs, J_cs = fi_cs(N, L, rho / TAU, nu)
    lo, hi, sd = bootstrap_ci(lambda idx: fi_mi(N[idx], L[idx], rho / TAU, nu)[0],
                              n_mc, 300, base_seed + 1)
    rel_mi = (I_mi - closed) / closed * 100
    rel_cs = (I_cs - closed) / closed * 100
    print(f"[toy] I_MI = {I_mi:.6f}  ({rel_mi:+.2f}% vs closed)  "
          f"boot95=[{lo:.6f},{hi:.6f}]  tail_mass<{TAIL_THRESH}={tail:.4f}", flush=True)
    print(f"[toy] I_CS = {I_cs:.6f}  ({rel_cs:+.2f}% vs closed)  "
          f"MI-CS agree: |{I_mi - I_cs:.2e}|", flush=True)

    # sign-unbiasedness: many independent seeds, mean relative error must straddle 0
    rels = []
    for s in range(n_seeds):
        Ns, Ls, _, _ = simulate(rho, cv, nu, n_mc, base_seed + 100 + s, want_hist=False)
        Is, _, _ = fi_mi(Ns, Ls, rho / TAU, nu)
        rels.append((Is - closed) / closed * 100)
    rels = np.array(rels)
    se = rels.std(ddof=1) / np.sqrt(len(rels))
    print(f"[toy] sign check over {n_seeds} seeds: mean rel err = "
          f"{rels.mean():+.3f}% +/- {se:.3f}% (1SE)  "
          f"[{'UNBIASED: straddles 0' if abs(rels.mean()) < 2 * se else 'check'}]",
          flush=True)
    print(f"[toy] per-seed rel errs (%): "
          f"{np.array2string(rels, precision=2, floatmode='fixed')}", flush=True)
    print(f"[toy] RESULT: closed={closed:.5f} I_MI={I_mi:.5f} rel={rel_mi:+.2f}% "
          f"N_MC={n_mc} sign_mean={rels.mean():+.3f}% se={se:.3f}%", flush=True)


def run_toy_poisson(n_mc, base_seed=20250716):
    """Weaker anchor: zero dead-hold => L=T, N~Poisson(rho*nu), I_N=rho*nu, J=rho."""
    rho, nu = 3.0, 50
    closed_J = rho                                    # J = I/nu = rho exactly
    print(f"[toyP] zero-dead-hold Poisson  rho={rho} nu={nu}  closed J=rho={rho}",
          flush=True)
    # B=0: reuse simulate with a monkeypatched zero hold via cv=0 but tau->0 is not
    # exposed; simulate directly here.
    T = nu * TAU
    rng = np.random.default_rng(base_seed)
    m_max = int(rho * nu + 12 * np.sqrt(rho * nu) + 50)
    N = np.empty(n_mc, dtype=np.int64)
    L = np.empty(n_mc, dtype=np.float64)
    pos = 0
    while pos < n_mc:
        nb = min(4000, n_mc - pos)
        A = np.cumsum(rng.exponential(1.0, size=(nb, m_max)), axis=1)
        n0 = (A <= (rho / TAU) * T).sum(axis=1)       # B=0 => D=0, slack=T
        N[pos:pos + nb] = n0
        L[pos:pos + nb] = T                           # live the whole window
        pos += nb
    I_mi, J_mi, _ = fi_mi(N, L, rho / TAU, nu)
    print(f"[toyP] I_MI={I_mi:.4f}  J_MI={J_mi:.5f}  ({(J_mi-closed_J)/closed_J*100:+.2f}% "
          f"vs closed J={closed_J})", flush=True)
    print(f"[toyP] RESULT: J_closed={closed_J} J_MI={J_mi:.5f}", flush=True)


def run_compare(rho, cv, nu, n_mc, base_seed=20250716, n_boot=200,
                dlog=DLOG, alpha=HIST_ALPHA):
    print(f"[cmp] identical-draw comparison  nu={nu} cv={cv} rho={rho} N_MC={n_mc}",
          flush=True)
    print(f"[cmp] audited histogram: dlog={dlog} floor=alpha/N with alpha={alpha}",
          flush=True)
    t0 = time.time()
    N, L, cm, cp = simulate(rho, cv, nu, n_mc, base_seed, want_hist=True, dlog=dlog)
    lam = rho / TAU
    I_mi, J_mi, tail = fi_mi(N, L, lam, nu)
    I_cs, J_cs = fi_cs(N, L, lam, nu)
    J_hist = fi_hist(cm, N, cp, nu, dlog=dlog, alpha=alpha)
    print(f"[cmp] sim+estimate wall={time.time()-t0:.1f}s  "
          f"mean(N)={np.mean(N):.2f} tail_mass<{TAIL_THRESH}={tail:.4f}", flush=True)

    lo_mi, hi_mi, _ = bootstrap_ci(lambda i: fi_mi(N[i], L[i], lam, nu)[1],
                                   n_mc, n_boot, base_seed + 1)
    lo_cs, hi_cs, _ = bootstrap_ci(lambda i: fi_cs(N[i], L[i], lam, nu)[1],
                                   n_mc, n_boot, base_seed + 2)
    lo_h, hi_h, _ = bootstrap_ci(lambda i: fi_hist(cm[i], N[i], cp[i], nu,
                                                   dlog=dlog, alpha=alpha),
                                 n_mc, n_boot, base_seed + 3)
    print(f"[cmp] J_MI   (score, no-hist,no-FD) = {J_mi:.6f}  boot95=[{lo_mi:.6f},{hi_mi:.6f}]",
          flush=True)
    print(f"[cmp] J_CS   (score cross-check)     = {J_cs:.6f}  boot95=[{lo_cs:.6f},{hi_cs:.6f}]",
          flush=True)
    print(f"[cmp] J_HIST (audited estimator)     = {J_hist:.6f}  boot95=[{lo_h:.6f},{hi_h:.6f}]",
          flush=True)
    print(f"[cmp] MI-CS gap = {J_mi - J_cs:+.2e} (must be ~0); "
          f"HIST-MI gap = {J_hist - J_mi:+.2e}", flush=True)
    print(f"[cmp] RESULT: J_MI={J_mi:.6f} J_CS={J_cs:.6f} J_HIST={J_hist:.6f} "
          f"nu={nu} cv={cv} rho={rho} N_MC={n_mc}", flush=True)
    print("[cmp] NOTE: single-point smoke at small N_MC -- machinery check only, "
          "NOT a scientific peak estimate.", flush=True)


def main():
    ap = argparse.ArgumentParser(description="R16 single-load score/MI Fisher diagnostic.")
    ap.add_argument("--toy", action="store_true",
                    help="single-count closed-form validation (primary)")
    ap.add_argument("--toy-poisson", action="store_true",
                    help="zero-dead-hold Poisson anchor (weaker; J=rho)")
    ap.add_argument("--compare", action="store_true",
                    help="score vs histogram on identical draws at one point")
    ap.add_argument("--nu", type=int, default=200)
    ap.add_argument("--cv", type=float, default=0.05)
    ap.add_argument("--rho", type=float, default=5.65)
    ap.add_argument("--nmc", type=int, default=20000)
    ap.add_argument("--nboot", type=int, default=200)
    args = ap.parse_args()

    if args.nmc > 20000:
        print(f"[warn] N_MC={args.nmc} exceeds the local 20k smoke contract; "
              "use jitter_score_diag_colab.py for the full run.", file=sys.stderr)

    ran = False
    if args.toy:
        run_toy(args.nmc); ran = True
    if args.toy_poisson:
        run_toy_poisson(args.nmc); ran = True
    if args.compare:
        run_compare(args.rho, args.cv, args.nu, args.nmc, n_boot=args.nboot); ran = True
    if not ran:
        ap.print_help()


if __name__ == "__main__":
    main()
