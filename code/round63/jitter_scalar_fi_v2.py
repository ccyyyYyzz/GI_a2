"""Candidate-B decisive test: SCALAR count Fisher information of the
jittered nonparalyzable renewal channel, by Monte Carlo.

Channel: active start; wait W_i ~ Exp(lam); hold B_i (deterministic tau,
or lognormal mean tau with given cv). Count N(lam) = #{m >= 1:
sum_{i<m}(W_i+B_i) + W_m <= T}. With W_i = A_i/lam (A = cumsum of unit
exponentials) and D_m = sum_{i<m} B_i, the condition is
A_m <= lam*(T - D_m), so counts for MANY lam values come from ONE (A, D)
draw (common random numbers: exact pathwise coupling in lam).

FI estimate at lam0: symmetric log-pmf finite difference at
lam0*exp(+-dlog) using the SAME draws; J = I/nu (per-slot, log-rate).
Grid: nu=2000; rho geomspace 2..40; cv in {0, 0.05, 0.1, 0.3}.
Prediction (R14 candidate B): peak rho stays ~22.9 at cv=0 and is CAPPED
~ c*cv^(-2/3) for cv > 0.
"""
import time

import numpy as np

NU = int(__import__('sys').argv[1]) if len(__import__('sys').argv) > 1 else 2000
TAU = 1.0
T = NU * TAU
RHOS = np.unique(np.round(np.geomspace(1.2, 40.0, 25), 3))
CVS = [0.0, 0.02, 0.05, 0.1, 0.2, 0.3]
N_MC = 100_000
CHUNK = 6_000
DLOG = 0.01

def counts_three(rho, cv, seed):
    """Counts at lam0*exp({-dlog,0,+dlog}) from common random numbers."""
    lam0 = rho / TAU
    lams = lam0 * np.exp(np.array([-DLOG, 0.0, DLOG]))
    m_max = int(min(NU + 1, NU * rho / (1 + rho) * (1 + 3 * DLOG)
                    + 8 * np.sqrt(NU) + 40))
    rng = np.random.default_rng(seed)
    out = [np.empty(N_MC, dtype=np.int32) for _ in range(3)]
    if cv > 0:
        sig2 = np.log1p(cv * cv)
        mu = -0.5 * sig2
    pos = 0
    while pos < N_MC:
        nb = min(CHUNK, N_MC - pos)
        A = np.cumsum(rng.exponential(1.0, size=(nb, m_max)), axis=1)
        if cv > 0:
            B = TAU * np.exp(mu + np.sqrt(sig2)
                             * rng.standard_normal((nb, m_max)))
        else:
            B = np.full((nb, m_max), TAU)
        D = np.concatenate([np.zeros((nb, 1)), np.cumsum(B, axis=1)[:, :-1]],
                           axis=1)
        slack = T - D                       # lam-independent
        for j, lam in enumerate(lams):
            reg = A <= lam * slack          # monotone rows -> count via sum
            out[j][pos:pos + nb] = reg.sum(axis=1, dtype=np.int32)
        pos += nb
    return out

def fi_at(rho, cv, seed):
    cm, c0, cp = counts_three(rho, cv, seed)
    lo = int(min(cm.min(), c0.min(), cp.min()))
    hi = int(max(cm.max(), c0.max(), cp.max())) + 1
    bins = np.arange(lo, hi + 1)
    eps = 0.5 / N_MC
    p0 = np.histogram(c0, bins=bins)[0] / N_MC
    pp = np.maximum(np.histogram(cp, bins=bins)[0] / N_MC, eps)
    pm = np.maximum(np.histogram(cm, bins=bins)[0] / N_MC, eps)
    score = (np.log(pp) - np.log(pm)) / (2 * DLOG)
    return float(np.sum(p0 * score ** 2)) / NU

print(f"[sfi] nu={NU} n_mc={N_MC} dlog={DLOG} rhos={list(RHOS)}", flush=True)
print("[sfi]   rho     " + "  ".join(f"cv={cv:<5g}" for cv in CVS), flush=True)
tab = {}
t0 = time.time()
for rho in RHOS:
    row = []
    for k, cv in enumerate(CVS):
        row.append(fi_at(rho, cv, seed=int(rho * 1000) + 977 * k))
    tab[rho] = row
    print(f"[sfi] {rho:7.3f}  " + "  ".join(f"{J:8.4f}" for J in row)
          + f"   ({time.time()-t0:.0f}s)", flush=True)

print("\n[sfi] === peaks ===", flush=True)
for k, cv in enumerate(CVS):
    peak_rho = max(tab, key=lambda r: tab[r][k])
    cap = "" if cv == 0 else f"; cap~c*{cv ** (-2.0/3.0):.1f}"
    print(f"[sfi] cv={cv:<5g}: peak J={max(t[k] for t in tab.values()):.4f} "
          f"at rho={peak_rho}  (cv=0 prediction 22.9{cap})", flush=True)
print("[sfi] DONE", flush=True)
