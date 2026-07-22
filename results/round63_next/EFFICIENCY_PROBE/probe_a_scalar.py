"""PROBE A -- SCALAR efficiency of the RQL rate estimator vs the jitter-capped
information ceiling (CRB), at the M1 operating points.

Reuses the FROZEN, validated score-identity machinery (fi_mi / fi_cs /
simulate_multi) from code/round63/jitter_score_diag_colab_frozen.py verbatim
(imported), and cross-checks fi_mi against physics.exact_fisher_analytic at c=0.

RQL scalar rate estimate (closed form; minimizer of the RQL per-frame objective
Q(lam;N)=(T-N*tau)*lam - N*log lam, identical to physics.precorrect_rates):
    lam_hat = N / max(T - N*tau, floor).

Efficiency (log-rate theta=log lam, the parameterization of J):
    I_frame = nu * J(rho;c)   [per-frame Fisher info about log-rate]
    CRB_theta = 1 / I_frame
    eff = CRB_theta / Var(theta_hat).
Also reports rate-level efficiency and MSE-based efficiency (incl. bias).
"""
import json, os, sys, time
import numpy as np

ROOT = r"D:\GI_another"
sys.path.insert(0, os.path.join(ROOT, "code", "round63"))
import jitter_score_diag_colab_frozen as jd   # frozen validated estimators
import physics

TAU = 1.0
NUS = [2000, 200]
RHOS = [1.0, 5.7, 22.25]
CS = [0.0, 0.05]
N_MC = 150_000
OUT = os.path.join(ROOT, "results", "round63_next", "EFFICIENCY_PROBE")

def rql_rate(N, T, tau):
    """RQL scalar rate estimate lam_hat = N/(T-N*tau), clipped at the
    physics.precorrect_rates saturation floor (denom >= 1/(1+1e4))."""
    r = np.maximum(N, 0.0) / T
    denom = np.maximum(1.0 - r * tau, 1.0 / (1.0 + 1e4))
    return r / denom

def run_point(nu, rho, c, seed):
    T = nu * TAU
    lam = rho / TAU
    N, L, _ = jd.simulate_multi(rho, c, nu, N_MC, seed, mults=[1.0])
    # ceiling: per-slot log-rate info J, per-frame info = nu*J
    J_mi, tail = jd.fi_mi(N, L, lam, nu)
    J_cs = jd.fi_cs(N, L, lam, nu)
    I_frame = nu * J_mi
    crb_theta = 1.0 / I_frame
    # RQL estimator
    lam_hat = rql_rate(N.astype(np.float64), T, TAU)
    n_sat = int(np.sum(N.astype(np.float64) * TAU >= T))
    theta_hat = np.log(lam_hat)
    var_theta = float(np.var(theta_hat, ddof=1))
    bias_theta = float(np.mean(theta_hat) - np.log(lam))
    mse_theta = var_theta + bias_theta ** 2
    eff_var = crb_theta / var_theta
    eff_mse = crb_theta / mse_theta
    # rate-level cross-check
    I_lam_frame = I_frame / (lam ** 2)     # I_lam = I_theta / lam^2
    crb_lam = 1.0 / I_lam_frame
    var_lam = float(np.var(lam_hat, ddof=1))
    bias_lam = float(np.mean(lam_hat) - lam)
    eff_lam = crb_lam / var_lam
    eff_lam_mse = crb_lam / (var_lam + bias_lam ** 2)
    # c=0 exact analytic cross-check
    J_exact = None
    if c == 0.0:
        J_exact = physics.exact_fisher_analytic(lam, T, TAU) / nu
    return dict(nu=nu, rho=rho, c=c, n_mc=N_MC, mean_N=float(N.mean()),
                fano_N=float(N.var() / N.mean()), n_sat=n_sat,
                J_mi=J_mi, J_cs=J_cs, J_exact_c0=J_exact,
                mi_cs_reldiff=abs(J_mi - J_cs) / J_mi,
                I_frame=I_frame, crb_theta=crb_theta, var_theta=var_theta,
                bias_theta=bias_theta, mse_theta=mse_theta,
                eff_theta_var=eff_var, eff_theta_mse=eff_mse,
                crb_lam=crb_lam, var_lam=var_lam, bias_lam=bias_lam,
                eff_lam_var=eff_lam, eff_lam_mse=eff_lam_mse,
                headroom_dB=10.0 * np.log10(1.0 / eff_var))

def main():
    t0 = time.time()
    results = []
    print("[scalar] nu   rho     c     meanN   J_mi      J_cs      Jexact   "
          "eff_var  eff_mse  head_dB", flush=True)
    for nu in NUS:
        for c in CS:
            for rho in RHOS:
                seed = int(rho * 1000) + 7919 * int(c * 1000) + 101 * nu
                r = run_point(nu, rho, c, seed)
                results.append(r)
                je = "  exact" if r["J_exact_c0"] is None else f"{r['J_exact_c0']:.5f}"
                print(f"[scalar] {nu:<4d} {rho:6.2f} {c:4.2f} {r['mean_N']:8.1f} "
                      f"{r['J_mi']:.5f}  {r['J_cs']:.5f}  {je}  "
                      f"{r['eff_theta_var']:.4f}  {r['eff_theta_mse']:.4f}  "
                      f"{r['headroom_dB']:+.3f}   ({time.time()-t0:.0f}s)", flush=True)
    with open(os.path.join(OUT, "scalar_efficiency.json"), "w") as f:
        json.dump({"n_mc": N_MC, "tau": TAU, "results": results}, f, indent=2)
    # c=0 validation summary
    print("\n[scalar] === c=0 fi_mi vs exact analytic (validation) ===", flush=True)
    for r in results:
        if r["c"] == 0.0:
            rel = abs(r["J_mi"] - r["J_exact_c0"]) / r["J_exact_c0"]
            print(f"[scalar] nu={r['nu']} rho={r['rho']}: J_mi={r['J_mi']:.6f} "
                  f"J_exact={r['J_exact_c0']:.6f} reldiff={rel*100:.3f}%", flush=True)
    print("[scalar] DONE", flush=True)

if __name__ == "__main__":
    main()
