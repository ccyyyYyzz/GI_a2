# E1 -- TRACKING BLIND RECOVERY grid.  Blind-from-cold null-err vs fixed-A baseline 1.000.
# Grid: d_w in {16,64}, sig_w in {0.15,0.30}, tau in {4,16,64}, noise in {clean,1e6,1e5}.
# Methods per config: ORACLE (upper bound), TRACK-EM cold (blind, the headline), JOINT-ADAM cold.
# SUCCESS SIGNAL: any blind-from-cold null-err <= 0.3.
import sys, time, json, argparse
import numpy as np
import torch
from fog_common import (make_scene, smooth_basis, make_patterns, projectors,
                        null_metric, ou_coeffs)
import fog_tracker as ft

def log(*a): print(*a, flush=True)

def build_case(n_side, M, c, sig_w, tau, T, photons, patt_kind, seed):
    """Generate one experiment instance. Returns everything needed for solves + metric."""
    rng = np.random.default_rng(seed)
    N = n_side * n_side
    x_true = make_scene(n_side)
    P = make_patterns(M, N, patt_kind, rng)
    Pi, rangeP = projectors(P)
    null_true = x_true - rangeP(x_true); nrm0 = np.linalg.norm(null_true)
    U = smooth_basis(c, n_side); d_w = U.shape[1]
    res = ou_coeffs(T, d_w, sig_w, tau, rng)
    Z, rho = res if isinstance(res, tuple) else (res, 0.0)
    W = np.clip(1.0 + Z @ U.T, 0.05, None)
    mu = np.einsum('mn,tn->tm', P, W * x_true)            # clean buckets (>=0 for binary P)
    if photons is None:
        B = mu.copy(); wts = np.ones_like(B); R_diag = np.ones_like(B) * 1e-8
    else:
        s = photons / mu.mean()                          # scale so mean bucket = photons
        Y = rng.poisson(np.clip(mu, 0, None) * s) / s
        B = Y
        wts = (s / np.clip(mu, 1e-9, None))              # Poisson weights 1/var
        R_diag = np.clip(mu, 1e-9, None) / s             # obs variance in bucket units
    return dict(P=P, U=U, W=W, B=B, wts=wts, R_diag=R_diag, rho=rho, x_true=x_true,
                rangeP=rangeP, null_true=null_true, nrm0=nrm0, d_w=d_w)

def metric(x, cs):
    return null_metric(x, cs['x_true'], cs['rangeP'], cs['null_true'], cs['nrm0'])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='E1_results.json')
    ap.add_argument('--T', type=int, default=256)
    ap.add_argument('--n_em', type=int, default=25)
    ap.add_argument('--adam', type=int, default=1)         # 1=run joint-adam too
    ap.add_argument('--adam_steps', type=int, default=5000)
    ap.add_argument('--quick', type=int, default=0)        # 1 = tiny grid for timing
    args = ap.parse_args()
    dev = 'cuda' if torch.cuda.is_available() else 'cpu'
    log("device:", dev, torch.cuda.get_device_name(0) if dev == 'cuda' else '')
    n_side = 32; M = 64
    lam_x = 1e-6
    d_w_list = [16] if args.quick else [16, 64]
    sig_list = [0.30] if args.quick else [0.15, 0.30]
    tau_list = [16] if args.quick else [4, 16, 64]
    noise_list = [None] if args.quick else [None, 1e6, 1e5]
    c_of = {16: 4, 64: 8}
    results = []
    t_start = time.time()
    for d_w in d_w_list:
        for sig_w in sig_list:
            for tau in tau_list:
                for photons in noise_list:
                    seed = 20260723 + d_w * 131 + int(sig_w * 100) * 17 + tau * 7 + \
                           (0 if photons is None else int(np.log10(photons)))
                    cs = build_case(n_side, M, c_of[d_w], sig_w, tau, args.T,
                                    photons, 'gauss' if photons is None else 'binary', seed)
                    # regularization: heavier x-ridge under noise
                    lam = lam_x if photons is None else 1e-3
                    tag = f"d_w={d_w} sig={sig_w} tau={tau} ph={'clean' if photons is None else f'{photons:.0e}'}"
                    row = dict(d_w=d_w, sig_w=sig_w, tau=tau,
                               photons=None if photons is None else float(photons), T=args.T)
                    # oracle
                    t0 = time.time()
                    x = ft.solve_oracle(cs['P'], cs['W'], cs['B'], cs['wts'], lam, dev)
                    ne, te = metric(x, cs); row['oracle_ne'] = ne; row['oracle_te'] = te
                    # track-EM cold (headline blind method)
                    x = ft.track_em(cs['P'], cs['U'], cs['B'], cs['wts'], cs['rho'], sig_w,
                                    cs['R_diag'], lam, n_em=args.n_em, anneal=True, dev=dev)
                    ne, te = metric(x, cs); row['em_ne'] = ne; row['em_te'] = te
                    # joint-adam cold
                    if args.adam:
                        x = ft.joint_adam(cs['P'], cs['U'], cs['B'], cs['wts'], cs['rho'],
                                          sig_w, lam, steps=args.adam_steps, lr=2e-2, dev=dev)
                        ne, te = metric(x, cs); row['adam_ne'] = ne; row['adam_te'] = te
                    row['sec'] = time.time() - t0
                    best_blind = min([row.get('em_ne', 9), row.get('adam_ne', 9)])
                    row['best_blind_ne'] = best_blind
                    results.append(row)
                    log(f"{tag:42s} | oracle {row['oracle_ne']:.3f} | EM {row['em_ne']:.3f}"
                        f" | adam {row.get('adam_ne', float('nan')):.3f} | best {best_blind:.3f}"
                        f" | {row['sec']:.1f}s")
                    with open(args.out, 'w') as f:
                        json.dump(dict(meta=dict(n_side=n_side, M=M, T=args.T,
                                   n_em=args.n_em, adam_steps=args.adam_steps,
                                   wall_s=time.time() - t_start), results=results), f, indent=2)
    best = min(r['best_blind_ne'] for r in results)
    log(f"\n=== E1 best blind-from-cold null-err = {best:.3f} "
        f"({'SUCCESS <=0.3' if best <= 0.3 else 'no success signal'}) ===")
    log(f"total wall {time.time() - t_start:.1f}s -> {args.out}")

if __name__ == '__main__':
    main()
