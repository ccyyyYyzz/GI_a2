#!/usr/bin/env python
# run_e6.py -- E6 driver. Five design arms x two estimators; three endpoints per cell.
#   endpoint 1 TRUE-SEED STABILITY (pathwise): seed medium from truth -> STAY(<=0.1) or DRIFT(~0.7)?
#   endpoint 2 BLIND COLD-START null-NMSE (5-seed median) + multi-start agreement.
#   endpoint 3 (A2/A3) certificate<->outcome correlation across seeds.
# Fixed cell: N=1024, M=96, d=48, sigma_f=0.15, OU tau=16, T=128 blocks, S=96 slots/block,
#            total exposures = T*S = 12288 (same across arms), clean + 1e5 photons/bucket.
import argparse, json, time, sys, os
import numpy as np
import fog_e6 as fe

CFG = dict(n_side=32, N=1024, M=96, d=48, sigma_f=0.15, tau=16.0, T=128, S=96,
           photons_list=[None, 1e5], lam_x=1e-4)


def build_arm(arm, seed, cfg, P_bank, U, x_true):
    """Return dict(P, idx, idx_flat, W, Z, rho, coeff_sd, ou_lambda) for this arm+seed."""
    T, S, M, N = cfg['T'], cfg['S'], cfg['M'], cfg['N']
    rng_med = np.random.default_rng(seed * 100 + 1)
    rng_sch = np.random.default_rng(seed * 100 + 3)
    # medium: A4 = multi-timescale; others single-tau (matches E4/E5 tau=16)
    if arm == 'A4':
        tau_arr = fe.tau_schedule(cfg['d'], 64.0, 4.0)
        W, Z, rho, csd = fe.lognormal_medium_mt(U, T, cfg['sigma_f'], tau_arr, rng_med)
    else:
        W, Z, rho, csd = fe.lognormal_medium(U, T, cfg['sigma_f'], cfg['tau'], rng_med)
    # schedule: slot (randomized chronology) for A2/A3/A4; cartesian for A0/A1
    if arm in ('A2', 'A3', 'A4'):
        idx, idx_flat = fe.slot_idx(T, S, M, rng_sch)
    else:
        idx = fe.cartesian_idx(T, M)
        idx_flat = idx.reshape(-1).copy()
    ou_lambda = 3.0 if arm == 'A4' else 0.0
    return dict(P=P_bank[arm], idx=idx, idx_flat=idx_flat, W=W, Z=Z, rho=rho,
                coeff_sd=csd, ou_lambda=ou_lambda)


def run_cell(arm, photons, seeds, cfg, P_bank, U, x_true, rangeP, null_true, nrm0,
             dev, kn):
    """One (arm, photon) cell over all seeds. Returns per-seed records + summary."""
    recs = []
    for sd in seeds:
        t0 = time.time()
        A = build_arm(arm, sd, cfg, P_bank, U, x_true)
        P, idx, W, Z, rho, csd = A['P'], A['idx'], A['W'], A['Z'], A['rho'], A['coeff_sd']
        rng_noise = np.random.default_rng(sd * 100 + 7)
        b_clean = fe.forward_clean(P, idx, W, x_true)
        b, Rdiag = fe.add_photons(b_clean, photons, rng_noise)

        # oracle (known medium)
        x_or = fe.oracle_solve(P, idx, W, b, cfg['lam_x'], dev)
        or_nmse = fe.null_metric(x_or, x_true, rangeP, null_true, nrm0)[0]

        # (i) pathwise TRUE-SEED STABILITY (full-Z refinement seeded at truth -> STAY or DRIFT)
        x_ts = fe.pathwise_fullz_refine(P, U, idx, b, csd, rho, Z, cfg['lam_x'], dev=dev,
                                        outer=kn['ts_outer'], zsteps=kn['ts_zsteps'],
                                        zlr=6e-3, ou_lambda=0.05)
        pw_ts = fe.null_metric(x_ts, x_true, rangeP, null_true, nrm0)[0]
        # (i) pathwise BLIND cold-start
        r_bl = fe.pathwise_solve(P, U, idx, b, csd, rho, cfg['lam_x'], dev=dev,
                                 z_seed=None, n_starts=kn['starts'], seed0=sd,
                                 outer=kn['outer'], qsteps=kn['qsteps'], ou_lambda=A['ou_lambda'])
        pw_bl_starts = [fe.null_metric(x, x_true, rangeP, null_true, nrm0)[0]
                        for x in r_bl['per_start_x']]
        pw_bl = fe.null_metric(r_bl['x_best'], x_true, rangeP, null_true, nrm0)[0]

        # (ii) covariance-domain BLIND (marginalizes medium -> no true-seed test)
        r_cv = fe.cov_solve(P, U, idx, b, csd, rho, dev=dev, n_starts=kn['starts'],
                            steps=kn['cov_steps'], seed0=sd)
        cv_starts = [fe.null_metric(x, x_true, rangeP, null_true, nrm0)[0]
                     for x in r_cv['per_start_x']]
        cv_bl = fe.null_metric(r_cv['x_best'], x_true, rangeP, null_true, nrm0)[0]

        cert = fe.carrier_certificate(P, A['idx_flat'], x_true)
        recs.append(dict(seed=int(sd), oracle_nmse=or_nmse,
                         pw_trueseed_nmse=pw_ts, pw_blind_nmse=pw_bl,
                         pw_blind_starts=pw_bl_starts,
                         cov_blind_nmse=cv_bl, cov_blind_starts=cv_starts,
                         cert=cert, secs=round(time.time() - t0, 1)))
        print(f"    [{arm} ph={photons} sd={sd}] or={or_nmse:.3f} "
              f"pw_TS={pw_ts:.3f} pw_BL={pw_bl:.3f} cv_BL={cv_bl:.3f} "
              f"sf={cert['sf_fourier']:.3f} ({recs[-1]['secs']}s)", flush=True)

    def med(key):
        return float(np.median([r[key] for r in recs]))

    def agree(key):  # median within-seed spread across starts
        return float(np.median([np.std(r[key]) for r in recs]))
    # certificate<->outcome correlation across seeds
    sf = np.array([r['cert']['sf_fourier'] for r in recs])
    bl = np.array([r['pw_blind_nmse'] for r in recs])
    corr = float(np.corrcoef(sf, bl)[0, 1]) if len(recs) > 2 and sf.std() > 1e-9 else None
    summ = dict(oracle_median=med('oracle_nmse'),
                pw_trueseed_median=med('pw_trueseed_nmse'),
                pw_blind_median=med('pw_blind_nmse'),
                pw_blind_agreement=agree('pw_blind_starts'),
                cov_blind_median=med('cov_blind_nmse'),
                cov_blind_agreement=agree('cov_blind_starts'),
                cert_sf_median=float(np.median(sf)),
                cert_vs_blind_corr=corr)
    return dict(seeds=recs, summary=summ)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--smoke', action='store_true')
    ap.add_argument('--arms', default='A0,A1,A2,A3,A4')
    ap.add_argument('--seeds', default='0,1,2,3,4')
    ap.add_argument('--dev', default='cuda')
    ap.add_argument('--out', default='E6_results.json')
    args = ap.parse_args()
    cfg = dict(CFG)
    if args.smoke:
        cfg['T'] = 48; cfg['S'] = 96
        kn = dict(outer=8, qsteps=6, starts=2, cov_steps=400, ts_outer=6, ts_zsteps=20)
        seeds = [0, 1]
    else:
        kn = dict(outer=40, qsteps=22, starts=4, cov_steps=4000, ts_outer=30, ts_zsteps=60)
        seeds = [int(s) for s in args.seeds.split(',')]
    arms = args.arms.split(',')

    import torch
    dev = args.dev if (args.dev != 'cuda' or torch.cuda.is_available()) else 'cpu'
    print(f"device={dev} torch={torch.__version__} np={np.__version__} smoke={args.smoke}", flush=True)

    n_side = cfg['n_side']
    x_true = fe.make_scene(n_side)
    U = fe.dct_basis(cfg['d'], n_side)
    # pattern banks (fixed per arm, seed-independent)
    rng_pat = np.random.default_rng(777)
    P_bin = fe.make_binary_bank(cfg['M'], cfg['N'], rng_pat)
    P_four, four_meta = fe.make_fourier_overlap_bank(cfg['M'], n_side, spacing=3, band=4.0)
    P_bank = {'A0': P_bin, 'A1': P_four, 'A2': P_bin, 'A3': P_four, 'A4': P_bin}
    arm_desc = {
        'A0': 'CONTROL: random-binary bank, Cartesian schedule (E4 config)',
        'A1': 'Fourier lattice + overlapping sidebands, Cartesian schedule',
        'A2': 'random-binary bank, CHRONOLOGY-RANDOMIZED slot schedule (GI_a1 port)',
        'A3': 'Fourier lattice + chronology-randomized slot schedule (A1+A2)',
        'A4': 'A2 base (random+slot) + multi-timescale medium + temporal-law penalty',
    }

    results = dict(config=dict(cfg, exposures=cfg['T'] * cfg['S'], four_meta=four_meta,
                               kn=kn, seeds=seeds, arm_desc=arm_desc),
                   null_fraction=fe.null_fraction(x_true, fe.projectors(P_bin)[1]),
                   arms={})
    outpath = args.out

    for arm in arms:
        P = P_bank[arm]
        Pi, rangeP = fe.projectors(P)
        null_true = x_true - rangeP(x_true)
        nrm0 = np.linalg.norm(null_true)
        results['arms'][arm] = dict(desc=arm_desc[arm], null_fraction=float(nrm0 ** 2 / np.linalg.norm(x_true) ** 2))
        for photons in cfg['photons_list']:
            pk = 'clean' if photons is None else f'{photons:g}'
            print(f"== ARM {arm} | photons={pk} ==", flush=True)
            cell = run_cell(arm, photons, seeds, cfg, P_bank, U, x_true, rangeP,
                            null_true, nrm0, dev, kn)
            results['arms'][arm][pk] = cell
            with open(outpath, 'w') as f:
                json.dump(results, f, indent=1)
    print("DONE ->", outpath, flush=True)


if __name__ == '__main__':
    main()
