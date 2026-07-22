#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
DLGI CONFIRMATORY CAMPAIGN -- the four arms and three schedules (R34 sec.2.3),
ported from the feasibility probe with the photon-SNR axis and sealed-bank seeds.

ARMS (R34 sec.2.3):
  1. pilot-free DLGI     : C.joint_dual_ledger on the ordinary bucket record (frozen).
  2. pilot-interleaved   : the frozen strongest 5% pilot design; pilots REPLACE scene
                           rows, identical total exposures/photons/duration; scene cost
                           charged (ported from C.pilot_interleaved, photon-scaled).
  3. oracle monitor      : C.oracle_monitor_estimate on the gain path observed at every
                           exposure (nondeployable ceiling; frozen).
  4. plain-linear camera : C.arm_A2 on the BYTE-IDENTICAL pilot-free record, no medium
                           correction (scene noninferiority endpoint; frozen).

SCHEDULES (R34 sec.2.3): {paired, random, ordered}.  Only temporal order changes; the
pattern multiset, photons, exposure count and duration are identical.  The pilot-free
schedule is what varies; the pilot baseline is always paired (the mandatory strongest
baseline).  NAMING RECONCILIATION (flagged): the probe implements schedules named
{paired, split, random}; R34 names {paired, random, ordered}.  "ordered" == the probe's
deterministic non-paired "split" order (all m+ then all m-, complementary members
maximally separated in time -> the large-cross-confusion ordered acquisition).  The
frozen make_schedule is reused verbatim through this fixed rename.

The raw record is preserved; no differencing step discards the common-mode medium
information (R34 sec.2.3).  Read-only on the frozen probe.  CPU.
"""
import numpy as np

import dlgi_common as K
import dl_common as C                                   # frozen estimator (via K's path)

SCHEDULES = ["paired", "random", "ordered"]
_SCHED_MAP = {"paired": "paired", "random": "random", "ordered": "split"}


def _arm_rngs(bank, ckey, scene_id, rep):
    """Shared gain path across arms (fair paired comparison on ONE gain realization);
    independent Poisson streams per arm; no seed crosses banks."""
    ss = np.random.SeedSequence([K._BANK_SALT[bank], K._u32(ckey), K._u32(scene_id),
                                 int(rep)])
    g, pf, pp, s = ss.spawn(4)
    return (np.random.default_rng(g), np.random.default_rng(pf),
            np.random.default_rng(pp), int(s.generate_state(1)[0]) % (2 ** 31))


def pilot_interleaved_scaled(x, a_time, poisson_rng, sigma_l, photon_scale=1.0,
                             pilot_frac=0.05):
    """Frozen pilot-interleaved baseline (C.pilot_interleaved), ported to a passed
    scene x, a SHARED gain path a_time, an explicit Poisson rng, and a photon scale.
    Pilots replace 5% of Hadamard rows by all-ones pilot pairs (paired schedule)."""
    x = np.asarray(x, np.float64); sumx = float(x.sum())
    n_pair = C.NPIX
    n_pilot = max(int(round(pilot_frac * n_pair)), 2)
    pilot_rows = np.unique(np.linspace(1, n_pair - 1, n_pilot).round().astype(int))
    is_pilot = np.zeros(n_pair, bool); is_pilot[pilot_rows] = True
    s_meas = C.s_of_x(x)
    s_meas[2 * pilot_rows] = sumx
    s_meas[2 * pilot_rows + 1] = sumx
    Y = poisson_rng.poisson(photon_scale * a_time * C.PHI * s_meas).astype(np.float64)
    pe = np.concatenate([2 * pilot_rows, 2 * pilot_rows + 1]); pe.sort()
    a_reads = Y[pe] / (photon_scale * C.PHI * sumx)
    l_reads = np.log(np.maximum(a_reads, 1e-9))
    tc_p, sig_p = C.ou_fit_irregular(pe, l_reads)
    a_interp = np.interp(np.arange(C.N_EXP), pe, a_reads)
    Yp = Y[0::2]; Ym = Y[1::2]
    c = (Yp / np.maximum(a_interp[0::2], 1e-9)
         - Ym / np.maximum(a_interp[1::2], 1e-9))
    c[pilot_rows] = 0.0
    x_hat = np.clip(C.rof_denoise(C.hadamard_inverse(c) / (photon_scale * C.PHI),
                                  C.LAM_TV), 0.0, None)
    med = dict(tc_hat=float(tc_p), cv_hat=C.cv_of_sigma_l(sig_p),
               sigma_l_hat=float(sig_p))
    return dict(med=med, x_hat=x_hat, n_pilot=int(len(pilot_rows)),
                a_interp=a_interp, is_pilot_row=is_pilot, Y=Y)


def run_arms(bank, ckey, scene_id, x, rep, tc, sigma_l, cv, ph, schedule="paired"):
    """All four arms on one (scene, cell, schedule, replicate).  The pilot-free
    schedule is `schedule`; the pilot baseline is always paired.  Returns medium
    estimates, scene PSNRs (vs the true x), gain-path correlation, and the identical-
    ledger accounting (bar C7)."""
    x = np.asarray(x, np.float64)
    gain_rng, pf_rng, pp_rng, sched_seed = _arm_rngs(bank, ckey, scene_id, rep)

    # shared gain realization
    a_time = (np.ones(C.N_EXP) if sigma_l <= 0
              else C.ou_path(gain_rng, sigma_l, tc, C.N_EXP))
    slot = C.make_schedule(_SCHED_MAP[schedule], sched_seed)
    a_exp = a_time[slot]
    s = C.s_of_x(x)

    # arm 1 + arm 4 share the byte-identical pilot-free record
    Y = pf_rng.poisson(ph * a_exp * C.PHI * s).astype(np.float64)
    r = C.joint_dual_ledger(Y, slot=slot, n_outer=2)
    m = r["med"]
    pf = dict(tc_hat=float(m["tc_hat"]), cv_hat=float(m["cv_hat"]),
              gaincorr=C.gain_path_corr(m["a_time"], a_time),
              scene_psnr=C.psnr(r["x_hat"], x))
    lin = dict(scene_psnr=C.psnr(C.arm_A2(Y[0::2], Y[1::2], C.LAM_A2), x))

    # arm 2 pilot (paired baseline, shared gain path)
    p = pilot_interleaved_scaled(x, a_time, pp_rng, sigma_l, photon_scale=ph)
    pilot = dict(tc_hat=float(p["med"]["tc_hat"]), cv_hat=float(p["med"]["cv_hat"]),
                 scene_psnr=C.psnr(p["x_hat"], x), n_pilot=p["n_pilot"])

    # arm 3 oracle monitor (nondeployable ceiling)
    om = C.oracle_monitor_estimate(a_time)
    oracle = dict(tc_hat=float(om["tc_hat"]), cv_hat=float(om["cv_hat"]))

    accounting = dict(
        pilot_free_exposures=int(Y.size), pilot_free_photons=float(Y.sum()),
        plain_linear_exposures=int(Y.size), plain_linear_photons=float(Y.sum()),
        identical_record_pilot_free_vs_linear=True,
        pilot_exposures=int(p["Y"].size), pilot_photons=float(p["Y"].sum()),
        pilot_extra_photons_or_exposures=0,
        note="pilot-free and plain-linear consume the byte-identical record; the "
             "medium ledger is post-processing at zero extra acquisition")
    return dict(cell=ckey, scene_id=scene_id, schedule=schedule, rep=int(rep),
                tc_true=tc, cv_true=cv, ph=ph,
                pilot_free=pf, plain_linear=lin, pilot=pilot, oracle=oracle,
                accounting=accounting)


def _smoke():
    scenes, _ = K.load_bank("confirmatory")
    sid = list(scenes)[0]; x = scenes[sid]
    tc, cv, ph = 64.0, 0.40, 1.0; sig = C.sigma_l_of_cv(cv)
    ckey = K.cell_key(tc, cv, ph)
    print(f"[arms smoke] scene={sid} cell={ckey}", flush=True)
    for sch in SCHEDULES:
        r = run_arms("confirmatory", ckey, sid, x, 0, tc, sig, cv, ph, schedule=sch)
        pf, pl, pi, orc = r["pilot_free"], r["plain_linear"], r["pilot"], r["oracle"]
        print(f"  {sch:7s}: PF tc={pf['tc_hat']:6.1f} cv={pf['cv_hat']:.3f} "
              f"PSNR={pf['scene_psnr']:.2f} corr={pf['gaincorr']:.3f} | "
              f"linPSNR={pl['scene_psnr']:.2f} | pilot tc={pi['tc_hat']:6.1f} "
              f"PSNR={pi['scene_psnr']:.2f} | oracle tc={orc['tc_hat']:6.1f}",
              flush=True)
    a = r["accounting"]
    print(f"  accounting: PF/linear identical record = "
          f"{a['identical_record_pilot_free_vs_linear']} "
          f"({a['pilot_free_exposures']} exp, {a['pilot_free_photons']:.0f} photons)",
          flush=True)


if __name__ == "__main__":
    _smoke()
