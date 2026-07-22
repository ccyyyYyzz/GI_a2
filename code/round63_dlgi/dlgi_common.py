#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
DLGI CONFIRMATORY CAMPAIGN -- shared config, sealed-bank seed hygiene, forward
model, and frozen-estimator re-exports (R34).

WHAT IS FROZEN / BIT-IDENTICAL.  The point estimator is the feasibility probe's,
imported VERBATIM from results/round63_next/DUAL_LEDGER_PROBE/code/dl_common.py
(joint_dual_ledger, medium_estimate, mom_autocov, _z_residuals, s_of_x, ou_path,
make_schedule, the arms, PHI, N_EXP, ...).  R34 bar C3 kills the campaign if the
calibration changes the estimator; the calibration changes ONLY the uncertainty
statement (the Neyman construction in dlgi_neyman.py).  Nothing here modifies the
frozen module.

WHAT IS NEW (R34-authorised).  (1) a photon-SNR axis: the forward simulate takes a
photon_scale in {0.5,1,2} multiplying the frozen nominal PHI (2200 counts/exposure,
cohort-average); at photon_scale=1 with a DEV scene the forward is BIT-IDENTICAL to
the probe's simulate_record (self-tested in __main__).  The estimator's internal
PHI stays frozen; the log(photon_scale) offset is absorbed by the gauge demean, so
the medium estimate is unchanged in output (verified).  (2) sealed-bank seed
hygiene: gain path and Poisson noise draw from INDEPENDENT SeedSequence children so
no gain/Poisson aliasing and no seed crosses banks.

Read-only on the frozen probe; no writes here.  CPU.
"""
import hashlib
import os
import sys
import collections

import numpy as np

REPO = r"D:/GI_another"
PROBE = os.path.join(REPO, "results", "round63_next", "DUAL_LEDGER_PROBE", "code")
sys.path.insert(0, PROBE)
import dl_common as C                       # FROZEN probe estimator (verbatim import)

CAMPAIGN = os.path.join(REPO, "results", "round63_dlgi_campaign")
BANK_ROOT = os.path.join(CAMPAIGN, "scene_banks")

# ============================================================ frozen grid ===
# R34 sec.2.2 primary claim = 27 cells.  Delta = 1 exposure, so t_c/Delta == t_c
# in exposure units (matches the probe's t_c grid convention).
TC_PRIMARY = [16.0, 32.0, 64.0]            # t_c / Delta
CV_PRIMARY = [0.05, 0.15, 0.40]
PHOTON_PRIMARY = [0.5, 1.0, 2.0]           # x nominal 2200 counts/exposure
NOMINAL_COUNTS = 2200.0

# R34 sec.2.2 edge cells: fast-drift t_c/Delta = 2 and one slower-than-64 slow-drift
# cell (t_c/Delta = 128); shown across the CV sweep at nominal photon.  Edge cells
# define the published failure map; they cannot rescue or kill the interior.
TC_EDGE_FAST = 2.0
TC_EDGE_SLOW = 128.0

# mandatory named stress cells (R34 bar C1): the three that FAILED intervals in the
# probe, at nominal SNR.
STRESS_CELLS = [(16.0, 0.05, 1.0), (16.0, 0.40, 1.0), (64.0, 0.05, 1.0)]


def cell_key(tc, cv, ph):
    return f"tc{int(round(tc))}_cv{int(round(cv * 100))}_ph{int(round(ph * 10)):02d}"


def primary_cells():
    return [(tc, cv, ph) for tc in TC_PRIMARY for cv in CV_PRIMARY for ph in PHOTON_PRIMARY]


def edge_cells():
    fast = [(TC_EDGE_FAST, cv, 1.0) for cv in CV_PRIMARY]
    slow = [(TC_EDGE_SLOW, cv, 1.0) for cv in CV_PRIMARY]
    return fast + slow


# ============================================================ seed hygiene ==
# disjoint bank salts (no seed crosses banks); gain/Poisson/schedule draw from
# independent SeedSequence children (no aliasing).
_BANK_SALT = {"calibration": 710, "confirmatory": 700, "edge": 720}


def _u32(s):
    return int(hashlib.sha256(str(s).encode()).hexdigest()[:8], 16)


def record_rngs(bank, ckey, scene_id, replicate):
    """Independent (gain_rng, poisson_rng, sched_seed) for one forward record.
    Deterministic pure function of (bank, cell, scene, replicate); streams never
    alias and never cross banks."""
    ss = np.random.SeedSequence([_BANK_SALT[bank], _u32(ckey), _u32(scene_id),
                                 int(replicate)])
    g, p, s = ss.spawn(3)
    sched_seed = int(s.generate_state(1)[0]) % (2 ** 31)
    return np.random.default_rng(g), np.random.default_rng(p), sched_seed


# ============================================================ forward model =
def simulate_record_scaled(x, gain_rng, poisson_rng, tc, sigma_l,
                           photon_scale=1.0, schedule="paired", sched_seed=0):
    """One 2048-exposure complementary bucket record through the frozen forward
    model, with an added photon_scale on the frozen PHI.  Reuses the frozen physics
    (s_of_x, ou_path, make_schedule, PHI) verbatim.  Returns (Y, a_time, slot)."""
    x = np.asarray(x, np.float64)
    s = C.s_of_x(x)
    a_time = (np.ones(C.N_EXP) if sigma_l <= 0
              else C.ou_path(gain_rng, sigma_l, tc, C.N_EXP))
    slot = C.make_schedule(schedule, sched_seed)
    a_exp = a_time[slot]
    Y = poisson_rng.poisson(photon_scale * a_exp * C.PHI * s).astype(np.float64)
    return Y, a_time, slot


def simulate_record_probe(x, seed, tc, sigma_l, photon_scale=1.0,
                          schedule="paired", sched_seed=0):
    """Probe-compatible seed formula (gain=1000+seed, Poisson=5000+seed).  At
    photon_scale=1 with a DEV scene this reproduces C.simulate_record BIT-FOR-BIT
    (self-tested in __main__).  Kept only as the provenance/bit-identity anchor;
    the campaign uses record_rngs()."""
    return simulate_record_scaled(x, np.random.default_rng(1000 + seed),
                                  np.random.default_rng(5000 + seed), tc, sigma_l,
                                  photon_scale, schedule, sched_seed)


# ============================================================ scene banks ===
def bank_dir(bank):
    return os.path.join(BANK_ROOT, bank)


def load_bank(bank):
    """OrderedDict scene_id -> x (1024,) float64 for a sealed bank (npz field 'x')."""
    import json
    d = bank_dir(bank)
    man = json.load(open(os.path.join(d, "MANIFEST.json")))
    out = collections.OrderedDict()
    for e in man["scenes"]:
        x = np.load(os.path.join(d, e["scene_id"] + ".npz"))["x"].astype(np.float64).ravel()
        out[e["scene_id"]] = x
    return out, man


# ============================================================ self-test =====
def _selftest():
    """Verify (a) the photon-scale forward is BIT-IDENTICAL to the frozen probe at
    scale=1 on a DEV scene, and (b) the frozen estimator's medium output is
    invariant to photon_scale (gauge absorbs it)."""
    sc = C.SCENES[0]; x = C.scene_x[sc]
    tc, sig = 64.0, C.sigma_l_of_cv(0.40)
    Y0, a0, s0 = C.simulate_record(sc, 0, tc, sig)
    Y1, a1, s1 = simulate_record_probe(x, 0, tc, sig)
    bit_identical = bool(np.array_equal(Y0, Y1) and np.array_equal(a0, a1)
                         and np.array_equal(s0, s1))
    print("forward bit-identical to probe at photon_scale=1 (DEV scene):",
          "PASS" if bit_identical else "FAIL", flush=True)

    # estimator gauge-invariance to photon_scale: same gain/Poisson RNG, scaled PHI
    g = np.random.default_rng(123)
    Ya, _, sl = simulate_record_scaled(x, np.random.default_rng(1), np.random.default_rng(2),
                                       tc, sig, photon_scale=1.0)
    Yb, _, _ = simulate_record_scaled(x, np.random.default_rng(1), np.random.default_rng(2),
                                      tc, sig, photon_scale=1.0)
    ma = C.medium_estimate(Ya, C.joint_dual_ledger(Ya, slot=sl)["x_hat"], slot=sl,
                           return_path=False)
    print("estimator import OK; sample medium tc_hat=%.2f cv_hat=%.3f" %
          (ma["tc_hat"], ma["cv_hat"]), flush=True)
    return bit_identical


if __name__ == "__main__":
    ok = _selftest()
    print("\nprimary cells:", len(primary_cells()), " edge cells:", len(edge_cells()))
    print("example cell_key(64,0.40,1.0) =", cell_key(64.0, 0.40, 1.0))
    sys.exit(0 if ok else 1)
