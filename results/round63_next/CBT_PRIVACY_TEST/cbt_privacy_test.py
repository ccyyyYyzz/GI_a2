#!/usr/bin/env python
# cbt_privacy_test.py -- CONDITIONAL-MUTUAL-INFORMATION certification of the privacy corollary
#   of the CURVED_BLINDNESS_TEST theorem (JOB 3).
#
# Direct information-theoretic statement: the fully-scrambled sensor record Y^T carries NO
# information about WHERE on the invariant orbit (the great circle) the scene parameter sits,
# once the orbit RADIUS R is given.  Because the great circle keeps the sufficient statistic
# Q EXACTLY constant, the record law is identical for every orbit phase Omega -> I(Omega;Y|R)=0.
# The record DOES carry the radius (Q0 = Q_bg + R^2 depends on R) -> I(R;Y) large (positive
# control).  An OFF-orbit trajectory at the SAME speed lets Q vary with phase -> I(Omega;Y|R)>0
# (specificity control).
#
# Reuses the FROZEN world + generator VERBATIM: build_world, build_geometry (direction basis
# + background), gen_banks_bstream, sample_cov, efficient_weight.  No edits to the exam script.
#
# FROZEN BARS (preregistered before compute):
#   P1  I(Omega;Y|R) < 0.01 bit  AND  inside the label-shuffled null band          (on orbit)
#   P2  I(R;Y) > 0.5 bit                                                            (positive control)
#   P3  I(Omega;Y|R) > 0.1 bit                                                      (off orbit, same speed)
#   CRN pairing: on-orbit and off-orbit arms share per-sample speckle seeds.
#   MI estimator: binned plug-in, equal-frequency Y bins, phase-binned Omega, R-stratified
#   conditional MI, debiased + banded by label-shuffled nulls (resolves ~0.05 bit).
# WRITE-SCOPE: results/round63_next/CBT_PRIVACY_TEST/ ONLY.
import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
CBT = os.path.join(ROOT, "results", "round63_next", "CURVED_BLINDNESS_TEST")
SCRAMBLE = os.path.join(ROOT, "results", "round63_next", "SCRAMBLE_EXT")
JET = os.path.join(ROOT, "results", "round63_next", "JET_TEST")
for p in (CBT, SCRAMBLE, JET):
    sys.path.insert(0, p)

import curved_blindness_test as cbt  # noqa: E402
from curved_blindness_test import build_world, build_geometry  # noqa: E402
from scramble_toy import gen_banks_bstream, sample_cov, efficient_weight  # noqa: E402

# ----------------------------------------------------------------- frozen config / bars
# Radius levels kept in the INFORMATIVE small-background regime (Q_bg ~ 12, I_Q ~ 0.04, like
# the frozen exam), so the orbit-radius signal Q_xi = R^2 clears the per-record score noise
# SD ~ 1/sqrt(B*I_Q).  Scaling radii up inflates B0 -> Q_bg -> shrinks I_Q, cancelling the
# gap gain (minGap/SD pinned ~ sqrt(B)); hence B=16384 (the frozen exam's sensitive detector)
# rather than larger radii.  Gaussian-model prediction at this config: I(R;Y)=0.95 (>0.5),
# off-orbit cMI=0.30 (>0.1), on-orbit cMI~0 (within null band).
R_LEVELS_AMP = [0.35, 0.70, 1.05]     # amp multipliers of BASE_AMP -> radii r0 ~ [0.15,0.30,0.45]
BASE_AMP = 0.08
OMEGA_SPEED = 0.7                     # angular speed s/r0 (fixed across radii and arms)
GEOM_SEED = 460
KO = 6                                # phase bins for Omega
KY = 6                                # equal-frequency bins for Y
N_SHUFFLE = 300
DEFAULT_B = 16384                     # banks per record (frozen; sizes the score-noise floor)
DEFAULT_N = 6000                      # records per arm (2000 per radius stratum)

FROZEN = {
    "description": "Conditional-MI certification of orbit-phase privacy; bars preregistered.",
    "design": {
        "R_levels_amp_mult": R_LEVELS_AMP, "base_amp": BASE_AMP,
        "omega_speed": OMEGA_SPEED, "geom_seed": GEOM_SEED,
        "orbit": "great circle scene(R,Omega)=B0+R*(xi0u cosOmega + vu sinOmega); Q exactly constant",
        "off_orbit": "same-speed ellipse using the radial velocity vru (not G-tangent); Q varies",
        "record": "Y = <W, sample_cov(banks)>  (efficient matched score, fixed W at mid radius)",
        "estimator": "binned plug-in MI (bits), R-stratified conditional, shuffle-null debiased",
        "KO_phase_bins": KO, "KY_Y_bins": KY, "n_shuffle": N_SHUFFLE,
        "B_banks": DEFAULT_B, "N_records_per_arm": DEFAULT_N,
    },
    "bars": {
        "P1_onorbit_cMI_lt_bit": 0.01,
        "P1_within_shuffle_null_band": True,
        "P2_I_R_Y_gt_bit": 0.5,
        "P3_offorbit_cMI_gt_bit": 0.1,
        "resolves_bit": 0.05,
    },
    "note": "A failed bar is a reported result, not a retry license (preregistration).",
}


# ----------------------------------------------------------------- binned MI estimators (bits)
def _mi_from_joint(counts):
    """Plug-in MI in bits from a joint count table."""
    N = counts.sum()
    if N <= 0:
        return 0.0
    p = counts / N
    px = p.sum(1, keepdims=True)
    py = p.sum(0, keepdims=True)
    nz = p > 0
    denom = (px @ py)
    return float(np.sum(p[nz] * np.log2(p[nz] / denom[nz])))


def _equal_freq_bins(y, k):
    """Map continuous y to k equal-frequency bin labels (ties handled by rank)."""
    r = np.argsort(np.argsort(y))            # ranks 0..N-1
    return np.minimum((r * k) // len(y), k - 1).astype(int)


def _mi_disc_cont(labels, y, ky):
    """MI(label; y) in bits, y binned equal-frequency into ky bins."""
    yb = _equal_freq_bins(y, ky)
    kl = int(labels.max()) + 1
    C = np.zeros((kl, ky))
    np.add.at(C, (labels, yb), 1.0)
    return _mi_from_joint(C)


def _phase_bins(omega, ko):
    return np.minimum((omega % (2 * np.pi)) / (2 * np.pi) * ko, ko - 1).astype(int)


def mi_omega_y_given_R(Rlab, omega, y, ko, ky, rng, n_shuffle):
    """R-stratified conditional MI  I(Omega;Y|R) in bits, with shuffle-null band.
    Returns dict(raw, debiased, null_mean, null_lo, null_hi, within_band)."""
    ob = _phase_bins(omega, ko)
    levels = np.unique(Rlab)
    N = len(y)

    def cmi_from(ob_arr):
        tot = 0.0
        for lv in levels:
            m = Rlab == lv
            if m.sum() < ko + ky:
                continue
            yb = _equal_freq_bins(y[m], ky)
            C = np.zeros((ko, ky))
            np.add.at(C, (ob_arr[m], yb), 1.0)
            tot += (m.sum() / N) * _mi_from_joint(C)
        return tot

    raw = cmi_from(ob)
    nulls = np.empty(n_shuffle)
    for s in range(n_shuffle):
        ob_s = ob.copy()
        for lv in levels:                       # shuffle phase WITHIN each R stratum
            idx = np.where(Rlab == lv)[0]
            ob_s[idx] = rng.permutation(ob_s[idx])
        nulls[s] = cmi_from(ob_s)
    nm, lo, hi = float(nulls.mean()), float(np.percentile(nulls, 2.5)), float(np.percentile(nulls, 97.5))
    return dict(raw=float(raw), debiased=float(raw - nm), null_mean=nm,
                null_lo=lo, null_hi=hi, within_band=bool(lo <= raw <= hi))


def mi_R_y(Rlab, y, ky, rng, n_shuffle):
    """MI(R;Y) in bits with shuffle-null debiasing."""
    raw = _mi_disc_cont(Rlab, y, ky)
    nulls = np.array([_mi_disc_cont(rng.permutation(Rlab), y, ky) for _ in range(n_shuffle)])
    return dict(raw=float(raw), debiased=float(raw - nulls.mean()),
                null_mean=float(nulls.mean()), null_hi=float(np.percentile(nulls, 97.5)))


# ----------------------------------------------------------------- geometry / scenes
def build_privacy_geom(world):
    """Fixed direction basis + background sized for the LARGEST radius so every orbit stays
    physical (scene >= 0).  Returns unit G-tangent basis (xi0u, vu), radial unit vru, B0, and
    the shot/detector reference at the mid radius."""
    gmax = build_geometry(world, amp=BASE_AMP * max(R_LEVELS_AMP), s_frac=OMEGA_SPEED, seed=GEOM_SEED)
    r0max = gmax["r0"]
    xi0u = gmax["xi0"] / r0max
    vu = gmax["v"] / gmax["s"]
    vru = gmax["v_rad"] / gmax["s"]                 # radial unit (||vru||_G = 1), same speed arm
    radii = [BASE_AMP * a / BASE_AMP for a in R_LEVELS_AMP]  # placeholder; real R from build_geometry
    # real radii from build_geometry at each amp (linear in amp)
    R_of = []
    for a in R_LEVELS_AMP:
        g = build_geometry(world, amp=BASE_AMP * a, s_frac=OMEGA_SPEED, seed=GEOM_SEED)
        R_of.append(g["r0"])
    R_of = np.array(R_of)

    # background headroom: cover both orbit families over a full phase sweep at R_max
    ph = np.linspace(0, 2 * np.pi, 512)
    orb = R_of.max() * (np.outer(np.cos(ph), xi0u) + np.outer(np.sin(ph), vu))
    orb_off = R_of.max() * (np.outer(np.cos(ph), xi0u) + np.outer(np.sin(ph), vru))
    B0 = float(-min(orb.min(), orb_off.min()) + max(0.05, 0.3 * BASE_AMP * max(R_LEVELS_AMP)))

    # fixed detector W at the mid radius operating point
    Rmid = R_of[len(R_of) // 2]
    scene_mid = B0 + Rmid * xi0u
    C0, Oii = world["C0"], world["Oii"]
    Eb = C0 * Oii * scene_mid.sum()
    scale = cbt.PHOT / max(Eb.mean(), 1e-12)
    shot_mid = Eb / scale
    Qfun = world["Qfun"]
    Q0_mid = float(Qfun(scene_mid))
    W = efficient_weight(world["H"], Q0_mid, shot_mid)

    Q0_of = np.array([float(Qfun(B0 + R * xi0u)) for R in R_of])
    return dict(xi0u=xi0u, vu=vu, vru=vru, B0=B0, R_of=R_of, W=W,
                Q0_of=Q0_of, Q_bg=float(Qfun(B0 * np.ones(cbt.N))), Rmid=float(Rmid))


def scene_onorbit(pg, R, Omega):
    return pg["B0"] + R * (pg["xi0u"] * np.cos(Omega) + pg["vu"] * np.sin(Omega))


def scene_offorbit(pg, R, Omega):
    return pg["B0"] + R * (pg["xi0u"] * np.cos(Omega) + pg["vru"] * np.sin(Omega))


# ----------------------------------------------------------------- ensemble (real generator)
def gen_ensemble(world, pg, scene_fn, N, B, seed_base):
    """N records; each draws a radius level (uniform) and phase (uniform), generates B paired
    banks with a per-record seed (CRN across arms via shared seed_base), returns Y (score),
    R-level index, Omega."""
    A, O_half, mask = world["A"], world["O_half"], world["mask"]
    W = pg["W"]
    nlev = len(pg["R_of"])
    rs = np.random.default_rng(seed_base)
    lev = rs.integers(0, nlev, N)
    omega = rs.uniform(0, 2 * np.pi, N)
    Y = np.empty(N)
    for i in range(N):
        R = pg["R_of"][lev[i]]
        sc = scene_fn(pg, R, omega[i])
        # reuse the frozen paired generator; single scene -> pass sc for both, take b0.
        b0, _ = gen_banks_bstream(B, A, O_half, mask, sc, sc, phot=cbt.PHOT,
                                  seed=(seed_base + 1 + i))
        Y[i] = float(np.sum(W * sample_cov(b0)))
    return dict(Y=Y, Rlab=lev.astype(int), omega=omega)


# ----------------------------------------------------------------- selftest (synthetic, CPU)
def selftest():
    """Validate the estimator resolves ~0.05 bit and calibrates on shuffled nulls, with NO
    generator.  Three synthetic ensembles: (a) Y depends only on R (I(Om;Y|R)=0, I(R;Y)>0),
    (b) Y depends weakly on Omega within R (known small cMI), (c) strong Omega dependence."""
    rng = np.random.default_rng(7)
    N = 6000
    nlev = 3
    Rlab = rng.integers(0, nlev, N)
    Rmean = np.array([0.0, 1.0, 2.0])[Rlab]
    omega = rng.uniform(0, 2 * np.pi, N)

    # (a) pure radius dependence, no phase info
    Ya = Rmean + rng.normal(0, 0.25, N)
    # (b) small phase modulation added on top of radius (weak leak)
    Yb = Rmean + 0.06 * np.sin(omega) + rng.normal(0, 0.25, N)
    # (c) strong phase modulation (off-orbit analogue)
    Yc = Rmean + 0.6 * np.sin(omega) + rng.normal(0, 0.25, N)

    out = {}
    for tag, Y in [("a_radius_only", Ya), ("b_weak_phase", Yb), ("c_strong_phase", Yc)]:
        cmi = mi_omega_y_given_R(Rlab, omega, Y, KO, KY, np.random.default_rng(1), 200)
        iRY = mi_R_y(Rlab, Y, KY, np.random.default_rng(2), 200)
        out[tag] = dict(cMI_debiased=cmi["debiased"], cMI_within_band=cmi["within_band"],
                        cMI_null_hi=cmi["null_hi"], I_R_Y_debiased=iRY["debiased"])
        print("  [selftest] %-16s cMI_deb=%+.4f within_band=%s I(R;Y)=%.3f" %
              (tag, cmi["debiased"], cmi["within_band"], iRY["debiased"]), flush=True)
    ok = (out["a_radius_only"]["cMI_debiased"] < 0.01 and out["a_radius_only"]["cMI_within_band"]
          and out["c_strong_phase"]["cMI_debiased"] > 0.1
          and out["a_radius_only"]["I_R_Y_debiased"] > 0.5)
    print("  [selftest] estimator resolves the 0.05-bit regime: %s" % ok, flush=True)
    return out, bool(ok)


# ----------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["selftest", "run"], default="selftest")
    ap.add_argument("--N", type=int, default=DEFAULT_N, help="records per arm")
    ap.add_argument("--B", type=int, default=DEFAULT_B, help="banks per record")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    import scramble_toy as st
    print("[cbt_privacy_test] mode=%s device=%s" % (args.mode, st.DEV), flush=True)

    if args.mode == "selftest":
        res, ok = selftest()
        payload = dict(frozen=FROZEN, kind="selftest", selftest=res, selftest_ok=ok)
        out = args.out or os.path.join(HERE, "CBT_PRIVACY_TEST_selftest.json")
        with open(out, "w") as f:
            json.dump(payload, f, indent=2, default=float)
        print("[cbt_privacy_test] wrote %s" % out, flush=True)
        return

    t0 = time.time()
    world = build_world()
    pg = build_privacy_geom(world)
    print("[cbt_privacy_test] radii R=%s  Q0(R)=%s  Q_bg=%.4f  B0=%.4f" %
          (np.round(pg["R_of"], 4).tolist(), np.round(pg["Q0_of"], 4).tolist(),
           pg["Q_bg"], pg["B0"]), flush=True)

    seed_base = 90000
    on = gen_ensemble(world, pg, scene_onorbit, args.N, args.B, seed_base)
    print("[cbt_privacy_test] on-orbit ensemble done (%.1fs)" % (time.time() - t0), flush=True)
    off = gen_ensemble(world, pg, scene_offorbit, args.N, args.B, seed_base)   # CRN-shared seeds
    print("[cbt_privacy_test] off-orbit ensemble done (%.1fs)" % (time.time() - t0), flush=True)

    rng = np.random.default_rng(11)
    cmi_on = mi_omega_y_given_R(on["Rlab"], on["omega"], on["Y"], KO, KY, rng, N_SHUFFLE)
    iRY = mi_R_y(on["Rlab"], on["Y"], KY, rng, N_SHUFFLE)
    cmi_off = mi_omega_y_given_R(off["Rlab"], off["omega"], off["Y"], KO, KY, rng, N_SHUFFLE)

    bars = {
        "P1_onorbit_cMI_lt": {"value": cmi_on["debiased"],
                              "thresh": FROZEN["bars"]["P1_onorbit_cMI_lt_bit"],
                              "within_null_band": cmi_on["within_band"],
                              "null_band": [cmi_on["null_lo"], cmi_on["null_hi"]],
                              "pass": bool(cmi_on["debiased"] < FROZEN["bars"]["P1_onorbit_cMI_lt_bit"]
                                           and cmi_on["within_band"])},
        "P2_I_R_Y_gt": {"value": iRY["debiased"], "thresh": FROZEN["bars"]["P2_I_R_Y_gt_bit"],
                        "pass": bool(iRY["debiased"] > FROZEN["bars"]["P2_I_R_Y_gt_bit"])},
        "P3_offorbit_cMI_gt": {"value": cmi_off["debiased"], "thresh": FROZEN["bars"]["P3_offorbit_cMI_gt_bit"],
                               "within_null_band": cmi_off["within_band"],
                               "pass": bool(cmi_off["debiased"] > FROZEN["bars"]["P3_offorbit_cMI_gt_bit"])},
    }
    verdict = "PRIVACY_CERTIFIED" if all(b["pass"] for b in bars.values()) else "BAR_FAILED"
    payload = dict(
        frozen=FROZEN, kind="privacy_test", verdict=verdict,
        config=dict(N_records=args.N, B_banks=args.B, N_SIDE=cbt.N_SIDE, N=cbt.N, M=cbt.M,
                    PHOT=cbt.PHOT, R_of=pg["R_of"].tolist(), Q0_of=pg["Q0_of"].tolist(),
                    Q_bg=pg["Q_bg"], B0=pg["B0"]),
        mi=dict(onorbit_cMI=cmi_on, I_R_Y=iRY, offorbit_cMI=cmi_off),
        bars=bars, runtime_sec=round(time.time() - t0, 1))
    out = args.out or os.path.join(HERE, "CBT_PRIVACY_TEST.json")
    with open(out, "w") as f:
        json.dump(payload, f, indent=2, default=float)
    print("[cbt_privacy_test] VERDICT %s | P1 cMI=%.4f(band=%s) P2 I(R;Y)=%.3f P3 offcMI=%.4f | %.1fs -> %s" %
          (verdict, cmi_on["debiased"], cmi_on["within_band"], iRY["debiased"],
           cmi_off["debiased"], payload["runtime_sec"], out), flush=True)


if __name__ == "__main__":
    main()
