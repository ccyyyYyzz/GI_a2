#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
FIVE SEALED BANKS for the R41 sealed detection probe (sec 4.3), k_p=5 detection geometry.

  calibration   : 6 scenes / independent seeds -> build score filters, nuisance bases, thresholds,
                  cross-fit baseline, and the fixed-vs-fresh forecast. (The ONLY bank opened now.)
  confirmatory  : 12 FRESH scenes (6 natural + 6 synthetic) + fresh code/medium/shot/anomaly seeds.
                  Supplies every primary endpoint. SEALED until the coordinator's freeze order.
  specificity   : base scenes + the SIX event classes (H0, in-band scene, beyond-band scene,
                  medium-amplitude, medium-timescale tau x{0.5,2}, mixed scene+medium). SEALED.
  mismatch      : base scenes + the D5 mismatch axes (basis rotation 10/20%, spectral slope +-1,
                  band-edge +-20%, static envelope, shot +-20%, convolutive blend 10/25/50%). SEALED.
  oracle        : disjoint scenes for the exact mean-wall / true-law ceiling only. SEALED.

NO scene, code seed, medium seed, anomaly seed, or event seed crosses banks (disjoint salt bases,
distinct from the FOG_DMD_PROBE64 banks 800000/810000/820000/830000). Every scene is byte-committed
by sha256; every stochastic ingredient is a committed integer seed. Thresholds (sealed_common.BARS)
are frozen BEFORE this generation. Writes only under SEALED_DET/banks/.
"""
import hashlib
import json
import os
import sys
from datetime import datetime, timezone

import numpy as np

import sealed_common as sc

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "banks")

# disjoint salt bases (SEALED_DET series; distinct from every prior repo cohort)
SALT = dict(calibration=840000, confirmatory=841000, specificity=842000,
            mismatch=843000, oracle=844000)


def sha(x):
    return hashlib.sha256(np.ascontiguousarray(x, np.float64).tobytes()).hexdigest()


def synth_scene(kind, seed):
    """Synthetic scenes in the k_p=5 Fourier convention (nonneg, [0,1])."""
    rng = np.random.default_rng(int(seed))
    if kind == "witness":
        x = sc.witness_scene(seed)
    elif kind == "onef":
        w = rng.standard_normal((sc.n, sc.n)); F = np.fft.fft2(w)
        fy = np.fft.fftfreq(sc.n)[:, None]; fx = np.fft.fftfreq(sc.n)[None, :]
        fr = np.sqrt(fy ** 2 + fx ** 2); fr[0, 0] = 1.0
        beta = float(rng.uniform(0.9, 1.4))
        img = np.real(np.fft.ifft2(F / fr ** beta)).ravel()
        x = (img - img.min()) / (np.ptp(img) + 1e-12)
    elif kind == "texture":
        # band-pass concentrated in the beyond-band annulus (stress stratum)
        Ube = sc.band_modes(sc.KP + 1, 9)
        a_in = sc.U_in_np @ (0.5 * rng.standard_normal(sc.U_in_np.shape[1]))
        a_be = Ube @ (1.2 * rng.standard_normal(Ube.shape[1]))
        x = a_in + a_be; x = x - x.min(); x = x / x.max()
    elif kind == "smooth":
        from scipy.ndimage import gaussian_filter
        f = rng.standard_normal((sc.n, sc.n))
        sm = gaussian_filter(f, float(rng.uniform(4, 8)), mode="reflect").ravel()
        x = (sm - sm.min()) / (np.ptp(sm) + 1e-12)
    else:
        raise ValueError(kind)
    return np.asarray(x, np.float64)


def nat_scene(name):
    return sc.natural_scene(name)


# --------- bank scene tables (globally-unique scenes; no scene crosses banks) ----------
# natural sources partitioned across banks so none repeats
NAT = dict(
    calibration=["camera", "coins"],
    confirmatory=["moon", "text", "clock", "page", "gravel", "brick"],
    specificity=["horse"],
    mismatch=["cell"],
    oracle=["checkerboard"],
)
# synthetic slots (kind, local_index) per bank
SYN = dict(
    calibration=[("witness", 0), ("onef", 0), ("texture", 0), ("smooth", 0)],
    confirmatory=[("witness", 1), ("witness", 2), ("witness", 3), ("onef", 1), ("texture", 1), ("smooth", 1)],
    specificity=[("witness", 4), ("texture", 2)],
    mismatch=[("witness", 5), ("onef", 2)],
    oracle=[("witness", 6)],
)


_KIND_OFFSET = {"witness": 1, "onef": 2, "texture": 3, "smooth": 4}   # deterministic (no PYTHONHASHSEED)


def scene_seed(bank, kind, idx):
    return SALT[bank] + 1000 * _KIND_OFFSET[kind] + idx


def build_scenes(bank):
    entries = []
    for name in NAT[bank]:
        x = nat_scene(name)
        entries.append(dict(scene_id=f"{bank}_nat_{name}", family="natural", source=name,
                            seed=None, x_sha256=sha(x)))
    for (kind, idx) in SYN[bank]:
        seed = scene_seed(bank, kind, idx)
        x = synth_scene(kind, seed)
        entries.append(dict(scene_id=f"{bank}_{kind}_{idx}", family=kind, source=None,
                            seed=int(seed), x_sha256=sha(x)))
    return entries


# --------- per-bank stochastic ingredient seeds (committed; used only when the bank is opened) -----
def code_seed(bank):
    return SALT[bank] + 10          # the sealed code bank seed for this partition


def medium_seed(bank, rep):
    return SALT[bank] + 20000 + int(rep)


def anomaly_seed(bank, rep):
    return SALT[bank] + 30000 + int(rep)


SPECIFICITY_EVENTS = [
    dict(cls="H0", desc="declared scene, no change"),
    dict(cls="inband_scene", desc="in-band scene change, matched norm (lights mean channel)"),
    dict(cls="beyond_scene", desc="beyond-band scene change (headline target)"),
    dict(cls="medium_amplitude", desc="sigma_f x {0.8,1.2} (amplitude/law direction)"),
    dict(cls="medium_timescale", desc="tau x {0.5,2} (lag/correlation-time direction)"),
    dict(cls="mixed_scene_medium", desc="simultaneous beyond-band scene + medium change"),
]

MISMATCH_AXES = [
    dict(axis="basis_rotation", levels=[0.10, 0.20], desc="declared medium subspace rotation"),
    dict(axis="spectral_slope", levels=[-1, +1], desc="declared power-law slope error"),
    dict(axis="band_edge", levels=[-0.20, +0.20], desc="declared medium band-edge error"),
    dict(axis="tau_scale", levels=[0.5, 2.0], desc="declared correlation-time error"),
    dict(axis="shot_level", levels=[0.8, 1.2], desc="declared photon/shot level error"),
    dict(axis="static_envelope", levels=[0.10], desc="10% RMS static multiplicative envelope"),
    dict(axis="convolutive_blend", levels=[0.10, 0.25, 0.50],
         desc="multiplicative->convolutive medium blend (0.50 = mapped edge, not a pass condition)"),
]


def build_bank(bank, write=True):
    scenes = build_scenes(bank)
    manifest = dict(
        campaign="R41 sealed detection probe (SEALED_DET)",
        bank=bank,
        role={
            "calibration": "build score filters/nuisance bases/thresholds/cross-fit baseline; ONLY bank opened in the dry run",
            "confirmatory": "all primary endpoints (12 fresh scenes: 6 natural + 6 synthetic); SEALED",
            "specificity": "six event classes incl tau x{0.5,2} and mixed scene+medium; SEALED",
            "mismatch": "D5 mismatch axes; SEALED",
            "oracle": "exact mean-wall / true-law ceiling anchor; SEALED",
        }[bank],
        frozen_ledger=dict(n=sc.n, N=sc.N, k_p=sc.KP, M=sc.M, exposures_per_bank=sc.N_COMPLEMENTARY,
                           ms_per_bank=sc.MS_PER_BANK, T_cap=sc.T_CAP, T_cap_sec=round(sc.T_CAP_SEC, 2),
                           photons_per_bucket=sc.PHOT, tau=sc.TAU, rho=round(sc.RHO, 4)),
        seed_salt_base=SALT[bank],
        seed_disjoint_from=[v for k, v in SALT.items() if k != bank] +
        [800000, 810000, 820000, 830000, "FOG_DMD_PROBE64 banks", "DLGI 700/710/720"],
        code_seed=code_seed(bank), medium_seed_base=SALT[bank] + 20000,
        anomaly_seed_base=SALT[bank] + 30000,
        n_scenes=len(scenes), scenes=scenes,
        thresholds_frozen_before_generation=True,
        bars_reference="sealed_common.BARS (frozen R41 sec 4.6)",
        utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
        python=sys.version.split()[0], numpy=np.__version__,
    )
    if bank == "specificity":
        manifest["event_classes"] = SPECIFICITY_EVENTS
        manifest["eps_levels"] = sc.EPS_LEVELS
        manifest["min_records_per_class"] = 250
    if bank == "mismatch":
        manifest["mismatch_axes"] = MISMATCH_AXES
        manifest["min_records_per_class"] = 250
    if write:
        os.makedirs(os.path.join(OUT, bank), exist_ok=True)
        json.dump(manifest, open(os.path.join(OUT, bank, "MANIFEST.json"), "w"), indent=2)
    return manifest


def main():
    print("=== SEALED_DET five sealed banks ===", flush=True)
    all_sha = {}
    bank_hashes = {}
    for bank in SALT:
        man = build_bank(bank, write=True)
        h = hashlib.sha256(json.dumps(man["scenes"], sort_keys=True).encode()).hexdigest()
        bank_hashes[bank] = h
        print(f"\n[{bank}] {man['n_scenes']} scenes  salt={man['seed_salt_base']}  "
              f"scene-set sha={h[:16]}", flush=True)
        for e in man["scenes"]:
            if e["x_sha256"] in all_sha:
                print(f"  !! DUPLICATE scene {e['scene_id']} == {all_sha[e['x_sha256']]}", flush=True)
            all_sha[e["x_sha256"]] = e["scene_id"]
            print(f"  {e['scene_id']:30s} {e['family']:9s} {e['x_sha256'][:16]}", flush=True)
    # determinism + cross-bank disjointness
    ok = True
    for bank in SALT:
        for e in build_bank(bank, write=False)["scenes"]:
            if e["x_sha256"] not in all_sha or all_sha[e["x_sha256"]] != e["scene_id"]:
                ok = False
    n_unique = len(all_sha)
    n_total = sum(len(build_bank(b, write=False)["scenes"]) for b in SALT)
    print(f"\nunique scenes = {n_unique} of {n_total} total  ->  "
          f"{'ALL DISJOINT' if n_unique == n_total else 'COLLISION'}", flush=True)
    print(f"determinism/disjointness: {'PASS' if ok and n_unique == n_total else 'FAIL'}", flush=True)
    # write a top-level manifest index with the committed bank hashes
    idx = dict(campaign="R41 sealed detection probe", banks=list(SALT.keys()),
               salt_bases=SALT, bank_scene_set_sha256=bank_hashes,
               n_unique_scenes=n_unique, n_total_scenes=n_total,
               utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"))
    os.makedirs(OUT, exist_ok=True)
    json.dump(idx, open(os.path.join(OUT, "BANKS_INDEX.json"), "w"), indent=2)
    return 0 if (ok and n_unique == n_total) else 1


if __name__ == "__main__":
    sys.exit(main())
