#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
SEALED SCENE + BANK GENERATOR for the beyond-modulator-band sealed probe (R39, R34 pattern).

FOUR DISJOINT, sealed partitions generated BLIND (no reconstruction / arm / bucket metric
touches these during generation; only image-domain statistics are computed):

  calibration  : construct the aperture-law prediction map + the coeff_sd gauge check ONLY.
  confirmatory : fresh scenes/seeds used ONLY for the final beyond-band endpoints.
  mismatch     : scenes used ONLY for the F5 declared-law / geometry mismatch degradation bars.
  oracle       : disjoint scenes used ONLY to anchor the oracle-ceiling reference (so the
                 ceiling is never computed on a confirmatory scene).

NO scene, and NO seed, crosses partitions.  Salts are disjoint from every prior repo cohort
and from the DLGI banks (700/710/720).  New disjoint salt bases (R39 probe):

  calibration  : 810000 + 100*stratum + instance
  confirmatory : 800000 + 100*stratum + instance
  mismatch     : 820000 + 100*stratum + instance
  oracle       : 830000 + 100*stratum + instance

STRATA (beyond-band-relevant):
  witness    synthetic in-band random + explicit super-resolution-zone spikes (the beyond-band
             witness -- the load-bearing stratum: content the modulator band CANNOT reach).
  natural    a fixed, deterministic skimage image (cameraman/coins/moon/text/...) downsampled
             to 64x64 and normalized -- honest natural-image statistics (sha256-recorded).
  onef       1/f^beta power-law random field (synthetic natural-spectrum proxy).
  smooth     low-pass Gaussian random field (in-band-heavy control; little beyond-band energy).
  texture    band-pass field concentrated in the super-resolution annulus (stress stratum).

Every image is SIDE x SIDE float64 in [0,1].  Determinism: every synthetic image is a pure
function of its integer seed; natural images are pure functions of their fixed source id.
Byte-reproducible, self-checked in __main__.  Writes ONLY under FOG_DMD_PROBE64/scene_banks/.
CPU.  [R39] tags mark thresholds/choices the external referee will freeze.
"""
import hashlib
import json
import os
import sys
from datetime import datetime, timezone

import numpy as np
from scipy.fftpack import idct
from scipy.ndimage import gaussian_filter

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_ROOT = os.path.join(HERE, "scene_banks")

SIDE = 64
N = SIDE * SIDE
PB = 8                                   # pattern band (must match the probe spec)
MED_LO, MED_HI = 1, 16                   # medium band annulus (must match the probe spec)

STRATA = ["witness", "natural", "onef", "smooth", "texture"]
BANKS = {
    "calibration": dict(base=810000, n_inst=2),
    "confirmatory": dict(base=800000, n_inst=3),
    "mismatch": dict(base=820000, n_inst=2),
    "oracle": dict(base=830000, n_inst=2),
}
# fixed, deterministic natural sources; each bank's natural slots draw GLOBALLY-UNIQUE sources
# (no natural image crosses banks) via NAT_OFFSET below. Recorded by sha256.
NATURAL_SOURCES = ["camera", "coins", "moon", "text", "clock", "gravel",
                   "brick", "page", "horse", "checkerboard", "cell", "shepp_logan_phantom"]


# cumulative natural-slot offset per bank (canonical BANKS order) -> globally-unique source per
# (bank, instance) natural slot, so no natural image ever crosses partitions.
NAT_OFFSET = {}
_acc = 0
for _bk, _sp in BANKS.items():
    NAT_OFFSET[_bk] = _acc
    _acc += _sp["n_inst"]
assert _acc <= len(NATURAL_SOURCES), "not enough distinct natural sources for all bank slots"


def I2(A):
    return idct(idct(A, axis=0, norm="ortho"), axis=1, norm="ortho")


def _superres_mask():
    """DCT frequencies in the Minkowski-sum coverage set but OUTSIDE the pattern box."""
    band = [(i, j) for i in range(SIDE) for j in range(SIDE) if MED_LO <= i + j <= MED_HI]
    pat_box = np.zeros((SIDE, SIDE), bool); pat_box[:PB + 1, :PB + 1] = True
    cov = np.zeros((SIDE, SIDE), bool)
    for pi in range(PB + 1):
        for pj in range(PB + 1):
            for (ki, kj) in band:
                for si in (pi + ki, abs(pi - ki)):
                    for sj in (pj + kj, abs(pj - kj)):
                        if si < SIDE and sj < SIDE:
                            cov[si, sj] = True
    return cov & (~pat_box), pat_box, cov


SUPERRES, PAT_BOX, COV = _superres_mask()


def img_stats(x):
    x = np.asarray(x, np.float64)
    from scipy.fftpack import dct
    D2 = dct(dct(x, axis=0, norm="ortho"), axis=1, norm="ortho") ** 2
    tot = D2.sum() + 1e-12
    return dict(mean=float(x.mean()), std=float(x.std()), min=float(x.min()), max=float(x.max()),
                superres_energy_frac=float(D2[SUPERRES].sum() / tot),
                inband_energy_frac=float(D2[PAT_BOX].sum() / tot),
                beyond_coverage_frac=float(D2[~COV].sum() / tot))


# --------------------------------------------------------------- generators ---
def gen_witness(seed):
    rs = np.random.default_rng(int(seed))
    C = np.zeros((SIDE, SIDE))
    C[:PB + 1, :PB + 1] = rs.standard_normal((PB + 1, PB + 1))
    zi, zj = np.where(SUPERRES)
    pick = rs.choice(len(zi), size=int(rs.integers(12, 24)), replace=False)
    for k in pick:
        C[zi[k], zj[k]] = 2.5 * rs.choice([-1, 1])
    x = I2(C).ravel(); x = x - x.min(); x = (x / x.max()).reshape(SIDE, SIDE)
    return x, {"generator": "bank_gen64.gen_witness (in-band random + super-res-zone spikes)",
               "family": "witness", "seed": int(seed), "n_spikes": int(len(pick))}


def gen_natural(global_idx):
    from skimage import data, transform
    name = NATURAL_SOURCES[global_idx]
    src = getattr(data, name)().astype(np.float64)
    if src.ndim == 3:
        src = src.mean(axis=2)
    img = transform.resize(src, (SIDE, SIDE), anti_aliasing=True, mode="reflect")
    img = (img - img.min()) / (np.ptp(img) + 1e-12)
    return img, {"generator": "bank_gen64.gen_natural (skimage.data, downsampled 64x64, "
                              "deterministic, sha256-recorded)", "family": "natural",
                 "source": name, "global_idx": int(global_idx)}


def gen_onef(seed):
    rng = np.random.default_rng(int(seed))
    w = rng.standard_normal((SIDE, SIDE)); F = np.fft.fft2(w)
    fy = np.fft.fftfreq(SIDE)[:, None]; fx = np.fft.fftfreq(SIDE)[None, :]
    fr = np.sqrt(fy ** 2 + fx ** 2); fr[0, 0] = 1.0
    beta = float(rng.uniform(0.9, 1.4))
    img = np.real(np.fft.ifft2(F / fr ** beta))
    img = (img - img.min()) / (np.ptp(img) + 1e-12)
    return np.clip(0.05 + 0.90 * img, 0, 1), {
        "generator": "bank_gen64.gen_onef (1/f^beta power-law field)", "family": "onef",
        "seed": int(seed), "beta": beta}


def gen_smooth(seed):
    rng = np.random.default_rng(int(seed))
    f = rng.standard_normal((SIDE, SIDE))
    sig = float(rng.uniform(4.0, 8.0))
    sm = gaussian_filter(f, sig, mode="reflect")
    sm = (sm - sm.min()) / (np.ptp(sm) + 1e-12)
    return np.clip(0.06 + 0.90 * sm, 0, 1), {
        "generator": "bank_gen64.gen_smooth (low-pass Gaussian field)", "family": "smooth",
        "seed": int(seed), "blur_sigma": sig}


def gen_texture(seed):
    """Band-pass field concentrated in the super-resolution annulus (a beyond-band stress)."""
    rng = np.random.default_rng(int(seed))
    C = np.zeros((SIDE, SIDE))
    C[:PB + 1, :PB + 1] = 0.5 * rng.standard_normal((PB + 1, PB + 1))
    zi, zj = np.where(SUPERRES)
    C[zi, zj] = rng.standard_normal(len(zi)) * 1.2
    x = I2(C).ravel(); x = x - x.min(); x = (x / x.max()).reshape(SIDE, SIDE)
    return x, {"generator": "bank_gen64.gen_texture (band-pass, super-res annulus)",
               "family": "texture", "seed": int(seed)}


def generate(stratum, seed, inst, bank):
    if stratum == "witness":
        return gen_witness(seed)
    if stratum == "natural":
        return gen_natural(NAT_OFFSET[bank] + inst)
    if stratum == "onef":
        return gen_onef(seed)
    if stratum == "smooth":
        return gen_smooth(seed)
    if stratum == "texture":
        return gen_texture(seed)
    raise ValueError(stratum)


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def bank_table(bank):
    spec = BANKS[bank]; base = spec["base"]; n_inst = spec["n_inst"]
    return [(f"{bank}_{st}_{i}", st, i, base + 100 * si + i)
            for si, st in enumerate(STRATA) for i in range(n_inst)]


def build_bank(bank, write=True):
    out_dir = os.path.join(OUT_ROOT, bank)
    if write:
        os.makedirs(out_dir, exist_ok=True)
    entries = []
    for scene_id, stratum, inst, seed in bank_table(bank):
        x, params = generate(stratum, seed, inst, bank)
        x = np.asarray(x, np.float64)
        assert x.shape == (SIDE, SIDE) and np.all(np.isfinite(x))
        assert x.min() >= 0.0 and x.max() <= 1.0 and x.var() > 0.0, scene_id
        x_sha = sha256_bytes(np.ascontiguousarray(x, np.float64).tobytes())
        params = dict(params); params["x_sha256"] = x_sha
        entry = dict(scene_id=scene_id, bank=bank, stratum=stratum, instance=inst,
                     seed=int(seed), gen_params=params, x_sha256=x_sha, stats=img_stats(x))
        if write:
            np.savez(os.path.join(out_dir, scene_id + ".npz"), x=x)
        entries.append(entry)
    manifest = dict(
        campaign="R39 beyond-modulator-band sealed probe (R34 pattern)",
        bank=bank, role={
            "calibration": "aperture-law prediction map + coeff_sd gauge check ONLY",
            "confirmatory": "final beyond-band endpoints ONLY (fresh scenes/seeds)",
            "mismatch": "F5 declared-law / geometry mismatch degradation bars ONLY",
            "oracle": "oracle-ceiling reference anchor ONLY (disjoint from confirmatory)"}[bank],
        side=SIDE, pattern_band=PB, medium_band=[MED_LO, MED_HI], value_range=[0.0, 1.0],
        seed_base=BANKS[bank]["base"], instances_per_stratum=BANKS[bank]["n_inst"],
        strata=STRATA, n_scenes=len(entries),
        seed_disjoint_from=[700000, 710000, 720000, "the other three R39 partitions"],
        blind_generation="no reconstruction, arm, or bucket metric computed during generation",
        superres_frac=float(SUPERRES.mean()), coverage_frac=float(COV.mean()),
        utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
        python=sys.version.split()[0], numpy=np.__version__, scenes=entries)
    if write:
        json.dump(manifest, open(os.path.join(out_dir, "MANIFEST.json"), "w"), indent=2)
    return manifest


def main():
    print("=== R39 beyond-band sealed scene banks ===", flush=True)
    all_sha = {}
    for bank in BANKS:
        man = build_bank(bank, write=True)
        print(f"\n[{bank}] {man['n_scenes']} scenes base={man['seed_base']} "
              f"(superres-frac {man['superres_frac']:.3f})", flush=True)
        for e in man["scenes"]:
            all_sha[e["scene_id"]] = e["x_sha256"]
            s = e["stats"]
            print(f"  {e['scene_id']:28s} {e['stratum']:8s} seed={e['seed']:6d}  "
                  f"mean {s['mean']:.3f}  superres-E {s['superres_energy_frac']:.3f}  "
                  f"{e['x_sha256'][:12]}", flush=True)
    # determinism + disjointness self-check
    print("\n=== determinism (regenerate, no write) + all-distinct ===", flush=True)
    ok = True; seen = {}
    for bank in BANKS:
        for e in build_bank(bank, write=False)["scenes"]:
            if e["x_sha256"] != all_sha[e["scene_id"]]:
                ok = False; print("  MISMATCH", e["scene_id"])
            if e["x_sha256"] in seen:
                ok = False; print("  DUPLICATE", e["scene_id"], "==", seen[e["x_sha256"]])
            seen[e["x_sha256"]] = e["scene_id"]
    print("determinism + all-distinct across banks:", "PASS" if ok else "FAIL", flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
