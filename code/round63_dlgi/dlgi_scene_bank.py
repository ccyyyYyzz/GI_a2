#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
DLGI CONFIRMATORY CAMPAIGN -- sealed scene banks (R34 sec.2.1).

Three DISJOINT, sealed scene sets, generated BLIND (no reconstruction, no arm, no
bucket metric touches these images; only image-domain statistics are computed):

  calibration  : used ONLY to construct the critical-value table c_0.95(eta) (R34 sec.1).
  confirmatory : fresh scenes/seeds used ONLY for the final endpoints.
  edge         : used ONLY to display the fast/slow-drift and weak-CV failure map.

NO scene, and NO seed, crosses banks.  Seed bases are disjoint from every prior
cohort in the repo (630900 detail24-dev, 631000 detail24-conf, 631900/632000
detail32, 632900/633000 M1, 650000 R63-bridge DEV -- the 6 DEV scenes used by the
feasibility probe).  New disjoint bases:

  confirmatory : 700000 + 100*stratum + instance
  calibration  : 710000 + 100*stratum + instance
  edge         : 720000 + 100*stratum + instance

STRATA (6; "2 each" for the confirmatory cohort per R34 -> >=12 scenes):
  contour   frozen detail24._gen_contour        (verbatim import)
  texture   frozen detail24._gen_microtexture   (verbatim import)
  twopop    frozen gen_bridge_scenes.gen_twopop (verbatim import; R23 T2-B witness)
  smooth    NEW here: low-pass Gaussian random field (low spatial frequency)
  binsparse NEW here: sparse bright features on a modest background (binary-like)
  natural   NEW here: 1/f power-law random field = synthetic natural-image-statistics
            proxy (NOT a photographic crop; labelled honestly as a 1/f proxy)

Every image is 32x32 float64 in [0,1] with values chosen to keep the cohort mean
intensity in the bridge-scene band (~0.15-0.55) so that under the FROZEN forward
PHI (2200 nominal counts/exposure, cohort-average) the per-scene detected-count
levels match the feasibility probe's regime (no per-scene brightness renormal;
frozen PHI is reused verbatim by the estimator -- see dlgi_forward.py).

Reused generators are IMPORTED verbatim (zero code drift); the three NEW generators
live here and are archived with the bank.  Determinism: every image is a pure
function of its integer seed (numpy default_rng only); regeneration is
byte-reproducible and self-checked in __main__.

Writes ONLY under results/round63_dlgi_campaign/scene_banks/<bank>/.  CPU.
"""
import hashlib
import json
import os
import sys
from datetime import datetime, timezone

import numpy as np
from scipy.ndimage import gaussian_filter

REPO = r"D:/GI_another"
CODE_R63 = os.path.join(REPO, "code", "round63")
sys.path.insert(0, CODE_R63)
import detail24 as d24                     # frozen family generators (verbatim, import-only)
import gen_bridge_scenes as gbs            # frozen twopop witness (verbatim, import-only)

SIDE = 32
OUT_ROOT = os.path.join(REPO, "results", "round63_dlgi_campaign", "scene_banks")

# canonical stratum order (index drives the seed offset)
STRATA = ["contour", "texture", "twopop", "smooth", "binsparse", "natural"]

# bank -> (seed_base, instances_per_stratum)
BANKS = {
    "calibration": dict(base=710000, n_inst=2),
    "confirmatory": dict(base=700000, n_inst=2),   # 6 x 2 = 12 (>= R34 minimum)
    "edge": dict(base=720000, n_inst=1),
}
# two-population brightness ratio R per instance index (brackets the nominal ~5-10)
TWOPOP_R = [8, 16, 4, 32]


# ------------------------------------------------------------------ statistics
def img_stats(x):
    x = np.asarray(x, np.float64)
    mean = float(x.mean()); std = float(x.std())
    mn, mx = float(x.min()), float(x.max())
    F = np.fft.fftshift(np.abs(np.fft.fft2(x - mean)))
    n = x.shape[0]
    cc, rr = np.meshgrid(np.arange(n) - n // 2, np.arange(n) - n // 2)
    rad = np.sqrt(rr ** 2 + cc ** 2); p = F ** 2
    spec = float((rad * p).sum() / p.sum() / (n / 2.0)) if p.sum() > 0 else 0.0
    return dict(mean=mean, std=std, min=mn, max=mx,
                rms_contrast=float(std / mean) if mean > 0 else 0.0,
                michelson=float((mx - mn) / (mx + mn)) if (mx + mn) > 0 else 0.0,
                spectral_centroid=spec)


# ------------------------------------------------- reused frozen generators ---
def gen_contour(seed):
    x = d24._generate01("contour", int(seed), SIDE)
    return np.asarray(x, np.float64), {
        "generator": "code/round63/detail24.py::_gen_contour via _generate01 "
                     "(FROZEN, verbatim import; fresh seed)",
        "family": "contour", "seed": int(seed),
        "params_snapshot": d24.PARAMS["contour"]}


def gen_texture(seed):
    x = d24._generate01("microtexture", int(seed), SIDE)
    # record realised DoG sigma-pair (family scale knob), replicating the draw order
    rng = np.random.default_rng(int(seed))
    rng.integers(0, 2, size=(SIDE, SIDE))
    pairs = d24.PARAMS["microtexture"]["sigma_pairs"]
    idx = int(rng.integers(0, len(pairs)))
    return np.asarray(x, np.float64), {
        "generator": "code/round63/detail24.py::_gen_microtexture via _generate01 "
                     "(FROZEN, verbatim import; fresh seed)",
        "family": "microtexture", "seed": int(seed),
        "realized_sigma_pair": pairs[idx], "realized_sigma_pair_index": idx,
        "params_snapshot": d24.PARAMS["microtexture"]}


def gen_twopop(seed, R):
    x, p = gbs.gen_twopop(SIDE, int(R), int(seed))   # frozen R23 T2-B witness, verbatim
    p = dict(p)
    p["generator"] = ("code/round63/gen_bridge_scenes.py::gen_twopop "
                      "(FROZEN, verbatim import; fresh seed)")
    return np.asarray(x, np.float64), p


# ----------------------------------------------------------- NEW generators ---
def gen_smooth(seed):
    """Low-pass Gaussian random field: energy concentrated at low spatial
    frequency (a smooth-scene stratum).  Pure function of seed."""
    rng = np.random.default_rng(int(seed))
    f = rng.standard_normal((SIDE, SIDE))
    sig = float(rng.uniform(3.0, 5.0))
    sm = gaussian_filter(f, sig, mode="reflect")
    sm = (sm - sm.min()) / (np.ptp(sm) + 1e-12)      # -> [0,1]
    x = 0.06 + 0.90 * sm                             # gentle floor; mean ~0.5
    x = np.clip(x, 0.0, 1.0)
    return x, {"generator": "dlgi_scene_bank.gen_smooth (NEW; low-pass Gaussian "
                            "random field)", "family": "smooth", "seed": int(seed),
               "blur_sigma": sig, "floor": 0.06, "gain": 0.90}


def gen_binsparse(seed):
    """Sparse bright features (1-3 px squares) on a modest background: a
    binary-like sparse stratum.  Pure function of seed."""
    rng = np.random.default_rng(int(seed))
    bg = 0.06
    x = np.full((SIDE, SIDE), bg, np.float64)
    n_feat = int(rng.integers(28, 60))
    feats = []
    for _ in range(n_feat):
        r0 = int(rng.integers(0, SIDE)); c0 = int(rng.integers(0, SIDE))
        w = int(rng.integers(1, 4))                  # 1..3 px square
        val = float(rng.uniform(0.85, 1.0))
        r1 = min(SIDE, r0 + w); c1 = min(SIDE, c0 + w)
        x[r0:r1, c0:c1] = val
        feats.append([r0, c0, w, round(val, 4)])
    x = np.clip(x, 0.0, 1.0)
    return x, {"generator": "dlgi_scene_bank.gen_binsparse (NEW; sparse bright "
                            "features on a modest background)", "family": "binsparse",
               "seed": int(seed), "bg_level": bg, "n_features": n_feat,
               "feature_fill_frac": float((x > bg).mean())}


def gen_natural(seed):
    """1/f^beta power-law random field: a synthetic proxy for natural-image
    spectral statistics (natural images have ~1/f amplitude spectra).  NOT a
    photographic crop -- an honest 1/f proxy, self-contained and reproducible.
    Pure function of seed."""
    rng = np.random.default_rng(int(seed))
    w = rng.standard_normal((SIDE, SIDE))
    F = np.fft.fft2(w)
    fy = np.fft.fftfreq(SIDE)[:, None]; fx = np.fft.fftfreq(SIDE)[None, :]
    fr = np.sqrt(fy ** 2 + fx ** 2); fr[0, 0] = 1.0
    beta = float(rng.uniform(0.9, 1.4))
    img = np.real(np.fft.ifft2(F / fr ** beta))
    img = (img - img.min()) / (np.ptp(img) + 1e-12)
    x = np.clip(0.05 + 0.90 * img, 0.0, 1.0)         # mean ~0.5
    return x, {"generator": "dlgi_scene_bank.gen_natural (NEW; 1/f^beta power-law "
                            "random field = synthetic natural-image-statistics proxy, "
                            "NOT a photographic crop)", "family": "natural",
               "seed": int(seed), "beta": beta, "floor": 0.05, "gain": 0.90}


def generate(stratum, seed, inst):
    if stratum == "contour":
        return gen_contour(seed)
    if stratum == "texture":
        return gen_texture(seed)
    if stratum == "twopop":
        return gen_twopop(seed, TWOPOP_R[inst % len(TWOPOP_R)])
    if stratum == "smooth":
        return gen_smooth(seed)
    if stratum == "binsparse":
        return gen_binsparse(seed)
    if stratum == "natural":
        return gen_natural(seed)
    raise ValueError(stratum)


# ------------------------------------------------------------------- writers
def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def bank_table(bank):
    """[(scene_id, stratum, instance, seed)] for one bank, canonical order."""
    spec = BANKS[bank]; base = spec["base"]; n_inst = spec["n_inst"]
    out = []
    for si, stratum in enumerate(STRATA):
        for i in range(n_inst):
            seed = base + 100 * si + i
            out.append((f"{bank}_{stratum}_{i}", stratum, i, seed))
    return out


def build_bank(bank, write=True):
    out_dir = os.path.join(OUT_ROOT, bank)
    if write:
        os.makedirs(out_dir, exist_ok=True)
    entries = []
    for scene_id, stratum, inst, seed in bank_table(bank):
        x, params = generate(stratum, seed, inst)
        x = np.asarray(x, np.float64)
        assert x.shape == (SIDE, SIDE), (scene_id, x.shape)
        assert np.all(np.isfinite(x)) and x.min() >= 0.0 and x.max() <= 1.0, scene_id
        assert x.var() > 0.0, ("degenerate scene", scene_id)
        x_sha = sha256_bytes(np.ascontiguousarray(x, np.float64).tobytes())
        params = dict(params); params["x_sha256"] = x_sha
        entry = dict(scene_id=scene_id, bank=bank, stratum=stratum, instance=inst,
                     seed=int(seed), gen_params=params, x_sha256=x_sha,
                     stats=img_stats(x))
        if write:
            npz = os.path.join(out_dir, scene_id + ".npz")
            np.savez(npz, x=x)
            with open(npz, "rb") as fh:
                entry["sha256_npz"] = sha256_bytes(fh.read())
            try:
                from imageio.v2 import imwrite
                imwrite(os.path.join(out_dir, scene_id + ".png"),
                        np.round(np.clip(x, 0, 1) * 255).astype(np.uint8))
            except Exception:
                pass
        entries.append(entry)
    manifest = dict(
        campaign="ROUND63 DLGI confirmatory campaign (R34)",
        bank=bank, role={
            "calibration": "construct the critical-value table c_0.95(eta) ONLY",
            "confirmatory": "final endpoints ONLY (fresh scenes/seeds)",
            "edge": "display fast/slow-drift + weak-CV failure map ONLY"}[bank],
        side=SIDE, value_range=[0.0, 1.0], seed_base=BANKS[bank]["base"],
        instances_per_stratum=BANKS[bank]["n_inst"], strata=STRATA,
        n_scenes=len(entries),
        seed_disjoint_from=[630900, 631000, 631900, 632000, 632900, 633000, 650000,
                            "the other two DLGI banks"],
        blind_generation="no reconstruction, arm, or bucket metric computed during generation",
        utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
        python=sys.version.split()[0], numpy=np.__version__,
        scenes=entries)
    if write:
        with open(os.path.join(out_dir, "MANIFEST.json"), "w") as fh:
            json.dump(manifest, fh, indent=2)
    return manifest


def main():
    print("=== DLGI campaign scene banks ===", flush=True)
    all_sha = {}
    for bank in ["calibration", "confirmatory", "edge"]:
        man = build_bank(bank, write=True)
        print(f"\n[{bank}]  {man['n_scenes']} scenes  base={man['seed_base']}", flush=True)
        print("%-26s %-9s %5s  mean   std    rmsC   specC   x_sha256[:12]" %
              ("scene_id", "stratum", "seed"), flush=True)
        for e in man["scenes"]:
            s = e["stats"]
            all_sha[e["scene_id"]] = e["x_sha256"]
            print("%-26s %-9s %6d  %.3f  %.3f  %.3f  %.3f  %s" %
                  (e["scene_id"], e["stratum"], e["seed"], s["mean"], s["std"],
                   s["rms_contrast"], s["spectral_centroid"], e["x_sha256"][:12]),
                  flush=True)

    # -------- determinism + disjointness self-check --------
    print("\n=== determinism (regenerate, no write) + disjointness ===", flush=True)
    ok = True
    seen = {}
    for bank in ["calibration", "confirmatory", "edge"]:
        man2 = build_bank(bank, write=False)
        for e in man2["scenes"]:
            if e["x_sha256"] != all_sha[e["scene_id"]]:
                ok = False; print("  MISMATCH", e["scene_id"], flush=True)
            if e["x_sha256"] in seen:
                ok = False
                print("  DUPLICATE scene content:", e["scene_id"], "==",
                      seen[e["x_sha256"]], flush=True)
            seen[e["x_sha256"]] = e["scene_id"]
    print("determinism + all-distinct across banks:", "PASS" if ok else "FAIL", flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
