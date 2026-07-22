"""T-A: generate + freeze the 16 R63 DEV bridge scenes.

BLIND generation: NO reconstruction, NO arm, NO bucket-signal quality metric is
computed here. Only image-domain statistics (mean/std/contrast/spectrum) are used
for parameter targeting and reporting -- these are simple image statistics, not
reconstruction or bucket evaluations.

Outputs (all inside data/r63_bridge_scenes/ only):
  bridge_{group}_{0..3}.npz   -> field x (32x32 float64, [0,1])
  bridge_{group}_{0..3}.png   -> 32x32 uint8 preview (native resolution)
  BRIDGE_SCENES.json          -> manifest {scene_id, group, gen_params, sha256}
  GENERATION_NOTES.md         -> (written by a separate step)

Generators for contour / microtex / control reuse the FROZEN detail24 family
generators verbatim (code/round63/detail24.py), fresh seeds only. twopop is an
analytic two-population witness per R23 T2-B (docs/ROUND63_GPT_ROUND23_RULING_RAW.md
section 3.2).
"""
import hashlib
import json
import os
import sys

import numpy as np
from scipy.ndimage import gaussian_filter

CODE = r"D:\GI_another\code\round63"
sys.path.insert(0, CODE)
import detail24 as d24  # noqa: E402  (frozen family generators; import only)

OUT = r"D:\GI_another\data\r63_bridge_scenes"
SIDE = 32

# ----------------------------------------------------------------- seeds --- #
# seed base 650000+, disjoint from 632900 DEV and 633000+ confirmatory.
SEEDS = {
    "contour":  [650000, 650001, 650002, 650003],
    "twopop":   [650100, 650101, 650102, 650103],   # R = 4, 8, 16, 32
    "microtex": [650200, 650201, 650202, 650203],
    "control":  [650300, 650301, 650302, 650303],   # maze, maze, glyph, glyph
}
TWOPOP_R = [4, 8, 16, 32]
CONTROL_FAMILY = ["maze", "maze", "glyph", "glyph"]


# ------------------------------------------------------- image statistics -- #
def img_stats(x):
    x = np.asarray(x, dtype=np.float64)
    mean = float(x.mean())
    std = float(x.std())
    mn, mx = float(x.min()), float(x.max())
    rms_contrast = float(std / mean) if mean > 0 else 0.0        # RMS contrast
    michelson = float((mx - mn) / (mx + mn)) if (mx + mn) > 0 else 0.0
    # radial spectral centroid (normalized), a pure image statistic
    F = np.fft.fftshift(np.abs(np.fft.fft2(x - mean)))
    n = x.shape[0]
    cc, rr = np.meshgrid(np.arange(n) - n // 2, np.arange(n) - n // 2)
    rad = np.sqrt(rr ** 2 + cc ** 2)
    p = F ** 2
    spec_centroid = float((rad * p).sum() / p.sum() / (n / 2.0)) if p.sum() > 0 else 0.0
    return dict(mean=mean, std=std, min=mn, max=mx,
                rms_contrast=rms_contrast, michelson=michelson,
                spectral_centroid=spec_centroid)


# ------------------------------------------- analytic two-population (T2-B) - #
def gen_twopop(side, R, seed):
    """R23 Theorem T2-B witness rendered as a 32x32 image in [0,1].

    Population H (bright, low task leverage): a large, spatially smooth bright
    plateau  ->  high mean intensity, negligible fine structure.
    Population L (dim, task-carrying): thin fine-structure strokes at a low
    intensity, spatially disjoint from the plateau so the constructed brightness
    ratio bright_level / dim_level == R exactly.
    """
    rng = np.random.default_rng(int(seed))
    bg = 0.02
    bright_level = 0.80
    dim_level = bright_level / float(R)                # amplitude = bright/R
    rr, cc = np.mgrid[0:side, 0:side]
    rr = rr.astype(np.float64); cc = cc.astype(np.float64)

    # --- bright smooth large region (low fine structure) ---
    cy = side * 0.50 + float(rng.uniform(-3.0, 3.0))
    cx = side * 0.40 + float(rng.uniform(-2.0, 2.0))
    r_disk = side * 0.34
    disk = (((cc - cx) ** 2 + (rr - cy) ** 2) <= r_disk ** 2).astype(np.float64)
    plateau = gaussian_filter(disk, sigma=2.5)         # smooth edges
    plateau = plateau / float(plateau.max())
    img = np.full((side, side), bg, np.float64)
    img = np.maximum(img, bright_level * plateau)      # smooth bright region

    # --- dim fine structure carrying the task (thin strokes, small period) ---
    plateau_mask = plateau > 0.5
    period = int(rng.choice([3, 4]))
    ph_v = int(rng.integers(0, period))
    ph_h = int(rng.integers(0, period))
    ph_d = int(rng.integers(0, period + 2))
    vlines = ((cc.astype(int) + ph_v) % period) == 0
    hlines = ((rr.astype(int) + ph_h) % period) == 0
    dlines = (((rr + cc).astype(int) + ph_d) % (period + 2)) == 0
    strokes = vlines | hlines | dlines
    task = strokes & (~plateau_mask)                   # disjoint from plateau
    img = np.where(task, dim_level, img)               # set strokes to dim level
    img = np.clip(img, 0.0, 1.0)

    core = plateau > 0.9
    bright_mean = float(img[core].mean()) if core.any() else float("nan")
    dim_mean = float(img[task].mean()) if task.any() else float("nan")
    params = {
        "recipe": "R23 Theorem T2-B two-population witness (docs/"
                  "ROUND63_GPT_ROUND23_RULING_RAW.md section 3.2)",
        "seed": int(seed),
        "side": int(side),
        "brightness_ratio_R": int(R),
        "bg_level": bg,
        "bright_level": bright_level,
        "dim_level": dim_level,
        "fine_amp_frac_of_bright": dim_level / bright_level,   # == 1/R
        "disk_radius_frac": 0.34,
        "plateau_blur_sigma": 2.5,
        "disk_center_row_col": [round(cy, 4), round(cx, 4)],
        "stroke_period_px": period,
        "stroke_phases_vhd": [ph_v, ph_h, ph_d],
        "stroke_region": "complement of plateau (plateau<=0.5); disjoint",
        "plateau_core_mean_measured": bright_mean,
        "dim_structure_mean_measured": dim_mean,
        "brightness_ratio_measured": (bright_mean / dim_mean
                                      if dim_mean and dim_mean > 0 else None),
    }
    return img, params


# ------------------------------------------- frozen-generator families ----- #
def gen_family(family, seed, side):
    """Reuse a FROZEN detail24 family generator verbatim (import-only), fresh seed.
    _generate01 is a pure function of (family, seed, side); it writes nothing."""
    img = d24._generate01(family, int(seed), int(side))
    params = {
        "generator": "code/round63/detail24.py :: _GEN[%r] via _generate01 "
                     "(FROZEN, verbatim; fresh seed)" % family,
        "family": family,
        "seed": int(seed),
        "side": int(side),
        "params_ref": "detail24.PARAMS[%r] (unchanged)" % family,
        "params_snapshot": d24.PARAMS[family],
    }
    if family == "microtexture":
        # replicate the generator's draw order to record the realized DoG sigma
        # pair (the family's low-amplitude / scale knob). Pure, deterministic.
        rng = np.random.default_rng(int(seed))
        rng.integers(0, 2, size=(int(side), int(side)))          # field draw
        pairs = d24.PARAMS["microtexture"]["sigma_pairs"]
        idx = int(rng.integers(0, len(pairs)))
        params["realized_sigma_pair"] = pairs[idx]
        params["realized_sigma_pair_index"] = idx
    return np.asarray(img, dtype=np.float64), params


# ---------------------------------------------------------------- writers -- #
def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def write_scene(scene_id, group, x, gen_params):
    from imageio.v2 import imwrite
    x = np.asarray(x, dtype=np.float64)
    assert x.shape == (SIDE, SIDE), (scene_id, x.shape)
    assert np.all(np.isfinite(x)) and x.min() >= 0.0 and x.max() <= 1.0, scene_id
    npz_path = os.path.join(OUT, scene_id + ".npz")
    png_path = os.path.join(OUT, scene_id + ".png")
    # reproducible content hash of the array (independent of zip timestamps)
    x_sha = sha256_bytes(np.ascontiguousarray(x, dtype=np.float64).tobytes())
    gen_params = dict(gen_params)
    gen_params["x_sha256"] = x_sha
    np.savez(npz_path, x=x)
    u8 = np.round(np.clip(x, 0.0, 1.0) * 255.0).astype(np.uint8)
    imwrite(png_path, u8)
    with open(npz_path, "rb") as fh:
        npz_sha = sha256_bytes(fh.read())
    stats = img_stats(x)
    return {
        "scene_id": scene_id,
        "group": group,
        "gen_params": gen_params,
        "sha256": npz_sha,             # sha256 of the .npz file bytes
        "stats": stats,
    }


def main():
    os.makedirs(OUT, exist_ok=True)
    entries = []

    # contour (frozen contour generator, fresh seeds)
    for i, sd in enumerate(SEEDS["contour"]):
        x, p = gen_family("contour", sd, SIDE)
        entries.append(write_scene("bridge_contour_%d" % i, "contour", x, p))

    # twopop (analytic T2-B; R sweep)
    for i, (sd, R) in enumerate(zip(SEEDS["twopop"], TWOPOP_R)):
        x, p = gen_twopop(SIDE, R, sd)
        entries.append(write_scene("bridge_twopop_%d" % i, "twopop", x, p))

    # microtex (frozen microtexture generator, fresh seeds)
    for i, sd in enumerate(SEEDS["microtex"]):
        x, p = gen_family("microtexture", sd, SIDE)
        entries.append(write_scene("bridge_microtex_%d" % i, "microtex", x, p))

    # control (2 maze + 2 glyph, frozen generators, fresh seeds, std params)
    for i, (sd, fam) in enumerate(zip(SEEDS["control"], CONTROL_FAMILY)):
        x, p = gen_family(fam, sd, SIDE)
        entries.append(write_scene("bridge_control_%d" % i, "control", x, p))

    manifest = {
        "cohort": "ROUND63 DEV bridge scenes (T-A)",
        "spec": ["docs/ROUND63_BRIDGE_BUILD_PLAN.md (interface contract)",
                 "docs/ROUND63_GPT_ROUND24_RULING_RAW.md section 3.1",
                 "docs/ROUND63_GPT_ROUND23_RULING_RAW.md section 3.2 (T2-B)"],
        "side": SIDE,
        "value_range": [0.0, 1.0],
        "n_scenes": len(entries),
        "seed_base": 650000,
        "seed_disjoint_from": [632900, 633000],
        "blind_generation": ("no reconstruction, no arm, no bucket-signal quality "
                             "metric was computed during generation"),
        "sha256_field": "sha256 of the .npz file bytes; gen_params.x_sha256 is the "
                        "reproducible content hash of x.tobytes()",
        "scenes": entries,
    }
    with open(os.path.join(OUT, "BRIDGE_SCENES.json"), "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=False)

    # -------- M1 contour-family contrast reference (in-memory, no disk writes) --
    # For the "contour contrast comparable to M1 contour family" verification.
    m1_contour_seeds = [633400, 633401, 633402, 633403, 632904]  # conf 0..3 + dev
    m1_ref = []
    for sd in m1_contour_seeds:
        xm = d24._generate01("contour", sd, SIDE)   # pure fn, writes nothing
        m1_ref.append((sd, img_stats(xm)))

    # -------- print report --------
    print("=== BRIDGE SCENES (16) ===", flush=True)
    hdr = ("%-20s %-9s  mean    std     min    max    rmsC    michel  specC" % ("scene_id", "group"))
    print(hdr, flush=True)
    for e in entries:
        s = e["stats"]
        print("%-20s %-9s  %.4f  %.4f  %.3f  %.3f  %.4f  %.4f  %.4f"
              % (e["scene_id"], e["group"], s["mean"], s["std"], s["min"],
                 s["max"], s["rms_contrast"], s["michelson"],
                 s["spectral_centroid"]), flush=True)

    print("\n=== twopop measured brightness ratios ===", flush=True)
    for e in entries:
        if e["group"] == "twopop":
            g = e["gen_params"]
            print("%-20s R=%2d  bright_core=%.4f dim=%.4f  measured_ratio=%.2f"
                  % (e["scene_id"], g["brightness_ratio_R"],
                     g["plateau_core_mean_measured"], g["dim_structure_mean_measured"],
                     g["brightness_ratio_measured"]), flush=True)

    print("\n=== microtex realized sigma pairs (low-amplitude scale knob) ===", flush=True)
    for e in entries:
        if e["group"] == "microtex":
            g = e["gen_params"]
            print("%-20s sigma_pair=%s (idx %d)  mean=%.4f"
                  % (e["scene_id"], g["realized_sigma_pair"],
                     g["realized_sigma_pair_index"], e["stats"]["mean"]), flush=True)

    print("\n=== M1 contour family reference (in-memory, float generator) ===", flush=True)
    print("%-10s  mean    std     rmsC    michel  specC" % "seed", flush=True)
    for sd, s in m1_ref:
        print("%-10d  %.4f  %.4f  %.4f  %.4f  %.4f"
              % (sd, s["mean"], s["std"], s["rms_contrast"], s["michelson"],
                 s["spectral_centroid"]), flush=True)
    bc = [e["stats"]["rms_contrast"] for e in entries if e["group"] == "contour"]
    mc = [s["rms_contrast"] for _, s in m1_ref]
    print("\nbridge contour rms_contrast range: [%.4f, %.4f]" % (min(bc), max(bc)), flush=True)
    print("M1     contour rms_contrast range: [%.4f, %.4f]" % (min(mc), max(mc)), flush=True)

    # emit reference for the notes file (scratchpad only; NOT in the frozen dir)
    scratch = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(scratch, "_m1_contour_ref.json"), "w") as fh:
        json.dump({"m1_contour_reference": [{"seed": sd, "stats": s}
                                            for sd, s in m1_ref]}, fh, indent=2)
    print("\n[gen] wrote %d scenes + manifest to %s" % (len(entries), OUT), flush=True)


if __name__ == "__main__":
    main()
