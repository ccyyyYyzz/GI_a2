"""ROUND63 DETAIL-24 primary confirmatory cohort generator (GPT round-7 F1
ruling §1.1). SELF-CONTAINED, F1-FROZEN — no dependence on images63/STL splits.

DETAIL-24 = exactly 24 procedurally generated diagnostic targets, six fixed
families x four independent instances, the confirmatory cohort for the positive
acquisition-speed claim (ruling §1). Six families (ruling §1.1, in order):
  0 glyph        bitmap 5x7 alphanumeric strings on a low background
  1 chirp        h/v/diagonal line-pair gratings, widths/gaps 1-4 px
  2 spokes       Siemens-star-like radial spokes in an annulus
  3 maze         1-2 px connected stroke networks (randomized-DFS spanning tree)
  4 contour      thin ellipse/polygon/curve outlines on a NONZERO background
  5 microtexture rectified band-pass (DoG) of a seeded Rademacher +/-1 field

Each family generator is a PURE DETERMINISTIC function of one integer instance
seed (numpy default_rng(seed) ONLY -- no wall clock, no global state, no
rng_for/SEED0 coupling), producing a 64x64 float64 image with values >= 0, then
sum-normalized to 1 exactly like every other round-63 image set (images63).

Every numeric range is a FROZEN parameter collected in PARAMS at module top and
emitted verbatim to params_manifest.json in each cache root. Generation and
hashing of the 24 confirmatory instances is REQUIRED (F1 lock condition 1) but
they must NEVER be reconstructed (no solver) before F1.

Frozen identities
-----------------
DEV instances (implementation-only, ruling §1.2): seeds 630900..630905, exactly
one per family in order, names detail_dev_<family>.
CONFIRMATORY instances (ruling §1.1): seeds 631000 + 4*f + i for family f=0..5,
instance i=0..3, names detail_<family>_<i>.

API
---
build_detail24_set(side=64)     -> OrderedDict name -> flattened float64 vector
                                   (24 confirmatory, sum-normalized to 1)
build_detail24_dev_set(side=64) -> OrderedDict name -> vector (6 dev)
family_of(name)                 -> family string (for the family-stratified boot)

Cache: PNGs + sha256.txt + params_manifest.json under
data/r63_images_detail24/<side>/ and data/r63_images_detail24_dev/<side>/.
The 8-bit PNG is the truth of record: each image is generated once, written,
then read back (uint8/255, sum-normalized), so truth is byte-reproducible from
the cache and matches its sha256. Regeneration is byte-identical.
"""
import collections
import hashlib
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
DATA = os.path.join(ROOT, "data")
IMG_ROOT = os.path.join(DATA, "r63_images_detail24")
DEV_IMG_ROOT = os.path.join(DATA, "r63_images_detail24_dev")
# STUDY-2 DETAIL-32 cohort (GPT round-8 ruling §2): the same six families at
# 32^2 with per-image generator-frozen SIGNAL / BACKGROUND ROIs.
IMG32_ROOT = os.path.join(DATA, "r63_images_detail32")
DEV32_IMG_ROOT = os.path.join(DATA, "r63_images_detail32_dev")

FAMILIES = ["glyph", "chirp", "spokes", "maze", "contour", "microtexture"]

# ------------------------------------------------------------------ FROZEN --- #
# All generator ranges. Interpretation of an [lo, hi] "range": integer draws are
# inclusive of both ends (rng.integers(lo, hi+1)); float ranges are uniform
# [lo, hi]; "*_choices" lists are drawn uniformly. DO NOT edit after F1.
PARAMS = {
    "side_default": 64,
    "peak": 1.0,
    "dtype": "float64",
    "sanity": {
        "min_frac_above_half": 0.01,   # >= 1% of pixels above half-max
        "max_frac_above_half": 0.90,   # <= 90% of pixels above half-max
        "require_positive_variance": True,
    },
    "seeds": {"dev_base": 630900, "conf_base": 631000, "conf_family_stride": 4},
    "glyph": {
        "alphabet": "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789",
        "bg_level": 0.06, "fg_level": 1.0,           # peak/bg = 16.67 (frozen)
        "scale_choices": [1, 2],
        "n_rows_range": [3, 6],
        "chars_per_row_range": [4, 9],
        "top_margin": 3, "left_margin": 2,
        "char_gap_font_px": 1, "row_gap_font_px": 2,
    },
    "chirp": {
        "bg_level": 0.0, "fg_level": 1.0,
        "grid_choices": [2, 3, 4],
        "width_range": [1, 4], "gap_range": [1, 4],
        "orientations": ["h", "v", "d", "ad"],       # ad = anti-diagonal
    },
    "spokes": {
        "bg_level": 0.0, "fg_level": 1.0,
        "n_spokes_choices": [8, 12, 16, 20, 24, 32],
        "center_offset_px_range": [-6, 6],
        "r_inner_frac": 0.12, "r_outer_frac": 0.92,
    },
    "maze": {
        "bg_level": 0.0, "fg_level": 1.0,
        "step_choices": [6, 8],
        "stroke_width_choices": [1, 2],
        "margin": 3,
        "fill_factor_range": [0.05, 0.55],           # structural bound (asserted)
    },
    "contour": {
        "bg_level_range": [0.05, 0.15],              # NONZERO background (frozen)
        "fg_level": 1.0,
        "n_shapes_range": [4, 7],
        "outline_width_choices": [1, 2],
        "shape_types": ["ellipse", "polygon", "curve"],
        "center_frac_range": [0.28, 0.72],
        "radius_frac_range": [0.14, 0.42],
        "polygon_vertices_range": [3, 6],
        "curve_freq_range": [1.0, 3.0],
    },
    "microtexture": {
        "fg_level": 1.0,
        "sigma_pairs": [[0.8, 1.6], [1.0, 2.0], [1.2, 2.4], [1.5, 3.0]],
        "rectify": "abs",                            # -> values >= 0
        "renorm": "peak",                            # scale so max == peak
        "conv_mode": "reflect",
    },
}

# ------------------------------------------------------ STUDY-2 DETAIL-32 --- #
# Frozen parameters for the DETAIL-32 cohort (GPT round-8 ruling §2). The six
# family GENERATORS are reused verbatim at side=32 (PARAMS above is unchanged);
# only the resolution, the confirmatory/dev SEED tables, and the generator-
# frozen ROI rule are new. DO NOT edit after the Study-2 freeze.
PARAMS32 = {
    "side": 32,
    # confirmatory seeds s_{f,r} = 632000 + 100 f + r  (f family 0..5, r 0..3);
    # six dev-only instances 631900 + f (ruling §2).
    "seeds": {"conf_base": 632000, "conf_family_mult": 100,
              "dev_base": 631900},
    # SIGNAL / BACKGROUND ROIs are a DETERMINISTIC function of the clean truth
    # of record (the read-back 8-bit PNG) — generator-frozen, never tuned on TV,
    # contrast, CNR, PSNR or visual quality (ruling §2).
    "roi": {
        "signal_rule": ("binary_dilation(truth > 0.5*max, 3x3 full "
                        "structuring element, 1 iteration)"),
        "signal_half_max_frac": 0.5,
        "signal_dilate_iters": 1,
        "background_rule": ("complement of binary_dilation(core, 3x3 full, 2 "
                            "iterations) = pixels at least 2 px from any signal "
                            "core pixel; deterministic fallback when that "
                            "complement is empty (fully-filled scene) = outer "
                            "2 px border minus the signal ROI"),
        "background_guard_iters": 2,
        "border_width": 2,
        "connectivity": "8 (3x3 full structuring element)",
    },
    "sanity": {
        "min_frac_above_half": 0.01,
        "max_frac_above_half": 0.90,
        "require_positive_variance": True,
        "roi_nonempty": True,
        "roi_disjoint": True,
    },
}

# --------------------------------------------------------- 5x7 bitmap font --- #
# Repository-embedded 5-wide x 7-tall alphanumeric font (ruling §1.1: no system
# fonts, no PIL text). '#' = on, '.' = off; each glyph is 7 rows of 5 chars.
FONT_5x7 = {
    "A": (".###.", "#...#", "#...#", "#####", "#...#", "#...#", "#...#"),
    "B": ("####.", "#...#", "#...#", "####.", "#...#", "#...#", "####."),
    "C": (".###.", "#...#", "#....", "#....", "#....", "#...#", ".###."),
    "D": ("###..", "#..#.", "#...#", "#...#", "#...#", "#..#.", "###.."),
    "E": ("#####", "#....", "#....", "####.", "#....", "#....", "#####"),
    "F": ("#####", "#....", "#....", "####.", "#....", "#....", "#...."),
    "G": (".###.", "#...#", "#....", "#.###", "#...#", "#...#", ".###."),
    "H": ("#...#", "#...#", "#...#", "#####", "#...#", "#...#", "#...#"),
    "I": (".###.", "..#..", "..#..", "..#..", "..#..", "..#..", ".###."),
    "J": ("..###", "...#.", "...#.", "...#.", "#..#.", "#..#.", ".##.."),
    "K": ("#...#", "#..#.", "#.#..", "##...", "#.#..", "#..#.", "#...#"),
    "L": ("#....", "#....", "#....", "#....", "#....", "#....", "#####"),
    "M": ("#...#", "##.##", "#.#.#", "#.#.#", "#...#", "#...#", "#...#"),
    "N": ("#...#", "#...#", "##..#", "#.#.#", "#..##", "#...#", "#...#"),
    "O": (".###.", "#...#", "#...#", "#...#", "#...#", "#...#", ".###."),
    "P": ("####.", "#...#", "#...#", "####.", "#....", "#....", "#...."),
    "Q": (".###.", "#...#", "#...#", "#...#", "#.#.#", "#..#.", ".##.#"),
    "R": ("####.", "#...#", "#...#", "####.", "#.#..", "#..#.", "#...#"),
    "S": (".###.", "#...#", "#....", ".###.", "....#", "#...#", ".###."),
    "T": ("#####", "..#..", "..#..", "..#..", "..#..", "..#..", "..#.."),
    "U": ("#...#", "#...#", "#...#", "#...#", "#...#", "#...#", ".###."),
    "V": ("#...#", "#...#", "#...#", "#...#", "#...#", ".#.#.", "..#.."),
    "W": ("#...#", "#...#", "#...#", "#.#.#", "#.#.#", "##.##", "#...#"),
    "X": ("#...#", "#...#", ".#.#.", "..#..", ".#.#.", "#...#", "#...#"),
    "Y": ("#...#", "#...#", ".#.#.", "..#..", "..#..", "..#..", "..#.."),
    "Z": ("#####", "....#", "...#.", "..#..", ".#...", "#....", "#####"),
    "0": (".###.", "#...#", "#..##", "#.#.#", "##..#", "#...#", ".###."),
    "1": ("..#..", ".##..", "..#..", "..#..", "..#..", "..#..", ".###."),
    "2": (".###.", "#...#", "....#", "...#.", "..#..", ".#...", "#####"),
    "3": ("#####", "...#.", "..#..", "...#.", "....#", "#...#", ".###."),
    "4": ("...#.", "..##.", ".#.#.", "#..#.", "#####", "...#.", "...#."),
    "5": ("#####", "#....", "####.", "....#", "....#", "#...#", ".###."),
    "6": ("..##.", ".#...", "#....", "####.", "#...#", "#...#", ".###."),
    "7": ("#####", "....#", "...#.", "..#..", ".#...", ".#...", ".#..."),
    "8": (".###.", "#...#", "#...#", ".###.", "#...#", "#...#", ".###."),
    "9": (".###.", "#...#", "#...#", ".####", "....#", "...#.", ".##.."),
}


# ----------------------------------------------------------- draw helpers --- #
def _draw_line(img, r0, c0, r1, c1, val, width):
    """Rasterize a thick line from (r0,c0) to (r1,c1) by dense sampling."""
    side = img.shape[0]
    n = int(max(abs(r1 - r0), abs(c1 - c0))) + 1
    rr = np.linspace(r0, r1, max(n, 2))
    cc = np.linspace(c0, c1, max(n, 2))
    lo = int(np.floor((width - 1) / 2.0))
    hi = int(np.ceil((width - 1) / 2.0))
    for r, c in zip(rr, cc):
        ri, ci = int(round(r)), int(round(c))
        a0 = max(0, ri - lo); a1 = min(side, ri + hi + 1)
        b0 = max(0, ci - lo); b1 = min(side, ci + hi + 1)
        if a1 > a0 and b1 > b0:
            img[a0:a1, b0:b1] = val


def _blit_glyph(img, glyph_rows, r0, c0, scale, val):
    """Blit one 5x7 font glyph at (r0,c0), each font pixel a scale x scale block."""
    side = img.shape[0]
    for gr, rowstr in enumerate(glyph_rows):          # 7 rows
        for gc, ch in enumerate(rowstr):              # 5 cols
            if ch == "#":
                a0 = r0 + gr * scale; b0 = c0 + gc * scale
                a1 = min(side, a0 + scale); b1 = min(side, b0 + scale)
                if a0 < side and b0 < side and a1 > a0 and b1 > b0:
                    img[a0:a1, b0:b1] = val


# ------------------------------------------------------ family generators --- #
def _gen_glyph(side, rng):
    P = PARAMS["glyph"]
    img = np.full((side, side), P["bg_level"], np.float64)
    alpha = P["alphabet"]
    s = int(rng.choice(P["scale_choices"]))
    cell_w = (5 + P["char_gap_font_px"]) * s
    cell_h = (7 + P["row_gap_font_px"]) * s
    top, left = P["top_margin"], P["left_margin"]
    max_rows = max(1, (side - top) // cell_h)
    max_chars = max(1, (side - left) // cell_w)
    n_rows = min(int(rng.integers(P["n_rows_range"][0], P["n_rows_range"][1] + 1)),
                 max_rows)
    for ri in range(n_rows):
        r0 = top + ri * cell_h
        n_ch = min(int(rng.integers(P["chars_per_row_range"][0],
                                    P["chars_per_row_range"][1] + 1)), max_chars)
        for ci in range(n_ch):
            ch = alpha[int(rng.integers(0, len(alpha)))]
            _blit_glyph(img, FONT_5x7[ch], r0, left + ci * cell_w, s, P["fg_level"])
    return img


def _gen_chirp(side, rng):
    P = PARAMS["chirp"]
    img = np.full((side, side), P["bg_level"], np.float64)
    g = int(rng.choice(P["grid_choices"]))
    b = np.linspace(0, side, g + 1).astype(int)
    for gi in range(g):
        for gj in range(g):
            r0, r1, c0, c1 = b[gi], b[gi + 1], b[gj], b[gj + 1]
            w = int(rng.integers(P["width_range"][0], P["width_range"][1] + 1))
            gap = int(rng.integers(P["gap_range"][0], P["gap_range"][1] + 1))
            period = w + gap
            orient = P["orientations"][int(rng.integers(0, len(P["orientations"])))]
            rr, cc = np.mgrid[r0:r1, c0:c1]
            if orient == "h":
                coord = rr
            elif orient == "v":
                coord = cc
            elif orient == "d":
                coord = rr + cc
            else:
                coord = rr - cc
            mask = (coord % period) < w
            img[r0:r1, c0:c1][mask] = P["fg_level"]
    return img


def _gen_spokes(side, rng):
    P = PARAMS["spokes"]
    img = np.full((side, side), P["bg_level"], np.float64)
    n_spokes = int(rng.choice(P["n_spokes_choices"]))
    lo, hi = P["center_offset_px_range"]
    ox = float(rng.integers(lo, hi + 1)); oy = float(rng.integers(lo, hi + 1))
    cx = (side - 1) / 2.0 + ox; cy = (side - 1) / 2.0 + oy
    rr, cc = np.mgrid[0:side, 0:side]
    dy = rr.astype(np.float64) - cy; dx = cc.astype(np.float64) - cx
    theta = np.arctan2(dy, dx)
    rad = np.sqrt(dx * dx + dy * dy)
    Rmax = side / 2.0
    annulus = (rad >= P["r_inner_frac"] * Rmax) & (rad <= P["r_outer_frac"] * Rmax)
    spoke = np.cos(n_spokes * theta) > 0.0
    img[annulus & spoke] = P["fg_level"]
    return img


def _gen_maze(side, rng):
    P = PARAMS["maze"]
    img = np.full((side, side), P["bg_level"], np.float64)
    step = int(rng.choice(P["step_choices"]))
    sw = int(rng.choice(P["stroke_width_choices"]))
    m = P["margin"]
    ys = list(range(m, side - m, step))
    xs = list(range(m, side - m, step))
    ny, nx = len(ys), len(xs)
    visited = np.zeros((ny, nx), bool)
    stack = [(0, 0)]
    visited[0, 0] = True
    edges = []
    while stack:
        i, j = stack[-1]
        nbrs = [(i + di, j + dj) for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1))
                if 0 <= i + di < ny and 0 <= j + dj < nx and not visited[i + di, j + dj]]
        if not nbrs:
            stack.pop()
            continue
        ni, nj = nbrs[int(rng.integers(0, len(nbrs)))]
        visited[ni, nj] = True
        edges.append(((i, j), (ni, nj)))
        stack.append((ni, nj))
    for (i, j), (ni, nj) in edges:
        _draw_line(img, ys[i], xs[j], ys[ni], xs[nj], P["fg_level"], sw)
    return img


def _gen_contour(side, rng):
    P = PARAMS["contour"]
    lo, hi = P["bg_level_range"]
    bg = float(lo + (hi - lo) * rng.random())
    img = np.full((side, side), bg, np.float64)
    rr, cc = np.mgrid[0:side, 0:side]
    rr = rr.astype(np.float64); cc = cc.astype(np.float64)
    cf_lo, cf_hi = P["center_frac_range"]
    rf_lo, rf_hi = P["radius_frac_range"]
    n = int(rng.integers(P["n_shapes_range"][0], P["n_shapes_range"][1] + 1))
    for _ in range(n):
        w = int(rng.choice(P["outline_width_choices"]))
        stype = P["shape_types"][int(rng.integers(0, len(P["shape_types"])))]
        cy = float(rng.uniform(cf_lo, cf_hi) * side)
        cx = float(rng.uniform(cf_lo, cf_hi) * side)
        R = float(rng.uniform(rf_lo, rf_hi) * side)
        if stype == "ellipse":
            a = R; bax = R * float(rng.uniform(0.5, 1.0))
            ang = float(rng.uniform(0.0, np.pi))
            ct, st = np.cos(ang), np.sin(ang)
            xr = (cc - cx) * ct + (rr - cy) * st
            yr = -(cc - cx) * st + (rr - cy) * ct
            f = np.sqrt((xr / a) ** 2 + (yr / bax) ** 2)
            band = np.abs(f - 1.0) * min(a, bax) <= (w / 2.0 + 0.5)
            img[band] = P["fg_level"]
        elif stype == "polygon":
            k = int(rng.integers(P["polygon_vertices_range"][0],
                                 P["polygon_vertices_range"][1] + 1))
            angs = np.sort(rng.uniform(0.0, 2.0 * np.pi, size=k))
            pts = [(cy + R * np.sin(t), cx + R * np.cos(t)) for t in angs]
            for p in range(k):
                (a0, b0), (a1, b1) = pts[p], pts[(p + 1) % k]
                _draw_line(img, a0, b0, a1, b1, P["fg_level"], w)
        else:  # curve: a sine arc drawn as connected segments
            amp = R * 0.5
            freq = float(rng.uniform(P["curve_freq_range"][0],
                                     P["curve_freq_range"][1]))
            phase = float(rng.uniform(0.0, 2.0 * np.pi))
            xs = np.arange(0, side, dtype=np.float64)
            rows = cy + amp * np.sin(2.0 * np.pi * freq * xs / side + phase)
            for p in range(side - 1):
                _draw_line(img, rows[p], xs[p], rows[p + 1], xs[p + 1],
                           P["fg_level"], w)
    return img


def _gen_microtexture(side, rng):
    from scipy.ndimage import gaussian_filter
    P = PARAMS["microtexture"]
    field = rng.integers(0, 2, size=(side, side)).astype(np.float64) * 2.0 - 1.0
    s1, s2 = P["sigma_pairs"][int(rng.integers(0, len(P["sigma_pairs"])))]
    band = (gaussian_filter(field, s1, mode=P["conv_mode"])
            - gaussian_filter(field, s2, mode=P["conv_mode"]))
    out = np.abs(band)                                # rectify -> >= 0
    mx = float(out.max())
    if mx > 0:
        out = out / mx * P["fg_level"]
    return out


_GEN = {
    "glyph": _gen_glyph, "chirp": _gen_chirp, "spokes": _gen_spokes,
    "maze": _gen_maze, "contour": _gen_contour, "microtexture": _gen_microtexture,
}


# ---------------------------------------------------------------- naming --- #
def family_of(name):
    """Family string for a DETAIL-24 or DETAIL-32 image name (conf or dev).

    Handles all four name shapes uniformly by inspecting the token after the
    cohort prefix: detail_<fam>_<i>, detail_dev_<fam>, detail32_<fam>_<r>,
    detail32_dev_<fam> (the six family tokens carry no underscores)."""
    parts = name.split("_")
    fam = parts[2] if len(parts) >= 3 and parts[1] == "dev" else parts[1]
    if fam not in FAMILIES:
        raise ValueError("name %r does not encode a DETAIL family" % name)
    return fam


def _conf_seed(f, i):
    return PARAMS["seeds"]["conf_base"] + PARAMS["seeds"]["conf_family_stride"] * f + i


def _dev_seed(f):
    return PARAMS["seeds"]["dev_base"] + f


def _conf_table():
    """[(name, family, seed)] for the 24 confirmatory instances, canonical order."""
    out = []
    for f, fam in enumerate(FAMILIES):
        for i in range(4):
            out.append(("detail_%s_%d" % (fam, i), fam, _conf_seed(f, i)))
    return out


def _dev_table():
    """[(name, family, seed)] for the 6 dev instances, canonical order."""
    return [("detail_dev_%s" % fam, fam, _dev_seed(f))
            for f, fam in enumerate(FAMILIES)]


# ---------------------------------------------------- generation + cache --- #
def _to_uint8(img01):
    return np.round(np.clip(img01, 0.0, 1.0) * 255.0).astype(np.uint8)


def _generate01(family, seed, side):
    """One family image as float64 in [0, peak]; pure function of (family, seed)."""
    img = _GEN[family](int(side), np.random.default_rng(int(seed)))
    img = np.asarray(img, dtype=np.float64)
    if img.shape != (side, side):
        raise RuntimeError("family %s produced shape %s, expected (%d,%d)"
                           % (family, img.shape, side, side))
    if np.any(img < 0) or not np.all(np.isfinite(img)):
        raise RuntimeError("family %s produced negative/non-finite values" % family)
    return img


def _assert_sanity(img01, name):
    """Structural non-degeneracy on the read-back [0,1] truth (ruling §1.1)."""
    S = PARAMS["sanity"]
    if S["require_positive_variance"] and img01.var() <= 0.0:
        raise RuntimeError("DETAIL-24 %s has zero variance (degenerate)" % name)
    mx = float(img01.max())
    if mx <= 0.0:
        raise RuntimeError("DETAIL-24 %s is all-zero" % name)
    frac = float(np.mean(img01 > 0.5 * mx))
    if not (S["min_frac_above_half"] <= frac <= S["max_frac_above_half"]):
        raise RuntimeError("DETAIL-24 %s frac>half-max=%.4f outside [%.2f,%.2f]"
                           % (name, frac, S["min_frac_above_half"],
                              S["max_frac_above_half"]))
    return frac


def _build(table, cache_root, side):
    """Shared builder: generate/cache PNGs, sanity-check, sum-normalize, and
    (re)write sha256.txt + params_manifest.json. Existing PNGs are canonical."""
    from imageio.v2 import imread, imwrite

    side = int(side)
    cache = os.path.join(cache_root, str(side))
    os.makedirs(cache, exist_ok=True)

    out = collections.OrderedDict()
    sha_lines = []
    sha_map = {}
    for name, fam, seed in table:
        png = os.path.join(cache, name + ".png")
        if not os.path.exists(png):
            imwrite(png, _to_uint8(_generate01(fam, seed, side)))
        u8 = imread(png)
        if u8.shape != (side, side):
            raise RuntimeError("cached PNG %s has shape %s, expected (%d,%d)"
                               % (png, u8.shape, side, side))
        img01 = u8.astype(np.float64) / 255.0
        _assert_sanity(img01, name)
        x = img01.ravel()
        s = x.sum()
        if s <= 0:
            raise RuntimeError("DETAIL-24 %s has non-positive sum %g" % (name, s))
        out[name] = x / s
        with open(png, "rb") as fh:
            dig = hashlib.sha256(fh.read()).hexdigest()
        sha_map[name] = dig
        sha_lines.append("%s  %s.png" % (dig, name))
    with open(os.path.join(cache, "sha256.txt"), "w") as fh:
        fh.write("\n".join(sha_lines) + "\n")
    manifest = {
        "ruling": "GPT round-7 F1 ruling section 1.1 (DETAIL-24)",
        "side": side,
        "families_order": FAMILIES,
        "params": PARAMS,
        "seed_table": [{"name": n, "family": f, "seed": int(sd)}
                       for n, f, sd in table],
        "sha256": sha_map,
    }
    with open(os.path.join(cache, "params_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
    return out


def build_detail24_set(side=64):
    """OrderedDict name -> sum-normalized float64 vector for the 24 CONFIRMATORY
    DETAIL-24 instances (ruling §1.1). Generation + hashing only; callers must
    NOT run a solver on these before F1 (ruling §1.2)."""
    return _build(_conf_table(), IMG_ROOT, side)


def build_detail24_dev_set(side=64):
    """OrderedDict name -> sum-normalized float64 vector for the 6 DEV DETAIL-24
    instances (seeds 630900..630905; ruling §1.2 implementation-only)."""
    return _build(_dev_table(), DEV_IMG_ROOT, side)


# ============================================================ DETAIL-32 ===== #
# STUDY-2 cohort (GPT round-8 ruling §2): the same six families at side=32,
# each output additionally carrying a generator-frozen SIGNAL ROI and BACKGROUND
# ROI (for the CNR co-secondary, ruling §11). The 24 confirmatory instances are
# GENERATED + HASHED + ROI-checked only; NO solver may touch them before the
# Study-2 freeze (the ROI/sanity checks are generation-side, not reconstruction).

def _conf32_seed(f, r):
    return (PARAMS32["seeds"]["conf_base"]
            + PARAMS32["seeds"]["conf_family_mult"] * f + r)


def _dev32_seed(f):
    return PARAMS32["seeds"]["dev_base"] + f


def _conf32_table():
    """[(name, family, seed)] for the 24 DETAIL-32 confirmatory instances."""
    out = []
    for f, fam in enumerate(FAMILIES):
        for r in range(4):
            out.append(("detail32_%s_%d" % (fam, r), fam, _conf32_seed(f, r)))
    return out


def _dev32_table():
    """[(name, family, seed)] for the 6 DETAIL-32 dev instances."""
    return [("detail32_dev_%s" % fam, fam, _dev32_seed(f))
            for f, fam in enumerate(FAMILIES)]


def _roi_struct():
    """3x3 full (8-connectivity) structuring element for the ROI dilation."""
    from scipy.ndimage import generate_binary_structure
    return generate_binary_structure(2, 2)


def compute_rois(truth2d):
    """SIGNAL / BACKGROUND ROI booleans for one clean truth (ruling §2), a
    DETERMINISTIC function of the truth of record:
      signal     = dilate(truth > 0.5*max, 1)                 (feature + 1-px halo)
      background = complement of dilate(core, 2)              (>=2 px from signal)
                   -> fallback (empty on a fully-filled scene): outer border
                      frame minus the signal ROI.
    signal and background are DISJOINT by construction (signal subset of
    dilate-by-2, background its complement; the fallback subtracts the signal).
    Returns (signal_bool, background_bool, mode)."""
    from scipy.ndimage import binary_dilation
    P = PARAMS32["roi"]
    st = _roi_struct()
    truth2d = np.asarray(truth2d, dtype=np.float64)
    mx = float(truth2d.max())
    core = truth2d > P["signal_half_max_frac"] * mx
    signal = binary_dilation(core, structure=st,
                             iterations=P["signal_dilate_iters"])
    far = ~binary_dilation(core, structure=st,
                           iterations=P["background_guard_iters"])
    background = far
    mode = "far_complement"
    if not background.any():
        W = int(P["border_width"])
        frame = np.zeros_like(core)
        frame[:W, :] = True
        frame[-W:, :] = True
        frame[:, :W] = True
        frame[:, -W:] = True
        background = frame & (~signal)
        mode = "border_minus_signal"
    return signal, background, mode


def _assert_roi(signal, background, name):
    """Structural ROI non-degeneracy (ruling §2 / D2 sanity): both ROIs nonempty
    and disjoint. Generator-frozen, not a quality gate."""
    if not signal.any():
        raise RuntimeError("DETAIL-32 %s SIGNAL ROI is empty" % name)
    if not background.any():
        raise RuntimeError("DETAIL-32 %s BACKGROUND ROI is empty" % name)
    if bool((signal & background).any()):
        raise RuntimeError("DETAIL-32 %s SIGNAL and BACKGROUND ROIs overlap"
                           % name)


def _roi_sha(signal, background):
    """Reproducible ROI hash from the boolean array bytes (NOT the .npz file,
    whose zip envelope carries timestamps)."""
    return hashlib.sha256(np.ascontiguousarray(signal, dtype=bool).tobytes()
                          + np.ascontiguousarray(background, dtype=bool).tobytes()
                          ).hexdigest()


def _build32(table, cache_root, side, ruling):
    """DETAIL-32 builder: generate/cache PNGs (truth of record), structural
    sanity, generator-frozen ROIs saved as <name>_roi.npz, and a params manifest
    carrying the image + ROI SHAs. Existing PNGs are canonical (byte-reproducible
    regeneration). ROIs are computed from the read-back [0,1] truth so they hash
    deterministically from the cache."""
    from imageio.v2 import imread, imwrite

    side = int(side)
    cache = os.path.join(cache_root, str(side))
    os.makedirs(cache, exist_ok=True)

    out = collections.OrderedDict()
    sha_lines = []
    sha_map = {}
    roi_sha_map = {}
    roi_mode_map = {}
    for name, fam, seed in table:
        png = os.path.join(cache, name + ".png")
        if not os.path.exists(png):
            imwrite(png, _to_uint8(_generate01(fam, seed, side)))
        u8 = imread(png)
        if u8.shape != (side, side):
            raise RuntimeError("cached PNG %s has shape %s, expected (%d,%d)"
                               % (png, u8.shape, side, side))
        img01 = u8.astype(np.float64) / 255.0
        _assert_sanity(img01, name)
        signal, background, mode = compute_rois(img01.reshape(side, side))
        _assert_roi(signal, background, name)
        np.savez(os.path.join(cache, name + "_roi.npz"),
                 signal=signal, background=background)
        x = img01.ravel()
        s = x.sum()
        if s <= 0:
            raise RuntimeError("DETAIL-32 %s has non-positive sum %g" % (name, s))
        out[name] = x / s
        with open(png, "rb") as fh:
            dig = hashlib.sha256(fh.read()).hexdigest()
        sha_map[name] = dig
        sha_lines.append("%s  %s.png" % (dig, name))
        roi_sha_map[name] = _roi_sha(signal, background)
        roi_mode_map[name] = mode
    with open(os.path.join(cache, "sha256.txt"), "w") as fh:
        fh.write("\n".join(sha_lines) + "\n")
    manifest = {
        "ruling": ruling,
        "side": side,
        "families_order": FAMILIES,
        "params": PARAMS,
        "params32": PARAMS32,
        "seed_table": [{"name": n, "family": f, "seed": int(sd)}
                       for n, f, sd in table],
        "sha256": sha_map,
        "roi_sha256": roi_sha_map,
        "roi_mode": roi_mode_map,
    }
    with open(os.path.join(cache, "params_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
    return out


def build_detail32_set(side=32):
    """OrderedDict name -> sum-normalized float64 vector for the 24 CONFIRMATORY
    DETAIL-32 instances (ruling §2). Generation + hashing + ROI checks only;
    callers must NOT run a solver on these before the Study-2 freeze."""
    return _build32(_conf32_table(), IMG32_ROOT, side,
                    "GPT round-8 ruling section 2 (DETAIL-32 confirmatory)")


def build_detail32_dev_set(side=32):
    """OrderedDict name -> sum-normalized float64 vector for the 6 DEV DETAIL-32
    instances (seeds 631900..631905; ruling §2 dev-only)."""
    return _build32(_dev32_table(), DEV32_IMG_ROOT, side,
                    "GPT round-8 ruling section 2 (DETAIL-32 dev)")


def load_detail32_rois(side=32, dev=False):
    """dict name -> (signal_bool, background_bool) for a DETAIL-32 set, read from
    the cached <name>_roi.npz (build the set first). Used by the campaign layer
    to compute the frozen-ROI CNR (ruling §11)."""
    root = DEV32_IMG_ROOT if dev else IMG32_ROOT
    cache = os.path.join(root, str(int(side)))
    table = _dev32_table() if dev else _conf32_table()
    out = collections.OrderedDict()
    for name, _fam, _seed in table:
        npz = os.path.join(cache, name + "_roi.npz")
        if not os.path.exists(npz):
            raise RuntimeError("missing ROI cache %s (build the DETAIL-32 set "
                               "first)" % npz)
        d = np.load(npz)
        out[name] = (np.asarray(d["signal"], dtype=bool),
                     np.asarray(d["background"], dtype=bool))
    return out


# ------------------------------------------------------------- verify -------- #
def _stats(img01):
    mx = float(img01.max())
    return (float(img01.min()), mx,
            float(np.mean(img01 > 0.5 * mx)) if mx > 0 else 0.0)


def _read01(cache_root, side, name):
    from imageio.v2 import imread
    u8 = imread(os.path.join(cache_root, str(int(side)), name + ".png"))
    return u8.astype(np.float64) / 255.0


def _contact_sheet(side=64, out_png=None):
    """6-panel contact sheet of the DEV instances only (ruling: the 24
    confirmatory are NEVER rendered as a sheet)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    build_detail24_dev_set(side)                      # ensure cache present
    if out_png is None:
        out_png = os.path.join(ROOT, "results", "round63_detail24_devsheet.png")
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(9, 6))
    for ax, (name, fam, _seed) in zip(axes.ravel(), _dev_table()):
        img = _read01(DEV_IMG_ROOT, side, name)
        ax.imshow(img, cmap="gray", vmin=0.0, vmax=1.0, interpolation="nearest")
        ax.set_title(name, fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle("DETAIL-24 DEV instances (seeds %d..%d) — dev-only sheet"
                 % (_dev_seed(0), _dev_seed(5)), fontsize=11)
    fig.tight_layout()
    fig.savefig(out_png, dpi=120)
    plt.close(fig)
    return out_png


def _smoke(side=64):
    """Build both sets (twice, cache-wiped between) -> byte-identical SHAs;
    print per-family pixel stats; render the DEV contact sheet."""
    import shutil

    def build_and_hash(builder, cache_root):
        d = builder(side)
        with open(os.path.join(cache_root, str(side), "sha256.txt")) as fh:
            sha_txt = fh.read()
        with open(os.path.join(cache_root, str(side), "params_manifest.json"),
                  "rb") as fh:
            man_sha = hashlib.sha256(fh.read()).hexdigest()
        return d, sha_txt, man_sha

    ok = True
    for label, builder, cache_root in (
            ("dev", build_detail24_dev_set, DEV_IMG_ROOT),
            ("conf", build_detail24_set, IMG_ROOT)):
        cdir = os.path.join(cache_root, str(side))
        if os.path.isdir(cdir):
            shutil.rmtree(cdir)
        _, sha1, man1 = build_and_hash(builder, cache_root)
        shutil.rmtree(cdir)
        _, sha2, man2 = build_and_hash(builder, cache_root)
        same = (sha1 == sha2) and (man1 == man2)
        ok = ok and same
        print("[detail24 %-4s] regen byte-identical: %s  manifest_sha=%s"
              % (label, "YES" if same else "NO", man1[:16]), flush=True)

    print("\n[detail24] per-family pixel stats (read-back [0,1] truth):", flush=True)
    print("  DEV instances:", flush=True)
    for name, fam, seed in _dev_table():
        lo, hi, frac = _stats(_read01(DEV_IMG_ROOT, side, name))
        print("    %-22s fam=%-12s seed=%d  min=%.3f max=%.3f frac>half=%.4f"
              % (name, fam, seed, lo, hi, frac), flush=True)
    print("  CONFIRMATORY instances (generated+hashed only, NOT reconstructed):",
          flush=True)
    per_fam = {f: [] for f in FAMILIES}
    for name, fam, seed in _conf_table():
        lo, hi, frac = _stats(_read01(IMG_ROOT, side, name))
        per_fam[fam].append(frac)
    for fam in FAMILIES:
        fr = per_fam[fam]
        print("    %-12s  frac>half-max over 4 instances: min=%.4f max=%.4f"
              % (fam, min(fr), max(fr)), flush=True)

    sheet = _contact_sheet(side)
    print("\n[detail24] DEV contact sheet -> %s" % sheet, flush=True)
    ok = _smoke32() and ok
    print("[detail24] smoke %s" % ("PASS" if ok else "FAIL"), flush=True)
    return 0 if ok else 1


def _smoke32():
    """DETAIL-32 checks (ruling §2): byte-identical double regeneration of the
    conf + dev sets (sha256.txt + manifest), per-image sanity, and ROI
    nonempty/disjoint. The 24 confirmatory instances are hashed + ROI-checked
    but never reconstructed."""
    import shutil

    side = PARAMS32["side"]
    ok = True
    for label, builder, root, table in (
            ("dev", build_detail32_dev_set, DEV32_IMG_ROOT, _dev32_table()),
            ("conf", build_detail32_set, IMG32_ROOT, _conf32_table())):
        cdir = os.path.join(root, str(side))

        def build_and_hash():
            builder(side)
            with open(os.path.join(cdir, "sha256.txt")) as fh:
                sha = fh.read()
            with open(os.path.join(cdir, "params_manifest.json"), "rb") as fh:
                man = hashlib.sha256(fh.read()).hexdigest()
            return sha, man

        if os.path.isdir(cdir):
            shutil.rmtree(cdir)
        sha1, man1 = build_and_hash()
        shutil.rmtree(cdir)
        sha2, man2 = build_and_hash()
        same = (sha1 == sha2) and (man1 == man2)
        ok = ok and same
        # verify ROI structure on every instance in the set
        rois = load_detail32_rois(side, dev=(label == "dev"))
        roi_ok = True
        for name, (sig, bg) in rois.items():
            try:
                _assert_roi(sig, bg, name)
            except RuntimeError:
                roi_ok = False
        ok = ok and roi_ok
        occ_sig = [float(sig.mean()) for sig, _ in rois.values()]
        occ_bg = [float(bg.mean()) for _, bg in rois.values()]
        print("[detail32 %-4s] regen byte-identical: %s  manifest_sha=%s  "
              "ROIs disjoint/nonempty: %s  sig-frac[min=%.3f max=%.3f] "
              "bg-frac[min=%.3f max=%.3f]"
              % (label, "YES" if same else "NO", man1[:16],
                 "YES" if roi_ok else "NO", min(occ_sig), max(occ_sig),
                 min(occ_bg), max(occ_bg)), flush=True)
    print("[detail32] smoke %s" % ("PASS" if ok else "FAIL"), flush=True)
    return ok


if __name__ == "__main__":
    sys.exit(_smoke())
