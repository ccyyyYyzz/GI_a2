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
    """Family string for a DETAIL-24 image name (conf or dev)."""
    parts = name.split("_")
    fam = parts[2] if len(parts) >= 3 and parts[1] == "dev" else parts[1]
    if fam not in FAMILIES:
        raise ValueError("name %r does not encode a DETAIL-24 family" % name)
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
    print("[detail24] smoke %s" % ("PASS" if ok else "FAIL"), flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_smoke())
