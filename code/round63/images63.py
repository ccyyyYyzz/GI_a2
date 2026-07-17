"""ROUND63 image set — spec §4 (>=24 自然图 + 6 结构靶 @ 64^2 / 128^2).

build_image_set(side) returns an ordered dict  name -> flattened float64
sum-normalized truth (sum_j x_j = 1).  Order: 24 natural (stl_00..stl_23) then
6 synthetic structural targets.

Canonical protocol (mirrors gi_core/images.py + data_fetch.py): the 8-bit PNG is
the ground truth of record.  Each image is generated once, written to
data/r63_images/<side>/<name>.png, and the returned truth is read back from that
PNG (uint8/255 then sum-normalized), so truth is exactly reproducible from the
cache and matches its sha256.  A sha256 manifest is written alongside.

Natural images: STL-10 test split indices 0..23, parsed directly from the raw
binary already on disk (data/stl10_raw/stl10_binary/test_X.bin; uint8,
column-major 3x96x96 per image -> transpose to HWC).  NO re-download.  Grayscale
via skimage rgb2gray (Rec.709 luma), bicubic (order=3) resize to side, clip [0,1].

Synthetic targets: generated procedurally with numpy/PIL, no downloads.  Any
randomness lives under rng_for(0, 63, 2, synth_id) (SEED0 system, spec §5;
sub-stream 2 = image layer, disjoint from patterns' sub-stream 1).
"""
import collections
import hashlib
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from gi_core.utils import rng_for

ROOT = os.path.dirname(os.path.dirname(HERE))
DATA = os.path.join(ROOT, "data")
STL_TEST = os.path.join(DATA, "stl10_raw", "stl10_binary", "test_X.bin")
STL_TRAIN = os.path.join(DATA, "stl10_raw", "stl10_binary", "train_X.bin")
IMG_ROOT = os.path.join(DATA, "r63_images")
DEV_IMG_ROOT = os.path.join(DATA, "r63_images_dev")

R63 = 63
IMG_SUBSTREAM = 2          # image layer sub-stream (patterns use sub-stream 1)
N_NATURAL = 24
N_DEV_NATURAL = 16         # STL-10 TRAIN split indices 0..15 (S1 development set)
SYNTH = ["usaf_bars", "fine_lines", "text", "sparse_dots", "low_contrast",
         "gray_staircase"]
SYNTH_IDS = {name: i for i, name in enumerate(SYNTH)}


# --------------------------------------------------------------------------- #
# natural images (STL-10 test split, parsed from raw binary)
# --------------------------------------------------------------------------- #
def _stl10_one(idx, side):
    """Grayscale side x side float64 [0,1] truth of STL-10 test image `idx`."""
    from skimage.color import rgb2gray
    from skimage.transform import resize

    need = (idx + 1) * 3 * 96 * 96
    raw = np.fromfile(STL_TEST, dtype=np.uint8, count=need)
    if raw.size < need:
        raise RuntimeError(
            "STL-10 test split too short/missing at %s (need >=%d bytes, got %d); "
            "expected raw binary already on disk (do NOT re-download)."
            % (STL_TEST, need, raw.size))
    # column-major 3x96x96 per image -> (n, 96, 96, 3) HWC (canonical STL-10 read)
    img = raw.reshape(-1, 3, 96, 96).transpose(0, 3, 2, 1)[idx]
    g = rgb2gray(img)                                       # float64 [0,1]
    g = resize(g, (side, side), order=3, anti_aliasing=True, mode="reflect")
    return np.clip(g, 0.0, 1.0)


def _stl10_train_one(idx, side):
    """Grayscale side x side float64 [0,1] truth of STL-10 TRAIN image `idx`
    (S1 development set). Identical parse/read to _stl10_one but from the TRAIN
    split, which is DISJOINT from the test split used by build_image_set — so the
    S1 pilot images never intersect the S2 confirmatory set (digest §1.5A)."""
    from skimage.color import rgb2gray
    from skimage.transform import resize

    need = (idx + 1) * 3 * 96 * 96
    raw = np.fromfile(STL_TRAIN, dtype=np.uint8, count=need)
    if raw.size < need:
        raise RuntimeError(
            "STL-10 train split too short/missing at %s (need >=%d bytes, got %d);"
            " expected raw binary already on disk (do NOT re-download)."
            % (STL_TRAIN, need, raw.size))
    img = raw.reshape(-1, 3, 96, 96).transpose(0, 3, 2, 1)[idx]
    g = rgb2gray(img)
    g = resize(g, (side, side), order=3, anti_aliasing=True, mode="reflect")
    return np.clip(g, 0.0, 1.0)


# --------------------------------------------------------------------------- #
# synthetic structural targets (procedural, side-agnostic)
# --------------------------------------------------------------------------- #
def _syn_usaf_bars(side):
    """USAF-1951-style three-bar groups at 3 scales; horizontal groups in the
    top half, vertical groups in the bottom half."""
    img = np.zeros((side, side), np.float64)
    scales = [max(3, side // 12), max(2, side // 20), max(1, side // 32)]
    third = max(1, side // 3)
    for gi, t in enumerate(scales):                        # horizontal 3-bar groups
        L = min(third - 1, 6 * t)
        c0 = gi * third + max(0, (third - L) // 2)
        r0 = max(1, (side // 2 - 5 * t) // 2)
        for b in range(3):
            r = r0 + b * 2 * t
            if r + t <= side // 2 and c0 + L <= side:
                img[r:r + t, c0:c0 + L] = 1.0
    for gi, t in enumerate(scales):                        # vertical 3-bar groups
        L = min(side - side // 2 - 1, 6 * t)
        c0 = gi * third + max(0, (third - 6 * t) // 2)
        r0 = side // 2 + max(1, (side // 2 - L) // 2)
        for b in range(3):
            c = c0 + b * 2 * t
            if c + t <= side and r0 + L <= side:
                img[r0:r0 + L, c:c + t] = 1.0
    return img


def _syn_fine_lines(side):
    """1-px and 2-px gratings at several orientations, one per quadrant."""
    img = np.zeros((side, side), np.float64)
    h = max(1, side // 2)
    img[0:h:2, 0:h] = 1.0                                  # TL: horizontal 1-px
    for c in range(h, side, 4):                            # TR: vertical 2-px
        img[0:h, c:min(c + 2, side)] = 1.0
    img[h:side, 0:h:2] = 1.0                               # BL: vertical 1-px
    rr, cc = np.mgrid[h:side, h:side]                      # BR: 1-px diagonals
    img[h:side, h:side] = ((rr + cc) % 4 == 0).astype(np.float64)
    return img


def _syn_text(side):
    """Binarized render of "OE 2026 GI" with PIL's default bitmap font, scaled to
    fit the frame width (aspect preserved) and vertically centered."""
    from PIL import Image, ImageDraw, ImageFont
    from skimage.transform import resize

    s = "OE 2026 GI"
    font = ImageFont.load_default()
    probe = ImageDraw.Draw(Image.new("L", (16, 16), 0))
    l, t, r, b = probe.textbbox((0, 0), s, font=font)
    tw, th = max(1, r - l), max(1, b - t)
    txt = Image.new("L", (tw, th), 0)
    ImageDraw.Draw(txt).text((-l, -t), s, fill=255, font=font)
    arr = np.asarray(txt, dtype=np.float64) / 255.0
    scale = min((side - 2) / tw, (side - 2) / th)          # fit within frame
    nw, nh = max(1, int(round(tw * scale))), max(1, int(round(th * scale)))
    small = resize(arr, (nh, nw), order=1, anti_aliasing=True, mode="constant")
    img = np.zeros((side, side), np.float64)
    r0, c0 = (side - nh) // 2, (side - nw) // 2
    img[r0:r0 + nh, c0:c0 + nw] = small
    return (img > 0.4).astype(np.float64)


def _syn_sparse_dots(side, rng):
    """12 isolated 1-2 px dots on a dark background."""
    img = np.zeros((side, side), np.float64)
    margin = max(2, side // 16)
    hi = max(margin + 1, side - margin)
    for _ in range(12):
        r = int(rng.integers(margin, hi))
        c = int(rng.integers(margin, hi))
        rad = int(rng.integers(1, 3))                      # 1 or 2 px
        img[r:min(r + rad, side), c:min(c + rad, side)] = 1.0
    return img


def _syn_low_contrast(side, rng):
    """Smooth 0.5 background with a few small +0.05-contrast structures."""
    img = np.full((side, side), 0.5, np.float64)
    lo, hi = max(2, side // 12), max(3, side // 6)
    for _ in range(6):
        h = int(rng.integers(lo, hi + 1))
        w = int(rng.integers(lo, hi + 1))
        r = int(rng.integers(0, max(1, side - h)))
        c = int(rng.integers(0, max(1, side - w)))
        img[r:r + h, c:c + w] += 0.05
    return np.clip(img, 0.0, 1.0)


def _syn_gray_staircase(side):
    """8-level horizontal staircase wedge, levels i/7 for i = 0..7."""
    img = np.zeros((side, side), np.float64)
    for i in range(8):
        r0, r1 = i * side // 8, (i + 1) * side // 8
        img[r0:r1, :] = i / 7.0
    return img


def _generate01(name, side):
    """Generate one image as float64 [0,1] side x side (pre-quantization)."""
    if name.startswith("stl_"):
        return _stl10_one(int(name.split("_")[1]), side)
    if name == "usaf_bars":
        return _syn_usaf_bars(side)
    if name == "fine_lines":
        return _syn_fine_lines(side)
    if name == "text":
        return _syn_text(side)
    rng = rng_for(0, R63, IMG_SUBSTREAM, SYNTH_IDS.get(name, 0))
    if name == "sparse_dots":
        return _syn_sparse_dots(side, rng)
    if name == "low_contrast":
        return _syn_low_contrast(side, rng)
    if name == "gray_staircase":
        return _syn_gray_staircase(side)
    raise ValueError("unknown image name %r" % name)


def _to_uint8(img01):
    return np.round(np.clip(img01, 0.0, 1.0) * 255.0).astype(np.uint8)


def build_image_set(side):
    """Ordered dict name -> flattened float64 sum-normalized truth (spec §4).

    PNGs + sha256.txt are cached under data/r63_images/<side>/; existing PNGs are
    treated as canonical and reused."""
    from imageio.v2 import imread, imwrite

    side = int(side)
    cache = os.path.join(IMG_ROOT, str(side))
    os.makedirs(cache, exist_ok=True)
    names = ["stl_%02d" % i for i in range(N_NATURAL)] + list(SYNTH)

    out = collections.OrderedDict()
    sha_lines = []
    for name in names:
        png = os.path.join(cache, name + ".png")
        if not os.path.exists(png):
            imwrite(png, _to_uint8(_generate01(name, side)))
        u8 = imread(png)
        if u8.shape != (side, side):
            raise RuntimeError("cached PNG %s has shape %s, expected (%d,%d)"
                               % (png, u8.shape, side, side))
        x = (u8.astype(np.float64) / 255.0).ravel()
        s = x.sum()
        if s <= 0:
            raise RuntimeError("image %s has non-positive sum %g" % (name, s))
        out[name] = x / s
        with open(png, "rb") as f:
            sha_lines.append("%s  %s.png"
                             % (hashlib.sha256(f.read()).hexdigest(), name))
    with open(os.path.join(cache, "sha256.txt"), "w") as f:
        f.write("\n".join(sha_lines) + "\n")
    return out


# --------------------------------------------------------------------------- #
# S1 development image set (STL-10 TRAIN split + reused synthetic targets)
# --------------------------------------------------------------------------- #
def _generate01_dev(name, side):
    """Generate one DEVELOPMENT image as float64 [0,1] side x side.

    dev_stl_NN -> STL-10 TRAIN image NN (disjoint from the test split).
    dev_<synth> -> the SAME procedural synthetic target as build_image_set
    (the 6 structural targets are not drawn from any split, so reusing them
    under dev_ names is exact and carries no test-set leakage)."""
    if name.startswith("dev_stl_"):
        return _stl10_train_one(int(name.split("_")[2]), side)
    if not name.startswith("dev_"):
        raise ValueError("dev image name %r must start with 'dev_'" % name)
    base = name[len("dev_"):]
    if base == "usaf_bars":
        return _syn_usaf_bars(side)
    if base == "fine_lines":
        return _syn_fine_lines(side)
    if base == "text":
        return _syn_text(side)
    rng = rng_for(0, R63, IMG_SUBSTREAM, SYNTH_IDS.get(base, 0))
    if base == "sparse_dots":
        return _syn_sparse_dots(side, rng)
    if base == "low_contrast":
        return _syn_low_contrast(side, rng)
    if base == "gray_staircase":
        return _syn_gray_staircase(side)
    raise ValueError("unknown dev image name %r" % name)


def build_dev_image_set(side):
    """Ordered dict name -> flattened float64 sum-normalized truth for the S1
    DEVELOPMENT set (digest §1.5A): 16 natural images from the STL-10 TRAIN split
    (indices 0..15) named dev_stl_00..dev_stl_15, then the 6 synthetic structural
    targets reused under dev_ names. These are S1-pilot-only and DISJOINT from the
    confirmatory build_image_set (which uses the STL-10 TEST split) — protocol
    tuning on this set never touches the frozen S2 test images.

    PNGs + sha256.txt are cached under data/r63_images_dev/<side>/; existing PNGs
    are treated as canonical and reused. build_image_set is not touched."""
    from imageio.v2 import imread, imwrite

    side = int(side)
    cache = os.path.join(DEV_IMG_ROOT, str(side))
    os.makedirs(cache, exist_ok=True)
    names = (["dev_stl_%02d" % i for i in range(N_DEV_NATURAL)]
             + ["dev_" + s for s in SYNTH])

    out = collections.OrderedDict()
    sha_lines = []
    for name in names:
        png = os.path.join(cache, name + ".png")
        if not os.path.exists(png):
            imwrite(png, _to_uint8(_generate01_dev(name, side)))
        u8 = imread(png)
        if u8.shape != (side, side):
            raise RuntimeError("cached dev PNG %s has shape %s, expected (%d,%d)"
                               % (png, u8.shape, side, side))
        x = (u8.astype(np.float64) / 255.0).ravel()
        s = x.sum()
        if s <= 0:
            raise RuntimeError("dev image %s has non-positive sum %g" % (name, s))
        out[name] = x / s
        with open(png, "rb") as f:
            sha_lines.append("%s  %s.png"
                             % (hashlib.sha256(f.read()).hexdigest(), name))
    with open(os.path.join(cache, "sha256.txt"), "w") as f:
        f.write("\n".join(sha_lines) + "\n")
    return out


def _smoke():
    """Build side=64, print names + sha count, verify sums == 1 within 1e-12."""
    side = 64
    d = build_image_set(side)
    print("[images63 smoke] side=%d built %d images:" % (side, len(d)), flush=True)
    print("  " + ", ".join(d.keys()), flush=True)
    with open(os.path.join(IMG_ROOT, str(side), "sha256.txt")) as f:
        n_sha = sum(1 for ln in f if ln.strip())
    max_err = 0.0
    for name, x in d.items():
        err = abs(float(x.sum()) - 1.0)
        max_err = max(max_err, err)
        if err >= 1e-12:
            print("  FAIL %s sum=%.16f" % (name, x.sum()), flush=True)
    ok = max_err < 1e-12
    print("[images63 smoke] sha lines=%d  max|sum-1|=%.2e  %s"
          % (n_sha, max_err, "PASS" if ok else "FAIL"), flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_smoke())
