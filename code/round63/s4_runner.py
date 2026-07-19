"""S4(ii): image-level EXACT-vs-RQL split (spec D2 S4; R9 Q5 ruling
"Run it. Do not drop it.").

Frozen grid (spec): rho in {.03,.1,.3,.6,1,2} x nu in {20,100,500,2000},
3 seeds, sides 8 and 16, EXACT-TV and RQL sharing the identical TV
functional and frozen analytic selection rule (both flow through
campaign.run_cell verbatim).

IMAGE-SET NOTE (documented BEFORE any S4 result existed): the conf/dev
builders degenerate at <=16 px (text-family images collapse to zero sum),
and the F1 freeze pinned no S4 image list (expected_cells has no S4
stage). S4 images are therefore the six detail32_dev family instances
(dev-legal), block-averaged to the target side and re-normalized. This is
a reference/no-gate stage (R9: "No new gate is needed").

M = n (square geometry, campaign convention), bern50 patterns, seed-frozen.
Resume-safe CSV; summary = EXACT-TV minus RQL PSNR_rad gap per (rho, nu).
"""
import os
import sys
import csv
import glob
import itertools

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.chdir(os.path.dirname(HERE) if os.path.basename(HERE) == "round63"
         else HERE)
os.chdir(r"D:\GI_another")
import campaign

RHOS = [0.03, 0.1, 0.3, 0.6, 1.0, 2.0]
NUS = [20, 100, 500, 2000]
SEEDS = [0, 1, 2]
SIDES = [8, 16]
ARMS = ["RQL", "EXACT-TV"]

def build_s4_images(side):
    imgs = {}
    for f in sorted(glob.glob(
            r"D:\GI_another\data\r63_images_detail32_dev\32\detail32_dev_*.png")):
        name = os.path.basename(f)[:-4]
        if name.endswith("_roi"):
            continue
        x = np.asarray(Image.open(f).convert("L"), float)
        blk = 32 // side
        x = x.reshape(side, blk, side, blk).mean(axis=(1, 3))
        s = x.sum()
        if s <= 0:
            continue
        imgs["s4_%d_%s" % (side, name.split("_")[-1])] = (x / s).ravel()
    return imgs

_S4_IMGS = {side: build_s4_images(side) for side in SIDES}
_orig_images = campaign._images
def routed(side, spec, dev=False, imageset="conf"):
    if isinstance(spec, list) and spec and str(spec[0]).startswith("s4_"):
        pool = _S4_IMGS[side]
        return {k: pool[k] for k in spec}
    return _orig_images(side, spec, dev=dev, imageset=imageset)
campaign._images = routed

OUT = r"D:\GI_another\results\round63_s4\s4_rows.csv"
os.makedirs(os.path.dirname(OUT), exist_ok=True)
done = set()
if os.path.exists(OUT):
    with open(OUT, newline="") as f:
        for r in csv.DictReader(f):
            done.add((r["image"], r["arm"], r["rho_bar"], r["nu"], r["seed"]))
fout = open(OUT, "a" if done else "w", newline="")
w = csv.writer(fout)
if not done:
    w.writerow(["side", "image", "arm", "rho_bar", "nu", "seed",
                "PSNR_rad", "lam_tv", "runtime_s"])

cells = [(side, img, rho, nu, seed)
         for side in SIDES
         for img in sorted(_S4_IMGS[side])
         for rho, nu, seed in itertools.product(RHOS, NUS, SEEDS)]
total = len(cells)
for i, (side, img, rho, nu, seed) in enumerate(cells, 1):
    if all((img, arm, str(float(rho)), str(float(nu)), str(seed)) in done
           for arm in ARMS):
        continue
    cell = {"side": side, "nu": nu, "rho_bar": rho, "seed": seed,
            "M": side * side, "pattern": "bern50", "arms": ARMS,
            "images": [img], "cell_id": "s4_%s_%g_%g_%d" % (img, rho, nu, seed),
            "audit": False}
    rows = campaign.run_cell(cell)
    for r in rows:
        w.writerow([side, r["image"], r["arm"], rho, nu, seed,
                    r["PSNR_rad"], r.get("lam_tv", ""), r.get("runtime_s", "")])
    fout.flush()
    print("[s4] cell %d/%d %s rho=%g nu=%g s%d done" % (i, total, img, rho, nu, seed),
          flush=True)
fout.close()
print("[s4] SWEEP DONE", flush=True)

# ---- summary: EXACT-TV minus RQL gap ----
import collections
vals = collections.defaultdict(list)
with open(OUT, newline="") as f:
    for r in csv.DictReader(f):
        if r["PSNR_rad"] in ("", "inf"):
            continue
        vals[(r["side"], r["image"], r["rho_bar"], r["nu"], r["seed"])
             ] = vals[(r["side"], r["image"], r["rho_bar"], r["nu"], r["seed"])]
with open(OUT, newline="") as f:
    per = collections.defaultdict(dict)
    for r in csv.DictReader(f):
        if r["PSNR_rad"] in ("", "inf"):
            continue
        per[(r["side"], r["image"], r["rho_bar"], r["nu"], r["seed"])][r["arm"]] = \
            float(r["PSNR_rad"])
gaps = collections.defaultdict(list)
for k, d in per.items():
    if "RQL" in d and "EXACT-TV" in d:
        side, _img, rho, nu, _s = k
        gaps[(side, float(rho), float(nu))].append(d["EXACT-TV"] - d["RQL"])
print("\n[s4] === EXACT-TV minus RQL PSNR_rad gap (median over images/seeds) ===")
print("side rho    " + "  ".join("nu=%-5g" % n for n in NUS))
for side in ("8", "16"):
    for rho in RHOS:
        row = []
        for nu in NUS:
            g = gaps.get((side, rho, float(nu)), [])
            row.append("%+.2f  " % float(np.median(g)) if g else "  --   ")
        print("%4s %-5g " % (side, rho) + "".join(row))
print("[s4] DONE", flush=True)
