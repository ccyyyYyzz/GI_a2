"""M1 unblinding: assemble fetched shard CSVs -> frozen analyzer -> verdicts.

Assembly rules (frozen by R17/R18; this script only maps columns):
  - per image (24), per arm, mean PSNR_rad over the 5 frozen seeds at each nu;
  - safe / ridge nine-dwell curves carry per-image mean achieved rho_bar(nu);
  - q060_terminal / ridge_terminal at nu=2000; any missing/non-finite cell
    stays NaN -> analyzer ITT treats it as nonpositive dQ / censored;
  - cert_cells: the image's 20 R18 status strings (5 seeds x 2 nu x 2 b),
    passed through untouched (branch-aware analyzer handles descriptive mode).
Run from repo root: python code/round63/m1_score.py
"""
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.getcwd(), "code"))
sys.path.insert(0, os.path.join(os.getcwd(), "code", "round63"))
import m1_runner as m1  # noqa: E402

SHARDS = os.path.join("results", "round63_m1", "shards")
OUT_JSON = os.path.join("results", "round63_m1", "M1_VERDICTS.json")
OUT_MD = os.path.join("results", "round63_m1", "M1_VERDICTS.md")


def read_csv(path):
    with open(path, encoding="utf-8") as f:
        hdr = f.readline().strip().split(",")
        rows = []
        for line in f:
            # cert solver-status fields may embed commas inside quotes is NOT
            # the case here (writer uses plain join); but LP_FAIL messages
            # contain commas -> split with maxsplit guard via column count
            parts = line.rstrip("\n").split(",")
            if len(parts) > len(hdr):
                # merge overflow back into the two solver_status columns
                i0 = hdr.index("solver_status_primary")
                extra = len(parts) - len(hdr)
                parts = (parts[:i0]
                         + [",".join(parts[i0:i0 + 1 + extra])]
                         + parts[i0 + 1 + extra:])
                if len(parts) > len(hdr):
                    i1 = hdr.index("solver_status_retry")
                    extra = len(parts) - len(hdr)
                    parts = (parts[:i1]
                             + [",".join(parts[i1:i1 + 1 + extra])]
                             + parts[i1 + 1 + extra:])
            rows.append(dict(zip(hdr, parts)))
    return rows


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


imaging, cert = [], []
for p in sorted(glob.glob(os.path.join(SHARDS, "*.csv"))):
    rows = read_csv(p)
    (cert if os.path.basename(p).startswith("M1_CERT") else imaging).extend(rows)

print("rows: imaging=%d cert=%d" % (len(imaging), len(cert)))
assert len(imaging) == 5400 and len(cert) == 480, "cell accounting mismatch"

images = sorted({r["image"] for r in imaging})
assert len(images) == 24, images
NUS = sorted({fnum(r["nu"]) for r in imaging})
print("images=%d nus=%s" % (len(images), NUS))

# arm -> image -> nu -> [per-seed PSNR_rad], and rho accumulator
# NOTE: the CSV "arm" column holds the estimator name (RQL); the M1 arm
# identity is the shard_id prefix (M1_<ARM>_<idx>).
def m1_arm(r):
    sid = r["shard_id"]
    return sid[len("M1_"):sid.rindex("_")]


acc, rho_acc = {}, {}
for r in imaging:
    key = (m1_arm(r), r["image"], fnum(r["nu"]))
    acc.setdefault(key, []).append(fnum(r["PSNR_rad"]))
    rho_acc.setdefault(key, []).append(fnum(r["rho_bar"]))


def curve(arm, img):
    qs, rhos = [], []
    for nu in NUS:
        v = acc.get((arm, img, nu), [])
        q = float(np.mean(v)) if len(v) == 5 and np.isfinite(v).all() \
            else float("nan")
        rh = rho_acc.get((arm, img, nu), [])
        rhov = float(np.mean(rh)) if rh else float("nan")
        qs.append(q)
        rhos.append(rhov)
    return NUS, rhos, np.array(qs)


def terminal(arm, img):
    v = acc.get((arm, img, 2000.0), [])
    return float(np.mean(v)) if len(v) == 5 and np.isfinite(v).all() \
        else float("nan")


cert_by_img = {}
for r in cert:
    cert_by_img.setdefault(r["image"], []).append(r["status"])

curves = {}
for img in images:
    curves[img] = {
        "family": m1._family_of(img),
        "safe": curve("SCAT32-SAFE", img),
        "ridge": curve("RIDGE-SCAT32", img),
        "q060_terminal": terminal("SCAT32-060", img),
        "ridge_terminal": terminal("RIDGE-SCAT32", img),
        "cert_cells": cert_by_img.get(img, []),
    }
n_cc = sum(len(c["cert_cells"]) for c in curves.values())
assert n_cc == 480, n_cc

res = m1.m1_analyze_r17(curves, n_images=24, boot_B=10000)

# descriptive extras for the paper (per-family dQ medians, cert status dist)
fam_dQ = {}
for img, c in curves.items():
    dq = (c["ridge_terminal"] - c["q060_terminal"]
          if np.isfinite(c["ridge_terminal"]) and np.isfinite(c["q060_terminal"])
          else 0.0)
    fam_dQ.setdefault(c["family"], []).append(dq)
res["per_family_median_dQ"] = {f: float(np.median(v))
                               for f, v in sorted(fam_dQ.items())}
from collections import Counter  # noqa: E402
res["cert_status_distribution"] = dict(Counter(r["status"] for r in cert))
res["cert_by_anchor"] = {}
for r in cert:
    k = "nu%g_b%g" % (fnum(r["nu"]), fnum(r["b"]))
    res["cert_by_anchor"].setdefault(k, Counter())[r["status"]] += 1
res["cert_by_anchor"] = {k: dict(v) for k, v in res["cert_by_anchor"].items()}

with open(OUT_JSON, "w") as f:
    json.dump(res, f, indent=1, default=str)

lines = ["# M1 verdicts (unblinded %s, tag m1-freeze @ 6f00932)" % "2026-07-19",
         "",
         "RIDGE_OPERATING_PASS = %s" % res["RIDGE_OPERATING_PASS"],
         "  median dQ = %.3f dB (bar >= 1.0)" % res["median_dQ_dB"],
         "  LB2.5 = %.3f (bar > 0)" % res["dQ_LB2.5"],
         "  n_pos = %d/24 (bar >= 18)" % res["n_dQ_pos"],
         "",
         "RIDGE_SPEED_PASS = %s" % res["RIDGE_SPEED_PASS"],
         "  median S = %.3f (bar >= 3)" % res["median_S"],
         "  S LB2.5 = %.3f (bar > 1)" % res["S_LB2.5"],
         "  n_S>1 = %d/24 (bar >= 18)" % res["n_S_gt1"],
         "",
         "cert branch = %s" % res["cert_branch"],
         "  descriptive certified fraction = %s" %
         res.get("full_stack_cert_descriptive_fraction"),
         "  status distribution: %s" % res["cert_status_distribution"],
         "  by anchor: %s" % res["cert_by_anchor"],
         "",
         "per-family median dQ: %s" % res["per_family_median_dQ"]]
with open(OUT_MD, "w") as f:
    f.write("\n".join(lines) + "\n")
print("\n".join(lines))
