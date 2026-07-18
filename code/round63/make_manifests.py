"""ROUND63 shard-manifest + expected-cell generator (spec D2.3 §0/§3, GPT
round-7 F1 ruling).

F1-FREEZE ARTIFACT BUILDER.  Pure, deterministic generation + validation — this
script NEVER runs a reconstruction (no campaign.run_cell in the default path).
Re-running it is a pure function of the spec constants + the (frozen) input PNGs
and the S1 Pass-B runtime calibration, so the outputs are byte-identical across
runs.  Everything is sorted; there is no wall-clock and no RNG.

It emits, under results/round63/manifests/ + results/round63/:
  * manifests/<shard_id>.json  — one shard manifest per shard, following
    shard_runner.py's documented schema EXACTLY (shard_id, stage, frozen_inputs
    with sha256, cells, output_csv, meta_json) plus the requested extra fields
    (priority, account_hint, est_hours, and blocked_by on blocked shards).
  * manifests/MANIFEST_SUMMARY.md — the human table (shard, stage, cells, est
    hours, account, blocked) with per-stage / per-account totals and notes.
  * manifests/MANIFEST_INDEX.json — machine index of default (runnable) vs
    blocked shards + totals, for the runner orchestration layer.
  * expected_cells.csv — the COMPLETE expected-cell table: one row per
    (shard, stage, side, pattern, rho_bar, nu, M, seed, image, arm, imageset,
    mismatch, blocked_by) with images expanded to real names (DETAIL rows expand
    to the 24 detail_<family>_<i> names).  This is the F1 "complete expected-cell
    table" and the input to the merge-time expected/duplicate/missing gate.

Grids encoded (spec D2.3 §3, GPT round-7 ruling §1):
  * S2A_DETAIL   (pri 1) — 64^2 bern50 M=4096; the DETAIL-24 primary confirmatory
                           positive-speed cohort (imageset detail24, all 24).
  * S2A_NATURAL  (pri 2) — same physics grid on the STL-test NATURAL-24 secondary
                           regime-boundary cohort (imageset conf, all 30; no gate).
  * S2B          (pri 3) — undersampling robustness M in {1024,2048} (M=4096 rows
                           dropped: they duplicate the S2A_NATURAL anchor).
  * S2C / S2C128 (pri 3) — hadpair/gam4 pattern anchors (M=2048) + the 128^2 block
                           (blocked on the matrix-free operator).
  * S3 / S3JIT / S3CONT  — OAT mismatch axes + interactions (M=2048, nu=500 anchor,
                           M/n=0.5 per spec); jitter is priority-4, continuous is
                           blocked on eta* inheritance.
S4 is EXCLUDED from manifests (scalar Fisher map = standalone script; the 8^2/16^2
image split is blocked on the EXACT+TV arm) and is only noted in the summary.

Every emitted cell dict satisfies the campaign.run_cell contract (side, pattern,
rho_bar, nu, M, seed, arms, images, imageset, tau, sigma_b, select_rule) by
construction; validate_cell() re-checks keys/types by inspection on every cell.
No cell carries C0/c_force — the frozen C0_FROZEN.json is authoritative at run
time; select_rule stays at the default "discrepancy" (the F1 analytic rule).
Run `--selftest` to additionally push ONE tiny cell through campaign.run_cell.

Usage:
  D:/Anacondar/anaconda3/python.exe code/round63/make_manifests.py [--target-hours 6.0]
  D:/Anacondar/anaconda3/python.exe code/round63/make_manifests.py --selftest
"""
import argparse
import csv
import glob
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))          # code/ (for gi_core)
ROOT = os.path.dirname(os.path.dirname(HERE))       # repo root (D:\GI_another)

import images63                                     # light: numpy + gi_core.utils
import detail24                                     # light: numpy (scipy lazy)

# --------------------------------------------------------------------------- #
# canonical constants
# --------------------------------------------------------------------------- #
TAU = 50e-9                                          # spec §2: tau = 50 ns fixed
SIGMA_B = 0.0                                        # spec §2: sigma_b = 0 main
SELECT_RULE = "discrepancy"                          # spec §4: F1 analytic rule

# ------- confirmatory NATURAL-24 image build order (mirrors images63) -------- #
CONF_NAMES = (["stl_%02d" % i for i in range(images63.N_NATURAL)]
              + list(images63.SYNTH))               # 24 STL test + 6 structural
CONF_INDEX = {name: i for i, name in enumerate(CONF_NAMES)}
STL12 = ["stl_%02d" % i for i in range(12)]
STL8_PLUS_4 = (["stl_%02d" % i for i in range(8)]    # first 8 STL ...
               + list(images63.SYNTH[:4]))           # + first 4 structural targets

# S1 development images (used only by the non-confirmatory selftest)
DEV_NAMES = (["dev_stl_%02d" % i for i in range(images63.N_DEV_NATURAL)]
             + ["dev_" + s for s in images63.SYNTH])
DEV_INDEX = {name: i for i, name in enumerate(DEV_NAMES)}

# ------- DETAIL-24 primary confirmatory build order (detail24 ruling §1) ----- #
DETAIL_NAMES = [name for name, _fam, _seed in detail24._conf_table()]   # 24
DETAIL_INDEX = {name: i for i, name in enumerate(DETAIL_NAMES)}
DETAIL_DEV_NAMES = [name for name, _f, _s in detail24._dev_table()]     # 6
DETAIL_DEV_INDEX = {name: i for i, name in enumerate(DETAIL_DEV_NAMES)}

# imageset -> ("all" expansion list, canonical index map, cache-root rel path)
IMAGESET_ALL = {
    "conf": CONF_NAMES, "dev": DEV_NAMES,
    "detail24": DETAIL_NAMES, "detail24_dev": DETAIL_DEV_NAMES,
}
IMAGESET_INDEX = {
    "conf": CONF_INDEX, "dev": DEV_INDEX,
    "detail24": DETAIL_INDEX, "detail24_dev": DETAIL_DEV_INDEX,
}
IMAGESET_ROOT = {
    "conf": "data/r63_images", "dev": "data/r63_images_dev",
    "detail24": "data/r63_images_detail24",
    "detail24_dev": "data/r63_images_detail24_dev",
}
ALLOWED_IMAGESET = set(IMAGESET_ALL)

# arm names known to solvers.run_arm (for validation) + a canonical order
KNOWN_ARMS = ["RQL", "POISSON-LIN", "SAT-POISSON", "PRECORRECT", "QMLE",
              "QMLE-FULLGAUSS", "GI", "DGI", "EXACT"]
ARM_RANK = {a: i for i, a in enumerate(KNOWN_ARMS)}
LINEAR_ARMS = {"GI", "DGI"}

KNOWN_PATTERNS = {"bern50", "hadpair", "gam4"}
ALLOWED_DET = {"p_ap", "ap_tau", "paralyzable", "tau_jitter_cv", "start_mode",
               "guard", "dark_frac"}
ALLOWED_EST = {"tau_err", "dark_known", "assume_paralyzable"}

# cost-model bases — analytic FALLBACK only.  The default cost model is the S1
# Pass-B runtime_calibration.json (REAL 64^2/M=4096 per-(nu,arm) wall seconds
# INCLUDING the analytic lam_TV selection + descriptive audit), M-scaled by
# (M/4096); the fallback below is exercised only if a (nu,arm) key is missing.
BASE_RQL = 300.0
BASE_ITER = 150.0
BASE_LINEAR = 5.0
NU_REF = 500.0
NU_EXP = 0.15
M_REF = 4096.0                                       # calibration reference M
# hadpair charges 2 physical exposures per signed row (spec §3 "按 2 次曝光计费"),
# so its count-likelihood fits run on 2M physical rows -> ~2x per-fit cost.
PAT_COST_FACTOR = {"bern50": 1.0, "gam4": 1.0, "hadpair": 2.0}


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def g(v):
    """Compact deterministic numeric string (1.0 -> '1', 0.05 -> '0.05')."""
    return "%g" % v


def sort_arms(arms):
    return sorted(dict.fromkeys(arms), key=lambda a: ARM_RANK[a])


def img_index(name):
    """Canonical build-order index for an image name, across every imageset
    (names are disjoint between imagesets: stl_* / <synth> / dev_* / detail_*)."""
    for idx in (CONF_INDEX, DETAIL_INDEX, DETAIL_DEV_INDEX, DEV_INDEX):
        if name in idx:
            return idx[name]
    raise KeyError(name)


def all_names(imageset):
    return list(IMAGESET_ALL[imageset])


def expand_images(cell):
    """'all' -> the imageset's canonical names; an explicit list -> itself."""
    spec = cell["images"]
    if spec == "all":
        return all_names(cell.get("imageset", "conf"))
    return list(spec)


_SHA_CACHE = {}


def sha256_png(imageset, side, name):
    """sha256 of a cached image PNG, or None if the file does not exist yet
    (e.g. the 128^2 anchor block, blocked on the matrix-free operator)."""
    rel = "%s/%d/%s.png" % (IMAGESET_ROOT[imageset], side, name)
    if rel in _SHA_CACHE:
        return _SHA_CACHE[rel]
    ap = os.path.join(ROOT, rel.replace("/", os.sep))
    val = None
    if os.path.exists(ap):
        h = hashlib.sha256()
        with open(ap, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        val = h.hexdigest()
    _SHA_CACHE[rel] = val
    return val


# --------------------------------------------------------------------------- #
# cell construction
# --------------------------------------------------------------------------- #
def make_cell(stage, side, pattern, rho, nu, M, seed, arms, images, imageset,
              det=None, est=None, blocked_by=None):
    """Build one run_cell-ready cell dict (cell_id assigned later).

    NB: no C0/c_force key — the frozen C0_FROZEN.json is authoritative at run
    time (spec §4: C0 is the single frozen calibration constant, Pass A chose
    C0*=inf).  select_rule stays at the default 'discrepancy' (F1 analytic rule)."""
    cell = {
        "stage": stage, "side": int(side), "pattern": pattern,
        "rho_bar": float(rho), "nu": float(nu), "M": int(M), "seed": int(seed),
        "arms": sort_arms(arms), "images": images, "imageset": imageset,
        "tau": TAU, "sigma_b": SIGMA_B, "select_rule": SELECT_RULE,
    }
    if det:
        cell["det"] = dict(det)
    if est:
        cell["est"] = dict(est)
    if blocked_by:
        cell["blocked_by"] = blocked_by
    return cell


# S2-A physics grid (shared by DETAIL and NATURAL; spec D2.3 §3)
S2A_RHO = (0.05, 0.3, 0.6, 1.0, 2.0)
S2A_NU = (5, 10, 20, 50, 100, 200, 500, 1000, 2000)   # nine dwell steps
S2A_SEEDS = range(5)
S2A_M = 4096                                            # M/n = 1 (round-6 ruling)


def s2a_arms(rho, nu):
    """S2-A per-cell arm set (spec §3 / ruling §1): RQL at every (rho,nu); the
    three misspecified physical arms at rho in {0.05,0.6,1} x all nu; GI (display
    reference) at rho in {0.05,0.6} x nu in {100,500,2000}."""
    rk = round(rho, 3)
    arms = ["RQL"]
    if rk in (0.05, 0.6, 1.0):
        arms += ["POISSON-LIN", "SAT-POISSON", "PRECORRECT"]
    if rk in (0.05, 0.6) and int(nu) in (100, 500, 2000):
        arms.append("GI")
    return arms


def build_s2a(stage, imageset, images):
    """One S2-A stage over the shared physics grid; DETAIL and NATURAL differ
    only by (stage, imageset, images)."""
    cells = []
    for rho in S2A_RHO:
        for nu in S2A_NU:
            for seed in S2A_SEEDS:
                cells.append(make_cell(
                    stage, 64, "bern50", rho, nu, S2A_M, seed,
                    s2a_arms(rho, nu), images, imageset))
    return cells


def build_s2a_detail():
    """PRIMARY confirmatory positive-speed cohort — the 24 DETAIL-24 targets."""
    return build_s2a("S2A_DETAIL", "detail24", "all")


def build_s2a_natural():
    """SECONDARY regime-boundary cohort (no positive gate) — STL-test NATURAL-24
    + the 6 structural display targets (30 names)."""
    return build_s2a("S2A_NATURAL", "conf", "all")


def build_s2b():
    """Undersampling robustness (spec §3).  M in {1024,2048,4096} x rho{.05,.6,1}
    x nu{100,500,2000}, 3 seeds, first-12 STL naturals.  The M=4096 rows are a
    strict subset of the S2A_NATURAL anchor (same physics, subset of images/arms/
    seeds), so they are DROPPED here to avoid duplicate expected-cell rows — the
    analogue of the historic 'drop the S2A-anchor M' rule (see LOUD DECISION)."""
    cells = []
    for M in (1024, 2048):                            # 4096 == S2A_NATURAL anchor
        for rho in (0.05, 0.6, 1.0):
            for nu in (100, 500, 2000):
                for seed in range(3):
                    cells.append(make_cell(
                        "S2B", 64, "bern50", rho, nu, M, seed,
                        ["RQL", "POISSON-LIN", "SAT-POISSON", "PRECORRECT"],
                        list(STL12), "conf"))
    return cells


def build_s2c():
    cells = []
    for pattern in ("hadpair", "gam4"):
        for rho in (0.05, 0.6, 1.0):
            for nu in (100, 500, 2000):
                for seed in range(3):
                    cells.append(make_cell(
                        "S2C", 64, pattern, rho, nu, 2048, seed,
                        ["RQL", "POISSON-LIN", "PRECORRECT", "GI"],
                        list(STL8_PLUS_4), "conf"))
    return cells


def build_s2c_128_blocked():
    # 128^2 anchor: Bernoulli, M/n=0.5 -> M=8192, same anchor grid, RQL/SAT/
    # PRECORRECT (PnP-BM3D deferred with the matrix-free operator — see NOTES).
    cells = []
    for rho in (0.05, 0.6, 1.0):
        for nu in (100, 500, 2000):
            for seed in range(3):
                cells.append(make_cell(
                    "S2C128", 128, "bern50", rho, nu, 8192, seed,
                    ["RQL", "SAT-POISSON", "PRECORRECT"],
                    list(STL8_PLUS_4), "conf", blocked_by="matrixfree_128"))
    return cells


def build_s3():
    """S3 OAT mismatch (spec §3).  OAT axes are perturbed from the two primary
    operating points rho in {0.05 (safe ref), 0.6 (high-flux)}; the rho x axis
    interactions add rho in {0.3, 1}.  Their union is exactly {0.05,0.3,0.6,1}
    (spec header 'rho in {0.3,0.6,1}+ref 0.05'), with no duplicate cell — see
    the LOUD DECISION note in main().  Fixed 64^2/bern50/M=2048/nu=500 anchor
    (M/n=0.5, the historic mismatch anchor per spec §3).  Jitter is priority-4
    (first to cut, §1); continuous acquisition is blocked (eta* inheritance)."""
    RHO_OAT = (0.05, 0.6)
    RHO_INT = (0.3, 1.0)
    SEEDS = range(3)
    cells, jit, cont = [], [], []

    def s3(stage, rho, seed, det=None, est=None, blocked_by=None):
        return make_cell(stage, 64, "bern50", rho, 500, 2048, seed,
                         ["RQL", "PRECORRECT"], list(STL12), "conf",
                         det=det, est=est, blocked_by=blocked_by)

    for rho in RHO_OAT:
        for seed in SEEDS:
            # estimator tau error
            for e in (-0.20, -0.10, 0.10, 0.20):
                cells.append(s3("S3", rho, seed, est={"tau_err": e}))
            # dark count, known
            for d in (0.05, 0.10, 0.25, 0.50):
                cells.append(s3("S3", rho, seed, det={"dark_frac": d},
                                est={"dark_known": True}))
            # dark count 0.1, unknown (estimator ignores dark)
            cells.append(s3("S3", rho, seed, det={"dark_frac": 0.10},
                            est={"dark_known": False}))
            # afterpulsing
            for p in (0.01, 0.02, 0.05, 0.10):
                cells.append(s3("S3", rho, seed, det={"p_ap": p}))
            # frame-start: delayed (active is the shared S2-A baseline)
            cells.append(s3("S3", rho, seed, det={"start_mode": "delayed"}))
            # inter-frame guard (continuous-only field; see LOUD DECISION note)
            for gk in (1, 5):
                cells.append(s3("S3", rho, seed, det={"guard": gk * TAU}))
            # dead-time jitter -> priority 4 (own stage S3JIT)
            for j in (0.05, 0.10):
                jit.append(s3("S3JIT", rho, seed, det={"tau_jitter_cv": j}))
            # frame-start: continuous -> BLOCKED (eta* inheritance)
            cont.append(s3("S3CONT", rho, seed, det={"start_mode": "continuous"},
                           blocked_by="continuous_eta_inheritance"))

    for rho in RHO_INT:
        for seed in SEEDS:
            for e in (-0.10, 0.10):                 # rho x tau-err interaction
                cells.append(s3("S3", rho, seed, est={"tau_err": e}))
            for p in (0.02, 0.10):                   # rho x afterpulse interaction
                cells.append(s3("S3", rho, seed, det={"p_ap": p}))
            for p in (0.02, 0.10):                   # continuous x afterpulse (BLOCKED)
                cont.append(s3("S3CONT", rho, seed,
                               det={"start_mode": "continuous", "p_ap": p},
                               blocked_by="continuous_eta_inheritance"))
    return cells, jit, cont


# --------------------------------------------------------------------------- #
# validation (code-inspection against the run_cell contract)
# --------------------------------------------------------------------------- #
def validate_cell(cell):
    """Inspection-only check that a cell dict satisfies campaign.run_cell's
    contract (keys present + typed correctly; det/est/arms/images/imageset legal).
    Raises AssertionError on any violation."""
    cid = cell.get("cell_id", "<no id>")
    for k in ("side", "pattern", "rho_bar", "nu", "M", "seed", "arms", "images",
              "imageset"):
        assert k in cell, "%s: missing key %r" % (cid, k)
    assert isinstance(cell["side"], int) and cell["side"] > 0, cid
    assert cell["pattern"] in KNOWN_PATTERNS, "%s: pattern %r" % (cid, cell["pattern"])
    assert isinstance(cell["rho_bar"], float) and cell["rho_bar"] >= 0.0, cid
    assert isinstance(cell["nu"], float) and cell["nu"] > 0.0, cid
    assert isinstance(cell["M"], int) and cell["M"] > 0, cid
    assert isinstance(cell["seed"], int) and cell["seed"] >= 0, cid
    assert isinstance(cell["arms"], list) and cell["arms"], "%s: empty arms" % cid
    for a in cell["arms"]:
        assert a in KNOWN_ARMS, "%s: unknown arm %r" % (cid, a)
    assert cell["arms"] == sort_arms(cell["arms"]), "%s: arms unsorted" % cid
    imageset = cell["imageset"]
    assert imageset in ALLOWED_IMAGESET, "%s: imageset %r" % (cid, imageset)
    valid_names = set(IMAGESET_ALL[imageset])
    imgs = cell["images"]
    if imgs == "all":
        pass
    else:
        assert isinstance(imgs, list) and imgs, "%s: bad images" % cid
        for nm in imgs:
            assert nm in valid_names, ("%s: image %r not in imageset %r"
                                       % (cid, nm, imageset))
    assert isinstance(cell["tau"], float) and cell["tau"] > 0, cid
    assert isinstance(cell["sigma_b"], float) and cell["sigma_b"] >= 0, cid
    assert cell["select_rule"] in ("discrepancy", "legacy"), cid
    det = cell.get("det", {})
    assert set(det).issubset(ALLOWED_DET), "%s: det keys %r" % (cid, set(det))
    if "start_mode" in det:
        assert det["start_mode"] in ("active", "delayed", "continuous"), cid
    est = cell.get("est", {})
    assert set(est).issubset(ALLOWED_EST), "%s: est keys %r" % (cid, set(est))
    if "dark_known" in est:
        assert isinstance(est["dark_known"], bool), cid
    # no frozen-provenance leak: C0/c_force must NOT be baked into a cell
    assert "C0" not in cell and "c_force" not in cell, \
        "%s: C0/c_force must not be in a cell (C0_FROZEN.json is authoritative)" % cid


# --------------------------------------------------------------------------- #
# cost model + sharding
# --------------------------------------------------------------------------- #
def load_calibration():
    """Per-(nu,arm) mean wall seconds from the S1 runtime calibration.  The F1
    cost model is the Pass-B calibration (REAL 64^2 / M=4096 walls, INCLUDING the
    analytic lam_TV selection + descriptive audit) — no smoke x16 extrapolation.
    Falls back to the older S1 / smoke calibrations (with the historic x16 32^2
    rescale) only if Pass B is absent.  Returns (calib, meta) where calib maps
    (nu, arm) -> per-(image,arm) wall seconds at M=4096."""
    candidates = ["results/round63_s1_passB/runtime_calibration.json",
                  "results/round63_s1/runtime_calibration.json",
                  "results/round63_s1_smoke/runtime_calibration.json"]
    for rel in candidates:
        ap = os.path.join(ROOT, rel.replace("/", os.sep))
        if not os.path.exists(ap):
            continue
        with open(ap) as f:
            raw = json.load(f)
        # detect smoke / side via the sibling pilot_summary.json config
        side, smoke = 64, False
        sib = os.path.join(os.path.dirname(ap), "pilot_summary.json")
        if os.path.exists(sib):
            with open(sib) as f:
                cfg = (json.load(f) or {}).get("config", {})
            side = int(cfg.get("side", 64))
            smoke = bool(cfg.get("smoke", side == 32))
        else:
            smoke = ("smoke" in rel)
            side = 32 if smoke else 64
        scale = 16.0 if (smoke or side == 32) else 1.0
        calib = {}
        for key, v in raw.items():
            nu_s, arm = key.split("|", 1)
            calib[(float(nu_s), arm)] = float(v["wall_mean_s"]) * scale
        return calib, {"path": rel, "smoke": smoke, "side": side,
                       "scale": scale, "n_keys": len(calib)}
    return {}, {"path": None, "smoke": None, "side": None, "scale": 1.0,
                "n_keys": 0}


def arm_cost_4096(nu, arm, calib):
    """Per-(image,arm) wall seconds at the M=4096 calibration reference."""
    if (float(nu), arm) in calib:
        return calib[(float(nu), arm)]
    if arm in LINEAR_ARMS:
        base = BASE_LINEAR
    elif arm == "RQL":
        base = BASE_RQL
    else:
        base = BASE_ITER
    return base * (float(nu) / NU_REF) ** NU_EXP


def cell_cost_s(cell, calib):
    """Wall-second cost estimate for one cell = n_images * sum_arms(arm_cost@4096)
    * pattern factor * (M / 4096) — count-likelihood fit wall scales ~linearly in
    the number of measurement rows M."""
    n_img = len(expand_images(cell))
    pf = PAT_COST_FACTOR[cell["pattern"]]
    mscale = cell["M"] / M_REF
    per = sum(arm_cost_4096(cell["nu"], a, calib) for a in cell["arms"])
    return n_img * per * pf * mscale


def cell_sort_key(cell):
    det = json.dumps(cell.get("det", {}), sort_keys=True)
    est = json.dumps(cell.get("est", {}), sort_keys=True)
    pat_rank = {"bern50": 0, "hadpair": 1, "gam4": 2}[cell["pattern"]]
    return (cell["side"], pat_rank, cell["M"], cell["rho_bar"], cell["nu"],
            cell["seed"], det, est)


def pack_shards(stage, cells, target_s, calib):
    """Deterministic greedy next-fit packing over the canonically-sorted cells:
    accumulate into the current shard until the next cell would exceed
    target_s (a single over-target cell still forms its own shard).  Assigns
    shard_id '{stage}_{index:02d}', account by index parity (even->pro1,
    odd->pro2), and cell_ids '{stage}_c{n:04d}' in sorted order."""
    cells = sorted(cells, key=cell_sort_key)
    for i, c in enumerate(cells):
        c["cell_id"] = "%s_c%04d" % (stage, i + 1)
    shards, cur, cur_cost = [], [], 0.0
    for c in cells:
        cc = cell_cost_s(c, calib)
        if cur and cur_cost + cc > target_s:
            shards.append(cur)
            cur, cur_cost = [], 0.0
        cur.append(c)
        cur_cost += cc
    if cur:
        shards.append(cur)
    out = []
    for idx, group in enumerate(shards):
        sid = "%s_%02d" % (stage, idx)
        est_s = sum(cell_cost_s(c, calib) for c in group)
        blocked = group[0].get("blocked_by")
        out.append({
            "shard_id": sid,
            "stage": stage,
            "priority": STAGE_PRIORITY[stage],
            "account_hint": "pro1" if idx % 2 == 0 else "pro2",
            "est_hours": round(est_s / 3600.0, 3),
            "blocked_by": blocked,
            "cells": group,
        })
    return out


STAGE_PRIORITY = {"S2A_DETAIL": 1, "S2A_NATURAL": 2, "S2B": 3,
                  "S2C": 3, "S2C128": 3, "S3": 2, "S3JIT": 4, "S3CONT": 2}
STAGE_BLOCKED = {"S2C128": "matrixfree_128",
                 "S3CONT": "continuous_eta_inheritance"}


# --------------------------------------------------------------------------- #
# manifest / summary / expected-cells writers
# --------------------------------------------------------------------------- #
def shard_frozen_inputs(shard):
    """Union of image PNGs consumed by the shard's cells, sorted, with sha256
    (None -> '' for not-yet-built inputs, e.g. the 128^2 block).  A shard's cells
    are all one stage, hence one imageset, so the cache root is uniform."""
    side = shard["cells"][0]["side"]
    imageset = shard["cells"][0]["imageset"]
    names = set()
    for c in shard["cells"]:
        for nm in expand_images(c):
            names.add(nm)
    fi = []
    for nm in sorted(names, key=img_index):
        sha = sha256_png(imageset, side, nm)
        fi.append({"path": "%s/%d/%s.png" % (IMAGESET_ROOT[imageset], side, nm),
                   "sha256": sha if sha is not None else ""})
    return fi


def write_manifest(shard, mdir):
    man = {
        "shard_id": shard["shard_id"],
        "stage": shard["stage"],
        "priority": shard["priority"],
        "account_hint": shard["account_hint"],
        "est_hours": shard["est_hours"],
        "frozen_inputs": shard_frozen_inputs(shard),
        "cells": shard["cells"],
        "output_csv": "results/round63/shards/%s.csv" % shard["shard_id"],
        "meta_json": "results/round63/shards/%s_meta.json" % shard["shard_id"],
    }
    if shard["blocked_by"]:
        man["blocked_by"] = shard["blocked_by"]
    path = os.path.join(mdir, shard["shard_id"] + ".json")
    with open(path, "w", newline="\n") as f:
        json.dump(man, f, indent=2, sort_keys=True, ensure_ascii=True)
        f.write("\n")


def mismatch_sig(cell):
    """Compact deterministic det/est signature (empty for the no-mismatch
    baseline).  REQUIRED in the expected-cell table: without it, S3 mismatch
    rows project onto the identical (side,pattern,rho,nu,M,seed,image,arm)
    tuple as their S2-A baseline and as each other, so the merge-time
    duplicate/missing gate cannot key them.  See LOUD DECISION note."""
    parts = []
    for pre, d in (("det", cell.get("det", {})), ("est", cell.get("est", {}))):
        for k in sorted(d):
            v = d[k]
            v = int(v) if isinstance(v, bool) else (g(v) if isinstance(v, float)
                                                    else v)
            parts.append("%s.%s=%s" % (pre, k, v))
    return ";".join(parts)


def expected_rows(shard):
    rows = []
    for c in shard["cells"]:
        ms = mismatch_sig(c)
        for nm in expand_images(c):
            for arm in c["arms"]:               # already in ARM_RANK order
                rows.append({
                    "shard": shard["shard_id"], "stage": c["stage"],
                    "side": c["side"], "pattern": c["pattern"],
                    "rho_bar": g(c["rho_bar"]), "nu": g(c["nu"]),
                    "M": c["M"], "seed": c["seed"],
                    "image": nm, "arm": arm, "imageset": c["imageset"],
                    "mismatch": ms,
                    "blocked_by": c.get("blocked_by", "") or "",
                })
    return rows


EXP_FIELDS = ["shard", "stage", "side", "pattern", "rho_bar", "nu", "M",
              "seed", "image", "arm", "imageset", "mismatch", "blocked_by"]


def write_expected_cells(all_shards, path):
    rows = []
    for sh in all_shards:
        rows.extend(expected_rows(sh))
    rows.sort(key=lambda r: (r["stage"], r["shard"], r["side"], r["pattern"],
                             float(r["rho_bar"]), float(r["nu"]), r["M"],
                             r["seed"], img_index(r["image"]), ARM_RANK[r["arm"]],
                             r["mismatch"]))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=EXP_FIELDS)
        w.writeheader()
        w.writerows(rows)
    return len(rows)


def write_summary(default_shards, blocked_shards, stage_stats, calib_meta,
                  target_hours, path):
    L = []
    L.append("# ROUND63 shard manifest summary (F1 freeze artifact, D2.3)\n")
    L.append("Generator: `code/round63/make_manifests.py` (deterministic; "
             "no wall-clock, no RNG). Target shard size: %.2f h.\n" % target_hours)
    if calib_meta["path"]:
        L.append("Cost model: S1 runtime calibration `%s` "
                 "(smoke=%s, side=%s, scale x%g, %d (nu,arm) keys) — REAL "
                 "64^2/M=4096 walls INCLUDING the analytic lam_TV selection + "
                 "descriptive audit; M-scaled by (M/4096); analytic fallback "
                 "(RQL %gs / iter %gs / linear %gs @ nu=500, x(nu/500)^%.2f), "
                 "hadpair x2. 128^2 est is a nominal (M-scaled, blocked) figure — "
                 "the matrix-free 128^2 operator is a prerequisite gate and its "
                 "real per-fit cost is substantially higher than this lower bound.\n"
                 % (calib_meta["path"], calib_meta["smoke"], calib_meta["side"],
                    calib_meta["scale"], calib_meta["n_keys"], BASE_RQL,
                    BASE_ITER, BASE_LINEAR, NU_EXP))
    else:
        L.append("Cost model: analytic fallback only (no S1 calibration found).\n")
    L.append("Cohorts (GPT round-7 ruling §1): **S2A_DETAIL** = DETAIL-24 primary "
             "confirmatory positive-speed cohort (imageset detail24, 24 targets); "
             "**S2A_NATURAL** = STL-test NATURAL-24 secondary regime-boundary "
             "cohort (imageset conf, 30 names, NO positive gate). Both run the "
             "identical 64^2/bern50/M=4096 physics grid (rho in {0.05,0.3,0.6,1,2} "
             "x nu in {5,10,20,50,100,200,500,1000,2000} x 5 seeds).\n")
    L.append("")

    L.append("## Default (runnable) shards\n")
    L.append("| shard | stage | pri | cells | est h | account | blocked |")
    L.append("|---|---|---|---|---|---|---|")
    for sh in default_shards:
        L.append("| %s | %s | %d | %d | %.2f | %s | %s |"
                 % (sh["shard_id"], sh["stage"], sh["priority"],
                    len(sh["cells"]), sh["est_hours"], sh["account_hint"],
                    sh["blocked_by"] or ""))
    L.append("")

    if blocked_shards:
        L.append("## Blocked shards (EXCLUDED from the default run set)\n")
        L.append("| shard | stage | pri | cells | est h | account | blocked_by |")
        L.append("|---|---|---|---|---|---|---|")
        for sh in blocked_shards:
            L.append("| %s | %s | %d | %d | %.2f | %s | %s |"
                     % (sh["shard_id"], sh["stage"], sh["priority"],
                        len(sh["cells"]), sh["est_hours"], sh["account_hint"],
                        sh["blocked_by"] or ""))
        L.append("")

    L.append("## Per-stage totals\n")
    L.append("| stage | pri | shards | cells | rows | est h | blocked |")
    L.append("|---|---|---|---|---|---|---|")
    for st in sorted(stage_stats, key=lambda s: (STAGE_PRIORITY[s], s)):
        d = stage_stats[st]
        L.append("| %s | %d | %d | %d | %d | %.2f | %s |"
                 % (st, STAGE_PRIORITY[st], d["shards"], d["cells"], d["rows"],
                    d["est_hours"], STAGE_BLOCKED.get(st, "")))
    L.append("")

    # account balance (default shards only)
    acc = {"pro1": 0.0, "pro2": 0.0}
    accn = {"pro1": 0, "pro2": 0}
    for sh in default_shards:
        acc[sh["account_hint"]] += sh["est_hours"]
        accn[sh["account_hint"]] += 1
    L.append("## Account parity (frozen: even shard index -> pro1, odd -> pro2; "
             "default shards)\n")
    L.append("| account | shards | est h |")
    L.append("|---|---|---|")
    for a in ("pro1", "pro2"):
        L.append("| %s | %d | %.2f |" % (a, accn[a], acc[a]))
    L.append("")

    L.append("## Notes\n")
    L.append("- S3 (+S3JIT/S3CONT) stays at the historic mismatch anchor "
             "**64^2 / bern50 / M=2048 / nu=500 (M/n=0.5)** per spec §3 — it is "
             "NOT lifted to the S2-A M=4096 grid.")
    L.append("- S2B drops its M=4096 rows: they are a strict subset of the "
             "S2A_NATURAL anchor (same physics, subset of images/arms/seeds); "
             "S2B therefore runs only M in {1024, 2048}.")
    L.append("- No cell carries C0/c_force; the frozen `C0_FROZEN.json` "
             "(C0*=inf, Pass A) is authoritative at run time.")
    L.append("")

    L.append("## S4 (EXCLUDED from manifests)\n")
    L.append("- S4(i) scalar exact-Fisher map (nu=20..2000, rho adaptive to "
             "max(64, 2 rho*(nu))): standalone analytic-derivative script, not a "
             "cell grid.")
    L.append("- S4(ii) 8^2/16^2 exact-vs-RQL image split (rho in "
             "{.03,.1,.3,.6,1,2} x nu in {20,100,500,2000}, 3 seeds): BLOCKED on "
             "the EXACT+TV arm (same TV + selector as RQL); emit once that arm "
             "exists.")
    L.append("- S4/Supplement optional: nu in {1,2} exact Bernoulli/trinomial vs "
             "RQL + ceiling-count fraction (outside the deployment zone; no nu<5 "
             "result enters either confirmatory endpoint).")
    L.append("")
    with open(path, "w", newline="\n") as f:
        f.write("\n".join(L))


# --------------------------------------------------------------------------- #
# selftest: ONE tiny cell through the REAL run_cell (proves dict shape accepted)
# --------------------------------------------------------------------------- #
def selftest():
    # side=32 (not 16): the "text" structural target renders empty at 16^2
    # (images63 raises before GI ever runs) — a pre-existing images63 limitation,
    # unrelated to this generator.  A NON-CONFIRMATORY image is used on purpose:
    # imageset "dev" + dev_stl_00 is the S1 development set (STL-10 TRAIN split),
    # disjoint from every S2 confirmatory cohort (NATURAL-24 and DETAIL-24), so
    # the selftest never reconstructs a frozen confirmatory target.  It proves the
    # exact same thing: a make_cell()-shaped dict is accepted by the real run_cell.
    from campaign import run_cell
    cell = {"cell_id": "selftest", "side": 32, "pattern": "bern50",
            "rho_bar": 0.05, "nu": 100.0, "M": 64, "seed": 0, "arms": ["GI"],
            "images": ["dev_stl_00"], "imageset": "dev", "tau": TAU,
            "sigma_b": 0.0, "select_rule": "legacy"}
    validate_cell(cell)
    rows = run_cell(cell)
    ok = (len(rows) == 1 and rows[0]["arm"] == "GI"
          and set(rows[0]) >= {"side", "pattern", "rho_bar", "nu", "M", "seed",
                               "image", "arm", "PSNR", "SSIM"})
    print("[selftest] run_cell accepted the cell dict: rows=%d arm=%s PSNR=%.3f  %s"
          % (len(rows), rows[0]["arm"], rows[0]["PSNR"], "PASS" if ok else "FAIL"),
          flush=True)
    return 0 if ok else 1


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--target-hours", type=float, default=6.0,
                    help="greedy shard target size in wall-hours (default 6.0)")
    ap.add_argument("--selftest", action="store_true",
                    help="also push ONE tiny GI cell through campaign.run_cell")
    args = ap.parse_args(argv)

    if args.selftest:
        return selftest()

    target_s = args.target_hours * 3600.0
    calib, calib_meta = load_calibration()

    # 1) build every stage's cells
    groups = [
        ("S2A_DETAIL", build_s2a_detail()),
        ("S2A_NATURAL", build_s2a_natural()),
        ("S2B", build_s2b()),
        ("S2C", build_s2c()),
        ("S2C128", build_s2c_128_blocked()),
    ]
    s3, s3jit, s3cont = build_s3()
    groups += [("S3", s3), ("S3JIT", s3jit), ("S3CONT", s3cont)]

    # 2) pack + validate
    all_shards = []
    stage_stats = {}
    for stage, cells in groups:
        shards = pack_shards(stage, cells, target_s, calib)
        for sh in shards:
            for c in sh["cells"]:
                validate_cell(c)
        n_rows = sum(len(expected_rows(sh)) for sh in shards)
        stage_stats[stage] = {
            "shards": len(shards), "cells": len(cells), "rows": n_rows,
            "est_hours": round(sum(s["est_hours"] for s in shards), 3),
        }
        all_shards.extend(shards)

    default_shards = [s for s in all_shards if not s["blocked_by"]]
    blocked_shards = [s for s in all_shards if s["blocked_by"]]

    # 3) write outputs (wipe stale shard manifests from a previous design first)
    mdir = os.path.join(ROOT, "results", "round63", "manifests")
    os.makedirs(mdir, exist_ok=True)
    for stale in glob.glob(os.path.join(mdir, "*.json")):
        os.remove(stale)
    for sh in all_shards:
        write_manifest(sh, mdir)

    n_exp = write_expected_cells(all_shards,
                                 os.path.join(ROOT, "results", "round63",
                                              "expected_cells.csv"))
    write_summary(default_shards, blocked_shards, stage_stats, calib_meta,
                  args.target_hours, os.path.join(mdir, "MANIFEST_SUMMARY.md"))

    index = {
        "target_hours": args.target_hours,
        "cost_model": calib_meta,
        "n_expected_rows": n_exp,
        "default_shards": [{"shard_id": s["shard_id"], "stage": s["stage"],
                            "priority": s["priority"],
                            "account_hint": s["account_hint"],
                            "est_hours": s["est_hours"],
                            "n_cells": len(s["cells"])} for s in default_shards],
        "blocked_shards": [{"shard_id": s["shard_id"], "stage": s["stage"],
                            "blocked_by": s["blocked_by"],
                            "account_hint": s["account_hint"],
                            "est_hours": s["est_hours"],
                            "n_cells": len(s["cells"])} for s in blocked_shards],
        "stage_stats": stage_stats,
    }
    with open(os.path.join(mdir, "MANIFEST_INDEX.json"), "w", newline="\n") as f:
        json.dump(index, f, indent=2, sort_keys=True, ensure_ascii=True)
        f.write("\n")

    # 4) stdout digest
    tot_cells = sum(len(s["cells"]) for s in all_shards)
    tot_def_h = sum(s["est_hours"] for s in default_shards)
    tot_blk_h = sum(s["est_hours"] for s in blocked_shards)
    print("[make_manifests] target=%.2fh  cost=%s"
          % (args.target_hours,
             calib_meta["path"] or "analytic-fallback-only"), flush=True)
    for st in sorted(stage_stats, key=lambda s: (STAGE_PRIORITY[s], s)):
        d = stage_stats[st]
        print("  %-12s pri%d  shards=%3d cells=%4d rows=%6d est=%.1fh %s"
              % (st, STAGE_PRIORITY[st], d["shards"], d["cells"], d["rows"],
                 d["est_hours"], "[BLOCKED %s]" % STAGE_BLOCKED[st]
                 if st in STAGE_BLOCKED else ""), flush=True)
    print("  ---- default shards=%d (%.1fh)  blocked shards=%d (%.1fh)  "
          "total cells=%d  expected rows=%d"
          % (len(default_shards), tot_def_h, len(blocked_shards), tot_blk_h,
             tot_cells, n_exp), flush=True)
    print("[make_manifests] wrote %d manifests + MANIFEST_SUMMARY.md + "
          "MANIFEST_INDEX.json + expected_cells.csv under results/round63/"
          % len(all_shards), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
