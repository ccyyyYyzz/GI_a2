"""M1 shard-manifest generator (docs/ROUND63_METHOD_SPEC_M1.md §10; modeled
on make_manifests.py, emitting the SAME manifest JSON schema shard_runner.py
consumes — shard_id, stage, priority, account_hint, est_hours, frozen_inputs
(sha256), cells, output_csv, meta_json, blocked_by — so the proven Colab
infra (session_driver.sh -> remote_lane.py -> shard_runner.py) runs the M1
grid unchanged (cells resolve through campaign.run_cell via the M1 appendix).

Sharding: the full grid is 6 arms x 24 m1 images x 5 seeds x {safe 0.05,
arm-fast} x 9 nu = 12,960 single-image cells. The atomic unit is the
(arm, image, seed) BLOCK (18 cells): blocks are never split, so every design
cache is computed exactly once and a shard is self-contained. Blocks are
packed per arm (spec §10 "shard by (arm x scene block)") by a deterministic
greedy into ~TARGET_SHARDS shards balanced by estimated cost; OED/RIDGE
design cost is amortized into the block's estimate (the design npz caches
are declared frozen_inputs, i.e. they ship with the bundle and are built
LOCALLY by `m1_runner.py --designs` before manifest generation).

Cost model (deterministic constants, calibrated from the m1_runner --smoke
timings on the local reference CPU; documented, not fitted at run time):
  cell cost  c(nu, rho) = C0_S + C1_S * nu * rho   [simulation + RQL solve]
  design amortization   per block: OED-DT ~ DT_S, OED-EQLOAD ~ EQ_S,
                        RIDGE-FIXED ~ RF_S (xhat + 9 dwell calibrations;
                        ridge targets are global, cached once)

Blocks whose design caches are missing are emitted with
blocked_by = "m1_design_cache_missing" (RIDGE-FIXED fast rho_bar then
carries the sentinel -1.0): they are NOT runnable until `--designs`
completes and this generator is re-run (A13).

Usage:
  python code/round63/m1_make_manifests.py [--target-shards 40]
"""
import argparse
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
ROOT = os.path.dirname(os.path.dirname(HERE))

import m1_runner as m1                      # noqa: E402

TARGET_SHARDS = 40
STAGE_PRIORITY = {"M1_OED-DT": 1, "M1_OED-EQLOAD": 2, "M1_RIDGE-FIXED": 2,
                  "M1_LBLOB16": 3, "M1_SCAT16": 3, "M1_SCAT32": 4,
                  "M1_MATCH1": 3}

# ---- deterministic cost constants (seconds) -------------------------------- #
# Calibrated from m1_runner --smoke (2026-07-19, local reference CPU):
# OED cells 0.2-0.3 s (nu 20 -> 2000), RIDGE cell 0.5 s, SCAT16 first cell
# 2.4 s (sparsek 2-switch build; cached per (image, seed) block thereafter);
# OED-DT design 95 s, RIDGE-FIXED design (9 dwell calibrations) 4.3 s.
# NOTE: at these rates the FULL campaign is ~6-7 h single-process, i.e.
# LOCAL-feasible with 2-3 processes without Colab; the shards remain valid
# for either backend.
C0_S = 0.30           # per-cell floor: RQL + select_eta at 32^2, M=1024
C1_S = 8.4e-5         # per-cell nu*rho slope (count simulation + likelihood)
DESIGN_AMORT_S = {"OED-DT": 100.0, "OED-EQLOAD": 60.0, "RIDGE-FIXED": 8.0,
                  "SCAT16": 2.5, "SCAT32": 4.0, "LBLOB16": 1.0,
                  "MATCH1": 30.0}

MDIR = os.path.join(ROOT, "results", "round63_m1", "manifests")
SHARD_CSV_DIR = "results/round63_m1/shards"
IMG_ROOT_REL = "data/r63_images_m1"


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def cell_cost_s(cell):
    return C0_S + C1_S * float(cell["nu"]) * max(float(cell["rho_bar"]), 0.0)


def block_cells(arm, image, seed, C0):
    """The 18 cells of one (arm, image, seed) block; RIDGE-FIXED fast rho
    from the design cache (sentinel -1.0 + blocked flag when missing)."""
    cells, blocked = [], None
    for nu in m1.NU_FULL:
        cells.append(m1.m1_cell(arm, image, seed, m1.RHO_SAFE, nu, None))
    for nu in m1.NU_FULL:
        rf = (m1.fast_rho_for(arm, nu, image, seed)
              if arm == "RIDGE-FIXED" else m1.RHO_FAST)
        if rf is None:
            rf = -1.0
            blocked = "m1_design_cache_missing"
        cells.append(m1.m1_cell(arm, image, seed, rf, nu, None))
    if arm in ("OED-DT", "OED-EQLOAD", "MATCH1") and not os.path.exists(
            m1._design_path(image, seed, arm)):
        blocked = "m1_design_cache_missing"
    for c in cells:
        c.pop("C0", None)               # frozen C0_FROZEN.json rules at run time
    return cells, blocked


def build_blocks():
    C0 = None
    blocks = []
    for arm in m1.ARMS_ALL:
        for image in m1.m1_images():
            for seed in m1.SEEDS5:
                cells, blocked = block_cells(arm, image, seed, C0)
                cost = sum(cell_cost_s(c) for c in cells) \
                    + DESIGN_AMORT_S[arm]
                blocks.append({"arm": arm, "image": image, "seed": seed,
                               "cells": cells, "cost_s": cost,
                               "blocked_by": blocked})
    return blocks


def pack(blocks, target_shards):
    """Per-arm deterministic greedy: shards allocated to arms proportional to
    arm cost; within an arm, blocks (sorted by cost desc, then image, seed)
    go to the lightest shard (ties -> lowest shard index)."""
    total = sum(b["cost_s"] for b in blocks)
    shards = []
    for arm in m1.ARMS_ALL:
        ab = [b for b in blocks if b["arm"] == arm]
        acost = sum(b["cost_s"] for b in ab)
        n_sh = max(1, int(round(target_shards * acost / total)))
        bins = [{"stage": "M1_%s" % arm, "blocks": [], "cost_s": 0.0}
                for _ in range(n_sh)]
        for b in sorted(ab, key=lambda b: (-b["cost_s"], b["image"], b["seed"])):
            j = min(range(n_sh), key=lambda i: (bins[i]["cost_s"], i))
            bins[j]["blocks"].append(b)
            bins[j]["cost_s"] += b["cost_s"]
        shards.extend(bins)
    out = []
    for idx, sh in enumerate(shards):
        blocked = sorted({b["blocked_by"] for b in sh["blocks"]
                          if b["blocked_by"]})
        cells = []
        for b in sorted(sh["blocks"], key=lambda b: (b["image"], b["seed"])):
            cells.extend(b["cells"])
        sid = "%s_%02d" % (sh["stage"], sum(1 for s in out
                                            if s["stage"] == sh["stage"]))
        out.append({"shard_id": sid, "stage": sh["stage"],
                    "priority": STAGE_PRIORITY[sh["stage"]],
                    "account_hint": "pro1" if idx % 2 == 0 else "pro2",
                    "est_hours": round(sh["cost_s"] / 3600.0, 3),
                    "blocked_by": (blocked[0] if blocked else None),
                    "blocks": sh["blocks"], "cells": cells})
    return out


def shard_frozen_inputs(shard):
    ents = []
    seen = set()
    for b in shard["blocks"]:
        png = "%s/%d/%s.png" % (IMG_ROOT_REL, m1.SIDE, b["image"])
        if png not in seen:
            seen.add(png)
            p = os.path.join(ROOT, png)
            ents.append({"path": png,
                         "sha256": _sha256(p) if os.path.exists(p) else ""})
        if b["arm"] in m1.DESIGN_ARMS:
            dp = m1._design_path(b["image"], b["seed"], b["arm"])
            rel = os.path.relpath(dp, ROOT).replace("\\", "/")
            if rel not in seen:
                seen.add(rel)
                ents.append({"path": rel,
                             "sha256": _sha256(dp) if os.path.exists(dp)
                             else ""})
    return sorted(ents, key=lambda e: e["path"])


def write_manifests(shards, mdir=MDIR):
    os.makedirs(mdir, exist_ok=True)
    for i, sh in enumerate(shards):
        for j, c in enumerate(sh["cells"]):
            c["cell_id"] = "%s_c%04d" % (sh["shard_id"], j + 1)
        man = {"shard_id": sh["shard_id"], "stage": sh["stage"],
               "priority": sh["priority"], "account_hint": sh["account_hint"],
               "est_hours": sh["est_hours"],
               "frozen_inputs": shard_frozen_inputs(sh),
               "cells": sh["cells"],
               "output_csv": "%s/%s.csv" % (SHARD_CSV_DIR, sh["shard_id"]),
               "meta_json": "%s/%s_meta.json" % (SHARD_CSV_DIR,
                                                 sh["shard_id"])}
        if sh["blocked_by"]:
            man["blocked_by"] = sh["blocked_by"]
        with open(os.path.join(mdir, sh["shard_id"] + ".json"), "w",
                  newline="\n") as f:
            json.dump(man, f, indent=2, sort_keys=True, ensure_ascii=True)
            f.write("\n")

    runnable = [s for s in shards if not s["blocked_by"]]
    blocked = [s for s in shards if s["blocked_by"]]
    index = {
        "cost_model": {"C0_S": C0_S, "C1_S": C1_S,
                       "design_amort_s": DESIGN_AMORT_S,
                       "provenance": "m1_runner --smoke local calibration"},
        "n_expected_rows": sum(len(s["cells"]) for s in shards),
        "default_shards": [{"shard_id": s["shard_id"], "stage": s["stage"],
                            "priority": s["priority"],
                            "account_hint": s["account_hint"],
                            "n_cells": len(s["cells"]),
                            "est_hours": s["est_hours"]} for s in runnable],
        "blocked_shards": [{"shard_id": s["shard_id"], "stage": s["stage"],
                            "account_hint": s["account_hint"],
                            "n_cells": len(s["cells"]),
                            "est_hours": s["est_hours"],
                            "blocked_by": s["blocked_by"]} for s in blocked],
    }
    with open(os.path.join(mdir, "MANIFEST_INDEX.json"), "w",
              newline="\n") as f:
        json.dump(index, f, indent=2, sort_keys=True)
        f.write("\n")

    lines = ["# M1 MANIFEST SUMMARY", "",
             "| shard | stage | pri | acct | cells | est_h | blocked |",
             "|---|---|---|---|---|---|---|"]
    for s in shards:
        lines.append("| %s | %s | %d | %s | %d | %.3f | %s |"
                     % (s["shard_id"], s["stage"], s["priority"],
                        s["account_hint"], len(s["cells"]), s["est_hours"],
                        s["blocked_by"] or ""))
    tot = sum(s["est_hours"] for s in shards)
    lines += ["", "total shards: %d  cells: %d  est: %.1f h (sum of shards)"
              % (len(shards), sum(len(s["cells"]) for s in shards), tot)]
    with open(os.path.join(mdir, "MANIFEST_SUMMARY.md"), "w",
              newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    return index, lines


def write_expected_cells(shards, mdir=MDIR):
    """R13 item 13: the complete expected-cell ledger (one row per cell:
    shard, stage, arm, image, seed, rho_bar, nu, pattern, blocked_by)."""
    import csv
    os.makedirs(mdir, exist_ok=True)
    path = os.path.join(mdir, "expected_cells.csv")
    n = 0
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["shard_id", "stage", "arm", "image", "seed", "rho_bar",
                    "nu", "pattern", "cell_id", "blocked_by"])
        for sh in shards:
            for b in sorted(sh["blocks"], key=lambda b: (b["image"],
                                                         b["seed"])):
                for c in b["cells"]:
                    w.writerow([sh["shard_id"], sh["stage"], b["arm"],
                                b["image"], b["seed"], c["rho_bar"], c["nu"],
                                c["pattern"], c.get("cell_id", ""),
                                b["blocked_by"] or ""])
                    n += 1
    return n


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-shards", type=int, default=TARGET_SHARDS)
    a = ap.parse_args(sys.argv[1:] if argv is None else argv)
    blocks = build_blocks()
    shards = pack(blocks, a.target_shards)
    index, lines = write_manifests(shards)
    print("\n".join(lines))
    print("\n[m1 manifests] wrote %d shard manifests -> %s"
          % (len(shards), MDIR))
    n_blk = len(index["blocked_shards"])
    if n_blk:
        print("[m1 manifests] %d shard(s) BLOCKED (design caches missing; "
              "run `m1_runner.py --designs` then regenerate)." % n_blk)
    return 0


if __name__ == "__main__":
    sys.exit(main())
