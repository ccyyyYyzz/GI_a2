"""M1 shard-manifest generator — R17 re-architecture.

Emits the shard_runner-schema manifests (shard_id, stage, priority,
account_hint, est_hours, frozen_inputs+sha256, cells, output_csv,
meta_json) for the REVISED architecture (amendment §B/§C):

  Imaging stages (campaign.run_cell cells, 24 img x 5 seeds x 9 nu each):
    M1_SCAT32-SAFE   load 0.05            (endpoint arm, priority 1)
    M1_SCAT32-060    load 0.60            (endpoint arm, priority 1)
    M1_RIDGE-SCAT32  runtime-calibrated   (endpoint arm, priority 1;
                     cells carry rho_bar=-1 sentinel + m1_ridge_dynamic —
                     the campaign M1 appendix resolves the achieved ridge
                     load on the VM from the common pre-scan estimate)
    M1_SCAT16 / M1_LBLOB16  load 0.60     (context, descriptive, priority 3)

  Certificate stage (amendment §C.3; run_cert_cell via the appendix):
    M1_CERT — 480 cells: 24 scenes x 5 seeds x nu{200,2000} x b{0.05,0.60},
    grouped 4-per-(image,seed) block, priority 2.

R17 retirements honored: NO designs/<image>_<seed>_<ARM>.npz frozen inputs
(the fleet audit found 36/40 old shards declaring them); the only design
artifact any shard needs is the scene-independent deployed SCAT32 multiset
(results/round63_m1/designs/scat32_deployed.npz, sha-pinned) plus the image
PNGs. Ridge calibration and certificate subspaces are per-cell runtime
computation.

Deterministic; pure generation (no reconstruction).
Usage:  python code/round63/m1_make_manifests.py [--target-shards 40]
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
STAGE_PRIORITY = {"M1_SCAT32-SAFE": 1, "M1_SCAT32-060": 1,
                  "M1_RIDGE-SCAT32": 1, "M1_CERT": 2,
                  "M1_SCAT16": 3, "M1_LBLOB16": 3}
DESCRIPTIVE_STAGES = {"M1_SCAT16", "M1_LBLOB16"}

# ---- deterministic cost constants (seconds) -------------------------------- #
# Imaging: smoke-calibrated (freeze-era measurements, 32^2 RQL cells).
# Cert: measured by the R17 selftest item 6 (per-cell wall time recorded in
# the ledger; constant below updated from that measurement).
C0_S = 0.35
C1_S = 8.4e-5
CERT_CELL_S = 240.0
RIDGE_CAL_S = 25.0        # once per (image, seed) block on the VM

MDIR = os.path.join(ROOT, "results", "round63_m1", "manifests")
SHARD_CSV_DIR = "results/round63_m1/shards"
IMG_ROOT_REL = "data/r63_images_m1"
SCAT32_REL = "results/round63_m1/designs/scat32_deployed.npz"


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def imaging_cell_cost(arm, nu):
    rho = m1.MODE_LOADS.get(arm)
    if rho is None:                        # RIDGE-SCAT32: ridge-load slope
        rt = m1.ridge_targets([nu])["%g" % nu]
        rho = rt["rho_R_production"]
    return C0_S + C1_S * float(nu) * float(rho)


def build_blocks():
    blocks = []
    for arm in m1.ARMS_ALL:
        for image in m1.m1_images():
            for seed in m1.SEEDS5:
                cells = [m1.m1_cell(arm, image, seed,
                                    m1.MODE_LOADS.get(arm, -1.0), nu)
                         for nu in m1.NU_FULL]
                for c in cells:
                    c.pop("C0", None)
                cost = sum(imaging_cell_cost(arm, nu) for nu in m1.NU_FULL)
                if arm == "RIDGE-SCAT32":
                    cost += RIDGE_CAL_S
                blocks.append({"kind": "imaging", "arm": arm, "image": image,
                               "seed": seed, "cells": cells, "cost_s": cost,
                               "stage": "M1_%s" % arm})
    for image in m1.m1_images():
        for seed in m1.SEEDS5:
            cells = [m1.cert_cell(image, seed, nu, b)
                     for nu in m1.CERT_NUS for b in m1.CERT_BUDGETS]
            blocks.append({"kind": "cert", "arm": "CERT", "image": image,
                           "seed": seed, "cells": cells,
                           "cost_s": len(cells) * CERT_CELL_S,
                           "stage": "M1_CERT"})
    return blocks


def pack(blocks, target_shards):
    total = sum(b["cost_s"] for b in blocks)
    stages = sorted({b["stage"] for b in blocks})
    shards = []
    for stage in stages:
        sb = [b for b in blocks if b["stage"] == stage]
        scost = sum(b["cost_s"] for b in sb)
        n_sh = max(1, int(round(target_shards * scost / total)))
        bins = [{"stage": stage, "blocks": [], "cost_s": 0.0}
                for _ in range(n_sh)]
        for b in sorted(sb, key=lambda b: (-b["cost_s"], b["image"],
                                           b["seed"])):
            j = min(range(n_sh), key=lambda i: (bins[i]["cost_s"], i))
            bins[j]["blocks"].append(b)
            bins[j]["cost_s"] += b["cost_s"]
        shards.extend(bins)
    out = []
    for idx, sh in enumerate(shards):
        cells = []
        for b in sorted(sh["blocks"], key=lambda b: (b["image"], b["seed"])):
            cells.extend(b["cells"])
        sid = "%s_%02d" % (sh["stage"], sum(1 for s in out
                                            if s["stage"] == sh["stage"]))
        out.append({"shard_id": sid, "stage": sh["stage"],
                    "priority": STAGE_PRIORITY[sh["stage"]],
                    "descriptive": sh["stage"] in DESCRIPTIVE_STAGES,
                    "account_hint": "pro1" if idx % 2 == 0 else "pro2",
                    "est_hours": round(sh["cost_s"] / 3600.0, 3),
                    "blocks": sh["blocks"], "cells": cells})
    return out


def shard_frozen_inputs(shard):
    """Image PNGs + the deployed SCAT32 multiset for SCAT32-mode/cert
    shards. NO per-(image,seed,arm) design caches (R17 retirement)."""
    ents = []
    seen = set()
    needs_deployed = shard["stage"] in ("M1_SCAT32-SAFE", "M1_SCAT32-060",
                                        "M1_RIDGE-SCAT32", "M1_CERT")
    if needs_deployed:
        p = os.path.join(ROOT, SCAT32_REL)
        ents.append({"path": SCAT32_REL,
                     "sha256": _sha256(p) if os.path.exists(p) else ""})
    for b in shard["blocks"]:
        png = "%s/%d/%s.png" % (IMG_ROOT_REL, m1.SIDE, b["image"])
        if png not in seen:
            seen.add(png)
            p = os.path.join(ROOT, png)
            ents.append({"path": png,
                         "sha256": _sha256(p) if os.path.exists(p) else ""})
    return sorted(ents, key=lambda e: e["path"])


def write_manifests(shards, mdir=MDIR):
    os.makedirs(mdir, exist_ok=True)
    for sh in shards:
        for j, c in enumerate(sh["cells"]):
            c["cell_id"] = "%s_c%04d" % (sh["shard_id"], j + 1)
        man = {"shard_id": sh["shard_id"], "stage": sh["stage"],
               "priority": sh["priority"],
               "descriptive": sh["descriptive"],
               "account_hint": sh["account_hint"],
               "est_hours": sh["est_hours"],
               "frozen_inputs": shard_frozen_inputs(sh),
               "cells": sh["cells"],
               "output_csv": "%s/%s.csv" % (SHARD_CSV_DIR, sh["shard_id"]),
               "meta_json": "%s/%s_meta.json" % (SHARD_CSV_DIR,
                                                 sh["shard_id"])}
        with open(os.path.join(mdir, sh["shard_id"] + ".json"), "w",
                  newline="\n") as f:
            json.dump(man, f, indent=2, sort_keys=True, ensure_ascii=True)
            f.write("\n")
    index = {
        "architecture": "R17 (issue #9): SCAT32 modes + context + 480 cert",
        "cost_model": {"C0_S": C0_S, "C1_S": C1_S,
                       "CERT_CELL_S": CERT_CELL_S,
                       "RIDGE_CAL_S": RIDGE_CAL_S},
        "n_expected_rows": sum(len(s["cells"]) for s in shards),
        "default_shards": [{"shard_id": s["shard_id"], "stage": s["stage"],
                            "priority": s["priority"],
                            "descriptive": s["descriptive"],
                            "account_hint": s["account_hint"],
                            "n_cells": len(s["cells"]),
                            "est_hours": s["est_hours"]} for s in shards],
        "blocked_shards": [],
    }
    with open(os.path.join(mdir, "MANIFEST_INDEX.json"), "w",
              newline="\n") as f:
        json.dump(index, f, indent=2, sort_keys=True)
        f.write("\n")
    lines = ["# M1 MANIFEST SUMMARY (R17 architecture)", "",
             "| shard | stage | pri | desc | acct | cells | est_h |",
             "|---|---|---|---|---|---|---|"]
    for s in shards:
        lines.append("| %s | %s | %d | %s | %s | %d | %.3f |"
                     % (s["shard_id"], s["stage"], s["priority"],
                        "Y" if s["descriptive"] else "",
                        s["account_hint"], len(s["cells"]), s["est_hours"]))
    tot = sum(s["est_hours"] for s in shards)
    lines += ["", "total shards: %d  cells: %d  est: %.1f h"
              % (len(shards), sum(len(s["cells"]) for s in shards), tot)]
    with open(os.path.join(mdir, "MANIFEST_SUMMARY.md"), "w",
              newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    return index, lines


def write_expected_cells(shards, mdir=MDIR):
    """R13/R17 expected-cell ledger: one row per cell (imaging + cert)."""
    import csv
    os.makedirs(mdir, exist_ok=True)
    path = os.path.join(mdir, "expected_cells.csv")
    n = 0
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["shard_id", "stage", "kind", "arm", "image", "seed",
                    "rho_bar_or_b", "nu", "pattern_or_cert", "cell_id",
                    "descriptive"])
        for sh in shards:
            for b in sorted(sh["blocks"], key=lambda b: (b["image"],
                                                         b["seed"])):
                for c in b["cells"]:
                    if c.get("m1_cert"):
                        w.writerow([sh["shard_id"], sh["stage"], "cert",
                                    "CERT", c["image"], c["seed"], c["b"],
                                    c["nu"], "m1_cert", c["cell_id"],
                                    ""])
                    else:
                        w.writerow([sh["shard_id"], sh["stage"], "imaging",
                                    b["arm"], b["image"], b["seed"],
                                    c["rho_bar"], c["nu"], c["pattern"],
                                    c["cell_id"],
                                    int(sh["stage"] in DESCRIPTIVE_STAGES)])
                    n += 1
    return n


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-shards", type=int, default=TARGET_SHARDS)
    a = ap.parse_args(sys.argv[1:] if argv is None else argv)
    blocks = build_blocks()
    shards = pack(blocks, a.target_shards)
    index, lines = write_manifests(shards)
    n_exp = write_expected_cells(shards)
    print("\n".join(lines[-6:]))
    print("[m1 manifests R17] %d shards, %d expected cells -> %s"
          % (len(shards), n_exp, MDIR))
    return 0


if __name__ == "__main__":
    sys.exit(main())
