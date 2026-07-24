#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
SEALED-PROBE SHARD PLANNER (R39 Phase-2 machinery, 5-shard Colab-fleet lane format).

Enumerates the sealed beyond-band confirmatory into CELLS and shards them across the 5-session
fleet with the FROZEN account parity (even shards -> pro1, odd shards -> pro2), writing one
manifest JSON per shard in the exact contract the proven templates consume
(code/round63/colab/{make_bundle,live_launch_all,session_driver,remote_lane}) so the confirmatory
launches on 5 VMs the moment R39 freezes the bars.

CELL TYPES (each a pure function of sealed banks + frozen [R39] budget):
  arms      : (partition in {confirmatory, oracle}, scene_id, replicate) -> all six arms,
              beyond-band NMSE + multi-start agreement std.
  mismatch  : (mismatch scene_id, axis, condition) -> beyond-band NMSE degradation vs matched.
  aperture  : (calibration scene_id) -> aperture-law in/out coverage separation.

This planner WRITES ONLY manifests + a lane plan; it does NOT run the confirmatory (thresholds
are [R39] placeholders until the referee ruling).  Frozen choices are tagged [R39].
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
BANKS_DIR = os.path.join(HERE, "scene_banks")
MAN_DIR = os.path.join(HERE, "manifests")

# ------------------------------------------------- frozen budget + fleet [R39]
FROZEN = dict(
    T=4096,               # [R39] epochs / independent medium realizations per record
    photons=1e5,          # [R39] detected photons per bucket
    replicates=3,         # [R39] Poisson/medium replicates per scene
    steps=4000,           # [R39] moment-matching Adam steps
    n_starts=6,           # [R39] covariance multi-start count (for the agreement bar B4)
    n_em=25,              # [R39] Kalman-EM iterations for the MLE arm
    mismatch_axes={       # [R39] the F5 mismatch axes the sealed probe must include
        "fine_rotation": ["rot10", "rot20"],
        "band_width": ["narrow", "wide"],
        "radial_profile": ["flat_vs_1overf"],
        "geometry_alpha": ["a0.1", "a0.25", "a0.5"],
        "wrong_tau": ["tau_half", "tau_double"],
        "wrong_sigma_f": ["sf_0.2", "sf_0.4"],
    },
)

# 5-session fleet with FROZEN shard<->account parity (even=pro1, odd=pro2) per COLAB guide
FLEET = [
    dict(shard=0, account="pro1", home="/var/tmp/codex-colab-accounts/pro1"),
    dict(shard=1, account="pro2", home="/var/tmp/codex-colab-accounts/pro2"),
    dict(shard=2, account="pro1", home="/var/tmp/codex-colab-accounts/pro1"),
    dict(shard=3, account="pro2", home="/var/tmp/codex-colab-accounts/pro2"),
    dict(shard=4, account="pro2", home="/var/tmp/codex-colab-accounts/pro2"),  # pro2 has 3xL4
]
N_SHARDS = len(FLEET)


def _scene_ids(partition):
    man = json.load(open(os.path.join(BANKS_DIR, partition, "MANIFEST.json")))
    return [e["scene_id"] for e in man["scenes"]]


def enumerate_cells():
    cells = []
    for part in ("confirmatory", "oracle"):
        for sid in _scene_ids(part):
            for rep in range(FROZEN["replicates"]):
                cells.append(dict(kind="arms", partition=part, scene_id=sid, replicate=rep,
                                  cell_id=f"arms::{sid}::r{rep}"))
    for sid in _scene_ids("mismatch"):
        for axis, conds in FROZEN["mismatch_axes"].items():
            for cond in conds:
                cells.append(dict(kind="mismatch", partition="mismatch", scene_id=sid,
                                  axis=axis, condition=cond,
                                  cell_id=f"mismatch::{sid}::{axis}::{cond}"))
    for sid in _scene_ids("calibration"):
        cells.append(dict(kind="aperture", partition="calibration", scene_id=sid,
                          cell_id=f"aperture::{sid}"))
    return cells


def plan(write=True):
    if not os.path.isdir(BANKS_DIR):
        raise SystemExit("scene_banks/ not found -- run bank_gen64.py first.")
    cells = enumerate_cells()
    shards = {i: [] for i in range(N_SHARDS)}
    for k, c in enumerate(cells):                      # round-robin -> balanced shards
        shards[k % N_SHARDS].append(c)
    if write:
        os.makedirs(MAN_DIR, exist_ok=True)
    manifests = []
    for f in FLEET:
        i = f["shard"]
        # frozen parity: even->pro1, odd->pro2; shard4 is the pro2 EXTRA (pro2 has 3 L4 slots)
        parity = ("pro2-extra (5th slot)" if i == 4
                  else ("even->pro1" if i % 2 == 0 else "odd->pro2"))
        man = dict(shard_id=f"r39bb_shard{i}", account=f["account"], home=f["home"],
                   parity=parity,
                   frozen=FROZEN, n_cells=len(shards[i]), cells=shards[i],
                   output_csv=f"results/round63_next/FOG_DMD_PROBE64/shards/r39bb_shard{i}.csv",
                   note="thresholds/budget [R39]-placeholder; confirmatory UNRUN until referee ruling")
        manifests.append(man)
        if write:
            json.dump(man, open(os.path.join(MAN_DIR, f"r39bb_shard{i}.json"), "w"), indent=2)
    plan_blob = dict(n_cells=len(cells), n_shards=N_SHARDS, fleet=FLEET, frozen=FROZEN,
                     cells_per_shard={i: len(shards[i]) for i in range(N_SHARDS)},
                     cell_kind_counts={kk: sum(c["kind"] == kk for c in cells)
                                       for kk in ("arms", "mismatch", "aperture")})
    if write:
        json.dump(plan_blob, open(os.path.join(MAN_DIR, "LANE_PLAN.json"), "w"), indent=2)
    return plan_blob, manifests


if __name__ == "__main__":
    blob, mans = plan(write=True)
    print("=== R39 beyond-band sealed-probe lane plan (5-shard fleet) ===")
    print(f"total cells {blob['n_cells']}  ->  {blob['cell_kind_counts']}")
    for f in FLEET:
        i = f["shard"]
        lab = "pro2-extra" if i == 4 else ("even->pro1" if i % 2 == 0 else "odd->pro2")
        print(f"  shard{i} [{f['account']}] {blob['cells_per_shard'][i]:3d} cells ({lab})")
    print(f"\nwrote {len(mans)} shard manifests + LANE_PLAN.json under manifests/")
    print("confirmatory UNRUN (thresholds [R39]-placeholder); launch only after the R39 ruling.")
