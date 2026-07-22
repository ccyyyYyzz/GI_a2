"""ROUND63 DEV-bridge shard runner (R28-authorized grid).

Consumes a bridge shard manifest (results/round63_bridge/manifests/{layout}/
bridge_{route}.json), verifies its frozen_inputs sha256, and runs each cell =
one (scene, nu, c, seed) group through bridge_harness.run_group (all six arms,
sharing the charged pre-scan; TRUE-X-FW cached per (scene,nu); XHAT-FW seed-0
diagnostic economy).  Resumable: a cell is DONE when its six per-arm JSONs exist
under results/round63_bridge/cells/.  Emits a heartbeat line per cell for the
watchdog (progress / stall detection).  No gate is evaluated here (outcome-blind
production); gate evaluation is a separate post-grid step (bridge_gates).

Usage:
  python code/round63/bridge_shard_runner.py --manifest <path> [--wall-budget-s N]
"""
import argparse
import hashlib
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
ROOT = os.path.dirname(os.path.dirname(HERE))

ARMS = ("SCAT32-060", "RIDGE-SCAT32", "TRUE-X-FW", "XHAT-FW", "RLMI", "ORACLE-LIB")


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _abspath(p):
    return p if os.path.isabs(p) else os.path.normpath(os.path.join(ROOT, p))


def _cell_done(cells_dir, scene, nu, c, seed):
    tag = "%s_%%s_nu%d_c%03d_s%d.json" % (scene, int(nu), int(round(c * 1000)), seed)
    return all(os.path.exists(os.path.join(cells_dir, tag % arm)) for arm in ARMS)


def verify_frozen(manifest):
    """sha256-verify every frozen_inputs entry (scenes + banks + code)."""
    bad = []
    for fi in manifest.get("frozen_inputs", []):
        p = _abspath(fi["path"])
        if not os.path.exists(p):
            bad.append((fi["path"], "MISSING"))
        elif _sha256(p) != fi["sha256"]:
            bad.append((fi["path"], "SHA_MISMATCH"))
    return bad


def run_shard(manifest_path, wall_budget_s=None, log=print, fw_iters=40):
    with open(manifest_path) as f:
        man = json.load(f)
    bad = verify_frozen(man)
    if bad:
        log("[FROZEN_FAIL] %s" % bad)
        return 3
    import bridge_harness as bh
    cells_dir = os.path.join(ROOT, "results", "round63_bridge", "cells")
    os.makedirs(cells_dir, exist_ok=True)
    cells = man["cells"]
    t0 = time.time()
    done = skipped = 0
    log("[shard %s] %d cells; layout=%s route=%s"
        % (man.get("shard_id"), len(cells), man.get("layout"), man.get("route")))
    for i, cell in enumerate(cells):
        sc, nu, c, seed = cell["scene"], cell["nu"], cell["c"], cell["seed"]
        if _cell_done(cells_dir, sc, nu, c, seed):
            skipped += 1
            continue
        if wall_budget_s is not None and (time.time() - t0) > wall_budget_s:
            log("[shard %s] wall budget reached; %d/%d done (resumable)"
                % (man.get("shard_id"), done + skipped, len(cells)))
            return 0
        tc = time.time()
        try:
            bh.run_group(sc, nu, c, seed, arms=ARMS, fw_iters=fw_iters, save=True)
            done += 1
            log("[hb] %s done=%d skip=%d cell=%s_nu%d_c%03d_s%d %.1fs elapsed=%.0fs"
                % (man.get("shard_id"), done, skipped, sc, int(nu),
                   int(round(c * 1000)), seed, time.time() - tc, time.time() - t0))
        except Exception as e:                        # noqa - keep the shard alive
            log("[CELL_ERROR] %s_nu%d_c%03d_s%d :: %r"
                % (sc, int(nu), int(round(c * 1000)), seed, e))
    log("[shard %s] COMPLETE done=%d skipped=%d total=%d wall=%.0fs"
        % (man.get("shard_id"), done, skipped, len(cells), time.time() - t0))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--wall-budget-s", type=float, default=None)
    ap.add_argument("--fw-iters", type=int, default=40)
    ap.add_argument("--logfile", default=None)
    args = ap.parse_args(argv)
    logf = open(args.logfile, "a") if args.logfile else None

    def log(m):
        line = "%s %s" % (time.strftime("%H:%M:%S"), m)
        print(line, flush=True)
        if logf:
            logf.write(line + "\n"); logf.flush()

    rc = run_shard(_abspath(args.manifest), wall_budget_s=args.wall_budget_s,
                   log=log, fw_iters=args.fw_iters)
    if logf:
        logf.close()
    return rc


if __name__ == "__main__":
    sys.exit(main())
