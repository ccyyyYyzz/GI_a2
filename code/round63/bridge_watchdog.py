"""DEV-bridge grid watchdog (R28 launch): 30-min digests + stall detection.

Runs detached alongside the 5 shard processes.  Every DIGEST_S seconds it writes
one digest line to results/round63_bridge/logs/WATCHDOG_DIGEST.log with per-shard
cells-done, total progress, throughput, ETA, and a STALL flag for any shard whose
log has had no [hb] heartbeat for > STALL_S while its process is presumed alive.
Exits when all 320 cells are done or all shard logs have gone quiet past STALL_S.
"""
import glob
import os
import re
import time

ROOT = r"D:/GI_another/results/round63_bridge"
LOGS = os.path.join(ROOT, "logs")
CELLS = os.path.join(ROOT, "cells")
DIGEST = os.path.join(LOGS, "WATCHDOG_DIGEST.log")
ARMS = ("SCAT32-060", "RIDGE-SCAT32", "TRUE-X-FW", "XHAT-FW", "RLMI", "ORACLE-LIB")
ROUTES = ("pro2_a", "pro2_b", "pro2_c", "pro1_a", "pro1_b")
TOTAL_CELLS = 320
DIGEST_S = 1800          # 30 min
STALL_S = 2400           # 40 min with no heartbeat -> flag


def cells_done():
    files = set(os.path.basename(p) for p in glob.glob(os.path.join(CELLS, "*.json")))
    stems = set()
    for f in files:
        m = re.match(r"(bridge_.+?)_(%s)_(nu\d+_c\d+_s\d+)\.json$"
                     % "|".join(re.escape(a) for a in ARMS), f)
        if m:
            stems.add((m.group(1), m.group(3)))
    n = 0
    for scene, tail in stems:
        if all("%s_%s_%s.json" % (scene, a, tail) in files for a in ARMS):
            n += 1
    return n


def shard_status():
    out = {}
    now = time.time()
    for r in ROUTES:
        lp = os.path.join(LOGS, "%s.log" % r)
        if not os.path.exists(lp):
            out[r] = ("no-log", None, 0)
            continue
        age = now - os.path.getmtime(lp)
        try:
            with open(lp) as f:
                lines = f.readlines()
        except Exception:
            lines = []
        hb = sum(1 for ln in lines if "[hb]" in ln)
        complete = any("COMPLETE" in ln for ln in lines[-5:])
        state = "COMPLETE" if complete else ("STALL" if age > STALL_S else "running")
        out[r] = (state, int(age), hb)
    return out


def main():
    open(DIGEST, "a").write("%s [watchdog START] 5 shards, %d cells\n"
                            % (time.strftime("%H:%M:%S"), TOTAL_CELLS))
    start = time.time()
    last = 0
    while True:
        time.sleep(min(DIGEST_S, 300))
        n = cells_done()
        ss = shard_status()
        el = time.time() - start
        rate = n / max(el / 3600.0, 1e-6)             # cells/hr
        remain = TOTAL_CELLS - n
        eta_h = remain / rate if rate > 0 else float("inf")
        # emit a digest at DIGEST_S cadence (or on completion)
        if time.time() - last >= DIGEST_S or n >= TOTAL_CELLS:
            last = time.time()
            per = " ".join("%s=%s(hb%d,%ds)" % (r, ss[r][0], ss[r][2], ss[r][1] or 0)
                           for r in ROUTES)
            line = ("%s [digest] done=%d/%d (%.0f%%) rate=%.1f/hr ETA=%.1fh | %s"
                    % (time.strftime("%H:%M:%S"), n, TOTAL_CELLS,
                       100.0 * n / TOTAL_CELLS, rate, eta_h, per))
            with open(DIGEST, "a") as f:
                f.write(line + "\n")
        if n >= TOTAL_CELLS:
            with open(DIGEST, "a") as f:
                f.write("%s [watchdog DONE] all %d cells materialized; wall=%.1fh\n"
                        % (time.strftime("%H:%M:%S"), n, el / 3600.0))
            break
        # exit if every shard is COMPLETE or STALL (nothing left running)
        if all(ss[r][0] in ("COMPLETE", "STALL", "no-log") for r in ROUTES):
            with open(DIGEST, "a") as f:
                f.write("%s [watchdog EXIT] no running shard; done=%d/%d %s\n"
                        % (time.strftime("%H:%M:%S"), n, TOTAL_CELLS,
                           {r: ss[r][0] for r in ROUTES}))
            break


if __name__ == "__main__":
    main()
