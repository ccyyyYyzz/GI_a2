# ROUND63 S2 — Colab launch runbook (operator sequence)

Run the 82 default S2/S3 shards across **five Colab sessions** (pro2 ×3 +
pro1 ×2), six lanes each, and merge their CSVs. This is the concrete command
sequence; the parallel/resume design lives in `README_COLAB.md`.

Real runtime values below are taken from the 2026-07-14 runtime handoff
(`F5_RUNTIME_ENVIRONMENT_HANDOFF_20260714.md §2`). A handful of account-local
transport values the handoff does **not** pin (Drive remote name, Drive dir,
keep-alive endpoint) are marked `<FILL_LIVE>` — pull them from the live account
/ the reusable launcher, do **not** invent them.

> Driver discipline (frozen program rules): **PowerShell → `wsl.exe` only** — no
> inline `python` in PowerShell; every `colab` CLI call is wrapped in
> `timeout 240`; large files move over Drive in `<=45 MB` chunks.

---

## 0. Constants (filled from the handoff)

```
COLAB_CLI   = /var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
AUTH        = --auth oauth2
TIMEOUT     = timeout 240              # wrap EVERY colab CLI call
HOME_pro2   = /var/tmp/codex-colab-accounts/pro2     # pro+ ; Drive RELAY account
HOME_pro1   = /var/tmp/codex-colab-accounts/pro1
GPU         = L4
REPO        = /mnt/d/GI_another        # == D:\GI_another
BUNDLE      = round63_s2_bundle.tar.gz
FROZEN ENV  = py3.11.15 numpy1.26.4 scipy1.13.1 torch2.2.1+cu121 tv0.17.1 lpips0.1.4
              4 BLAS thread-vars = 1 ; RSS<2GiB ; avail>1.5GiB ; LPIPS batch 16

# transport values the handoff does NOT pin — confirm from the live account /
# the reusable launcher .../conditional_fiberwise_pareto_controller/launch/reentry_v22/
DRIVE_REMOTE_pro2 = <FILL_LIVE>        # e.g. the rclone/gdrive remote for pro2
DRIVE_DIR         = <FILL_LIVE>        # e.g. GAN_FCC_WORK/round63/s2
KEEPALIVE_ENDPOINT= <FILL_LIVE>        # arg to `colab keep-alive ENDPOINT SESSION`
VM_DRIVE_MOUNT    = /content/drive/MyDrive     # confirm on the VM
```

### Session ↔ account ↔ plan map (frozen even→pro1 / odd→pro2 parity)

| session | account | HOME | plan file | shards | est_hours |
|---------|---------|------|-----------|-------:|----------:|
| `s2_pro1a` | pro1 | `HOME_pro1` | `plans/session_pro1a.txt` | 21 | 114.16 |
| `s2_pro1b` | pro1 | `HOME_pro1` | `plans/session_pro1b.txt` | 21 | 116.97 |
| `s2_pro2a` | pro2 | `HOME_pro2` | `plans/session_pro2a.txt` | 14 | 71.38 |
| `s2_pro2b` | pro2 | `HOME_pro2` | `plans/session_pro2b.txt` | 13 | 70.06 |
| `s2_pro2c` | pro2 | `HOME_pro2` | `plans/session_pro2c.txt` | 13 | 71.40 |

(82 shards; regenerate the table any time with `python lane_plan.py`.)

---

## 1. Build the bundle (local, WSL)

```powershell
# PowerShell:
wsl.exe -e bash -lc 'cd /mnt/d/GI_another/code/round63/colab && python3 lane_plan.py && bash make_bundle.sh'
```

`lane_plan.py` (re)writes the 5 `plans/session_*.txt`; `make_bundle.sh` prints
the bundle **sha256 + size** and embeds the plans + drivers. Record the sha —
it is the transport integrity gate for every session. The bundle unpacks into a
single `round63_s2_bundle/` directory laid out exactly like the repo root, so
`shard_runner.py` resolves `frozen_inputs`, `output_csv`, the data caches and
`C0_FROZEN.json` unchanged.

If the bundle exceeds 50 MB (it is ~0.25 MB today, so this will not fire) the
script emits `<=45 MB` chunks and prints the reassembly line.

---

## 2. Upload the bundle to Drive (pro2 relay, chunk protocol)

pro2 is the **only** Drive relay account. Stage the canonical copy there; if a
pro1 VM cannot read the pro2 relay, also push the identical (sha-matched) bundle
to pro1's own Drive dir.

```powershell
# PowerShell -> WSL ; chunked, resume-aware (weak local net: 5x retry / 12 s)
wsl.exe -e bash -lc '
  set -e
  cd /mnt/d/tmp/round63_s2_build
  SZ=$(stat -c%s round63_s2_bundle.tar.gz)
  if [ "$SZ" -gt 52428800 ]; then split -b 45M -d -a3 round63_s2_bundle.tar.gz round63_s2_bundle.tar.gz.part_; fi
  # upload via the reentry_v22 chunked transport (remote hash gate over every byte):
  HOME=/var/tmp/codex-colab-accounts/pro2 timeout 240 \
    <FILL_LIVE: colab/rclone upload> round63_s2_bundle.tar.gz* <DRIVE_REMOTE_pro2>:<DRIVE_DIR>/
'
```

Confirm the remote sha equals the local sha before launching anything.

---

## 3. Create the 5 sessions, prove liveness, attach keep-alive

Repeat per row of the map (§0). Example for `s2_pro2a`:

```powershell
# (a) create the session
wsl.exe -e bash -lc 'HOME=/var/tmp/codex-colab-accounts/pro2 timeout 240 /var/tmp/codex-colab-tools/colab-cli-venv/bin/colab new --session s2_pro2a --gpu L4 --auth oauth2'

# (b) PROVE liveness by a nonce exec — never trust a cached sessions row
wsl.exe -e bash -lc 'N=$RANDOM; HOME=/var/tmp/codex-colab-accounts/pro2 timeout 240 /var/tmp/codex-colab-tools/colab-cli-venv/bin/colab exec s2_pro2a "echo NONCE_$N" | grep -q "NONCE_$N" && echo LIVE || echo DEAD'

# (c) attach keep-alive IMMEDIATELY, detached (setsid) — an unwatched VM
#     idle-recycles in ~1-2 h and wipes /content
wsl.exe -e bash -lc 'HOME=/var/tmp/codex-colab-accounts/pro2 setsid timeout 240 /var/tmp/codex-colab-tools/colab-cli-venv/bin/colab keep-alive <KEEPALIVE_ENDPOINT> s2_pro2a &'
```

pro1 sessions are identical with `HOME=/var/tmp/codex-colab-accounts/pro1` and
session names `s2_pro1a` / `s2_pro1b`. (Confirm the exact `exec` / `keep-alive`
subcommand names with `colab --help` and the reentry_v22 launcher if the CLI
differs.)

---

## 4. On each VM: mount Drive, stage + verify bundle, launch the driver

Run once per session (as an `exec` payload or first cell). Keep-alive is already
attached (step 3c); `session_driver.sh` adds its own touch loop as a second
guard. The driver itself is launched **detached via `setsid`** so the exec call
returns while compute continues.

```bash
# on the VM (payload for `colab exec <session> 'bash -lc "<this>"'`)
set -euo pipefail
SESSION=s2_pro2a                       # <- per session
PLAN=session_pro2a.txt                 # <- per session (see map §0)
DRIVE=/content/drive/MyDrive/<DRIVE_DIR>
cd /content

# (a) reassemble (if chunked) + verify sha, then unpack
cat "$DRIVE"/round63_s2_bundle.tar.gz.part_* > round63_s2_bundle.tar.gz 2>/dev/null \
  || cp "$DRIVE"/round63_s2_bundle.tar.gz .
echo "<BUNDLE_SHA256>  round63_s2_bundle.tar.gz" | sha256sum -c -
tar -xzf round63_s2_bundle.tar.gz
BDIR=/content/round63_s2_bundle

# (b) self-check every bundled file, then (optional resume) restore prior CSV+meta
( cd "$BDIR" && sha256sum -c MANIFEST_SHA256 --quiet )
mkdir -p "$BDIR/results/round63/shards"
cp -f "$DRIVE"/out/"$SESSION"/shards/* "$BDIR/results/round63/shards/" 2>/dev/null || true

# (c) checkpoint loop: mirror results back to Drive every 2 min (recycle-safe)
( while true; do
    mkdir -p "$DRIVE/out/$SESSION"
    cp -rf "$BDIR/results/round63/shards" "$DRIVE/out/$SESSION/" 2>/dev/null || true
    cp -f  "$BDIR/results/round63/status/"*.json "$DRIVE/out/$SESSION/" 2>/dev/null || true
    sleep 120
  done ) & echo "checkpoint pid $!"

# (d) launch the session driver, detached; 6 lanes; BLAS pinned inside each lane
setsid bash "$BDIR/code/round63/colab/session_driver.sh" \
    "$BDIR/code/round63/colab/plans/$PLAN" \
    --bundle-root "$BDIR" --session "$SESSION" --n-lanes 6 \
    --status-json "$BDIR/results/round63/status/session_$SESSION.json" \
    --python python3 --wall-budget-s 21600 \
    > "$BDIR/results/round63/driver_$SESSION.out" 2>&1 &
echo "BACKGROUND_LAUNCHED driver_$SESSION pid $!"
```

`--n-lanes 6` assumes ≥6 vCPUs on the VM (each lane pins BLAS to 1 thread).
Drop it to the VM's vCPU count if the L4 runtime is CPU-thin; six side-64 lanes
fit comfortably in L4 RAM (the 128-px S2C128 shards are **not** in this pack).

---

## 5. Monitor (local): pull the 5 status JSONs from Drive

Each driver rewrites `results/round63/status/session_<name>.json` every 2 min;
the checkpoint loop (4c) mirrors it to `Drive/out/<session>/`. Pull + read:

```powershell
wsl.exe -e bash -lc '
  for S in s2_pro1a s2_pro1b s2_pro2a s2_pro2b s2_pro2c; do
    ACCT=pro2; case "$S" in *pro1*) ACCT=pro1;; esac
    HOME=/var/tmp/codex-colab-accounts/$ACCT timeout 240 \
      <FILL_LIVE: colab/rclone download> <DRIVE_REMOTE>:<DRIVE_DIR>/out/$S/session_$S.json /mnt/d/tmp/status/ 2>/dev/null || true
  done
  python3 - <<PY
import glob, json, time
now=int(time.time())
for p in sorted(glob.glob("/mnt/d/tmp/status/session_*.json")):
    d=json.load(open(p))
    act=[l for l in d["lanes"] if l["state"]=="running"]
    print("%-9s done %d/%d  running %d  status_age %ds"%(
        d["session"], d["completed_count"], d["total_shards"],
        len(act), now-d["ts"]))
    for l in act:
        print("     lane %d %-16s hb_age %ss cells %s rows %s"%(
            l["lane"], l["shard"], l.get("heartbeat_age_s"),
            l.get("cells_done"), l.get("rows")))
PY
'
```

Reading it: `cells_done` can sit still for ~40–50 min (one cell) and that is
**normal**; liveness is `heartbeat_age_s` (a live lane ticks every 60 s
regardless of cell progress). `status_age` is how stale the whole session file
is (a healthy driver rewrites every 120 s).

---

## 6. Retrieve + merge results (local, WSL)

```powershell
wsl.exe -e bash -lc '
  set -e
  DST=/mnt/d/GI_another/results/round63/shards
  mkdir -p "$DST"
  for S in s2_pro1a s2_pro1b s2_pro2a s2_pro2b s2_pro2c; do
    ACCT=pro2; case "$S" in *pro1*) ACCT=pro1;; esac
    HOME=/var/tmp/codex-colab-accounts/$ACCT timeout 240 \
      <FILL_LIVE: download> <DRIVE_REMOTE>:<DRIVE_DIR>/out/$S/shards/ "$DST"/
  done
  # each shard writes an independent CSV + *_meta.json; concat CSVs per stage
  # (identical headers within a stage). Example for S2A_DETAIL:
  head -1 "$DST"/S2A_DETAIL_00.csv > "$DST"/../S2A_DETAIL_merged.csv
  for f in "$DST"/S2A_DETAIL_*.csv; do tail -n +2 "$f"; done >> "$DST"/../S2A_DETAIL_merged.csv
  wc -l "$DST"/../S2A_DETAIL_merged.csv
'
```

Expected total rows across all default shards: **45396** (per
`MANIFEST_INDEX.json:n_expected_rows`). A short shard count means resume is still
pending — see §7.

---

## 7. VM recycle → re-upload + relaunch (shards resume)

An idle recycle wipes `/content`. Because `shard_runner.py` is META-as-truth and
the checkpoint loop (4c) mirrored each shard's `*.csv` + `*_meta.json` to Drive:

1. Recreate the session (§3) with the **same** name; prove liveness; keep-alive.
2. Re-run the §4 payload with the **same** `PLAN`. Step 4b restores the last
   checkpointed `shards/` into the bundle; step 4d relaunches the driver.
3. Each shard's runner sanitizes any orphan CSV rows from the interrupted cell,
   skips completed cells, and continues. `--wall-budget-s` is the **accumulated**
   budget across resumes — raise it if a shard needs several sessions.

Nothing is recomputed that already landed in a checkpointed CSV.

---

## 8. WATCHDOG (local) — poll every 15 min, print stalls

Save as `s2_watchdog.sh` and run via `wsl.exe -e bash -lc 'bash s2_watchdog.sh'`
(or from PowerShell in a loop). It pulls the 5 status JSONs and flags: a **lane
stall** (live lane whose `heartbeat_age_s` exceeds `HB_STALL`, default 300 s) and
a **dead session** (`status_age` exceeds `SESS_STALL`, default 600 s → the driver
stopped writing, i.e. the VM likely recycled → go to §7).

```bash
#!/usr/bin/env bash
set -uo pipefail
HB_STALL=${HB_STALL:-300}          # a live lane heartbeats every 60 s
SESS_STALL=${SESS_STALL:-600}      # a live driver rewrites status every 120 s
STATUS=/mnt/d/tmp/status
SESSIONS="s2_pro1a s2_pro1b s2_pro2a s2_pro2b s2_pro2c"
while true; do
  mkdir -p "$STATUS"
  for S in $SESSIONS; do
    ACCT=pro2; case "$S" in *pro1*) ACCT=pro1;; esac
    HOME=/var/tmp/codex-colab-accounts/$ACCT timeout 240 \
      <FILL_LIVE: download> <DRIVE_REMOTE>:<DRIVE_DIR>/out/$S/session_$S.json "$STATUS"/ 2>/dev/null || true
  done
  echo "===== watchdog $(date -u +%H:%M:%SZ) ====="
  HB_STALL=$HB_STALL SESS_STALL=$SESS_STALL python3 - "$STATUS" <<'PY'
import glob, json, os, sys, time
status_dir=sys.argv[1]; now=int(time.time())
hb_stall=int(os.environ["HB_STALL"]); sess_stall=int(os.environ["SESS_STALL"])
seen=set()
for p in sorted(glob.glob(os.path.join(status_dir,"session_*.json"))):
    try: d=json.load(open(p))
    except Exception: continue
    S=d["session"]; seen.add(S); age=now-d["ts"]
    tag="OK"
    if age>sess_stall: tag="DEAD-SESSION(status_age=%ds) -> recycle? see runbook §7"%age
    print("%-9s %-40s done %d/%d"%(S, tag, d["completed_count"], d["total_shards"]))
    for l in d["lanes"]:
        if l["state"]=="running":
            ha=l.get("heartbeat_age_s",-1)
            if ha is not None and ha>hb_stall:
                print("   STALL lane %d %s heartbeat_age=%ss (>%ds)"%(
                    l["lane"], l["shard"], ha, hb_stall))
for S in "s2_pro1a s2_pro1b s2_pro2a s2_pro2b s2_pro2c".split():
    if S not in seen: print("%-9s NO-STATUS-FILE (never started / not mirrored)"%S)
PY
  sleep 900
done
```

---

## Exit codes (per-shard, surfaced by `remote_lane` → `session_driver` completed[].rc)

| rc | meaning | action |
|----|---------|--------|
| 0 | all cells done | shard complete |
| 2 | wall-budget abort | resumable — relaunch (§7) with larger `--wall-budget-s` |
| 3 | `campaign.run_cell` unavailable | code bug — should not happen (campaign is bundled) |
| 4 | frozen-input sha256 mismatch | bundle corruption — rebuild + re-upload (§1–2) |

A non-zero session exit means ≥1 shard did not reach rc 0.

## Placeholder checklist (fill before first launch)

- [ ] `DRIVE_REMOTE_pro2` (+ pro1 Drive dir if pro1 VMs can't read the relay)
- [ ] `DRIVE_DIR`, `VM_DRIVE_MOUNT`
- [ ] `KEEPALIVE_ENDPOINT`, and confirm `colab exec` / `keep-alive` subcommand names
- [ ] the chunked upload/download command (reentry_v22 transport)
- [ ] `<BUNDLE_SHA256>` from step 1 pasted into the step-4 verify line
- [ ] per-session `--wall-budget-s` (default 21600 s = 6 h; raise for multi-session shards)
