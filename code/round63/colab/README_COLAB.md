# ROUND63 — launching ONE campaign shard on Colab (from WSL)

This is the operator runbook for running a single `shard_runner.py` shard on one
Colab session and getting its `results/round63/*.csv` back. It follows the
program's standing runtime conventions:

- **Driver is WSL, not PowerShell.** From Windows, only ever do
  `powershell -> wsl.exe -- <bash>`; never run inline `python` in PowerShell.
- **Transport is Google Drive.** Bundles up, result CSV + meta back. Anything
  `>50 MB` is uploaded in `<=45 MB` chunks and reassembled on the far side.
- **BLAS threads pinned to 1 BEFORE numpy is imported** (`OMP/MKL/OPENBLAS/
  NUMEXPR = 1`) — set in the launch script's environment, not in Python.
- **Keep-alive is attached inside the launch script itself**, plus a guardian
  watchdog, because an idle Colab VM recycles and wipes `/content`. The runner is
  resumable (META-as-truth), so a recycle costs at most the in-flight cell.

Anything that is machine-/account-specific (the exact `colab` CLI binary, its
launch flags, the Drive remote name, the shard↔account map) is marked
`<FILL_FROM_RUNTIME_HANDOFF>` — pull the real value from the current runtime
handoff (see the "Local + Colab runtime envs" and "Colab account policy" notes)
rather than inventing it.

---

## 0. Pick the account/session for this shard

Formal shard↔account parity is frozen: **even shard index -> pro1, odd -> pro2**
(pro2 is HOME). Confirm the live mapping and the session slot before launching.

```
ACCOUNT            = <FILL_FROM_RUNTIME_HANDOFF>   # e.g. pro2 (odd shards) / pro1 (even)
COLAB_CLI          = <FILL_FROM_RUNTIME_HANDOFF>   # path to the colab launcher binary
COLAB_LAUNCH_FLAGS = <FILL_FROM_RUNTIME_HANDOFF>   # runtime/GPU/session flags
DRIVE_REMOTE       = <FILL_FROM_RUNTIME_HANDOFF>   # rclone/gdrive remote for this account
DRIVE_DIR          = <FILL_FROM_RUNTIME_HANDOFF>   # e.g. GAN_FCC_WORK/round63/shards
```

---

## 1. Assemble the shard bundle (local, in WSL)

A shard bundle is fully self-contained: the code it imports, the frozen inputs
its cells consume, and the manifest. Build it from the repo root
(`D:\GI_another` == `/mnt/d/GI_another` in WSL).

```bash
# from Windows:  powershell -> wsl.exe -- bash -lc '<the block below>'
set -euo pipefail
REPO=/mnt/d/GI_another
SHARD=<FILL_FROM_RUNTIME_HANDOFF>          # e.g. S2_rho0.05_a  (matches manifest shard_id)
STAGE=<manifest-derived>
WORK=/mnt/d/tmp/round63_bundles/$SHARD
rm -rf "$WORK"; mkdir -p "$WORK/code" "$WORK/data"

# code: the round63 tools + the gi_core package they import (flat-import layout)
cp -r "$REPO/code/round63"  "$WORK/code/round63"
cp -r "$REPO/code/gi_core"  "$WORK/code/gi_core"
# the manifest for THIS shard (author it per the schema in shard_runner.py's docstring)
cp "$REPO/results/round63/manifests/$SHARD.json" "$WORK/$SHARD.json"

# frozen inputs: copy exactly the files the manifest lists (images, patterns, etc.)
#   -> keep the same repo-relative paths so manifest paths resolve unchanged
#   e.g. cp --parents data/r63_images/64/stl_00.png ...   (run from $REPO)
( cd "$REPO" && python code/round63/_collect_frozen.py "$WORK/$SHARD.json" "$WORK" ) \
  2>/dev/null || echo "NOTE: copy frozen inputs listed in the manifest into $WORK/ (same rel paths)"

# record a bundle-level manifest of sha256 so the far side can self-check
( cd "$WORK" && find . -type f -print0 | xargs -0 sha256sum > BUNDLE_SHA256.txt )
tar -C "$WORK/.." -czf "$WORK.tar.gz" "$(basename "$WORK")"
ls -l "$WORK.tar.gz"
```

The manifest's `frozen_inputs[*].sha256` are the authoritative gate — the runner
re-verifies them on the VM before computing anything.

---

## 2. Upload the bundle to Drive (chunked if >50 MB)

```bash
BUNDLE=$WORK.tar.gz
SZ=$(stat -c%s "$BUNDLE")
if [ "$SZ" -gt 52428800 ]; then           # >50 MB -> split into <=45 MB parts
  split -b 45M -d "$BUNDLE" "$BUNDLE.part_"
  for p in "$BUNDLE".part_*; do
    <FILL_FROM_RUNTIME_HANDOFF: DRIVE upload cmd> "$p" "$DRIVE_REMOTE:$DRIVE_DIR/"
  done
  # reassemble on the VM in step 4:  cat <bundle>.part_* > <bundle>.tar.gz
else
  <FILL_FROM_RUNTIME_HANDOFF: DRIVE upload cmd> "$BUNDLE" "$DRIVE_REMOTE:$DRIVE_DIR/"
fi
```

---

## 3. Launch the Colab session (keep-alive attached in the launch script)

Start the session with the account/flags from step 0. The launch payload is the
`run_shard.sh` in step 4 — pass it as the session's entry command so keep-alive
and the runner come up together.

```bash
<COLAB_CLI> <COLAB_LAUNCH_FLAGS> \
    --account "$ACCOUNT" \
    --run "bash /content/run_shard.sh $SHARD"     # <FILL_FROM_RUNTIME_HANDOFF: exact flag names>
```

If the launcher cannot take an entry command, open the session and paste the
step-4 script into the first cell. Keep-alive MUST be part of that same script
(do not rely on a separate manual keep-alive).

---

## 4. On the VM: pin BLAS, mount Drive, verify, run with keep-alive + watchdog

Save this as `run_shard.sh` (upload alongside the bundle, or paste as cell 1).
Note the BLAS env is exported **before** any Python/numpy touches the process.

```bash
#!/usr/bin/env bash
set -euo pipefail
SHARD="${1:?shard id}"

# (a) BLAS threads pinned to 1 BEFORE numpy import — export in the env, not in py
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

# (b) keep-alive: attached here, not manual. Survives idle detection.
( while true; do echo "[keepalive $(date -u +%H:%M:%S)]"; sleep 60; done ) &
KEEPALIVE_PID=$!
trap 'kill $KEEPALIVE_PID 2>/dev/null || true' EXIT

# (c) mount Drive and stage the bundle into /content (wiped on recycle)
DRIVE=<FILL_FROM_RUNTIME_HANDOFF: Drive mount path on VM, e.g. /content/drive/MyDrive>
DDIR="$DRIVE/<FILL_FROM_RUNTIME_HANDOFF: DRIVE_DIR>"
cd /content
cat "$DDIR/$SHARD.tar.gz".part_* > "$SHARD.tar.gz" 2>/dev/null || \
    cp "$DDIR/$SHARD.tar.gz" "$SHARD.tar.gz"
tar -xzf "$SHARD.tar.gz"
BDIR="/content/$SHARD"

# (d) self-check the bundle sha256, then let the runner re-verify frozen inputs
( cd "$BDIR" && sha256sum -c BUNDLE_SHA256.txt --quiet )

# (e) guardian watchdog: checkpoint results back to Drive every 2 min so an idle
#     recycle (which wipes /content) never loses more than the in-flight cell.
RESULTS="$BDIR/results/round63"
( while true; do
    mkdir -p "$DDIR/out/$SHARD"
    [ -d "$RESULTS" ] && cp -f "$RESULTS"/* "$DDIR/out/$SHARD/" 2>/dev/null || true
    sleep 120
  done ) &
WATCHDOG_PID=$!
trap 'kill $KEEPALIVE_PID $WATCHDOG_PID 2>/dev/null || true' EXIT

# (f) run the shard. cwd = bundle/code so the flat gi_core / round63 imports resolve.
cd "$BDIR/code"
python round63/shard_runner.py \
    --manifest "$BDIR/$SHARD.json" \
    --wall-budget-s <FILL_FROM_RUNTIME_HANDOFF: session budget, e.g. 39600>
RC=$?

# (g) final push of CSV + meta to Drive
mkdir -p "$DDIR/out/$SHARD"; cp -f "$RESULTS"/* "$DDIR/out/$SHARD/" 2>/dev/null || true
echo "[run_shard] shard $SHARD exit code $RC"
exit $RC
```

Runner exit codes: `0` all cells done · `2` budget abort (resume later) ·
`3` `campaign.py` not yet present · `4` frozen-input sha mismatch.

---

## 5. Retrieve + merge results (local, WSL)

```bash
<FILL_FROM_RUNTIME_HANDOFF: DRIVE download cmd> \
    "$DRIVE_REMOTE:$DRIVE_DIR/out/$SHARD/" /mnt/d/GI_another/results/round63/
# each shard writes an independent CSV + *_meta.json; merge CSVs by concatenation
# (headers identical within a stage) once all shards for a stage are back.
```

---

## 6. Resume after an idle recycle

The runner is idempotent. If the VM recycled mid-shard:

1. Re-launch the session (step 3) with the **same** `SHARD`.
2. `run_shard.sh` re-stages the bundle AND pulls back the last-checkpointed
   `results/round63/<shard>.csv` + `_meta.json` from `$DDIR/out/$SHARD`
   into the bundle's `results/round63/` before invoking the runner
   (add that `cp` in step 4(c) when resuming).
3. `shard_runner.py` reads the meta, sanitizes any orphan CSV rows from the
   interrupted cell, skips completed cells, and continues. `--wall-budget-s` is
   the **accumulated** budget across resumes, so raise it if a shard needs
   several sessions.

---

## Placeholder checklist (fill before first launch)

- [ ] `COLAB_CLI`, `COLAB_LAUNCH_FLAGS`, entry-command flag name
- [ ] `ACCOUNT` + shard→account parity for this shard index
- [ ] `DRIVE_REMOTE`, `DRIVE_DIR`, VM Drive mount path, up/download commands
- [ ] per-session `--wall-budget-s`
- [ ] frozen-input collection step (or copy the manifest's files by hand)
