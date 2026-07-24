# Sealed-Probe Fleet Runbook (R39, 5-shard Colab lane) — LAUNCH ONLY AFTER THE R39 RULING

Build-only. Do **not** run until `bars64.py::THRESHOLDS` are frozen (the `[R39]` placeholders → the
referee's numbers). Every step below complies with `docs/COLAB_USAGE_GUIDE.md` to the letter.

## Preconditions
- `bank_gen64.py` has been run (sealed banks under `scene_banks/`, determinism+all-distinct PASS).
- `sealed_probe_planner.py` has been run (5 shard manifests + `LANE_PLAN.json` under `manifests/`).
- Thresholds frozen in `bars64.py` (one edit).

## Shard↔account parity (FROZEN for the official run)
`{shard0:pro1, shard1:pro2, shard2:pro1, shard3:pro2, shard4:pro2}` — even→pro1 (2×L4),
odd→pro2 (3×L4). pro2 carries the extra shard (3 slots).

## Colab discipline (mandatory, every call)
1. Every CLI call: `HOME=<account-dir> timeout 240 $C --auth oauth2 …`; long execs get their own
   `--timeout <s>`. `C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab`.
2. **All** wsl-side commands are FILE-CARRIED `.sh` written to `D:\`, then `sed -i 's/\r$//'`, then
   `wsl.exe -- bash <script>`. Never inline `$`-vars through PowerShell.
3. On 401 / vanished session rows: run `code/round63/colab/live_rebind.sh` — **never recreate
   sessions** (the VM is still computing).
4. Keep-alives are `setsid`-detached (or `Start-Process wsl.exe` from Windows), re-attached every
   watchdog cycle.
5. VM-side driver is autonomous: task list on the VM (shard manifest), heartbeat JSON, per-cell
   checkpoints (`<shard>_meta.json` + CSV) — token/network loss never kills computation.
6. Watchdog cycle order: **rebind → heartbeat check → idempotent fetch** (`live_watch.sh` /
   `live_fetch_all.sh`).
7. VM stack numpy 1.26.4 / scipy 1.13.1 / torch preinstalled — **never upgrade packages**.
8. On completion: RELEASE every VM (`state.client.unassign`, `live_release_*.sh`), kill all
   keep-alives, report `assignments=NONE`.

## Launch sequence (maps to the proven templates)
1. **Bundle** — `make_bundle.sh`-style: pack `results/round63_next/FOG_DMD_PROBE64/` (arms64,
   fog_tracker, sealed_shard_runner, scene_banks/, manifests/) + `../FOG_DMD_PROBE/fog_tracker.py`.
2. **Create sessions** — `live_create_sessions.sh` (5 sessions, GPU L4; pro1×2, pro2×3).
3. **Launch all** — `live_launch_all.sh`: upload bundle, unpack, and per lane run
   `remote_lane.py --manifest manifests/r39bb_shard<i>.json --bundle-root <root> --python <py>`
   with `--python` pointing at the VM interpreter and `remote_lane` invoking
   `sealed_shard_runner.py` (swap the runner path in the lane launcher — the R39 runner is a drop-in
   for `shard_runner.py`: same manifest/heartbeat/CSV/meta contract).
4. **Keep-alives** — `live_keepalive_attach.sh` (setsid), one per session.
5. **Watch** — `live_watch.sh` (rebind → heartbeat → fetch), hourly summary line.
6. **Fetch** — `live_fetch_all.sh`: idempotent pull of `shards/r39bb_shard<i>.csv` + `_meta.json`.
7. **Merge + adjudicate** — merge the 5 shard CSVs into `confirmatory/CONFIRMATORY_RESULTS.json`
   (shape per `bars64._synthetic`), then `python bars64.py` → `B_BARS_VERDICT.json`.
8. **Release** — `live_release_pro1.sh` + `live_release_pro2.sh`; kill keep-alives; confirm
   `assignments=NONE`.

## Local-vs-Colab decision (this session)
The dev-grade gates G1/G2/G3 ran **locally** on the free RTX 4060 (jobs ≪ 1 h), which is the
coordinator's own stated condition for local execution ("if contended OR >~1h → Colab"). No live
Colab calls were made; items (1)–(8) above bind the *confirmatory* launch, not the dev gates.
