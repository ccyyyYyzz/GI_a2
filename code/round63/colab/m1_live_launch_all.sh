#!/bin/bash
# ROUND63 **M1** live launch -- per session: nonce liveness, bundle upload,
# env-pin bootstrap, per-session lane-plan upload, detached driver launch,
# keep-alive. Modeled on live_launch_all.sh (S2) + session_driver.sh, retargeted
# at the M1 manifests/bundle. File-carried per the runtime discipline; every CLI
# call wrapped in timeout.
#
# ==========================================================================
#  SAFETY: this launches the M1 CONFIRMATORY campaign. It MUST NOT run before
#  the `m1-freeze` git tag exists (no confirmatory cell may run pre-freeze).
#  A hard guard below refuses to launch until the tag is present.
#    --plan-only         generate + print the 5 lane plans and exit (no VM,
#                        no upload, no compute) -- safe to run anytime.
#    --force-no-freeze   override the freeze guard (NOT recommended).
# ==========================================================================
#
# Fleet (see results/round63_m1/COLAB_FLEET.json): 5 x L4, 6 lanes each.
#   pro2: r63m1_a  r63m1_b  r63m1_c        pro1: r63m1_p1a  r63m1_p1b
# Shards: all 40 from results/round63_m1/manifests/MANIFEST_INDEX.json,
# assigned by greedy-LPT to balance total est_hours across ALL FIVE sessions
# (NOT the S2 even/odd parity -- the task calls for a 5-way cost balance).
set -u

REPO=/mnt/d/GI_another
STAGE=/var/tmp/r63_m1_stage
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
BUNDLE_ROOT_VM=/content/round63_m1_bundle
BUNDLE_TGZ_VM=/content/round63_m1_bundle.tar.gz
BUNDLE="$STAGE/round63_m1_bundle.tar.gz"
INDEX="$REPO/results/round63_m1/manifests/MANIFEST_INDEX.json"
PLAN_DIR="$STAGE/plans"
N_LANES=6
WALL_BUDGET_S=43200

PLAN_ONLY=0
FORCE_NO_FREEZE=0
while [ $# -gt 0 ]; do
  case "$1" in
    --plan-only)       PLAN_ONLY=1; shift ;;
    --force-no-freeze) FORCE_NO_FREEZE=1; shift ;;
    -h|--help) sed -n '2,25p' "$0"; exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done

mkdir -p "$STAGE" "$PLAN_DIR"

# ---- (0) lane plan: greedy-LPT balance of est_hours across all five sessions - #
#   Reads MANIFEST_INDEX.json (default_shards + blocked_shards, dedup by id).
#   Emits session_m1_<tag>.txt listing one REPO-RELATIVE manifest path per line
#   (remote_lane resolves it as bundle_root/<path>; a bare shard-id would wrongly
#   resolve under results/round63/manifests, not results/round63_m1/manifests).
python3 - "$INDEX" "$PLAN_DIR" <<'PYEOF'
import json, os, sys
index_path, out_dir = sys.argv[1], sys.argv[2]
idx = json.load(open(index_path))
shards, seen = [], set()
for key in ("default_shards", "blocked_shards"):
    for s in idx.get(key, []):
        sid = s["shard_id"]
        if sid in seen:
            continue
        seen.add(sid)
        shards.append({"id": sid, "est": float(s.get("est_hours", 0.0)),
                       "manifest": "results/round63_m1/manifests/%s.json" % sid})
SESSIONS = ["a", "b", "c", "p1a", "p1b"]      # 3 pro2 + 2 pro1
ACCT = {"a": "pro2", "b": "pro2", "c": "pro2", "p1a": "pro1", "p1b": "pro1"}
plan = {s: [] for s in SESSIONS}
load = {s: 0.0 for s in SESSIONS}
for sh in sorted(shards, key=lambda c: (-c["est"], c["id"])):   # LPT: descending
    tgt = min(SESSIONS, key=lambda s: (load[s], SESSIONS.index(s)))
    plan[tgt].append(sh); load[tgt] += sh["est"]
os.makedirs(out_dir, exist_ok=True)
for s in SESSIONS:
    p = os.path.join(out_dir, "session_m1_%s.txt" % s)
    lines = ["# M1 lane plan -- session m1_%s (account %s)" % (s, ACCT[s]),
             "# %d shards, %.3f est_hours (greedy-LPT balanced across 5 sessions)"
             % (len(plan[s]), load[s]),
             "# one repo-relative manifest path per line; '#'/blank ignored."]
    for sh in sorted(plan[s], key=lambda c: c["id"]):
        lines.append(sh["manifest"])
    with open(p, "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
print("M1 lane plan (%d shards, balance total est_hours across 5 sessions):"
      % len(shards))
tot = 0.0; ns = 0
for s in SESSIONS:
    print("  m1_%-4s (%s): %2d shards  %6.3f est_hours"
          % (s, ACCT[s], len(plan[s]), load[s]))
    tot += load[s]; ns += len(plan[s])
print("  TOTAL: %d shards  %.3f est_hours  |  spread(max-min) %.3f"
      % (ns, tot, max(load.values()) - min(load.values())))
if ns != len(shards):
    raise SystemExit("BUG: assigned %d != %d shards" % (ns, len(shards)))
PYEOF
[ $? -eq 0 ] || { echo "PLAN_GEN_FAILED" >&2; exit 1; }
echo "== plans written to $PLAN_DIR"

if [ "$PLAN_ONLY" = 1 ]; then
  echo "== --plan-only: no VM interaction, no launch. Done."
  exit 0
fi

# ---- (0b) FREEZE GUARD: refuse to launch before the m1-freeze tag exists ----- #
if [ "$FORCE_NO_FREEZE" != 1 ]; then
  if ! git -C "$REPO" rev-parse -q --verify "refs/tags/m1-freeze" >/dev/null 2>&1; then
    echo "REFUSING TO LAUNCH: git tag 'm1-freeze' does not exist in $REPO." >&2
    echo "The M1 confirmatory campaign must not run before the freeze is tagged." >&2
    echo "Tag the freeze first, or pass --force-no-freeze to override (NOT recommended)." >&2
    exit 3
  fi
  echo "== freeze guard OK: m1-freeze tag present."
fi

plain() { local home="$1"; shift; HOME="/var/tmp/codex-colab-accounts/$home" timeout 240 "$C" --auth oauth2 "$@"; }

# ---- rebuild the bundle fresh (deterministic; prints sha) -------------------- #
bash "$REPO/code/round63/colab/m1_make_bundle.sh" --out "$STAGE" >"$STAGE/bundle_build.log" 2>&1 \
  || { echo "BUNDLE_BUILD_FAILED"; tail -8 "$STAGE/bundle_build.log"; exit 1; }
sha256sum "$BUNDLE" 2>/dev/null || { echo "BUNDLE_MISSING $BUNDLE"; exit 1; }
SPLIT=0
compgen -G "$BUNDLE.part_*" >/dev/null 2>&1 && SPLIT=1

launch_one() {   # account session tag
  local acct="$1" sess="$2" tag="$3"
  echo "==== $sess ($acct, plan m1_$tag) ===="

  # 1. nonce liveness (never trust a cached sessions row; trivial print only)
  local nonce="r63m1nonce_${sess}_$$"
  printf 'print("%s")\n' "$nonce" > "$STAGE/nonce_$sess.py"
  local got
  got=$(plain "$acct" exec --session "$sess" --file "$STAGE/nonce_$sess.py" --timeout 120 2>&1 | tr -d '\r')
  echo "$got" | grep -q "$nonce" || { echo "LIVENESS_FAIL $sess: $got"; return 1; }
  echo "  liveness OK"

  # 2. upload bundle (whole, or split parts)
  if [ "$SPLIT" = 1 ]; then
    local i=0
    for p in "$BUNDLE".part_*; do
      plain "$acct" upload --session "$sess" "$p" "$BUNDLE_TGZ_VM.part_$(printf '%03d' "$i")" >/dev/null \
        || { echo "UPLOAD_FAIL $sess part $i"; return 1; }
      i=$((i+1))
    done
    echo "  upload OK ($i parts)"
  else
    plain "$acct" upload --session "$sess" "$BUNDLE" "$BUNDLE_TGZ_VM" >/dev/null \
      || { echo "UPLOAD_FAIL $sess"; return 1; }
    echo "  upload OK"
  fi

  # 3. bootstrap: (reassemble if split) untar + pin numpy/scipy + make dirs
  cat > "$STAGE/boot_$sess.py" <<PYEOF
import subprocess, sys, os
os.chdir('/content')
if ${SPLIT}:
    import glob
    parts = sorted(glob.glob('${BUNDLE_TGZ_VM}.part_*'))
    with open('${BUNDLE_TGZ_VM}', 'wb') as out:
        for p in parts:
            with open(p, 'rb') as fh:
                out.write(fh.read())
subprocess.run(['tar','xzf','${BUNDLE_TGZ_VM}'], check=True)
r = subprocess.run([sys.executable,'-m','pip','install','-q',
                    'numpy==1.26.4','scipy==1.13.1'], capture_output=True, text=True)
print('PIP_RC', r.returncode)
if r.returncode: print(r.stderr[-800:])
chk = subprocess.run([sys.executable,'-c',
    'import numpy,scipy;print("VERSIONS",numpy.__version__,scipy.__version__)'],
    capture_output=True, text=True)
print(chk.stdout.strip() or chk.stderr[-300:])
for d in ('${BUNDLE_ROOT_VM}/code/round63/colab/plans',
          '${BUNDLE_ROOT_VM}/results/round63_m1/status',
          '${BUNDLE_ROOT_VM}/results/round63_m1/shards'):
    os.makedirs(d, exist_ok=True)
print('BOOT_OK')
PYEOF
  local boot
  boot=$(plain "$acct" exec --session "$sess" --file "$STAGE/boot_$sess.py" --timeout 420 2>&1)
  echo "$boot" | grep -q "BOOT_OK" || { echo "BOOT_FAIL $sess: $boot"; return 1; }
  echo "$boot" | grep -E "VERSIONS|PIP_RC"
  echo "$boot" | grep -q "VERSIONS 1.26.4 1.13.1" || { echo "ENV_PIN_FAIL $sess"; return 1; }

  # 4. upload this session's lane plan into the unpacked bundle tree
  plain "$acct" upload --session "$sess" "$PLAN_DIR/session_m1_$tag.txt" \
    "$BUNDLE_ROOT_VM/code/round63/colab/plans/session_m1_$tag.txt" >/dev/null \
    || { echo "PLAN_UPLOAD_FAIL $sess"; return 1; }
  echo "  plan uploaded"

  # 5. detached driver launch (BACKGROUND_LAUNCHED pattern), 6 lanes
  cat > "$STAGE/launch_$sess.py" <<PYEOF
import os, subprocess
B='${BUNDLE_ROOT_VM}'
tag='${tag}'
plan=f'{B}/code/round63/colab/plans/session_m1_{tag}.txt'
os.makedirs(f'{B}/results/round63_m1/status', exist_ok=True)
cmd=['setsid','nohup','bash',f'{B}/code/round63/colab/session_driver.sh',plan,
     '--bundle-root',B,'--session',f'm1_{tag}',
     '--status-json',f'{B}/results/round63_m1/status/m1_{tag}.json',
     '--n-lanes','${N_LANES}','--wall-budget-s','${WALL_BUDGET_S}']
p=subprocess.Popen(cmd, stdout=open(f'{B}/driver_m1_{tag}.log','w'),
                   stderr=subprocess.STDOUT, start_new_session=True)
print('BACKGROUND_LAUNCHED', p.pid)
PYEOF
  local lo
  lo=$(plain "$acct" exec --session "$sess" --file "$STAGE/launch_$sess.py" --timeout 120 2>&1)
  echo "$lo" | grep -q "BACKGROUND_LAUNCHED" || { echo "LAUNCH_FAIL $sess: $lo"; return 1; }
  echo "  driver launched: $(echo "$lo" | grep BACKGROUND_LAUNCHED)"

  # 6. keep-alive (endpoint parsed from the live sessions receipt; setsid)
  local line ep
  line=$(plain "$acct" sessions 2>&1 | grep -F "[$sess]")
  ep=$(printf '%s' "$line" | sed -nE 's/^\[[^]]+\] ([A-Za-z0-9._-]+) \|.*/\1/p')
  [ -n "$ep" ] || { echo "ENDPOINT_PARSE_FAIL $sess: $line"; return 1; }
  if pgrep -f "keep-alive $ep $sess" >/dev/null; then
    echo "  keep-alive already running ($ep)"
  else
    HOME="/var/tmp/codex-colab-accounts/$acct" setsid nohup timeout 86400 \
      "$C" --auth oauth2 keep-alive "$ep" "$sess" \
      >"$STAGE/ka_$sess.log" 2>&1 < /dev/null &
    sleep 3
    pgrep -f "keep-alive $ep $sess" >/dev/null \
      && echo "  keep-alive attached ($ep)" || { echo "KEEPALIVE_FAIL $sess"; return 1; }
  fi
  return 0
}

fails=0
launch_one pro2 r63m1_a   a   || fails=$((fails+1))
launch_one pro2 r63m1_b   b   || fails=$((fails+1))
launch_one pro2 r63m1_c   c   || fails=$((fails+1))
launch_one pro1 r63m1_p1a p1a || fails=$((fails+1))
launch_one pro1 r63m1_p1b p1b || fails=$((fails+1))
echo "==== M1 LAUNCH SUMMARY: failures=$fails ===="
exit "$fails"
