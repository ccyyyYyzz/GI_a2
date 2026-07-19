#!/bin/bash
# One-shot M1 fleet status pull: rebind-first (self-healing, guide 3.2), then
# download each session's driver status JSON and print counts ONLY (outcome-
# blind: cell/row counts + rc codes, never endpoint metric values).
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
VENVPY=/var/tmp/codex-colab-tools/colab-cli-venv/bin/python
STAGE=/var/tmp/r63_m1_stage/status_once
mkdir -p "$STAGE"

# ---- rebind (endpoint map = generation-2 fleet, 2026-07-19 15:0x UTC) ------- #
cat > /tmp/r63m1_map.json <<'EOF'
{"gpu-l4-s-kkb-ass1a0-31p0z6j43zm2y": "r63m1_a",
 "gpu-l4-s-kkb-ass1a0-1qoqoql79dk3e": "r63m1_b",
 "gpu-l4-s-kkb-ass1a1-2qniwt6j3y0pq": "r63m1_c",
 "gpu-l4-s-kkb-ass1a2-2ink89ozfh59f": "r63m1_p1a",
 "gpu-l4-s-kkb-usw4a0-64dgi267p0qa": "r63m1_p1b"}
EOF
cat > /tmp/r63m1_rebind_all.py <<'PYEOF'
import json
mapping = json.load(open('/tmp/r63m1_map.json'))
from colab_cli.common import state
from colab_cli.state import SessionState
n = 0
for a in state.client.list_assignments():
    name = mapping.get(a.endpoint)
    if not name:
        continue
    state.store.add(SessionState(name=name, token=a.runtime_proxy_info.token,
                                 url=a.runtime_proxy_info.url, endpoint=a.endpoint,
                                 variant=a.variant.name, accelerator=a.accelerator.name))
    n += 1
print('TOTAL_REBOUND', n)
PYEOF
for acct in pro1 pro2; do
  HOME="/var/tmp/codex-colab-accounts/$acct" timeout 240 "$VENVPY" /tmp/r63m1_rebind_all.py 2>&1 | tail -1
done

# ---- per-session VM ground-truth summary (exec) ----------------------------- #
# NOTE: the driver status JSON's lane table is BLIND for M1 (plans carry
# manifest paths; session_driver v1 joins hb filenames from the raw plan line).
# Ground truth = heartbeat files written by remote_lane under
# results/round63/heartbeats (bare shard_id names) + the shards output dir.
cat > "$STAGE/summ_remote.py" <<'PYEOF'
import glob, json, time
B = '/content/round63_m1_bundle'
now = int(time.time())
done = 0; fails = []; inflight = 0
cells = rows = 0; ctot = 0; ages = []
for h in sorted(glob.glob(B + '/results/round63/heartbeats/*.hb.json')):
    try:
        d = json.load(open(h))
    except Exception:
        continue
    rc = d.get("returncode")
    if rc == 0:
        done += 1
    elif rc is None:
        inflight += 1
        cells += int(d.get("cells_done") or 0)
        ctot += int(d.get("cells_total") or 0)
        rows += int(d.get("rows") or 0)
        ages.append(now - int(d.get("epoch", now)))
    else:
        fails.append("%s:rc=%s" % (d.get("shard"), rc))
ncsv = len(glob.glob(B + '/results/round63_m1/shards/*.csv'))
print("M1SUMM done=%d fail=%s inflight=%d cells=%d/%d rows=%d max_hb=%s csv=%d"
      % (done, ";".join(fails) if fails else "-", inflight, cells, ctot, rows,
         max(ages) if ages else -1, ncsv))
PYEOF
declare -A ACCT=( [a]=pro2 [b]=pro2 [c]=pro2 [p1a]=pro1 [p1b]=pro1 )
declare -A SESS=( [a]=r63m1_a [b]=r63m1_b [c]=r63m1_c [p1a]=r63m1_p1a [p1b]=r63m1_p1b )
declare -A TOTAL=( [a]=11 [b]=8 [c]=8 [p1a]=8 [p1b]=8 )
for tag in a b c p1a p1b; do
  OUT=$(HOME="/var/tmp/codex-colab-accounts/${ACCT[$tag]}" timeout 120 "$C" --auth oauth2 \
        exec --session "${SESS[$tag]}" --file "$STAGE/summ_remote.py" --timeout 60 2>&1 \
        | grep "^M1SUMM" | tail -1)
  if [ -n "$OUT" ]; then
    echo "m1_$tag  ${OUT#M1SUMM } total_shards=${TOTAL[$tag]}"
  else
    echo "m1_$tag  SUMMARIZER_EXEC_FAILED"
  fi
done
