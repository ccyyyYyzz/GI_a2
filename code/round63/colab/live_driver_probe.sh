#!/bin/bash
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
cat > /tmp/r63_dprobe.py <<'PYEOF'
import os, subprocess, json
B = '/content/round63_s2_bundle'
print('BUNDLE', os.path.isdir(B))
p = subprocess.run(['pgrep', '-c', '-f', 'remote_lane|shard_runner'],
                   capture_output=True, text=True)
print('WORKER_PROCS', p.stdout.strip() or '0')
st = None
sd = os.path.join(B, 'results', 'round63', 'status')
if os.path.isdir(sd):
    for f in os.listdir(sd):
        d = json.load(open(os.path.join(sd, f)))
        print('STATUS', f, 'done', d.get('completed_count'), '/',
              d.get('total_shards'), 'ts', d.get('ts_utc'))
csvs = 0
rd = os.path.join(B, 'results', 'round63')
if os.path.isdir(rd):
    csvs = [f for f in os.listdir(rd) if f.endswith('.csv')]
print('CSVS', len(csvs) if isinstance(csvs, list) else 0, csvs[:4] if csvs else [])
PYEOF
probe() {
  echo "== $2"
  HOME="/var/tmp/codex-colab-accounts/$1" timeout 240 "$C" --auth oauth2 \
    exec --session "$2" --file /tmp/r63_dprobe.py --timeout 120 2>&1 \
    | grep -E "BUNDLE|WORKER|STATUS|CSVS" | head -8
}
probe pro1 r63_pro1a
probe pro1 r63_pro1b
probe pro2 r63_pro2a
probe pro2 r63_pro2b
probe pro2 r63_pro2c
