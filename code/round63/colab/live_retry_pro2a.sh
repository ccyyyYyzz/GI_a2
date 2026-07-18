#!/bin/bash
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
STAGE=/var/tmp/r63_s2_stage
H=/var/tmp/codex-colab-accounts/pro2
S=r63_pro2a

# sanity: bundle + env really are in place from the successful boot step
cat > "$STAGE/verify_pro2a.py" <<'PYEOF'
import os, subprocess, sys
ok = os.path.isdir('/content/round63_s2_bundle/code/round63')
chk = subprocess.run([sys.executable,'-c',
    'import numpy,scipy;print("VERSIONS",numpy.__version__,scipy.__version__)'],
    capture_output=True, text=True)
print('BUNDLE_PRESENT', ok)
print(chk.stdout.strip())
PYEOF
out=$(HOME=$H timeout 240 "$C" --auth oauth2 exec --session "$S" --file "$STAGE/verify_pro2a.py" --timeout 120 2>&1)
echo "$out" | grep -E "BUNDLE_PRESENT|VERSIONS"
echo "$out" | grep -q "BUNDLE_PRESENT True" || { echo "BUNDLE_MISSING -> need full relaunch"; exit 2; }
echo "$out" | grep -q "VERSIONS 1.26.4 1.13.1" || { echo "ENV_PIN_MISSING -> need re-boot step"; exit 3; }

# retry the detached driver launch
lo=$(HOME=$H timeout 240 "$C" --auth oauth2 exec --session "$S" --file "$STAGE/launch_r63_pro2a.py" --timeout 120 2>&1)
echo "$lo" | grep -q "BACKGROUND_LAUNCHED" || { echo "LAUNCH_FAIL_AGAIN: $lo" | head -5; exit 1; }
echo "driver launched: $(echo "$lo" | grep BACKGROUND_LAUNCHED)"

# keep-alive
line=$(HOME=$H timeout 240 "$C" --auth oauth2 sessions 2>&1 | grep -F "[$S]")
ep=$(printf '%s' "$line" | sed -nE 's/^\[[^]]+\] ([A-Za-z0-9._-]+) \|.*/\1/p')
[ -n "$ep" ] || { echo "ENDPOINT_PARSE_FAIL: $line"; exit 1; }
HOME=$H setsid nohup timeout 86400 "$C" --auth oauth2 keep-alive "$ep" "$S" \
  >"$STAGE/ka_$S.log" 2>&1 < /dev/null &
sleep 3
pgrep -af "keep-alive $ep $S" >/dev/null && echo "keep-alive attached ($ep)" || { echo "KEEPALIVE_FAIL"; exit 1; }
echo "PRO2A_RETRY_OK"
