#!/bin/bash
# ROUND63 S2 live ops step 2: per session — nonce liveness, bundle upload,
# env-pin bootstrap, detached driver launch, keep-alive. File-carried per the
# runtime discipline; every CLI call wrapped in timeout.
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
REPO=/mnt/d/GI_another
STAGE=/var/tmp/r63_s2_stage
BUNDLE="$STAGE/round63_s2_bundle.tar.gz"
mkdir -p "$STAGE"

plain() { local home="$1"; shift; HOME="/var/tmp/codex-colab-accounts/$home" timeout 240 "$C" --auth oauth2 "$@"; }

# --- rebuild the bundle fresh (deterministic; prints sha) ---
bash "$REPO/code/round63/colab/make_bundle.sh" --out "$STAGE" >"$STAGE/bundle_build.log" 2>&1 \
  || { echo "BUNDLE_BUILD_FAILED"; tail -5 "$STAGE/bundle_build.log"; exit 1; }
sha256sum "$BUNDLE"

launch_one() {  # account session plan_tag
  local acct="$1" sess="$2" tag="$3"
  echo "==== $sess ($acct, plan $tag) ===="

  # 1. nonce liveness (never trust a cached sessions row)
  local nonce="r63nonce_${sess}_$$"
  printf 'print("%s")\n' "$nonce" > "$STAGE/nonce_$sess.py"
  local got
  got=$(plain "$acct" exec --session "$sess" --file "$STAGE/nonce_$sess.py" --timeout 120 2>&1 | tr -d '\r')
  echo "$got" | grep -q "$nonce" || { echo "LIVENESS_FAIL $sess: $got"; return 1; }
  echo "  liveness OK"

  # 2. upload bundle
  plain "$acct" upload --session "$sess" "$BUNDLE" /content/round63_s2_bundle.tar.gz >/dev/null \
    || { echo "UPLOAD_FAIL $sess"; return 1; }
  echo "  upload OK"

  # 3. bootstrap: untar + pin the frozen numpy/scipy + fresh-python version check
  cat > "$STAGE/boot_$sess.py" <<'PYEOF'
import subprocess, sys, os
os.chdir('/content')
subprocess.run(['tar','xzf','round63_s2_bundle.tar.gz'], check=True)
r = subprocess.run([sys.executable,'-m','pip','install','-q',
                    'numpy==1.26.4','scipy==1.13.1'],
                   capture_output=True, text=True)
print('PIP_RC', r.returncode)
if r.returncode: print(r.stderr[-800:])
chk = subprocess.run([sys.executable,'-c',
    'import numpy,scipy;print("VERSIONS",numpy.__version__,scipy.__version__)'],
    capture_output=True, text=True)
print(chk.stdout.strip() or chk.stderr[-300:])
print('BOOT_OK')
PYEOF
  local boot
  boot=$(plain "$acct" exec --session "$sess" --file "$STAGE/boot_$sess.py" --timeout 420 2>&1)
  echo "$boot" | grep -q "BOOT_OK" || { echo "BOOT_FAIL $sess: $boot"; return 1; }
  echo "$boot" | grep -E "VERSIONS|PIP_RC"
  echo "$boot" | grep -q "VERSIONS 1.26.4 1.13.1" || { echo "ENV_PIN_FAIL $sess"; return 1; }

  # 4. detached driver launch (v22 BACKGROUND_LAUNCHED pattern)
  cat > "$STAGE/launch_$sess.py" <<PYEOF
import os, subprocess
B='/content/round63_s2_bundle'
plan=f'{B}/code/round63/colab/plans/session_${tag}.txt'
os.makedirs(f'{B}/results/round63/status', exist_ok=True)
cmd=['setsid','nohup','bash',f'{B}/code/round63/colab/session_driver.sh',plan,
     '--bundle-root',B,'--session','${tag}',
     '--status-json',f'{B}/results/round63/status/${tag}.json',
     '--n-lanes','6','--wall-budget-s','43200']
p=subprocess.Popen(cmd, stdout=open(f'{B}/driver_${tag}.log','w'),
                   stderr=subprocess.STDOUT, start_new_session=True)
print('BACKGROUND_LAUNCHED', p.pid)
PYEOF
  local lo
  lo=$(plain "$acct" exec --session "$sess" --file "$STAGE/launch_$sess.py" --timeout 120 2>&1)
  echo "$lo" | grep -q "BACKGROUND_LAUNCHED" || { echo "LAUNCH_FAIL $sess: $lo"; return 1; }
  echo "  driver launched: $(echo "$lo" | grep BACKGROUND_LAUNCHED)"

  # 5. keep-alive (endpoint parsed from the live sessions receipt)
  local line ep
  line=$(plain "$acct" sessions 2>&1 | grep -F "[$sess]")
  ep=$(printf '%s' "$line" | sed -nE 's/^\[[^]]+\] ([A-Za-z0-9._-]+) \|.*/\1/p')
  [ -n "$ep" ] || { echo "ENDPOINT_PARSE_FAIL $sess: $line"; return 1; }
  HOME="/var/tmp/codex-colab-accounts/$acct" setsid nohup timeout 86400 \
    "$C" --auth oauth2 keep-alive "$ep" "$sess" \
    >"$STAGE/ka_$sess.log" 2>&1 < /dev/null &
  sleep 3
  pgrep -af "keep-alive $ep $sess" >/dev/null \
    && echo "  keep-alive attached ($ep)" \
    || { echo "KEEPALIVE_FAIL $sess"; return 1; }
  return 0
}

fails=0
launch_one pro1 r63_pro1a pro1a || fails=$((fails+1))
launch_one pro1 r63_pro1b pro1b || fails=$((fails+1))
launch_one pro2 r63_pro2a pro2a || fails=$((fails+1))
launch_one pro2 r63_pro2b pro2b || fails=$((fails+1))
launch_one pro2 r63_pro2c pro2c || fails=$((fails+1))
echo "==== LAUNCH SUMMARY: failures=$fails ===="
exit "$fails"
