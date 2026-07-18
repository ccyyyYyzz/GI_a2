#!/bin/bash
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
printf 'print("alive-check-ok")\n' > /tmp/alv.py
chk() { # account session
  echo "== $2"
  HOME="/var/tmp/codex-colab-accounts/$1" timeout 240 "$C" --auth oauth2 \
    exec --session "$2" --file /tmp/alv.py --timeout 90 2>&1 | tail -2
}
chk pro1 r63_pro1a
chk pro1 r63_pro1b
chk pro2 r63_pro2a
chk pro2 r63_pro2b
chk pro2 r63_pro2c
