#!/bin/bash
# ROUND63 S2 live ops step 1: create the five sessions (2x pro1, 3x pro2).
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
create() {  # account session
  local home="/var/tmp/codex-colab-accounts/$1" s="$2"
  echo "== creating $s on $1"
  HOME="$home" timeout 240 "$C" --auth oauth2 new --session "$s" --gpu L4 2>&1 | tail -3
  echo "== rc=$?"
}
create pro1 r63_pro1a
create pro1 r63_pro1b
create pro2 r63_pro2a
create pro2 r63_pro2b
create pro2 r63_pro2c
echo "== sessions receipts =="
HOME=/var/tmp/codex-colab-accounts/pro1 timeout 240 "$C" --auth oauth2 sessions 2>&1
HOME=/var/tmp/codex-colab-accounts/pro2 timeout 240 "$C" --auth oauth2 sessions 2>&1
