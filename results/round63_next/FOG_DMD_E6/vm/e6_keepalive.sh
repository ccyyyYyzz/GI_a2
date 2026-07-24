#!/bin/bash
# Persistent keep-alive daemons for both e6 sessions. Launch via Start-Process wsl.exe.
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
H=/var/tmp/codex-colab-accounts/pro2
setsid bash -c "HOME=$H $C --auth oauth2 keep-alive gpu-l4-s-kkb-ass1b2-876uvi64s25f e6s1 >/tmp/e6_ka1.log 2>&1" &
setsid bash -c "HOME=$H $C --auth oauth2 keep-alive gpu-l4-s-kkb-ass1b0-32rhhwq9z4gca e6s2 >/tmp/e6_ka2.log 2>&1" &
wait
