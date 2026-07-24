#!/bin/bash
set -e
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
H=/var/tmp/codex-colab-accounts/pro2
S=cbt_replica
D=/mnt/d/GI_another/results/round63_next
echo "=== upload 3 frozen files (flat to /content) ==="
HOME=$H timeout 240 "$C" --auth oauth2 upload --session "$S" "$D/CURVED_BLINDNESS_TEST/curved_blindness_test.py" /content/curved_blindness_test.py 2>&1 | tail -3
HOME=$H timeout 240 "$C" --auth oauth2 upload --session "$S" "$D/SCRAMBLE_EXT/scramble_toy.py" /content/scramble_toy.py 2>&1 | tail -3
HOME=$H timeout 240 "$C" --auth oauth2 upload --session "$S" "$D/JET_TEST/jet_test.py" /content/jet_test.py 2>&1 | tail -3
HOME=$H timeout 240 "$C" --auth oauth2 upload --session "$S" "$D/CURVED_BLINDNESS_TEST/colab_a100/launcher.py" /content/launcher.py 2>&1 | tail -3
echo "=== exec launcher (places files, verifies md5, background-launches run) ==="
HOME=$H timeout 300 "$C" --auth oauth2 exec --session "$S" --file "$D/CURVED_BLINDNESS_TEST/colab_a100/launcher.py" --timeout 600 2>&1 | tail -25
