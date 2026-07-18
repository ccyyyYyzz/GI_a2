#!/bin/bash
# Re-adoption probe: call `new` once on pro2; if the returned endpoint matches
# one of the three orphaned endpoints, assignment-pool adoption works and the
# running VMs (with our drivers) are recoverable with fresh tokens.
set -u
C=/var/tmp/codex-colab-tools/colab-cli-venv/bin/colab
H=/var/tmp/codex-colab-accounts/pro2
echo "== orphans before:"
HOME=$H timeout 240 "$C" --auth oauth2 sessions 2>&1
echo "== adopting as r63_pro2x:"
HOME=$H timeout 240 "$C" --auth oauth2 new --session r63_pro2x --gpu L4 2>&1 | tail -3
echo "== state row:"
python3 -c "import json;d=json.load(open('$H/.config/colab-cli/sessions.json'));print({k:{'endpoint':v.get('endpoint')} for k,v in d.items()})"
echo "== sessions after:"
HOME=$H timeout 240 "$C" --auth oauth2 sessions 2>&1
