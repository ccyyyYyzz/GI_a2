#!/bin/bash
LOG=/d/GI_another/results/round63_study2_primary.log
while true; do
  if grep -q "cell 90/90" "$LOG" 2>/dev/null; then echo "PRIMARY_COMPLETE"; exit 0; fi
  if ! tasklist.exe 2>/dev/null | grep -qi python; then echo "PRIMARY_PROCESS_GONE_BEFORE_90"; exit 1; fi
  sleep 30
done
