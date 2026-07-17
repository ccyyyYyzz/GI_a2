# GOF coverage probe (outcome-blind, round-4 adjudication)

Correct-model false-alarm rate (MODEL_FAIL on data whose
generator IS the assumed renewal model) and misspecified-arm
detection rate (POISSON-LIN at rho=0.6, dead time ignored):

- **RQL_M1500**: {"false_alarm_fixed": "5/5", "false_alarm_refit": "4/5"}
- **PLIN_M1500**: {"detect_fixed": "2/2", "detect_refit": "2/2"}
- **RQL_M512**: {"false_alarm_fixed": "5/5", "false_alarm_refit": "4/5"}
- **PLIN_M512**: {"detect_fixed": "2/2", "detect_refit": "2/2"}

Full rows in probe_results.json. Verdict logic: the honest
gate must have LOW false alarm on the first block and HIGH
detection on the second, in BOTH M regimes.
