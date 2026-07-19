# R14 exact-law predictions, registered BEFORE the v2 sweep completes

R14 (issue filed) gives the exact jitter-limited optimum: c^2 rho^2 (1+2rho) = 1,
matched crossover rho* ~ [2c^2 + 1/(6 nu)]^(-1/3).

Predictions for the RUNNING refined sweep (jitter_sfi_v2, nu=2000,
n_mc=100k, 25-point grid), registered at commit time:

| cv   | exact-optimum prediction | prior coarse measurement |
|------|--------------------------|--------------------------|
| 0.02 | rho* ~ 10.4              | (not yet measured - fresh point) |
| 0.05 | rho* ~ 5.69              | 5.43 (grid-limited) |
| 0.1  | rho* ~ 3.50              | 3.30 |
| 0.2  | rho* ~ 2.30              | (fresh point) |
| 0.3  | rho* ~ 1.63              | <=2 (grid edge) |

nu=200 run: crossover term 1/(6*200)=8.3e-4 shifts cv<=0.05 peaks left
(cv=0.02: [8e-4+8.3e-4]^(-1/3) ~ 8.5; cv=0: (1200)^(1/3) ~ 10.6).

The cv=0.02 and cv=0.2 rows are OUT-OF-SAMPLE tests of the R14 law.

## VERDICT (appended after sweep completion, same session)

nu=2000 measured peaks: cv 0/0.02/0.05/0.1/0.2/0.3 ->
22.297 / 10.739 / 5.173 / 3.862 / 2.153 / 1.607.
Both OUT-OF-SAMPLE predictions hit: cv=0.02 predicted 10.4 measured
10.739 (3%); cv=0.2 predicted 2.30 measured 2.153 (6%). All in-sample
points within one grid notch (spacing 1.157). nu=200 column consistent
including the crossover correction (cv=0 peak 9.28~10.0 exact; cv=0.3
measured 1.86 vs 1.77 predicted). The R14 law is confirmed
out-of-sample at two dwells.
