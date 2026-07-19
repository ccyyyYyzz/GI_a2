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
