# E6 -- Does measurement-design redundancy kill the bilinear collusion?

Cell: N=1024, M=96, d=48, sigma_f=0.15, OU tau=16.0, T=128 blocks x S=96 slots = 12288 exposures (fixed across arms). Photons: clean + 1e5/bucket. Seeds: [0, 1, 2, 3, 4].
Fourier bank (A1/A3/A4): lattice spacing=3, medium band=4.0, achieved sideband overlap ~= 0.625.
Null fraction of scene (fraction of energy in ker(P)) ~= 0.680.

## Verdict table (arm x estimator x endpoint)

null-NMSE (fixed-averaging baseline = 1.0; oracle ~0 = data contain the null). TRUE-SEED: refinement seeded from the true medium path, medium re-estimated -> STAY(<=0.10)=collusion killed, DRIFT(~0.7)=collusion present. BLIND=data-only cold start (5-seed median); agree=median within-seed std across starts.

| arm | photons | oracle | **PW true-seed** | PW true-seed verdict | PW blind | PW blind agree | COV blind | COV blind agree | cert SF(Fourier) |
|---|---|---|---|---|---|---|---|---|---|
| A2 | clean | 0.001 | **0.694** | DRIFT | 1.073 | 0.006 | 0.894 | 0.018 | 0.560 |
| A4 | clean | 0.001 | **0.802** | DRIFT | 1.074 | 0.006 | 0.883 | 0.015 | 0.560 |

Arm legend:
- A2: random-binary bank, CHRONOLOGY-RANDOMIZED slot schedule (GI_a1 port)
- A4: A2 base (random+slot) + multi-timescale medium + temporal-law penalty

## Endpoint 3 -- whitening certificate vs outcome (slot arms)

| arm | photons | cert SF median | corr(SF, blind null-NMSE) across seeds |
|---|---|---|---|
| A2 | clean | 0.560 | -0.971 |
| A4 | clean | 0.560 | -0.955 |

## Automatic verdict summary

NO ARM KILLED THE COLLUSION: every arm true-seed DRIFTS (>0.10) at matched refinement. Measurement-design redundancy did not convexify the pathwise medium re-estimation.
