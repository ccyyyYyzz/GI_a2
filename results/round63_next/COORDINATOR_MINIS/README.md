# Coordinator mini-experiments (session scratchpad archive, 2026-07-23/24)

Dev-grade oxygen/prediction tests run by the coordinator during the fog-DMD →
covariance → beyond-band → detection arc. These are the primary sources for
several numbers quoted in R38–R40 round documents that predate the agent-run
64×64 replications. All 16×16 (N=256) unless noted; conventions match
FOG_DMD_PROBE/fog_common.py lineage. Outputs (.out) included where captured;
scripts are deterministic (seeds inline) and re-runnable from repo root with
the base Anaconda python.

| script | what it established | headline numbers |
|---|---|---|
| fog_dmd_oxygen.py / fog_dmd_oxygen2.py | first oxygen test: oracle full rank, blind ALS cold fails | rank 1024/1024; oracle 0.033–0.000; cold ~1.0 |
| fog_noise_oxygen.py | shot-noise oracle curve (unregularized) | 0.15 @ 1e6 photons |
| mini_cov_predict.py | independent covariance-route replication + T-scaling + mismatch + lattice predictions | 0.10 @ T=2048; rot10% +23%; lattice ~1.12 (poison) |
| mini_lat_fix.py | lattice failure not a rank artifact (quadrature pairs, rank 46/48) | still 1.12 |
| mini_a2_sched.py | ordered beats randomized slot schedule on covariance route (concentration) | 0.196 vs 0.335 @ 49k slots |
| mini_p1_cumulant.py | third-order cumulants do NOT unlock small-M (tower floor 3 starved) | M=12: 0.725→0.712; M=20: 0.616→0.673 |
| mini_aperture.py | aperture-law test v1 — CONFOUNDED (structured-cosine estimator collapse); superseded by v2 | in 1.04 / out 0.62 (invalid) |
| mini_aperture2.py | aperture law verified with band-limited incoherent random patterns | in 0.073 / out 1.213 (17×) — the R39 "17x separation" source |
| mini_shell_resolved.py | beyond-band recovery is shell-resolved; depth ≈ 1 shell | Δk=1: 0.318; Δk=2: 0.91; Δk=4: 2.1 |
| mini_f4_duel.py | fresh-pattern arm dominates in-band at equal budget | 0.037 vs 0.285 (8×) |
| mini_f4b_medium.py | fresh arm also reads the medium law via cross-pattern regression | σ² 7% err; τ 7.0–7.8 vs 8.0 |
| mini_f4c_crossover.py | no crossover: fresh wins all 9 (σ_f × photons) cells | 9/9 fresh |
| mini_superres_duel.py | beyond-modulator-band: fresh+mean wall exact, fixed+cov recovers | 1.000 vs 0.47–0.79 |
| mini_fresh_cov.py (+ _eq.out) | fresh+covariance fails beyond-band even at equal budget (first-pass estimator) | 0.955/0.957 |

Caveats: coordinator-written, single-implementation (except where agent
replications exist at 64×64 in FOG_DMD_PROBE64/); the shot-model bug found in
the P1/map engine (CRB_RECONCILIATION) does NOT affect these Monte Carlo
scripts (they simulate records directly), but analytic comparisons quoted
alongside them in round docs should defer to the CORRECTED engine artifacts.
