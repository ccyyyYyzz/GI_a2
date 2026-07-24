# KT6 — Mean-Covariance Channel-Ratio Ranging (R44 §5.5, innovation #5, score 80)

**Verdict: ORACLE confirms the z2⁴ ratio law; BLIND ranging is ESTIMATOR_KILLED —
exactly the coordinator's binding caveat.** `KT6_channel_ratio_ranging.py` +
`.json`. Committed result is the M=8 run (complete, 4 cells); the M=16 run confirms
(first cell: exponent 4.08, signal/floor 5.5e-3) and larger M only *strengthens* the
kill (the Wishart floor grows as M²/T).

## Two required versions (binding: oracle-only is NOT a pass)
Dense z2 sweep {0.05…5}mm, z1∈{5,10}mm, screens {weak σ=0.3, developed σ=2π}.

| cell (z1, screen) | ORACLE near-contact exponent | BLIND median signal/floor (near) |
|---|---|---|
| 5, weak | 4.08 | 4.6e-3 |
| 5, developed | 3.95 | 4.5e-3 |
| 10, weak | 4.02 | 4.4e-3 |
| 10, developed | 3.93 | 4.3e-3 |

- **ORACLE (P_oracle_z4_law: PASS)** — J_cov/J_mean ∝ z2⁴ near contact confirmed
  (exponents 3.93–4.08 ≈ 4), instantiating R44 eq 5.3.
- **BLIND (P_blind_estimable: ESTIMATOR_KILLED)** — the finite-bank Wishart null floor
  `λ_null = 0.5 tr[(V⁻¹(Ĉ₀ₐ−Ĉ₀ᵦ))²]` (two independent H0 half-samples) sits ~220×
  **above** the near-contact covariance signal (signal/floor ~4.4e-3 ≪ 1). Blind
  λ_cov cannot be separated from sampling noise near contact.

## Reading (honors the coordinator's binding caveat)
The channel-ratio depth observable is a **valid theory relationship** (oracle z2⁴ law
holds cleanly) but is **NOT bench-observable via naive blind estimation**: near
contact the covariance Fisher is far below the Wishart M²/T sampling floor. This
matches the internal estimator-layer finding that blind λ_cov sits 600–9000× under
the floor (our conservative M=8 shows ~220×; larger M pushes it further under). KT6 is
therefore an ESTIMATOR_KILLED negative on the blind bench observable, with the oracle
law recorded as a theory result. Nothing here touches the Letter or sealed artifacts.
