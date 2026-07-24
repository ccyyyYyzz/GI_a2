# CBT cross-venue replication — local_full vs colab_a100

The frozen `curved_blindness_test.py` (CURVED_BLINDNESS_TEST, the invariant-orbit
statistic-invariance exam) was run as **two independent arms of the same byte-identical
script**, on different hardware and with different RNG streams, as a cross-venue
replication. Both arms return the same verdict; the deterministic float64 stages match to
the last printed digit, and the stochastic Monte-Carlo stage agrees within Monte-Carlo
error (agreement is statistical, not bitwise, by design — the two arms draw disjoint
speckle streams).

| arm | venue | device | runtime | verdict | predictions |
|---|---|---|---|---|---|
| `local_full` | local GPU (PID 19876) | cuda (RTX 4060) | 7079 s | MOONSHOT_SURVIVES | 22/22 |
| `colab_a100` | Colab (pro2, `cbt_replica`) | A100 | 998 s | MOONSHOT_SURVIVES | 22/22 |

## Deterministic stages — identical to printed precision

The float64 analytic divergences (stage 1), the block-CUSUM geometry (stage 3), and the
non-optical second model (stage 4) carry no RNG and reproduce bit-for-bit:

| quantity | local_full | colab_a100 |
|---|---|---|
| efficient Fisher `I_Q` | 0.044923 | 0.044923 |
| STRAIGHT KL log-log slope | 3.99999 | 3.99999 |
| acceleration-cap kink `A_kink` | 0.209343 | 0.209343 (`A* = 0.209343`) |
| exact-cloak `max\|ΔQ\|/Q` | 1.5371575669e-16 | 1.5371575669e-16 |
| CUSUM partial-cloak delay slope | −3.947 | −3.947 |
| CUSUM exact-cloak detect_frac | 0.0000 | 0.0000 |

## Monte-Carlo stage — agreement within MC error

The exact-cloak (on-orbit great-circle) paired-detector duel uses CRN-paired records with
disjoint RNG streams between the two arms. Every cell stays inside the frozen razor band
`paired AUC ∈ [0.49, 0.51]` in both venues; the paired d′ is consistent with 0 to within
its 2σ null SE (0.0365) in every cell. The largest cross-venue difference in paired AUC is
**0.013** (th0.15, B=16384), the size expected from independent finite-MC sampling.

| cell | paired AUC local | paired AUC a100 | \|Δ\| | paired d′ local | paired d′ a100 |
|---|---|---|---|---|---|
| th0.15 B4096 | 0.4997 | 0.4937 | 0.0060 | −0.0010 | −0.0223 |
| th0.15 B16384 | 0.4924 | 0.5058 | 0.0134 | −0.0271 | +0.0205 |
| th0.35 B4096 | 0.4945 | 0.5015 | 0.0070 | −0.0195 | +0.0055 |
| th0.35 B16384 | 0.4992 | 0.4966 | 0.0026 | −0.0029 | −0.0121 |
| th0.70 B4096 | 0.4952 | 0.5051 | 0.0099 | −0.0170 | +0.0182 |
| th0.70 B16384 | 0.5013 | 0.4978 | 0.0035 | +0.0046 | −0.0078 |

Positive controls fire in both arms (unpaired AUC: STRAIGHT 0.9973 / 0.9956; PARTIAL α=0.5
0.9338 / 0.9262), so the null result on the orbit is not a power failure.

## Conclusion

The exam replicates across venues: identical verdict (MOONSHOT_SURVIVES, 22/22, no kill
triggers), bit-identical deterministic geometry, and statistically indistinguishable
Monte-Carlo detectability. The `colab_a100/` bundle is the venue-stamped replica of this
`local_full` arm.
