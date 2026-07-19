# Analysis correction (R19)

This additive correction preserves the immutable `m1-freeze` tag, raw cells,
original analyzer outputs, and `M1_VERDICTS.json`. After unblinding, an
independent audit found that the frozen scorer used `log(nu*rho)` instead of
the preregistered elapsed optical-time coordinate, resampled family labels,
and did not reproduce the frozen PAVA/censoring machinery. The promoted
independent implementation and in-project rerun agreed on 18/18 audited
numbers with zero deep-diff leaves.

The corrected source of record is `M1_VERDICTS_SPEC_CORRECTED_R19.json`.
The `nu*rho` result is retained only as a post-hoc incident-exposure
sensitivity and carries no verdict. Original outputs remain available in
`../M1_VERDICTS.json` and `../AUDIT_2026-07-19/`.
