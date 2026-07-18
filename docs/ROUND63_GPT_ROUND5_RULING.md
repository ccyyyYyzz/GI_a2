# GPT round-5 final ruling (2026-07-18, 3m53s reply) — verbatim extract

Chat 候选方案评审与建议, 5th assistant message (2,819 chars; one ~140-char
DLP-blocked hole around the corr_obs/q_corr definition — reconstructed by
symmetry with q_mean, marked [inferred]). THE F1-frozen text for the λ_TV/
audit layer, superseding the round-4 binary gate.

---

## Final frozen ruling: no binary adequacy gate

Delete MODEL_FAIL_PREDICTIVE entirely. Do not replace it with
MEAN_RESIDUAL_WARN, and do not roughen the plug-in generator.

The probe invalidates the binary gate on both size and power: it falsely
rejects 7/20 matched-model cells while detecting 0/5 paralyzable and 0/5
τ+20% cells. The underlying non-identifiability is structural: at one
count-only operating point, detector mean mismatch can be absorbed into
scene scale.

The DEV-only K=5, one-SE η⋆ selection remains unchanged. No audit statistic
may affect η⋆, reconstruction, cell inclusion, campaign gates, or
confirmatory inference.

## Frozen descriptive audit

Run the audit once per RQL cell and record:

1. AUDIT_STATUS ∈ {OK, UNDERPOWERED, NA_DEPENDENT}.
2. D_obs, D_star_mean, D_star_sd, and D_ratio = D_obs / D̄⋆.
3. With all B_diag = 39 replicates completed, record
   q_D = (1 + #{D⋆_b ≥ D_obs}) / 40 as plugin_upper_rank, NOT as a calibrated
   p-value.
4. mean_r_obs and q_mean = (1 + #{|r̄⋆_b| ≥ |r̄_obs|}) / 40.
   Keep MEAN_RESIDUAL_WARN = (q_mean == 0.025) as a descriptive marker only.
5. corr_obs [inferred by symmetry: and q_corr = (1 + #{|c⋆_b| ≥ |c_obs|})/40,
   LOAD_CORR_WARN = (q_corr == 0.025) descriptive only; NA if undefined].
   For the fixed-λ̂, B_leak = 199 lower-tail diagnostic:
   q_low = (1 + #{D(0)_b ≤ D_obs}) / 200, with
   LEAKAGE_SUSPECT = (q_low <= 0.01).
6. Produce the residual-vs-load plot.

Disable early stopping. Once the output is continuous ranks rather than a
binary decision, early-stopped values are not comparable. Current paralyzable
rows often stop after one replicate, so their apparent warning discrimination
is not statistically interpretable; the code currently computes warnings from
that truncated replicate set.

Remove these fields: MODEL_FAIL, MODEL_FAIL_PREDICTIVE, P_FAIL, GOF_OK as a
model-acceptance claim.

## Preregistered wording

> The measurement audit is descriptive and does not constitute a
> model-adequacy test. All audit statistics are excluded from regularization
> selection, reconstruction, cell inclusion, and confirmatory pass/fail
> decisions. The refit-bootstrap reference distribution is conditional on a
> smoothed plug-in scene; its empirical tail ranks are therefore not
> calibrated model-class p-values for arbitrary scenes and can be
> anti-conservative in the upper tail when the true scene is rougher than the
> plug-in. Moreover, count-only observations collected at a single operating
> point cannot generally distinguish detector mean-response mismatch from
> radiometric rescaling of the scene. Detector-mismatch consequences are
> therefore quantified directly through the preregistered radiometric metrics
> in S3.

This is the honest and review-safe endpoint. F1 may freeze with the binary
gate deleted.
