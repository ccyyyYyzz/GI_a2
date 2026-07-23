# R37 — DLGI v2: instrument-fix campaign under UNCHANGED bars (legitimacy + prereg ruling)

The confirmatory campaign returned FLAGSHIP_KILL (@ 2be9f65): C2 (width
2.2–2.4× t_c / 3.8–4.6× CV vs the 1.5× bar; bounded fraction ≈0),
C3 (single cell tc64_cv5_ph05 at rmse ratio 3.53; 26/27 pass, mostly
dominant), C4 (three low-SNR cells at task-fidelity 0.065–0.083 vs
0.05; PSNR passed everywhere, median +0.5 dB). C1 over-covered in all
27 cells (≥0.98 vs the 0.95 target). C5/C6/C7 passed — reciprocity
confirmed at 3.66e-15 with the preregistered schedule ordering and zero
contradictions. The R34 no-repair rule was honored: the campaign's
verdict stands permanently; the fail branch is the default.

The operator demands a genuine positive method. We propose DLGI v2 —
NOT a repair of the closed campaign but a new method version with
mechanistically-diagnosed instrument fixes, preregistered fresh, under
the IDENTICAL C1–C7 bars. Your ruling is requested on legitimacy and
the v2 prereg.

## The three diagnosed mechanisms and the proposed fixes

1. **C2 (width/boundedness)** — the calibration surface was built
   deliberately conservative: monotone UPPER-ENVELOPE interpolation +
   conservative MC order statistic + HOLD-FLAT extrapolation outside
   the calibration box (which makes Neyman regions unbounded by
   construction). Evidence it is the instrument, not the statistic:
   coverage landed ≥0.98 in every cell — 3–6% coverage slack per cell
   is being wasted on width. V2 fix: per-cell exact MC quantiles (no
   envelope), trilinear interpolation WITHOUT the conservative box-max,
   and a BOUNDED declared search domain for the inversion (the region
   is clipped to the declared (t_c, CV) domain — regions touching the
   boundary are flagged, not unbounded).
2. **C3 (one corner cell)** — tc64_cv5_ph05 is the identifiability
   corner your own R33 §4.2 predicted (slow drift = few effective
   samples; weak CV; low photons). V2 fix: honest domain scoping — the
   primary grid excludes cells your identifiability analysis marks as
   near-edge; they move to the published edge bank. No bar change: C3
   still requires EVERY primary cell to pass.
3. **C4 (three low-SNR task-margin cells)** — the scene estimator
   carried a diagnosed gauge bug (estimated gain not renormalized to
   the E[a]=1 gauge; fixing it recovers +2.19 dB at high CV with zero
   learning — results/round63_next/DL_TOOL_EXPLORATORY/ @ afbdac8).
   V2 fix: the v2 estimator includes the gauge renormalization,
   declared upfront (this is a new estimator version, hash-frozen
   before the v2 calibration bank; the v1 bit-identity rule then binds
   v2's calibration to v2's estimator).

## R37 asks

1. **Legitimacy ruling**: does v2-as-proposed constitute a legitimate
   new preregistered campaign (instrument fixes under unchanged bars),
   or bar-shopping? If any element crosses the line, name it and strike
   it. The v1 verdict remains permanently in the record either way and
   will be disclosed in any manuscript.
2. **The v2 prereg**: freeze the campaign spec — new seeds/banks
   (fresh salts), the primary-domain definition (which cells are
   in-domain given the identifiability analysis; state the rule, not a
   hand-picked list), calibration protocol (exact per-cell quantiles;
   replicates per cell), the identical C1–C7 bars verbatim, and the
   disclosure obligations (v1 verdict cited; domain narrowing stated in
   the abstract-level claim).
3. **Width prognosis**: from the calibration table's c-values and the
   confirmatory width data, estimate whether the de-conservatized
   construction plausibly clears the 1.5× width bar or whether the
   width is intrinsic (in which case say so now and we stop before
   spending the compute).
4. **Manuscript consequence**: if v2 passes, does the flagship revert
   to the DLGI-protagonist form (R36) with the v1 kill disclosed? State
   the required disclosure wording.

Deliver as a GitHub issue titled R37. Fast — the machinery is warm and
the operator is waiting.
