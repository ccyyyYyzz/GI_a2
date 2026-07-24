# R43 — Final architecture: two validated pillars, two honest dents, one manuscript decision

Both decisive exams have concluded. This is the last architecture round
before writing begins.

## Exam 1 — the R42 moonshot: MOONSHOT_SURVIVES (JET_TEST/ @ 1bf29f1)

Your Curvature-Rescued Detection / quotient-jet theorem passed its own
falsification verbatim: 17/17 checks, 0/4 kill conditions. Exact-KL
contact orders 2.038 / 4.000 (coefficients matched to 1.07% / 0.00%);
mixed-direction crossover at ratio 1.000 of eq (6.23); Monte-Carlo
detector m = 0.95 / 2.05; CUSUM delay slopes −2.16 / −3.92; the
amplitude-gauge collapse, lag-persistence, anchor restoration, iso-Q
cancellation (AUC→0.5 at ε*=−2c/d), and monotone-cone positivity all
verified; no broad blind set; the ε⁻⁴ class is expensive-not-vanishing.
The program now holds a validated PRL-shaped theorem with its log-log
figure (jet_slopes.pdf).

## Exam 2 — the sealed detection probe: PARTIAL (SEALED_DET/ @ b37c841)

Frozen at 5910277 before unblinding; 2.76/6.0 GPU-h; no repair round.

**PASS — D0 (wall 5.7e-16; two engines 1.7%; shot 1.0%), D1 (sealed
calibration median d'_emp/d'_ana = 0.998 [0.909, 1.099], n=108),
D2 (THE CAPABILITY: 2% headline 77/81 analytic, best cell T_det = 453
banks ≤ 600, MC power LCB 0.990; 5%: 81/81 spread, 75/81 worst-mode;
1% best cell 1811 ≤ 2048; 0.5% edge audit max d' 3.76), D4 (fixed-bank
concentration retained: best-cell latency 459 vs production-fresh 1013;
mid-cell 1271 vs 1488), D6 (online 0.020 ms/bank = 0.16% of the bank
interval; 0.26 MB).**

**FAIL — D3 (specificity): target TPR 0.988, in-band FA 0.020, balanced
accuracy 0.916 (≥0.90 ok), beyond-score |d'| on non-targets 0.02 — but
medium-amplitude FA 0.096 and medium-timescale FA 0.084 vs the 0.05
bar, and the intended per-class scores reached min d' 1.79 vs the
declared 5 at the frozen event magnitudes. FAIL — D5 (mismatch): AUC
1.000 and T_det inflation within bounds on EVERY axis (rot10% passes
everything: FA 0.032, T_det +19%), but non-target FA leaks at
rot20% (0.052) and spectral-slope −1 (0.076). Kill-tree node: D5
("narrow the physical domain or kill robustness wording"); D3's node
kills the strong "specific sentinel" wording.**

Reading: the sealed exam certified the DETECTION POWER completely and
deflated the SPECIFICITY CALIBRATION by a few percent — the alarm
detects exactly as predicted everywhere; its false-alarm rate under
medium events and law mismatch runs 5–10% instead of the frozen 5%.

Integrity disclosure: after the freeze commit but before unblinding,
the record generator was refactored for GPU memory chunking (audited
diff: chunking only; no thresholds, banks, scenes, or estimator
changes; realizations differ from the unchunked stream by RNG order).
Frozen wording for this is requested below.

## R43 asks

1. **The center of gravity.** Two candidate cuts: (a) the JETS-centered
   manuscript — the validated observability-order theorem leads (PRL
   register: "statistical observability has an order"), with the
   sealed detector as its severe optical instantiation (its certified
   D2/D4/D6 as the demonstration; its D3/D5 dents as measured
   boundaries); (b) the R41 sentinel-centered cut with claims narrowed
   per the kill nodes (Optica register). Rule which, or a third cut.
   Note the jets theorem needs NO specificity claim — the D3/D5 dents
   do not touch it.
2. **One-paper discipline vs PRL Letter + companion.** The operator's
   standing rule is ONE good paper. A PRL Letter (jets + scrambling
   instantiation + the latency-slope figure) with a long companion
   (the full atlas/sentinel) is a classic legitimate pattern — rule
   whether it violates the discipline or fulfills it, and if a Letter:
   exactly what goes in its 3.5 pages.
3. **Final frozen wording set** under the kill nodes: (i) the detection
   capability sentence (power certified; specificity stated as measured
   numbers — balanced accuracy 0.916, medium-event FA ≤ 0.10 at the 1%
   threshold — never as "zero false alarm"); (ii) the robustness domain
   sentence (primary domain rot ≤ 10% where ALL bars pass; rot 20% and
   slope −1 as the mapped FA boundary); (iii) the D3/D5 disclosure
   sentences; (iv) the generator-chunking integrity sentence; (v) the
   jets theorem statement for the main text.
4. **Venue ruling, final.** Given: a validated PRL-shaped theorem + a
   sealed-certified capability with honest dents + sim-only. PRL first
   with the jets cut? Optica with the sentinel cut? Sequencing if a
   Letter+companion is ruled legitimate.
5. **Final page map + figures** for whichever cut you rule: the jets
   log-log figure and the scrambling three-row inversion (mean support
   → {0}; rank → 1; jets survive) are new assets; the R41 4-figure plan
   predates both exams.
6. **Writing green-light protocol**: confirm that after this ruling the
   manuscript build proceeds with NO further campaigns or probes — all
   remaining work is prose, figures from frozen tables, and the
   supplement (S-notes per R40/R41 maps + JET_TEST + SCRAMBLE_EXT +
   SEALED_DET as new S-notes).

Deliver as a GitHub issue titled R43. Both exam artifact sets are
committed (JET_TEST @ 1bf29f1, SEALED_DET @ b37c841). As standing: add
your own innovations section if the two exams suggest a sharper cut we
have not seen.
