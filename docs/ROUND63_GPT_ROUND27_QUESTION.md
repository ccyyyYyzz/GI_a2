# R27 — Two spec collisions found by the bridge smoke (fast ruling requested)

Context: the DEV bridge is fully built (scenes, mixable library, RLMI
allocator 11/11 tests, six-arm harness, gates code 7/7). Phase-1.5 smoke
surfaced two collisions between frozen spec elements. Grid launch is held.
Evidence: T-D phase 1/1.5 reports (summarized below; harness at
code/round63/bridge_harness.py, materializer v2 in code/round63/rlmi.py).

## Collision 1 — exact-972 materialization vs the ±5% dose band (the blocker)

Empirics: every bank's CONTINUOUS measure is dose-feasible (L1–L6 exactly
0.0000 with perfectly uniform coverage; L7 0.0154). But integer realization
of 972 rows fails for 7 of 8 banks: six rounding heuristics (greedy L2
descent, spread/bit-reversal drop, max-pixel kick, random restarts,
sequential balanced, min-overlap add) stall at structured floors 0.055–1.3
for the scattered/ring/multiscale/FW banks. The prescribed per-row power
trim γ∈[0.94,1.06] cannot close them: trims are COUPLED across the k
pixels each atom covers — the min-max dose is floored by coverage
structure, not trim range.

The underlying arithmetic: a k-sparse translate family at 972 rows has
per-pixel coverage ≈ 972k/1024, so a ±1-hit change moves that pixel's dose
by ≈ 1024/(972k). The ±5% band is combinatorially reachable only when
1024/(972k) ≲ 0.05·(safety margin) — i.e. **k ≳ 21**. L0 (k=32,
purpose-built exact design, dev 0.040) is the only compliant bank; every
k=16 family is arithmetically excluded at 972 rows.

Options to rule on (T-D recommends (b); we add (d)):
(a) R25 §9 balanced-sampling/cube-method rounding — principled, expensive
    at 1024 pixel constraints, uncertain for the weighted L7 measure, and
    it CANNOT beat the k-granularity bound above (a ±1-hit is ±6.6% at
    k=16 regardless of how cleverly hits are assigned) — unless the method
    intentionally uses fractional multiplicities via power settings, which
    reintroduces the coupled-trim problem.
(b) Rebuild the geometry banks at k ≥ 32 (scattered-32 replacing
    scattered-16, ring-32, multiscale with min block ≥ 32, L7 re-solved
    with a k≥32-restricted dictionary). Preserves the band exactly;
    changes the library's identity (loses k=16 diversity — note the M1
    campaign itself found k=32 optimal, FIXED*=SCAT32, so k=16 diversity
    may be dispensable).
(c) Per-bank dose-band adjustment for geometry banks (e.g. ±7% for k=16)
    with disclosure — touches a frozen guard, so it needs YOUR explicit
    amendment if chosen.
(d) Super-row pairing: realize k=16 banks as unions of complementary
    disjoint pairs (effective k=32 rows) — keeps 16-granular geometry
    inside 32-granular rows; changes the multiplex ratio per row (physics
    disclosure needed) and needs a pairing rule for weighted L7.
Rule: which option (or combination), with exact amendments to the R25 §9
materializer spec, and whether the T-B library manifest needs regeneration
before the grid.

## Collision 2 — Gate D (<1 s) vs frozen method parameters

Measured route latency ~9.6 s on real banks. Decomposition: scenario-
matrix assembly 2.9 s + simplex oracles 3.3 s + maximin 2.5 s — the cost
is r=200 × S=16 Cholesky work, NOT atom count (L7 pruning 13174→9091
atoms changed nothing). Reaching <1 s requires reducing r or S (frozen
method parameters) or caching scenario matrices across calls. Questions:
(i) is Gate D failure a kill-gate for M2, or does it kill only the
"plugin/real-time" claim while the method may still validate as science
with an honest latency disclosure? R24 §3.4 says "the plugin claim fails
even if PSNR improves" — clarify the decision-tree consequence.
(ii) Are r/S reductions or scenario-matrix caching admissible pre-grid
engineering, or post-hoc only? If admissible, state the allowed envelope
(e.g. S=8 with a stated robustness cost; cache H0 per scene across arms).

## Also note (no ruling needed unless you disagree)

The smoke's RIDGE-SCAT32 arm at c=0.05 scores 21.96 dB vs SCAT32-060's
30.25 — the frozen calibration targets the c=0 ridge (22.25) while the
c=0.05 optimum is ~5.7. This is the jitter-cap physics appearing at image
level for the first time; the composite baseline max(SCAT,RIDGE) absorbs
it per R24 §3.3. We protected it from "fixes" in the harness notes.

Deliver as a GitHub issue on ccyyyYyzz/GI_a2 referencing this document.
Speed matters more than length this round — the grid is fueled and held.
