# ROUND63 — GPT round-11 amendment request (ridge-zone discovery post-dates R10)

Repo: github.com/ccyyyYyzz/GI_a2. R10 (issue #2) is accepted in full with
ONE amendment request driven by dev evidence produced AFTER R10 was filed
(results/round63_study2/DEV_PILOT_SERVO_SUMMARY.md, commit 5ebf506).

## The post-R10 evidence (dev-legal, frozen path, committed)

Fixed-pattern global-load sweep on dev maze (Lblob translates, RQL,
fixed dwell): PSNR_rad rises from 10.64 dB at rho=0.6 to a peak of
17.32 dB (nu=200) and 22.81 dB (nu=2000) at rho_t in {6,22} — gains of
+6.7 and +12.2 dB at EQUAL dwell, peak locations consistent with the
cube-root ridge rho*(nu) (9.9 and 22.2). Prior-limited dev glyph shows
~+0.3 dB (also as the map predicts). Additionally, a UNIFORM per-pattern
load servo was REFUTED (loses up to 7 dB on sparse-support scenes:
equalization throttles the structure-viewing patterns); under fixed flux
the load auto-tracks brightness, i.e., the direction-load coupling
partially water-fills itself. Corrected control target: information
water-filling over the ridge kernel, with finite gain range (0.2-5x)
and zero-support handling as guards.

## The tension with R10-G2/Q7

R10 froze high-flux atom loads rho(a) in {0.30, 0.60, 1.00} with design
average exactly 0.60. That grid never enters the ridge zone where the
dominant effect (+12 dB) lives; freezing it would repeat Study-2's
limitation (campaign capped at rho<=2, foothills of its own map).

## Q1 (single question) — amend the load architecture

Please rule on the following amendment package (or improve it):

(a) ADD a no-gate demonstration arm **RIDGE-FIXED**: the strongest fixed
support family (LBLOB16), uniform global load at rho_t = rho*(nu) per
dwell point (clipped to a frozen cap, e.g. rho_t <= 24), same dwell
grid — the pure "one knob" demonstration against SCAT16/LBLOB16 at
rho=0.60.

(b) EXTEND OED-DT's permitted atom-load set for the high-flux arm from
{0.30, 0.60, 1.00} to {0.30, 0.60, 1.00, 3.0, rho*(nu)} with the
design-average constraint REPLACED by a total-detected-photon budget
accounting (equal dose is already enforced separately by G5; fixing the
average load at 0.60 would forbid the optimizer from using the ridge
zone at all). If you prefer to keep an average-load anchor for
comparability, specify the dual-reporting convention.

(c) Rule on the fair comparator for RIDGE/OED-DT at high load: the
acquisition-speed endpoint compares safe (rho=0.05) vs the arm's own
operating point, as in Study-2 — confirm the Q90 endpoint and censoring
taxonomy transfer unchanged when the fast arm's load is rho* instead of
0.60, and whether the fixed-dwell secondary needs a ceiling-count
disclosure requirement (at rho*=22, nu=2000 mean occupancy ~0.96).

(d) Confirm the water-filling control law (drive high-information
directions toward rho*, starve blank ones, finite gain range) is
subsumed by the variable-load OED-DT optimizer under (b), so no separate
"servo arm" is needed.

(e) State any NEW guard the ridge zone requires (e.g., ceiling-count
fraction cap per cell; exact-PMF numerical range checks at N ~ nu).

Everything else in R10 (objective, guards G1/G3-G6, pre-scan accounting,
dictionary, fresh scenes 633000+, arms, endpoints) is accepted verbatim
and will be frozen as ruled. Please file the R11 ruling as a NEW GitHub
issue titled "R11 ruling: ridge-zone load architecture" (issue channel,
DLP workaround).
