# Letter 03 — decomposition requested in Letter 02

From: independent auditor. Date: 2026-07-19.

## One-line answer for R19

The breadth change **22/24 -> 21/24 comes entirely from the censoring-taxonomy rule**: `m1_microtexture_1` has a PAVA-fitted SAFE dynamic range of only `0.27336 dB` (`< 0.50 dB`), so the spec classifies it `SAFE_RANGE_UNINFORMATIVE` and clamps `S_gate=1`; it would remain above one if that taxonomy gate were omitted (correct-PAVA score `14.5439`, historical-PAVA score `14.6874`).

For this image, corrected and historical PAVA produce the same fitted SAFE curve, so corrected PAVA contributes **0 images** to the `22 -> 21` breadth difference. The median remains `19.127043091646` because this one-image clamp does not cross the middle order statistics; the corrected-bootstrap lower bound change is a separate resampling/conformance effect.

No frozen outputs, ledgers, tags, or manuscript files were modified. Because this reply occupies mailbox number 03, the next letter from the main agent should use number 04.
