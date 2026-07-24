# IT8 — Dither-render leak (internal divergence finer #5)

**Verdict: PWM CONSTRAINT NEEDED — ordered-dither rendering leaks 1.38e-2 into the
beyond-band annulus at z1=0 (≥1e-2), nearly doubling the analog floor.**
`IT8_dither_render_leak.py` + `.json`. float64, M=16, Bayer 8×8 ordered dither, ~2 s.

## Measured (beyond-band mean leak)
| render | z1=0 | z1=10mm |
|---|---|---|
| analog (grayscale hold) | 7.02e-3 | 5.59e-2 |
| **dither (Bayer 8×8)** | **1.38e-2** | 5.73e-2 |

- At z1=0 the ordered-dither render **doubles** the beyond-band leak (7.0e-3 → 1.38e-2),
  crossing the 1e-2 threshold: the binary dither pattern adds mirror-scale spatial
  structure that the complementary differencing does not cancel.
- At z1=10mm the dither excess is small (5.73e-2 vs 5.59e-2) because the propagation
  pedestal already dominates.

## Reading
Binary spatial dithering of the DMD grayscale amplitude is **not** a non-issue: it adds a
~1.4e-2 beyond-band leak at the conjugate plane (z1=0), on the order of the whole
pixelation floor. **PWM temporal dithering with exact whole-period bucket integration**
(so the time-averaged amplitude is realized without a static spatial dither pattern)
should be a hard line in the apparatus spec. Nothing here touches the Letter or sealed
artifacts.
