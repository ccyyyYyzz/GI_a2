# Application wounds for the beyond-band change sentinel (hostile-scouted, 2026-07-24)

Scout deliverable (full text in session task archive). Capability scored
against: static target + sub-resolution change + thin near-field
multiplicative fluctuating screen, stationary over tens of seconds; one
photodiode + band-limited structured illumination; 2–25 s; theorem-grade
specificity (mean/intensity instruments have exactly zero response;
in-band/beyond-band/medium events separate on orthogonal statistics).

## Ranked wounds

1. **Laser welding / laser-DED metal AM defect-onset QC (fit 7.5).** The
   incumbent sensor is LITERALLY a single photodiode reading mean/RMS
   emission, and it is documented that porosity and cracks cannot be
   sensed from those intensity signals (while keyhole/penetration can);
   fume/spatter fouls coaxial cameras. Stakes: metal AM ≈ $6.7B (2025);
   qualification "millions of dollars and 5–15 years"; ~10% scrap;
   porosity is the qualification-blocking defect. CAVEAT (honesty flag):
   LPBF melt pools are microsecond/moving — the honest scope is slower
   DED/WAAM/continuous-seam, reframed as a per-track process-regime-change
   sentinel. Sources: MDPI Processes 13(1):121; ScienceDirect
   S003039922502136X; DED reviews PMC12472029, S1526612523010964;
   10.1080/17452759.2021.2018938. VERIFY primary PDFs before verbatim
   quoting the photodiode-blindness sentence (publisher 403s blocked the
   scout).
2. **Thermal/plasma spray coating monitoring (fit 7.5 — regime-cleanest
   twin).** Coating integrity known only post-process by destructive
   metallography; in-flight diagnostics see particles, not the growing
   coating; substrate quasi-static over seconds; spray cloud = thin
   fluctuating screen. Stakes: $11–15B, aerospace ~35%, no repair path
   (strip+recoat). Lead the DEMONSTRATION here; headline the wound with
   welding. Sources: Struers metallography practice; JTST
   10.1007/s11666-017-0543-8, -0559-0.
3. **Multimode-fiber endoscopy drift (fit 3 now / ~7 after scrambling-
   channel theory extension).** Real, funded wound (TM decalibration in
   vivo); WRONG CHANNEL for v1 (random-unitary, not multiplicative).
   Future-extension target only. Sources: Nat. Photon.
   s41566-023-01240-x; Comms Phys s42005-023-01410-x.
4. Bioreactor/cell-culture QC (fit 5, geometry awkward — medium is the
   target). 5. Furnace/boiler (fit 4 — blackbody self-emission swamps a
   passive diode; refractive turbulence dominates). 6. Packaging tamper
   (fit 4 — film is static; our USP absent). 7. Deep tissue (fit 2 —
   diffusive regime; motivation only, never a claimed operating point).

## Explicit misfits (the paper must not claim)

Underwater/turbid-flow inspection (thick diffusive + sonar already
solves it); deep-tissue targets; furnace interiors; corrosion-under-
insulation (static opaque medium); translucent-packaging (static film);
MMF before the scrambling-channel extension.

## The one-liner (adoption pitch)

Industry's photodiodes have run in mean-mode for decades; a covariance
firmware update turns the same diode into a sub-resolution sentinel —
the instrument that cannot see the pore can hear its birth.

## Chen-bench mapping (post-publication pitch)

Rotating ground-glass diffuser = the fluctuating screen (a second,
independently-driven diffuser emulates a medium-change event for the
three-way separation demo); DMD = the band-limited projector (coarser
than the target feature); bucket photodiode = readout; static resolution
chart with a toggleable sub-DMD-resolution feature = the change.
Deliverables: detection in 2–25 s; clean {in-band, beyond-band, medium}
separation; graceful degradation vs a controlled convolutive blend to
~50%. Minimal new hardware on their existing single-pixel-through-
scattering apparatus, operated as a sentinel rather than an imager.
