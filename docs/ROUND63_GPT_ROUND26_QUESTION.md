# R26 — New theory mainline: the information theory of ghost imaging through dynamic scattering media

Program pivot (operator directive): the primary THEORY thread now moves
into the dynamic-scattering GI problem domain (the Peng–Chen line:
per-frame multiplicative gain from dynamic/complex media, supervised
corrections, single-pixel readout). This is theory-only work; the existing
two ROUND63 papers and the RLMI bridge continue independently. This round
architects the new program.

## 1. The channel and the two assets we hold

Channel (general form of Peng & Chen APL 124, 181104 (2024), Eq. 2):

    y_n = a_n · (m_n^T x) + noise_n ,   n = 1..N

with a_n = exp(l_n) a hidden stochastic gain path of correlation time t_c
(their instance: slow deterministic exponential; general: OU/low-pass
class), patterns m_n programmable, and readout either Gaussian/analog
(their bench) or photon-counting (our extension).

Asset A — an existing identifiability theory (operator's separate
manuscript, summarized): square designs are algebraically non-identifiable
mod gauge (invertible M: T' = M^{-1} diag(a) M T reproduces data); tall
thresholds N >= K+p-1 (p = drift-class dimension) for local
identifiability, N >= 2K+p-1 uniform; a stationarity anchor: windowed-log
estimators are consistent iff a stationary-carrier condition holds;
temporal randomization (SRHT/permutation) makes the object's carrier
stationary ("the object stops pretending to be time"); a master
finite-noise relMSE identity with a residual-gain "static-correction
ceiling" term v·B_L; a design flip boundary rho* between paired-Hadamard
and randomized designs. All simulation-verified.

Asset B — the ROUND63 hidden-state information machinery: the exact
missing-information identity (loss = E[Var(hidden|observed)]), the
controlled conditional-score Gramian architecture (latent score -> channel
contraction -> Gramian -> variational design), long-window renewal
evaluations, KKT/design machinery, and the standardized-maximin allocation
theory (R22-R25 in this repo's docs/).

## 2. The five objects to architect (rank, refine, or refute)

- **O1 — the scattering missing-information identity.** For the gain
  channel with a known drift-class prior on l = (l_1..l_N):
  I_observed(x) = I_complete(x) − (a conditional-covariance functional of
  the hidden gain path given the record). Make it exact and matrix-valued
  on the object parameters; identify what plays the role of "live time".
  Expected corollary: the Asset-A static-correction ceiling as the
  long-window scalarization.
- **O2 — the temporal ridge / operating map.** Information rate about x as
  a function of the schedule: frame time T_f, pattern ordering, and revisit
  structure relative to t_c and the photon budget. Conjecture: a ridge in
  T_f/t_c — too fast wastes photons per frame (readout/shot floor), too
  slow makes gain effectively per-frame-free (identifiability cost p -> N).
  Deliver the correct axes of the "operating map for GI through dynamic
  media" and the law along the ridge, with the same rigor ladder we used
  for the dead-time ridge (exact identity -> long-window law -> matched
  asymptotics, honestly labeled).
- **O3 — the correction-ceiling theorem.** A Fisher/van-Trees bound on ANY
  estimator of x (network or not) given the gain-class prior — the
  theoretical yardstick for supervised corrections (how far is SCGI's
  U-Net from the limit; when is more training pointless). State what extra
  assumptions (Gaussianity, known drift class, known pattern intensities)
  are load-bearing.
- **O4 — optimal temporal design.** Efficiency companion to Asset A's
  identifiability: which schedules maximize the O1 information —
  derive interleaving/chopping/randomization as OPTIMAL solutions of the
  design problem (measure-valued, KKT certificates as in R25) rather than
  heuristics; connect to Asset A's rho* flip boundary.
- **O5 — the composite channel (reach probe, one section only).** Gain
  drift + photon-counting/dead-time readout: does the Gramian architecture
  compose cleanly (R22's chain rule), and which NEW phenomenon appears
  first (e.g. gain drift aliasing into the count nonlinearity)?

## 3. Questions

- **Q1.** Prior-art sweep, aggressive: information limits for imaging
  through dynamic/turbid media, speckle-decorrelation bounds, blind-gain
  estimation (astronomy self-calibration, MRI Nyquist ghosts, microscopy
  flat-fielding), random-gain channels in communications (fading channel
  capacity analogues!), time-varying system identification. What of O1–O4
  already exists under other names? Citations, not rediscovery.
- **Q2.** Rank O1–O4 by (theorem probability × novelty × depth ×
  Chen-domain relevance). Identify THE first theorem to prove and sketch
  its proof at R23 rigor (object, statement, strategy, first
  counterexample risk).
- **Q3.** The honest-scope ruling: what should this program NOT claim
  (image-domain fences, network-specific statements, media physics beyond
  the multiplicative model)? Where does the multiplicative-gain model
  itself break (speckle is multiplicative only after bucket integration
  under which conditions?) — state the physical validity boundary as a
  hypothesis to carry.
- **Q4.** Program architecture: is this ONE paper or two (identifiability
  [exists] + efficiency/design [new])? Where does it live relative to the
  operator's existing identifiability manuscript — extension, companion,
  or merge? (Structure recommendation only; no writing this round.)

Deliver as a GitHub issue on ccyyyYyzz/GI_a2 referencing this document.
House rules as R22–R25: brutal candor, citations aggressive, conjectures
labeled with obstacles named.
