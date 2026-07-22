# R31 — NEW MAINLINE: saturation-coded bucket detection (theory round; GPT-led proofs invited)

Operator direction pick (2026-07-22, overrides the R30 package as
flagship; R30's certified-imaging fence is downgraded to material bank):
the mainline is now the **defect-inversion idea** — SiPM per-cell
saturation as built-in nonlinear optics that breaks the linear-GI null
space. Operator's stated aesthetic: one exquisite theory-discovered
sledgehammer with large, robust, universal gains — not a web of
conservative micro-methods. Backup lane (no work now): dual-ledger
scene+medium imaging (the E[Var(H|R)] complementarity).

## The mechanism (candidate theorem seed — please prove or break it)

SiPM bucket detector: C microcells (10³–10⁴), each fires BINARY per gate
window; output S = number of fired cells. The nonlinearity acts BEFORE
spatial summation on per-cell intensities u_c = ρ·w_cᵀ(m⊙x), where w_c
is the collection speckle weight vector. Under fully-developed-speckle
collection (w_cp i.i.d. exponential across cells; one grain per cell):

1 − E[S]/C = E_w[e^{−ρ wᵀv}] = ∏_p (1 + ρ̃ v_p)^{−1},  v = m⊙x.

**The saturated response curve in ρ is the reciprocal generating
polynomial of the masked image; its coefficients are the elementary
symmetric polynomials e_k(v).** A power sweep therefore measures e_1
(the ordinary bucket) AND e_2 (hence Σv², the masked second moment) and
in principle higher e_k — genuinely nonlinear functionals of the scene,
outside the span of any linear bucket set. Two scenes identical to every
linear bucket (null-space pair) are separated by the saturation channel.
A fixed detector self-averages over its C cells (LLN), so under the
i.i.d.-speckle model NO per-cell calibration is needed — the ensemble
statistics are the instrument. Gamma-speckle generalization
∏(1+ρv_p/k)^{−k} preserves coefficient identifiability.

Hardware status: zero new hardware — SiPMs are the most common bucket
detectors; power sweeps are standard; the method OPERATES the detector
in the deep-saturation regime everyone currently avoids and corrects
away. The "defect" literature (SiPM saturation correction, where
nonuniform illumination is a known SYSTEMATIC) becomes our signal.

A numerical life-or-death probe is running (T1 product-formula
verification; T2 Fisher effect size at C=3600 and realistic photon
budgets; T3 null-space-pair discriminability d′ vs budget; T4 image
teaser). This round is the THEORY half.

## R31 asks (theory first — you may DO the main proofs, R11 pattern)

1. **Prove or refute the theorem seed properly**: exact conditions on
   the speckle ensemble for the product/MGF formula; quenched
   self-averaging rates at finite C (relative fluctuation of the fixed-W
   response around the ensemble curve — we need the C-scaling and
   whether C=3600 suffices); which e_k are Fisher-identifiable at
   realistic photon budgets before noise drowns them (we expect e_1, e_2
   robust; e_3+ questionable — quantify).
2. **The information geometry**: cast the fired-cell record in the HSGI
   frame — what is the exact Fisher information matrix of (e_1,…,e_k)
   from Poisson-binomial S at swept ρ; what is the optimal sweep design
   (levels, allocation) for e_2 at fixed total photons; is there a
   dead-time-style operating ridge for the saturation channel?
3. **Identifiability downstream**: K masks × {e_1,e_2} constraints =
   linear tomography of x plus linear tomography of x² (pointwise
   square). Characterize the joint feasible set: when do the quadratic
   constraints shrink the linear null space substantially (generic-x
   arguments, dimension counting, any exact recovery statements at rate
   5–10%)? This is the jailbreak's mathematical teeth.
4. **Ruthless prior-art fence** (kill it now if it's a renaming): one-bit
   / quantized compressed sensing (per-measurement quantization — ours
   is per-cell WITHIN a measurement, aggregated); occupancy/balls-into-
   bins sketches and distinct-count estimators (Flajolet-Martin lineage
   — cells as random bins!); SiPM saturation-correction literature
   (1−exp(−μ/C) response models; nonuniformity corrections); speckle
   correlometry / HBT intensity interferometry / classical thermal-light
   GI (second moments via intensity correlations — the CLOSEST classical
   relative: state exactly what is new vs measuring g² with a fast
   detector); "single-pixel imaging with detector nonlinearity" if it
   exists; random-feature moment estimation. For each: DOI, what
   survives, what must be cited, what must not be claimed.
5. **Name the object** (the operator values elegance): propose 2-3
   candidate names/acronyms for the method and the theorem.
6. If the mechanism is a renaming or the effect is provably too thin at
   physical budgets, say so plainly — before we spend a campaign.

Deliver as a GitHub issue on ccyyyYyzz/GI_a2. Reference commit stated in
chat. Depth over speed; you may write the full derivations in the issue.
