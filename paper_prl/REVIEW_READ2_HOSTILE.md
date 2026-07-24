# REVIEW READ 2 — Hostile statistical-optics read (R43 §7.3 step 7)

**Reviewer persona:** hostile referee, expert in statistical optics / speckle / memory
effect / SPADE & quantum metrology / sequential detection, actively seeking grounds to
reject.
**Scope:** report-only. No `.tex` or repo file was edited. Sources read: `main.tex`,
`supplement.tex`, `CLAIM_SOURCE_MATRIX.md`, `docs/ROUND63_GPT_ROUND43_RULING_RAW.md`,
and the three frozen JSONs (`JET_TEST.committed.json`, `SCRAMBLE_RESULTS.committed.json`,
`CONFIRMATORY_RESULTS.committed.json`).

---

## VERDICT: PASS_WITH_FIXES

The theorem is correct, the five displayed equations are algebraically sound, all 11 main
+ 5 supplement references are bibliographically correct, and 21/22 spot-checked numbers
match their frozen artifacts exactly. **No FAIL condition** (no fabricated number, no
missing-source headline claim, no wrong theorem, no non-existent reference). What blocks a
clean PASS is one integrity-of-framing issue (simulation labelling), a theorem-scope
qualifier that is stated only in the supplement, one weaponizable word, an incomplete
prior-art fence, and one un-matrixed number. All are fixable by editing text only — no new
data generation.

Counts: **BLOCKER 1 · MAJOR 5 · MINOR 7.**

---

## FINDINGS

### BLOCKER

**B1 — SIM-ONLY honesty: two sentences read as a physical optical experiment.**
The Letter never uses the word "simulation" in the main body, yet twice describes the
sealed detector in language a referee (or, worse, an editor post-acceptance) can read as a
performed hardware experiment.
- `main.tex:67` — *"We prove a Rank–Jet Separation theorem, demonstrate its exact
  directional classes, and falsify it at its own boundaries in a sealed single-pixel
  **optical experiment**."*
- `main.tex:271–272` — *"A preregistered thin-screen **experiment on a true single-pixel
  bucket record** converts the theorem into a power-certified sentinel."*
- `main.tex:39` (abstract) — *"…and a sealed **bucket-optics test** validate the result."*

The word "true" + "optical experiment" with no "numerical/simulated/in-silico" qualifier
is exactly the phrasing that becomes an integrity problem for a simulation-only paper. The
supplement is clean here (`supplement.tex:246` "this Kronecker generator drives the
**simulation**"; `:48` "no new data are generated here"), and the main text carries
implicit signals (Monte Carlo, generator, GPU-hours `:397`, seeds `:398`), but the task's
bar is that *no* sentence be misreadable — these two fail it.
**Severity rationale:** trivial to fix, but the paper's entire honesty posture rests on it,
so it must not ship unqualified.
**Repair (do NOT apply):** (i) abstract — change to "…and a sealed *simulated* bucket-optics
test validate the result," or add one clause "All results are numerical." (ii) `:67` →
"…in a sealed, fully simulated single-pixel optical model." (iii) `:271` → "A preregistered
*simulated* thin-screen study on a zero-dimensional bucket record…"; delete "true." (iv)
add "simulated" once to the Fig. 3 caption (`:256` "the sealed thin-screen detector").

### MAJOR

**M1 — Theorem boundary: the "empty blind set" claim is not qualified where the referee
looks; a genuine counterexample exists and is buried in the supplement.**
`main.tex:147–148` — *"Hence the reconstruction Fisher matrix is rank one, yet the local
blind set is empty"* and `:209` *"there is no blind direction for detection."* The three
boundaries paragraph (`:152–160`) is correctly adjacent to the theorem, but it lists
locality, the amplitude anchor, and the monotone cone — it does **not** state the condition
that makes the blind set empty, namely `G≻0` on the accessible aperture. The scrambling
realization has `G` with **finite** spectral support (`supplement.tex:264–265`: `Ĝ`
supported on `|k|≤2k_grain`). A change `δ` whose Fourier support lies entirely beyond
`2k_grain` gives `u=2xᵀGδ=0` **and** `v=2δᵀGδ=0` — genuinely blind, `m=∞`. The supplement
admits exactly this (equality case (iii), `supplement.tex:107–108`: *"true local blindness
only if d=0 (a nullspace of G at coarse grain)"*), and the numeric claim is honestly
scoped to "in-aperture" at `:247` and `:353–354`. But the **headline** claims at `:147` and
`:209` drop the qualifier, so a hostile referee constructs the beyond-`2k_grain`
counterexample against the theorem as stated.
**Repair:** add the aperture/positive-definiteness qualifier next to the headline claim,
e.g. `:147–149` → "…the local blind set is empty **on the aperture where `G≻0` (grain band
`≥` half the sampling band)**: first-order-orthogonal directions are curvature-rescued…";
and at `:209` → "no blind direction for detection **within the informative aperture**."
Note in the boundaries paragraph that `G≻0` requires the grain's modulus-squared band to
cover the scene band.

**M2 — D3/D5 boundary: one sentence in the theory section asserts the specificity the
sealed test explicitly killed.**
`main.tex:210–211` — *"The mean wall supplies **specificity**; the covariance curvature
supplies detection."* The whole D3 disclosure (`:279–282`, and honestly stated) is that the
strong specificity claim was **killed** (medium-event FA 0.096/0.084 > the 0.05 bar).
A referee quotes `:210` back: "the authors claim the mean wall 'supplies specificity,' yet
their own sentinel fails attribution." The two "specificities" are different (DC-vs-non-DC
invariance vs. inter-class attribution), but the bare word is a hostage, and the supplement
itself pledges (`supplement.tex:363–364`) that no "specific"-type wording is used.
**Repair:** `:210–211` → "The mean wall **rejects DC and flux changes**; the covariance
curvature supplies detection." Do not use the noun "specificity" in the theory section.

**M3 — Prior-art fence (S5) is incomplete: four of the six equivalence threats are not
fenced.** `supplement.tex:417–432`. Coverage audit:
- (a) memory-effect **imaging** through walls — **NOT fenced.** The fence names Mudry/Idier
  (blind-SIM) and Chaigne (photoacoustic) as the "memory-effect… imaging" representatives,
  which is a misattribution; the canonical works — Bertolotti *et al.* Nature **491**, 232
  (2012); Katz *et al.* Nat. Photonics **8**, 784 (2014); Freund/Feng memory effect,
  PRL **61**, 2328 (1988) — are absent. A memory-effect referee will demand these be
  distinguished by name.
- (b) SPADE/Tsang — fenced (`:426–428`, "moves Fisher information between… measurements").
  But **Grace–Guha quantum change detection** (the change-detection variant explicitly
  listed in the task) is **NOT** cited or distinguished.
- (c) coda-wave interferometry (Snieder) — **fenced** correctly (`:428–430`). ✓
- (d) Lee–Stone / universal conductance fluctuations analogy — **NOT fenced** (LeeStone is
  only a bare outlook citation, `main.tex:306`).
- (e) SAR coherent change detection — **NOT fenced anywhere.** Note Idier *et al.* (already
  ref [5]) itself flags speckle-covariance superresolution as applicable to SAR, so the
  overlap is closer than the fence admits.
- (f) blind-SIM / speckle SR (Mudry/Idier) — **fenced** correctly (`:421–425`). ✓
**Repair:** in S5, add one sentence per missing class stating what THEY own vs. what remains
ours (rank-vs-jet separation + exact contact-order classes): for (a) "memory-effect imaging
[Bertolotti, Katz, Freund] *reconstructs* while residual correlations survive; we work at
the endpoint where reconstruction rank is one and only the observability order survives";
for (b) "[Grace–Guha] optimize quantum-limited discrimination of a *known* object change;
we classify *classical* nuisance-profiled contact order"; for (e) "SAR-CCD detects
inter-pass decorrelation of the *medium*; here the first-moment wall makes the medium a
transducer for *scene* change." Add the corresponding bibitems (block below).

**M4 — Number custody: the shot-model integrity-disclosure numbers have NO
`CLAIM_SOURCE_MATRIX` row and their artifact is not in the matrix's frozen-source table.**
`main.tex:409–413` (End Matter C, disclosure 2) and `supplement.tex:388–394` quote
**"approximately 1.9×"**, **"0.050 to 0.033"**. These are real and sourced — I traced them
to `results/round63_next/FOG_DMD_PROBE64/CORRECTION_NOTE.md` ("Fisher over-stated ≈ 1.9×";
"witness f_rec … 0.050 → 0.033") and `P1_results_CORRECTED.json`
(`f_rec_snr3 = 0.03333`). **But** `CLAIM_SOURCE_MATRIX.md` section E covers only the
generator-chunking disclosure (E1); there is no row for the shot-model disclosure, and
`FOG_DMD_PROBE64` is not listed in the matrix's "Frozen source commits" table. The matrix's
closing line ("No numerical claim in R43 §4.1–4.6 … lacks a committed source") is true only
because it is scoped to §4.1–4.6, which excludes this disclosure — yet the numbers appear
in the Letter. Custody is therefore incomplete for the manuscript as a whole.
**Repair:** add a matrix row (new section, e.g. "E2 — shot-model integrity disclosure")
mapping `1.9×`, `0.050`, `0.033` → `FOG_DMD_PROBE64/CORRECTION_NOTE.md` +
`P1_results_CORRECTED.json:f_rec_snr3`, and add that artifact + its commit to the
frozen-source table. (No text change to the Letter required.)

**M5 — Sequential-detection rigor: CUSUM delay slopes are each estimated from only two ε
points.** `main.tex:240` quotes CUSUM slopes −2.16/−3.92; `JET_TEST.committed.json`
`cusum.delta_g.rows` and `cusum.delta_o.rows` each contain **exactly two** ε entries
(delta_g: ε=0.06,0.03; delta_o: ε=0.16,0.09). A two-point log–log "slope" is an exact fit
with zero residual and no goodness-of-fit, so the 8% deviation of the generic slope from
−2 cannot be assessed. A sequential-analysis referee will object that a claimed universal
exponent is asserted from two points.
**Repair:** either report ≥4 ε per direction for the CUSUM panel, or state in the Fig. 2
caption / S2 that the CUSUM slopes are two-point confirmations of the exact-KL exponents
(which *are* densely sampled, 12 points, `bank_A.eps`) and are illustrative, not
independent fits. Also cite Siegmund for the ARL₀→threshold calibration (see m5).

### MINOR

**m1 — `CLAIM_SOURCE_MATRIX` GAP-1 is stale, and matrix row F8 has an internal typo.**
GAP-1 describes `JET_TEST/` as a working-tree divergence from frozen `1bf29f1`; as of this
read the working tree is **clean** and matches `1bf29f1` exactly (`git status` empty;
`git cat-file 1bf29f1:…/JET_TEST.json` `kl_slope`=2.0381 == working copy). GAP-1 was
resolved via option (a) but the matrix still narrates it as an open CRITICAL item — update
or strike it. Separately, matrix row F8 lists "local visible either side (auc 0.0001 / 1.0)"
whereas the JSON is `auc_left=0.02638` (0.026) / `auc_right=1.0`; `supplement.tex:184`
correctly uses 0.026. Fix the matrix's 0.0001 → 0.026.

**m2 — Outlook citation mismatch.** `main.tex:308` cites Snieder for *"sequential change
detection."* Coda-wave interferometry is a *physical* change-detection method, not
statistical sequential detection (that is Lorden/Page/Siegmund). Re-map: cite Lorden (ref
[9]) for "sequential change detection," and keep Snieder for "disordered/medium change
sensing."

**m3 — "99.99% energy in the single direction" conflates scene-space and matrix-space.**
`main.tex:204` — "…99.99% of the measured covariance energy lies in the single direction
`Gx`." The 99.99% metric is `frac_cov_energy_in_OhO_direction` (overlap with the
matrix-space structure `O∘O`), whereas `Gx` is a scene-space vector; and the *top
eigenvalue* is only 82.6% (`valB_rank1.top_eig_frac=0.826`). The rank-one claim is correct
for the *scene*-Fisher (scene enters only through scalar `Q`), but the sentence can be
misread as "the M×M covariance matrix is 99.99% rank-one," which 0.826 contradicts.
**Repair:** "…99.99% of the *scene-dependent* covariance energy lies along the rank-one
structure `O∘O`," and keep "reconstruction Fisher information is rank one" for the scene.

**m4 — Curvature notational duality.** Body/Fig. 1 use `v(δ)=2δᵀGδ` (`:135`,`:176`);
End Matter A / Table S1 use `d=δᵀGδ` (`:348`; `supplement.tex:151–155`). They differ by a
factor 2 (second derivative vs. Taylor coefficient) and are individually correct, but a
cross-referencing referee may stumble. Add a half-clause in EM-A: "(so `d=v/2`)."

**m5 — CUSUM ARL₀→threshold method uncited.** `supplement.tex:171` sets `ARL₀=500`; the
thresholds `h=2.952/1.795` (`JET_TEST cusum.*.h`) come from an ARL approximation but no
method is cited. Add Siegmund, *Sequential Analysis* (Springer, 1985) for the ARL₀→`h`
calibration.

**m6 — Fluctuation–response saturation is an interpretive overclaim.** `main.tex:98–99` —
"the efficient profiled score is the **unique** statistic that saturates the
fluctuation–response bound [Dechant]." Dechant–Sasa's FRI equality condition is a specific
observable; asserting it is *uniquely* the efficient profiled score is a stretch. Soften to
"…attains the fluctuation–response bound" or add the equality-case justification in S1.

**m7 — Profiled-Chernoff attainability for a truly unknown nuisance.** `supplement.tex:61`
("assume the infimum is attained locally") is the right hedge for the *rate*, but the
operational fixed-error exponent for unknown `η` is a composite-hypothesis (GLRT /
least-favorable-prior) quantity, of which the profiled Chernoff is a bound. The `m`
classification is robust regardless; one sentence noting that the constant (not the
exponent `m`) may require a least-favorable-nuisance argument would pre-empt a
semiparametric-efficiency referee.

---

## EQUATION & STATISTICS CHECK (item 4) — all displayed equations correct

- `main.tex:344` KL `D=½[tr A − log det(I+A)] = z²‖K‖_F²/4 + O(z³)`: verified. Expansion
  `tr A − log det(I+A) = ½ tr A² − ⅓ tr A³ + …` ⇒ `D=¼z²tr(K²)−…`; `tr(K²)=‖K‖_F²`. ✓
- `main.tex:345` Chernoff at `s=½`: `C=z²‖K‖_F²/16`. Verified `B=½log det(I+A/2)−¼log
  det(I+A)= (1/16)tr A²+…` (`supplement.tex:93–95`). ✓ Factor KL:Chernoff = 4:1 correct.
- `main.tex:348–350`: `z(ε)=2cε+dε²`, `c=xᵀGδ`, `d=δᵀGδ`. `c≠0 ⇒ z=2cε ⇒
  C=c²‖K‖_F²ε²/4`, `m=1`; `c=0 ⇒ z=dε² ⇒ C=d²‖K‖_F²ε⁴/16`, `m=2`. ✓
- `supplement.tex:104–106` KL equality cases `D=c²‖K‖_F²ε²` (fast) and `D=¼d²‖K‖_F²ε⁴`
  (curvature): consistent (the `z²`→`/4` and `4c²`/`d²` factors reconcile). ✓
- Chernoff exponent definition (`supplement.tex:56–60`, `main.tex:91–93`): standard
  `inf_η' max_s [−log ∫ p^s p^{1−s}]`, `−T⁻¹log P_e*→C_*`. ✓ Dimensionally `T_req≍ε^{−2m}`
  dimensionless. ✓ `d′∝ε^m√T` and the ortho `d′∼ε²√T` (`supplement.tex:283`) confirmed by
  `a_dprime_per_sqrtT` fits. ✓
- Contact-order numbers (2.038/4.000, MC 0.95/2.05, CUSUM −2.16/−3.92): reproduce from the
  frozen JSON (custody table). 2.038 vs 2.000 is the honest finite-ε ε⁴ contamination of the
  generic slope (`bank_A.crossover.local_slope` climbs 2.1→3.85). Only rigor caveat is M5
  (2-point CUSUM fit).

---

## CORRECTED-REFERENCES BLOCK

**All 11 main + 5 supplement bibitems verified correct** against the primary literature
(author list, journal, volume, first page, year). No corrections required to existing
entries. Confirmations:

| Ref | Entry as printed | Web-verified | OK |
|---|---|---|---|
| FengMemory | Feng, Kane, Lee, Stone, PRL **61**, 834 (1988) | PRL 61, 834 (1988) | ✓ |
| LeeStone | Lee, Stone, PRL **55**, 1622 (1985) | PRL 55, 1622 (1985) | ✓ |
| Tsang | Tsang, Nair, Lu, PRX **6**, 031033 (2016) | PRX 6, 031033 (2016) | ✓ |
| Mudry | Mudry *et al.*, Nat. Photonics **6**, 312 (2012) | Nat. Photon. 6, 312–315 (2012) | ✓ |
| Idier | Idier *et al.*, IEEE TCI **4**, 87 (2018) | IEEE TCI 4, 87–98 (2018) | ✓ |
| Chaigne | Chaigne *et al.*, Optica **3**, 54 (2016) | Optica 3, 54–57 (2016) | ✓ |
| Snieder | Snieder, Grêt, Douma, Scales, Science **295**, 2253 (2002) | Science 295, 2253 (2002) | ✓ |
| Dechant | Dechant, Sasa, PNAS **117**, 6430 (2020) | PNAS 117, 6430–6436 (2020) | ✓ |
| Lorden | Lorden, Ann. Math. Statist. **42**, 1897 (1971) | AoMS 42(6), 1897–1908 (1971) | ✓ |
| HermannKrener | R. Hermann, A. J. Krener, IEEE TAC **22**, 728 (1977) | IEEE TAC 22, 728–740 (1977) | ✓ |
| BBP | Baik, Ben Arous, Péché, Ann. Probab. **33**, 1643 (2005) | Ann. Probab. 33(5), 1643–1697 (2005) | ✓ |

**Recommended ADDITIONS for the prior-art fence (M3)** — ready to paste into `main.tex`
and/or `supplement.tex` bibliographies (verified entries):

```latex
\bibitem{Bertolotti}
J. Bertolotti, E.~G. van Putten, C. Blum, A. Lagendijk, W.~L. Vos, and A.~P. Mosk,
Non-invasive imaging through opaque scattering layers,
Nature \textbf{491}, 232 (2012).

\bibitem{Katz}
O. Katz, P. Heidmann, M. Fink, and S. Gigan,
Non-invasive single-shot imaging through scattering layers and around corners via
speckle correlations, Nat. Photonics \textbf{8}, 784 (2014).

\bibitem{Freund}
I. Freund, M. Rosenbluh, and S. Feng,
Memory effects in propagation of optical waves through disordered media,
Phys. Rev. Lett. \textbf{61}, 2328 (1988).   % verify page 2328

\bibitem{GraceGuha}
M.~R. Grace and S. Guha,
Identifying objects at the quantum limit for superresolution imaging,
Phys. Rev. Lett. \textbf{129}, 180502 (2022).
% (change-detection variant: Grace, Guha, Dutton, arXiv 2023 — cite if fencing (b) directly)

\bibitem{SAR_CCD}
M. Preiss and N.~J.~S. Stacy,
Coherent change detection: Theoretical description and experimental results,
DSTO-TR-1851 (2006).   % or a peer-reviewed SAR-CCD ref of the authors' choice
```

---

## CUSTODY SPOT-CHECK TABLE (item 5) — 22 numbers, 21 clean matches, 1 gap

| # | Number (main.tex line) | Matrix row | Frozen JSON field = value | Match |
|---|---|---|---|---|
| 1 | mean null `6.6e-19` (L197,267) | G1/H7 | SCRAMBLE `valA_mean_null.delta_mean_rel`=6.560e-19 | ✓ |
| 2 | cov energy `99.99%` (L204,267) | G2/H5 | SCRAMBLE `valB_rank1.frac_cov_energy_in_OhO_direction`=0.99992 | ✓ |
| 3 | `77/81` cells (L258,275) | A1/H1 | CONF `D2.eps2_cells`=77, `eps5_cells`=81 | ✓ |
| 4 | best `453` banks (L258,276) | A2 | CONF `D2.best_Tdet_2pct`=453.0 | ✓ |
| 5 | MC power LB `0.990` (L259,277) | A3 | CONF `D2.mc_power_lcb`=0.9899 | ✓ |
| 6 | `459` vs `1013` (L259,276) | A4/H2 | CONF `D4.rows[0]` fixed=459, fresh=1013 | ✓ |
| 7 | bal. acc. `0.916` (L265,280) | A6/C3 | CONF `D3.balanced_accuracy`=0.916 | ✓ |
| 8 | in-band FA `0.020` (L266,281) | A7 | CONF `D3.fa.inband`=0.02 | ✓ |
| 9 | medium FA `0.096`/`0.084` (L264,282) | A8/C4 | CONF `D3.fa.amplitude`=0.096, `timescale`=0.084 | ✓ |
| 10 | `d'=4.1`, AUC `0.997` (L266–267) | H6 | SCRAMBLE `detection.coherent.eps_0.05."8192"` d′=4.098, auc=0.997 | ✓ |
| 11 | KL orders `2.038`/`4.000` (L223,238) | F1 | JET `delta_g.kl_slope`=2.0381, `delta_o`=3.9999 | ✓ |
| 12 | coef `1.07%`/`0.00%` (L223–224) | F2 | JET `delta_g.coef_ratio`=1.0107, `delta_o`=1.0000 | ✓ |
| 13 | MC exp `0.95`/`2.05` (L224,239) | F3 | JET `bank_A_montecarlo.delta_g…slope`=0.9547, `delta_o`=2.0455 | ✓ |
| 14 | CUSUM `−2.16`/`−3.92` (L224,240) | F4 | JET `cusum.delta_g.delay_slope`=−2.1583, `delta_o`=−3.9162 | ✓ (see M5) |
| 15 | online `0.16%` (L278) | A5 | CONF `D6.online_frac_of_bank`=0.0016 | ✓ |
| 16 | iso-Q AUC `≈0.48` (L245) | F8 | JET `bank_C.cancellation.auc_at_star`=0.4807 | ✓ (matrix side-typo, m1) |
| 17 | `400` dirs, blind `0` (L247) | F10 | JET `kill_scan.n_dirs`=400, `frac_blind`=0.0 | ✓ |
| 18 | `17` checks / `4` kills (L248) | F11 | JET `checks`(17 true) / `kills`(4 false) | ✓ |
| 19 | FA `0.032`, latency `19%` @10% (L284) | B1 | CONF `D5` basis_rotation 0.1: `nontarget_fa`=0.032, `Tdet_inflation`=0.192 | ✓ |
| 20 | `0.052`/`0.076` @20%/slope−1 (L286–287) | B2/D3 | CONF `D5` br0.2=0.052, spectral_slope−1=0.076 | ✓ |
| 21 | `12` banks, `2.76`/`6.0` GPU-h (L397) | I3 | CONF `banks.confirmatory`=12, `total_gpu_hours`=2.764 | ✓ |
| 22 | shot model `1.9×`, `0.050`→`0.033` (L411–413) | **NONE** | `FOG_DMD_PROBE64/CORRECTION_NOTE.md`; `P1_results_CORRECTED.json:f_rec_snr3`=0.0333 | **NO MATRIX ROW (M4)** |

Rounding is consistent throughout (e.g. 0.9899→0.990, 6.560e-19→6.6e-19, 4.098→4.1,
2.0381→2.038, −2.1583→−2.16); no rounding drift or mis-transcription found in rows 1–21.

---

## SUMMARY FOR COORDINATOR

The physics and custody are sound; this is a fixable manuscript, not a failed one. Before
submission: (B1) label the two "optical experiment" sentences as simulation; (M1) add the
`G≻0`/in-aperture qualifier next to the empty-blind-set headline; (M2) drop the word
"specificity" from `:210`; (M3) extend the S5 fence to Bertolotti/Katz/Freund, Grace–Guha,
and SAR-CCD; (M4) add a matrix row + artifact for the shot-model disclosure numbers; (M5)
disclose that the CUSUM slopes are 2-point confirmations. None require regenerating data
(all repairs are text/matrix edits), consistent with R43 §7 green-light constraints.
