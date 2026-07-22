# paper3_scatter — scaffolding notes (Part II)

Companion **Part II** to the existing identifiability manuscript (**Part I**).
Scaffolded from the **NORMATIVE** R26 ruling,
`docs/ROUND63_GPT_ROUND26_RULING_RAW.md`. Spine follows ruling §24 (eight
sections). This file records, per block, what is **verbatim-from-ruling**, what
is **TODO**, and what is an **open lemma** — plus the numerical-verification
plan sketch.

Working title (ruling §24): *Information limits and temporal design for ghost
imaging through dynamic scattering media*. O2 does **not** appear in the title:
the ruling forbids it until a concrete resource-conditioned law is derived.

## Build

```
cd paper3_scatter
pdflatex main_scatter.tex && bibtex main_scatter && pdflatex main_scatter.tex && pdflatex main_scatter.tex
```

Status at scaffold time: all passes exit 0; **zero** undefined references, zero
undefined citations, zero overfull hboxes; **11 pages**.

## Deliberate deviation from paper2 conventions

paper2/main_m1.tex is `[10pt,twocolumn]`. This file is single-column `[11pt]`.
Reason: the O3/O4 matrix and Schur-complement displays (O3.5, O3.6, O4.1, O4.2)
are too wide for a two-column measure and would produce overfull boxes. All
other conventions are inherited from paper2: `amsmath/amssymb/bm`, `geometry`,
`xcolor`, `hyperref`, `placeins`, the `\SPH{...}` red-placeholder macro,
`\newtheorem`, and `\bibliographystyle{unsrt}` + `\bibliography{refs}`. Added
here: `amsthm` (for `proof` and the extra theorem-like environments) and a
`\TODO{...}` blue macro.

## The fence (ruling §23)

The seven claims the program must **not** make are reproduced as a comment block
at the top of `main_scatter.tex`. Do not let any of them re-enter the abstract,
title, or body.

## Section-by-section provenance

Equation tags below (O1.1, O1.2, …) are the ruling's own tags, reproduced via
`\tag{...}` so cross-references print the ruling numbering.

| Spine § | File section | Provenance |
|---|---|---|
| Front matter | title / author / abstract | Title **verbatim** from §24. Author = `\SPH` placeholder. Abstract = **STUB**, newly written (5 sentences, honest, no O2 ridge claim). Part I cited via `companion` placeholder. |
| — | Introduction | New scaffolding prose; the clean hierarchy line is **verbatim** from §23; framing of Part I vs Part II is **verbatim-in-substance** from §24–25. |
| 1 | Physical hypothesis MG | (MG) box **verbatim** from §21; the six sufficient conditions **verbatim** from §21; the break list **verbatim** from §22; Remark (scope + gauge singularity) condenses §12 + §22. |
| 2 | O1 missing-information theorem | Statement O1.1 / O1.2, the four assumptions, the proof, and the Louis-attribution sentence ("This proof is classical Louis theory. The following programmable-GI specialization is the important part.") all **verbatim** from §9. |
| 3 | O1-P + Gaussian readout | O1.3 / O1.4 / O1.5, the proof, the gain-weighted illumination-vector discussion, and the three consequences **verbatim** from §10. Gaussian readout counterpart **verbatim** from §11. |
| 4 | O3 correction ceiling (van Trees) | O3.1 / O3.2 / O3.3 **verbatim** from §13; may-claim / may-not-claim lists **verbatim** from §13 (as scoping Remarks). §14 connection: O3.4 / O3.5 / O3.6 **verbatim** from §14. **The `vB_L` scalarization is a TODO** — the ruling only sketches it (see below). |
| 5 | O4-A nuisance-orthogonal design | O4.1 model + O4.2 Schur complement **verbatim** from §17; O4-A statement, O4.3, O4.4, and the proof **verbatim** from §18; the chopping/interleaving/randomization meaning paragraph **verbatim** from §18. |
| 6 | Randomization + flip boundary | Randomization proposition + concentration bounds **verbatim** from §19; the **named obstacle is carried as an explicit open-lemma Remark** (see below). Flip-boundary rederivation O4.5 **verbatim** from §20. |
| 7 | Model-specific operating surface | **STUB** — one paragraph + one `\begin{conjecture}` (O2). Resource preconditions (r, η, R, ω, N/K, read-noise ratio) **verbatim** from §16; the named obstacle **verbatim** from §16. Labeled FUTURE/conditional; **no ridge claim**. |
| 8 | Reach: gain drift + dead time | O5 coupling discussion and the four consequence bullets **verbatim** from §5 (the "O5" section); the explicit **no-scalar-product-law** statement **verbatim** from that section. |

### Minor transcription fix (not a content change)

- Ruling (O3.5) renders as `+rac12\operatorname{tr}(...)` — a dropped-backslash
  artifact of the source markdown. Transcribed as the intended
  `+\tfrac12\operatorname{tr}(...)`. No mathematical change.

## TODO items (explicitly marked in the source)

1. **`vB_L` derivation (§4, `\TODO` + Remark `rem:o3-todo`).** The ruling §14
   only *sketches* that a scalarized/relMSE expansion of (O3.6) produces a
   `v B_L` term matching Part I's static-correction ceiling. Task: fix the task
   scalarization `W`, expand `tr{W [I_pi + E I_Y]^{-1}}` from (O3.6) to first
   order in `v`, and identify `B_L` with Part I's coefficient. This is a theorem
   path, not an automatic identity between manuscripts.
2. **Front-matter `\SPH` placeholders:** author block, repo URL, funding/ack.
3. **Companion (Part I) citation:** `companion` entry title/author are
   placeholders pending Part I's final front matter.

## Open lemma (must remain until closed)

- **Randomization proposition, §6, Remark `rem:rand-obstacle`.**
  Proposition `prop:rand` is a theorem only under bounded/centered
  permutation-sampling assumptions. For arbitrary nonnegative GI patterns the
  scene-dependent carrier `s_n m_n` can destroy centering. **Closing the lemma
  requires the stationary-carrier hypothesis of the companion Part I.** The
  Remark states this explicitly and must not be deleted until the lemma is
  proved.

## Bibliography provenance (`refs.bib`)

Follows paper2/refs.bib discipline (`%UNVERIFIED` / `VERIFY at submission`).

- **DOIs transcribed verbatim from the ruling's hyperlinks** (Q1 §§1–8):
  `xiao2022temporal`, `pengchen2023`, `zhou2023dual`, `pengchen2024apl`,
  `hao2025photonlimited`, `louis1982`, `tichavsky1998`, `lilee2015`,
  `lingstrohmer2018`, `cambareri2019`, `kechkrahmer2017`, `fixsen2000calib`,
  `arendt2000dither`, `uecker2008`, `cidre2015`, `basic2017`, `pine1988`,
  `bandyopadhyay2005`, `yuan2005`, `scirep2023speckle`, `hassibi2003`,
  `qureshi_ng`, `mehra1974`. For these, the **DOI/arXiv ID is from the ruling**;
  author/volume/page metadata beyond what the ruling names is marked
  `VERIFY at submission` and needs a live Crossref pass before submission.
- **Task-named, NOT in the ruling raw** (seeded as standard OED references, all
  fields `VERIFY at submission`): `vantrees` (van Trees), `dette1997`,
  `imhofwong2000`, `pukelsheim2006`.
- **Reused verbatim from the Crossref-verified paper2/refs.bib**:
  `duarte2008singlepixel`, `edgar2019principles`, `shapiro2008computational`,
  `hao2025photonlimited`.
- `\nocite{*}` is used so every seed entry appears in the bibliography.

## Numerical-verification plan sketch (DO NOT RUN — plan only)

Goal: independently confirm the two load-bearing linear-algebra claims of §5
before they anchor any design recommendation. Both are exact finite-dimensional
identities, so verification is deterministic (no Monte Carlo needed for the core
check; MC only for the O1-P posterior-covariance sanity check).

1. **O4.2 Schur-complement formula.**
   - Draw small random `M_pi` (N×K), diagonal `R`, drift basis `H` (N×p),
     prior precision `Lambda_beta`, object `x`; form `s_pi = M_pi x`,
     `D_{s_pi}`.
   - Assemble the joint information `I(pi)` block-by-block from the §17
     definitions, and independently compute the object sub-block by the Schur
     complement (O4.2).
   - Check `I_{x|gain}` equals `(I(pi)^{-1})_{xx-block}^{-1}` (Schur identity)
     to machine precision across many random seeds and a sweep of `p`, `N`, `K`.

2. **O4.4 orthogonality conditions.**
   - Construct schedules that *satisfy* (O4.3): pick `M_pi`, `H`, `R` with
     `M_pi^T R^{-1} D_{s_pi} H = 0` by projecting `H` onto the null space of
     `(D_{s_pi} R^{-1} M_pi)`; verify column-wise moment sums (O4.4) vanish and
     that `I_{x|gain} = I_{xx}` exactly.
   - Construct schedules that *violate* it and confirm `I_{x|gain} \prec I_{xx}`
     (strict, by checking the smallest eigenvalue of `I_{xx} - I_{x|gain} > 0`).

3. **O1-P posterior-covariance sanity check (OU gain + schedules).**
   - Simulate an OU log-gain path `a_n = exp(l_n)` on a schedule grid; generate
     Poisson bucket counts `Y_n ~ Poisson(a_n s_n)`.
   - Estimate `Cov(a | Y, x)` (small K, e.g. particle/grid posterior) and check
     the exact hidden-exposure identity (O1.4):
     `I_Y = M^T diag(E a_n / s_n) M - M^T E_Y[Cov(a|Y,x)] M`, comparing against
     a direct finite-difference Fisher of the marginal log-likelihood
     `log p(Y|x)`.
   - Sweep OU correlation time `t_c` and pattern order to confirm consequence
     (1)–(2) of §3: the loss depends on the *full* posterior covariance and on
     pattern ordering, not on a scalar variance.

Deliverables when run: a single script + a table of max abs errors per identity;
target < 1e-10 for the exact linear-algebra checks (1,2) and < a few % Monte
Carlo agreement for (3). **Not executed at scaffold time.**


## Numerical verification results (2026-07-22, results/scatter_verify/)

All three blocks PASS (O1.4 identity via exact HMM posterior, paired-CRN
max z=2.19 vs FD Fisher; O4.2 Schur 9.1e-16; exact-zero chop + delta^2
perturbation slope 1.984; ordered 1770 >> random 63.6 > paired 0.0;
concentration slope -0.497; obstacle demo: uncentered carrier bias flat
~31-33 under permutation). Two findings to fold into the manuscript:

1. **DC-mode gauge remark (add near Theorem O4-A):** the constant drift
   mode is NOT orthogonalisable by any schedule — its moment equals
   sum s_n^2/sigma_n^2 > 0 identically (verified to 16 digits). The DC
   gain mode IS the scale gauge, degenerate with object amplitude; O4.4
   is satisfiable only for centered/AC modes. This is Part I's gauge
   invariance (a,T)->(ca,T/c) seen from the design side — a genuine
   Part I/Part II interlock sentence.
2. **Methods note (verification section):** compute the missing term as
   the small posterior covariance Cov(G|Y) DIRECTLY; the difference form
   E[aa^T]-E[a|Y-hat outer] suffers catastrophic cancellation along the
   gauge direction (factor-30 error observed).
