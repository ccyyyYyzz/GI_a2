# Fig. 1 — mechanism / schematic figure for `paper/main_oe.tex`

**Files**
- Figure: `paper/figs/fig_mechanism_p1.pdf` (vector; `.png` at 220 dpi for preview)
- Generator: `code/round63/figs/fig_mechanism_p1.py`
  (self-contained; numpy + matplotlib only; reads no campaign artifact)
- Regenerate (cwd = repo root `D:\GI_another`):
  `python code\round63\figs\fig_mechanism_p1.py`

**Placement.** Designed for the full two-column text width, i.e. a `figure*`
float (the same width as `fig:study1` and `fig:ladder`). Base font 8 pt; the
Okabe–Ito/Wong colourblind-safe palette matches the other round-63 figures.

**What it shows (three connected panels).**
- **(a)** the single-pixel imaging chain — light source → DMD (binary
  illumination pattern, inset) → object (low- vs high-contrast pair) → SPD →
  RQL reconstruction, with the forward model along the bottom;
- **(b)** the dead-time mechanism — a photon-arrival timeline with the blind
  interval τ censoring counts (× = lost), and the mean detected-vs-incident
  saturation ρ/(1+ρ) bending away from the ideal, with the count-information
  ridge ρ*(ν) marked;
- **(c)** the contrast–dead-time operating map — bucket contrast C_u vs
  detector load ρ, the photon-limited onset Γ = 1, and the six Study-2
  families placed by measured contrast and coloured by the frozen gate outcome.

---

## Ready-to-paste LaTeX

Add near the top of `main_oe.tex` (e.g. just after the Introduction). **This
block is provided for you to paste; no `.tex` file was modified.**

```latex
\begin{figure*}[t]
  \centering
  \includegraphics[width=\textwidth]{figs/fig_mechanism_p1.pdf}
  \caption{Contrast--dead-time operating map for high-flux single-pixel
  imaging. DMD: digital micromirror device; SPD: single-pixel (bucket)
  photon-counting detector; RQL: renewal quasi-likelihood reconstruction;
  $a_i$: binary illumination pattern; $x$: scene; $N_i$: integrated bucket
  count; $\hat{x}$: reconstruction; $\Phi$: incident-flux scale; $d$: dark
  rate; $\tau$: detector dead time; $\rho=\lambda\tau$: per-slot detector
  load; $\nu=T/\tau$: dead-time slots per exposure $T$; $\Cu=\mathrm{sd}_i(u_i)
  /\mathrm{mean}_i(u_i)$: measured bucket-energy contrast of $u_i=a_i^{\top}x$;
  $\Gamma=\Cu\sqrt{\nu\rho/(1+\rho)}$: contrast-to-count-noise index.
  (a)~Single-pixel imaging chain: structured illumination from a DMD probes a
  low- or high-contrast scene $x$, and the bucket counts $N_i$ feed the RQL
  reconstruction $\hat{x}$ under the forward model
  $\lambda_i=\Phi\,a_i^{\top}x+d$, exposure $T=\nu\tau$.
  (b)~Dead-time mechanism: after each registered count the detector is blind
  for $\tau$ (shaded), censoring later arrivals ($\times$); the mean detected
  counts per slot $\nu_{\det}/\nu=\rho/(1+\rho)$ bend away from the ideal
  (dashed) as the load rises, while the exact count-information peaks on the
  ridge $\rho^{*}(\nu)=(6\nu)^{1/3}-2/3$ ($\rho^{*}\approx7.8$ and $22.2$ at
  $\nu=100$ and $2000$, respectively), well beyond the conventional
  $\rho=0.05$ and the high-flux $\rhobar=0.6$ operating points.
  (c)~Operating map in the contrast--load plane: the photon-limited onset
  $\Gamma=1$ is a descriptive photon-limited-onset proxy, marking bucket
  modulation comparable to count noise; it is not a predicted or validated
  boundary for high-flux image benefit. The six Study-2
  structural families are placed at the high-flux operating load
  $\rhobar=0.6$ by their measured bucket contrast and coloured by the frozen
  acquisition-speed gate outcome, with the number of accelerating instances
  shown as a green wedge: maze ($\Cu=0.53$), glyph ($0.41$) and spokes
  ($0.38$) accelerate $4/4$, chirp ($0.28$) $3/4$, contour ($0.27$) $1/4$, and
  microtexture ($0.17$) $0/4$, respectively, while the dense $k=512$ rung
  ($\Cu=0.040$) sits at the $\Gamma=1$ multiplex-limited boundary. All six
  family markers share the single operating load $\rhobar=0.6$; the
  horizontal offsets are cosmetic, applied only to separate near-coincident
  markers (chirp and contour differ by $0.01$ in $\Cu$), and thin leader
  lines tie each family label to its marker. The contrast-controlled
  transition band is descriptive and does not estimate a universal critical
  contrast. Grayscale object icons are for display purposes.}
  \label{fig:mechanism}
\end{figure*}
```

The caption uses the `\Cu` and `\rhobar` macros already defined in
`main_oe.tex` (lines 20–21).

---

## Numbers pulled from `main_oe.tex` (no value invented)

| quantity | value(s) | source (main_oe.tex) |
|---|---|---|
| conventional load | ρ ≲ 0.05 | L98, L167 |
| high-flux operating point | ρ̄ = 0.6 | L368, L479 |
| exposure / slots | T = ντ, ν = T/τ | L51, L158 |
| forward model | λ_i = Φ a_iᵀx + d | L156–157, L189 |
| cube-root ridge | ρ*(ν) = (6ν)^{1/3} − 2/3 | Eq. (6), L265 |
| ridge values | ρ* ≈ 7.8 (ν=100), 22.2 (ν=2000) | closed form, cf. L275–278 |
| detected fraction (non-paralyzable mean) | ρ/(1+ρ) | ridge kernel `J_CLT`, `fig_a_ridge_map.py` L87 |
| contrast-to-count-noise index | Γ = C_u √(νρ/(1+ρ)); onset Γ=1 | L465–466 |
| family contrast & gate | maze 0.53 4/4; glyph 0.41 4/4; spokes 0.38 4/4; chirp 0.28 3/4; contour 0.27 1/4; microtexture 0.17 0/4 | L505–508 |
| dense rung contrast | C_u = 0.040 (k=512), Γ = 1.09 | L476–477, L480 |

Quantities **deliberately omitted** because they are not in the manuscript
(per the "if a number is not in the tex, omit it" rule): the raster-rung
contrast 1.270 and per-occupancy Γ (6.07/8.66/34.79) are Study-2 ladder
numbers not needed by this schematic; exact-PMF ρ* values (4.53 … 22.16) are
replaced by the closed-form law the panel actually draws.

## Housekeeping / honesty notes baked into the design
- The family markers all sit at the single operating load ρ̄ = 0.6 (the
  Study-2 high-flux point); the horizontal spread is cosmetic, only to keep
  the six pass-fraction pies from overlapping (chirp/contour are 0.01 apart
  in C_u and need the largest offsets), and thin leader lines run from each
  family label to the centre of its pie — both stated in the caption.
- The contrast-controlled transition is drawn as a *band* and labelled
  descriptive, matching the manuscript's "consistent with — but not
  establishing — a contrast-controlled transition" and "does not estimate a
  universal critical contrast" (L71–73, L508–510). No "predicted boundary"
  wording is used.
- Γ = 1 is the manuscript's photon-limited *onset*, not a pass/fail gate.
