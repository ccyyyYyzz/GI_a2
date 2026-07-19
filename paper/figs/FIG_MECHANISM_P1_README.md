# Paper-1 figures — R20 presentation update (GitHub issue #12)

R20 (M11/M12) reorganized the paper-1 (`paper/main_oe.tex`) figures. This README
gives the new figure sequence, the ready-to-paste journal captions, and the
provenance. **No `.tex` file was modified**; the caption blocks are for the text
agent to install.

**Style (M10 / SF7 / SF8).** All owned paper-1 generators were converted to a
clean sans-serif family; the reduced palette is neutral gray = safe/baseline,
one saturated blue = ridge/high-flux, one secondary vermillion = dead-time /
negative. Final absolute sizes and single/double-column widths are approximate
(built at 8.6 cm / 17.8 cm guesses); a resize + font-match pass follows the
Optica template port (M9).

## New figure sequence

| # | file(s) | content | column |
|---|---|---|---|
| Fig. 1 | `fig_mechanism_p1.pdf` | SPI chain + dead-time/ridge + **neutral** $\Gamma=1$ schematic (M12) | double |
| Fig. 2 | `fig_a_ridge_map.pdf` + `fig_b_masks.pdf` + `fig_c_ladder.pdf` | mechanism / pattern ladder (M11) | double |
| Fig. 3 | `fig_d_dwell.pdf` + `fig_e_gate.pdf` | image outcomes (M11) | double |
| results | `fig_p1_families.pdf` | six-family empirical outcomes (moved out of Fig. 1) (M12) | single |
| suppl. | `fig_f_scattering.pdf` | dynamic scattering (moved out of the ladder) (M11) | single |

The old six-panel `fig:ladder` (ridge + masks + $C_u/\Gamma$ + dwell +
scattering + gate) is **split**: mechanism/pattern panels → Fig. 2, image
outcomes → Fig. 3, dynamic scattering → supplement. Compose Fig. 2 and Fig. 3
as `figure*` floats with three / two `minipage` panels as before.

Regenerate (cwd = repo root `D:\GI_another`):

```
python code/round63/figs/fig_mechanism_p1.py     # Fig. 1
python code/round63/figs/fig_p1_families.py       # results
python code/round63/figs/fig_b_masks.py fig_c_ladder.py fig_d_dwell.py fig_e_gate.py fig_f_scattering.py
```

(`fig_a_ridge_map` lives under `results/round63_theory/` and is **not** owned
here — see template-port flags.)

---

## Ready-to-paste journal captions

### Fig. 1 — `fig_mechanism_p1.pdf` (double column, ~150 words)

> **Contrast–dead-time operating map for high-flux single-pixel imaging.**
> (a) Imaging chain: a DMD imposes a binary pattern $a_i$ on a low- or
> high-contrast scene $x$, and the bucket counts $N_i$ feed the renewal
> quasi-likelihood (RQL) reconstruction $\hat{x}$ under
> $\lambda_i=\Phi\,a_i^{\top}x+d$, exposure $T=\nu\tau$.
> (b) Dead-time mechanism: after each registered count the detector is blind for
> $\tau$ (shaded), censoring later arrivals ($\times$); the mean detected counts
> per slot $\rho/(1+\rho)$ bend away from the ideal, while the exact
> count-information peaks on the ridge $\rho^{*}(\nu)=(6\nu)^{1/3}-2/3$
> ($\rho^{*}\approx7.8$ and $22.2$ at $\nu=100$ and $2000$), far beyond the
> conventional $\rho=0.05$ and high-flux $\bar\rho=0.6$ points.
> (c) The photon-limited onset $\Gamma=1$ in the contrast–load plane is a
> **descriptive onset proxy**, not a predicted or validated boundary; the
> conventional and high-flux operating loads are marked. Grayscale object icons
> are illustrative.

### Fig. 2 — mechanism / pattern ladder (double column, ~140 words)

> **Exact information ridge and the illumination ladder (Study 2).**
> (a) Exact count-Fisher information per unit time versus detector load,
> peaking on the ridge $\rho^{*}(\nu)=(6\nu)^{1/3}-2/3$.
> (b) Representative SCAT illumination masks across the occupancy ladder
> $k\in\{512,32,16,1\}$; the per-pixel dose is held fixed by the amplitude
> normalization $A=(n/k)\,B$, $\Phi=\bar\rho/\tau$ (equal mean load, equal
> incident dose per pixel).
> (c) Measured bucket-energy contrast $C_u$ (left) and the photon-limited flux
> figure $\Gamma=C_u\sqrt{\nu\rho/(1+\rho)}$ at $\bar\rho=0.6,\ \nu=2000$
> (right) versus occupancy; $C_u$ falls as $\sim1/\sqrt{k}$ and the dense
> $k=512$ rung sits at the descriptive photon-limited onset $\Gamma=1$. Error
> bars are $\pm1$ s.d. over image–seed pairs.

### Fig. 3 — image outcomes (double column, ~135 words)

> **Dwell curves, fixed-dwell gain, and acquisition speed (Study 2).**
> (a) Radiometric PSNR and CNR versus dwell budget $\nu$ for the preregistered
> subject at the safe ($\bar\rho=0.05$) and high-flux ($\bar\rho=0.6$)
> operating points.
> (b) Per-image fixed-dwell quality gain $\Delta$PSNR$_\mathrm{rad}$
> ($\bar\rho=0.6$ vs $0.05$, $\nu=2000$) across the occupancy ladder, with the
> per-rung median (diamond) and the $+1$\,dB preregistered secondary bar.
> (c) The frozen $k=16$ acquisition-speedup gate $S$: per-scene speedups with
> the median (diamond) and the $S=3$ gate; the other rungs carry no gate by the
> frozen ruling, and censored ($S=0$) scenes are drawn as open markers. The
> speedup gate is evaluated at $k=16$ only.

### Results — `fig_p1_families.pdf` (single column, ~70 words)

> **High-flux acceleration by scene family (empirical).** For the six Study-2
> structural families, the number of accelerating instances (out of four) rises
> with measured bucket contrast $C_u$: maze, glyph and spokes accelerate $4/4$,
> chirp $3/4$, contour $1/4$, and microtexture $0/4$. These are empirical
> outcomes at the single operating load $\bar\rho=0.6$; no critical-contrast
> boundary or transition band is claimed.

### Supplement — `fig_f_scattering.pdf` (single column, ~70 words)

> **Dynamic-scattering robustness.** Radiometric PSNR versus measured detected
> signal $S_\mathrm{det}$ under a fluctuating scattering medium
> ($\mathrm{CV}(\alpha)=0.20$) at the safe and high-flux operating points; the
> high-flux point retains its advantage across the experimentally reported
> $S_\mathrm{det}$ range (shaded). Error bars span the dwell replicates.

---

## PROVENANCE — data sources and selection rules

- **`fig_mechanism_p1` (M12).** Self-contained schematic (numpy + matplotlib);
  every constant is transcribed from `paper/main_oe.tex` (loads $\rho=0.05$,
  $\bar\rho=0.6$; ridge $(6\nu)^{1/3}-2/3$; $\Gamma=C_u\sqrt{\nu\rho/(1+\rho)}$).
  **Change vs. R19:** the six-family pass-fraction pies, the "high-flux
  helps" / "multiplex-limited" region labels, the transition band, and the
  cosmetic horizontal offsets were **removed** from panel (c); it is now a
  neutral $\Gamma=1$ schematic labelled *descriptive onset proxy* (R9 language
  preserved). The removed family outcomes moved to `fig_p1_families`.
- **`fig_p1_families` (M12).** Family contrast/gate transcribed from
  `main_oe.tex` L505–508: maze $0.53$ $4/4$; glyph $0.41$ $4/4$; spokes $0.38$
  $4/4$; chirp $0.28$ $3/4$; contour $0.27$ $1/4$; microtexture $0.17$ $0/4$.
  Points at true measured $C_u$ (no offsets); coloured gray→blue by pass
  fraction. No regions, no band, no $\Gamma=1$ curve.
- **Fig. 2 components.** `fig_a_ridge_map` (exact ridge, `results/round63_theory`);
  `fig_b_masks` (illustrative masks + dose normalization, self-contained);
  `fig_c_ladder` reads `results/round63_study2/fluxmap_rows.csv` ($C_u$ mean over
  all rows at $k$; $\Gamma$ at $\bar\rho=0.6,\nu=2000$; $\pm1$ s.d. over
  image–seed pairs).
- **Fig. 3 components.** `fig_d_dwell` (PSNR/CNR vs $\nu$, prereg subject);
  `fig_e_gate` reads `study2_summary.json` (`endpoint_analysis.per_image`
  `DeltaQ_budget`, `S_gate` for $k=16$) and `controls_rows.csv` /
  `robustness_rows.csv` for the $k\in\{512,32,1\}$ fixed-dwell $\Delta Q$ rungs.
  Empty / `nan` cells are dropped, never imputed; horizontal jitter is
  deterministic and perturbs position only.
- **`fig_f_scattering`.** `results/round63_study2/scattering_rows.csv`,
  $\mathrm{CV}(\alpha)=0.20$.

---

## Flags for the Optica template port (M9)

- **`fig_a_ridge_map` is not owned here** (`results/round63_theory/`,
  referenced by `main_oe.tex` as `../results/round63_theory/fig_a_ridge_map.pdf`)
  and is still in its original style. Convert it to the shared sans-serif /
  reduced palette in the template-port pass so Fig. 2 is internally consistent.
- **Resize pass required.** Widths are 8.6 cm / 17.8 cm guesses; rebuild at the
  true single/double-column widths after the port, re-run overfull/underfull
  checks, and repeat a 100%-scale legibility check (min in-panel text ≥ 7.5–8 pt).
- **Supplement generators** (`supp_s*.py`) were not restyled in this pass; apply
  the same sans-serif conversion at port time.
- **`fig_nat_pair`, `fig_s1b/c/d`, `fig_s2_gallery`** were converted to
  sans-serif here (SF7) but their layout/size (M14 — median+worst gallery
  scale, dense-study median-plus-band) is out of the R20 figure-owner scope and
  should be revisited at port.
