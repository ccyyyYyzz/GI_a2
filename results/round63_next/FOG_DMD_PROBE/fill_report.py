# Populate FOG_DMD_PROBE_REPORT.md tables + verdict from the result JSONs.
import json

def f3(x):
    return "—" if x is None else f"{x:.3f}"
def pct(x):
    return "—" if x is None or x != x else f"{x*100:.0f}%"
def phs(p):
    return "clean" if p is None else f"{p:.0e}"

E1 = json.load(open('E1_results.json'))

def prim_table(rows, warm=True):
    hdr = "| cell | χ_full | photons | oracle | lin-approx oracle | **blind cold (median)** | captured | warm-Q ceiling |\n"
    hdr += "|:--|:--:|:--:|:--:|:--:|:--:|:--:|:--:|\n"
    body = ""
    for r in rows:
        body += (f"| M{r['M']} d{r['d']} σ_f{r['sigma_f']} | {r['chi']:.2f} | {phs(r['photons'])} "
                 f"| {f3(r['oracle_nmse'])} | {f3(r['linapprox_oracle_nmse'])} "
                 f"| **{f3(r['blind_cold_nmse_median'])}** | {pct(r['frac_oracle_improvement'])} "
                 f"| {f3(r.get('warm_Q_ceiling_nmse'))} |\n")
    return hdr + body

prim = E1['primary']
linref = E1.get('linear_reference', [])
ou = E1.get('ou_sensitivity', [])

primary_md = "**Lognormal medium (ruling-mandated), physical 0/1 patterns, T=128, 5-seed medians.**\n\n"
primary_md += prim_table(prim)
best = min(prim, key=lambda r: r['blind_cold_nmse_median'])
primary_md += (f"\nBest primary blind cold-start: **NMSE {best['blind_cold_nmse_median']:.3f}** "
               f"at {best['tag']} ph={phs(best['photons'])} "
               f"(captures {pct(best['frac_oracle_improvement'])} of the oracle improvement).\n")

primary_md += "\n### 5a. Model-matched linear-medium reference (isolates cold-start vs mismatch)\n\n"
primary_md += ("Same geometry, affine medium w=1+Uz that the estimator models exactly. The warm-Q "
               "column is now the true **profiling-tax ceiling** (no model mismatch).\n\n")
primary_md += prim_table(linref)
primary_md += ("\nReading: with the model matched, the warm-start ceiling falls to "
               f"{min(r['warm_Q_ceiling_nmse'] for r in linref if r.get('warm_Q_ceiling_nmse')):.3f}"
               "–"
               f"{max(r['warm_Q_ceiling_nmse'] for r in linref if r.get('warm_Q_ceiling_nmse')):.3f}"
               " (vs 0.86–1.0 at M=64) — the capacity margin is real — yet blind **cold-start** stays "
               f"{min(r['blind_cold_nmse_median'] for r in linref):.2f}"
               f"–{max(r['blind_cold_nmse_median'] for r in linref):.2f}. Cold-start is the wall.\n")

primary_md += "\n### 5b. OU-prior sensitivity (medium OU τ=16; prior weakened 4×)\n\n"
primary_md += prim_table(ou, warm=False)
if len(ou) == 2:
    d = ou[1]['blind_cold_nmse_median'] - ou[0]['blind_cold_nmse_median']
    primary_md += (f"\nWeakening the OU precision 4× changes blind null-NMSE by {d:+.3f} "
                   f"({ou[0]['blind_cold_nmse_median']:.3f} → {ou[1]['blind_cold_nmse_median']:.3f}). "
                   "The (small) recovery is **not** prior-imputed — it is likelihood-driven — but it "
                   "remains far above the bar either way.\n")

# E2
try:
    E2 = json.load(open('E2_noise_curve.json'))
    e2_md = (f"Cell M={E2['cell']['M']}, d={E2['cell']['d']}, T={E2['cell']['T']}, "
             f"σ_f={E2['cell']['sigma_f']}, lognormal; {E2['cell']['seeds']}-seed medians. "
             "null-NMSE vs mean detected photons/bucket (baseline fixed-averaging = 1.0):\n\n")
    e2_md += "| photons/bucket | oracle NMSE | blind cold NMSE | baseline |\n|:--:|:--:|:--:|:--:|\n"
    for r in E2['curve']:
        e2_md += f"| {r['photons']:.0e} | {r['oracle_nmse_median']:.3f} | {r['blind_nmse_median']:.3f} | 1.000 |\n"
    e2_md += ("\nThe oracle degrades **gracefully** toward 1.0 as photons fall (regularized, no explosion): "
              f"{E2['curve'][-1]['oracle_nmse_median']:.3f} at 1e7 → "
              f"{E2['curve'][0]['oracle_nmse_median']:.3f} at 1e3. The blind arm never approaches the F3 "
              "bar at any dose. See `E2_noise_curve.png`.\n")
except FileNotFoundError:
    e2_md = "_E2_noise_curve.json not found._\n"

# Verdict
passF3 = any((r['blind_cold_nmse_median'] <= 0.25 and (r['frac_oracle_improvement'] or 0) >= 0.5) for r in prim)
verdict = f"""**E1 (blind cold-start recovery) — SUCCESS SIGNAL: NOT MET.**
Best primary blind cold-start null-NMSE = **{best['blind_cold_nmse_median']:.3f}** (at {best['tag']},
ph={phs(best['photons'])}), capturing **{pct(best['frac_oracle_improvement'])}** of the oracle
improvement. The F3-shadow bar requires ≤0.25 **and** ≥50% captured; no cell meets it
(F3 pass = {passF3}). Across the whole M=96 grid, blind cold-start sits at NMSE ≈ 0.54–1.0.

- **Oracle viability holds** (null-NMSE ≈ 0 clean; ≤0.02 at ≤1e5 photons at d≤48): the data *do*
  contain the null information, and the M=96 capacity margin (χ≈2.0–2.3) is real — the
  *model-matched warm-start ceiling* drops to 0.20–0.26 (vs 0.86–1.0 at M=64).
- **The wall is the bilinear cold-start optimizer.** Even where the warm ceiling reaches the bar
  (d=64, σ_f=0.15: warm-Q 0.203), the data-only staged cold start stalls at 0.58. This is R38's
  predicted kill-tree node: *F1/F2-type geometry OK, F3 cold-start fails → software-only blind
  deployment defeated.*

**E2 (noise scaling):** oracle degrades gracefully toward 1.0 under regularization; blind never
approaches the bar at any dose (1e3–1e7). Delivered as table above + `E2_noise_curve.png`.

**E3 (scale):** **not run** — gated on E1 showing life, which it does not.

**Bottom line for the direction:** consistent with the R38 ruling. The oracle "second-DMD" effect
is real and the M=96 geometry is healthy, but **reference-free blind cold-start recovery of the
fixed null space fails** at the dev-grade F3 bar. Under the ruling's kill tree this is a **kill of
the software-only blind deployment**; any surviving value is the oracle/known-medium-law theorem
(materials-bank) or a much lower-dimensional scene class — not the arbitrary-dense-scene flagship.
"""

txt = open('FOG_DMD_PROBE_REPORT.md', encoding='utf-8').read()
txt = txt.replace('{{PRIMARY_TABLE}}', primary_md).replace('{{E2_TABLE}}', e2_md).replace('{{VERDICT}}', verdict)
open('FOG_DMD_PROBE_REPORT.md', 'w', encoding='utf-8').write(txt)
print("report filled. F3 pass =", passF3, "best blind =", best['blind_cold_nmse_median'])
