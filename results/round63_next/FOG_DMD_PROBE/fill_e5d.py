# Fill {{E5D}} in FOG_DMD_PROBE_REPORT.md from E5d_results.json.
import json
R = json.load(open('E5d_results.json'))
axes = [('basis_rotation', 'basis rotation (estimator U rotated; simulator U_true)'),
        ('wrong_tau', 'wrong τ (estimator assumes 16)'),
        ('wrong_sigma_f', 'wrong σ_f (estimator assumes 0.15 → pure scale/gauge)'),
        ('dim_mismatch', 'subspace dim (estimator d=48)')]
out = []
for T in ['T2048', 'T1024']:
    b0 = R[T]['matched']
    out.append(f"**{T[1:]} epochs** — matched baseline null-NMSE = **{b0:.3f}**:\n")
    out.append("| mismatch axis | condition | null-NMSE | degradation |")
    out.append("|:--|:--|:--:|:--:|")
    for key, desc in axes:
        for i, r in enumerate(R[T][key]):
            ax = desc if i == 0 else ""
            ang = f" (∠{r['measured_mismatch']*100:.0f}%)" if r.get('measured_mismatch') else ""
            out.append(f"| {ax} | {r['label']}{ang} | {r['nmse']:.3f} | {r['degradation_pct']:+.0f}% |")
    out.append("")

# verdict
def find(T, ax, sub): return next((r for r in R[T][ax] if sub in r['label']), None)
rot10 = find('T2048', 'basis_rotation', '10%')
tau8 = find('T2048', 'wrong_tau', 'tau=8'); tau32 = find('T2048', 'wrong_tau', 'tau=32')
worst = max(abs(rot10['degradation_pct']), abs(tau8['degradation_pct']), abs(tau32['degradation_pct']))
robust = worst < 25
verdict = ("MISMATCH_ROBUST" if robust else "EXACT_LAW_LOAD_BEARING")
out.append(f"**E5d verdict: {verdict}.** In the F5 key cells @T=2048 the degradation is "
           f"rot-10% {rot10['degradation_pct']:+.0f}%, τ-mismatch (true 8 / 32) "
           f"{tau8['degradation_pct']:+.0f}% / {tau32['degradation_pct']:+.0f}% "
           f"(worst {worst:.0f}% {'<' if robust else '≥'} 25%). ")
if robust:
    out.append("The estimator does **not** require the exact simulator law: τ and σ_f mismatches are "
               "absorbed by the scalar gauge (every lag stays ∝ G_P, so a wrong τ/σ_f only rescales the "
               "target), and moderate basis rotation degrades gracefully. **E5a graduates from signal to "
               "candidate.** The load-bearing assumption is the *subspace* (U's span and its dimension), "
               "not its exact realization — consistent with a declared-medium-class model, not an "
               "exact-simulator artifact.")
else:
    out.append("An exact simulator-matched law is load-bearing → per R38 F5 this is a simulation "
               "artifact and the covariance arm is killed.")
txt = open('FOG_DMD_PROBE_REPORT.md', encoding='utf-8').read()
txt = txt.replace('{{E5D}}', "\n".join(out))
open('FOG_DMD_PROBE_REPORT.md', 'w', encoding='utf-8').write(txt)
print("E5d filled. verdict:", verdict, "| placeholders:", txt.count('{{'))
