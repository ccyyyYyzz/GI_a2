# Fill §12 (E5) placeholders in FOG_DMD_PROBE_REPORT.md from E5_results.json / E5b_results.json.
import json, os
R = json.load(open('E5_results.json'))
def ph(p): return 'clean' if p is None else f'{p:.0e}'

# agreement sentence
a = R['cell_T128'][0]
agree = (f"Across 6 data-only random starts the null-NMSE agrees to std **{a['agreement_std']:.3f}** "
         f"(at T=128, clean) — every start lands in the same basin. There is no realization vector to "
         f"trade against, so the E4 drift is structurally impossible here.")

# T-scaling table
ts = "| T (epochs) | " + " | ".join(str(r['T']) for r in R['T_scaling']) + " |\n"
ts += "|:--|" + ":--:|" * len(R['T_scaling']) + "\n"
ts += "| blind null-NMSE | " + " | ".join(f"**{r['nmse_median']:.3f}**" for r in R['T_scaling']) + " |\n"
Tcurve = {r['T']: r['nmse_median'] for r in R['T_scaling']}
Tstar = next((T for T in sorted(Tcurve) if Tcurve[T] <= 0.25), None)

# cell table
cell = "Blind cell (T=128, 5-seed medians; multi-start agreement std in parentheses):\n\n"
cell += "| photons/bucket | blind null-NMSE | agreement std |\n|:--:|:--:|:--:|\n"
for r in R['cell_T128']:
    cell += f"| {ph(r['photons'])} | {r['nmse_median']:.3f} | {r['agreement_std']:.3f} |\n"

# photon robustness at T=1024
pht = "Photon-robustness at T=1024 (lags ℓ≥1 are shot-free):\n\n"
pht += "| photons/bucket | blind null-NMSE |\n|:--:|:--:|\n"
for r in R['photon_at_T1024']:
    pht += f"| {ph(r['photons'])} | {r['nmse_median']:.3f} |\n"

txt = open('FOG_DMD_PROBE_REPORT.md', encoding='utf-8').read()
txt = txt.replace('{{E5A_AGREE}}', agree).replace('{{E5A_TSCALE}}', ts)
txt = txt.replace('{{E5A_CELL}}', cell).replace('{{E5A_PHOTON}}', pht)
txt = txt.replace('{{E5A_TSTAR}}', str(Tstar) if Tstar else '>4096')

# E5b
if os.path.exists('E5b_results.json'):
    B = json.load(open('E5b_results.json'))['rows']
    sp = [r for r in B if r['kind'] == 'sparse'][0]
    dn = [r for r in B if r['kind'] == 'dense_ref'][0]
    broke = sp['trueZ_sparse'] < 0.5 * sp['trueZ_noSparse'] and sp['trueZ_sparse'] < 0.35
    e5b = ("Sparse scene (few bright squares; T=128, single-τ OU). §11c-style table, null-NMSE:\n\n"
           "| scene | oracle | frozen-known | true-Z-seed (no L1) | true-Z-seed (+L1) | baseline (+L1) |\n"
           "|:--|:--:|:--:|:--:|:--:|:--:|\n"
           f"| dense (control) | {dn['oracle']:.3f} | {dn['frozen_known']:.3f} | {dn['trueZ_noSparse']:.3f} | {dn['trueZ_sparse']:.3f} | {dn['baseline_sparse']:.3f} |\n"
           f"| **sparse** | {sp['oracle']:.3f} | {sp['frozen_known']:.3f} | {sp['trueZ_noSparse']:.3f} | **{sp['trueZ_sparse']:.3f}** | {sp['baseline_sparse']:.3f} |\n\n"
           f"**SPARSITY_BREAKS_COLLUSION: {'YES' if broke else 'NO'}** — with the L1/sparsity constraint the "
           f"true-Z-seeded pathwise refinement goes {sp['trueZ_noSparse']:.3f} → {sp['trueZ_sparse']:.3f} "
           f"(no-L1 → +L1). " + ("The structural prior penalizes the dense wrong-null collusion partner, so truth "
           "becomes the preferred pathwise minimum." if broke else "The sparsity prior does not by itself "
           "restore truth as the pathwise minimum here."))
    txt = txt.replace('{{E5B}}', e5b)

txt = txt.replace('{{E5C}}', open('E5c_note.txt').read() if os.path.exists('E5c_note.txt')
                  else "_E5c (anchor frames) not run — E5a already resurrects the direction; deprioritized per instruction._")
open('FOG_DMD_PROBE_REPORT.md', 'w', encoding='utf-8').write(txt)
print("E5 filled. Tstar:", Tstar, "| remaining placeholders:", txt.count('{{'))
