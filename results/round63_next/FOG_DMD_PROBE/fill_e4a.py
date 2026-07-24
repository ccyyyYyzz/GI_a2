# Fill the {{E4A_TABLE}} placeholder in FOG_DMD_PROBE_REPORT.md from E4a_results.json.
import json
R = json.load(open('E4a_results.json'))
def ph(p): return 'clean' if p is None else f'{p:.0e}'
hdr = ("**Medium generated (ruling-declared), physical 0/1 patterns, 5-seed medians; null-NMSE. "
       "'frozen-known' = solve x with the medium held at truth (no re-estimation).**\n\n"
       "| medium | photons | oracle | frozen-known | true-Z-seed refine | baseline (random) | **SOBI-fix** | SOBI quality |\n"
       "|:--|:--:|:--:|:--:|:--:|:--:|:--:|:--:|\n")
body = ""
for r in R['results']:
    body += (f"| {r['medium']} | {ph(r['photons'])} | {r['oracle_nmse']:.3f} | "
             f"{r['frozen_knownmedium_nmse']:.3f} | {r['trueZ_seed_nmse']:.3f} | "
             f"{r['baseline_nmse_median']:.3f} | **{r['sobifix_nmse_median']:.3f}** | "
             f"{r['sobi_quality']:.2f} |\n")
tbl = hdr + body
mm = [r for r in R['results'] if r['medium'] == 'multi']
best = min(mm, key=lambda r: r['sobifix_nmse_median'])
tbl += (f"\nBest multi-τ SOBI-fix null-NMSE = **{best['sobifix_nmse_median']:.3f}** (ph={ph(best['photons'])}); "
        f"the SOBI init never beats the random baseline by more than noise, and even the frozen-known "
        f"medium ({best['frozen_knownmedium_nmse']:.3f}) is unreachable once the medium is re-estimated "
        f"(true-Z-seed {best['trueZ_seed_nmse']:.3f}).\n")
txt = open('FOG_DMD_PROBE_REPORT.md', encoding='utf-8').read()
txt = txt.replace('{{E4A_TABLE}}', tbl)
open('FOG_DMD_PROBE_REPORT.md', 'w', encoding='utf-8').write(txt)
print("E4a table filled. multi-tau SOBI-fix:", [round(r['sobifix_nmse_median'],3) for r in mm])
