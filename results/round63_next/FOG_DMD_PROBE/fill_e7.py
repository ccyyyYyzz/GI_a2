# Fill §13 (E7) placeholders from E7_results.json.
import json
R = json.load(open('E7_results.json'))
def cross(rows, key):
    b = [r['T'] for r in rows if r[key] <= 0.25]; return min(b) if b else None

# E7a tables
e7a = []
for tag, title in [('mini', 'Mini (16×16, M=48, d=24)'), ('full32', 'Full (32×32, M=96, d=48)')]:
    if tag not in R.get('E7a', {}): continue
    rows = R['E7a'][tag]
    e7a.append(f"**{title}** — null-NMSE, 1e5 photons (wall-clock/fit in parentheses):\n")
    e7a.append("| T | moment-matching (E5a) | **two-stage MLE** | t_moment | t_MLE |")
    e7a.append("|:--:|:--:|:--:|:--:|:--:|")
    for r in rows:
        e7a.append(f"| {r['T']} | {r['moment_nmse']:.3f} | **{r['mle_nmse']:.3f}** | "
                   f"{r['sec_moment']:.1f}s | {r['sec_mle']:.1f}s |")
    cm, cl = cross(rows, 'moment_nmse'), cross(rows, 'mle_nmse')
    e7a.append(f"\nCrossing 0.25: moment at T≈{cm}, **MLE at T≈{cl}**"
               + (f" (**{cm/cl:.1f}× reduction**)." if (cm and cl and cl < cm) else ".") + "\n")

# E7b table
e7b = []
if 'E7b' in R:
    for tag, rows in R['E7b'].items():
        e7b.append(f"**Mini d=24**, fixed total exposures m·T = {rows[0]['mT']}:\n")
        e7b.append("| m (patterns/epoch) | T (epochs) | χ_full | MLE null-NMSE |")
        e7b.append("|:--:|:--:|:--:|:--:|")
        for r in rows:
            deg = "degenerate (m≤d)" if r['chi'] <= 0 else f"{r['mle_nmse']:.3f}"
            e7b.append(f"| {r['m']} | {r['T']} | {r['chi']:.2f} | {deg} |")
        viable = [r for r in rows if r['chi'] > 0]
        if viable:
            best = min(viable, key=lambda r: r['mle_nmse'])
            e7b.append(f"\nExposure-optimal viable point: **m={best['m']}, T={best['T']}** "
                       f"(null-NMSE {best['mle_nmse']:.3f}). ")
            if len(viable) >= 2:
                v = sorted(viable, key=lambda r: r['m'])
                trend = ("smaller m (more, cheaper epochs) HELPS" if v[0]['mle_nmse'] < v[-1]['mle_nmse'] - 0.02
                         else "smaller m HURTS" if v[0]['mle_nmse'] > v[-1]['mle_nmse'] + 0.02
                         else "m is roughly neutral")
                e7b.append(f"At fixed exposure, **{trend}** ({v[0]['m']}→{v[-1]['m']}: "
                           f"{v[0]['mle_nmse']:.3f}→{v[-1]['mle_nmse']:.3f}). ")
            e7b.append("Below m=d the estimator is degenerate (no blind capacity), so subsampling "
                       "cannot go arbitrarily cheap — the medium dimension d floors patterns-per-epoch.\n")

# verdict
mi = R['E7a']['mini']; cm, cl = cross(mi, 'moment_nmse'), cross(mi, 'mle_nmse')
f32 = R['E7a'].get('full32', [])
cl32 = cross(f32, 'mle_nmse') if f32 else None
red = (cm / cl) if (cm and cl) else None
verdict = (f"**E7 verdict.** The marginal-likelihood MLE (moment-matching init → EM) is the "
           f"efficiency win: on the mini cell it crosses 0.25 at T≈{cl} vs T≈{cm} for moment-matching"
           + (f" (**{red:.1f}× sample-complexity reduction**)" if red and red >= 1.5 else "")
           + (f"; on the 32×32 cell the MLE crosses at T≈{cl32}" if cl32 else "")
           + ". Exposure economics: patterns-per-epoch is floored by the medium dimension (m>d), so "
           "the exposure budget is best spent on epochs at m just above d, not on the full bank if "
           "d≪M. These numbers anchor the sample-complexity (E7a) and optimal-design (E7b) questions "
           "for R39.")

txt = open('FOG_DMD_PROBE_REPORT.md', encoding='utf-8').read()
txt = txt.replace('{{E7A}}', "\n".join(e7a)).replace('{{E7B}}', "\n".join(e7b)).replace('{{E7VERDICT}}', verdict)
open('FOG_DMD_PROBE_REPORT.md', 'w', encoding='utf-8').write(txt)
print("E7 filled. mini crossings moment/MLE:", cm, cl, "| placeholders:", txt.count('{{'))
