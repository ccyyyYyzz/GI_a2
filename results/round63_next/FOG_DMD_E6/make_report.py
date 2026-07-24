#!/usr/bin/env python
# make_report.py -- merge E6_s1.json + E6_s2.json, build the verdict table, write
# E6_results.json (merged) + E6_REPORT.md. Robust to partial results.
import json, os, sys, numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
FILES = [os.path.join(BASE, 'E6_s1.json'), os.path.join(BASE, 'E6_s2.json')]
ARM_ORDER = ['A0', 'A1', 'A2', 'A3', 'A4']


def load_merged():
    merged = None
    arms = {}
    for f in FILES:
        if not os.path.exists(f):
            continue
        d = json.load(open(f))
        if merged is None:
            merged = dict(config=d.get('config'), null_fraction=d.get('null_fraction'))
        arms.update(d.get('arms', {}))
    if merged is None:
        merged = dict(config=None)
    merged['arms'] = arms
    return merged


def cell(arm, pk):
    a = M['arms'].get(arm, {})
    return a.get(pk, {}).get('summary', None)


def fmt(v, nd=3):
    return '  -  ' if v is None else (f'{v:.{nd}f}' if isinstance(v, float) else str(v))


def stay_verdict(ts):
    if ts is None:
        return '?'
    if ts <= 0.10:
        return 'STAY'
    if ts <= 0.35:
        return 'partial-drift'
    return 'DRIFT'


M = load_merged()
present = [a for a in ARM_ORDER if a in M['arms']]

# ---- merged results json ----
json.dump(M, open(os.path.join(BASE, 'E6_results.json'), 'w'), indent=1)

# ---- certificate table (per-seed, for the correlation endpoint) ----
cert_rows = []
for arm in present:
    for pk in ('clean', '100000'):
        cd = M['arms'][arm].get(pk, {})
        for r in cd.get('seeds', []):
            cert_rows.append(dict(arm=arm, photons=pk, seed=r['seed'],
                                  sf_fourier=r['cert']['sf_fourier'],
                                  sf_walsh=r['cert']['sf_walsh'],
                                  pw_blind=r['pw_blind_nmse'],
                                  pw_trueseed=r['pw_trueseed_nmse']))
json.dump(cert_rows, open(os.path.join(BASE, 'E6_whitening_certificate.json'), 'w'), indent=1)

# ---- build report ----
L = []
L.append('# E6 -- Does measurement-design redundancy kill the bilinear collusion?')
L.append('')
cfg = M.get('config') or {}
L.append(f"Cell: N={cfg.get('N')}, M={cfg.get('M')}, d={cfg.get('d')}, sigma_f={cfg.get('sigma_f')}, "
         f"OU tau={cfg.get('tau')}, T={cfg.get('T')} blocks x S={cfg.get('S')} slots "
         f"= {cfg.get('exposures')} exposures (fixed across arms). "
         f"Photons: clean + 1e5/bucket. Seeds: {cfg.get('seeds')}.")
fm = cfg.get('four_meta', {})
L.append(f"Fourier bank (A1/A3/A4): lattice spacing={fm.get('spacing')}, medium band={fm.get('band')}, "
         f"achieved sideband overlap ~= {fm.get('overlap_frac')}.")
L.append(f"Null fraction of scene (fraction of energy in ker(P)) ~= {fmt(M.get('null_fraction'))}.")
L.append('')
L.append('## Verdict table (arm x estimator x endpoint)')
L.append('')
L.append('null-NMSE (fixed-averaging baseline = 1.0; oracle ~0 = data contain the null). '
         'TRUE-SEED: refinement seeded from the true medium path, medium re-estimated -> '
         'STAY(<=0.10)=collusion killed, DRIFT(~0.7)=collusion present. '
         'BLIND=data-only cold start (5-seed median); agree=median within-seed std across starts.')
L.append('')
hdr = ('| arm | photons | oracle | **PW true-seed** | PW true-seed verdict | PW blind | PW blind agree '
       '| COV blind | COV blind agree | cert SF(Fourier) |')
L.append(hdr)
L.append('|' + '---|' * 10)
for arm in present:
    desc = M['arms'][arm].get('desc', '')
    for pk in ('clean', '100000'):
        s = cell(arm, pk)
        if s is None:
            continue
        L.append('| {a} | {p} | {orc} | **{ts}** | {tsv} | {bl} | {ba} | {cv} | {ca} | {sf} |'.format(
            a=arm, p=pk, orc=fmt(s['oracle_median']),
            ts=fmt(s['pw_trueseed_median']), tsv=stay_verdict(s['pw_trueseed_median']),
            bl=fmt(s['pw_blind_median']), ba=fmt(s['pw_blind_agreement']),
            cv=fmt(s['cov_blind_median']), ca=fmt(s['cov_blind_agreement']),
            sf=fmt(s['cert_sf_median'])))
L.append('')
L.append('Arm legend:')
for arm in present:
    L.append(f"- {arm}: {M['arms'][arm].get('desc','')}")
L.append('')

# ---- endpoint 3: certificate<->outcome correlation (A2/A3/A4 slot arms) ----
L.append('## Endpoint 3 -- whitening certificate vs outcome (slot arms)')
L.append('')
L.append('| arm | photons | cert SF median | corr(SF, blind null-NMSE) across seeds |')
L.append('|---|---|---|---|')
for arm in present:
    for pk in ('clean', '100000'):
        s = cell(arm, pk)
        if s is None:
            continue
        if s.get('cert_sf_median', 0) and s['cert_sf_median'] > 0.05:  # only slot arms have nonzero SF
            L.append(f"| {arm} | {pk} | {fmt(s['cert_sf_median'])} | {fmt(s['cert_vs_blind_corr'])} |")
L.append('')

# ---- automatic verdict summary ----
L.append('## Automatic verdict summary')
L.append('')
any_stay = []
any_full = []
for arm in present:
    for pk in ('clean', '100000'):
        s = cell(arm, pk)
        if s is None:
            continue
        ts = s['pw_trueseed_median']; bl = s['pw_blind_median']
        if ts is not None and ts <= 0.10:
            any_stay.append((arm, pk, ts, bl))
            if bl is not None and bl <= 0.30:
                any_full.append((arm, pk, ts, bl))
if any_full:
    L.append('SUCCESS SIGNAL MET: the following arm x estimator kill the collusion (true-seed STAY) '
             'AND recover blind (<=0.30):')
    for a, p, ts, bl in any_full:
        L.append(f"- {a} ({p}): true-seed {ts:.3f}, blind {bl:.3f}")
elif any_stay:
    L.append('PARTIAL SIGNAL: collusion KILLED (true-seed STAY) but blind cold-start remains >0.30 '
             '("collusion killed, optimization remains"):')
    for a, p, ts, bl in any_stay:
        L.append(f"- {a} ({p}): true-seed {ts:.3f}, blind {bl:.3f}")
else:
    L.append('NO ARM KILLED THE COLLUSION: every arm true-seed DRIFTS (>0.10) at matched refinement. '
             'Measurement-design redundancy did not convexify the pathwise medium re-estimation.')
L.append('')
open(os.path.join(BASE, 'E6_REPORT_TABLE.md'), 'w').write('\n'.join(L))
print('\n'.join(L))
print('\n[wrote E6_results.json, E6_whitening_certificate.json, E6_REPORT_TABLE.md]')
