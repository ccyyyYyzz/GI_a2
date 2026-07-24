# GATE G1 -- FINE-BAND MISMATCH (the sharpest risk).
#
# The beyond-band claim leans on the estimator KNOWING the medium's fine spectral band.
# Data are generated with the TRUE fine band (1<=i+j<=6, flat radial profile). The estimator
# is fed a PERTURBED declared band:
#   (a) rotated 10% / 20% within the fine annulus (true modes tilted into the just-beyond 7-8 shell)
#   (b) too narrow : declares i+j<=5 (truth <=6)
#   (c) too wide   : declares i+j<=7 (truth <=6)
#   (d) wrong radial profile: estimator declares FLAT while the true medium has a 1/f radial
#       power profile across the fine annulus.
# Degradation is reported vs the MATCHED fixed+cov moment baseline (R39's 0.47-0.79).
# KILL SIGNAL: degradation > 50% at the 10-20% perturbations.
import json
import time

import numpy as np

import bb_common as bb

T = 2048
N_SEED = 3
STEPS = 4000
N_STARTS = 3
t0 = time.time()


def rotate_within_annulus(Ub_true, eps, lo=bb.MED_LO, hi=bb.MED_HI, shell_hi=8,
                          n=bb.N_SIDE, seed=0):
    """Declared subspace = true fine band tilted by angle arcsin(eps) toward the just-beyond
    fine shell (i+j in (hi, shell_hi]) -- the medium's own neighboring fine modes, so the
    perturbation stays 'within the fine annulus'. eps is the tilted-energy fraction."""
    shell = [(i, j) for i in range(n) for j in range(n) if hi < i + j <= shell_hi]
    Up = np.linalg.qr(np.array([bb.dct_spike(i, j, n) for (i, j) in shell]).T)[0]
    db = Ub_true.shape[1]
    k = min(db, Up.shape[1])                       # rotate the first k declared directions
    Urot = Ub_true.copy().astype(float)
    Urot[:, :k] = np.sqrt(1 - eps ** 2) * Ub_true[:, :k] + eps * Up[:, :k]
    return np.linalg.qr(Urot)[0]


def run_condition(tag, gen_gen_sd, declared_Ub, declared_sd):
    """gen_gen_sd: per-mode sd vector for the TRUE medium (None -> flat coeff_sd on the TRUE band).
    declared_Ub/declared_sd: what the estimator is told."""
    errs = []
    for s in range(N_SEED):
        W = bb.medium_field(T, 600 + s, Ub_true, coeff_sd_true, gen_sd=gen_gen_sd)
        B, Rd, scp = bb.buckets(Pfix, W, x, 650 + s)
        xm, _ = bb.est_moment(Pfix, declared_Ub, B, declared_sd,
                              n_starts=N_STARTS, steps=STEPS)
        errs.append(bb.hi_err(xm, x, hi_mask))
    return float(np.median(errs)), errs


# --- true cell (flat medium, band 1<=i+j<=6) ---
Ub_true, db_true, coeff_sd_true, Kg_true = bb.medium_subspace(bb.MED_LO, bb.MED_HI)
x, hi_mask = bb.make_scene()
Pfix = bb.bl_patterns(bb.MED_PAT if hasattr(bb, "MED_PAT") else bb.M_PAT, 10)

rows = []
print(f"G1 FINE-BAND MISMATCH  T={T} seeds={N_SEED} true db={db_true}", flush=True)

# matched baseline
base_med, base_all = run_condition("matched", None, Ub_true, coeff_sd_true)
rows.append(dict(cond="matched (band 1..6, flat)", err=base_med, all=base_all, degr=0.0))
print(f"  matched baseline           err {base_med:.3f}  {[round(q,3) for q in base_all]}", flush=True)


def degr(m):
    return (m - base_med) / base_med


# (a) rotation within the fine annulus
for eps in (0.10, 0.20):
    Ud = rotate_within_annulus(Ub_true, eps)
    m, a = run_condition(f"rot{int(eps*100)}", None, Ud, coeff_sd_true)
    rows.append(dict(cond=f"(a) rotated {int(eps*100)}% within fine annulus",
                     err=m, all=a, degr=degr(m)))
    print(f"  (a) rot {int(eps*100):2d}%              err {m:.3f}  degr {degr(m)*100:+.0f}%  {[round(q,3) for q in a]}", flush=True)

# (b) too narrow: declare 1..5
Ub_n, db_n, sd_n, _ = bb.medium_subspace(1, 5)
m, a = run_condition("narrow", None, Ub_n, sd_n)
rows.append(dict(cond=f"(b) too narrow: declare i+j<=5 (db={db_n})", err=m, all=a, degr=degr(m)))
print(f"  (b) too narrow (<=5)        err {m:.3f}  degr {degr(m)*100:+.0f}%  {[round(q,3) for q in a]}", flush=True)

# (c) too wide: declare 1..7
Ub_w, db_w, sd_w, _ = bb.medium_subspace(1, 7)
m, a = run_condition("wide", None, Ub_w, sd_w)
rows.append(dict(cond=f"(c) too wide: declare i+j<=7 (db={db_w})", err=m, all=a, degr=degr(m)))
print(f"  (c) too wide (<=7)          err {m:.3f}  degr {degr(m)*100:+.0f}%  {[round(q,3) for q in a]}", flush=True)

# (d) wrong radial profile: TRUE medium has 1/f radial profile; estimator declares flat.
modes = [(i, j) for i in range(bb.N_SIDE) for j in range(bb.N_SIDE) if bb.MED_LO <= i + j <= bb.MED_HI]
prof = np.array([1.0 / (1.0 + np.hypot(i, j)) for (i, j) in modes])
gen_sd = coeff_sd_true * prof / np.sqrt((prof ** 2).mean())   # keeps pixelwise RMS = SIG_F
m, a = run_condition("profile", gen_sd, Ub_true, coeff_sd_true)   # declared = flat
rows.append(dict(cond="(d) wrong radial profile (true 1/f, declared flat)", err=m, all=a, degr=degr(m)))
print(f"  (d) wrong radial profile    err {m:.3f}  degr {degr(m)*100:+.0f}%  {[round(q,3) for q in a]}", flush=True)

# verdict: kill if any 10-20% perturbation degrades > 50%
key = [r for r in rows if r["cond"].startswith("(a)")]
worst_ab = max(r["degr"] for r in rows[1:])
kill = any(r["degr"] > 0.50 for r in key)
verdict = "KILL" if kill else ("WATCH" if worst_ab > 0.25 else "PASS")
out = dict(gate="G1_fine_band_mismatch", T=T, seeds=N_SEED, matched_baseline=base_med,
           rows=rows, worst_degradation=worst_ab,
           kill_rule="degradation > 50% at 10-20% (rotation) perturbations",
           verdict=verdict, wall_s=time.time() - t0)
json.dump(out, open("G1_results.json", "w"), indent=2)
print(f"\nG1 VERDICT: {verdict}  (worst degradation across all axes {worst_ab*100:+.0f}%; "
      f"rotation-axis kill rule fired: {kill})  [{time.time()-t0:.0f}s]", flush=True)
