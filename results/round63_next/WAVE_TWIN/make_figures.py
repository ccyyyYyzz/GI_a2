#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Render the WAVE_TWIN figures from the T1-T5 JSON results."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIG = "figs"
os.makedirs(FIG, exist_ok=True)
def load(p):
    with open(p) as f: return json.load(f)

# ---------- Fig 1: T1 wall-leak (z1 sweep + bar) ----------
try:
    t1 = load("T1_WALL_LEAK.json")
    A = t1["partA_z1_sweep_z2_0_deterministic"]
    z1 = [r["z1_mm"] for r in A]; lk = [r["leak_rel_inband"] for r in A]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].semilogy(z1, np.maximum(lk, 1e-18), 'o-', color="#c1121f", lw=2, ms=7)
    ax[0].axhline(1e-3, ls='--', color='k', label='flag threshold 1e-3')
    ax[0].axhline(5.68e-16, ls=':', color='gray', label='sealed wall 5.7e-16')
    ax[0].set_xlabel("DMD->diffuser distance z1 (mm)"); ax[0].set_ylabel("mean-wall leak (rel. in-band)")
    ax[0].set_title("T1: deterministic mean-wall leak vs z1 (z2=0)")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=.3, which='both')
    # right panel: operational mean d'/cov d' from T5 (the leak that actually matters)
    try:
        t5 = load("T5_SENTINEL.json")
        labels=[]; ratios=[]
        for g in t5["geometries"]:
            for r in g["rows"]:
                labels.append(f"{g['geometry'].split('_z2')[0][:9]}\nz2={g['z2_mm']} {r['change']}")
                ratios.append(r.get("mean_dp_over_cov_dp", 0))
        ax[1].bar(range(len(labels)), ratios, color="#457b9d")
        ax[1].axhline(1e-3, ls='--', color='k', label='1e-3')
        ax[1].set_yscale('log'); ax[1].set_xticks(range(len(labels)))
        ax[1].set_xticklabels(labels, fontsize=6.5)
        ax[1].set_ylabel("operational mean d' / cov d'")
        ax[1].set_title("T1/T5: operational wall leak (mean vs cov channel)")
        ax[1].legend(fontsize=8); ax[1].grid(alpha=.3, axis='y')
    except Exception as ee:
        ax[1].text(0.5,0.5,"T5 pending", ha='center')
    plt.tight_layout(); plt.savefig(f"{FIG}/fig1_wall_leak.png", dpi=130); plt.close()
    print("fig1 ok")
except Exception as e:
    print("fig1 skip:", e)

# ---------- Fig 2: T2 grain curves (grain & contrast vs z2, l_c) ----------
try:
    t2 = load("T2_SPECKLE_STATS.json")
    rows = t2["rows"]
    dev = [r for r in rows if r["regime"]=="developed_2pi"]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for l_c in [5,15,50]:
        sub = [r for r in dev if r["l_c_um"]==l_c]
        z = [r["z2_mm"] for r in sub]
        ax[0].plot(z, [r["grain_um"] for r in sub], 'o-', label=f"l_c={l_c}um")
        ax[1].plot(z, [r["contrast"] for r in sub], 's-', label=f"l_c={l_c}um")
    ax[0].set_xlabel("z2 (mm)"); ax[0].set_ylabel("speckle grain (um)")
    ax[0].set_title("T2: grain size vs z2 (developed)"); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
    ax[1].axhline(1.0, ls=':', color='gray', label='developed contrast=1')
    ax[1].set_xlabel("z2 (mm)"); ax[1].set_ylabel("speckle contrast")
    ax[1].set_title("T2: contrast vs z2 (developed)"); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
    plt.tight_layout(); plt.savefig(f"{FIG}/fig2_grain_curves.png", dpi=130); plt.close()
    print("fig2 ok")
except Exception as e:
    print("fig2 skip:", e)

# ---------- Fig 3: T3 aperture curves + edge vs grain ----------
try:
    t3 = load("T3_APERTURE_LAW.json")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    kp = t3["kp"]
    for cur in t3["curves"]:
        ks = cur["ks"]; r = np.array(cur["response"]); r = r/ (r[1:6].mean()+1e-30)
        ax[0].plot(ks, r, 'o-', ms=4, label=f"l_c={cur['l_c_um']}um (k_w~{cur['k_w_mapped']})")
        ax[0].axvline(cur["edge_pred_kp_plus_kw"], ls='--', alpha=.5)
    ax[0].axvline(kp, ls=':', color='k', label=f'k_p={kp}')
    ax[0].set_xlabel("scene Chebyshev freq k"); ax[0].set_ylabel("cov response (norm.)")
    ax[0].set_title("T3: aperture curve vs grain (dashed = k_p+k_w)"); ax[0].legend(fontsize=7); ax[0].grid(alpha=.3)
    gw = [c["grain_um"] for c in t3["curves"]]
    em = [c.get("edge_measured_40pct", c["edge_measured"]) for c in t3["curves"]]
    ax[1].plot(gw, em, 's-', color="#e63946", label="measured usable edge (40% peak)")
    ax[1].axhline(t3.get("nyquist",32), ls=':', color='gray', label='64-grid Nyquist')
    ax[1].set_xlabel("grain size (um)"); ax[1].set_ylabel("usable scene-freq edge")
    ax[1].set_title("T3: usable edge vs grain size"); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
    plt.tight_layout(); plt.savefig(f"{FIG}/fig3_edge_vs_grain.png", dpi=130); plt.close()
    print("fig3 ok")
except Exception as e:
    print("fig3 skip:", e)

# ---------- Fig 4: T4 multiplicative->convolutive ----------
try:
    t4 = load("T4_MULT_TO_CONV.json")
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for reg in t4["regimes"]:
        z = [r["z2_mm"] for r in reg["rows"]]
        ic = [r["intensity_obj_corr"] for r in reg["rows"]]
        ax.plot(z, ic, 'o-', label=reg["regime"])
    for lvl in [0.9, 0.75, 0.5]:
        ax.axhline(lvl, ls=':', alpha=.4)
    ax.set_xlabel("z2 (mm)  [screen-object separation]"); ax.set_ylabel("field/intensity correlation")
    ax.set_title("T4: multiplicative->convolutive (corr vs z2)")
    ax.legend(fontsize=8); ax.grid(alpha=.3)
    plt.tight_layout(); plt.savefig(f"{FIG}/fig4_mult_to_conv.png", dpi=130); plt.close()
    print("fig4 ok")
except Exception as e:
    print("fig4 skip:", e)

# ---------- Fig 5: T5 ROC / d' bars ----------
try:
    t5 = load("T5_SENTINEL.json")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    labels=[]; tdet=[]; auc=[]
    for g in t5["geometries"]:
        for r in g["rows"]:
            labels.append(f"{g['geometry'].split('_z2')[0][:9]}\nz2={g['z2_mm']} {r['change']}")
            tdet.append(r["T_det_dprime5"]); auc.append(r["auc_at_T453"])
    x = np.arange(len(labels))
    ax[0].bar(x, tdet, color="#2a9d8f")
    ax[0].axhline(453, ls='--', color='k', label='sealed best 453 banks')
    ax[0].set_yscale('log')
    ax[0].set_xticks(x); ax[0].set_xticklabels(labels, fontsize=6.5, rotation=0)
    ax[0].set_ylabel("T_det (banks) for d'=5  [log]"); ax[0].set_title("T5: detection time vs sealed 453")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=.3, axis='y', which='both')
    ax[1].bar(x, auc, color="#e76f51")
    ax[1].axhline(1.0, ls=':', color='gray')
    ax[1].set_ylim(0.4, 1.02); ax[1].set_xticks(x); ax[1].set_xticklabels(labels, fontsize=6.5)
    ax[1].set_ylabel("AUC @ T_eff=453"); ax[1].set_title("T5: AUC at sealed budget")
    ax[1].grid(alpha=.3, axis='y')
    plt.tight_layout(); plt.savefig(f"{FIG}/fig5_roc_dprime.png", dpi=130); plt.close()
    print("fig5 ok")
except Exception as e:
    print("fig5 skip:", e)

# ---------- Fig 6: T5b M-scaling (comparability correction) ----------
try:
    t5b = load("T5b_MSCALING.json")
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    Ms = t5b["M_values"]
    for g in t5b["geometries"]:
        for r in g["rows"]:
            lam = [r["scaling"][str(m)]["cov_lam"] for m in Ms]
            ax.loglog(Ms, lam, 'o-', label=f"{g['geometry'].split('_z2')[0]} {r['change']} (p={r['scaling_exponent_p']})")
    # reference M^2 and M^1 slopes
    m0 = np.array(Ms, float)
    ax.loglog(m0, 3e-6*(m0/m0[0])**2, 'k--', alpha=.4, label='M^2 (naive cov)')
    ax.loglog(m0, 3e-6*(m0/m0[0])**1, 'k:', alpha=.4, label='M^1 (mean)')
    ax.set_xlabel("number of codes M"); ax.set_ylabel("per-bank cov noncentrality λ")
    ax.set_title("T5b: λ ∝ M^1.8 (cov NOT saturated at finite z2)")
    ax.legend(fontsize=7); ax.grid(alpha=.3, which='both')
    plt.tight_layout(); plt.savefig(f"{FIG}/fig6_mscaling.png", dpi=130); plt.close()
    print("fig6 ok")
except Exception as e:
    print("fig6 skip:", e)

print("figures done")
