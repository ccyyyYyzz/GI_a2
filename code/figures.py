"""Publication-portable figures — spec §9: A1 phase diagram, B1 double-gate
heatmaps, GAM2 divergence curve (UT6), DT30 calibration curve.
Written to results/figures/ (PNG, 200 dpi)."""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gi_core import config as C

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
FIGS = os.path.join(RESULTS, "figures")


def fig_a1_phase_diagram():
    df = pd.read_csv(os.path.join(RESULTS, "phaseA_witness.csv"))
    fams = ["GAM1", "GAM2", "GAM3", "GAM4", "GAM8", "GAUSS", "CORR-LOGN", "MIX-LOGN"]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    meds, p90s = [], []
    for f in fams:
        e = df[df.family == f]["err"].astype(float)
        meds.append(np.median(e))
        p90s.append(np.percentile(e, 90))
    xs = np.arange(len(fams))
    ax.semilogy(xs, meds, "o-", label="median witness err")
    ax.semilogy(xs, p90s, "s--", label="P90 witness err")
    with open(os.path.join(RESULTS, "phaseA_gates.json")) as f:
        gates = json.load(f)["A1"]
    qs = []
    for f2 in fams:
        sub = df[df.family == f2]
        qs.append(np.median(sub["q_theory"].astype(float)))
    ax.semilogy(xs, qs, "k:", label="q (theory/interp SE)")
    ax.axhline(1.0, color="r", lw=0.8, alpha=0.5)
    ax.set_xticks(xs)
    ax.set_xticklabels(fams, rotation=30)
    ax.set_ylabel("relative witness error")
    ax.set_title("A1 Stein-validity phase diagram (Mw=2e5, n=256)\n"
                 "k=1 boundary failure | 1<k<=2 infinite variance | k>2 regular")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "A1_phase_diagram.png"), dpi=200)
    plt.close(fig)


def fig_gam2_divergence():
    with open(os.path.join(RESULTS, "unit_tests.json")) as f:
        ut = json.load(f)["UT6"]
    eps = np.array(ut["eps"])
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.semilogx(eps, ut["mc"], "o-", label="MC truncated E[s^2 1{a>eps}]")
    ax.semilogx(eps, ut["analytic"], "k--", label="analytic 4E1(2e)+(8e-4)e^{-2e}")
    ax.invert_xaxis()
    ax.set_xlabel("lower truncation eps")
    ax.set_ylabel("truncated second moment of GAM2 score")
    ax.set_title("UT6: GAM2 score second moment diverges (no saturation)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "GAM2_divergence.png"), dpi=200)
    plt.close(fig)

    # A1 GAM2 convergence-rate exhibit
    with open(os.path.join(RESULTS, "phaseA_gates.json")) as f:
        g = json.load(f)["A1"]["family_stats"]
    fig, ax = plt.subplots(figsize=(6, 4))
    for fam, style in (("GAM2", "o-"), ("GAM4", "s--")):
        d = g[fam]["median_err_by_Mw"]
        ms = sorted(int(k) for k in d)
        ax.loglog(ms, [d[str(m)] for m in ms], style, label=fam)
    ms = np.array(sorted(int(k) for k in g["GAM4"]["median_err_by_Mw"]))
    ref = g["GAM4"]["median_err_by_Mw"][str(ms[0])] * (ms / ms[0]) ** -0.5
    ax.loglog(ms, ref, "k:", label="M^{-1/2} reference")
    ax.set_xlabel("Mw")
    ax.set_ylabel("median witness err")
    ax.set_title("A1: GAM2 non-M^{-1/2} convergence signature")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "GAM2_convergence.png"), dpi=200)
    plt.close(fig)


def fig_dt30_calibration():
    from gi_core import links as L
    from gi_core.utils import rng_for

    u = np.linspace(0.05, 3.0, 30)
    fig, ax = plt.subplots(figsize=(6, 4))
    for s, c in ((1e3, "tab:blue"), (1e4, "tab:orange")):
        rng = rng_for(0, 91, int(s))
        sim = []
        for uv in u:
            lam = np.full(400, s * uv)
            N = L._dt30_counts(lam, s, rng)
            sim.append(N.mean() / s)
        ax.plot(u, sim, ".", color=c, label="renewal sim s=%g" % s)
    ax.plot(u, u / (1.0 + C.DT30_TAU_COEF * u), "k-", label="u/(1+(3/7)u)")
    ax.plot(u, u, "k:", alpha=0.5, label="identity")
    ax.set_xlabel("u")
    ax.set_ylabel("E[b|u]")
    ax.set_title("DT30 dead-time calibration (mean compression 0.7 at U0=1)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "DT30_calibration.png"), dpi=200)
    plt.close(fig)


def fig_b1_heatmaps():
    path = os.path.join(RESULTS, "phaseB_gates.json")
    if not os.path.exists(path):
        return
    with open(path) as f:
        gates = json.load(f)
    for gname in ("RELATION_HEADROOM", "PRACTICAL_HEADROOM"):
        detail = gates["combos"][gname]
        if not detail:
            continue
        keys = sorted(detail.keys())
        illums = sorted({k.split("|")[0] for k in keys})
        links = sorted({"%s|%s" % (k.split("|")[1], k.split("|")[2]) for k in keys})
        Z = np.full((len(illums), len(links)), np.nan)
        for k in keys:
            il, lk, ph = k.split("|")
            Z[illums.index(il), links.index("%s|%s" % (lk, ph))] = \
                detail[k]["mean_delta"]["PSNR"]
        fig, ax = plt.subplots(figsize=(1.1 * len(links) + 3, 0.8 * len(illums) + 2))
        im = ax.imshow(Z, cmap="RdBu_r", vmin=-2, vmax=2)
        for i in range(len(illums)):
            for j in range(len(links)):
                if np.isfinite(Z[i, j]):
                    k = "%s|%s" % (illums[i], links[j].replace("|", "|"))
                    kk = "%s|%s" % (illums[i], links[j])
                    passed = detail.get(kk, {}).get("pass", False)
                    ax.text(j, i, "%.2f%s" % (Z[i, j], "*" if passed else ""),
                            ha="center", va="center", fontsize=8)
        ax.set_xticks(range(len(links)))
        ax.set_xticklabels(links, rotation=45, ha="right")
        ax.set_yticks(range(len(illums)))
        ax.set_yticklabels(illums)
        ax.set_title("%s: mean dPSNR (dB), * = combo passes all gates" % gname)
        fig.colorbar(im, ax=ax, shrink=0.8)
        fig.tight_layout()
        fig.savefig(os.path.join(FIGS, "B1_%s.png" % gname.lower()), dpi=200)
        plt.close(fig)


def main():
    os.makedirs(FIGS, exist_ok=True)
    for fn in (fig_a1_phase_diagram, fig_gam2_divergence, fig_dt30_calibration,
               fig_b1_heatmaps):
        try:
            fn()
            print("figure ok:", fn.__name__, flush=True)
        except Exception as e:
            print("figure FAILED:", fn.__name__, repr(e), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
