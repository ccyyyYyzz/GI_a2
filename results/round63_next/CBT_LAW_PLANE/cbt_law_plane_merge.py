#!/usr/bin/env python
# cbt_law_plane_merge.py -- deterministic merge of the CBT_LAW_PLANE shards (JOB 2).
# Combines the analytic sweep + the 4 MC corner shards into CBT_LAW_PLANE.json, re-evaluates
# the frozen bars + the MC on-orbit-AUC bar, stamps venue metadata, and renders the hero figure.
import glob
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    p = os.path.join(HERE, name)
    return json.load(open(p)) if os.path.exists(p) else None


def main():
    analytic = load("CBT_LAW_PLANE_analytic.json")
    assert analytic is not None, "missing CBT_LAW_PLANE_analytic.json"
    corners = []
    for k in range(4):
        c = load("CBT_LAW_PLANE_mc_corner%d.json" % k)
        if c is not None:
            corners.append(c["corner"])
    corners = sorted(corners, key=lambda r: r["corner_idx"])

    band = analytic["frozen"]["bars"]["MC_paired_auc_on_orbit_in"]
    mc_bar = {
        "band": band,
        "per_corner": [{"corner_idx": c["corner_idx"], "r0_scale": c["r0_scale"],
                        "omega": c["omega"], "paired_auc": c["paired_auc"],
                        "in_band": c["paired_auc_in_band"]} for c in corners],
        "all_in_band": bool(len(corners) == 4 and all(c["paired_auc_in_band"] for c in corners)),
        "n_corners": len(corners),
    }

    payload = {
        "meta": {"venue": "colab_a100", "job": "JOB2_CBT_LAW_PLANE",
                 "analytic_arm": "float64 deterministic (venue-independent; reproduced bitwise vs local)",
                 "mc_arm": "colab A100 GPU generator, CRN-paired on-orbit AUC",
                 "sessions": ["cbt_camp1"]},
        "frozen": analytic["frozen"],
        "config": analytic["config"],
        "analytic_bars": analytic["bars"],
        "mc_onorbit_auc_bar": mc_bar,
        "operating_points": analytic["operating_points"],
        "mc_corners": corners,
        "verdict_bars": {
            "B1": analytic["bars"]["B1_median_kink_rel_err"]["pass"],
            "B2": analytic["bars"]["B2_median_supercritical_C_ratio"]["pass"],
            "B3": analytic["bars"]["B3_min_supercritical_quartic_R2"]["pass"],
            "B4": analytic["bars"]["B4_max_subcritical_C_fit_over_C0"]["pass"],
            "MC_on_orbit_auc": mc_bar["all_in_band"],
        },
    }
    payload["all_bars_pass"] = bool(all(payload["verdict_bars"].values()))
    out = os.path.join(HERE, "CBT_LAW_PLANE.json")
    with open(out, "w") as f:
        json.dump(payload, f, indent=2, default=float)
    print("wrote", out, "| bars:", payload["verdict_bars"], "| all_pass:", payload["all_bars_pass"])

    make_figure(payload)


def make_figure(payload):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pts = payload["operating_points"]
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    cmap = plt.get_cmap("viridis")

    # Panel 1: universal collapse  C_fit/C0 vs A/A*  onto  [1-A/A*]_+^2
    xx = np.linspace(0, 2.2, 200)
    ax[0].plot(xx, np.clip(1 - xx, 0, None) ** 2, "k-", lw=1.6, label=r"$[1-A/A_*]_+^2$")
    for i, p in enumerate(pts):
        col = cmap(i / max(len(pts) - 1, 1))
        x = [1.0 / r["chi"] for r in p["chi_rows"]]           # A_traj/A* = 1/chi
        y = [(r["C_fit_over_C0"] if r["C_fit_over_C0"] is not None else np.nan) for r in p["chi_rows"]]
        ax[0].plot(x, y, "o", color=col, ms=7, alpha=0.85,
                   label="r0x%.1f w%.1f" % (p["r0_scale"], p["omega"]))
    ax[0].axvline(1.0, color="red", ls="--", lw=1)
    ax[0].set_xlabel(r"budget ratio $A/A_*$"); ax[0].set_ylabel(r"$C_{\rm fit}/C_0$")
    ax[0].set_title("Universal curvature-tax collapse\n(all operating points)")
    ax[0].legend(fontsize=6, ncol=2); ax[0].grid(True, alpha=.15); ax[0].set_ylim(-0.05, 1.1)

    # Panel 2: supercritical C_fit vs C_pred (on y=x)
    cf, cp = [], []
    for p in pts:
        for r in p["chi_rows"]:
            if r["supercritical"] and r["C_pred"] > 0:
                cf.append(r["C_fit"]); cp.append(r["C_pred"])
    cf, cp = np.array(cf), np.array(cp)
    lim = [min(cp.min(), cf.min()) * 0.7, max(cp.max(), cf.max()) * 1.3]
    ax[1].plot(lim, lim, "k-", lw=1, alpha=.6, label="y=x")
    ax[1].loglog(cp, cf, "o", color="#1b7837", ms=7, alpha=.8)
    ax[1].set_xlabel(r"predicted $C_{\rm pred}=\frac{\|K\|_F^2}{16}(s^2-Ar_0)^2$")
    ax[1].set_ylabel(r"fitted $C_{\rm fit}$ (Chernoff)")
    ax[1].set_title("Coefficient law verified\n(supercritical branch, %d pts)" % len(cf))
    ax[1].legend(fontsize=8); ax[1].grid(True, which="both", alpha=.15)

    # Panel 3: kink A_kink vs A*  +  MC corners on-orbit AUC inset text
    As = [p["A_star"] for p in pts]; Ak = [p["A_kink"] for p in pts]
    lim2 = [min(As) * 0.7, max(As) * 1.3]
    ax[2].plot(lim2, lim2, "k-", lw=1, alpha=.6, label="y=x")
    ax[2].loglog(As, Ak, "s", color="#2166ac", ms=8, alpha=.85)
    ax[2].set_xlabel(r"critical budget $A_*=s^2/r_0$"); ax[2].set_ylabel(r"measured kink $A_{\rm kink}$")
    ax[2].set_title("Kink locates $A_*$\n(6 operating points)")
    mcb = payload.get("mc_onorbit_auc_bar", {})
    if mcb.get("per_corner"):
        txt = "MC on-orbit AUC (B=%d):\n" % payload["frozen"]["bars"]["MC_B"]
        for c in mcb["per_corner"]:
            txt += "  r0x%.1f w%.1f: %.4f%s\n" % (c["r0_scale"], c["omega"], c["paired_auc"],
                                                  " OK" if c["in_band"] else " !")
        ax[2].text(0.03, 0.97, txt, transform=ax[2].transAxes, fontsize=7, va="top",
                   family="monospace", bbox=dict(fc="white", alpha=.7, ec="grey"))
    ax[2].legend(fontsize=8, loc="lower right"); ax[2].grid(True, which="both", alpha=.15)

    fig.suptitle("CBT_LAW_PLANE -- calibration plane of the curvature-tax law "
                 "(kink at $A_*$, coefficient $\\propto[1-A/A_*]_+^2$)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = os.path.join(HERE, "cbt_law_plane_hero.png")
    fig.savefig(out, dpi=130); plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
