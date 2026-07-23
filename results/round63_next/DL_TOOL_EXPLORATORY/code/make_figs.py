#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""EXPLORATORY / DEV-ONLY figures for E1 (analytic-fix head-to-head) and E2
(learned temporal prior, honest recon loop).  Reads the json outputs; no compute."""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "D:/GI_another/results/round63_next/DL_TOOL_EXPLORATORY"
J = os.path.join(OUT, "json"); F = os.path.join(OUT, "figs")


def load(n): return json.load(open(os.path.join(J, n)))


# ---------------- E1: 4-cell analytic-fix head-to-head ----------------
def fig_e1():
    hh = load("e1_fix_headtohead.json"); cc = load("e1c_corrector.json")
    cells = hh["cells"]; corr = {(c["tc"], c["cv"]): c for c in cc["cells"]}
    arms = ["arm0", "D0", "D_arith", "corr", "L_mix", "oracle"]
    labels = ["no-corr\nfloor", "DLGI\nfrozen", "DLGI\n+E[a]=1 fix", "  +learned\n  corrector",
              "SCGI\n(learned)", "oracle\n(ceiling)"]
    colors = ["#8c8c8c", "#d1495b", "#edae49", "#66a182", "#2e86ab", "#111111"]
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.4), sharey=True)
    for ax, c in zip(axes, cells):
        p = c["psnr"]; key = (c["tc"], c["cv"])
        vals = [p["arm0"], p["D0"], p["D_arith"], corr[key]["psnr_corr"], p["L_mix"], p["oracle"]]
        ax.bar(range(len(arms)), vals, color=colors, edgecolor="k", width=0.72)
        for i, v in enumerate(vals):
            ax.text(i, v + 0.15, f"{v:.1f}", ha="center", fontsize=8)
        ax.axhline(p["oracle"], ls="--", c="k", lw=0.9, alpha=0.5)
        ax.set_xticks(range(len(arms))); ax.set_xticklabels(labels, fontsize=7.5)
        strong = " (fast-strong corner)" if (c["tc"] == 16 and c["cv"] == 0.40) else ""
        ax.set_title(f"t_c={c['tc']:.0f}, CV={c['cv']:.2f}{strong}", fontsize=9.5,
                     fontweight=("bold" if strong else "normal"))
        ax.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("scene PSNR (dB)")
    fig.suptitle("E1  Analytic E[a]=1 gauge fix closes the fast-strong-drift corner deficit with ZERO learning; "
                 "the learned corrector adds only marginal, non-robust extra   (EXPLORATORY, DEV-only)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fp = os.path.join(F, "e1_headtohead.png"); fig.savefig(fp, dpi=140); print("saved", fp)


# ---------------- E2: learned temporal prior, recon loop + adequacy ----------------
def fig_e2():
    rec = load("e2_recon_psnr.json")["rows"]
    mis = load("e2_ou_misfit.json")["rows"]
    lp = load("e2_learned_prior.json")["rows"]
    misd = {r["cls"]: r for r in mis}; lpd = {r["cls"]: r for r in lp}
    classes = [r["cls"] for r in rec]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 4.8))
    # left: recon-loop scene PSNR, OU-model vs learned vs oracle
    x = np.arange(len(classes)); w = 0.26
    ou = [r["psnr_ou"] for r in rec]; le = [r["psnr_lp"] for r in rec]; orc = [r["psnr_or"] for r in rec]
    ax1.bar(x - w, ou, w, label="OU-model DLGI", color="#d1495b", edgecolor="k")
    ax1.bar(x, le, w, label="learned-prior DLGI", color="#2e86ab", edgecolor="k")
    ax1.bar(x + w, orc, w, label="oracle (ceiling)", color="#111111", edgecolor="k")
    for i, r in enumerate(rec):
        d = r["dpsnr"]; ax1.text(i, max(r["psnr_ou"], r["psnr_lp"]) + 0.2,
                                 f"{d:+.2f} dB", ha="center", fontsize=8.5,
                                 color=("#1a7f3c" if d > 0.2 else ("#b03030" if d < -0.2 else "#666")),
                                 fontweight=("bold" if abs(d) > 0.5 else "normal"))
    ax1.set_xticks(x); ax1.set_xticklabels(classes); ax1.set_ylabel("scene PSNR (dB)")
    ax1.set_title("(A) recon-loop scene PSNR (learned − OU annotated)\nlearned prior wins where OU autocovariance is structurally wrong", fontsize=9.5)
    ax1.legend(fontsize=8, loc="lower right"); ax1.grid(axis="y", alpha=0.25); ax1.set_ylim(0, 21)
    # right: OU-model adequacy (innovation whiteness pass rate) OU vs learned predictor
    po = [misd[c]["frac_pass"] * 100 for c in classes]
    pl = [lpd[c]["pass_lp"] * 100 for c in classes]
    ax2.bar(x - w / 2, po, w, label="OU-Kalman innovations", color="#d1495b", edgecolor="k")
    ax2.bar(x + w / 2, pl, w, label="learned predictor innovations", color="#2e86ab", edgecolor="k")
    ax2.axhline(95, ls=":", c="k", alpha=0.5); ax2.text(3.35, 96, "adequate", fontsize=7, ha="right")
    for i, c in enumerate(classes):
        ax2.text(i - w / 2, po[i] + 1.5, f"{po[i]:.0f}", ha="center", fontsize=8)
        ax2.text(i + w / 2, pl[i] + 1.5, f"{pl[i]:.0f}", ha="center", fontsize=8)
    ax2.set_xticks(x); ax2.set_xticklabels(classes); ax2.set_ylabel("innovation-whiteness pass rate (%)")
    ax2.set_title("(B) model adequacy: learned prior restores whiteness on REGIME\n(where OU is structurally wrong); heavy-tail tell is kurtosis, not whiteness", fontsize=9.5)
    ax2.legend(fontsize=8, loc="lower center"); ax2.grid(axis="y", alpha=0.25); ax2.set_ylim(0, 108)
    fig.suptitle("E2  Learned temporal prior = domain extension beyond scalar OU   (EXPLORATORY, DEV-only)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(F, "e2_learned_prior.png"); fig.savefig(fp, dpi=140); print("saved", fp)


if __name__ == "__main__":
    fig_e1(); fig_e2()
