"""
FIG 3 - covariance-layer law and boundary map (the ONLY main figure for the
entire fog-DMD lifecycle).

Panels (2x2):
  (a) spectral law: B_p vs B_p (+) B_w next to the exact first-moment wall
  (b) empirical aperture validation: 18.76x in/out separation
  (c) the 81-cell Fisher atlas (REGION_EMPTY, mode-coverage bottleneck, pocket)
  (d) branch inset placeholder (both variants prepared)

Sources: T_fig3_aperture.json, T_fig3_p1.json, T_fig3_map.csv,
         T_fig3_map_meta.json, T_fig3_pocket.json
"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import common as K
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

AP = K.load_json("T_fig3_aperture.json")
P1 = K.load_json("T_fig3_p1.json")
PK = K.load_json("T_fig3_pocket.json")
MM = K.load_json("T_fig3_map_meta.json")
_, MAP = K.load_csv("T_fig3_map.csv")


# ---------------------------------------------------------------- (a) law
def panel_law(ax):
    wall = AP["mean_wall_duel_16x16"]
    T = wall["T_eff"]
    x = np.arange(len(T))
    w = 0.36
    ax.bar(x - w/2, wall["fresh_bandlimited_wall"], w, color=K.C["gray"],
           label="mean channel (fresh band-limited patterns)", zorder=3)
    ax.bar(x + w/2, wall["fixed_bank_mle"], w, color=K.C["blue"],
           label="covariance channel (fixed-bank)", zorder=3)
    ax.plot(x + w/2, wall["oracle_medium_known"], "o--", color=K.C["green"], ms=4,
            lw=1.0, label="oracle (medium known)", zorder=4)
    ax.axhline(1.0, color=K.C["red"], lw=1.0, ls=(0, (4, 2)), zorder=2)
    ax.text(x[-1] + 0.05, 1.0, "exact\nmean wall", ha="right", va="top",
            fontsize=6.2, color=K.C["red"], fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([str(t) for t in T])
    ax.set_xlabel("independent medium states  $T_{\\mathrm{eff}}$", fontsize=7)
    ax.set_ylabel("beyond-band error\n(1.0 = recovers nothing)", fontsize=7)
    ax.set_ylim(0, 1.15)
    ax.set_title("(a)  spectral law: covariance crosses the mean wall",
                 fontsize=7.6, loc="left", fontweight="bold")
    ax.legend(loc="lower left", frameon=False, fontsize=5.7, handlelength=1.2,
              labelspacing=0.25, borderpad=0.15)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# ------------------------------------------------------- (b) aperture 18.76x
def panel_separation(ax):
    sep = AP["separation_64x64"]
    ina = sep["in_all"]; outa = sep["out_all"]
    ax.scatter([0]*len(ina), ina, s=26, color=K.C["blue"], zorder=4,
               label="in-coverage modes")
    ax.scatter([1]*len(outa), outa, s=26, color=K.C["red"], zorder=4,
               label="out-of-coverage modes")
    ax.plot([-0.28, 0.28], [sep["in_coverage_err"]]*2, color=K.C["blue"], lw=1.6)
    ax.plot([0.72, 1.28], [sep["out_coverage_err"]]*2, color=K.C["red"], lw=1.6)
    ax.annotate("", (0.5, sep["out_coverage_err"]), (0.5, sep["in_coverage_err"]),
                arrowprops=dict(arrowstyle="<->", color=K.C["ink"], lw=1.0))
    ax.text(0.56, np.sqrt(sep["in_coverage_err"]*sep["out_coverage_err"]),
            "%.1f$\\times$\nseparation" % sep["separation"], ha="left", va="center",
            fontsize=8, fontweight="bold", color=K.C["ink"])
    ax.text(0, sep["in_coverage_err"]-0.35, "%.3f" % sep["in_coverage_err"],
            ha="center", va="top", fontsize=6.3, color=K.C["blue"])
    ax.text(1, sep["out_coverage_err"]+0.28, "%.3f" % sep["out_coverage_err"],
            ha="center", va="bottom", fontsize=6.3, color=K.C["red"])
    ax.set_xlim(-0.6, 1.9)
    ax.set_ylim(0, 5.3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["inside\naperture", "outside\naperture"], fontsize=6.8)
    ax.set_ylabel("reconstruction error", fontsize=7)
    ax.set_title("(b)  the aperture is real (64$\\times$64)", fontsize=7.6,
                 loc="left", fontweight="bold")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# ------------------------------------------------------ (c) 81-cell atlas
def panel_atlas(ax, cax):
    # rows = (sigma_f, shape) 9 combos ; cols = (k_w, claim) 9 combos
    sig = MM["axes"]["sigma_f"]; shp = MM["axes"]["shape"]
    kw = MM["axes"]["k_w_over_kp"]; cl = MM["axes"]["claim"]
    row_key = [(s, sh) for s in sig for sh in shp]      # 9
    col_key = [(k, c) for k in kw for c in cl]          # 9
    G = np.full((9, 9), np.nan)
    lut = {}
    for r in MAP:
        rk = (float(r["sigma_f"]), r["shape"])
        ck = (int(r["k_w_over_kp"]), float(r["claim"]))
        lut[(rk, ck)] = r
    for i, rk in enumerate(row_key):
        for j, ck in enumerate(col_key):
            G[i, j] = float(lut[(rk, ck)]["f_rec_witness"])
    P1_BAR = P1["bar"]["f_rec_snr3_required"]           # 0.70

    cmap = LinearSegmentedColormap.from_list(
        "frec", ["#F2F2F2", "#CDE3EF", "#7FB2D4", K.C["blue"], "#023E5A"])
    im = ax.imshow(G, cmap=cmap, vmin=0, vmax=P1_BAR, aspect="auto", origin="upper")
    # every cell fails: overlay a light 'empty' wash is implicit (none >= 0.70)

    # gridlines separating the 3x3 macro-blocks
    for p in (2.5, 5.5):
        ax.axhline(p, color="white", lw=1.4); ax.axvline(p, color="white", lw=1.4)

    # highlight the pocket cell (sigma_f=0.3, k^-2, k_w=1, claim=1.25)
    pr = row_key.index((0.3, "k^-2")); pc = col_key.index((1, 1.25))
    ax.add_patch(Rectangle((pc-0.5, pr-0.5), 1, 1, fill=False,
                           edgecolor=K.C["amber"], lw=2.4, zorder=5))
    ax.annotate("$1.25\\times$ pocket\n(sole residue)", (pc+0.4, pr-0.4),
                textcoords="offset points", xytext=(10, 12), ha="left", va="bottom",
                fontsize=6.0, color="#8A6200", fontweight="bold",
                bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.4),
                arrowprops=dict(arrowstyle="-", color=K.C["amber"], lw=1.0), zorder=7)

    ax.set_xticks(range(9))
    ax.set_xticklabels([("%g,\n%s" % (c, {1: "$k_w$1", 2: "$k_w$2", 4: "$k_w$4"}[k]))
                        for (k, c) in col_key], fontsize=5.0)
    ax.set_yticks(range(9))
    ax.set_yticklabels([("%.1f %s" % (s, sh)) for (s, sh) in row_key], fontsize=5.2)
    ax.set_xlabel("claim ($\\times k_p$)  $\\times$  medium band $k_w$", fontsize=6.6)
    ax.set_ylabel("$\\sigma_f$  $\\times$  spectrum shape", fontsize=6.6)
    ax.set_title("(c)  81-cell Fisher atlas: REGION_EMPTY, 0/81 pass",
                 fontsize=7.4, loc="left", fontweight="bold")
    ax.tick_params(length=0)

    cb = plt.colorbar(im, cax=cax)
    cb.set_label("mode coverage $f_{\\mathrm{rec}}$  (witness)", fontsize=6.3)
    cb.ax.tick_params(labelsize=5.6)
    cb.ax.axhline(P1_BAR, color=K.C["red"], lw=1.2)   # bar sits at top of scale
    cb.ax.text(1.7, P1_BAR, "0.70 bar\n(none pass)", transform=cb.ax.get_yaxis_transform(),
               fontsize=5.4, color=K.C["red"], va="center", ha="left")


# ------------------------------------------------------ (d) branch inset
def panel_branch(ax, variant="pending"):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.03, 0.06), 0.94, 0.88,
                 boxstyle="round,pad=0.01,rounding_size=0.03",
                 fill=False, edgecolor=K.C["gray"], lw=1.0, ls=(0, (4, 2))))
    ax.text(0.5, 0.88, "(d)  pocket branch", ha="center", va="top",
            fontsize=7.6, fontweight="bold")
    px, py = 0.24, 0.55
    if variant == "A":
        ax.plot(px, py, marker="D", ms=13, mfc=K.C["amber"], mec="white", mew=1.4)
        head = "POCKET_ACHIEVED"
        body = ("blind lifted-GLS realized the\npredeclared $1.25\\times$ pocket:\n"
                "one modest positive check\n(one reconstruction inset).")
    elif variant == "B":
        ax.plot(px, py, marker="D", ms=13, mfc="white", mec=K.C["amber"], mew=2.4)
        head = "ESTIMATOR_GAP_PERSISTS"
        body = ("Fisher-feasible but the frozen\nblind reconstruction did not\n"
                "realize it — an estimator-reach\ngap, not a contradiction.")
    else:
        ax.plot(px, py, marker="D", ms=13, mfc="white", mec=K.C["amber"], mew=2.4,
                alpha=0.85)
        head = "BRANCH  =  PENDING"
        body = ("blind demo verdict arrives\nseparately.  A = achieved (filled);\n"
                "B = unrealized (hollow).\nBoth variants prepared.")
    ax.text(0.44, py+0.16, head, ha="left", va="center", fontsize=6.8,
            fontweight="bold", color="#8A6200")
    ax.text(0.44, py-0.09, body, ha="left", va="top", fontsize=5.9,
            color=K.C["ink"], linespacing=1.15)
    ax.text(0.5, 0.14, "pocket: $f_{\\mathrm{rec}}$(nat) %.2f  ·  NRMSE %.2f  ·  "
            "$T_{\\mathrm{req}}\\!\\approx\\!%d$"
            % (PK["f_rec_nat"], PK["nmse_crb_nat_med"], int(PK["t_req_0p30"])),
            ha="center", va="center", fontsize=5.9, color=K.C["line"])


def build(variant="pending"):
    K.apply_style()
    fig = plt.figure(figsize=(K.cm(17.8), K.cm(12.6)))
    gs = GridSpec(2, 2, width_ratios=[1.0, 1.05], height_ratios=[1.0, 1.12],
                  wspace=0.30, hspace=0.46, left=0.085, right=0.945,
                  top=0.895, bottom=0.075)
    panel_law(fig.add_subplot(gs[0, 0]))
    panel_separation(fig.add_subplot(gs[0, 1]))
    # atlas with an attached colorbar axis
    axc = fig.add_subplot(gs[1, 0])
    pos = axc.get_position()
    cax = fig.add_axes([pos.x1 + 0.012, pos.y0 + 0.02, 0.012, pos.height - 0.14])
    panel_atlas(axc, cax)
    panel_branch(fig.add_subplot(gs[1, 1]), variant)
    fig.suptitle("The covariance layer: broad formal support, thin usable aperture",
                 x=0.085, y=0.965, ha="left", fontsize=9.5, fontweight="bold")
    return fig


if __name__ == "__main__":
    for v, stem in [("pending", "fig3"), ("A", "fig3_branchA"), ("B", "fig3_branchB")]:
        fig = build(v)
        K.save(fig, stem)
        plt.close(fig)
    print("FIG3 done.")
