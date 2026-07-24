"""
FIG 1 (THE GATE) - master information atlas of a bucket record.

Two rows (mean channel / covariance channel) on a shared spatial-frequency axis.
  pale fill  = formal support        (mean: B_p ; covariance: B_p (+) B_w = 2x)
  dark fill  = Fisher-usable region  (mean: full band ; covariance: thin rim)
  green dot  = certified mean operating point
  red cross  = killed 1.8x covariance claim
  amber mark = 1.25x natural-scene pocket  (filled=achieved / hollow=unrealized)
Legend line: "support is not supply".

Must be intelligible in 10 s WITHOUT a caption: a MAP, not a method anthology.
All geometry read from source_tables/T_fig1_atlas.json.
"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import common as K
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.lines import Line2D

A = K.load_json("T_fig1_atlas.json")
M = A["rows"]["mean"]
V = A["rows"]["covariance"]
mk = {m["id"]: m for m in A["markers"]}

XMAX = 2.18
ROW = {"mean": 1.0, "covariance": 0.0}   # y-centres
BH  = 0.34                                # half band height


def band(ax, y, x0, x1, face, z, alpha=1.0, edge="none", lw=0):
    ax.add_patch(Rectangle((x0, y - BH), x1 - x0, 2 * BH, facecolor=face,
                           edgecolor=edge, linewidth=lw, alpha=alpha, zorder=z))


def build(variant="pending"):
    K.apply_style()
    fig, ax = plt.subplots(figsize=(K.cm(18.0), K.cm(8.2)))

    # ---- row bands -------------------------------------------------------
    for key in ("mean", "covariance"):
        y = ROW[key]
        r = A["rows"][key]
        band(ax, y, 0, r["formal_support_extent"], K.C["grayfill"], z=1)   # pale support
        band(ax, y, 0, r["fisher_usable_extent"], K.C["darkfill"], z=2)    # dark usable
        # bright edge on the dark boundary so the usable extent is legible
        ax.plot([r["fisher_usable_extent"]] * 2, [y - BH, y + BH],
                color="white", lw=1.4, zorder=3)

    # covariance: hatch the pale-only "support without supply" gap
    band(ax, ROW["covariance"], V["fisher_usable_extent"], V["formal_support_extent"],
         "none", z=3)
    ax.add_patch(Rectangle((V["fisher_usable_extent"], ROW["covariance"] - BH),
                           V["formal_support_extent"] - V["fisher_usable_extent"], 2 * BH,
                           facecolor="none", edgecolor=K.C["gray"], linewidth=0.9,
                           linestyle=(0, (4, 2)), hatch="////", zorder=3))

    # ---- the exact mean wall at k_p = 1 ---------------------------------
    wall_top = ROW["mean"] + BH + 0.16
    wall_bot = ROW["covariance"] - BH - 0.16
    ax.plot([1.0, 1.0], [wall_bot, wall_top], color=K.C["ink"], lw=1.1,
            ls=(0, (5, 2)), zorder=6)

    # ---- per-row verdict tags (the 10-second read) ----------------------
    ax.text(0.5, ROW["mean"] + BH + 0.10, "support $\\approx$ supply",
            ha="center", va="bottom", fontsize=8.4, fontweight="bold",
            color=K.C["green"])
    ax.text(1.55, ROW["covariance"] + BH + 0.12, "support $\\gg$ supply",
            ha="center", va="bottom", fontsize=8.4, fontweight="bold",
            color=K.C["red"])

    # ---- in-band fill labels --------------------------------------------
    ax.text(0.5, ROW["mean"], "Fisher-usable\n(full band)", ha="center", va="center",
            color="white", fontsize=7.4, fontweight="bold", zorder=5, linespacing=0.95)
    ax.text(0.5, ROW["covariance"], "usable", ha="center", va="center",
            color="white", fontsize=7.2, fontweight="bold", zorder=5)
    ax.text((V["fisher_usable_extent"] + V["formal_support_extent"]) / 2,
            ROW["covariance"], "formal support,\nno usable supply", ha="center",
            va="center", color=K.C["line"], fontsize=6.9, style="italic",
            zorder=5, linespacing=0.95)

    # ---- support-band captions ------------------------------------------
    ax.text(0.5, ROW["mean"] - BH - 0.10, "formal support = modulator band $B_p$",
            ha="center", va="top", fontsize=6.8, color=K.C["line"])
    ybar = ROW["covariance"] - BH - 1.02
    ax.annotate("", (V["formal_support_extent"], ybar), (0.0, ybar),
                arrowprops=dict(arrowstyle="<->", color=K.C["line"], lw=0.8))
    ax.text(V["formal_support_extent"] / 2, ybar - 0.06,
            "formal support = $B_p \\oplus B_w$  (2$\\times$ aperture)", ha="center",
            va="top", fontsize=6.8, color=K.C["line"])

    # ---- markers with leader callouts -----------------------------------
    gx = 0.70
    ax.plot(gx, ROW["mean"], "o", ms=8.5, mfc=K.C["green"], mec="white", mew=1.1, zorder=8)
    ax.annotate("certified mean\noperating point", (gx, ROW["mean"]),
                (0.20, ROW["mean"] + BH + 0.42), ha="center", va="bottom",
                fontsize=6.9, color=K.C["green"], fontweight="bold",
                arrowprops=dict(arrowstyle="-", color=K.C["green"], lw=0.8,
                                connectionstyle="arc3,rad=-0.2"), zorder=9)

    # both covariance markers callout BELOW the row (clear of everything)
    rx = mk["killed_1p8"]["x_kp"]
    ax.plot(rx, ROW["covariance"], marker="X", ms=10, mfc=K.C["red"], mec="white",
            mew=1.1, zorder=8, linestyle="none")
    ax.annotate("killed $1.8\\times$\ncovariance claim", (rx, ROW["covariance"] - BH),
                (rx, ROW["covariance"] - BH - 0.18), ha="center", va="top",
                fontsize=6.9, color=K.C["red"], fontweight="bold",
                arrowprops=dict(arrowstyle="-", color=K.C["red"], lw=0.8), zorder=9)

    px = mk["pocket_1p25"]["x_kp"]
    if variant == "A":
        ax.plot(px, ROW["covariance"], marker="D", ms=8, mfc=K.C["amber"],
                mec="white", mew=1.0, zorder=9, linestyle="none")
        ptxt = "$1.25\\times$ pocket\n(achieved)"
    elif variant == "B":
        ax.plot(px, ROW["covariance"], marker="D", ms=8, mfc="white",
                mec=K.C["amber"], mew=1.8, zorder=9, linestyle="none")
        ptxt = "$1.25\\times$ pocket\n(feasible, unrealized)"
    else:
        ax.plot(px, ROW["covariance"], marker="D", ms=8, mfc="white",
                mec=K.C["amber"], mew=1.8, zorder=9, linestyle="none")
        ptxt = "$1.25\\times$ pocket\n(pending)"
    ax.annotate(ptxt, (px, ROW["covariance"] - BH),
                (px - 0.18, ROW["covariance"] - BH - 0.18), ha="center", va="top",
                fontsize=6.7, color="#8A6200", fontweight="bold",
                arrowprops=dict(arrowstyle="-", color=K.C["amber"], lw=0.9), zorder=9)

    # ---- row labels ------------------------------------------------------
    ax.text(-0.07, ROW["mean"], "mean\nchannel", ha="right", va="center",
            fontsize=9, fontweight="bold", color=K.C["ink"], linespacing=0.95)
    ax.text(-0.07, ROW["covariance"], "covariance\nchannel", ha="right", va="center",
            fontsize=9, fontweight="bold", color=K.C["ink"], linespacing=0.95)

    # ---- axes cosmetics --------------------------------------------------
    ax.set_xlim(-0.46, XMAX)
    ax.set_ylim(ROW["covariance"] - BH - 1.42, ROW["mean"] + BH + 0.98)
    ax.set_xticks([0, 0.5, 1.0, 1.25, 1.5, 1.8, 2.0])
    ax.set_xticklabels(["0", "0.5", "$k_p$", "1.25", "1.5", "1.8", "$2k_p$"])
    ax.set_yticks([])
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_position(("data", ROW["covariance"] - BH - 1.28))
    ax.set_xlabel("spatial frequency  /  task mode   (units of modulator band edge $k_p$)",
                  fontsize=7.8, labelpad=2)
    ax.tick_params(axis="x", length=3)
    ax.text(1.0, wall_top + 0.02, "exact first-moment wall", ha="center", va="bottom",
            fontsize=6.6, color=K.C["ink"], style="italic")

    # ---- title + punchline ----------------------------------------------
    ax.text(-0.46, ROW["mean"] + BH + 0.74,
            "One bucket record, two information apertures",
            ha="left", va="bottom", fontsize=11, fontweight="bold")
    ax.text(XMAX, ROW["mean"] + BH + 0.74, "support is not supply",
            ha="right", va="bottom", fontsize=11.5, fontstyle="italic",
            fontweight="bold", color=K.C["ink"])

    # ---- legend: horizontal strip below the axis ------------------------
    handles = [
        Rectangle((0, 0), 1, 1, fc=K.C["grayfill"], ec=K.C["gray"], lw=0.6),
        Rectangle((0, 0), 1, 1, fc=K.C["darkfill"], ec="none"),
        Line2D([0], [0], marker="o", color="none", mfc=K.C["green"], mec="white", ms=8, mew=1.0),
        Line2D([0], [0], marker="X", color="none", mfc=K.C["red"], mec="white", ms=8.5, mew=1.0),
        Line2D([0], [0], marker="D", color="none", mfc="white", mec=K.C["amber"], ms=7, mew=1.5),
    ]
    labels = ["formal support", "Fisher-usable region", "certified operating point",
              "killed $1.8\\times$ claim", "$1.25\\times$ pocket"]
    leg = ax.legend(handles, labels, loc="upper center",
                    bbox_to_anchor=(0.5, -0.14), ncol=5, frameon=False,
                    handlelength=1.3, columnspacing=1.4, handletextpad=0.5)

    fig.subplots_adjust(left=0.11, right=0.99, top=0.88, bottom=0.20)
    return fig


if __name__ == "__main__":
    for v, stem in [("pending", "fig1"), ("A", "fig1_pocketA"), ("B", "fig1_pocketB")]:
        fig = build(v)
        K.save(fig, stem)
        plt.close(fig)
    print("FIG1 done.")
