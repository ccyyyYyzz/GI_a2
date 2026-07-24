"""
FIGURE 1 (R43 §6.2) — HERO / THE GATE: rank collapses; jets survive.

A schematic (no fitted data): a horizontal disorder axis from thin multiplicative
screen to complete scrambling, three rows —
  mean spatial support        B_p  ->  {0}
  reconstruction score rank   many ->   1
  quotient jet order          finite m=1/2  ->  finite m=1/2   (SURVIVES)
— with the scrambling-endpoint law  Cov(b)=Q(x)(O o O)+R, Q=x^T G x, the two
directional classes, and the amplitude-anchor condition as a visible icon.

Ten-second conclusion: disorder erases coordinates, not local change observability.

The few numeric anchors printed (rank-one 99.99%, jet orders 2/4) are read from the
COMMITTED FROZEN artifacts SCRAMBLE_RESULTS.json@ed7a1e0 and JET_TEST.json@1bf29f1
(extracted to _frozen_sources/); no data generation.
"""
import json, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, Ellipse, FancyBboxPatch, Arc
from prl_style import apply_style, save, COL2, C_GENERIC, C_ORTHO, C_MIXED, C_ANCHOR, C_WALL, OI

HERE = os.path.dirname(os.path.abspath(__file__))
FS = os.path.join(HERE, "_frozen_sources")
scr = json.load(open(os.path.join(FS, "SCRAMBLE_RESULTS.committed.json")))
jet = json.load(open(os.path.join(FS, "JET_TEST.committed.json")))
rank1_pct = 100.0 * scr["valB_rank1"]["frac_cov_energy_in_OhO_direction"]   # 99.99
apply_style()

fig = plt.figure(figsize=(COL2, 3.5))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

# regime background bands (thin screen vs complete scrambling) -- drawn first
ax.add_patch(FancyBboxPatch((0.225, 0.135), 0.185, 0.685, boxstyle="round,pad=0.006",
             facecolor="0.955", edgecolor="none", zorder=0))
ax.add_patch(FancyBboxPatch((0.445, 0.135), 0.205, 0.685, boxstyle="round,pad=0.006",
             facecolor=OI["blue"], alpha=0.065, edgecolor="none", zorder=0))

# lane geometry
Y = {"mean": 0.685, "rank": 0.470, "jet": 0.255}
XL, XR = 0.315, 0.545            # left-state / right-state centers
LAB = 0.015                      # row-label left edge
DEAD = "0.25"                    # neutral-dark = collapsed / lost channel

def lane_label(y, txt):
    ax.text(LAB, y, txt, ha="left", va="center", fontsize=7.3, fontweight="bold")

def arrow(x0, x1, y, color=C_WALL, lw=1.6, style="-|>"):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle=style,
                 mutation_scale=11, lw=lw, color=color, shrinkA=0, shrinkB=0))

# ---------------- top disorder axis ----------------
ax.add_patch(FancyArrowPatch((0.17, 0.935), (0.62, 0.935), arrowstyle="-|>",
             mutation_scale=16, lw=3.2, color="0.35"))
ax.text(0.395, 0.965, "OPTICAL DISORDER", ha="center", va="bottom",
        fontsize=8.2, fontweight="bold", color="0.25")
ax.text(XL, 0.885, "thin multiplicative\nscreen", ha="center", va="center", fontsize=7.0)
ax.text(XR, 0.885, "complete\nscrambling", ha="center", va="center", fontsize=7.0,
        fontweight="bold")

# ================= ROW 1 : mean spatial support =================
lane_label(Y["mean"], "mean channel\nspatial support")
# left: band-limited support disk B_p
ax.add_patch(Circle((XL, Y["mean"]), 0.052, facecolor=OI["skyblue"], alpha=0.55,
             edgecolor=OI["blue"], lw=1.0))
ax.text(XL, Y["mean"], r"$B_p$", ha="center", va="center", fontsize=7.5, color=OI["blue"])
arrow(0.385, 0.475, Y["mean"])
# right: collapses to {0} (single DC point) inside an X'd-out ghost disk
ax.add_patch(Circle((XR, Y["mean"]), 0.052, facecolor="none",
             edgecolor="0.7", lw=0.8, ls=(0, (2, 2))))
ax.plot([XR - 0.037, XR + 0.037], [Y["mean"] - 0.037, Y["mean"] + 0.037], color="0.7", lw=0.8)
ax.plot([XR - 0.037, XR + 0.037], [Y["mean"] + 0.037, Y["mean"] - 0.037], color="0.7", lw=0.8)
ax.add_patch(Circle((XR, Y["mean"]), 0.009, facecolor=DEAD, edgecolor="none"))
ax.text(XR, Y["mean"] + 0.075, r"$\to\{0\}$", ha="center", va="bottom",
        fontsize=8.0, color=DEAD, fontweight="bold")
ax.text(XR, Y["mean"] - 0.075, "DC only", ha="center", va="top",
        fontsize=6.6, color="0.4")

# ================= ROW 2 : reconstruction Fisher rank =================
lane_label(Y["rank"], "reconstruction\nFisher rank")
# left: many-mode ellipsoid with several principal axes
ax.add_patch(Ellipse((XL, Y["rank"]), 0.115, 0.075, angle=25,
             facecolor=OI["green"], alpha=0.28, edgecolor=OI["green"], lw=1.0))
for ang in (25, 70, 115, 160):
    dx, dy = 0.052 * np.cos(np.radians(ang)), 0.040 * np.sin(np.radians(ang))
    ax.plot([XL - dx, XL + dx], [Y["rank"] - dy, Y["rank"] + dy], color=OI["green"], lw=0.9)
ax.text(XL, Y["rank"] - 0.078, "rank many", ha="center", va="top", fontsize=6.8,
        color=OI["green"])
arrow(0.385, 0.475, Y["rank"])
# right: rank-one degenerate line (single direction Gx) -- imaging collapsed (dark)
ax.plot([XR - 0.060, XR + 0.060], [Y["rank"] - 0.030, Y["rank"] + 0.030],
        color=DEAD, lw=3.0, solid_capstyle="round")
ax.text(XR, Y["rank"] - 0.070, r"rank $1$  ($Gx$)", ha="center", va="top",
        fontsize=7.0, color=DEAD, fontweight="bold")
ax.text(XR, Y["rank"] - 0.115, f"{rank1_pct:.2f}% energy in one mode",
        ha="center", va="top", fontsize=6.6, color="0.45")

# ================= ROW 3 : quotient jet order (SURVIVES) =================
lane_label(Y["jet"], "local change\nquotient jet order")

def slope_glyph(cx, cy, boxed=False):
    w, h = 0.058, 0.058
    x0, y0 = cx - w / 2, cy - h / 2
    # m=1 (blue) and m=2 (vermillion) mini log-log lines
    ax.plot([x0 + 0.006, x0 + w - 0.006], [y0 + 0.010, y0 + 0.028],
            color=C_GENERIC, lw=1.6)
    ax.plot([x0 + 0.006, x0 + w - 0.006], [y0 + 0.006, y0 + h - 0.006],
            color=C_ORTHO, lw=1.6)
    ax.text(x0 + w + 0.006, y0 + 0.014, r"$m{=}1$", fontsize=6.6, color=C_GENERIC,
            va="center")
    ax.text(x0 + w + 0.006, y0 + h - 0.012, r"$m{=}2$", fontsize=6.6, color=C_ORTHO,
            va="center")
    if boxed:
        ax.add_patch(FancyBboxPatch((x0 - 0.012, y0 - 0.014), w + 0.062, h + 0.028,
                     boxstyle="round,pad=0.004", facecolor=OI["green"], alpha=0.14,
                     edgecolor=OI["green"], lw=1.3))

slope_glyph(XL, Y["jet"])
ax.text(XL, Y["jet"] - 0.058, "finite", ha="center", va="top", fontsize=6.8, color="0.3")
arrow(0.385, 0.475, Y["jet"], color=OI["green"], lw=2.2)
slope_glyph(XR, Y["jet"], boxed=True)
ax.text(XR + 0.016, Y["jet"] - 0.058, "finite — SURVIVES", ha="center", va="top",
        fontsize=6.9, color=OI["green"], fontweight="bold")

# ================= right endpoint law callout =================
bx0, by0, bw, bh = 0.715, 0.135, 0.278, 0.72
ax.add_patch(FancyBboxPatch((bx0, by0), bw, bh, boxstyle="round,pad=0.008",
             facecolor="white", edgecolor=OI["blue"], lw=1.2, alpha=0.95))
cx = bx0 + 0.014
ax.text(bx0 + bw / 2, by0 + bh - 0.030, "complete-scrambling law",
        ha="center", va="center", fontsize=7.4, fontweight="bold", color=OI["blue"])
ax.plot([cx, bx0 + bw - 0.014], [by0 + bh - 0.058] * 2, color=OI["blue"], lw=1.0, alpha=0.6)
ax.text(cx, by0 + bh - 0.10, r"$\mathrm{Cov}(b)=Q(x)\,(O\!\circ\!O)+R$",
        ha="left", va="center", fontsize=7.6)
ax.text(cx, by0 + bh - 0.155, r"$Q(x)=x^{\top} G\, x,\quad G>0$",
        ha="left", va="center", fontsize=7.6)
ax.plot([cx, bx0 + bw - 0.014], [by0 + bh - 0.195] * 2, color="0.7", lw=0.6)
ax.text(cx, by0 + bh - 0.235, "two directional classes:", ha="left", va="center",
        fontsize=6.6, style="italic", color="0.3")
ax.text(cx, by0 + bh - 0.295, r"$x^{\top}\!G\delta\neq0\;\Rightarrow\; m{=}1$",
        ha="left", va="center", fontsize=7.2, color=C_GENERIC)
ax.text(cx + 0.020, by0 + bh - 0.335, "first-order observable", ha="left", va="center",
        fontsize=6.6, color=C_GENERIC)
ax.text(cx, by0 + bh - 0.395, r"$x^{\top}\!G\delta=0\;\Rightarrow\; m{=}2$",
        ha="left", va="center", fontsize=7.2, color=C_ORTHO)
ax.text(cx + 0.020, by0 + bh - 0.435, r"curvature-rescued ($\delta^{\top}\!G\delta>0$)",
        ha="left", va="center", fontsize=6.6, color=C_ORTHO)
ax.add_patch(FancyBboxPatch((cx + 0.028, by0 + bh - 0.560), bw - 0.070, 0.052,
             boxstyle="round,pad=0.004", facecolor=OI["green"], alpha=0.14,
             edgecolor=OI["green"], lw=1.0))
ax.text(bx0 + bw / 2, by0 + bh - 0.534, r"$\Rightarrow$ blind set empty",
        ha="center", va="center", fontsize=7.4, fontweight="bold", color=OI["green"])

# amplitude-anchor icon (drawn) + label
ax_c = bx0 + 0.045; ax_y = by0 + 0.080
ax.add_patch(Circle((ax_c, ax_y + 0.030), 0.010, facecolor="none",
             edgecolor=OI["orange"], lw=1.6))                       # ring
ax.plot([ax_c, ax_c], [ax_y - 0.028, ax_y + 0.020], color=OI["orange"], lw=1.8)  # shank
ax.plot([ax_c - 0.016, ax_c + 0.016], [ax_y + 0.006, ax_y + 0.006],
        color=OI["orange"], lw=1.8)                                  # stock
ax.add_patch(Arc((ax_c, ax_y - 0.020), 0.044, 0.030, angle=0, theta1=200, theta2=340,
             color=OI["orange"], lw=1.8))                            # flukes
ax.text(ax_c + 0.028, ax_y, "amplitude anchor required\n(else quotient info $J\\to0$)",
        ha="left", va="center", fontsize=6.6, color=OI["orange"])

# ================= ten-second conclusion banner =================
ax.add_patch(FancyBboxPatch((0.02, 0.02), 0.96, 0.072, boxstyle="round,pad=0.006",
             facecolor=OI["blue"], alpha=0.08, edgecolor=OI["blue"], lw=0.9))
ax.text(0.5, 0.056,
        r"Strong disorder erases $\mathbf{coordinates}$ (support $\to\{0\}$, rank $\to 1$) "
        r"but not $\mathbf{local\ change\ observability}$: the quotient jet order stays finite.",
        ha="center", va="center", fontsize=7.4, color="0.1")

save(fig, os.path.join(HERE, "fig1_hero"))
