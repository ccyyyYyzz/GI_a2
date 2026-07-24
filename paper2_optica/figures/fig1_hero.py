"""
PAPER 2 (Optica) — FIGURE 1 / HERO: "What the sensor cannot see".

Architecture FROZEN by R45 §C3 (TLSG_PASS). Three linked panels along the causal
spine of the paper:

  SYMMETRY CREATES A WALL  ->  CALIBRATION MEASURES ITS BLINDNESS CAPACITY
                           ->  FINITE TANGENT DIMENSION SETS WHAT CAN BE CERTIFIED.

  LEFT  — two exact symmetry walls (score-null orbits) + their one controlled breaker each:
            support/field wall  6.6e-16  (breaker: finite-window leak law)
            energy/unitary wall 2.18e-16 (breaker: NA clipping, eps^1)
  CENTER— blindness capacity: generalized leak sigma_d vs code-subspace dimension d
            (IT4b joint spectrum), with the single-plane (80 dims <=1e-5), joint
            (35 dims <=1e-4) and DPSS (<=2.4e-4, 80->24, noncomposable) marks.
  RIGHT — certifiability ladder: required T*lambda ~ 1 : sqrt(p) : p (known/blind/estimation),
            TLSG-measured p-exponents -0.004/0.444/1.010; KT6 below the blind line;
            the projected jet (restored m=2) above its known-direction channel line.

Ten-second message: an optical sensor has EXACT things it cannot see, that blindness
has a MEASURED capacity that does not automatically compose, and a nonzero physical
signal is only CERTIFIABLE above a dimension-dependent threshold.

DATA: committed frozen artifacts only, extracted verbatim to _frozen_sources/.
  KT1 @19338a4, IT1 @f13a323, IT4b @5666d1e, IT5 @c8a0487, KT6 @635ff05, TLSG @35d858e.
No data generation.
"""
import json
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle
from o2_style import (apply_style, save, COL2, OI, C_WALL, C_BREAK, C_CAP, C_DPSS,
                      C_KNOWN, C_BLIND, C_EST, C_KT6, C_JET, C_GREY)

HERE = os.path.dirname(os.path.abspath(__file__))
FS = os.path.join(HERE, "_frozen_sources")


def load(name):
    return json.load(open(os.path.join(FS, name)))


kt1 = load("KT1_shell_sweep.committed.json")
it1 = load("IT1_two_wall_jet_invariance.committed.json")
it4b = load("IT4b_spectrum_trade.committed.json")
it5 = load("IT5_dpss_concentration.committed.json")
kt6 = load("KT6_channel_ratio_ranging.committed.json")
tlsg = load("TLSG_RESULTS.committed.json")

# ---- pull the frozen numbers (verbatim from committed JSONs) ----
FIELD_WALL = kt1["verdicts"]["P2a_field_wall"]["per_variant"]["V6_kR5"]["max_field_leak_beyond_2kR"]  # 6.59e-16
SUPPORT_BREAK = kt1["verdicts"]["P3"]["value"]                     # 5.3e-2 finite-window detector leak
ENERGY_WALL = it1["verdict"]["max_full_aperture_null"]            # 2.18e-16
CLIP_SLOPE = float(np.mean(it1["verdict"]["clip_slopes"]))        # ~0.996 (eps^1)

spec = it4b["results"]["joint_spectrum"]["rel_leak_vs_d"]
d_vals = np.array(sorted(int(k) for k in spec))
leak_vals = np.array([spec[str(d)] for d in d_vals])
D_JOINT_1E4 = it4b["results"]["joint_spectrum"]["d_at_1e4"]       # 35
D_JOINT_1E5 = it4b["results"]["joint_spectrum"]["d_at_1e5"]       # 8
D_SINGLE_1E5 = it4b["results"]["per_z1_vs_joint"]["z1_0_only_d_at_1e5"]  # 80
DPSS_LEAK = it5["verdict"]["dpss_max_witness_leak"]              # 2.38e-4
DPSS_NOWIN = it5["combined_single_z1_SVD"]["no_window_d_at_1e4"]  # 80
DPSS_WIN = it5["combined_single_z1_SVD"]["dpss_d_at_1e4"]         # 24

# TLSG ledger exponents (parse "measured" string -0.004 / 0.444 / 1.010)
exp_str = next(b for b in tlsg["bars"] if b["bar"] == 4)["measured"]
EXPS = [float(x) for x in re.findall(r"-?\d+\.\d+", exp_str)]       # [-0.004, 0.444, 1.010]
bar5 = next(b for b in tlsg["bars"] if b["bar"] == 5)["measured"]
KT6_TLAM = float(re.search(r"T\*lambda=([0-9.]+)", bar5).group(1))  # 0.0866
KT6_P = int(re.search(r"p=(\d+)", bar5).group(1))                   # 136
KT6_SIGFLOOR = kt6["cells"][0]["median_signal_to_floor_near"]      # ~5.4e-3

apply_style()
fig = plt.figure(figsize=(COL2, 3.45))
# top spine banner + three panels
axS = fig.add_axes([0.0, 0.90, 1.0, 0.10]); axS.axis("off")
axL = fig.add_axes([0.045, 0.11, 0.265, 0.70])
axC = fig.add_axes([0.400, 0.11, 0.255, 0.70])
axR = fig.add_axes([0.740, 0.11, 0.245, 0.70])

# ================= TOP SPINE =================
axS.set_xlim(0, 1); axS.set_ylim(0, 1)
spine = [("symmetry\ncreates a wall", 0.135, C_WALL),
         ("calibration measures\nits blindness capacity", 0.500, C_CAP),
         ("tangent dimension sets\nwhat is certifiable", 0.865, C_EST)]
for i, (txt, x, c) in enumerate(spine):
    axS.text(x, 0.55, txt, ha="center", va="center", fontsize=6.6, fontweight="bold",
             color=c, linespacing=0.95)
    if i < 2:
        axS.add_patch(FancyArrowPatch((x + 0.115, 0.55), (spine[i + 1][1] - 0.135, 0.55),
                      arrowstyle="-|>", mutation_scale=9, lw=1.4, color="0.45"))
axS.text(0.5, 0.02, "WHAT AN OPTICAL SENSOR CANNOT SEE", ha="center", va="bottom",
         fontsize=8.6, fontweight="bold", color="0.12")

# ================= LEFT : symmetry walls =================
axL.set_xlim(0, 1); axL.set_ylim(0, 1); axL.axis("off")
axL.text(0.5, 1.05, "Two exact walls", ha="center", va="bottom", fontsize=7.8,
         fontweight="bold", color=C_WALL, transform=axL.transAxes)


def wall_row(y, name, invar, wall_val, brk_label, brk_val):
    # WALL block (machine-zero, blue)
    axL.add_patch(FancyBboxPatch((0.02, y - 0.075), 0.44, 0.15,
                  boxstyle="round,pad=0.006", facecolor=OI["skyblue"], alpha=0.22,
                  edgecolor=C_WALL, lw=1.1))
    axL.text(0.24, y + 0.038, name, ha="center", va="center", fontsize=6.9,
             fontweight="bold", color=C_WALL)
    axL.text(0.24, y - 0.006, invar, ha="center", va="center", fontsize=5.7, color="0.35")
    axL.text(0.24, y - 0.050, f"null {wall_val:.1e}".replace("e-16", r"$\times10^{-16}$"),
             ha="center", va="center", fontsize=6.4, fontweight="bold", color="0.15")
    # breaker arrow -> vermillion
    axL.add_patch(FancyArrowPatch((0.47, y), (0.565, y), arrowstyle="-|>",
                  mutation_scale=8, lw=1.3, color=C_BREAK))
    axL.text(0.515, y + 0.052, "breaker", ha="center", va="bottom", fontsize=5.3,
             style="italic", color=C_BREAK)
    axL.add_patch(FancyBboxPatch((0.575, y - 0.075), 0.40, 0.15,
                  boxstyle="round,pad=0.006", facecolor=OI["vermil"], alpha=0.14,
                  edgecolor=C_BREAK, lw=1.0))
    axL.text(0.775, y + 0.030, brk_label, ha="center", va="center", fontsize=5.9,
             color=C_BREAK, linespacing=0.9)
    axL.text(0.775, y - 0.048, brk_val, ha="center", va="center", fontsize=6.2,
             fontweight="bold", color="0.15")


wall_row(0.80, "SUPPORT wall", r"hard-pupil field, $\ker A$", FIELD_WALL,
         "finite window\n(leak law)", f"{SUPPORT_BREAK*100:.0f}%")
wall_row(0.52, "ENERGY wall", r"$E\!\to\!UE$, full bucket", ENERGY_WALL,
         "NA clipping", rf"$\varepsilon^{{{CLIP_SLOPE:.2f}}}$")
# tie: score-null orbit label
axL.text(0.24, 0.335, "score-null orbits", ha="center", va="center", fontsize=5.8,
         style="italic", color="0.4")
axL.text(0.5, 0.245,
         r"$\mathbf{P_{g\cdot x}\!=\!P_x\;\Rightarrow\;I_x X_\xi\!=\!0}$",
         ha="center", va="center", fontsize=6.8, color="0.12", transform=axL.transAxes)
axL.add_patch(FancyBboxPatch((0.03, 0.02), 0.94, 0.135, boxstyle="round,pad=0.006",
              facecolor="0.955", edgecolor="none"))
axL.text(0.5, 0.088, "no data volume, moment order, or\nprior recovers an orbit direction",
         ha="center", va="center", fontsize=5.7, color="0.3", linespacing=0.95)

# ================= CENTER : blindness capacity =================
axC.set_yscale("log")
axC.plot(d_vals, leak_vals, "-o", color=C_CAP, lw=1.4, ms=2.6, zorder=5,
         label="joint spectrum")
# threshold guides
for thr, lab in [(1e-4, r"$10^{-4}$"), (1e-5, r"$10^{-5}$")]:
    axC.axhline(thr, color="0.6", lw=0.7, ls=(0, (3, 2)))
# joint marks
axC.plot([D_JOINT_1E4], [1e-4], marker="v", color=C_CAP, ms=6, zorder=7)
axC.annotate(f"joint {D_JOINT_1E4}\n@$10^{{-4}}$", (D_JOINT_1E4, 1e-4),
             xytext=(D_JOINT_1E4 + 6, 3.0e-4), fontsize=5.9, color=C_CAP, ha="left")
axC.plot([D_JOINT_1E5], [1e-5], marker="v", color=C_CAP, ms=5, zorder=7)
# single-plane capacity band (80 dims <= 1e-5) — annotated reference
axC.axvline(D_SINGLE_1E5, color=OI["blue"], lw=1.0, ls=(0, (4, 2)), alpha=0.8)
axC.annotate(f"single conjugate\nplane: {D_SINGLE_1E5} dims\n"r"$\leq10^{-5}$",
             (D_SINGLE_1E5, 1.3e-6), xytext=(D_SINGLE_1E5 - 2, 1.5e-6),
             fontsize=5.9, color=OI["blue"], ha="right", va="bottom", fontweight="bold")
# DPSS route (noncomposable): 80->24, leak <= 2.4e-4
axC.annotate("", xy=(DPSS_WIN, DPSS_LEAK), xytext=(DPSS_NOWIN, DPSS_LEAK),
             arrowprops=dict(arrowstyle="-|>", color=C_DPSS, lw=1.3))
axC.plot([DPSS_WIN], [DPSS_LEAK], marker="s", color=C_DPSS, ms=4.5, zorder=7)
axC.text(6, 1.4e-2,
         f"DPSS stack {DPSS_NOWIN}$\\to${DPSS_WIN}\n"rf"leak $\leq{DPSS_LEAK*1e4:.1f}\times10^{{-4}}$"
         "\n(noncomposable)",
         fontsize=5.7, color=C_DPSS, ha="left", va="center", fontweight="bold")
axC.set_xlim(0, 122); axC.set_ylim(6e-7, 3e-1)
axC.set_xlabel("code-subspace dimension  $d$", fontsize=7)
axC.set_ylabel(r"guaranteed leak  $\sigma_d=\beta_d(L)$", fontsize=7)
axC.set_title("Blindness has capacity", fontsize=7.8, color=C_CAP, fontweight="bold", pad=3)
axC.tick_params(labelsize=6)
axC.text(0.5, -0.30, r"$\beta_d(L)=\sigma_d(L)$ — and does not automatically compose",
         ha="center", va="top", fontsize=6.0, color="0.25", transform=axC.transAxes)

# ================= RIGHT : certifiability ladder =================
p = np.logspace(0, np.log10(560), 200)
KNOWN = 2.6 * np.ones_like(p)         # Ledger I  : T*lambda ~ 1     (50%-power ~2.6)
BLIND = 2.85 * np.sqrt(p)             # Ledger II : T*lambda ~ sqrt p (50%-power ~2.85)
EST = 1.0 * p                         # Ledger III: T*lambda ~ p     (rel. risk O(1) = CRB)
axR.set_xscale("log"); axR.set_yscale("log")
axR.plot(p, EST, color=C_EST, lw=1.5, label=rf"estimate  $\sim p$  ({EXPS[2]:.3f})")
axR.plot(p, BLIND, color=C_BLIND, lw=1.5, label=rf"blind  $\sim\sqrt{{p}}$  ({EXPS[1]:.3f})")
axR.plot(p, KNOWN, color=C_KNOWN, lw=1.5, label=rf"known  $\sim1$  ({EXPS[0]:.3f})")
axR.fill_between(p, KNOWN, BLIND, color=C_BLIND, alpha=0.06)
axR.fill_between(p, BLIND, EST, color=C_EST, alpha=0.06)
# KT6 operating point (below the blind line)
axR.plot([KT6_P], [KT6_TLAM], marker="X", color=C_KT6, ms=8, zorder=8,
         markeredgecolor="white", markeredgewidth=0.6)
axR.annotate(f"KT6 ranging\n"rf"$p={KT6_P}$, $T\lambda={KT6_TLAM:.2g}$"
             "\nbelow blind line",
             (KT6_P, KT6_TLAM), xytext=(KT6_P * 0.35, KT6_TLAM * 0.16),
             fontsize=5.7, color=C_KT6, ha="center", fontweight="bold",
             arrowprops=dict(arrowstyle="-|>", color=C_KT6, lw=0.9))
# projected jet (restored m=2) — certifiable, known-direction channel
axR.plot([1.0], [KNOWN[0] * 1.9], marker="*", color=C_JET, ms=11, zorder=8,
         markeredgecolor="white", markeredgewidth=0.5)
axR.annotate("projected jet\n(restored $m{=}2$)\nabove channel line", (1.0, KNOWN[0] * 1.9),
             xytext=(1.4, 40), fontsize=5.7, color=C_JET, ha="left", fontweight="bold",
             arrowprops=dict(arrowstyle="-|>", color=C_JET, lw=0.9))
axR.set_xlim(0.9, 560); axR.set_ylim(2e-2, 8e2)
axR.set_xlabel("efficient tangent dimension  $p$", fontsize=7)
axR.set_ylabel(r"required  $T\lambda$  (50% power)", fontsize=7)
axR.set_title("Three ledgers of certifiability", fontsize=7.8, color=C_EST,
              fontweight="bold", pad=3)
axR.tick_params(labelsize=6)
axR.legend(loc="upper left", fontsize=5.4, handlelength=1.2, borderpad=0.3,
           labelspacing=0.25, frameon=True, framealpha=0.9)
axR.text(0.5, -0.30, r"$1:\sqrt{p}:p$ — a nonzero signal is not a certifiable signal",
         ha="center", va="top", fontsize=6.0, color="0.25", transform=axR.transAxes)

save(fig, os.path.join(HERE, "fig1_hero"))
