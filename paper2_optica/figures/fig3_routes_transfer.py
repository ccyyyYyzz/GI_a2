"""
PAPER 2 (Optica) — FIGURE 3: three restoration routes + calibration transfer
(R45 §C3 fig3 spec).

Panel L — the descending leak ladder of the three physically distinct restoration
routes (plus the exact statistical projection), all reducing the forbidden mean-channel
response:
    raw random codes   5.6e-2   (IT4  meanleak_random_codes)
      -> DPSS window    2.4e-4   (IT5  dpss_max_witness_leak, hardware concentration)
      -> certified codes <=1e-5  (IT4b single-plane 80 dims, code-space capacity)
      -> statistical projection -> 0  (KT3 clean nuisance overlap after = 0, exact)
  The routes are ALTERNATIVES, not factors (noncomposable, §3 / hero center).

Panel R — calibration transfer (TLSG arm B, Eq. 4.4): a fresh calibration/test split
validates sigma_d(Lhat)+delta as a simultaneous cannot-see envelope: envelope 1.07e-2,
worst fresh bottom-64 leak 2.0e-6 (six orders below), 100% coverage on 24 fresh draws,
and >=80 code dims still certified <1e-4 at the conjugate plane.

DATA: IT4 @33b8ef2, IT5 @c8a0487, IT4b @5666d1e, KT3 @2430166, TLSG @35d858e (committed,
extracted to _frozen_sources/). No data generation.
"""
import json
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from o2_style import (apply_style, save, COL2, OI, C_BREAK, C_DPSS, C_CAP, C_KNOWN,
                      C_EST, C_GREY)

HERE = os.path.dirname(os.path.abspath(__file__))
FS = os.path.join(HERE, "_frozen_sources")
it4 = json.load(open(os.path.join(FS, "IT4_svd_certified_blind_codes.committed.json"))) \
    if os.path.exists(os.path.join(FS, "IT4_svd_certified_blind_codes.committed.json")) else None
it5 = json.load(open(os.path.join(FS, "IT5_dpss_concentration.committed.json")))
it4b = json.load(open(os.path.join(FS, "IT4b_spectrum_trade.committed.json")))
kt3 = json.load(open(os.path.join(FS, "KT3_leak_orthogonalized_sentinel.committed.json")))
tlsg = json.load(open(os.path.join(FS, "TLSG_RESULTS.committed.json")))

# route levels (committed)
RAW = 5.6365251833210714e-02          # IT4 meanleak_random_codes (extracted below if present)
if it4 is not None:
    RAW = it4["measured"]["meanleak_random_codes"]
DPSS = it5["verdict"]["dpss_max_witness_leak"]                       # 2.38e-4
CERT = 1e-5                                                          # IT4b 80 dims <=1e-5 (threshold cleared)
CERT_FLOOR = it4b["results"]["joint_spectrum"]["rel_leak_vs_d"]["1"]  # 1.43e-6 best-code floor

apply_style()
fig, (axL, axR) = plt.subplots(1, 2, figsize=(COL2, 3.0))

# ---- Panel L: descending route ladder ----
routes = [("raw random\ncodes", RAW, C_BREAK, "IT4"),
          ("DPSS window\n(hardware)", DPSS, C_DPSS, "IT5"),
          ("certified codes\n(80 dims, code-space)", CERT, C_CAP, "IT4b"),
          ("statistical\nprojection", None, C_KNOWN, "KT3")]
xs = np.arange(len(routes))
axL.set_yscale("log")
for x, (lab, val, col, src) in zip(xs, routes):
    if val is not None:
        axL.bar(x, val, width=0.62, color=col, alpha=0.85, zorder=3)
        axL.text(x, val * 1.5, f"{val:.1e}".replace("e-0", r"$\times10^{-").replace("e-", r"$\times10^{-") + "}$",
                 ha="center", va="bottom", fontsize=6.0, fontweight="bold", color=col)
    else:
        # exact projection -> 0 : draw an arrow diving to the floor
        axL.annotate("", xy=(x, 3e-6), xytext=(x, 3e-3),
                     arrowprops=dict(arrowstyle="-|>", color=col, lw=2.2))
        axL.text(x, 5e-3, r"overlap $\to 0$", ha="center", va="bottom", fontsize=6.2,
                 fontweight="bold", color=col)
        axL.text(x, 2.2e-6, "(exact)", ha="center", va="top", fontsize=5.6, color=col)
    axL.text(x, 1.3e-6, src, ha="center", va="top", fontsize=5.3, color="0.45", style="italic")
# connecting descent arrows
for x in xs[:-1]:
    axL.annotate("", xy=(x + 0.72, routes[x + 1][1] if routes[x + 1][1] else 2e-4),
                 xytext=(x + 0.28, routes[x][1]),
                 arrowprops=dict(arrowstyle="-|>", color="0.55", lw=1.0,
                                 connectionstyle="arc3,rad=-0.15"))
axL.axhline(1e-5, color="0.6", lw=0.7, ls=(0, (3, 2)))
axL.text(3.35, 1.15e-5, r"$10^{-5}$", fontsize=5.6, color="0.5", va="bottom")
axL.set_xticks(xs)
axL.set_xticklabels([r[0] for r in routes], fontsize=5.9)
axL.set_ylim(1e-6, 3e-1)
axL.set_ylabel("forbidden mean-channel leak")
axL.set_title("Three routes reduce the same leak — noncomposably", fontsize=7.0, pad=3)
kt3_ret = kt3["verdict"]["retained_info_frac"]
axL.text(0.03, 0.045,
         f"KT3 target info retained {min(kt3_ret)*100:.0f}–{max(kt3_ret)*100:.0f}%",
         transform=axL.transAxes, fontsize=5.7, color=C_KNOWN, va="bottom")

# ---- Panel R: calibration-transfer coverage (TLSG arm B) ----
bar6 = next(b for b in tlsg["bars"] if b["bar"] == 6)["measured"]
import re
sig64 = float(re.search(r"sigma_64=([0-9.eE+-]+)", bar6).group(1))    # 3.41e-8
delta = float(re.search(r"delta=([0-9.eE+-]+)", bar6).group(1))       # 1.07e-2
envelope = float(re.search(r"envelope=([0-9.eE+-]+)", bar6).group(1)) # 1.07e-2
worst = float(re.search(r"worst S64 leak on fresh=([0-9.eE+-]+)", bar6).group(1))  # 1.99e-6
bar7 = next(b for b in tlsg["bars"] if b["bar"] == 7)["measured"]
fresh_d = int(re.search(r"min = (\d+)", bar7).group(1))               # 80

def sci(x, p=1):
    m, e = f"{x:.{p}e}".split("e")
    return rf"{m}\times10^{{{int(e)}}}"

axR.set_yscale("log")
# envelope band (certified-blind region below envelope)
axR.axhspan(1e-9, envelope, color=C_CAP, alpha=0.09, zorder=0)
axR.axhline(envelope, color=C_EST, lw=1.6, zorder=4)
axR.text(0.5, envelope * 1.45,
         rf"cannot-see envelope  $\sigma_{{64}}(\hat L)+\delta={sci(envelope,2)}$",
         ha="center", va="bottom", fontsize=6.0, color=C_EST, fontweight="bold")
# calibration sigma_64 and worst fresh leak
axR.plot([0.28], [sig64], marker="o", color=C_KNOWN, ms=6, zorder=5)
axR.text(0.28, sig64 * 2.2, rf"calib $\sigma_{{64}}={sci(sig64,1)}$",
         ha="center", va="bottom", fontsize=5.8, color=C_KNOWN)
axR.plot([0.72], [worst], marker="D", color=C_DPSS, ms=6, zorder=5)
axR.text(0.72, worst * 2.4, "worst fresh\nbottom-64 leak\n"rf"${sci(worst,1)}$",
         ha="center", va="bottom", fontsize=5.8, color=C_DPSS, fontweight="bold")
axR.set_xlim(0, 1); axR.set_ylim(1e-9, 5e-1)
axR.set_xticks([])
axR.set_ylabel("bottom-64 subspace leak")
axR.set_title("Calibration transfer: a fresh-data cannot-see certificate", fontsize=6.8, pad=3)
axR.text(0.5, 0.045,
         f"coverage = 1.000 (24/24 fresh draws)\n$\\geq${fresh_d} code dims certified $<10^{{-4}}$",
         transform=axR.transAxes, ha="center", va="bottom", fontsize=6.2,
         color=C_CAP, fontweight="bold",
         bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=C_CAP, lw=0.8, alpha=0.9))

save(fig, os.path.join(HERE, "fig3_routes_transfer"))
