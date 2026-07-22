#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""DLGI probe -- figure generator (3 panels): T1 medium precision vs oracle floor
and pilot; the no-trade schedule result; and the canonical-confusion reciprocity
identity.  Writes figs/fig_dual_ledger.png."""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "D:/GI_another/results/round63_next/DUAL_LEDGER_PROBE"
G = json.load(open(os.path.join(OUT, "dual_ledger_results.json")))["grid"]
T3 = json.load(open(os.path.join(OUT, "t3_reciprocity.json")))

cells = ["tc2_cv5", "tc2_cv15", "tc2_cv40", "tc16_cv5", "tc16_cv15",
         "tc16_cv40", "tc64_cv5", "tc64_cv15", "tc64_cv40"]
labels = ["2/5", "2/15", "2/40", "16/5", "16/15", "16/40", "64/5", "64/15", "64/40"]

fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))

# ---- Panel A: T1 t_c precision, pilot-free vs pilot vs oracle floor ----
xa = np.arange(len(cells))
pf = np.array([abs(G[c]["pf_tc_relerr_med"]) * 100 for c in cells])
pil = np.array([min(abs(G[c]["pilot_tc_relerr_med"]) * 100, 200) for c in cells])
flo = np.array([G[c]["floor_tc_rel"] * 100 for c in cells])
w = 0.38
ax[0].bar(xa - w/2, pf, w, label="pilot-free (DLGI)", color="#2c7fb8")
ax[0].bar(xa + w/2, pil, w, label="pilot 5% (clipped 200%)", color="#d95f0e", alpha=0.8)
ax[0].plot(xa, flo, "k_", ms=14, mew=2.5, label="oracle CRB floor (R33 eq4.6)")
ax[0].axhline(20, ls="--", color="gray", lw=1)
ax[0].text(0.1, 21.5, "20% bar", color="gray", fontsize=8)
ax[0].set_xticks(xa); ax[0].set_xticklabels(labels, rotation=45, fontsize=8)
ax[0].set_xlabel("cell  (t_c / CV%)"); ax[0].set_ylabel("|t_c relative error|  (%)")
ax[0].set_title("A. T1 medium channel: correlation-time precision\n"
                "pilot-free tracks the oracle floor; pilot fails at fast drift",
                fontsize=9)
ax[0].set_ylim(0, 60); ax[0].legend(fontsize=7, loc="upper right")

# ---- Panel B: no-trade schedule result (paired best for BOTH ledgers) ----
sc = T3["empirical_schedule_scan"]
scheds = ["paired", "random", "split"]
xb = np.arange(3)
for i, (tc, mk) in enumerate([(16, "o"), (64, "s")]):
    psnr = [sc[f"tc{tc}_cv15_{s}"]["scene_psnr_med"] for s in scheds]
    terr = [sc[f"tc{tc}_cv15_{s}"]["medium_tc_absrelerr_med"] * 100 for s in scheds]
    ax[1].plot(xb, psnr, mk + "-", color="#2c7fb8", label=f"scene PSNR (t_c={tc})")
    ax[1].plot(xb, terr, mk + "--", color="#d95f0e", label=f"medium t_c err% (t_c={tc})")
ax[1].set_xticks(xb); ax[1].set_xticklabels(scheds)
ax[1].set_xlabel("schedule (same photons/patterns; order only)")
ax[1].set_ylabel("scene PSNR (dB)  /  medium t_c err (%)")
ax[1].set_title("B. No forced trade: PAIRED is best for BOTH ledgers\n"
                "(scene sensitivity = differential data-processing, R33 sec.3)",
                fontsize=9)
ax[1].legend(fontsize=7)

# ---- Panel C: canonical-confusion reciprocity  eff_scene == eff_medium ----
abk = T3["marginal_fisher_ABK"]
es = [abk[k]["eff_scene"] for k in abk]
em = [abk[k]["eff_medium"] for k in abk]
ax[2].plot([0.87, 0.91], [0.87, 0.91], "k-", lw=1, label="y = x (reciprocity, Thm 2)")
ax[2].scatter(es, em, s=70, color="#7b3294", zorder=3, label="schedules x cells")
ax[2].set_xlabel("scene ledger efficiency  det(J_x)/det(A)")
ax[2].set_ylabel("medium ledger efficiency  det(J_theta)/det(B)")
ax[2].set_title("C. Canonical-Confusion Ledger Reciprocity\n"
                "same normalized loss taxes both ledgers (exact to 1e-15)",
                fontsize=9)
ax[2].legend(fontsize=8)
ax[2].set_aspect("equal"); ax[2].grid(alpha=0.3)

plt.tight_layout()
fn = os.path.join(OUT, "figs", "fig_dual_ledger.png")
os.makedirs(os.path.dirname(fn), exist_ok=True)
plt.savefig(fn, dpi=130, bbox_inches="tight")
print("wrote", fn)
