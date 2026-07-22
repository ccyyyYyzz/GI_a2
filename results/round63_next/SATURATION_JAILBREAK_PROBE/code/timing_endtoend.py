"""
R32 BAR 4 -- END-TO-END WALL-CLOCK OVERHEAD (not MHz gate arithmetic).

R32 forbids "negligible time based only on MHz gate arithmetic".  The equal-time
frame is legitimate ONLY if the whole MOLT sweep fits the measured system
schedule with <= 10% end-to-end overhead INCLUDING mask changes, settling, power
switches and readout.  We model every hardware cost:

  DMD mask switch : tau_switch = 1/R_dmd,  R_dmd in {1, 10, 22} kHz
                    (10 kHz DMD -> 100 us/pattern switch, per the coordinator note)
  settling        : tau_settle after each switch
  microcell gate  : tau_gate = 1/R_gate,   R_gate in {100 kHz, 1 MHz}
  power switch    : tau_pow, 3-level sweep = 2 extra power switches/sparse mask
  readout         : tau_read per gate-block (per level)

Comparison (equal DMD-slot budget N = K_D + K_S so switching is common-mode):
  LINEAR arm : every slot is one linear dwell (G_lin gates).
  MOLT arm   : K_D linear dwells + K_S three-level sweeps (ridge-heavy gates +
               2 power switches + 2 extra readouts each).
  overhead % = (T_MOLT - T_LINEAR) / T_LINEAR  -- reported for all
  {R_gate} x {R_dmd}, and for two sparse-bank sizes K_S in {K_D, K_D/2}.

Gate counts from the PRODUCTION dose design: linear dwell at t_lin=0.1; sparse
sweep at the three-level design for a 10% p2 map (bar-2 budget).
"""
import json
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import saturation_core as sc
import gi_operator as op

OUT = os.path.dirname(HERE)
FIGS = os.path.join(OUT, "figs")
DATA = "data/r63_bridge_scenes"

C = 3600
B_LIN = 1.0e4
K_D = 51                                   # dense linear bank (original GI system)

# hardware timing constants (microseconds)
TAU_SETTLE = 10.0                          # DMD settling after a switch
TAU_POW = 10.0                             # power-level switch (per coordinator)
TAU_READ = 5.0                             # readout per gate-block/level
R_DMD = {"1kHz": 1.0e3, "10kHz": 1.0e4, "22kHz": 2.2e4}      # -> us switch = 1e6/R
R_GATE = {"100kHz": 1.0e5, "1MHz": 1.0e6}                    # -> us/gate = 1e6/R
T_LIN = 0.10                               # linear operating occupancy

# ---- gate counts from the production design -------------------------------- #
G_lin = B_LIN / (C * T_LIN)                # linear dwell gates (~27.8)
# 10% p2 map sweep on a representative sparse mask (n_eff~19 -> bar-2 ~x11)
x = op.load_scene(os.path.join(DATA, "bridge_contour_2.npz"))
mrng = np.random.default_rng(651_605)
m = np.zeros(op.P); m[mrng.choice(op.P, size=32, replace=False)] = 1.0
v = (m * x); v = v[v > 0]
p1 = float(v.sum()); n_eff = p1 ** 2 / float((v ** 2).sum())
B_map = 11.3 * B_LIN                        # bar-2 production multiple (max over masks)
a3, G3 = sc.three_level_design(p1, budget_incident=B_map, C=C)
G_sweep = float(np.sum(G3))                 # total gates over the 3 levels
n_levels = len(G3)
print("gate counts: G_lin=%.1f  G_sweep(3-level, 10%% p2 map)=%.1f (levels=%d, per-level %s)"
      % (G_lin, G_sweep, n_levels, G3.tolist()))


def slot_time_linear(tau_switch, tau_gate):
    return tau_switch + TAU_SETTLE + G_lin * tau_gate + TAU_READ


def slot_time_sweep(tau_switch, tau_gate):
    # one DMD switch, then 3 power levels each = gates + readout, + 2 power switches
    return (tau_switch + TAU_SETTLE
            + sum(g * tau_gate + TAU_READ for g in G3)
            + (n_levels - 1) * TAU_POW)


def campaign_overhead(K_S, R_g, R_d):
    tau_gate = 1.0e6 / R_g
    tau_switch = 1.0e6 / R_d
    N = K_D + K_S
    T_linear = N * slot_time_linear(tau_switch, tau_gate)
    T_molt = (K_D * slot_time_linear(tau_switch, tau_gate)
              + K_S * slot_time_sweep(tau_switch, tau_gate))
    return T_linear, T_molt, 100.0 * (T_molt - T_linear) / T_linear


report = dict(K_D=K_D, G_lin=float(G_lin), G_sweep=float(G_sweep),
    levels_gates=G3.tolist(), tau_settle=TAU_SETTLE, tau_pow=TAU_POW,
    tau_read=TAU_READ, t_lin=T_LIN, B_map_multiple=11.3,
    R_dmd=R_DMD, R_gate=R_GATE, overhead_pct={})

print("\nEND-TO-END OVERHEAD %% (equal DMD-slot budget; MOLT sweep vs linear dwell):")
grid = {}
for K_S in (K_D, K_D // 2):
    print("  --- K_S = %d sparse masks (K_D=%d dense) ---" % (K_S, K_D))
    print("  %-8s" % "gate\\dmd" + "".join("%10s" % d for d in R_DMD))
    tbl = {}
    for gname, R_g in R_GATE.items():
        row = {}
        line = "  %-8s" % gname
        for dname, R_d in R_DMD.items():
            Tl, Tm, ov = campaign_overhead(K_S, R_g, R_d)
            row[dname] = dict(overhead_pct=float(ov), T_linear_us=float(Tl),
                              T_molt_us=float(Tm))
            line += "%9.1f%%" % ov
        tbl[gname] = row
        print(line)
    grid[f"K_S_{K_S}"] = tbl
report["overhead_pct"] = grid

# pass/fail: <=10% requires the switching-limited regime (slow DMD or fast gates)
passing = {}
for ks, tbl in grid.items():
    for gname, row in tbl.items():
        for dname, cell in row.items():
            passing[f"{ks}|{gname}|{dname}"] = bool(cell["overhead_pct"] <= 10.0)
n_pass = sum(passing.values())
report["bar4_summary"] = dict(
    n_regimes=len(passing), n_pass_10pct=n_pass,
    passing_regimes=[k for k, ok in passing.items() if ok],
    verdict="MOLT meets <=10% end-to-end overhead ONLY in switching-limited "
            "regimes (fast gates and/or slow DMD); at 100 kHz gates with a 22 kHz "
            "DMD the ridge-heavy sweep exceeds 10%. The equal-time framing is "
            "hardware-contingent, NOT universal. 'Negligible time' is false in the "
            "gate-limited corner and must be stated per-regime.",
    BAR="<=10% end-to-end wall-clock overhead including mask/settle/power/readout")
print("\nBAR 4: %d/%d (R_gate x R_dmd x K_S) regimes meet <=10%% overhead" % (n_pass, len(passing)))

# ---- figure: overhead heatmap-ish grouped bars (K_S=K_D) ----
fig, ax = plt.subplots(figsize=(6.6, 4.2))
dnames = list(R_DMD); gnames = list(R_GATE)
width = 0.38
xpos = np.arange(len(dnames))
for gi, gname in enumerate(gnames):
    ov = [grid[f"K_S_{K_D}"][gname][d]["overhead_pct"] for d in dnames]
    ax.bar(xpos + (gi - 0.5) * width, ov, width, label="gates @ " + gname)
ax.axhline(10, ls="--", color="k", alpha=0.6, label="10% budget")
ax.set_xticks(xpos); ax.set_xticklabels(["DMD @ " + d for d in dnames])
ax.set_ylabel("end-to-end overhead  %"); ax.set_yscale("log")
ax.set_title("MOLT sweep wall-clock overhead vs linear (K_S=K_D=%d)" % K_D)
ax.legend(fontsize=8); ax.grid(True, axis="y", which="both", alpha=0.25)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "fig_timing_overhead.png"), dpi=130)
print("wrote figs/fig_timing_overhead.png")

with open(os.path.join(OUT, "timing_endtoend.json"), "w") as fh:
    json.dump(report, fh, indent=2, default=float)
print("wrote timing_endtoend.json")
