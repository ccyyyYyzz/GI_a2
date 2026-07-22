#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Consolidate Part 1 (M1_cont, M0_const) + Part 2 into ladder_results.json and
print the summary tables used in the report.  Read-only on the run JSONs."""
import json, os, numpy as np
OUT = "D:/GI_another/results/round63_next/LADDER_PROBE"
def L(f): return json.load(open(os.path.join(OUT, f)))
m1 = L("ladder_results_M1_cont_full.json")
m0 = L("ladder_results_M0_const_full.json")
p2 = L("part2_jitter_results.json")
SCENES = m1["meta"]["scenes"]
RUNGS = m1["meta"]["rungs"]
CELLS = m1["meta"]["cells"]

def cell_table(res, cellname):
    ps = res["results"][cellname]
    rows = {}
    for sc in SCENES:
        d = ps[sc]
        rows[sc] = {k: d[k] for k in ["A2", "A3", "B1", "B4", "B16", "B64", "oracle"]}
    def med(k): return float(np.median([rows[sc][k] for sc in SCENES]))
    med_row = {k: med(k) for k in ["A2", "A3", "B1", "B4", "B16", "B64", "oracle"]}
    best_sw = [max(rows[sc]["A2"], rows[sc]["A3"]) for sc in SCENES]
    ts_vs_sw = [rows[sc]["B64"] - max(rows[sc]["A2"], rows[sc]["A3"]) for sc in SCENES]
    ts_vs_b1 = [rows[sc]["B64"] - rows[sc]["B1"] for sc in SCENES]
    ofrac = []
    for sc in SCENES:
        base = max(rows[sc]["A2"], rows[sc]["A3"]); gap = rows[sc]["oracle"] - base
        ofrac.append((rows[sc]["B64"] - base) / gap if gap > 1e-6 else 0.0)
    # saturation: rung increments on the median ladder
    lad = [med_row[f"B{B}"] for B in RUNGS]
    incs = [lad[i + 1] - lad[i] for i in range(len(lad) - 1)]
    return dict(per_scene=rows, median=med_row,
                headline_ts_vs_bestCPLsoftware_dB=float(np.median(ts_vs_sw)),
                refine_only_ts_vs_B1smoother_dB=float(np.median(ts_vs_b1)),
                oracle_fraction=float(np.median(ofrac)),
                ladder_increments_dB=[float(x) for x in incs])

consol = {"meta": {"models": {"M1_cont": "continuous within-exposure drift (honest variant)",
                              "M0_const": "frozen within-exposure-constant gain (CPL model)"},
                   "rungs": RUNGS, "cells": CELLS, "scenes": SCENES,
                   "B_MAX_timestamp_proxy": m1["meta"]["B_MAX"],
                   "PHI": m1["meta"]["PHI"], "N_EXP": m1["meta"]["N_EXP"]},
           "part1_M1_cont": {c: cell_table(m1, c) for c in CELLS},
           "part1_M0_const": {c: cell_table(m0, c) for c in CELLS},
           "part2_jitter": {k: {kk: v[kk] for kk in
                    ["rho", "J_count", "J_c0", "jitter_loss_dB",
                     "recovered_fraction_analytic", "recovered_dB_analytic",
                     "EG_v"]} | {"recovered_fraction_mc": v["mc"]["recovered_fraction_mc"]}
                    for k, v in p2["results"].items()}}
json.dump(consol, open(os.path.join(OUT, "ladder_results.json"), "w"), indent=2)

# ---- print summary tables ----
def pt(model_key, title):
    print(f"\n===== {title} =====")
    print(f"{'cell':20s} {'A2':>6} {'A3':>6} | {'B1':>6} {'B4':>6} {'B16':>6} {'B64':>6} | "
          f"{'oracle':>6} | {'TSvsSW':>7} {'refine':>7} {'ofrac':>6}")
    for c in CELLS:
        t = consol[model_key][c]; md = t["median"]
        print(f"{c:20s} {md['A2']:6.2f} {md['A3']:6.2f} | {md['B1']:6.2f} {md['B4']:6.2f} "
              f"{md['B16']:6.2f} {md['B64']:6.2f} | {md['oracle']:6.2f} | "
              f"{t['headline_ts_vs_bestCPLsoftware_dB']:+7.2f} "
              f"{t['refine_only_ts_vs_B1smoother_dB']:+7.2f} {t['oracle_fraction']:6.2f}")
pt("part1_M1_cont", "PART 1  M1_cont (continuous drift, honest PRIMARY variant)")
pt("part1_M0_const", "PART 1  M0_const (frozen within-constant gain -> refinement VACUOUS)")
print("\n===== PART 2  jitter/dead-time timestamp recovery (nu=2000, c=0.05) =====")
for k, v in consol["part2_jitter"].items():
    print(f"  rho={v['rho']:6.2f}  J_count={v['J_count']:.4f} J_c0={v['J_c0']:.4f}  "
          f"jitter_loss={v['jitter_loss_dB']:.3f}dB  recovered_frac={v['recovered_fraction_analytic']:.3f}"
          f"(mc {v['recovered_fraction_mc']:.3f})  recovered={v['recovered_dB_analytic']:.3f}dB")
print("\nsaved ladder_results.json")
