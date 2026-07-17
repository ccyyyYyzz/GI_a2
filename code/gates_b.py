"""B1/B2 gate evaluation + paired hierarchical bootstrap — spec §8.

Reads phaseB_core_metrics.csv / phaseB_ext_metrics.csv, writes phaseB_gates.json.

Gate-eligible combos (nonlinear links only): DT30 (core, both photon levels)
and SAT30/SAT50/GAMMA07/LOG (extended, s=1e4). LIN and FGAIN are honesty arms:
reported, never gated. MLE-OR is an info upper bound: never in a control set.

Per-combo pass (each of RELATION / PRACTICAL):
  mean dPSNR >= +0.50 dB, mean dSSIM >= +0.010, mean dLPIPS_benefit >= +0.010,
  all three simultaneously improved on >= 6/8 images (per-image seed-means),
  and hierarchical-bootstrap LB90 of the PSNR gain > 0 (primary; SSIM/LPIPS
  LB90 reported). Baseline = per-(combo, image, metric) best of the gate's
  control set, frozen on seed-means for bootstrap pairing.

Gate PASS: passing combos span >=2 illuminations x >=2 nonlinear links and
include DT30 (spec: >=2 non-Gaussian illum x >=2 nonlinear links, must contain
DT30; all Phase-B illuminations are non-Gaussian).
"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gi_core import config as C
from gi_core.utils import rng_for

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")

REL_SET = ["WHITEN-OR", "SIR-10", "SIR-20", "L-ISOTRON"]
PRAC_SET = ["WHITEN-LW", "DGI", "CORR", "SIR-10", "SIR-20", "L-ISOTRON"]
METRICS = ["PSNR", "SSIM", "LPIPS"]


def load_tables():
    core = pd.read_csv(os.path.join(RESULTS, "phaseB_core_metrics.csv"))
    ext_path = os.path.join(RESULTS, "phaseB_ext_metrics.csv")
    ext = pd.read_csv(ext_path) if os.path.exists(ext_path) else pd.DataFrame(columns=core.columns)
    return core, ext


def eligible_combos(core, ext):
    combos = []
    for (il, lk, ph), _ in core[core.link == "DT30"].groupby(["illum", "link", "photons"]):
        combos.append((il, lk, float(ph), "core"))
    for (il, lk, ph), _ in ext.groupby(["illum", "link", "photons"]):
        combos.append((il, lk, float(ph), "ext"))
    return combos


def combo_eval(df, illum, link, photons, base_set, rng):
    sub = df[(df.illum == illum) & (df.link == link) & (df.photons == photons)]
    if sub.empty or "SCORE-OR" not in set(sub.method):
        return None
    out = {"images": {}, "per_seed_mean_delta_psnr": {}}
    lpips_ok = np.isfinite(sub["LPIPS"]).all()
    # per-image seed-means
    piv = {m: sub.pivot_table(index="method", columns="image", values=m,
                              aggfunc="mean") for m in METRICS}
    base_choice = {}
    deltas = {m: {} for m in METRICS}
    for img in C.IMAGES:
        for m in METRICS:
            col = piv[m][img]
            base_vals = col.loc[[b for b in base_set if b in col.index]]
            if m == "LPIPS":
                bmeth = base_vals.idxmin()
                deltas[m][img] = float(base_vals.min() - col.loc["SCORE-OR"])  # benefit
            else:
                bmeth = base_vals.idxmax()
                deltas[m][img] = float(col.loc["SCORE-OR"] - base_vals.max())
            base_choice[(m, img)] = bmeth
    mean_d = {m: float(np.mean(list(deltas[m].values()))) for m in METRICS}
    n_all3 = sum(1 for img in C.IMAGES
                 if all(deltas[m][img] > 0 for m in METRICS))
    # paired per-(image, seed) gains with frozen per-image baseline method
    seeds = sorted(sub.seed.unique())
    gains = {m: np.full((len(C.IMAGES), len(seeds)), np.nan) for m in METRICS}
    for i, img in enumerate(C.IMAGES):
        for j, sd in enumerate(seeds):
            cell = sub[(sub.image == img) & (sub.seed == sd)].set_index("method")
            for m in METRICS:
                b = base_choice[(m, img)]
                if m == "LPIPS":
                    gains[m][i, j] = cell.loc[b, m] - cell.loc["SCORE-OR", m]
                else:
                    gains[m][i, j] = cell.loc["SCORE-OR", m] - cell.loc[b, m]
    # per-seed directions (report)
    out["per_seed_mean_delta_psnr"] = {
        str(sd): float(np.nanmean(gains["PSNR"][:, j])) for j, sd in enumerate(seeds)}
    # hierarchical bootstrap: resample images, then seeds within image
    lb90 = {}
    nI, nS = len(C.IMAGES), len(seeds)
    for m in METRICS:
        stats = np.empty(C.BOOT_N)
        img_idx = rng.integers(0, nI, size=(C.BOOT_N, nI))
        seed_idx = rng.integers(0, nS, size=(C.BOOT_N, nI, nS))
        g = gains[m]
        for t in range(C.BOOT_N):
            sel = g[img_idx[t][:, None], seed_idx[t]]
            stats[t] = sel.mean()
        lb90[m] = float(np.quantile(stats, C.BOOT_LB_Q))
    thresholds = {"PSNR": C.B1_PSNR_DB, "SSIM": C.B1_SSIM, "LPIPS": C.B1_LPIPS}
    mean_gate = all(mean_d[m] >= thresholds[m] for m in METRICS)
    combo_pass = bool(mean_gate and n_all3 >= C.B1_MIN_POS and lb90["PSNR"] > 0)
    out.update({
        "mean_delta": mean_d, "n_images_all3_improved": int(n_all3),
        "lb90": lb90, "lpips_available": bool(lpips_ok),
        "delta_per_image": {m: deltas[m] for m in METRICS},
        "baseline_choice": {"%s|%s" % k: v for k, v in base_choice.items()},
        "pass": combo_pass,
    })
    return out


def main():
    core, ext = load_tables()
    combos = eligible_combos(core, ext)
    rng = rng_for(0, 77, 1)
    gates = {"combos": {}, "honesty_arms_note":
             "LIN and FGAIN reported in phaseB_core_metrics.csv, never gated; "
             "MLE-OR reported as information upper bound, never a control."}
    for gate_name, base_set in (("RELATION_HEADROOM", REL_SET),
                                ("PRACTICAL_HEADROOM", PRAC_SET)):
        passing = []
        detail = {}
        for il, lk, ph, tab in combos:
            df = core if tab == "core" else ext
            ev = combo_eval(df, il, lk, ph, base_set, rng)
            if ev is None:
                continue
            detail["%s|%s|%g" % (il, lk, ph)] = ev
            if ev["pass"]:
                passing.append((il, lk, ph))
        ills = {p[0] for p in passing}
        lks = {p[1] for p in passing}
        has_dt30 = any(p[1] == "DT30" for p in passing)
        full_pass = bool(len(ills) >= 2 and len(lks) >= 2 and has_dt30)
        gates["combos"][gate_name] = detail
        gates[gate_name] = {
            "passing_combos": ["%s|%s|%g" % p for p in passing],
            "n_passing": len(passing), "n_evaluated": len(detail),
            "n_illum": len(ills), "n_links": len(lks), "contains_DT30": has_dt30,
            "KILL": len(passing) == 0,          # spec: no passing combo anywhere
            "FULL_PASS": full_pass,             # >=2 illum x >=2 links, incl. DT30
        }
    rel, prac = gates["RELATION_HEADROOM"], gates["PRACTICAL_HEADROOM"]
    stage0_kill = bool(rel["KILL"] or prac["KILL"])
    dt30_both = bool(rel["contains_DT30"] and prac["contains_DT30"])
    if stage0_kill:
        regime = "KILL"
    elif dt30_both:
        regime = "HONEST_NONLINEARITY"
    else:
        regime = "THEORY_ONLY"                  # only strong math arms pass
    if stage0_kill:
        verdict = "KILL"
    elif regime == "HONEST_NONLINEARITY" and rel["FULL_PASS"] and prac["FULL_PASS"]:
        verdict = "PROCEED_TO_PHASE_C"
    else:
        verdict = "THEORY_PAPER_ONLY"
    gates["B2_REGIME"] = regime
    gates["STAGE0_KILL"] = stage0_kill
    gates["VERDICT"] = verdict
    with open(os.path.join(RESULTS, "phaseB_gates.json"), "w") as f:
        json.dump(gates, f, indent=2)
    print("B1 RELATION:", "KILL" if rel["KILL"] else
          "PASS(%d/%d, DT30=%s)" % (rel["n_passing"], rel["n_evaluated"], rel["contains_DT30"]),
          "| B1 PRACTICAL:", "KILL" if prac["KILL"] else
          "PASS(%d/%d, DT30=%s)" % (prac["n_passing"], prac["n_evaluated"], prac["contains_DT30"]),
          "| B2:", regime, "| VERDICT:", verdict, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
