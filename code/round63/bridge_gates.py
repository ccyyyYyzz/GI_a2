"""ROUND63 DEV bridge — kill-gate evaluation (T-D deliverable d).

Decision rules ONLY (R24 §3.4 verbatim; Gate C replaced by R25 §12; the four-gap
decomposition R24 §3.5).  No constant is tunable: every threshold below is quoted
from the ruling.  These functions consume a flat list of per-cell result records
and return verdicts; they are unit-tested on SYNTHETIC outcomes (test_bridge_
gates.py) and are NOT run on real bridge data in phase 1.

Record schema (one per scene x arm x nu x c x seed)
  {"scene": str, "group": "contour|twopop|microtex|control",
   "arm": "SCAT32-060|RIDGE-SCAT32|TRUE-X-FW|XHAT-FW|RLMI|ORACLE-LIB",
   "nu": float, "c": float, "seed": int, "PSNR_rad": float,
   # RLMI only:
   "A_realized": float,            # R25.4  A_j = 1 - w_hat_{j,0}
   "w0_realized": float,           # w_hat_{j,0}
   "route_latency_s": float,       # Gate D (rlmi alloc_time_s, ex-recon)
   # ORACLE-LIB only:
   "oracle_best_psnr": float,      # best bank PSNR_rad (the arm's reported value)
   "oracle_knob_psnr": float}      # the L0 (knob) bank sub-acquisition PSNR_rad

The hard bridge decision cell is (c=0.05, nu=2000) (R24 §3.2); other cells are
diagnostic.  Per-scene values aggregate the 5 paired seeds by MEAN (paired
radiometric PSNR); Q_base,j = max{Q_SCAT32-060,j, Q_RIDGE-SCAT32,j} (R24 §3.3).
"""
from __future__ import annotations

import numpy as np

STRESS_GROUPS = ("contour", "twopop", "microtex")   # 12 misalignment-stress scenes
CONTROL_GROUP = "control"                            # 4 aligned controls
DECISION_NU = 2000.0
DECISION_C = 0.05
BOOT_CONF = 0.90                                     # one-sided lower-bound level
BOOT_B = 2000


# --------------------------------------------------------------------------- #
#  aggregation
# --------------------------------------------------------------------------- #
def _sel(records, nu, c):
    return [r for r in records if r["nu"] == nu and r["c"] == c]


def scene_means(records, arm, nu, c, field="PSNR_rad"):
    """{scene: mean over the paired seeds of `field` for `arm`} at (nu,c)."""
    acc = {}
    for r in _sel(records, nu, c):
        if r["arm"] != arm:
            continue
        acc.setdefault(r["scene"], []).append(float(r[field]))
    return {sc: float(np.mean(v)) for sc, v in acc.items()}


def scene_groups_map(records):
    return {r["scene"]: r["group"] for r in records}


def q_base(records, nu, c):
    """Q_base,j = max{Q_SCAT32-060,j, Q_RIDGE-SCAT32,j} per scene (R24 §3.3)."""
    a = scene_means(records, "SCAT32-060", nu, c)
    b = scene_means(records, "RIDGE-SCAT32", nu, c)
    return {sc: max(a[sc], b[sc]) for sc in a if sc in b}


def paired_gains(records, arm, nu, c, groups=STRESS_GROUPS):
    """Per-scene gain Q_arm,j - Q_base,j restricted to `groups`.  Returns
    {scene: gain} (scene means; the pairing is the shared Q_base per scene)."""
    gmap = scene_groups_map(records)
    qa = scene_means(records, arm, nu, c)
    qb = q_base(records, nu, c)
    return {sc: qa[sc] - qb[sc] for sc in qa
            if sc in qb and gmap.get(sc) in groups}


# --------------------------------------------------------------------------- #
#  scene bootstrap (lower bound on the median)
# --------------------------------------------------------------------------- #
def scene_bootstrap_median_lb(values, conf=BOOT_CONF, B=BOOT_B, seed=650000):
    """One-sided `conf` lower bound on the median by resampling the SCENES with
    replacement (R24 §3.4 "90% scene-bootstrap lower bound on the median").
    LB = the (1-conf) quantile of the bootstrap median distribution."""
    v = np.asarray(values, dtype=np.float64)
    if v.size == 0:
        return float("nan")
    rng = np.random.default_rng(seed)
    meds = np.empty(B)
    for b in range(B):
        idx = rng.integers(0, v.size, size=v.size)
        meds[b] = np.median(v[idx])
    return float(np.quantile(meds, 1.0 - conf))


# --------------------------------------------------------------------------- #
#  Gate A — Fisher-guided geometry improves images (R24 §3.4)
# --------------------------------------------------------------------------- #
def gate_A(records, nu=DECISION_NU, c=DECISION_C, boot_seed=650000):
    g = paired_gains(records, "TRUE-X-FW", nu, c, STRESS_GROUPS)
    scenes = sorted(g)
    vals = np.array([g[s] for s in scenes])
    med = float(np.median(vals)) if vals.size else float("nan")
    n_pos = int((vals > 0).sum())
    lb = scene_bootstrap_median_lb(vals, seed=boot_seed)
    c1 = med >= 1.0
    c2 = n_pos >= 9
    c3 = lb > 0.50
    return {"gate": "A", "nu": nu, "c": c, "n_scenes": len(scenes),
            "median_gain_dB": med, "n_positive": n_pos, "boot_median_lb_dB": lb,
            "c1_median>=1.0": bool(c1), "c2_>=9of12_positive": bool(c2),
            "c3_bootlb>0.50": bool(c3),
            "PASS": bool(c1 and c2 and c3), "per_scene": g}


# --------------------------------------------------------------------------- #
#  Gate B — the deployable router retains the gain (R24 §3.4)
# --------------------------------------------------------------------------- #
def gate_B(records, nu=DECISION_NU, c=DECISION_C, boot_seed=650001):
    gr = paired_gains(records, "RLMI", nu, c, STRESS_GROUPS)
    scenes = sorted(gr)
    vals = np.array([gr[s] for s in scenes])
    med = float(np.median(vals)) if vals.size else float("nan")
    n_pos = int((vals > 0).sum())
    lb = scene_bootstrap_median_lb(vals, seed=boot_seed)
    gt = paired_gains(records, "TRUE-X-FW", nu, c, STRESS_GROUPS)
    med_truex = float(np.median([gt[s] for s in sorted(gt)])) if gt else float("nan")
    # controls
    gc = paired_gains(records, "RLMI", nu, c, (CONTROL_GROUP,))
    cvals = np.array([gc[s] for s in sorted(gc)])
    med_ctrl = float(np.median(cvals)) if cvals.size else float("nan")
    n_below1 = int((cvals < -1.0).sum())
    c1 = med >= 1.0
    c2 = n_pos >= 9
    c3 = lb > 0.0
    c4 = med >= 0.60 * med_truex if np.isfinite(med_truex) else False
    c5 = (med_ctrl >= -0.25) and (n_below1 <= 1)
    return {"gate": "B", "nu": nu, "c": c, "n_stress": len(scenes),
            "median_gain_dB": med, "n_positive": n_pos, "boot_median_lb_dB": lb,
            "median_truex_dB": med_truex, "60pct_truex_dB": 0.60 * med_truex,
            "control_median_dB": med_ctrl, "n_controls_below_-1dB": n_below1,
            "c1_median>=1.0": bool(c1), "c2_>=9of12_positive": bool(c2),
            "c3_bootlb>0": bool(c3), "c4_>=60pct_truex": bool(c4),
            "c5_noharm_controls": bool(c5),
            "PASS": bool(c1 and c2 and c3 and c4 and c5),
            "per_scene_stress": gr, "per_scene_control": gc}


# --------------------------------------------------------------------------- #
#  Gate C — allocation informativeness (R25 §12; A_j from materialized weights)
# --------------------------------------------------------------------------- #
def _oracle_knob_advantage(records, nu, c):
    """Per-scene ORACLE-LIB - knob(L0) advantage (DEV-only label; R25 §12).
    Uses the ORACLE-LIB record's oracle_best_psnr / oracle_knob_psnr (means over
    seeds)."""
    best = scene_means(records, "ORACLE-LIB", nu, c, field="oracle_best_psnr")
    knob = scene_means(records, "ORACLE-LIB", nu, c, field="oracle_knob_psnr")
    return {sc: best[sc] - knob[sc] for sc in best if sc in knob}


def gate_C(records, nu=DECISION_NU, c=DECISION_C):
    gmap = scene_groups_map(records)
    adv = _oracle_knob_advantage(records, nu, c)
    A = scene_means(records, "RLMI", nu, c, field="A_realized")
    w0 = scene_means(records, "RLMI", nu, c, field="w0_realized")
    # C1: rescue-needed scenes (ORACLE-LIB beats knob by > 1 dB)
    resc = [sc for sc in adv if adv[sc] > 1.0]
    A_resc = [A[sc] for sc in resc if sc in A]
    meanA = float(np.mean(A_resc)) if A_resc else float("nan")
    c1 = (meanA >= 0.80) if A_resc else False
    # C2: aligned controls where ORACLE-LIB advantage < 0.25 dB
    ctrl = [sc for sc in adv if gmap.get(sc) == CONTROL_GROUP and adv[sc] < 0.25]
    w0_ctrl = [w0[sc] for sc in ctrl if sc in w0]
    meanW0 = float(np.mean(w0_ctrl)) if w0_ctrl else float("nan")
    c2 = (meanW0 >= 0.75) if w0_ctrl else False
    return {"gate": "C", "nu": nu, "c": c,
            "n_rescue_needed": len(resc), "mean_A_j": meanA,
            "n_aligned_lowadv": len(ctrl), "mean_w0_j": meanW0,
            "c1_meanA>=0.80": bool(c1), "c2_meanw0>=0.75": bool(c2),
            "PASS": bool(c1 and c2),
            "rescue_scenes": resc, "aligned_control_scenes": ctrl}


# --------------------------------------------------------------------------- #
#  Gate D — route latency (R24 §3.4)
# --------------------------------------------------------------------------- #
def gate_D(records, field="route_latency_s"):
    lat = [float(r[field]) for r in records
           if r["arm"] == "RLMI" and r.get(field) is not None]
    if not lat:
        return {"gate": "D", "PASS": False, "reason": "no RLMI latencies",
                "median_s": float("nan"), "p95_s": float("nan")}
    med = float(np.median(lat))
    p95 = float(np.percentile(lat, 95))
    c1 = med < 1.0
    c2 = p95 < 2.0
    return {"gate": "D", "n": len(lat), "median_s": med, "p95_s": p95,
            "c1_median<1.0s": bool(c1), "c2_p95<2.0s": bool(c2),
            "PASS": bool(c1 and c2)}


# --------------------------------------------------------------------------- #
#  four-gap decomposition (R24 §3.5) — reported, not gated
# --------------------------------------------------------------------------- #
def _median_paired_diff(records, arm_hi, arm_lo, nu, c, groups=STRESS_GROUPS):
    a = scene_means(records, arm_hi, nu, c)
    b = scene_means(records, arm_lo, nu, c)
    gmap = scene_groups_map(records)
    diffs = [a[sc] - b[sc] for sc in a
             if sc in b and gmap.get(sc) in groups]
    return float(np.median(diffs)) if diffs else float("nan")


def four_gap(records, nu=DECISION_NU, c=DECISION_C):
    """Four separate gaps (medians over the 12 stress scenes; R24 §3.5):
      1 TRUE-X FW vs fixed baseline (surrogate-to-image),
      2 XHAT FW vs TRUE-X FW      (pre-scan / localization loss),
      3 ORACLE-LIB vs TRUE-X FW   (finite-bank approximation loss),
      4 RLMI vs ORACLE-LIB        (routing loss).
    Never averaged into one 'adaptive gain'."""
    gap1 = float(np.median([v for v in paired_gains(
        records, "TRUE-X-FW", nu, c, STRESS_GROUPS).values()]))
    gap2 = _median_paired_diff(records, "XHAT-FW", "TRUE-X-FW", nu, c)
    gap3 = _median_paired_diff(records, "ORACLE-LIB", "TRUE-X-FW", nu, c)
    gap4 = _median_paired_diff(records, "RLMI", "ORACLE-LIB", nu, c)
    return {"nu": nu, "c": c,
            "gap1_truex_vs_fixed_dB": gap1,
            "gap2_xhat_minus_truex_dB": gap2,
            "gap3_oraclelib_minus_truex_dB": gap3,
            "gap4_rlmi_minus_oraclelib_dB": gap4}


def evaluate_all(records, nu=DECISION_NU, c=DECISION_C):
    """Full gate report at the hard decision cell (R24 decision tree)."""
    A = gate_A(records, nu, c)
    B = gate_B(records, nu, c)
    C = gate_C(records, nu, c)
    D = gate_D(records)
    gaps = four_gap(records, nu, c)
    # decision tree (R24 §4 / build plan): Gate A fail -> method + M2 dead.
    if not A["PASS"]:
        verdict = "GATE_A_FAIL -> method line + M2 DEAD (paper 2 = jitter/knob)"
    elif A["PASS"] and B["PASS"] and D["PASS"] and C["PASS"]:
        verdict = "A-D PASS -> M2 preregistration; RLMI leads paper 2"
    elif A["PASS"] and B["PASS"] and D["PASS"] and not C["PASS"]:
        verdict = ("A,B,D PASS but C FAIL -> may be a static composite design; "
                   "MUST NOT lead M2 as adaptive (R25 §12)")
    else:
        verdict = "A pass but B/D fail -> do not launch M2 (no rescue redesign)"
    return {"decision_cell": {"nu": nu, "c": c},
            "gate_A": A, "gate_B": B, "gate_C": C, "gate_D": D,
            "four_gap": gaps, "verdict": verdict}
