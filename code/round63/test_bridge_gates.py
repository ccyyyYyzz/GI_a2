"""Unit tests for bridge_gates.py — DEV-bridge kill gates on SYNTHETIC outcomes.

The gates are decision rules (R24 §3.4 / R25 §12 / R24 §3.5); these tests build
synthetic per-cell records engineered to PASS and to FAIL each gate, and assert
the verdicts + the four-gap decomposition + the decision tree.  No real bridge
data is touched (phase-1 rule).  Runnable directly or under pytest.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bridge_gates as G

SCENES = {
    "contour": ["bridge_contour_%d" % i for i in range(4)],
    "twopop": ["bridge_twopop_%d" % i for i in range(4)],
    "microtex": ["bridge_microtex_%d" % i for i in range(4)],
    "control": ["bridge_control_%d" % i for i in range(4)],
}
STRESS = SCENES["contour"] + SCENES["twopop"] + SCENES["microtex"]
CTRL = SCENES["control"]
NU, C, SEEDS = 2000.0, 0.05, range(5)


def _rec(scene, group, arm, psnr, seed=0, **extra):
    return {"scene": scene, "group": group, "arm": arm, "nu": NU, "c": C,
            "seed": seed, "PSNR_rad": float(psnr), **extra}


def _base_records(base=20.0, ridge_delta=0.5, noise=0.0, rng=None):
    """SCAT32-060 and RIDGE-SCAT32 records for all 16 scenes x 5 seeds.
    Q_base,j = max = base + ridge_delta (RIDGE the stronger comparator)."""
    recs = []
    for grp, scs in SCENES.items():
        for sc in scs:
            for s in SEEDS:
                n = (rng.normal(0, noise) if (rng is not None and noise) else 0.0)
                recs.append(_rec(sc, grp, "SCAT32-060", base + n, s))
                recs.append(_rec(sc, grp, "RIDGE-SCAT32", base + ridge_delta + n, s))
    return recs


def _arm_records(arm, gain_stress, gain_ctrl, base=20.0, ridge_delta=0.5,
                 **extra):
    """arm PSNR = Q_base + gain (constant per group); `extra` (flat) applied to
    every record (e.g. A_realized, oracle_best_psnr)."""
    qbase = base + ridge_delta
    recs = []
    for grp, scs in SCENES.items():
        g = gain_ctrl if grp == "control" else gain_stress
        for sc in scs:
            for s in SEEDS:
                recs.append(_rec(sc, grp, arm, qbase + g, s, **extra))
    return recs


# =============================================================== bootstrap === #
def test_bootstrap_lb():
    v = np.full(12, 1.0)                              # zero-variance -> LB == 1.0
    lb = G.scene_bootstrap_median_lb(v, conf=0.90, B=1000, seed=1)
    assert abs(lb - 1.0) < 1e-9
    v2 = np.concatenate([np.full(6, 2.0), np.full(6, 0.0)])
    lb2 = G.scene_bootstrap_median_lb(v2, conf=0.90, B=4000, seed=1)
    assert 0.0 <= lb2 <= 2.0                          # median straddles 0..2
    print("[boot] scene-bootstrap median LB PASS (lb=%.3f, lb2=%.3f)" % (lb, lb2))


# =============================================================== Gate A ======= #
def test_gate_A_pass_and_fail():
    # PASS: TRUE-X gain = 1.5 dB on all stress scenes
    recs = _base_records() + _arm_records("TRUE-X-FW", 1.5, 0.0)
    a = G.gate_A(recs)
    assert a["PASS"] and a["c1_median>=1.0"] and a["c2_>=9of12_positive"] and a["c3_bootlb>0.50"]
    assert abs(a["median_gain_dB"] - 1.5) < 1e-9 and a["n_positive"] == 12
    # FAIL median: gain 0.6 dB (median < 1.0, and LB < 0.5)
    recs2 = _base_records() + _arm_records("TRUE-X-FW", 0.6, 0.0)
    a2 = G.gate_A(recs2)
    assert not a2["PASS"] and not a2["c1_median>=1.0"]
    # FAIL count: 8/12 positive (4 scenes at -1 dB), median still ~ high
    per = {sc: {} for sc in STRESS}
    recs3 = _base_records()
    for i, sc in enumerate(STRESS):
        grp = ("contour" if sc in SCENES["contour"] else "twopop" if sc in SCENES["twopop"] else "microtex")
        g = 2.0 if i >= 4 else -1.0                  # 8 positive, 4 negative
        for s in SEEDS:
            recs3.append(_rec(sc, grp, "TRUE-X-FW", 20.5 + g, s))
    a3 = G.gate_A(recs3)
    assert a3["n_positive"] == 8 and not a3["c2_>=9of12_positive"] and not a3["PASS"]
    print("[A] Gate A pass + median/count fail modes PASS")


# =============================================================== Gate B ======= #
def test_gate_B_pass_and_noharm_fail():
    # PASS: TRUE-X 2.0, RLMI 1.5 (>= 60% of 2.0 = 1.2), controls +0.0 (no harm)
    recs = (_base_records() + _arm_records("TRUE-X-FW", 2.0, 0.0)
            + _arm_records("RLMI", 1.5, 0.0, A_realized=0.9, w0_realized=0.1,
                           route_latency_s=0.5))
    b = G.gate_B(recs)
    assert b["PASS"] and b["c4_>=60pct_truex"] and b["c5_noharm_controls"]
    # FAIL c4: RLMI gain 1.0 < 0.6*2.0=1.2
    recs2 = (_base_records() + _arm_records("TRUE-X-FW", 2.0, 0.0)
             + _arm_records("RLMI", 1.05, 0.0))
    b2 = G.gate_B(recs2)
    assert not b2["c4_>=60pct_truex"] and not b2["PASS"]
    # FAIL c5 no-harm: controls lose 0.5 dB (median < -0.25)
    recs3 = (_base_records() + _arm_records("TRUE-X-FW", 2.0, 0.0)
             + _arm_records("RLMI", 1.5, -0.5))
    b3 = G.gate_B(recs3)
    assert not b3["c5_noharm_controls"] and not b3["PASS"]
    print("[B] Gate B pass + 60%%-of-TRUEX + no-harm fail modes PASS")


# =============================================================== Gate C ======= #
def test_gate_C_pass_and_fail():
    # rescue-needed scenes: give 6 stress scenes ORACLE-LIB advantage 2 dB with
    # RLMI A_j = 0.9 (>=0.80).  Controls: advantage 0.1 (<0.25) with w0 = 0.9.
    recs = _base_records()
    resc = STRESS[:6]
    for sc in STRESS:
        grp = ("contour" if sc in SCENES["contour"] else "twopop" if sc in SCENES["twopop"] else "microtex")
        adv = 2.0 if sc in resc else 0.1
        A = 0.9 if sc in resc else 0.2
        for s in SEEDS:
            recs.append(_rec(sc, grp, "ORACLE-LIB", 22.0, s,
                             oracle_best_psnr=22.0, oracle_knob_psnr=22.0 - adv))
            recs.append(_rec(sc, grp, "RLMI", 21.0, s, A_realized=A,
                             w0_realized=1.0 - A, route_latency_s=0.5))
    for sc in CTRL:
        for s in SEEDS:
            recs.append(_rec(sc, "control", "ORACLE-LIB", 25.0, s,
                             oracle_best_psnr=25.0, oracle_knob_psnr=24.9))
            recs.append(_rec(sc, "control", "RLMI", 25.0, s, A_realized=0.1,
                             w0_realized=0.9, route_latency_s=0.5))
    cc = G.gate_C(recs)
    assert cc["c1_meanA>=0.80"] and cc["c2_meanw0>=0.75"] and cc["PASS"]
    assert cc["n_rescue_needed"] == 6 and cc["n_aligned_lowadv"] == 4
    # FAIL C1: rescue scenes but RLMI stays on the knob (A_j = 0.3)
    recs2 = []
    for r in recs:
        r = dict(r)
        if r["arm"] == "RLMI" and r["scene"] in resc:
            r["A_realized"] = 0.3; r["w0_realized"] = 0.7
        recs2.append(r)
    cc2 = G.gate_C(recs2)
    assert not cc2["c1_meanA>=0.80"] and not cc2["PASS"]
    print("[C] Gate C (A_j rescue / w0 controls) pass + C1 fail PASS")


# =============================================================== Gate D ======= #
def test_gate_D():
    fast = [_rec("s", "contour", "RLMI", 20, s, route_latency_s=0.4) for s in range(20)]
    d = G.gate_D(fast)
    assert d["PASS"] and d["c1_median<1.0s"] and d["c2_p95<2.0s"]
    slow = [_rec("s", "contour", "RLMI", 20, s, route_latency_s=3.0) for s in range(20)]
    d2 = G.gate_D(slow)
    assert not d2["PASS"] and not d2["c1_median<1.0s"]
    # median ok but p95 spike
    mix = ([_rec("s", "contour", "RLMI", 20, s, route_latency_s=0.4) for s in range(18)]
           + [_rec("s", "contour", "RLMI", 20, s, route_latency_s=5.0) for s in range(2)])
    d3 = G.gate_D(mix)
    assert d3["c1_median<1.0s"] and not d3["c2_p95<2.0s"] and not d3["PASS"]
    print("[D] Gate D latency (median + p95) PASS")


# =============================================================== four-gap ==== #
def test_four_gap_and_verdict():
    recs = (_base_records()
            + _arm_records("TRUE-X-FW", 2.0, 0.0)
            + _arm_records("XHAT-FW", 1.7, 0.0)          # gap2 = 1.7-2.0 = -0.3
            + _arm_records("ORACLE-LIB", 1.8, 0.0, oracle_best_psnr=22.3,
                           oracle_knob_psnr=20.3)         # gap3 = 1.8-2.0 = -0.2
            + _arm_records("RLMI", 1.5, 0.0, A_realized=0.9, w0_realized=0.1,
                           route_latency_s=0.5))          # gap4 = 1.5-1.8 = -0.3
    fg = G.four_gap(recs)
    assert abs(fg["gap1_truex_vs_fixed_dB"] - 2.0) < 1e-9
    assert abs(fg["gap2_xhat_minus_truex_dB"] + 0.3) < 1e-9
    assert abs(fg["gap3_oraclelib_minus_truex_dB"] + 0.2) < 1e-9
    assert abs(fg["gap4_rlmi_minus_oraclelib_dB"] + 0.3) < 1e-9
    rep = G.evaluate_all(recs)
    assert "PASS" in rep["verdict"] or "M2" in rep["verdict"]
    print("[gap] four-gap decomposition + verdict PASS: %s" % rep["verdict"])


def test_decision_tree_gateA_fail():
    # Gate A fails -> method + M2 dead, regardless of the others
    recs = (_base_records() + _arm_records("TRUE-X-FW", 0.2, 0.0)
            + _arm_records("RLMI", 0.2, 0.0, A_realized=0.9, w0_realized=0.1,
                           route_latency_s=0.5)
            + _arm_records("XHAT-FW", 0.2, 0.0)
            + _arm_records("ORACLE-LIB", 0.2, 0.0, oracle_best_psnr=20,
                           oracle_knob_psnr=20))
    rep = G.evaluate_all(recs)
    assert not rep["gate_A"]["PASS"]
    assert "DEAD" in rep["verdict"]
    print("[tree] Gate-A-fail decision tree PASS: %s" % rep["verdict"])


ALL = [test_bootstrap_lb, test_gate_A_pass_and_fail, test_gate_B_pass_and_noharm_fail,
       test_gate_C_pass_and_fail, test_gate_D, test_four_gap_and_verdict,
       test_decision_tree_gateA_fail]

if __name__ == "__main__":
    for fn in ALL:
        fn()
    print("\nALL BRIDGE-GATE TESTS PASS (%d)" % len(ALL))
