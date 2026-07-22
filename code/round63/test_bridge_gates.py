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
    every record (e.g. A_realized)."""
    qbase = base + ridge_delta
    recs = []
    for grp, scs in SCENES.items():
        g = gain_ctrl if grp == "control" else gain_stress
        for sc in scs:
            for s in SEEDS:
                recs.append(_rec(sc, grp, arm, qbase + g, s, **extra))
    return recs


def _oracle_records(gain_stress, gain_ctrl, knob_stress=0.0, knob_ctrl=0.0,
                    base=20.0, ridge_delta=0.5, best_bank="L3"):
    """ORACLE-LIB records with per_bank (R28 §1).  best_bank reaches
    Q_base+gain (so the five-seed-mean argmax picks it), L0 reaches Q_base+knob;
    all other banks sit below both.  gain >= knob so the oracle gain == gain and
    the ORACLE-knob advantage == gain-knob."""
    qbase = base + ridge_delta
    recs = []
    for grp, scs in SCENES.items():
        g = gain_ctrl if grp == "control" else gain_stress
        kg = knob_ctrl if grp == "control" else knob_stress
        for sc in scs:
            for s in SEEDS:
                pb = []
                for k in range(8):
                    name = "L%d" % k
                    v = (qbase + g if name == best_bank
                         else qbase + kg if name == "L0"
                         else qbase + min(g, kg) - 1.0)
                    pb.append({"bank": name, "PSNR_rad": v})
                best = max(p["PSNR_rad"] for p in pb)
                recs.append(_rec(sc, grp, "ORACLE-LIB", best, s, per_bank=pb))
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


# ============================= LIBRARY_REACHABILITY_PASS (R28 §1, former Gate A) #
def test_library_reachability_pass_and_fail():
    # PASS: ORACLE-LIB gain = 1.5 dB on all stress scenes (best bank L3)
    recs = _base_records() + _oracle_records(1.5, 0.0)
    a = G.gate_library_reachability(recs)
    assert a["gate"] == "LIBRARY_REACHABILITY_PASS"
    assert a["PASS"] and a["c1_median>=1.0"] and a["c2_>=9of12_positive"] and a["c3_bootlb>0.50"]
    assert abs(a["median_gain_dB"] - 1.5) < 1e-9 and a["n_positive"] == 12
    assert all(v == "L3" for v in a["kstar_per_scene"].values())    # one bank/scene
    assert "class ceiling" not in a["reference_terminology"].lower().replace("not a class ceiling", "")
    # FAIL median: 0.6 dB
    a2 = G.gate_library_reachability(_base_records() + _oracle_records(0.6, 0.0))
    assert not a2["PASS"] and not a2["c1_median>=1.0"]
    # FAIL count: 8/12 positive
    recs3 = _base_records()
    for i, sc in enumerate(STRESS):
        grp = ("contour" if sc in SCENES["contour"] else "twopop" if sc in SCENES["twopop"] else "microtex")
        g = 2.0 if i >= 4 else -1.0
        for s in SEEDS:
            pb = [{"bank": "L%d" % k, "PSNR_rad": (20.5 + g if k == 3 else 20.5 + g - 1.0)}
                  for k in range(8)]
            recs3.append(_rec(sc, grp, "ORACLE-LIB", 20.5 + g, s, per_bank=pb))
    a3 = G.gate_library_reachability(recs3)
    assert a3["n_positive"] == 8 and not a3["c2_>=9of12_positive"] and not a3["PASS"]
    # min-index tie-break: L2 and L5 tie -> L2 chosen
    tie = _base_records()
    for sc in STRESS[:1]:
        for s in SEEDS:
            pb = [{"bank": "L%d" % k, "PSNR_rad": (22.0 if k in (2, 5) else 20.0)} for k in range(8)]
            tie.append(_rec(sc, "contour", "ORACLE-LIB", 22.0, s, per_bank=pb))
    assert G.gate_library_reachability(tie)["kstar_per_scene"][STRESS[0]] == "L2"
    print("[A] LIBRARY_REACHABILITY_PASS (ORACLE-LIB, five-seed argmax, min-index) PASS")


# =============================================================== Gate B ======= #
def test_gate_B_pass_and_noharm_fail():
    # PASS: ORACLE-LIB 2.0, RLMI 1.5 (>= 60% of 2.0 = 1.2), controls +0.0
    recs = (_base_records() + _oracle_records(2.0, 0.0)
            + _arm_records("RLMI", 1.5, 0.0, A_realized=0.9, w0_realized=0.1,
                           route_latency_s=0.5))
    b = G.gate_B(recs)
    assert b["PASS"] and b["c4_>=60pct_oracle"] and b["c5_noharm_controls"]
    assert abs(b["median_oracle_dB"] - 2.0) < 1e-9        # denominator is ORACLE-LIB
    # FAIL c4: RLMI 1.05 < 0.6*2.0 = 1.2
    b2 = G.gate_B(_base_records() + _oracle_records(2.0, 0.0)
                  + _arm_records("RLMI", 1.05, 0.0))
    assert not b2["c4_>=60pct_oracle"] and not b2["PASS"]
    # FAIL c5 no-harm: controls lose 0.5 dB
    b3 = G.gate_B(_base_records() + _oracle_records(2.0, 0.0)
                  + _arm_records("RLMI", 1.5, -0.5))
    assert not b3["c5_noharm_controls"] and not b3["PASS"]
    print("[B] Gate B pass + 60%%-of-ORACLE-LIB denominator + no-harm fail PASS")


# =============================================================== Gate C ======= #
def _orec(sc, grp, seed, best, knob, best_bank="L3"):
    pb = [{"bank": "L%d" % k, "PSNR_rad": (best if k == int(best_bank[1])
                                           else knob if k == 0 else min(best, knob) - 1.0)}
          for k in range(8)]
    return _rec(sc, grp, "ORACLE-LIB", best, seed, per_bank=pb)


def test_gate_C_pass_and_fail():
    # rescue scenes: ORACLE-LIB advantage 2 dB (best-knob), RLMI A_j 0.9; controls
    # advantage 0.1 (<0.25), w0 0.9.  Advantage read from per_bank (max - L0).
    recs = _base_records()
    resc = STRESS[:6]
    for sc in STRESS:
        grp = ("contour" if sc in SCENES["contour"] else "twopop" if sc in SCENES["twopop"] else "microtex")
        adv = 2.0 if sc in resc else 0.1
        A = 0.9 if sc in resc else 0.2
        for s in SEEDS:
            recs.append(_orec(sc, grp, s, 22.0, 22.0 - adv))
            recs.append(_rec(sc, grp, "RLMI", 21.0, s, A_realized=A,
                             w0_realized=1.0 - A, route_latency_s=0.5))
    for sc in CTRL:
        for s in SEEDS:
            recs.append(_orec(sc, "control", s, 25.0, 24.9))
            recs.append(_rec(sc, "control", "RLMI", 25.0, s, A_realized=0.1,
                             w0_realized=0.9, route_latency_s=0.5))
    cc = G.gate_C(recs)
    assert cc["c1_meanA>=0.80"] and cc["c2_meanw0>=0.75"] and cc["PASS"]
    assert cc["n_rescue_needed"] == 6 and cc["n_aligned_lowadv"] == 4
    recs2 = [dict(r) for r in recs]
    for r in recs2:
        if r["arm"] == "RLMI" and r["scene"] in resc:
            r["A_realized"] = 0.3; r["w0_realized"] = 0.7
    assert not G.gate_C(recs2)["c1_meanA>=0.80"] and not G.gate_C(recs2)["PASS"]
    print("[C] Gate C (per_bank advantage; A_j rescue / w0 controls) pass + C1 fail PASS")


# ================================================= PLUGIN_LATENCY_PASS (R27) == #
def test_plugin_latency():
    fast = [_rec("s", "contour", "RLMI", 20, s, route_latency_s=0.4) for s in range(20)]
    d = G.gate_plugin_latency(fast, cpu="ref-cpu", threads=1, cache_state="cold")
    assert d["gate"] == "PLUGIN_LATENCY_PASS" and d["PASS"]
    assert d["gates_M2"] is False                      # R27: never kills M2
    assert d["claim_label"].startswith("subsecond")
    assert d["cpu"] == "ref-cpu" and d["threads"] == 1 and d["cache_state"] == "cold"
    slow = [_rec("s", "contour", "RLMI", 20, s, route_latency_s=9.6) for s in range(20)]
    d2 = G.gate_plugin_latency(slow)
    assert not d2["PASS"] and d2["gates_M2"] is False
    assert "batch" in d2["claim_label"] and d2["max_s"] == 9.6
    # median ok but p95 spike -> plugin FAIL, still no M2 kill
    mix = ([_rec("s", "contour", "RLMI", 20, s, route_latency_s=0.4) for s in range(18)]
           + [_rec("s", "contour", "RLMI", 20, s, route_latency_s=5.0) for s in range(2)])
    d3 = G.gate_plugin_latency(mix)
    assert d3["c1_median<1.0s"] and not d3["c2_p95<2.0s"] and not d3["PASS"]
    assert G.gate_D(fast)["PASS"]                       # back-compat alias
    print("[D] PLUGIN_LATENCY_PASS claim-classification (never gates M2) PASS")


# ================================================= four-quantity decomp (R28 §4) #
def _full_cohort(oracle_gain, rlmi_gain, truex, xhat, oracle_knob=0.0):
    """Full A-C cohort: ORACLE-LIB (per_bank, gain=oracle_gain over stress, best-
    L0=oracle_gain-oracle_knob), RLMI (rlmi_gain), and dose-relaxed FW arms."""
    recs = _base_records()
    for grp, scs in SCENES.items():
        st = grp != "control"
        og = oracle_gain if st else 0.1
        rg = rlmi_gain if st else 0.0
        for sc in scs:
            for s in SEEDS:
                recs.append(_orec(sc, grp, s, 20.5 + og, 20.5 + og - (oracle_knob if st else 0.05)))
                recs.append(_rec(sc, grp, "RLMI", 20.5 + rg, s,
                                 A_realized=(0.9 if st else 0.1),
                                 w0_realized=(0.1 if st else 0.9), route_latency_s=0.5))
                recs.append(_rec(sc, grp, "TRUE-X-FW", 20.5 + truex, s,
                                 fw_label="NONDEPLOYABLE_DOSE_RELAXED_FW_DIAGNOSTIC"))
                recs.append(_rec(sc, grp, "XHAT-FW", 20.5 + xhat, s,
                                 fw_label="NONDEPLOYABLE_DOSE_RELAXED_FW_DIAGNOSTIC"))
    return recs


def test_four_quantity_decomposition():
    # oracle 1.8, rlmi 1.5, truex 2.0, xhat 1.7
    recs = _full_cohort(1.8, 1.5, 2.0, 1.7, oracle_knob=1.0)
    fg = G.four_gap(recs)
    assert abs(fg["q1_inclass_library_reachability_dB"] - 1.8) < 1e-9   # ORACLE-base
    assert abs(fg["q2_deployable_allocation_loss_dB"] - (1.5 - 1.8)) < 1e-9   # RLMI-ORACLE
    assert abs(fg["q3_plugin_fw_loss_dB_descriptive"] - (1.7 - 2.0)) < 1e-9   # XHAT-TRUEX
    assert abs(fg["q4_dose_relax_library_contrast_dB_descriptive"] - (2.0 - 1.8)) < 1e-9  # TRUEX-ORACLE
    assert "not gated" in fg["note"].lower()
    print("[gap] R28 four-quantity decomposition PASS")


# ===================================== R28 §3: no Gate A-C statistic reads an FW ==
def test_no_gate_reads_fw():
    recs = _full_cohort(1.8, 1.5, 2.0, 1.7, oracle_knob=1.0)
    assert G.gates_read_no_fw(recs) is True
    # sanity: the invariant would catch a gate that read an FW arm (flip via a
    # deliberately FW-reading probe)
    a = G.gate_library_reachability(recs); b = G.gate_B(recs); c = G.gate_C(recs)
    perturbed = [dict(r) for r in recs]
    for r in perturbed:
        if r["arm"] in ("TRUE-X-FW", "XHAT-FW"):
            r["PSNR_rad"] = 999.0
    assert G.gate_library_reachability(perturbed)["median_gain_dB"] == a["median_gain_dB"]
    assert G.gate_B(perturbed)["median_oracle_dB"] == b["median_oracle_dB"]
    assert G.gate_C(perturbed)["PASS"] == c["PASS"]
    print("[R28] no Gate A-C statistic reads an FW arm PASS")


def test_leave_smoke_out_sensitivity():
    recs = _full_cohort(1.8, 1.5, 2.0, 1.7, oracle_knob=1.0)
    s = G.leave_smoke_out_sensitivity(recs)
    assert s["disclosure_only"] and set(s["excluded_smoke_exposed"]) == {
        "bridge_twopop_0", "bridge_control_0"}
    # official gate uses the full 12-stress cohort; leave-out drops the one
    # smoke-exposed STRESS scene (twopop_0; control_0 is a control) -> 11 stress
    assert s["gate_A_leaveout"]["n_scenes"] == 11
    print("[R28] leave-smoke-exposed-out sensitivity (no-gate) PASS")


# ================================================= PLUGIN_LATENCY_PASS (R27) == #
def test_plugin_latency():
    fast = [_rec("s", "contour", "RLMI", 20, s, route_latency_s=0.4) for s in range(20)]
    d = G.gate_plugin_latency(fast, cpu="ref-cpu", threads=1, cache_state="cold")
    assert d["gate"] == "PLUGIN_LATENCY_PASS" and d["PASS"]
    assert d["gates_M2"] is False and d["claim_label"].startswith("subsecond")
    assert d["cpu"] == "ref-cpu" and d["threads"] == 1 and d["cache_state"] == "cold"
    slow = [_rec("s", "contour", "RLMI", 20, s, route_latency_s=9.6) for s in range(20)]
    d2 = G.gate_plugin_latency(slow)
    assert not d2["PASS"] and d2["gates_M2"] is False
    assert "batch" in d2["claim_label"] and d2["max_s"] == 9.6
    mix = ([_rec("s", "contour", "RLMI", 20, s, route_latency_s=0.4) for s in range(18)]
           + [_rec("s", "contour", "RLMI", 20, s, route_latency_s=5.0) for s in range(2)])
    d3 = G.gate_plugin_latency(mix)
    assert d3["c1_median<1.0s"] and not d3["c2_p95<2.0s"] and not d3["PASS"]
    assert G.gate_D(fast)["PASS"]
    print("[D] PLUGIN_LATENCY_PASS claim-classification (never gates M2) PASS")


def test_decision_tree_and_latency_verdict():
    # LIBRARY_REACHABILITY fail -> M2 stop
    fail = _full_cohort(0.2, 0.2, 0.2, 0.2)
    rep = G.evaluate_all(fail)
    assert not rep["gate_library_reachability"]["PASS"]
    assert "LIBRARY_REACHABILITY_FAIL" in rep["m2_verdict"]
    assert rep["no_gate_reads_fw"] is True
    # A-C pass but SLOW latency -> M2 still launches, claim -> batch (oracle_knob
    # 1.5 > 1.0 so stress scenes are rescue-needed and Gate C's C1 applies)
    recs = _full_cohort(2.0, 1.5, 2.0, 1.7, oracle_knob=1.5)
    for r in recs:
        if r["arm"] == "RLMI":
            r["route_latency_s"] = 9.6
    rep2 = G.evaluate_all(recs)
    assert (rep2["gate_library_reachability"]["PASS"] and rep2["gate_B"]["PASS"]
            and rep2["gate_C"]["PASS"])
    assert not rep2["plugin_latency"]["PASS"]
    assert "M2 preregistration" in rep2["m2_verdict"]
    assert "batch" in rep2["claim_label"].lower()
    print("[tree] R28 decision tree (reachability-fail stop; latency-fail keeps M2) PASS")


ALL = [test_bootstrap_lb, test_library_reachability_pass_and_fail,
       test_gate_B_pass_and_noharm_fail, test_gate_C_pass_and_fail,
       test_four_quantity_decomposition, test_no_gate_reads_fw,
       test_leave_smoke_out_sensitivity, test_plugin_latency,
       test_decision_tree_and_latency_verdict]

if __name__ == "__main__":
    for fn in ALL:
        fn()
    print("\nALL BRIDGE-GATE TESTS PASS (%d)" % len(ALL))
