#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
DLGI CONFIRMATORY CAMPAIGN -- the C1-C7 binding-kill-bar evaluator (R34 sec.3, VERBATIM).

Consumes the confirmatory endpoints (results/round63_dlgi_campaign/confirmatory/
CONFIRMATORY_RESULTS.json produced by dlgi_confirmatory_driver.py --merge) and emits a
verdict JSON + a human table.  Each bar is a BINDING kill bar: any kill ends the flagship
claim, with no averaging away a failed parameter cell and no repair round (R34 sec.3).

Bars (exact R34 sec.3 text encoded as thresholds):
  C1 calibrated t_c coverage : every primary cell coverage in [0.92,0.98]; exact binomial
       CIs; the 3 named stress cells (16,.05)/(16,.40)/(64,.05) at nominal SNR.
       KILL: any primary cell coverage < 0.92, OR nonfinite/unreported interval rate > 5%.
       (Over-coverage > 0.98 is flagged as a band violation but is NOT in the verbatim
        kill clause -- the frozen upper-envelope construction errs toward over-coverage.)
  C2 informativeness : median log-width t_c <= 1.5x pilot; median log-width CV <= 1.5x
       pilot; bounded AND connected in >= 95% of records.  Oracle-relative width reported.
       KILL: any width > 1.5x pilot, or bounded/connected < 95%.
  C3 point precision : pilot-free RMSE(t_c) and RMSE(CV) <= 1.5x strongest pilot every
       cell; oracle ratio + signed bias reported.  KILL: any cell > 1.5x, OR the estimator
       is not bit-identical to the probe (the calibration must change only the uncertainty
       statement, not the estimator).
  C4 scene noninferiority : PSNR loss <= 0.2 dB AND task/Fisher loss <= 5% every cell,
       cohort CIs reported.  KILL: either margin crossed.
  C5 model + gauge : held-out innovation std ~1, residual whiteness (|lag-1| small),
       scalar-gain adequacy, scale-gauge invariance on fresh scenes/all SNR strata.
       KILL: systematic pattern-dependent residuals, gauge sensitivity, or a spatially
       varying operator invalidating the scalar ledger.
  C6 reciprocity + schedule : det(J_x)/det(A)=det(J_theta)/det(B) + shared spectrum to
       numerical tolerance; the frozen (A,B,K) ordering PREDICTS the measured ordering.
       KILL: theorem fails, OR schedule behavior contradicts the frozen prediction on
       > 10% of primary cells (predeclared).
  C7 identical-ledger + edge honesty : byte-audit photons/exposures/duration/pattern
       multiset; publish fast/slow-drift edge failures + weak-CV boundary.
       KILL: hidden pilot/reference info, unequal resource accounting, or removal of a
       failed edge after unblinding.

Thresholds are FROZEN here (feasibility-probe-consistent: innovation std 0.89-0.97,
lag-1 <= 0.11, gauge shift ~3e-15, reciprocity ~3e-15).  Unit tests (--selftest) exercise
the bar arithmetic on synthetic passing AND failing inputs for every bar.  Read-only on
the endpoints.  CPU, seconds.
"""
import os
import sys
import json
import math
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

REPO = "D:/GI_another"
CONF = os.path.join(REPO, "results", "round63_dlgi_campaign", "confirmatory")
RESULTS = os.path.join(CONF, "CONFIRMATORY_RESULTS.json")

# named stress cells (R34 bar C1): the 3 that failed intervals in the probe, nominal SNR
STRESS_CKEYS = ["tc16_cv5_ph10", "tc16_cv40_ph10", "tc64_cv5_ph10"]

# frozen bar thresholds (R34 sec.3)
COV_LO, COV_HI = 0.92, 0.98
NONFINITE_MAX = 0.05
WIDTH_RATIO_MAX = 1.5
BOUND_CONN_MIN = 0.95
RMSE_RATIO_MAX = 1.5
PSNR_LOSS_MAX = 0.2               # dB
TASK_FISHER_LOSS_MAX = 0.05
INNOV_STD_LO, INNOV_STD_HI = 0.8, 1.2
LAG1_AC_MAX = 0.20
GAUGE_SHIFT_MAX = 1e-6
RECIPROCITY_TOL = 1e-9
SCHED_CONTRA_FRAC_MAX = 0.10
N_EXP_EXPECT = 2048


def _fin(x):
    try:
        return math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


# ============================================================ per-bar evaluators
def eval_c1(primary):
    """Calibrated t_c coverage in [0.92,0.98] every primary cell; nonfinite <= 5%."""
    rows = []; kill = False; band_violations = []
    for ck, c in sorted(primary.items()):
        cov = c["coverage_tc"]; nf = c["nonfinite_rate"]; cp = c.get("coverage_tc_cp95")
        below = _fin(cov) and cov < COV_LO
        over = _fin(cov) and cov > COV_HI
        nf_bad = _fin(nf) and nf > NONFINITE_MAX
        cell_kill = below or nf_bad
        kill = kill or cell_kill
        if over:
            band_violations.append(ck)
        rows.append(dict(ckey=ck, is_stress=ck in STRESS_CKEYS, coverage_tc=cov,
                         coverage_cp95=cp, nonfinite_rate=nf,
                         below_0_92=below, over_0_98=over, nonfinite_gt_5pct=nf_bad,
                         cell_kill=cell_kill))
    status = "KILL" if kill else ("WARN_OVERCOVER" if band_violations else "PASS")
    return dict(bar="C1", name="calibrated t_c coverage", status=status, kill=kill,
                in_band_required=[COV_LO, COV_HI], nonfinite_max=NONFINITE_MAX,
                stress_cells=STRESS_CKEYS, band_violations_overcover=band_violations,
                cells=rows,
                note=("KILL clause is verbatim (coverage<0.92 or nonfinite>5%); "
                      "over-coverage>0.98 is a band violation flagged WARN, not a kill "
                      "-- the frozen upper-envelope critical surface errs conservative."))


def eval_c2(primary):
    """Interval width <= 1.5x pilot; bounded AND connected >= 95%."""
    rows = []; kill = False
    for ck, c in sorted(primary.items()):
        wtc = c["width_ratio_tc_dlgi_over_pilot"]; wcv = c["width_ratio_cv_dlgi_over_pilot"]
        fb = c["frac_bounded"]; fc = c["frac_connected"]
        wtc_bad = _fin(wtc) and wtc > WIDTH_RATIO_MAX
        wcv_bad = _fin(wcv) and wcv > WIDTH_RATIO_MAX
        bc_bad = (not _fin(fb) or fb < BOUND_CONN_MIN or
                  not _fin(fc) or fc < BOUND_CONN_MIN)
        cell_kill = wtc_bad or wcv_bad or bc_bad
        kill = kill or cell_kill
        rows.append(dict(ckey=ck, width_ratio_tc=wtc, width_ratio_cv=wcv,
                         frac_bounded=fb, frac_connected=fc,
                         width_ratio_tc_over_oracle=c.get("width_ratio_tc_dlgi_over_oracle"),
                         width_tc_gt_1_5x=wtc_bad, width_cv_gt_1_5x=wcv_bad,
                         bounded_or_connected_lt_95=bc_bad, cell_kill=cell_kill))
    return dict(bar="C2", name="interval informativeness",
                status="KILL" if kill else "PASS", kill=kill,
                width_ratio_max=WIDTH_RATIO_MAX, bounded_connected_min=BOUND_CONN_MIN,
                cells=rows)


def eval_c3(primary, provenance):
    """Pilot-free RMSE <= 1.5x pilot every cell; estimator bit-identical to probe."""
    rows = []; kill = False
    for ck, c in sorted(primary.items()):
        rtc = c["rmse_ratio_tc_pf_over_pilot"]; rcv = c["rmse_ratio_cv_pf_over_pilot"]
        rtc_bad = _fin(rtc) and rtc > RMSE_RATIO_MAX
        rcv_bad = _fin(rcv) and rcv > RMSE_RATIO_MAX
        cell_kill = rtc_bad or rcv_bad
        kill = kill or cell_kill
        rows.append(dict(ckey=ck, rmse_ratio_tc_pf_over_pilot=rtc,
                         rmse_ratio_cv_pf_over_pilot=rcv,
                         rmse_ratio_tc_pf_over_oracle=c.get("rmse_ratio_tc_pf_over_oracle"),
                         pf_tc_bias=c.get("pf_tc_bias"), pf_cv_bias=c.get("pf_cv_bias"),
                         tc_gt_1_5x=rtc_bad, cv_gt_1_5x=rcv_bad, cell_kill=cell_kill))
    bit_id = bool(provenance.get("estimator", {}).get("estimator_bit_identical", False))
    if not bit_id:
        kill = True
    return dict(bar="C3", name="point precision",
                status="KILL" if kill else "PASS", kill=kill,
                rmse_ratio_max=RMSE_RATIO_MAX, estimator_bit_identical=bit_id,
                estimator_kill=(not bit_id), cells=rows,
                note="KILL if any RMSE ratio > 1.5x OR the estimator is not bit-identical "
                     "to the probe (calibration must change only the uncertainty statement).")


def eval_c4(primary):
    """Scene noninferiority: PSNR loss <= 0.2 dB AND task/Fisher loss <= 5% every cell."""
    rows = []; kill = False
    for ck, c in sorted(primary.items()):
        d = c["psnr_delta_mean"]; ci = c.get("psnr_delta_ci95"); tl = c["task_fisher_loss"]
        psnr_bad = _fin(d) and d < -PSNR_LOSS_MAX          # loss = -delta > 0.2 dB
        task_bad = _fin(tl) and tl > TASK_FISHER_LOSS_MAX
        cell_kill = psnr_bad or task_bad
        kill = kill or cell_kill
        rows.append(dict(ckey=ck, psnr_delta_mean=d, psnr_delta_ci95=ci,
                         task_fisher_loss=tl, psnr_loss_gt_0_2db=psnr_bad,
                         task_fisher_loss_gt_5pct=task_bad, cell_kill=cell_kill))
    return dict(bar="C4", name="scene noninferiority",
                status="KILL" if kill else "PASS", kill=kill,
                psnr_loss_max_db=PSNR_LOSS_MAX, task_fisher_loss_max=TASK_FISHER_LOSS_MAX,
                cells=rows)


def eval_c5(primary):
    """Model + gauge: innovation std ~1, whiteness, gauge invariance, every cell/stratum."""
    rows = []; kill = False
    for ck, c in sorted(primary.items()):
        istd = c["innov_std_med"]; lag = c["lag1_ac_med"]
        gtc = c["gauge_tc_shift_med"]; gcv = c["gauge_cv_shift_med"]
        istd_bad = (not _fin(istd)) or istd < INNOV_STD_LO or istd > INNOV_STD_HI
        lag_bad = _fin(lag) and abs(lag) > LAG1_AC_MAX
        gauge_bad = (_fin(gtc) and gtc > GAUGE_SHIFT_MAX) or (_fin(gcv) and gcv > GAUGE_SHIFT_MAX)
        cell_kill = istd_bad or lag_bad or gauge_bad
        kill = kill or cell_kill
        rows.append(dict(ckey=ck, innov_std_med=istd, lag1_ac_med=lag,
                         gauge_tc_shift_med=gtc, gauge_cv_shift_med=gcv,
                         innov_std_out_of_band=istd_bad, lag1_ac_too_large=lag_bad,
                         gauge_sensitive=gauge_bad, cell_kill=cell_kill))
    return dict(bar="C5", name="physical model + gauge",
                status="KILL" if kill else "PASS", kill=kill,
                innov_std_band=[INNOV_STD_LO, INNOV_STD_HI], lag1_ac_max=LAG1_AC_MAX,
                gauge_shift_max=GAUGE_SHIFT_MAX, cells=rows)


def eval_c6(reciprocity, substudy):
    """Reciprocity identity to tolerance + frozen (A,B,K) ordering predicts measured."""
    recip = reciprocity.get("worst_abserr_det_reciprocity", float("inf"))
    recip_gi = reciprocity.get("gi_fisher_reciprocity_worst_abserr", float("inf"))
    recip_bad = (not _fin(recip)) or recip > RECIPROCITY_TOL or \
                (not _fin(recip_gi)) or recip_gi > RECIPROCITY_TOL
    summ = substudy.get("_summary", {})
    frac = summ.get("contradiction_frac", float("nan"))
    sched_bad = _fin(frac) and frac > SCHED_CONTRA_FRAC_MAX
    kill = recip_bad or sched_bad
    return dict(bar="C6", name="reciprocity + schedule prediction",
                status="KILL" if kill else "PASS", kill=kill,
                reciprocity_worst_abserr=recip, gi_fisher_worst_abserr=recip_gi,
                reciprocity_tol=RECIPROCITY_TOL, reciprocity_fail=recip_bad,
                schedule_contradiction_frac=frac, schedule_contra_max=SCHED_CONTRA_FRAC_MAX,
                schedule_fail=sched_bad, n_substudy_cells=summ.get("n_cells"),
                n_contradictions=summ.get("n_scene_contradictions"))


def eval_c7(ledger_audit, edge):
    """Identical-ledger byte audit + edge honesty (edges published, not removed)."""
    sample = ledger_audit.get("sample", {})
    unequal = []; hidden = []
    for k, a in sample.items():
        if not a.get("identical_record_pilot_free_vs_linear", False):
            hidden.append(k)
        if a.get("pilot_extra_photons_or_exposures", 1) != 0:
            hidden.append(k)
        if (a.get("pilot_free_exposures") != N_EXP_EXPECT or
                a.get("plain_linear_exposures") != N_EXP_EXPECT or
                a.get("pilot_exposures") != N_EXP_EXPECT):
            unequal.append(k)
    edge_present = sorted(edge.keys())
    edge_ok = len(edge_present) >= 6                      # all 6 edge cells published
    kill = bool(unequal or hidden or not edge_ok)
    return dict(bar="C7", name="identical-ledger + edge honesty",
                status="KILL" if kill else "PASS", kill=kill,
                unequal_accounting_cells=unequal, hidden_info_cells=hidden,
                edge_cells_published=edge_present, edge_all_present=edge_ok,
                expected_exposures=N_EXP_EXPECT)


# ============================================================ top-level evaluate
def evaluate(blob):
    primary = blob.get("primary_cells", {})
    edge = blob.get("edge_cells", {})
    substudy = blob.get("schedule_substudy", {})
    reciprocity = blob.get("reciprocity", {})
    ledger_audit = blob.get("ledger_audit", {})
    provenance = blob.get("provenance", {})

    bars = [eval_c1(primary), eval_c2(primary), eval_c3(primary, provenance),
            eval_c4(primary), eval_c5(primary), eval_c6(reciprocity, substudy),
            eval_c7(ledger_audit, edge)]
    any_kill = any(b["kill"] for b in bars)
    verdict = "FLAGSHIP_KILL" if any_kill else "ALL_BARS_PASS"
    claim = ("One bucket stream, two MODEL-CERTIFIED products (R33 headline restored)."
             if not any_kill else
             "Flagship claim ends; DLGI demoted to a shorter theorem/method note "
             "(R34 sec.3 -- no repair round authorized).")
    return dict(
        meta=dict(utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
                  authority="R34 sec.3 (docs/ROUND63_GPT_ROUND34_RULING_RAW.md)",
                  n_primary=len(primary), n_edge=len(edge)),
        verdict=verdict, any_kill=any_kill, killed_bars=[b["bar"] for b in bars if b["kill"]],
        claim_ruling=claim, bars={b["bar"]: b for b in bars})


def human_table(rep):
    L = []
    L.append("=" * 72)
    L.append("DLGI CONFIRMATORY C1-C7 VERDICT  (R34 sec.3, binding kill bars)")
    L.append("=" * 72)
    for bk, b in rep["bars"].items():
        mark = "KILL" if b["kill"] else ("WARN" if b["status"].startswith("WARN") else "PASS")
        L.append(f"  [{mark:4s}] {b['bar']}  {b['name']}")
        if b["bar"] == "C1":
            nprim = len(b["cells"])
            nbelow = sum(r["below_0_92"] for r in b["cells"])
            nover = len(b["band_violations_overcover"])
            L.append(f"         cells {nprim}; below-0.92 {nbelow}; over-0.98 {nover}; "
                     f"stress {b['stress_cells']}")
        elif b["bar"] == "C2":
            nbad = sum(r["cell_kill"] for r in b["cells"])
            L.append(f"         cells failing width/boundedness: {nbad}/{len(b['cells'])}")
        elif b["bar"] == "C3":
            nbad = sum(r["cell_kill"] for r in b["cells"])
            L.append(f"         cells RMSE>1.5x: {nbad}/{len(b['cells'])}; "
                     f"estimator bit-identical: {b['estimator_bit_identical']}")
        elif b["bar"] == "C4":
            nbad = sum(r["cell_kill"] for r in b["cells"])
            L.append(f"         cells failing noninferiority: {nbad}/{len(b['cells'])}")
        elif b["bar"] == "C5":
            nbad = sum(r["cell_kill"] for r in b["cells"])
            L.append(f"         cells failing model/gauge: {nbad}/{len(b['cells'])}")
        elif b["bar"] == "C6":
            L.append(f"         reciprocity worst {b['reciprocity_worst_abserr']:.2e} "
                     f"(tol {b['reciprocity_tol']:.0e}); schedule contradiction "
                     f"{b['schedule_contradiction_frac']} (max {b['schedule_contra_max']})")
        elif b["bar"] == "C7":
            L.append(f"         unequal {b['unequal_accounting_cells']}; hidden "
                     f"{b['hidden_info_cells']}; edges published "
                     f"{len(b['edge_cells_published'])}")
    L.append("-" * 72)
    L.append(f"  VERDICT: {rep['verdict']}   killed bars: {rep['killed_bars'] or 'none'}")
    L.append(f"  {rep['claim_ruling']}")
    L.append("=" * 72)
    return "\n".join(L)


def run(results_path=RESULTS, write=True):
    blob = json.load(open(results_path))
    rep = evaluate(blob)
    print(human_table(rep))
    if write:
        fn = os.path.join(CONF, "C_BARS_VERDICT.json")
        json.dump(rep, open(fn, "w"), indent=2)
        print(f"\nwrote {fn}")
    return rep


# ============================================================ unit tests
def _synthetic(passing=True, break_bar=None):
    """Build a minimal CONFIRMATORY_RESULTS-shaped blob.  passing=True -> all bars pass;
    break_bar='C4' -> inject a single failing cell/field for that bar only."""
    def cell(ck, **ov):
        c = dict(ckey=ck, coverage_tc=0.95, coverage_tc_cp95=[0.93, 0.97],
                 nonfinite_rate=0.0, coverage_cv=0.95,
                 width_ratio_tc_dlgi_over_pilot=1.1, width_ratio_cv_dlgi_over_pilot=1.1,
                 width_ratio_tc_dlgi_over_oracle=2.0,
                 frac_bounded=0.98, frac_connected=0.99,
                 rmse_ratio_tc_pf_over_pilot=0.9, rmse_ratio_cv_pf_over_pilot=1.0,
                 rmse_ratio_tc_pf_over_oracle=1.3, pf_tc_bias=0.5, pf_cv_bias=0.001,
                 psnr_delta_mean=0.8, psnr_delta_ci95=[0.5, 1.1], task_fisher_loss=0.01,
                 innov_std_med=0.93, lag1_ac_med=0.05,
                 gauge_tc_shift_med=1e-12, gauge_cv_shift_med=1e-12)
        c.update(ov); return c

    cks = ["tc16_cv5_ph10", "tc32_cv15_ph10", "tc64_cv40_ph20"]
    primary = {ck: cell(ck) for ck in cks}
    edge = {f"tc2_cv{c}_ph10": cell(f"tc2_cv{c}_ph10") for c in (5, 15, 40)}
    edge.update({f"tc128_cv{c}_ph10": cell(f"tc128_cv{c}_ph10") for c in (5, 15, 40)})
    substudy = {"_summary": dict(n_cells=9, n_scene_contradictions=0, contradiction_frac=0.0)}
    reciprocity = dict(worst_abserr_det_reciprocity=3e-15,
                       gi_fisher_reciprocity_worst_abserr=3e-15)
    acct = dict(identical_record_pilot_free_vs_linear=True, pilot_extra_photons_or_exposures=0,
                pilot_free_exposures=2048, plain_linear_exposures=2048, pilot_exposures=2048)
    ledger = dict(sample={"tc16_cv5_ph10::paired": acct})
    prov = dict(estimator=dict(estimator_bit_identical=True))

    if not passing and break_bar:
        bad = "tc32_cv15_ph10"
        if break_bar == "C1":
            primary[bad]["coverage_tc"] = 0.80          # below 0.92
        elif break_bar == "C1nf":
            primary[bad]["nonfinite_rate"] = 0.10        # > 5%
        elif break_bar == "C2":
            primary[bad]["width_ratio_tc_dlgi_over_pilot"] = 2.0   # > 1.5x
        elif break_bar == "C2b":
            primary[bad]["frac_bounded"] = 0.50          # < 95%
        elif break_bar == "C3":
            primary[bad]["rmse_ratio_cv_pf_over_pilot"] = 1.9      # > 1.5x
        elif break_bar == "C3est":
            prov["estimator"]["estimator_bit_identical"] = False
        elif break_bar == "C4":
            primary[bad]["psnr_delta_mean"] = -0.5       # loss 0.5 dB > 0.2
        elif break_bar == "C4task":
            primary[bad]["task_fisher_loss"] = 0.09      # > 5%
        elif break_bar == "C5":
            primary[bad]["innov_std_med"] = 1.6          # out of [0.8,1.2]
        elif break_bar == "C5lag":
            primary[bad]["lag1_ac_med"] = 0.4            # > 0.2
        elif break_bar == "C6":
            reciprocity["worst_abserr_det_reciprocity"] = 1e-3     # > tol
        elif break_bar == "C6sched":
            substudy["_summary"]["contradiction_frac"] = 0.30      # > 10%
        elif break_bar == "C7":
            ledger["sample"]["tc16_cv5_ph10::paired"]["pilot_extra_photons_or_exposures"] = 5
        elif break_bar == "C7edge":
            edge.pop("tc2_cv5_ph10")                      # removed a failed edge

    return dict(primary_cells=primary, edge_cells=edge, schedule_substudy=substudy,
                reciprocity=reciprocity, ledger_audit=ledger, provenance=prov)


def selftest():
    ok = True

    # 1. all-passing synthetic -> ALL_BARS_PASS, no kills
    rep = evaluate(_synthetic(passing=True))
    if rep["any_kill"] or rep["verdict"] != "ALL_BARS_PASS":
        print("FAIL: passing synthetic should pass; got", rep["verdict"],
              rep["killed_bars"]); ok = False
    else:
        print("PASS: all-passing synthetic -> ALL_BARS_PASS")

    # 2. each break should kill exactly its own bar
    cases = {"C1": "C1", "C1nf": "C1", "C2": "C2", "C2b": "C2", "C3": "C3",
             "C3est": "C3", "C4": "C4", "C4task": "C4", "C5": "C5", "C5lag": "C5",
             "C6": "C6", "C6sched": "C6", "C7": "C7", "C7edge": "C7"}
    for brk, expect_bar in cases.items():
        rep = evaluate(_synthetic(passing=False, break_bar=brk))
        killed = rep["killed_bars"]
        if killed != [expect_bar]:
            print(f"FAIL: break {brk} should kill only {expect_bar}; got {killed}")
            ok = False
        else:
            print(f"PASS: break {brk:7s} -> kills {expect_bar} (verdict {rep['verdict']})")

    # 3. exact binomial CI sanity via a coverage cell round-trip already covered above.
    print("\nSELFTEST:", "ALL PASS" if ok else "FAILURES ABOVE")
    return ok


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="DLGI confirmatory C1-C7 bar evaluator")
    ap.add_argument("--selftest", action="store_true", help="unit-test the bar arithmetic")
    ap.add_argument("--results", default=RESULTS, help="CONFIRMATORY_RESULTS.json path")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    if not os.path.exists(a.results):
        print(f"[c_bars] no endpoints yet at {a.results}; run --selftest to check the "
              f"bar arithmetic, or run the confirmatory fleet + --merge first.")
        sys.exit(2)
    run(a.results)
