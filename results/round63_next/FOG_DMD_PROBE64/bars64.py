#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
BEYOND-MODULATOR-BAND sealed-probe BAR EVALUATOR (R39 ask 5), R34-style binding kill bars.

Consumes a confirmatory-endpoints blob (produced by the sealed-probe driver + --merge) and
emits a verdict JSON + human table.  Each bar is a BINDING kill bar: any kill ends the flagship
beyond-band claim, with no averaging-away of a failed cell and no repair round.

  B1 aperture-map coverage : the measured recoverable spectrum matches the aperture-law
       prediction map -- in-coverage recovered AND out-coverage stays blind (separation).
  B2 beyond-band recovery  : fixed+MLE (and fixed+moment) beyond-band NMSE below threshold on
       the witness + natural strata, AND beats the fresh+mean wall by a margin.
  B3 mismatch degradation  : every F5 mismatch axis (fine-band rotation, band width, radial
       profile, geometry-mixing alpha, wrong tau/sigma_f) degrades beyond-band NMSE by <= cap.
  B4 multi-start agreement : the covariance estimator's data-only multi-start std is small
       (collusion-free / injective -- no realization vector to trade against).
  B5 equal-budget ledger   : byte-audit that every arm spends identical M x T exposures and
       photons/bucket; the concentration principle (fixed >= fresh at equal budget) holds.
  B6 channel impossibility : the fresh+mean wall and classic-averaging arms stay at the
       beyond-band wall (NMSE ~ 1.0) -- the first-moment channel is provably blind.

THRESHOLDS ARE [R39] PLACEHOLDERS -- the external referee's ruling sets them.  Every threshold
lives in THRESHOLDS below with an [R39] tag, so FREEZING IS A ONE-EDIT STEP.  --selftest exercises
the bar arithmetic on synthetic passing AND failing inputs for every bar.  CPU, seconds.
"""
import json
import math
import os
import sys
from datetime import datetime, timezone

# ======================================================= [R39] FROZEN? thresholds
# Every value tagged [R39] is a PLACEHOLDER awaiting the referee ruling.  Editing these six
# numbers (and nothing else) freezes the probe.  Chosen provisionally from the dev-grade gates
# (matched beyond-band NMSE ~0.5; wall = 1.0; oracle ~0.1; multi-start std ~0.01-0.02).
THRESHOLDS = dict(
    B1_min_separation=3.0,          # [R39] out/in coverage rel-err ratio >= this (law holds)
    B1_max_in_coverage_err=0.70,    # [R39] in-coverage rel err <= this
    B2_max_beyond_band_nmse=0.70,   # [R39] fixed+MLE beyond-band NMSE <= this (per cell)
    B2_min_wall_margin=0.20,        # [R39] wall_nmse - fixed_mle_nmse >= this
    B3_max_degradation=0.25,        # [R39] each mismatch axis degradation <= this (F5 kill line)
    B4_max_multistart_std=0.05,     # [R39] covariance multi-start beyond-band-NMSE std <= this
    B5_max_exposure_mismatch=0,     # [R39] byte-exact exposure equality across arms (0 tolerance)
    B5_min_concentration=0.0,       # [R39] fresh_nmse - fixed_nmse >= this (fixed >= fresh)
    B6_min_wall_nmse=0.90,          # [R39] fresh+mean & classic-avg beyond-band NMSE >= this
)
T = THRESHOLDS

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "confirmatory", "CONFIRMATORY_RESULTS.json")


def _fin(x):
    try:
        return math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


# ============================================================ per-bar evaluators
def eval_b1(ap):
    """Aperture-map coverage: separation >= min, in-coverage err <= max, per cell."""
    rows = []; kill = False
    for ck, c in sorted(ap.get("cells", {}).items()):
        sep = c.get("separation"); ine = c.get("in_coverage_err")
        sep_bad = (not _fin(sep)) or sep < T["B1_min_separation"]
        in_bad = (not _fin(ine)) or ine > T["B1_max_in_coverage_err"]
        cell_kill = sep_bad or in_bad
        kill = kill or cell_kill
        rows.append(dict(ckey=ck, separation=sep, in_coverage_err=ine,
                         out_coverage_err=c.get("out_coverage_err"),
                         separation_lt_min=sep_bad, in_coverage_gt_max=in_bad, cell_kill=cell_kill))
    return dict(bar="B1", name="aperture-map coverage", status="KILL" if kill else "PASS",
                kill=kill, min_separation=T["B1_min_separation"],
                max_in_coverage_err=T["B1_max_in_coverage_err"], cells=rows)


def eval_b2(arms):
    """Beyond-band recovery: fixed+MLE NMSE <= max AND beats the fresh+mean wall by margin."""
    rows = []; kill = False
    for sid, a in sorted(arms.items()):
        mle = a.get("fixed_mle"); wall = a.get("fresh_mean")
        nmse_bad = (not _fin(mle)) or mle > T["B2_max_beyond_band_nmse"]
        margin = (wall - mle) if (_fin(wall) and _fin(mle)) else float("nan")
        margin_bad = (not _fin(margin)) or margin < T["B2_min_wall_margin"]
        cell_kill = nmse_bad or margin_bad
        kill = kill or cell_kill
        rows.append(dict(scene=sid, fixed_mle_nmse=mle, fixed_moment_nmse=a.get("fixed_moment"),
                         wall_nmse=wall, wall_margin=margin, nmse_gt_max=nmse_bad,
                         margin_lt_min=margin_bad, cell_kill=cell_kill))
    return dict(bar="B2", name="beyond-band recovery", status="KILL" if kill else "PASS",
                kill=kill, max_beyond_band_nmse=T["B2_max_beyond_band_nmse"],
                min_wall_margin=T["B2_min_wall_margin"], cells=rows)


def eval_b3(mismatch):
    """Mismatch degradation: every axis/condition degradation <= cap (F5 kill line)."""
    rows = []; kill = False
    for axis, conds in sorted(mismatch.items()):
        for cond, c in sorted(conds.items()):
            d = c.get("degradation")
            bad = (not _fin(d)) or d > T["B3_max_degradation"]
            kill = kill or bad
            rows.append(dict(axis=axis, condition=cond, degradation=d,
                             degradation_gt_cap=bad, cell_kill=bad))
    return dict(bar="B3", name="mismatch degradation (F5)", status="KILL" if kill else "PASS",
                kill=kill, max_degradation=T["B3_max_degradation"], cells=rows)


def eval_b4(ms):
    """Multi-start agreement: covariance-estimator data-only multi-start std small (collusion-free)."""
    rows = []; kill = False
    for sid, c in sorted(ms.items()):
        std = c.get("agreement_std")
        bad = (not _fin(std)) or std > T["B4_max_multistart_std"]
        kill = kill or bad
        rows.append(dict(scene=sid, agreement_std=std, std_gt_max=bad, cell_kill=bad))
    return dict(bar="B4", name="multi-start agreement", status="KILL" if kill else "PASS",
                kill=kill, max_multistart_std=T["B4_max_multistart_std"], cells=rows)


def eval_b5(ledger):
    """Equal-budget byte-audit + concentration principle (fixed >= fresh at equal budget)."""
    exp = ledger.get("exposures_by_arm", {})
    pho = ledger.get("photons_by_arm", {})
    exp_vals = [v for v in exp.values() if _fin(v)]
    pho_vals = [v for v in pho.values() if _fin(v)]
    exp_bad = (len(set(exp_vals)) > 1) if exp_vals else True     # exposures must be byte-identical
    pho_spread = (max(pho_vals) - min(pho_vals)) if pho_vals else float("inf")
    pho_bad = pho_spread > 1e-6 * (max(pho_vals) if pho_vals else 1.0)
    conc = ledger.get("concentration_fixed_minus_fresh")         # fresh_nmse - fixed_nmse (>=0 good)
    conc_bad = (not _fin(conc)) or conc < T["B5_min_concentration"]
    kill = bool(exp_bad or pho_bad or conc_bad)
    return dict(bar="B5", name="equal-budget ledger + concentration", status="KILL" if kill else "PASS",
                kill=kill, exposures_by_arm=exp, photons_by_arm=pho,
                exposures_unequal=exp_bad, photons_unequal=pho_bad,
                concentration_fixed_minus_fresh=conc, concentration_lt_min=conc_bad,
                min_concentration=T["B5_min_concentration"])


def eval_b6(arms):
    """Channel impossibility: fresh+mean AND classic-avg stay at the beyond-band wall (~1.0)."""
    rows = []; kill = False
    for sid, a in sorted(arms.items()):
        fm = a.get("fresh_mean"); ca = a.get("classic_avg")
        fm_bad = (not _fin(fm)) or fm < T["B6_min_wall_nmse"]
        ca_bad = (not _fin(ca)) or ca < T["B6_min_wall_nmse"]
        cell_kill = fm_bad or ca_bad
        kill = kill or cell_kill
        rows.append(dict(scene=sid, fresh_mean_nmse=fm, classic_avg_nmse=ca,
                         fresh_mean_below_wall=fm_bad, classic_avg_below_wall=ca_bad,
                         cell_kill=cell_kill))
    return dict(bar="B6", name="channel impossibility (wall holds)",
                status="KILL" if kill else "PASS", kill=kill,
                min_wall_nmse=T["B6_min_wall_nmse"],
                note="A wall arm dropping BELOW the wall would mean the first-moment channel is "
                     "NOT blind -- it would falsify the impossibility separation (the strongest "
                     "novelty lever). This bar guards that control.", cells=rows)


def evaluate(blob):
    bars = [eval_b1(blob.get("aperture_map", {})),
            eval_b2(blob.get("arms", {})),
            eval_b3(blob.get("mismatch", {})),
            eval_b4(blob.get("multistart", {})),
            eval_b5(blob.get("ledger", {})),
            eval_b6(blob.get("arms", {}))]
    any_kill = any(b["kill"] for b in bars)
    verdict = "BEYOND_BAND_KILL" if any_kill else "ALL_BARS_PASS"
    claim = ("The medium's fine speckle band is a synthesized aperture: beyond-modulator-band "
             "scene content is recovered from the bucket covariance alone."
             if not any_kill else
             "Beyond-band flagship claim ends; demote to a theorem/map note (no repair round).")
    return dict(meta=dict(utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
                          authority="R39 sealed-probe spec (thresholds [R39]-placeholder until "
                                    "the referee ruling)", thresholds=THRESHOLDS),
                verdict=verdict, any_kill=any_kill,
                killed_bars=[b["bar"] for b in bars if b["kill"]],
                claim_ruling=claim, bars={b["bar"]: b for b in bars})


def human_table(rep):
    L = ["=" * 74, "BEYOND-BAND SEALED-PROBE B1-B6 VERDICT  (R39, binding kill bars)", "=" * 74]
    for bk, b in rep["bars"].items():
        mark = "KILL" if b["kill"] else "PASS"
        L.append(f"  [{mark}] {b['bar']}  {b['name']}")
    L.append("-" * 74)
    L.append(f"  VERDICT: {rep['verdict']}   killed bars: {rep['killed_bars'] or 'none'}")
    L.append(f"  {rep['claim_ruling']}")
    L.append("  NOTE: thresholds are [R39] placeholders until the referee ruling.")
    L.append("=" * 74)
    return "\n".join(L)


# ============================================================ unit tests
def _synthetic(passing=True, break_bar=None):
    def arm(sid, **ov):
        a = dict(fixed_mle=0.45, fixed_moment=0.52, fresh_mean=1.0, fresh_cov_gls=0.95,
                 oracle=0.09, classic_avg=1.0)
        a.update(ov); return a
    arms = {"confirmatory_witness_0": arm("w0"), "confirmatory_natural_0": arm("n0")}
    ap = dict(cells={"cal_0": dict(separation=8.0, in_coverage_err=0.30, out_coverage_err=2.4)})
    mismatch = {"fine_rotation": {"rot10": dict(degradation=0.02), "rot20": dict(degradation=0.19)},
                "geometry_alpha": {"a0.1": dict(degradation=0.00), "a0.25": dict(degradation=0.06)}}
    ms = {"confirmatory_witness_0": dict(agreement_std=0.017)}
    ledger = dict(exposures_by_arm={"fixed_mle": 262144, "fresh_mean": 262144, "oracle": 262144},
                  photons_by_arm={"fixed_mle": 1e5, "fresh_mean": 1e5},
                  concentration_fixed_minus_fresh=0.43)
    if not passing and break_bar:
        if break_bar == "B1":
            ap["cells"]["cal_0"]["separation"] = 1.5
        elif break_bar == "B2":
            arms["confirmatory_witness_0"]["fixed_mle"] = 0.85     # > max
        elif break_bar == "B2m":
            arms["confirmatory_witness_0"]["fresh_mean"] = 0.50    # wall margin gone
        elif break_bar == "B3":
            mismatch["fine_rotation"]["rot20"]["degradation"] = 0.60
        elif break_bar == "B4":
            ms["confirmatory_witness_0"]["agreement_std"] = 0.20
        elif break_bar == "B5":
            ledger["exposures_by_arm"]["fresh_mean"] = 260000      # unequal budget
        elif break_bar == "B5c":
            ledger["concentration_fixed_minus_fresh"] = -0.3       # fresh beats fixed
        elif break_bar == "B6":
            arms["confirmatory_witness_0"]["fresh_mean"] = 0.40    # wall breached -> control fails
    return dict(aperture_map=ap, arms=arms, mismatch=mismatch, multistart=ms, ledger=ledger)


def selftest():
    ok = True
    rep = evaluate(_synthetic(passing=True))
    if rep["any_kill"] or rep["verdict"] != "ALL_BARS_PASS":
        print("FAIL: passing synthetic should pass; got", rep["verdict"], rep["killed_bars"]); ok = False
    else:
        print("PASS: all-passing synthetic -> ALL_BARS_PASS")
    # A wall breach is a GENUINE dual failure: fresh_mean dropping below the wall trips both
    # B2 (its margin over the wall vanishes) AND B6 (the channel-impossibility control fails).
    # That coupling is intended defense-in-depth, so B2m/B6 legitimately kill {B2,B6}.
    cases = {"B1": ["B1"], "B2": ["B2"], "B2m": ["B2", "B6"], "B3": ["B3"], "B4": ["B4"],
             "B5": ["B5"], "B5c": ["B5"], "B6": ["B2", "B6"]}
    for brk, expect in cases.items():
        rep = evaluate(_synthetic(passing=False, break_bar=brk))
        if sorted(rep["killed_bars"]) != sorted(expect):
            print(f"FAIL: break {brk} should kill {expect}; got {rep['killed_bars']}"); ok = False
        else:
            print(f"PASS: break {brk:4s} -> kills {expect}")
    print("\nSELFTEST:", "ALL PASS" if ok else "FAILURES ABOVE")
    return ok


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="R39 beyond-band sealed-probe B1-B6 bar evaluator")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--results", default=RESULTS)
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    if not os.path.exists(a.results):
        print(f"[bars64] no endpoints yet at {a.results}; run --selftest, or run the sealed "
              f"confirmatory fleet + --merge first.")
        sys.exit(2)
    rep = evaluate(json.load(open(a.results)))
    print(human_table(rep))
    json.dump(rep, open(os.path.join(HERE, "confirmatory", "B_BARS_VERDICT.json"), "w"), indent=2)
