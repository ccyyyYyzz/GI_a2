#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
BINDING BARS D0-D7 + KILL TREE for the sealed detection probe (R41 sec 4.6-4.7).

Thresholds are FROZEN in sealed_common.BARS (R41 froze them; no placeholders). This module is the
evaluator: each eval_Dk() takes MEASURED quantities and returns pass/fail + the failing sub-checks;
kill_tree() walks the R41 sec 4.7 branch logic (incl. the sec 4.5 fresh-vs-fixed 1.20x rule).

D0 and the calibration-side D1 are run now (dry run); D2-D7 are evaluated on SEALED banks after the
coordinator's freeze order. eval_* are pure functions so the same code scores the sealed run verbatim
(D7 = no repair: one sealed run, no threshold retuning after unblinding).
"""
import numpy as np

import sealed_common as sc

B = sc.BARS


def _wilson_lb(k, nrec, z=1.645):
    """One-sided Wilson lower confidence bound on a pass fraction (z=1.645 -> 95%)."""
    if nrec == 0:
        return 0.0
    p = k / nrec
    d = 1 + z ** 2 / nrec
    c = p + z ** 2 / (2 * nrec)
    hw = z * np.sqrt(p * (1 - p) / nrec + z ** 2 / (4 * nrec ** 2))
    return float((c - hw) / d)


# ================================================================== D0
def eval_D0(mean_deriv_rel, two_engine_rel, mc_shot_rel):
    """mean_deriv_rel: max over beyond-band events of ||dm||/||.||; two_engine_rel: max rel diff on
    {lam_mean,lam_max,tr,d'}; mc_shot_rel: max rel diff of MC shot/cov vs |P|x ledger."""
    c = dict(
        mean_derivative=(mean_deriv_rel <= B["D0"]["mean_deriv_rel_max"]),
        two_engine_fisher=(two_engine_rel <= B["D0"]["two_engine_rel_max"]),
        mc_shot_ledger=(mc_shot_rel <= B["D0"]["mc_shot_rel_max"]),
    )
    return dict(bar="D0", checks=c, passed=all(c.values()),
                measured=dict(mean_deriv_rel=mean_deriv_rel, two_engine_rel=two_engine_rel,
                              mc_shot_rel=mc_shot_rel))


# ================================================================== D1
def eval_D1(ratios):
    """ratios: list of d'_emp/d'_analytic across the 27 primary cells x anomaly sizes."""
    r = np.asarray(ratios, float)
    are = np.abs(r - 1.0)
    c = dict(
        band_0p80_1p20=bool(np.all((r >= B["D1"]["ratio_lo"]) & (r <= B["D1"]["ratio_hi"]))),
        median_are_le_10pct=bool(np.median(are) <= B["D1"]["median_are_max"]),
        none_outside_30pct=bool(np.all((r >= B["D1"]["hard_lo"]) & (r <= B["D1"]["hard_hi"]))),
    )
    return dict(bar="D1", checks=c, passed=all(c.values()),
                measured=dict(n=len(r), median_ratio=float(np.median(r)),
                              min_ratio=float(r.min()), max_ratio=float(r.max()),
                              median_abs_rel_err=float(np.median(are))))


# ================================================================== D2
def eval_D2(eps2_pass_cells, eps2_mc_pass_frac, eps2_mc_nrec, best_cell_Tdet,
            eps5_spread_all81, eps5_worstmode_frac, eps1_best_Tdet, eps05_max_dprime):
    d = B["D2"]
    mc_lb = _wilson_lb(int(round(eps2_mc_pass_frac * eps2_mc_nrec)), eps2_mc_nrec)
    c = dict(
        eps2_77of81=(eps2_pass_cells >= d["eps2_min_cells"]),
        eps2_mc_lb_gt_0p90=(mc_lb > d["eps2_mc_lb_min"]),
        best_cell_Tdet_le_600=(best_cell_Tdet <= d["best_cell_Tdet_max"]),
        eps5_spread_all81=(eps5_spread_all81 >= d["eps2_total_cells"]),
        eps5_worstmode_ge_90pct=(eps5_worstmode_frac >= d["eps5_worstmode_frac"]),
        eps1_best_le_2048=(eps1_best_Tdet <= d["eps1_best_Tdet_max"]),
    )
    edge = dict(eps05_audit_only=True, eps05_max_dprime=eps05_max_dprime,
                eps05_promotable=(eps05_max_dprime >= d["eps05_dp_floor"]))
    # kill = failure of the 2% headline OR the 5% robustness bar
    headline_ok = c["eps2_77of81"] and c["eps2_mc_lb_gt_0p90"]
    robust_ok = c["eps5_spread_all81"] and c["eps5_worstmode_ge_90pct"]
    return dict(bar="D2", checks=c, edge_audit=edge, passed=all(c.values()),
                kill=(not (headline_ok and robust_ok)),
                measured=dict(eps2_pass_cells=eps2_pass_cells, eps2_mc_lb=round(mc_lb, 4),
                              best_cell_Tdet=best_cell_Tdet, eps1_best_Tdet=eps1_best_Tdet))


# ================================================================== D3
def eval_D3(target_tpr, fa_inband, fa_amp, fa_tau, bal_acc, nontarget_beyond_dp,
            intended_dp, simplex_max_offdiag_cc, tau_aliases_beyond):
    d = B["D3"]
    c = dict(
        target_tpr_ge_90=(target_tpr >= d["target_tpr_min"]),
        fa_inband_le_5=(fa_inband <= d["offtarget_fa_max"]),
        fa_amplitude_le_5=(fa_amp <= d["offtarget_fa_max"]),
        fa_timescale_le_5=(fa_tau <= d["offtarget_fa_max"]),
        balanced_accuracy_ge_90=(bal_acc >= d["bal_acc_min"]),
        beyond_nontarget_dp_le_0p5=(max(np.abs(nontarget_beyond_dp)) <= d["nontarget_dp_max"]),
        intended_scores_dp_ge_5=(min(intended_dp) >= d["intended_dp_min"]),
        simplex_offdiag_cc_lt_0p10=(simplex_max_offdiag_cc < d["simplex_canon_corr_max"]),
    )
    return dict(bar="D3", checks=c, passed=all(c.values()),
                kill=bool(tau_aliases_beyond),   # tau/mixed aliasing as beyond kills the sentinel
                measured=dict(target_tpr=target_tpr, bal_acc=bal_acc,
                              simplex_max_offdiag_cc=simplex_max_offdiag_cc))


# ================================================================== D4 (fresh-cov branch)
def eval_D4(fixed_latency, fresh_latency, fixed_passes_bars, fresh_passes_bars):
    """R41 sec 4.5 frozen branch: compare FIXED-COV vs FRESH-COV-OPT on identical photons/exposures/
    duration/anomaly/law. latency = detection banks (T_det) over the primary grid (use medians)."""
    ratio = fixed_latency / max(fresh_latency, 1e-9)
    if not (fixed_passes_bars or fresh_passes_bars):
        branch = "BOTH_FAIL_kill_constructive_covariance"
    elif ratio <= B["D4"]["fixed_retain_ratio"]:
        branch = "RETAIN_fixed_concentration"        # fixed no worse than 1.20x fresh
    else:
        branch = "SWITCH_to_fresh_production"         # fresh materially better -> production design
    return dict(bar="D4", latency_ratio_fixed_over_fresh=round(ratio, 3), branch=branch,
                passed=(branch != "BOTH_FAIL_kill_constructive_covariance"),
                measured=dict(fixed_latency=fixed_latency, fresh_latency=fresh_latency))


# ================================================================== D5 (mismatch)
def eval_D5(auc_loss, Tdet_inflation, nontarget_fa, level_label):
    d = B["D5"]
    c = dict(
        auc_loss_le_0p05=(auc_loss <= d["auc_loss_max"]),
        Tdet_inflation_le_25=(Tdet_inflation <= d["Tdet_inflation_max"]),
        nontarget_fa_le_5=(nontarget_fa <= d["nontarget_fa_max"]),
    )
    return dict(bar="D5", level=level_label, checks=c, passed=all(c.values()),
                measured=dict(auc_loss=auc_loss, Tdet_inflation=Tdet_inflation, nontarget_fa=nontarget_fa))


# ================================================================== D6 (online feasibility)
def eval_D6(online_ms, mem_mb):
    d = B["D6"]
    frac = online_ms / sc.MS_PER_BANK
    c = dict(
        online_le_10pct_bank=(frac <= d["online_frac_of_bank"]),
        memory_le_10mb=(mem_mb <= d["mem_mb_max"]),
    )
    return dict(bar="D6", checks=c, passed=all(c.values()),
                measured=dict(online_ms=online_ms, online_frac_of_bank=round(frac, 4), mem_mb=mem_mb))


# ================================================================== KILL TREE (R41 sec 4.7)
def kill_tree(results):
    """Walk the R41 sec 4.7 branch logic; return the terminal verdict + the stopping node."""
    D = {r["bar"]: r for r in results if "bar" in r}

    def failed(bar):
        return (bar in D) and (not D[bar]["passed"])

    if failed("D0"):
        return dict(node="D0", verdict="STOP -- correction/theory note only (mechanism/engine integrity)")
    if failed("D1"):
        return dict(node="D1", verdict="STOP -- Fisher map not trustworthy")
    if "D2" in D and D["D2"].get("kill"):
        return dict(node="D2", verdict="remove detection tier; return to R40")
    if "D3" in D and D["D3"].get("kill"):
        return dict(node="D3", verdict="binary detector may survive in supplement; flagship sentinel killed")
    if "D4" in D:
        br = D["D4"]["branch"]
        if br == "BOTH_FAIL_kill_constructive_covariance":
            return dict(node="D4", verdict="kill constructive covariance capability")
        if br == "SWITCH_to_fresh_production":
            return dict(node="D4", verdict="switch frozen design branch to FRESH-COV-OPT; kill fixed-code concentration only (wall+sentinel survive)")
    if failed("D5"):
        return dict(node="D5", verdict="narrow the physical domain or kill robustness wording")
    if failed("D6"):
        return dict(node="D6", verdict="retain offline method only; delete real-time claim")
    return dict(node="all_pass", verdict="DETECTION ADMITTED as the flagship positive")


if __name__ == "__main__":
    # demonstrate the kill-tree logic on (a) all-pass and (b) a D2-fail mock
    allpass = [
        eval_D0(2e-16, 0.015, 0.06),
        eval_D1([1.03, 0.98, 1.10, 0.92, 1.05] * 6),
        eval_D2(78, 0.94, 500, 467, 81, 0.92, 1869, 2.1),
        eval_D3(0.94, 0.03, 0.02, 0.04, 0.93, [0.3, -0.2, 0.4], [7, 6, 5.5], 0.035, False),
        eval_D4(480, 505, True, True),
        eval_D5(0.03, 0.20, 0.04, "rot10%"),
        eval_D6(0.8, 4.2),
    ]
    print("[all-pass mock] kill tree ->", kill_tree(allpass)["verdict"])
    for r in allpass:
        print(f"  {r['bar']}: passed={r['passed']}" +
              (f" branch={r.get('branch')}" if r["bar"] == "D4" else ""))
    d2fail = list(allpass)
    d2fail[2] = eval_D2(60, 0.70, 500, 900, 70, 0.5, 3000, 1.0)   # headline+robustness fail
    print("\n[D2-fail mock] kill tree ->", kill_tree(d2fail)["node"], "::", kill_tree(d2fail)["verdict"])
    d0fail = list(allpass)
    d0fail[0] = eval_D0(1e-3, 0.02, 0.05)
    print("[D0-fail mock] kill tree ->", kill_tree(d0fail)["node"], "::", kill_tree(d0fail)["verdict"])
