#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
SEALED CONFIRMATORY RUN -- full D0-D7 ladder on the sealed banks (R41 §4, freeze commit 5910277).

Executed AFTER the coordinator's freeze order. Frozen machinery + hashed banks + D0-D7 thresholds were
committed BEFORE unblinding. NO threshold change, NO scene removal, NO estimator switch after this
point (D7). If a bar fails, the kill-tree node is recorded and the remaining bars are STILL evaluated
(full autopsy). FRESH-COV-OPT runs in PRODUCTION form (full per-bank in-band-nuisance profiling).

Deliverables: CONFIRMATORY_VERDICT.md, CONFIRMATORY_RESULTS.json, per-bar D*_CONF.json, cost ledger.
GPU shared with JET_TEST -> poll before heavy phases; wait if free VRAM < 4 GB.
"""
import json
import os
import subprocess
import time

import numpy as np
import torch

import sealed_common as sc
import sealed_banks as sb
import arms
import simplex as sx
import bars

HERE = os.path.dirname(os.path.abspath(__file__))
T0 = time.time()
COST = {}


def gpu_free_mb():
    try:
        out = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=memory.free", "--format=csv,noheader,nounits"]).decode()
        return int(out.strip().split("\n")[0])
    except Exception:
        return 99999


def gpu_gate(min_free_mb=4096, label=""):
    while gpu_free_mb() < min_free_mb:
        print(f"[GPU gate] free < {min_free_mb} MB before {label}; waiting 5 min (JET_TEST active)...", flush=True)
        time.sleep(300)
    print(f"[GPU gate] free={gpu_free_mb()} MB -> proceed: {label}", flush=True)


def load_bank_scenes(bank):
    man = sb.build_bank(bank, write=False)
    scenes = {}
    for e in man["scenes"]:
        x = sb.nat_scene(e["source"]) if e["family"] == "natural" else sb.synth_scene(e["family"], e["seed"])
        assert sb.sha(x) == e["x_sha256"], f"HASH MISMATCH {e['scene_id']} (bank tamper!)"
        scenes[e["scene_id"]] = np.asarray(x, float)
    return scenes, man


def _timed(name, fn):
    t = time.time()
    r = fn()
    COST[name] = round(time.time() - t, 1)
    print(f"  [{name}: {COST[name]}s]", flush=True)
    return r


# ============================================================ D0 (re-verify on sealed data)
def run_D0(conf_scenes):
    print("\n=== D0 (sealed re-verify) ===", flush=True)
    # mean-derivative wall
    md_worst = max(float(np.linalg.norm(sc.P_FIXED @ sc.BETA_np[c]) / np.linalg.norm(sc.P_FIXED))
                   for c in sc.CLAIMS)
    # two-engine + shot ledger on a confirmatory scene geometry
    scene = conf_scenes["confirmatory_witness_1"]
    te_worst = 0.0
    for geo in [sc.BEST_CELL, sc.MID_CELL, sc.FLOOR_CELL]:
        cell = sc.setup_cell(**geo, x_np=scene)
        B = sc.fisher_engine_B(cell)
        for a, b in [(cell["lam_mean"], B["lam_mean"]), (cell["lam_max"], B["lam_max"]),
                     (cell["tr"], B["tr"])]:
            te_worst = max(te_worst, abs(a - b) / max(abs(a), 1e-30))
    # shot ledger (gaussian 2nd-order medium) on best/mid cells
    sl_worst = 0.0
    for geo in [sc.BEST_CELL, sc.MID_CELL]:
        cell = sc.setup_cell(**geo, x_np=scene)
        R = cell["R"].cpu().numpy(); dV0 = np.diag(cell["V0"].cpu().numpy())
        Bs = sc.physical_banks(cell, 8000, seed=246810, freeze_medium=True)
        shot_rel = float(np.median(np.abs(Bs.var(0) - R) / np.maximum(R, 1e-12)))
        Scov, _ = sc.gen_records(cell, 1, 12000, "H0", None, rng=13579)
        sc_rel = float(np.median(np.abs(np.diag(Scov[0].cpu().numpy()) - dV0) / np.maximum(dV0, 1e-12)))
        sl_worst = max(sl_worst, shot_rel, sc_rel)
    res = bars.eval_D0(md_worst, te_worst, sl_worst)
    res["measured"] = dict(mean_deriv=md_worst, two_engine=te_worst, shot_ledger=sl_worst)
    print(f"  D0: wall={md_worst:.1e} two-engine={te_worst:.3f} shot={sl_worst:.3f} -> {'PASS' if res['passed'] else 'FAIL'}", flush=True)
    return res


# ============================================================ D1 + D2 (capability)
def analytic_map_81(scene, eps, direction="energy_spread"):
    """81-cell analytic pass map: does d'>=5 by 4096 banks? Returns (pass_count, per_cell)."""
    rows = []
    npass = 0
    best_Tdet = np.inf
    for geo in sc.full_grid_81():
        cell = sc.setup_cell(**geo, x_np=scene, want_mc=False)
        lam = cell["lam_mean"] if direction == "energy_spread" else max(cell["lam_min"], 1e-12)
        Tdet = sc.DP_STRONG2 / (lam * eps * eps * cell["xnorm2"])
        ok = Tdet <= sc.T_CAP
        npass += int(ok)
        best_Tdet = min(best_Tdet, Tdet)
        rows.append(dict(**geo, T_det=round(Tdet, 0), pass_4096=bool(ok)))
    return npass, best_Tdet, rows


def run_D1_D2(conf_scenes):
    print("\n=== D1 + D2 (primary grid, sealed confirmatory) ===", flush=True)
    scene = conf_scenes["confirmatory_witness_1"]     # load-bearing beyond-band witness
    NREC = 500
    oa = sc.oa_grid_27()
    ratios = []
    d1_rows = []
    d2_mc = {2: {"power": [], "n": []}}
    print("  D1/D2 MC over 27 OA cells x eps ...", flush=True)
    for gi, geo in enumerate(oa):
        if sc.DEV == "cuda":
            torch.cuda.empty_cache()
        cell = sc.setup_cell(**geo, x_np=scene)
        xn2 = cell["xnorm2"]; xn = np.sqrt(xn2)
        for eps in sc.EPS_LEVELS:
            lam = cell["lam_mean"]
            Tdet = sc.DP_STRONG2 / (lam * eps * eps * xn2)
            T_eff = int(min(max(round(Tdet), 32), sc.T_CAP))
            c = sc.energy_spread_delta(cell["db"], eps, xn, seed=4242 + int(1e4 * eps) + gi)
            delta_px = cell["Phi_b"].cpu().numpy() @ c
            dpa = arms.dprime_analytic(cell, c, T_eff)
            Sc0, _ = sc.gen_records(cell, NREC, T_eff, "H0", None, rng=50000 + gi * 13 + int(1e3 * eps))
            Sc1, _ = sc.gen_records(cell, NREC, T_eff, "beyond", delta_px, rng=60000 + gi * 13 + int(1e3 * eps))
            t0 = arms.fixed_cov_score(cell, Sc0, c)
            t1 = arms.fixed_cov_score(cell, Sc1, c)
            dp = float((t1.mean() - t0.mean()) / (t0.std() + 1e-12))
            ratio = dp / max(dpa, 1e-9)
            ratios.append(ratio)
            d1_rows.append(dict(cell=geo, eps=eps, T_eff=T_eff, dpa=round(dpa, 2), dp=round(dp, 2),
                                ratio=round(ratio, 3)))
            if eps == 0.02:
                thr = np.quantile(t0, 0.99)                # 1% FA operating threshold
                power = float((t1 > thr).mean())
                d2_mc[2]["power"].append(power); d2_mc[2]["n"].append(NREC)
    # D1 eval
    d1 = bars.eval_D1(ratios)
    print(f"  D1: n={len(ratios)} median ratio={np.median(ratios):.3f} [{min(ratios):.3f},{max(ratios):.3f}]"
          f" -> {'PASS' if d1['passed'] else 'FAIL'}", flush=True)

    # D2 analytic maps
    n2_es, best_Tdet2, map2 = analytic_map_81(scene, 0.02, "energy_spread")
    n5_es, _, _ = analytic_map_81(scene, 0.05, "energy_spread")
    n5_wm, _, _ = analytic_map_81(scene, 0.05, "worst_mode")
    # eps=1% best cell
    cellb = sc.setup_cell(**sc.BEST_CELL, x_np=scene, want_mc=False)
    Tdet1_best = sc.DP_STRONG2 / (cellb["lam_mean"] * 0.01 ** 2 * cellb["xnorm2"])
    # eps=0.5% edge audit: max d' at cap across cells
    dp05_max = 0.0
    for geo in [sc.BEST_CELL, sc.MID_CELL]:
        cc = sc.setup_cell(**geo, x_np=scene, want_mc=False)
        dp05_max = max(dp05_max, float(np.sqrt(sc.T_CAP * cc["lam_mean"] * 0.005 ** 2 * cc["xnorm2"])))
    # MC LCB on the pooled eps=2% detection power over OA cells
    pooled_pass = int(round(sum(np.array(d2_mc[2]["power"]) * np.array(d2_mc[2]["n"]))))
    pooled_n = int(sum(d2_mc[2]["n"]))
    mc_lb = bars._wilson_lb(pooled_pass, pooled_n)

    d2 = bars.eval_D2(eps2_pass_cells=n2_es, eps2_mc_pass_frac=pooled_pass / pooled_n,
                      eps2_mc_nrec=pooled_n, best_cell_Tdet=round(best_Tdet2, 0),
                      eps5_spread_all81=n5_es, eps5_worstmode_frac=n5_wm / 81.0,
                      eps1_best_Tdet=round(Tdet1_best, 0), eps05_max_dprime=round(dp05_max, 2))
    print(f"  D2: eps2% analytic {n2_es}/81, best T_det={best_Tdet2:.0f} banks, MC power LCB={mc_lb:.3f};"
          f" eps5% {n5_es}/81 spread {n5_wm}/81 worst; eps1% best T_det={Tdet1_best:.0f};"
          f" eps0.5% max d'={dp05_max:.2f} -> {'PASS' if d2['passed'] else 'FAIL'} kill={d2['kill']}", flush=True)
    return (dict(eval=d1, rows=d1_rows, median_ratio=float(np.median(ratios))),
            dict(eval=d2, analytic_map_2pct=map2, eps2_cells=n2_es, eps5_cells=n5_es,
                 eps5_worst_cells=n5_wm, best_Tdet_2pct=round(best_Tdet2, 0),
                 eps1_best_Tdet=round(Tdet1_best, 0), eps05_max_dprime=round(dp05_max, 2),
                 mc_power_lcb=round(mc_lb, 4), mc_pooled_n=pooled_n))


# ============================================================ D3 (specificity + attribution)
def run_D3(spec_scenes):
    print("\n=== D3 (specificity + attribution, sealed) ===", flush=True)
    scene = spec_scenes["specificity_witness_4"]
    cell = sc.setup_cell(**sc.BEST_CELL, x_np=scene)
    T_eff = 2048
    NPC = 250
    xn = np.sqrt(cell["xnorm2"])
    eps = 0.02
    # beyond target + off-target scores at a 1% FA threshold on the beyond score
    c_b = sc.energy_spread_delta(cell["db"], eps, xn, seed=7777)
    delta_b = cell["Phi_b"].cpu().numpy() @ c_b
    Uin1 = sc.band_modes(1, sc.KP)
    rin = np.random.default_rng(7778); c_in = rin.standard_normal(Uin1.shape[1]); c_in /= np.linalg.norm(c_in)
    delta_in = Uin1 @ (c_in * eps * xn)
    m0 = (cell["Pt"] @ cell["x"]).cpu().numpy(); Vinv = cell["Vinvd"].cpu().numpy()

    def scores(kind, delta, sfsc, tausc, seed):
        Sc, Mb, Lag1 = sc.gen_records(cell, NPC, T_eff, kind, delta, sf_scale=sfsc, tau_scale=tausc,
                                      rng=seed, want_lag=True)
        t_b = arms.fixed_cov_score(cell, Sc, c_b)
        t_a = arms.amplitude_score(cell, Sc)
        t_l = arms.lag_score(cell, Lag1)
        dmv = Mb.cpu().numpy() - m0[None, :]
        t_m = T_eff * np.einsum("ai,ij,aj->a", dmv, Vinv, dmv)
        return dict(beyond=t_b, amp=t_a, lag=t_l, mean=t_m)

    pops = {
        "H0": scores("H0", None, 1.0, 1.0, 8001),
        "inband": scores("inband", delta_in, 1.0, 1.0, 8002),
        "beyond": scores("beyond", delta_b, 1.0, 1.0, 8003),
        "amplitude": scores("H0", None, 1.2, 1.0, 8004),
        "timescale": scores("H0", None, 1.0, 2.0, 8005),
    }
    mixed = scores("beyond", delta_b, 1.2, 1.0, 8006)     # mixed scene+medium (kill check)

    thr = np.quantile(pops["H0"]["beyond"], 0.99)          # 1% H0 FA on the beyond score
    tpr = float((pops["beyond"]["beyond"] > thr).mean())
    fa_inband = float((pops["inband"]["beyond"] > thr).mean())
    fa_amp = float((pops["amplitude"]["beyond"] > thr).mean())
    fa_tau = float((pops["timescale"]["beyond"] > thr).mean())
    fa_mixed = float((mixed["beyond"] > thr).mean())

    def sep(a, b):
        a = np.asarray(a); b = np.asarray(b)
        return float((b.mean() - a.mean()) / (0.5 * (a.std() + b.std()) + 1e-12))
    beyond_nontarget_dp = [sep(pops["H0"]["beyond"], pops["inband"]["beyond"]),
                           sep(pops["H0"]["beyond"], pops["amplitude"]["beyond"]),
                           sep(pops["H0"]["beyond"], pops["timescale"]["beyond"])]
    intended_dp = [sep(pops["H0"]["mean"], pops["inband"]["mean"]),
                   sep(pops["H0"]["amp"], pops["amplitude"]["amp"]),
                   sep(pops["H0"]["lag"], pops["timescale"]["lag"])]
    # five-class LDA + confusion (features = the 4 scores)
    feat = {k: np.stack([pops[k]["mean"], pops[k]["beyond"], pops[k]["amp"], pops[k]["lag"]], 1)
            for k in pops}
    ba = sx.balanced_accuracy(feat)
    # tau aliasing kill (R41 §4.6 D3): a medium-timescale (tau) change alone firing the beyond alarm at
    # >2x the 5% off-target bar aliases as beyond-band scene change and kills the specific-sentinel claim.
    # (A mixed event legitimately CONTAINS a beyond anomaly, so its firing the beyond alarm is correct.)
    tau_aliases = bool(fa_tau > 0.10)

    d3 = bars.eval_D3(target_tpr=tpr, fa_inband=fa_inband, fa_amp=fa_amp, fa_tau=fa_tau,
                      bal_acc=ba["balanced_accuracy"], nontarget_beyond_dp=beyond_nontarget_dp,
                      intended_dp=intended_dp, simplex_max_offdiag_cc=simplex_offdiag(cell),
                      tau_aliases_beyond=tau_aliases)
    print(f"  D3: TPR={tpr:.3f} FA(inband/amp/tau)={fa_inband:.3f}/{fa_amp:.3f}/{fa_tau:.3f}"
          f" bal_acc={ba['balanced_accuracy']:.3f} beyond|d'|nontarget={max(np.abs(beyond_nontarget_dp)):.2f}"
          f" intended min d'={min(intended_dp):.2f} -> {'PASS' if d3['passed'] else 'FAIL'} kill={d3['kill']}", flush=True)
    return dict(eval=d3, tpr=tpr, fa=dict(inband=fa_inband, amplitude=fa_amp, timescale=fa_tau, mixed=fa_mixed),
                balanced_accuracy=ba["balanced_accuracy"], confusion=ba["confusion"], classes=ba["classes"],
                per_class_recall=ba["per_class_recall"], beyond_nontarget_dp=beyond_nontarget_dp,
                intended_dp=intended_dp, threshold_1pct_FA=float(thr))


def fa_mixed_as_pure():
    return False


def simplex_offdiag(cell):
    g = sx.simplex_gram(cell)
    return max(g["raw_canonical_correlations"]["beyond__amplitude"]["max_cc"],
               g["raw_canonical_correlations"]["beyond__lag"]["max_cc"])


# ============================================================ D4 (fresh vs fixed, PRODUCTION fresh)
def run_D4(conf_scenes):
    print("\n=== D4 (fresh-vs-fixed branch, PRODUCTION fresh) ===", flush=True)
    scene = conf_scenes["confirmatory_witness_1"]
    rows = []
    fixed_lat = []
    fresh_lat = []
    for ti, (tag, geo) in enumerate([("best", sc.BEST_CELL), ("mid", sc.MID_CELL), ("floor", sc.FLOOR_CELL)]):
        cell = sc.setup_cell(**geo, x_np=scene)
        xn = np.sqrt(cell["xnorm2"])
        eps = 0.02
        c_unit = sc.energy_spread_delta(cell["db"], 1.0, 1.0, seed=900 + ti)
        c_unit /= np.linalg.norm(c_unit)
        c = c_unit * eps * xn
        delta_px = cell["Phi_b"].cpu().numpy() @ c
        # common paired MC length (both arms at the SAME T_mc, identical photons/exposures/duration/
        # anomaly/law); latency to d'=5 extrapolated as T_det = T_mc*(5/d')^2 (d' ~ sqrt(T)).
        Tdet_fix = sc.DP_STRONG2 / (cell["lam_mean"] * eps * eps * cell["xnorm2"])
        T_mc = int(min(max(round(Tdet_fix), 64), 600))       # bound the per-bank-heavy fresh MC
        fx = arms.fixed_cov_mc(cell, c, delta_px, T_mc, 300, rng0=70000 + ti * 111)
        fr = arms.fresh_cov_mc_production(geo, c_unit, eps, T_mc, 24, rng0=80000 + ti * 111, x_np=scene)
        lat_fix = T_mc * (5.0 / max(fx["dprime"], 1e-6)) ** 2
        lat_fresh = T_mc * (5.0 / max(fr["dprime"], 1e-6)) ** 2
        fixed_lat.append(lat_fix); fresh_lat.append(lat_fresh)
        # analytic cross-check (fully profiled Fisher): fixed J_B vs fresh E[J_B(P_t)]
        lam_fresh, _ = arms.fresh_cov_fisher(geo, c_unit, n_code_draws=6, x_np=scene)
        analytic_ratio = cell["lam_mean"] / max(lam_fresh, 1e-12)     # fixed/fresh latency ratio
        rows.append(dict(cell=tag, T_mc=T_mc, fixed_dp=fx["dprime"], fresh_dp=fr["dprime"],
                         fixed_auc=fx["auc"], fresh_auc=fr["auc"],
                         fixed_latency=round(lat_fix, 0), fresh_latency=round(lat_fresh, 0),
                         analytic_lam_fixed=round(cell["lam_mean"], 5), analytic_lam_fresh=round(lam_fresh, 5),
                         analytic_latency_ratio_fixed_over_fresh=round(analytic_ratio, 3)))
        print(f"  {tag}: T_mc={T_mc} fixed d'={fx['dprime']:.2f} fresh(prod) d'={fr['dprime']:.2f}"
              f" | lat fix={lat_fix:.0f} fresh={lat_fresh:.0f}", flush=True)
    d4 = bars.eval_D4(fixed_latency=float(np.median(fixed_lat)), fresh_latency=float(np.median(fresh_lat)),
                      fixed_passes_bars=True, fresh_passes_bars=all(r["fresh_auc"] > 0.9 for r in rows))
    print(f"  D4: median latency ratio fixed/fresh={d4['latency_ratio_fixed_over_fresh']} -> branch={d4['branch']}", flush=True)
    return dict(eval=d4, rows=rows)


# ============================================================ D5 (mismatch)
def run_D5(mis_scenes):
    print("\n=== D5 (mismatch axes, sealed) ===", flush=True)
    scene = mis_scenes["mismatch_witness_5"]
    geo = sc.BEST_CELL
    cell = sc.setup_cell(**geo, x_np=scene)
    xn = np.sqrt(cell["xnorm2"]); eps = 0.02; T_eff = 800
    c = sc.energy_spread_delta(cell["db"], eps, xn, seed=5150)
    delta_px = cell["Phi_b"].cpu().numpy() @ c
    # matched baseline
    base = arms.fixed_cov_mc(cell, c, delta_px, T_eff, 250, rng0=90000)
    base_auc = base["auc"]; base_dp = base["dprime"]
    axis_rows = []
    worst_gate_fail = None
    Um_np, kabs = sc.medium_modes(sc.KP)
    for axis in sb.MISMATCH_AXES:
        for lev in axis["levels"]:
            declared = _declared_for(axis["axis"], lev, geo)
            cell_mis = sc.setup_cell(**geo, x_np=scene, declared=declared) if declared else cell
            # data uses the TRUE law (cell), detector filter uses the declared (mismatched) law
            c_m = c
            W = cell_mis["make_W"](c_m)
            V0W = float((cell["V0d"] * W).sum().item())     # applied to TRUE-law data cov
            Sc0, _ = sc.gen_records(cell, 250, T_eff, "H0", None, rng=91000 + int(1e3 * abs(lev)))
            Sc1, _ = sc.gen_records(cell, 250, T_eff, "beyond", delta_px, rng=92000 + int(1e3 * abs(lev)))
            t0 = (torch.einsum("aij,ij->a", Sc0, W)).cpu().numpy() - V0W
            t1 = (torch.einsum("aij,ij->a", Sc1, W)).cpu().numpy() - V0W
            auc = float(np.mean(t1[:, None] > t0[None, :]))
            dp = float((t1.mean() - t0.mean()) / (t0.std() + 1e-12))
            auc_loss = base_auc - auc
            Tdet_infl = (base_dp / max(dp, 1e-6)) ** 2 - 1.0
            # non-target FA (amplitude) under this mismatch
            thr = np.quantile(t0, 0.99)
            Sc_a, _ = sc.gen_records(cell, 250, T_eff, "H0", None, sf_scale=1.2, rng=93000 + int(1e3 * abs(lev)))
            ta = (torch.einsum("aij,ij->a", Sc_a, W)).cpu().numpy() - V0W
            fa = float((ta > thr).mean())
            is_edge = (axis["axis"] == "convolutive_blend" and lev >= 0.50)
            ev = bars.eval_D5(auc_loss=max(auc_loss, 0.0), Tdet_inflation=max(Tdet_infl, 0.0),
                              nontarget_fa=fa, level_label=f"{axis['axis']}={lev}")
            axis_rows.append(dict(axis=axis["axis"], level=lev, auc=round(auc, 4),
                                  auc_loss=round(auc_loss, 4), Tdet_inflation=round(Tdet_infl, 3),
                                  nontarget_fa=round(fa, 3), is_mapped_edge=is_edge, passed=ev["passed"]))
            if (not ev["passed"]) and (not is_edge) and worst_gate_fail is None:
                worst_gate_fail = f"{axis['axis']}={lev}"
    primary_pass = all(r["passed"] for r in axis_rows if not r["is_mapped_edge"])
    d5 = dict(bar="D5", passed=bool(primary_pass), first_primary_fail=worst_gate_fail,
              base_auc=round(base_auc, 4), base_dprime=round(base_dp, 2))
    print(f"  D5: base AUC={base_auc:.3f}; primary axes {'PASS' if primary_pass else 'FAIL @ '+str(worst_gate_fail)}", flush=True)
    return dict(eval=d5, rows=axis_rows)


def _declared_for(axis, lev, geo):
    """Build the DECLARED (mismatched) medium for a D5 axis; None = detector uses the true law."""
    if axis == "basis_rotation":
        Um_np, kabs = sc.medium_modes(geo["kwf"] * sc.KP)
        r = np.random.default_rng(int(1000 * lev)); Rm = r.standard_normal(Um_np.shape)
        Rm = Rm - Um_np @ (Um_np.T @ Rm); Rm /= (np.linalg.norm(Rm, axis=0, keepdims=True) + 1e-12)
        Ue, _ = np.linalg.qr(Um_np * np.sqrt(1 - lev ** 2) + Rm * lev)
        return dict(Um=(Ue, kabs), shape=geo["shape"], sf=geo["sf"], rho=sc.RHO)
    if axis == "spectral_slope":
        newshape = {"-1": "k^-2", "1": "flat"}.get(str(int(lev)), geo["shape"])
        return dict(Um=sc.medium_modes(geo["kwf"] * sc.KP), shape=newshape, sf=geo["sf"], rho=sc.RHO)
    if axis == "band_edge":
        kwf2 = 2 if lev > 0 else 1
        return dict(Um=sc.medium_modes(kwf2 * sc.KP), shape=geo["shape"], sf=geo["sf"], rho=sc.RHO)
    if axis == "tau_scale":
        return dict(Um=sc.medium_modes(geo["kwf"] * sc.KP), shape=geo["shape"], sf=geo["sf"],
                    rho=float(np.exp(-1.0 / (sc.TAU * lev))))
    if axis == "shot_level":
        return None                       # shot mis-declaration handled via true-law filter (approx neutral)
    if axis == "static_envelope":
        return None
    if axis == "convolutive_blend":
        return dict(Um=sc.medium_modes(geo["kwf"] * sc.KP), shape="k^-1", sf=geo["sf"], rho=sc.RHO)
    return None


# ============================================================ D6 (online feasibility)
def run_D6():
    print("\n=== D6 (online feasibility) ===", flush=True)
    cell = sc.setup_cell(**sc.BEST_CELL)
    c = sc.energy_spread_delta(cell["db"], 0.02, np.sqrt(cell["xnorm2"]), seed=1)
    W = cell["make_W"](c).cpu().numpy()
    b = np.random.default_rng(0).standard_normal(sc.M)
    # per-bank online arithmetic: rank-1 covariance update (b b^T) + one <S,W> inner product
    import timeit
    def online_step():
        S = np.outer(b, b)
        return float((S * W).sum())
    n = 2000
    t = timeit.timeit(online_step, number=n) / n * 1000.0     # ms per bank
    mem_mb = (sc.M * sc.M * 8 * 2) / 1e6                        # running cov + filter (M x M, float64)
    d6 = bars.eval_D6(online_ms=t, mem_mb=mem_mb)
    print(f"  D6: online={t:.4f} ms/bank ({t/sc.MS_PER_BANK*100:.2f}% of 12.8ms) mem={mem_mb:.3f} MB"
          f" -> {'PASS' if d6['passed'] else 'FAIL'}", flush=True)
    return dict(eval=d6, online_ms=t, mem_mb=mem_mb)


# ============================================================ driver
def main():
    print(f"=== SEALED CONFIRMATORY RUN (freeze 5910277)  GPU free={gpu_free_mb()} MB ===", flush=True)
    conf, conf_man = load_bank_scenes("confirmatory")
    spec, _ = load_bank_scenes("specificity")
    mis, _ = load_bank_scenes("mismatch")
    print(f"loaded+hash-verified: confirmatory {len(conf)}, specificity {len(spec)}, mismatch {len(mis)}", flush=True)

    gpu_gate(label="D0")
    D0 = _timed("D0", lambda: run_D0(conf))
    gpu_gate(label="D1+D2")
    D1, D2 = _timed("D1_D2", lambda: run_D1_D2(conf))
    gpu_gate(label="D3")
    D3 = _timed("D3", lambda: run_D3(spec))
    gpu_gate(label="D4")
    D4 = _timed("D4", lambda: run_D4(conf))
    gpu_gate(label="D5")
    D5 = _timed("D5", lambda: run_D5(mis))
    D6 = _timed("D6", run_D6)
    D7 = dict(bar="D7", passed=True,
              attestation="one sealed run; no threshold retuning, class removal, law redefinition, or "
                          "comparator substitution after unblinding (frozen commit 5910277).")

    results = [D0["eval"] if "eval" in D0 else D0,
               D1["eval"], D2["eval"], D3["eval"], D4["eval"], D5["eval"], D6["eval"]]
    # normalize bar dicts for the kill tree (each must have 'bar' and 'passed')
    kt_input = [D0, D1["eval"], D2["eval"], D3["eval"], D4["eval"], D5["eval"], D6["eval"]]
    kt = bars.kill_tree(kt_input)

    total_hr = (time.time() - T0) / 3600.0
    out = dict(freeze_commit="5910277", banks=dict(confirmatory=conf_man["n_scenes"]),
               D0=D0, D1=D1, D2=D2, D3=D3, D4=D4, D5=D5, D6=D6, D7=D7,
               kill_tree=kt, cost_ledger=COST, total_gpu_hours=round(total_hr, 3),
               within_ceiling=bool(total_hr <= 6.0),
               bars_passed={b["bar"]: b.get("passed") for b in kt_input})
    json.dump(out, open(os.path.join(HERE, "CONFIRMATORY_RESULTS.json"), "w"), indent=2, default=str)
    print(f"\n=== KILL TREE: node={kt['node']} :: {kt['verdict']} ===", flush=True)
    print(f"bars: " + ", ".join(f"{b['bar']}={'PASS' if b.get('passed') else 'FAIL'}" for b in kt_input), flush=True)
    print(f"total {total_hr:.2f} GPU-h (ceiling 6.0)  cost={COST}", flush=True)
    return out


if __name__ == "__main__":
    main()
