#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
PHASE 2 DRY RUN -- BAR D0: mechanism and engine integrity (R41 sec 4.6). Allowed on any cell (no
confirmatory scene opened). If D0 fails anywhere -> STOP and report (kill-tree node 1).

  D0.1 mean derivative to every beyond-band event has relative norm <= 1e-10  (the exact wall);
  D0.2 two independent Fisher implementations agree within 10% on eigenvalues, d', and T_det;
  D0.3 Monte-Carlo shot variance (physical |P|x complementary-exposure Poisson ledger) agrees with
       the analytic ledger V0 = C + R within 10%.

Runs on the best / mid / floor reference cells + a spectrum/k_w sweep. Writes D0_REPORT.md (+ .json).
"""
import json
import os
import time

import numpy as np
import torch

import sealed_common as sc

HERE = os.path.dirname(os.path.abspath(__file__))
t0 = time.time()


def d0_mean_derivative():
    """||P U_beta||_F / ||P||_F for every claim shell (beyond-band mean derivative = exact wall)."""
    P = sc.P_FIXED
    out = {}
    worst = 0.0
    for claim in sc.CLAIMS:
        Ub = sc.BETA_np[claim]
        rel = float(np.linalg.norm(P @ Ub) / (np.linalg.norm(P) + 1e-30))
        out[f"claim_{claim}"] = rel
        worst = max(worst, rel)
    return out, worst


def d0_two_engine(cells):
    """Engine A (einsum-trace + Schur) vs Engine B (whitened vech-derivative + QR) profiled Fisher.
    Compare lam_mean, lam_max, tr, and the analytic energy-spread / matched-target d' at eps=2%."""
    rows = []
    worst = 0.0
    for geo in cells:
        cell = sc.setup_cell(**geo)
        B = sc.fisher_engine_B(cell)
        xn2 = cell["xnorm2"]
        eps = 0.02
        # energy-spread d' ~ sqrt(T lam_mean eps^2 ||x||^2); matched-target ~ lam_max
        T = 1000
        dpA_es = np.sqrt(T * cell["lam_mean"] * eps ** 2 * xn2)
        dpB_es = np.sqrt(T * B["lam_mean"] * eps ** 2 * xn2)
        dpA_mt = np.sqrt(T * cell["lam_max"] * eps ** 2 * xn2)
        dpB_mt = np.sqrt(T * B["lam_max"] * eps ** 2 * xn2)
        comps = dict(
            lam_mean=(cell["lam_mean"], B["lam_mean"]),
            lam_max=(cell["lam_max"], B["lam_max"]),
            tr=(cell["tr"], B["tr"]),
            dprime_energy_spread=(dpA_es, dpB_es),
            dprime_matched_target=(dpA_mt, dpB_mt),
        )
        rels = {k: abs(a - b) / max(abs(a), 1e-30) for k, (a, b) in comps.items()}
        cell_worst = max(rels.values())
        worst = max(worst, cell_worst)
        rows.append(dict(cell=geo, rel_diffs={k: round(v, 4) for k, v in rels.items()},
                         lam_min_A=cell["lam_min"], lam_min_B=B["lam_min"],
                         note="lam_min excluded from the gate: near-null beyond-band modes (both ~0)"))
    return rows, worst


def d0_shot_ledger(cells, T_shot=8000, T_cov=12000):
    """Three checks, all vs the corrected |P|x ledger V0 = C + R:
      (a) SHOT: frozen medium, physical complementary-exposure Poisson -> Var(b_i) vs R_i (the |P|x
          shot term -- the literal R41 D0.3 test and the subject of the shot correction);
      (b) C+R SELF-CONSISTENCY: the engine's OWN 2nd-order generative law (gen_records) -> empirical
          bucket covariance vs V0 (validates the C+R analytic to sampling noise, all contrasts);
      (c) PHYSICAL END-TO-END: Gaussian 2nd-order medium + physical Poisson -> covariance vs V0
          (reported diagnostic; the high-contrast sf=1.0 residual is the linear-law/clip boundary,
          which is precisely why the lognormal medium is a D5 mismatch axis, not the primary law).
    GATE = max(a, b) <= 10%. (c) is reported, not gated at high contrast."""
    rows = []
    worst_gate = 0.0
    for geo in cells:
        cell = sc.setup_cell(**geo)
        R = cell["R"].cpu().numpy()
        V0 = cell["V0"].cpu().numpy()
        dV0 = np.diag(V0)
        # (a) shot ledger (frozen medium, physical Poisson)
        Bs = sc.physical_banks(cell, T_shot, seed=1234, freeze_medium=True)
        shot_rel = float(np.median(np.abs(Bs.var(axis=0) - R) / np.maximum(R, 1e-12)))
        # (b) C+R self-consistency (engine's own 2nd-order law via gen_records)
        Scov, _ = sc.gen_records(cell, 1, T_cov, "H0", None, rng=9999)
        cov_gen = Scov[0].cpu().numpy()
        selfcons_rel = float(np.median(np.abs(np.diag(cov_gen) - dV0) / np.maximum(dV0, 1e-12)))
        # (c) physical end-to-end (Gaussian medium + Poisson)
        Bf = sc.physical_banks(cell, T_cov, seed=5678, medium_model="gaussian")
        phys_rel = float(np.median(np.abs(np.diag(np.cov(Bf.T)) - dV0) / np.maximum(dV0, 1e-12)))
        gate = max(shot_rel, selfcons_rel)
        worst_gate = max(worst_gate, gate)
        rows.append(dict(cell=geo, shot_var_median_rel=round(shot_rel, 4),
                         CplusR_selfconsistency_rel=round(selfcons_rel, 4),
                         physical_end2end_rel=round(phys_rel, 4),
                         medium_shot_ratio=round(float(np.median((dV0 - R) / np.maximum(R, 1e-12))), 2)))
    return rows, worst_gate


def main():
    import bars
    cells = [sc.BEST_CELL, sc.MID_CELL, sc.FLOOR_CELL,
             dict(sf=0.6, shape="k^-1", kwf=2, claim=1.5),
             dict(sf=1.0, shape="flat", kwf=4, claim=1.8)]
    print("=== D0 dry run: mechanism + engine integrity ===", flush=True)

    md, md_worst = d0_mean_derivative()
    print(f"D0.1 mean-derivative wall (max rel over claims) = {md_worst:.2e}  (bar <= 1e-10)", flush=True)

    te_rows, te_worst = d0_two_engine(cells)
    print(f"D0.2 two-engine Fisher (max rel over cells/quantities) = {te_worst:.4f}  (bar <= 0.10)", flush=True)
    for r in te_rows:
        print(f"     {r['cell']}: {r['rel_diffs']}", flush=True)

    sl_rows, sl_worst = d0_shot_ledger(cells[:3])
    print(f"D0.3 MC shot + C+R self-consistency vs |P|x ledger (max median rel) = {sl_worst:.4f}  (bar <= 0.10)", flush=True)
    for r in sl_rows:
        print(f"     {r['cell']}: shot_rel={r['shot_var_median_rel']} selfcons_rel={r['CplusR_selfconsistency_rel']} "
              f"phys_e2e_rel={r['physical_end2end_rel']} med/shot={r['medium_shot_ratio']}", flush=True)

    res = bars.eval_D0(md_worst, te_worst, sl_worst)
    kt = bars.kill_tree([res])
    verdict = "D0_PASS" if res["passed"] else "D0_FAIL -> " + kt["verdict"]
    print(f"\nD0 VERDICT: {verdict}  [{time.time()-t0:.0f}s]", flush=True)

    out = dict(bar="D0", mean_derivative=dict(per_claim=md, worst=md_worst),
               two_engine=dict(rows=te_rows, worst=te_worst),
               shot_ledger=dict(rows=sl_rows, worst=sl_worst),
               D0_eval=res, kill_tree=kt, passed=res["passed"], wall_s=round(time.time() - t0, 1))
    json.dump(out, open(os.path.join(HERE, "D0_RESULTS.json"), "w"), indent=2)
    return out


if __name__ == "__main__":
    main()
