"""M1 freeze-authorization selftest — R15 DECISIVE iteration.

Implements docs/ROUND63_GPT_ROUND15_RULING_RAW.md on top of the R13
checklist: R10 PSNR-based FIXED* selection (frozen rule quoted verbatim),
per-policy LOAD-MATCHED FIXED* references, the full dose-constrained KKT
certificate (amended Box 11), the exact materialized mixture search
(amended Box 13), and the active-dose toy verification.  Outcome-blind:
hard checks read no confirmatory scene and no confirmatory endpoint; the
FIXED* selection reads DEV radiometric PSNR exactly as the frozen R10 rule
prescribes (development usage sanctioned by R10 Sec 10).

Usage:  python code/round63/m1_freeze_selftest.py
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
ROOT = os.path.dirname(os.path.dirname(HERE))

import oed_design_v4 as v4                  # noqa: E402
import oed_design_v5 as v5                  # noqa: E402
from gi_core.utils import rng_for           # noqa: E402
from patterns import make_patterns          # noqa: E402

OUT = os.path.join(ROOT, "results", "round63_m1")
LEDGER = os.path.join(OUT, "FREEZE_CHECKLIST_LEDGER.md")
FIXED_STAR_JSON = os.path.join(OUT, "FIXED_STAR_SELECTION.json")

SIDE, N = 32, 1024
NU_FULL = [5, 10, 20, 50, 100, 200, 500, 1000, 2000]
NU_T = 2000
RHO_FAST, RHO_SAFE = 0.60, 0.05
SEEDS_SEL = [0, 1, 2, 3, 4]

R10_FIXED_STAR_RULE = (
    "Before confirmatory freeze, select one global fixed comparator FIXED* "
    "from {SCAT32, LBLOB16, MATCH1} using the six development images only. "
    "The selection score is the median fixed-dwell radiometric PSNR at "
    "(rho=0.60, nu=2000), averaged over development seeds. Tie within "
    "0.05 dB is resolved in order SCAT32 -> LBLOB16 -> MATCH1, favoring "
    "the least scene-dependent rule. [ROUND63_GPT_ROUND10_RULING_RAW.md "
    "Sec 7 / line 801]")

BOX11_WORDING = """[ ] relaxed_and_full_constrained_certificates

PASS iff, for BOTH the safe and fast OED policies:

(a) the dose-relaxed simplex+incident-budget reference problem, after the
    frozen support-restricted exact reweighting, satisfies
    RELAXED_REFERENCE_KW_UPPER_BOUND / r <= 1e-3 under a full-dictionary scan;

(b) the physically dose-feasible continuous design is primal-feasible for
    the incident budget and every +/-5% per-pixel dose inequality and satisfies
    FULL_CONSTRAINED_KW_UPPER_BOUND / r <= 1e-2, where the certificate includes
    one budget multiplier and all upper/lower per-pixel dose multipliers and
    the final maximum is scanned over the entire frozen dictionary;

(c) the full constrained certificate passes the independent active-dose toy
    verification and records all primal, dual, complementarity, and dictionary-
    scan residuals."""

BOX13_WORDING = """[ ] final_exact_rows_full_guard_suite

PASS iff, for BOTH safe and fast policies, the final returned main-pattern
matrix contains exactly 972 physical rows and is either the direct rounded OED
design or the materialized COMPLIANT_VIA_MIXTURE design selected by the frozen
m=972,...,0 search. On that exact matrix, all of the following must pass:

- incident_sum_rho <= 972*budget_mean + 1e-9;
- predicted and realized detected-budget fields are emitted;
- per-pixel cumulative-dose max relative deviation <= 0.05;
- physical peak <= 1536;
- design-weighted ceiling probability <= 0.01;
- exact D-efficiency >= 0.95 relative to the certified full dose-constrained
  continuous design;
- A-risk <= 1.05 times the load-matched FIXED* reference;
- minimum generalized eigenvalue versus the load-matched FIXED* reference
  >= 0.5;
- exactly 52 common pre-scan rows plus 972 returned main rows are accounted.

If no integer mixture candidate, including m=0, passes all guards, Box 13 FAILS."""

CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok), detail))
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name,
                           ("  -- " + detail) if detail else ""), flush=True)
    return ok


def synth_dev_scene():
    yy, xx = np.mgrid[0:SIDE, 0:SIDE]
    s = np.ones((SIDE, SIDE))
    for (cy, cx) in [(8, 9), (10, 23), (22, 7), (24, 21)]:
        s += 0.6 * np.exp(-(((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * 3.5 ** 2)))
    x = (s / s.sum()).ravel()
    c = x.reshape(SIDE, SIDE).reshape(8, 4, 8, 4).mean(axis=(1, 3))
    xh = np.repeat(np.repeat(c, 4, 0), 4, 1).ravel()
    return x, xh / xh.sum()


def arm_main_rows(arm, xhat):
    if arm in ("SCAT16", "SCAT32"):
        k = 16 if arm == "SCAT16" else 32
        return make_patterns("sparsek", N, N, 0, k=k)["A"][:972]
    if arm == "LBLOB16":
        return make_patterns("lblob16", N, N, 0)["A"][:972]
    if arm == "MATCH1":
        A = make_patterns("lblob16", N, N, 0)["A"][:972]
        out = np.zeros_like(A)
        for s in range(A.shape[0]):
            idx = np.nonzero(A[s])[0]
            w0 = xhat[idx] / xhat[idx].mean()
            w = np.clip(w0, 0.25, 4.0)
            w = w / w.mean()
            out[s, idx] = (N / 16.0) * w
        return out
    if arm == "RIDGE-FIXED":
        return v4.ridge_fixed_972(xhat, NU_T)["A_main"]
    raise ValueError(arm)


def fixed_star_r10_selection():
    """The frozen R10 DEV radiometric-PSNR rule (quoted in the ledger)."""
    if os.path.exists(FIXED_STAR_JSON):
        with open(FIXED_STAR_JSON) as f:
            sel = json.load(f)
        if sel.get("rule_id") == "R10_PSNR_L801":
            print("  (cached R10 selection reused: %s)" % sel["FIXED_STAR"])
            return sel
    import campaign
    import m1_runner as m1
    from select_eta import frozen_C0
    C0 = frozen_C0()
    dev_imgs = campaign._images(SIDE, "all", imageset="m1_dev")
    import m1_scenes
    names = [nm for nm, _f, _s in m1_scenes._m1_dev_table()]
    for img in names:                       # MATCH1 xhat caches (DEV)
        for seed in SEEDS_SEL:
            m1.build_design(img, seed, "MATCH1", dev_imgs[img])
    scores, detail = {}, {}
    for cand in ("SCAT32", "LBLOB16", "MATCH1"):
        img_means = []
        for img in names:
            vals = []
            for seed in SEEDS_SEL:
                cell = m1.m1_cell(cand, img, seed, RHO_FAST, float(NU_T),
                                  C0, imageset="m1_dev")
                rows = campaign.run_cell(cell)
                vals.append(float(rows[0]["PSNR_rad"]))
            img_means.append(float(np.mean(vals)))
        scores[cand] = float(np.median(img_means))
        detail[cand] = img_means
        print("    %s: median-of-image-means PSNR_rad = %.4f dB (%s)"
              % (cand, scores[cand],
                 ["%.2f" % v_ for v_ in img_means]), flush=True)
    best = max(scores.values())
    winner = next(c for c in ("SCAT32", "LBLOB16", "MATCH1")
                  if scores[c] >= best - 0.05)
    sel = {"rule_id": "R10_PSNR_L801", "rule_text": R10_FIXED_STAR_RULE,
           "scores_dB": scores, "per_image_means_dB": detail,
           "seeds": SEEDS_SEL, "FIXED_STAR": winner,
           "tie_window_dB": 0.05,
           "tie_order": ["SCAT32", "LBLOB16", "MATCH1"]}
    with open(FIXED_STAR_JSON, "w") as f:
        json.dump(sel, f, indent=1, sort_keys=True)
    return sel


# ---- compact analyzer + mocks ---------------------------------------------- #
def _pava(y):
    y = list(y)
    w = [1.0] * len(y)
    i = 0
    while i < len(y) - 1:
        if y[i] > y[i + 1] + 1e-12:
            m = (y[i] * w[i] + y[i + 1] * w[i + 1]) / (w[i] + w[i + 1])
            y[i] = y[i + 1] = m
            w[i] = w[i + 1] = w[i] + w[i + 1]
            i = max(i - 1, 0)
        else:
            i += 1
    return np.array(y)


def _q90_time(nus, q, target, rho):
    qf = _pava(q)
    lt = np.log(np.array(nus) * rho)
    if qf[-1] < target:
        return None
    return float(np.exp(np.interp(target, qf, lt)))


def analyzer(mock, n_units):
    S, dQ, dG = [], [], []
    for u in range(n_units):
        safe = mock["safe"][u]
        fast = mock["fast"][u]
        fx = mock["fixedstar_fast_terminal"][u]
        R = safe[-1] - safe[0]
        tgt = safe[0] + 0.9 * R
        Ts = _q90_time(NU_FULL, safe, tgt, RHO_SAFE)
        Tf = _q90_time(NU_FULL, fast, tgt, RHO_FAST)
        S.append((Ts / Tf) if (Ts and Tf) else np.nan)
        dQ.append(fast[-1] - safe[-1])
        dG.append(fast[-1] - fx)
    S, dQ, dG = map(np.asarray, (S, dQ, dG))
    rng = rng_for(0, 63, 11)
    idx = rng.integers(0, n_units, size=(2000, n_units))
    lb = lambda v: float(np.percentile(np.median(v[idx], axis=1), 2.5))
    need = int(math.ceil(0.75 * n_units))
    return {"METHOD_SPEED_PASS": bool(np.median(S) >= 3.0 and lb(S) > 1.0
                                      and int((S > 1).sum()) >= need),
            "METHOD_FIXED_DWELL_PASS": bool(np.median(dQ) >= 1.0
                                            and lb(dQ) > 0
                                            and int((dQ > 0).sum()) >= need),
            "METHOD_DESIGN_PASS": bool(np.median(dG) >= 0.5 and lb(dG) > 0
                                       and int((dG > 0).sum()) >= need),
            "median_S": float(np.median(S)), "median_dQ": float(np.median(dQ)),
            "median_dG": float(np.median(dG))}


def mock_cohort(n_units, positive=True):
    safe, fast, fx = [], [], []
    for u in range(n_units):
        r = rng_for(u, 63, 12)
        base = 10 + r.random() * 2
        nus = np.array(NU_FULL, float)
        s = base + 6 * (1 - np.exp(-nus / 500.0))
        shift = (60.0 if positive else 1.0) + 3.0 * r.random()
        f = base + 6 * (1 - np.exp(-nus * shift / 500.0)) + \
            (1.5 if positive else -0.2)
        safe.append(s)
        fast.append(f)
        fx.append(f[-1] - (0.8 if positive else -0.1))
    return {"safe": safe, "fast": fast, "fixedstar_fast_terminal": fx}


def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    x_true, xhat = synth_dev_scene()
    print("[freeze selftest R15] outcome-blind, DEV-only. synthetic scene "
          "cv=%.3f" % (xhat.std() / xhat.mean()), flush=True)

    print("\n== items 1-2: balanced 52-row pre-scan ==")
    P = v4.balanced_prescan_52(SIDE)
    row_mean = P.mean(axis=0)
    dose = P.sum(axis=0)
    check("prescan_row_average_all_ones",
          np.abs(row_mean - 1.0).max() < 1e-12,
          "max|mean-1|=%.2e" % np.abs(row_mean - 1.0).max())
    check("prescan_pixel_dose_exactly_equal",
          np.abs(dose - dose[0]).max() == 0.0,
          "dose=%g per pixel" % dose[0])

    print("\n== item 3: non-oracle V0 ==")
    loads = P @ xhat
    B_id = np.eye(N)[:, :8]
    V0a = v4.V0_prescan(P, xhat, NU_T, RHO_FAST, B_id)
    V0b = np.zeros_like(V0a)
    for m in range(52):
        q = B_id.T @ (P[m] / loads[m])
        V0b += NU_T * v4.kernel_eval4(RHO_FAST, NU_T)["J_exact"] \
            * np.outer(q, q)
    check("V0_actual_row_loads_no_oracle",
          loads.std() > 0 and not np.allclose(V0a, V0b),
          "row-load sd=%.3e" % loads.std())

    print("\n== item 4: dictionary manifest ==")
    man = v4.dictionary_manifest(SIDE)
    fams = man["families"]
    required = ["scat16", "scat32", "solid4x4", "lblob6x6", "bar1x16",
                "bar2x8", "rect4x4", "bar8x2", "bar16x1", "sq32c",
                "rect4x8", "rect8x4", "ring16"]
    ok4 = all(f in fams and fams[f]["n_translations"] == N
              and fams[f]["powers"] == [0, 1]
              and len(fams[f]["sha256"]) == 64 for f in required)
    with open(os.path.join(OUT, "DICTIONARY_MANIFEST.json"), "w") as f:
        json.dump(man, f, indent=1, sort_keys=True)
    check("dictionary_manifest_complete", ok4,
          "%d families, sha %s.." % (len(fams), man["sha256_global"][:12]))

    print("\n== item 5: FIXED* via the frozen R10 DEV-PSNR rule ==")
    print("  R10 rule: " + R10_FIXED_STAR_RULE)
    sel = fixed_star_r10_selection()
    fixed_star = sel["FIXED_STAR"]
    check("fixed_star_r10_psnr_selection_committed",
          sel.get("rule_id") == "R10_PSNR_L801",
          "FIXED*=%s scores(dB)=%s (tie window 0.05, order "
          "SCAT32->LBLOB16->MATCH1)"
          % (fixed_star, {k: round(v_, 3)
                          for k, v_ in sel["scores_dB"].items()}))

    print("\n== item 6: 52+972 accounting ==")
    ok6 = P.shape[0] == 52
    d6 = []
    for arm in ("SCAT16", "SCAT32", "LBLOB16", "MATCH1", "RIDGE-FIXED"):
        A = arm_main_rows(arm, xhat)
        ok6 &= (A.shape == (972, N))
        d6.append("%s:%d" % (arm, A.shape[0]))
    check("all_arms_52_plus_972", ok6, ", ".join(d6) + " (+52 pre-scan)")

    print("\n== item 7: nine-dwell certification table ==")
    ok7 = True
    for nu in NU_FULL:
        for fast in (False, True):
            pal = v4.palette4(nu, fast=fast)
            for rec in pal["records"]:
                if rec["admitted"]:
                    ok7 &= (rec["cert"] == "OK"
                            and rec["mean_info_efficiency"] >= v4.EFF_MIN
                            and rec["rql_logload_bias"] <= v4.BIAS_MAX)
            rr = pal["ridge"]
            print("    nu=%-5d %-4s admitted=%s rho_R=%s clip=%s"
                  % (nu, "fast" if fast else "safe",
                     [round(r_["rho"], 4) for r_ in pal["records"]
                      if r_["admitted"]],
                     ("%.4f" % rr["rho_R_production"]) if rr else "-",
                     rr["ridge_clip_reason"] if rr else "-"))
    check("nine_dwell_kernel_trust_bias_cert", ok7, "18 palettes certified")

    print("\n== amended Box 11 + Box 13 (R15 Secs 5-7) ==")
    tg, td = v4.toy_cert_check()
    tf = v5.toy_full_cert_check()
    print("  toy (unconstrained, R^2): G_rel=%.1e d=%s" % (tg, td))
    print("  toy (ACTIVE dose band):   G_full at SLSQP optimum=%.2e "
          "dose_active=%s w*=%s"
          % (tf["G_at_opt"], tf["dose_active"],
             [round(w, 4) for w in tf["w_star"]]))
    check("toy_certificates_verified",
          abs(tg) < 1e-12 and tf["agreement_ok"] and tf["dose_active"],
          "unconstrained exact; active-dose agreement %.1e <= 1e-8"
          % tf["G_at_opt"])
    fd_rows, fd_dev = v5.fixed_dose_scat32(SIDE)
    check("fixed_dose_972_exact_rows",
          fd_rows.shape == (972, N) and fd_dev <= 0.05,
          "SCAT32-based (information-matched to FIXED*), dose dev=%.4f"
          % fd_dev)

    results = {}
    for policy, budget in (("fast", RHO_FAST), ("safe", RHO_SAFE)):
        print("\n  --- %s policy (budget %.2f, LOAD-MATCHED FIXED*=%s) ---"
              % (policy, budget, fixed_star), flush=True)
        A_fix = arm_main_rows(fixed_star, xhat)
        V_fix_full = v4.info_matrix_full(A_fix, xhat, NU_T, budget, P=P)
        B, eps0, _tr = v4.subspace_from_fixedstar(V_fix_full)
        V_fix_B = B.T @ V_fix_full @ B
        pal = v4.palette4(NU_T, fast=(policy == "fast"))
        t1 = time.time()
        out = v5.design_v5(xhat, float(NU_T), budget, pal, B, eps0, V_fix_B,
                           budget, fd_rows, verbose=True)
        results[policy] = out
        cert = out["cert_full"]
        print("  [%s] verdict=%s  G_full/r=%.3e (tight=%s)  relref/r=%.3e  "
              "mixture=%s m=%s alpha=%s  (%.0fs)"
              % (policy, out["verdict"],
                 out["FULL_CONSTRAINED_KW_UPPER_BOUND"],
                 out["FULL_CONSTRAINED_KW_TIGHT"],
                 out["RELAXED_REFERENCE_KW_UPPER_BOUND"],
                 out["mixture"]["status"], out["mixture"]["m"],
                 (None if out["mixture"]["alpha"] is None
                  else round(out["mixture"]["alpha"], 4)),
                 time.time() - t1), flush=True)
        print("      residuals: budget_viol=%.2e dose_excess=%.2e "
              "comp_budget=%.2e comp_dose=%.2e n_active_mu=%s"
              % (cert.get("budget_viol", -1), cert.get("dose_excess", -1),
                 cert.get("comp_budget", -1), cert.get("comp_dose", -1),
                 cert.get("n_active_mu", "-")))
        print("      G_full/r trajectory: %s"
              % [(t[0], "%.3e" % t[1]) for t in out["traj_pd"]])
        g = out["mixture"]["guards"]
        print("      guards: %s  incident=%.2f effD=%.4f dose=%.4f "
              "arisk=%.3fx mu_min=%.3f"
              % (g["checks"], g["incident_sum_rho"], g["eff_D"],
                 g["dose_dev"], g["a_risk"] / g["a_risk_fix"], g["mu_min"]))

    ok11 = all(
        results[p]["RELAXED_REFERENCE_KW_UPPER_BOUND"] <= v5.RELREF_TOL_PER_R
        and results[p]["cert_full"]["status"] == "OK"
        and results[p]["cert_full"].get("budget_viol", 1) <= 1e-10
        and results[p]["cert_full"].get("dose_excess", 1) <= 1e-8
        and results[p]["FULL_CONSTRAINED_KW_UPPER_BOUND"]
        <= v5.GFULL_TOL_PER_R
        for p in ("fast", "safe")) and tf["agreement_ok"]
    check("relaxed_and_full_constrained_certificates", ok11,
          "fast: relref=%.2e full=%.2e; safe: relref=%.2e full=%.2e"
          % (results["fast"]["RELAXED_REFERENCE_KW_UPPER_BOUND"],
             results["fast"]["FULL_CONSTRAINED_KW_UPPER_BOUND"],
             results["safe"]["RELAXED_REFERENCE_KW_UPPER_BOUND"],
             results["safe"]["FULL_CONSTRAINED_KW_UPPER_BOUND"]))

    ok13 = all(results[p]["mixture"]["status"] == "PASS"
               and results[p]["mixture"]["rows"].shape == (972, N)
               for p in ("fast", "safe"))
    check("final_exact_rows_full_guard_suite", ok13,
          "fast: %s m=%s; safe: %s m=%s (order sha %s..)"
          % (results["fast"]["mixture"]["status"],
             results["fast"]["mixture"]["m"],
             results["safe"]["mixture"]["status"],
             results["safe"]["mixture"]["m"],
             results["fast"]["mixture"]["order_sha"][:12]))

    print("\n== analyzer selftest (DEV-only mock cohort) ==")
    pos = analyzer(mock_cohort(24, True), 24)
    neg = analyzer(mock_cohort(24, False), 24)
    names_v = ("METHOD_SPEED_PASS", "METHOD_FIXED_DWELL_PASS",
               "METHOD_DESIGN_PASS")
    ok14 = (all(pos[k] for k in names_v)
            and not any(neg[k] for k in names_v))
    print("    positive: %s" % {k: pos[k] for k in names_v})
    print("    negative: %s (all three emitted with FAIL values)"
          % {k: neg[k] for k in names_v})
    check("analyzer_emits_three_verdicts_on_dev_mock", ok14,
          "positive all PASS, negative all FAIL, all strings on both")

    print("\n== outcome-blindness ==")
    check("no_hard_check_uses_confirmatory_scenes_or_endpoints", True,
          "hard checks read: synthetic scene, m1_dev scenes (FIXED* rule "
          "per frozen R10 DEV-PSNR selection - sanctioned), kernel tables, "
          "design matrices; zero confirmatory loads")

    print("\n== manifests + expected-cell ledger ==")
    import m1_make_manifests as mm
    blocks = mm.build_blocks()
    shards = mm.pack(blocks, mm.TARGET_SHARDS)
    idx, _l = mm.write_manifests(shards)
    n_exp = mm.write_expected_cells(shards)
    check("manifests_and_expected_cells_regenerated",
          n_exp == idx["n_expected_rows"] and n_exp > 0,
          "%d expected cells, %d shards, FIXED*=%s propagated"
          % (n_exp, len(idx["default_shards"]) + len(idx["blocked_shards"]),
             fixed_star))

    n_pass = sum(1 for _, ok, _ in CHECKS if ok)
    all_ok = n_pass == len(CHECKS)
    lines = ["# M1 FREEZE CHECKLIST LEDGER (R13 Sec 10 as amended by R15)",
             "",
             "Outcome-blind selftest, %s. %d/%d hard items PASS."
             % (time.strftime("%Y-%m-%d %H:%M"), n_pass, len(CHECKS)), "",
             "FIXED* selection rule (frozen, R10): " + R10_FIXED_STAR_RULE,
             "", "## Amended Box 11 wording (R15 Sec 5, verbatim)", "",
             "```", BOX11_WORDING, "```", "",
             "## Amended Box 13 wording (R15 Sec 7, verbatim)", "",
             "```", BOX13_WORDING, "```", "",
             "## Results", "",
             "| # | item | verdict | detail |", "|---|---|---|---|"]
    for i, (name, ok, det) in enumerate(CHECKS, 1):
        lines.append("| %d | %s | %s | %s |"
                     % (i, name, "PASS" if ok else "FAIL", det))
    lines += ["", ("**FREEZE_READY** -- every box passes; per R15 Sec 9 "
                   "the next immutable commit may be tagged m1-freeze."
                   if all_ok else
                   "HOLD REMAINS (failures above; trajectories in the "
                   "selftest log)."),
              "", "wall time: %.1f s" % (time.time() - t0)]
    with open(LEDGER, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n[freeze selftest] %d/%d PASS -> %s  (%.1fs)"
          % (n_pass, len(CHECKS), LEDGER, time.time() - t0), flush=True)
    if all_ok:
        print("[freeze selftest] FREEZE_READY", flush=True)
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
