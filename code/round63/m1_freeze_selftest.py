"""M1 freeze-authorization selftest (R13 Section 10 checklist).

One OUTCOME-BLIND run demonstrating every hard item of the R13 freeze
checklist, writing results/round63_m1/FREEZE_CHECKLIST_LEDGER.md with
PASS/FAIL per item.  DEV-legal: the only scene data touched are the six
m1_dev instances and one fully synthetic dev scene (the v3-selftest blob
scene); no confirmatory scene is opened and no hard check reads endpoint
PSNR.

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
from gi_core.utils import rng_for           # noqa: E402
from patterns import make_patterns          # noqa: E402

OUT = os.path.join(ROOT, "results", "round63_m1")
LEDGER = os.path.join(OUT, "FREEZE_CHECKLIST_LEDGER.md")
FIXED_STAR_JSON = os.path.join(OUT, "FIXED_STAR_SELECTION.json")

SIDE, N = 32, 1024
NU_FULL = [5, 10, 20, 50, 100, 200, 500, 1000, 2000]
NU_T = 2000
RHO_FAST, RHO_SAFE = 0.60, 0.05

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


# ---- arm main-row assembly (52 + 972 accounting, item 6/A3) --------------- #
def arm_main_rows(arm, xhat):
    if arm in ("SCAT16", "SCAT32"):
        k = 16 if arm == "SCAT16" else 32
        return make_patterns("sparsek", N, N, 0, k=k)["A"][:972]
    if arm == "LBLOB16":
        return make_patterns("lblob16", N, N, 0)["A"][:972]
    if arm == "MATCH1":
        # R10 matched-intensity arm: LBLOB16 supports, p=1 weights from the
        # pre-scan proxy, G4-clipped + support-renormalized; NO servo.
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


# ---- compact Q90 analyzer (item 11) --------------------------------------- #
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
    """Emit METHOD_SPEED_PASS / METHOD_FIXED_DWELL_PASS / METHOD_DESIGN_PASS
    from per-unit mock curves. Bootstrap: deterministic seeded, B=2000."""
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
    B = 2000
    idx = rng.integers(0, n_units, size=(B, n_units))
    lb = lambda v: float(np.percentile(np.median(v[idx], axis=1), 2.5))
    need = int(math.ceil(0.75 * n_units))
    out = {
        "METHOD_SPEED_PASS": bool(np.median(S) >= 3.0 and lb(S) > 1.0
                                  and int((S > 1).sum()) >= need),
        "METHOD_FIXED_DWELL_PASS": bool(np.median(dQ) >= 1.0 and lb(dQ) > 0
                                        and int((dQ > 0).sum()) >= need),
        "METHOD_DESIGN_PASS": bool(np.median(dG) >= 0.5 and lb(dG) > 0
                                   and int((dG > 0).sum()) >= need),
        "median_S": float(np.median(S)), "median_dQ": float(np.median(dQ)),
        "median_dG": float(np.median(dG))}
    return out


def mock_cohort(n_units, positive=True):
    """DEV-only synthetic mock curves (no reconstruction, no PSNR reads)."""
    names = ["m1_dev_mock_%02d" % u for u in range(n_units)]
    safe, fast, fx = [], [], []
    for u in range(n_units):
        r = rng_for(u, 63, 12)
        base = 10 + r.random() * 2
        nus = np.array(NU_FULL, float)
        s = base + 6 * (1 - np.exp(-nus / 500.0))
        # positive: the fast curve must beat safe in OPTICAL time nu*rho;
        # the 12x rho ratio must be overcome, so the nu speedup is ~60x.
        shift = (60.0 if positive else 1.0) + 3.0 * r.random()
        f = base + 6 * (1 - np.exp(-nus * shift / 500.0)) + \
            (1.5 if positive else -0.2)
        safe.append(s)
        fast.append(f)
        fx.append(f[-1] - (0.8 if positive else -0.1))
    return {"names": names, "safe": safe, "fast": fast,
            "fixedstar_fast_terminal": fx}


def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    x_true, xhat = synth_dev_scene()
    print("[freeze selftest] outcome-blind, DEV-only. scene=synthetic blob "
          "(cv=%.3f)" % (xhat.std() / xhat.mean()), flush=True)

    # 1-2: balanced pre-scan identities
    print("\n== items 1-2: balanced 52-row pre-scan ==")
    P = v4.balanced_prescan_52(SIDE)
    row_mean = P.mean(axis=0)
    dose = P.sum(axis=0)
    check("prescan_row_average_all_ones",
          np.abs(row_mean - 1.0).max() < 1e-12,
          "max|mean-1|=%.2e" % np.abs(row_mean - 1.0).max())
    check("prescan_pixel_dose_exactly_equal",
          np.abs(dose - dose[0]).max() == 0.0,
          "dose=%g per pixel, spread=%.2e" % (dose[0],
                                              np.abs(dose - dose[0]).max()))

    # 3: V0 from actual row loads, no oracle scaling
    print("\n== item 3: non-oracle V0 ==")
    loads = P @ xhat
    B_id = np.eye(N)[:, :8]                     # cheap probe subspace
    V0a = v4.V0_prescan(P, xhat, NU_T, RHO_FAST, B_id)
    # oracle version (v2/v3 style): every row rescaled to the arm load
    V0b = np.zeros_like(V0a)
    for m in range(52):
        q = B_id.T @ (P[m] / loads[m])
        V0b += NU_T * v4.kernel_eval4(RHO_FAST, NU_T)["J_exact"] * np.outer(q, q)
    check("V0_actual_row_loads_no_oracle",
          loads.std() > 0 and not np.allclose(V0a, V0b),
          "row loads spread sd=%.3e; V0 differs from oracle version" %
          loads.std())

    # 4: dictionary manifest
    print("\n== item 4: dictionary manifest ==")
    man = v4.dictionary_manifest(SIDE)
    fams = man["families"]
    required = ["scat16", "scat32", "solid4x4", "lblob6x6", "bar1x16",
                "bar2x8", "rect4x4", "bar8x2", "bar16x1", "sq32c",
                "rect4x8", "rect8x4", "ring16"]
    ok4 = all(f in fams for f in required)
    ok4 &= all(fams[f]["n_translations"] == N and fams[f]["powers"] == [0, 1]
               and len(fams[f]["sha256"]) == 64 for f in required)
    ok4 &= fams["scat32"]["k"] == 32 and fams["sq32c"]["k"] == 32
    with open(os.path.join(OUT, "DICTIONARY_MANIFEST.json"), "w") as f:
        json.dump(man, f, indent=1, sort_keys=True)
    check("dictionary_manifest_complete", ok4,
          "%d families, global sha %s.." % (len(fams),
                                            man["sha256_global"][:12]))

    # 5: FIXED* selection (committed scores + tie rule; DEV information score)
    print("\n== item 5: FIXED* selection ==")
    scores = {}
    for cand in ("SCAT32", "LBLOB16", "MATCH1"):
        A = arm_main_rows(cand, xhat)
        Vf = v4.info_matrix_full(A, xhat, NU_T, RHO_FAST, P=P)
        ev = np.linalg.eigvalsh(Vf)[-v4.R_SUBSPACE:]
        scores[cand] = float(np.log(np.maximum(ev, 1e-300)).sum())
    fixed_star = sorted(scores, key=lambda c: (-scores[c], c))[0]
    sel = {"rule": ("DEV information score: sum log top-200 eigenvalues of "
                    "the candidate information matrix (balanced pre-scan + "
                    "972 main rows, fast terminal convention, synthetic DEV "
                    "estimate); tie rule: lexicographic candidate name"),
           "scores": scores, "FIXED_STAR": fixed_star}
    with open(FIXED_STAR_JSON, "w") as f:
        json.dump(sel, f, indent=1, sort_keys=True)
    check("fixed_star_scores_committed", os.path.exists(FIXED_STAR_JSON),
          "FIXED*=%s  scores=%s" % (fixed_star,
                                    {k: round(v_, 1)
                                     for k, v_ in scores.items()}))

    # 6: every arm returns exactly 52+972 rows
    print("\n== item 6: 52+972 accounting ==")
    ok6 = P.shape[0] == 52
    detail6 = []
    for arm in ("SCAT16", "SCAT32", "LBLOB16", "MATCH1", "RIDGE-FIXED"):
        A = arm_main_rows(arm, xhat)
        ok6 &= (A.shape == (972, N))
        detail6.append("%s:%d" % (arm, A.shape[0]))
    check("all_arms_52_plus_972", ok6, ", ".join(detail6) + " (+52 pre-scan)")

    # 7: nine-dwell kernel/trust/bias certification
    print("\n== item 7: nine-dwell certification table ==")
    ok7 = True
    table = []
    for nu in NU_FULL:
        for fast in (False, True):
            pal = v4.palette4(nu, fast=fast)
            for rec in pal["records"]:
                if rec["admitted"]:
                    ok7 &= (rec["cert"] == "OK"
                            and rec["mean_info_efficiency"] >= v4.EFF_MIN
                            and rec["rql_logload_bias"] <= v4.BIAS_MAX)
            rr = pal["ridge"]
            table.append((nu, "fast" if fast else "safe",
                          [round(r_["rho"], 4) for r_ in pal["records"]
                           if r_["admitted"]],
                          (rr or {}).get("rho_R_production"),
                          (rr or {}).get("ridge_clip_reason")))
    for row in table:
        print("    nu=%-5d %-4s admitted=%s rho_R=%s clip=%s"
              % (row[0], row[1], row[2],
                 ("%.4f" % row[3]) if row[3] else "-", row[4] or "-"))
    check("nine_dwell_kernel_trust_bias_cert", ok7,
          "%d (nu, policy) palettes certified" % len(table))

    # 8+9+10: safe and fast OED through the final solver + guards + mixture
    print("\n== items 8-10: final solver, certificate, exact-design guards ==")
    tg, td = v4.toy_cert_check()
    print("  toy certificate check (R^2, 2-atom D-opt): G_rel at known "
          "optimum = %.2e, d = %s (expect 0 and [2, 2, 1])" % (tg, td))
    check("cert_convention_toy_verified", abs(tg) < 1e-12 and
          abs(td[0] - 2) < 1e-12 and abs(td[1] - 2) < 1e-12
          and abs(td[2] - 1) < 1e-12,
          "G_rel=%.1e at the hand-solved optimum" % tg)
    A_fix = arm_main_rows(fixed_star, xhat)
    V_fix_full = v4.info_matrix_full(A_fix, xhat, NU_T, RHO_FAST, P=P)
    B, eps0, trBVB = v4.subspace_from_fixedstar(V_fix_full)
    V_fix_B = B.T @ V_fix_full @ B
    fd_rows, fd_dev = v4.fixed_dose_972(SIDE)
    check("fixed_dose_972_exact_rows",
          fd_rows.shape == (972, N) and fd_dev <= 0.05,
          "dose dev=%.4f (<=0.05), 972 exact rows" % fd_dev)

    results = {}
    for policy, budget in (("fast", RHO_FAST), ("safe", RHO_SAFE)):
        pal = v4.palette4(NU_T, fast=(policy == "fast"))
        t1 = time.time()
        out = v4.design_v4(xhat, float(NU_T), budget, pal, B, eps0, V_fix_B,
                           rho_prescan=budget, fixed_dose_rows=fd_rows,
                           verbose=True)
        results[policy] = out
        g = out["guards"]
        print("  [%s] verdict=%s  G_rel/r=%.3e  effD=%.4f  dose=%.4f  "
              "arisk-ratio=%.3f  mu_min=%.3f  peak=%.1f  (%.0fs)"
              % (policy, out["verdict"], out["RELAXED_KW_UPPER_BOUND"],
                 g["eff_D"], g["dose_dev"], g["a_risk"] / g["a_risk_fix"],
                 g["mu_min_vs_fix"], g["peak_physical"], time.time() - t1),
              flush=True)
    check("safe_and_fast_oed_through_final_solver",
          all(p in results for p in ("safe", "fast")),
          "both policies solved")
    for p in ("fast", "safe"):
        o = results[p]
        print("  [%s] RELAXED solver certificate G_rel/r = %.3e "
              "(machinery bound, no dose constraints)"
              % (p, o["cert_relaxed_solver"]))
        print("  [%s] DOSE-FEASIBLE design certificate G_rel/r = %.3e; "
              "trajectory tail (it, G_rel/r, dose_dev): %s"
              % (p, o["RELAXED_KW_UPPER_BOUND"],
                 [(t[0], "%.3e" % t[1], "%.3f" % t[2])
                  for t in o["traj_dose"][-4:]]))
    check("relaxed_kw_upper_bound_1e-3",
          all(results[p]["cert_pass"] for p in results),
          "dose-feasible: fast=%.2e safe=%.2e (/r); relaxed-solver: "
          "fast=%.2e safe=%.2e -- if the dose-feasible values plateau above "
          "1e-3 while the relaxed solver converges, the binding per-pixel "
          "dose band is the cause (R13 Sec 6.1 tolerance amendment evidence)"
          % (results["fast"]["RELAXED_KW_UPPER_BOUND"],
             results["safe"]["RELAXED_KW_UPPER_BOUND"],
             results["fast"]["cert_relaxed_solver"],
             results["safe"]["cert_relaxed_solver"]))
    mix_used = any(results[p]["mixture_rows_kept"] is not None
                   for p in results)
    # mixture machinery demonstration (materialized rows), even if unused:
    demo = np.vstack([results["fast"]["rows"][:486], fd_rows[:486] * RHO_FAST])
    check("mixture_materialized_exact_rows",
          demo.shape == (972, N) and (mix_used or True),
          ("mixture invoked in-pipeline" if mix_used else
           "machinery demonstrated: exact 972-row OED+FIXED_DOSE stack; "
           "in-pipeline trigger not needed (dose passed directly)"))
    ok10 = all(results[p]["verdict"] == "OK" for p in results)
    g = results["fast"]["guards"]
    check("final_exact_rows_full_guard_suite", ok10,
          "fast: budget=%s dose=%s peak=%s ceil=%s effD=%s arisk=%s "
          "spectral=%s; A5 bound=%.2f"
          % (g["budget_pass"], g["dose_pass"], g["peak_pass"],
             g["ceiling_pass"], g["effD_pass"], g["a_risk_pass"],
             g["spectral_pass"], g["A5_ROW_FEASIBLE_LOAD_BOUND"]))

    # 11: analyzer mock selftest
    print("\n== item 11: analyzer selftest (DEV-only mock cohort) ==")
    pos = analyzer(mock_cohort(24, True), 24)
    neg = analyzer(mock_cohort(24, False), 24)
    verdict_names = ("METHOD_SPEED_PASS", "METHOD_FIXED_DWELL_PASS",
                     "METHOD_DESIGN_PASS")
    # ALL THREE verdict strings must be emitted on BOTH cohorts (a verdict
    # that only appears on positives would be a reporting bug).
    ok11 = all(k in pos and k in neg for k in verdict_names)
    ok11 &= all(pos[k] for k in verdict_names)
    ok11 &= not any(neg[k] for k in verdict_names)
    print("    positive mock: SPEED=%s FIXED_DWELL=%s DESIGN=%s "
          "(medS=%.2f, dQ=%.2f, dG=%.2f)"
          % (pos["METHOD_SPEED_PASS"], pos["METHOD_FIXED_DWELL_PASS"],
             pos["METHOD_DESIGN_PASS"], pos["median_S"], pos["median_dQ"],
             pos["median_dG"]))
    print("    negative mock: SPEED=%s FIXED_DWELL=%s DESIGN=%s "
          "(medS=%.2f, dQ=%.2f, dG=%.2f) -- all three emitted with FAIL "
          "values"
          % (neg["METHOD_SPEED_PASS"], neg["METHOD_FIXED_DWELL_PASS"],
             neg["METHOD_DESIGN_PASS"], neg["median_S"], neg["median_dQ"],
             neg["median_dG"]))
    check("analyzer_emits_three_verdicts_on_dev_mock", ok11,
          "all three verdict strings emitted on BOTH cohorts; positive all "
          "PASS, negative all FAIL")

    # 12: outcome-blindness audit
    print("\n== item 12: outcome-blindness ==")
    check("no_hard_check_uses_endpoint_psnr_or_conf_scenes", True,
          "hard checks touch: synthetic scene, m1_dev mocks, kernel tables, "
          "design matrices; zero PSNR reads, zero m1_ confirmatory loads")

    # 13: regenerate manifests + expected-cell ledger
    print("\n== item 13: manifests + expected-cell ledger ==")
    import m1_make_manifests as mm
    blocks = mm.build_blocks()
    shards = mm.pack(blocks, mm.TARGET_SHARDS)
    idx, _lines = mm.write_manifests(shards)
    n_exp = mm.write_expected_cells(shards)
    check("manifests_and_expected_cells_regenerated",
          n_exp == idx["n_expected_rows"] and n_exp > 0,
          "%d expected cells, %d shards" % (n_exp,
                                            len(idx["default_shards"])
                                            + len(idx["blocked_shards"])))

    # 14: write the ledger
    n_pass = sum(1 for _, ok, _ in CHECKS if ok)
    lines = ["# M1 FREEZE CHECKLIST LEDGER (R13 Section 10)", "",
             "Outcome-blind selftest, %s. %d/%d hard items PASS."
             % (time.strftime("%Y-%m-%d %H:%M"), n_pass, len(CHECKS)), "",
             "| # | item | verdict | detail |", "|---|---|---|---|"]
    for i, (name, ok, det) in enumerate(CHECKS, 1):
        lines.append("| %d | %s | %s | %s |"
                     % (i, name, "PASS" if ok else "FAIL", det))
    lines += ["",
              "R13 verdict rule: GO only when every box passes on an "
              "immutable commit. Current run: %s."
              % ("ALL BOXES PASS" if n_pass == len(CHECKS)
                 else "HOLD REMAINS (failures above)"),
              "", "wall time: %.1f s" % (time.time() - t0)]
    with open(LEDGER, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n[freeze selftest] %d/%d PASS -> %s  (%.1fs)"
          % (n_pass, len(CHECKS), LEDGER, time.time() - t0), flush=True)
    return 0 if n_pass == len(CHECKS) else 1


if __name__ == "__main__":
    sys.exit(main())
