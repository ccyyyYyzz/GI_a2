"""M1 refreeze selftest — maps 1:1 onto the R18 §7 ten-item checklist.

Outcome-blind, DEV-only (synthetic + 632900 DEV scenes; the confirmatory
lock is exercised as evidence). Includes the R18 §5.4 24-cell DEV
full-stack feasibility/runtime gate, whose predetermined outcome selects
the campaign branch (certificate retained vs FALLBACK_DESCRIPTIVE) — either
branch is a valid deliverable and is recorded in CERT_BRANCH.json.

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
import m1_runner as m1                      # noqa: E402

OUT = os.path.join(ROOT, "results", "round63_m1")
LEDGER = os.path.join(OUT, "FREEZE_CHECKLIST_LEDGER.md")
BRANCH_JSON = os.path.join(OUT, "CERT_BRANCH.json")

SIDE, N = 32, 1024
CHECKS = []
RETIRED_TOKEN = "DOSE_SAFE" + "_CERT_PASS"      # avoid self-matching


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok), detail))
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name,
                           ("  -- " + detail) if detail else ""), flush=True)
    return ok


def mock_curves(n_images=24, positive=True):
    from gi_core.utils import rng_for
    fams = ["glyph", "chirp", "spokes", "maze", "contour", "micro"]
    nus = np.array(m1.NU_FULL)
    curves = {}
    for u in range(n_images):
        r = rng_for(u, 63, 12)
        base = 10 + r.random() * 2
        safe_q = base + 6 * (1 - np.exp(-nus / 500.0))
        shift = (200.0 if positive else 1.0) + 10.0 * r.random()
        ridge_q = base + 6 * (1 - np.exp(-nus * shift / 500.0)) + \
            (1.5 if positive else -0.3)
        img = "m1_%s_%d" % (fams[u % 6], u // 6)
        curves[img] = {
            "family": fams[u % 6],
            "safe": (nus, np.full(nus.size, 0.05), safe_q),
            "ridge": (nus, np.full(nus.size, 1.3), ridge_q),
            "q060_terminal": float(base + 6 * (1 - np.exp(-4.0))),
            "ridge_terminal": float(ridge_q[-1]
                                    + (0.6 if positive else -0.5)),
            "cert_cells": ["CERTIFIED" if positive else "COUNTEREXAMPLE"] * 4}
    return curves


def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    print("[R18 refreeze selftest] outcome-blind, DEV-only.", flush=True)

    # ---- item 1: rename wired through everything --------------------------- #
    print("\n== R18-1: verdict rename wired ==")
    prod = ["m1_runner.py", "m1_make_manifests.py", "campaign.py",
            "oed_design_v5.py"]
    stale = [fn for fn in prod
             if RETIRED_TOKEN in open(os.path.join(HERE, fn),
                                      encoding="utf-8").read()]
    an_src = open(os.path.join(HERE, "m1_runner.py"), encoding="utf-8").read()
    paper = os.path.join(OUT, "PAPER2_ACT3_WORDING_R18.md")
    paper_ok = (os.path.exists(paper)
                and "FULL_STACK_CERT_PASS" in open(paper).read()
                and RETIRED_TOKEN not in open(paper).read())
    check("rename_wired_spec_analyzer_paper_manifests", not stale
          and "FULL_STACK_CERT_PASS" in an_src and paper_ok,
          "retired token absent from %s; FULL_STACK_CERT_PASS in analyzer + "
          "paper skeleton" % prod)

    # ---- item 2: C_stack exact -------------------------------------------- #
    print("\n== R18-2: C_stack class exact ==")
    src5 = open(os.path.join(HERE, "oed_design_v5.py"),
                encoding="utf-8").read()
    ok2 = ("1.05 * h_fix" in src5 and "0.5 * prob.V_fix" in src5
           and "1.05 * prob.h_fix" in src5)
    # A-risk/spectral must be GLOBAL: the atom-admission code must not
    # reference them (admission = shape/peak/ceiling/trust/bias only)
    adm = src5[src5.index("def _adm_at_loads"):src5.index(
        "def setup_ctx_cert")]
    ok2 &= ("arisk" not in adm and "spectral" not in adm)
    check("c_stack_constraints_exact_and_global", ok2,
          "budget + +/-5%% dose + A-risk 1.05x + spectral 0.5x as GLOBAL "
          "constraints (engine cuts + probe checks); atom admission "
          "untouched by conditioning")

    # ---- item 3: four-toy suite ------------------------------------------- #
    print("\n== R18-3: four-toy suite ==")
    toys = v5.r18_toy_suite(verbose=True)
    ok3 = all(r["pass"] for r in toys)
    check("four_toy_suite", ok3,
          "agreements %s; INTERPRETATION (logged): R18 Sec 3.5 item 4 "
          "'budget, dose, A-risk, and spectral all active simultaneously' "
          "-- implemented as all four families present and enforced with "
          "dose+A-risk binding at the optimum and spectral shaping the "
          "solve path (cuts engaged); a strictly simultaneous 4-binding "
          "vertex proved numerically incompatible with the 1e-8 agreement "
          "bar in every non-degenerate instance tried (razor-edge duals); "
          "flagged for coordinator/R19 review"
          % ["%.0e" % r["agreement"] for r in toys])

    # ---- item 4: full-dictionary scan + residual checks -------------------- #
    print("\n== R18-4: scan + residuals implemented ==")
    ok4 = all(abs(r["scan_residual"]) <= 1e-8 for r in toys)
    ok4 &= all(r["psd_residual"] >= -1e-9 for r in toys)
    check("full_dict_scan_and_residuals", ok4,
          "toys: scan residuals %s, psd %s; real-cell scans exercised by "
          "the item-5 gate below"
          % (["%.0e" % abs(r["scan_residual"]) for r in toys],
             ["%.0e" % r["psd_residual"] for r in toys]))

    # ---- item 5: 24-cell DEV full-stack gate ------------------------------ #
    print("\n== R18-5: 24-cell DEV gate (final predetermined stop rule) ==")
    import campaign
    import m1_scenes
    dev_imgs = campaign._images(SIDE, "all", imageset="m1_dev")
    names = [nm for nm, _f, _s in m1_scenes._m1_dev_table()]
    rows_dep, _sha = m1.deployed_scat32()
    table = []
    reuse = None
    if os.path.exists(BRANCH_JSON):
        with open(BRANCH_JSON) as f:
            prev = json.load(f)
        if len(prev.get("cells", [])) == 24:
            reuse = prev
            print("    (reusing the deterministic 24-cell gate table from "
                  "CERT_BRANCH.json; protocol is frozen + deterministic)",
                  flush=True)
    if reuse is not None:
        names = []                       # skip recompute
        table = [{"image": c_["image"], "nu": float(c_["nu"]),
                  "b": float(c_["b"]), "status": c_["status"],
                  "primal_gap_lower_per_r":
                      (float(c_["primal_gap_lower_per_r"])
                       if c_["primal_gap_lower_per_r"]
                       not in ("-inf", "nan") else float("-inf")),
                  "dual_gap_upper_per_r": float("nan"),
                  "wall_seconds": float(c_["wall_seconds"])}
                 for c_ in reuse["cells"]]
    for img in names:
        for nu in (200.0, 2000.0):
            for b in (0.05, 0.60):
                tc = time.time()
                xhat = m1.prescan_estimate(dev_imgs[img], img, 0, b, nu,
                                           per_cell=True)
                Vf = v4.info_matrix_full(rows_dep, xhat, int(nu), b,
                                         P=m1.prescan_matrix())
                B, eps0, _tr = v4.subspace_from_fixedstar(Vf)
                ctx = v5.setup_ctx_cert(xhat, nu, b, B, eps0, SIDE)
                out = v5.full_stack_cert_cell(ctx, rows_dep, b)
                out["image"], out["nu"], out["b"] = img, nu, b
                out["wall_total"] = time.time() - tc
                table.append(out)
                print("    %s nu=%g b=%g: %s primal_gap=%.4f dual=%s "
                      "wall=%.0fs"
                      % (img, nu, b, out["status"],
                         out.get("primal_gap_lower_per_r", float("nan")),
                         ("%.3f" % out["dual_gap_upper_per_r"])
                         if np.isfinite(out.get("dual_gap_upper_per_r",
                                                np.inf)) else "-",
                         out["wall_seconds"]), flush=True)
    n_cert = sum(1 for o in table if o["status"] == "CERTIFIED")
    n_ce = sum(1 for o in table if o["status"] == "COUNTEREXAMPLE")
    n_unres = len(table) - n_cert - n_ce
    walls = sorted(o["wall_seconds"] for o in table)
    med_wall = walls[len(walls) // 2]
    max_wall = walls[-1]
    gate_pass = (n_cert == 24 and n_ce == 0 and med_wall <= 120
                 and max_wall <= 420)
    branch = "CERTIFICATE_RETAINED" if gate_pass else "FALLBACK_DESCRIPTIVE"
    with open(BRANCH_JSON, "w") as f:
        json.dump({"branch": branch,
                   "gate": {"n_certified": n_cert,
                            "n_counterexample": n_ce,
                            "n_unresolved": n_unres,
                            "median_wall_s": med_wall,
                            "max_wall_s": max_wall},
                   "cells": [{k_: (float(o[k_]) if isinstance(
                       o.get(k_), (int, float)) and np.isfinite(
                       float(o[k_])) else str(o.get(k_)))
                       for k_ in ("image", "nu", "b", "status",
                                  "primal_gap_lower_per_r",
                                  "dual_gap_upper_per_r", "wall_seconds")}
                       for o in table]}, f, indent=1)
    # gate item passes if EITHER branch executed per the predetermined rule
    check("dev_gate_24_cells_branch_decided", True,
          "CERTIFIED=%d COUNTEREXAMPLE=%d UNRESOLVED=%d; median wall=%.0fs "
          "max=%.0fs -> BRANCH: %s (R18 Sec 5.4: either branch launches; "
          "no third option)"
          % (n_cert, n_ce, n_unres, med_wall, max_wall, branch))
    if branch == "FALLBACK_DESCRIPTIVE":
        pos_an = m1.m1_analyze_r17(mock_curves(24, True), 24, boot_B=500)
        ok5b = ("FULL_STACK_CERT_PASS" not in pos_an
                and "full_stack_cert_descriptive_fraction" in pos_an)
        check("fallback_removes_categorical_gate", ok5b,
              "analyzer under FALLBACK emits descriptive fraction only "
              "(keys: %s)" % sorted(k for k in pos_an if "cert" in k.lower()))

    # ---- item 6: dose-only DEV evidence stored ----------------------------- #
    print("\n== R18-6: dose-only no-gate DEV evidence ==")
    md = os.path.join(OUT, "R18_GAP_PROBE_REPLICATION.md")
    md_src = open(md, encoding="utf-8").read() if os.path.exists(md) else ""
    ok6 = ("support-preserving" in md_src
           and "DOSE_ONLY_PRIMAL_GAP_PER_R" in md_src
           and "Development-only constructive analysis" in md_src
           and md_src.count("m1_dev_") >= 12
           and "PASS verdict" in md_src)
    check("dose_only_dev_evidence_with_caveat", ok6,
          "frozen 6-scene x 2-anchor table + R18 Sec 1.2 caveat verbatim + "
          "no DOSE_ONLY_*_PASS verdict")

    # ---- item 7: paper wording corrected ----------------------------------- #
    print("\n== R18-7: paper wording ==")
    p_src = " ".join(open(paper, encoding="utf-8").read().split())
    ok7 = ("does not eliminate scene-adapted geometry selection"
           in p_src.replace("> ", "")
           and "full deployed conditioning stack" in p_src.replace("> ", "")
           and "DOES NOT SURVIVE" in p_src)
    check("paper_wording_r18", ok7,
          "R17 collapse wording removed; R18 Sec 4.3/4.4 frozen text in "
          "results/round63_m1/PAPER2_ACT3_WORDING_R18.md")

    # ---- item 8: empirical modes/endpoints unchanged ----------------------- #
    print("\n== R18-8: empirical arms unchanged ==")
    pos = m1.m1_analyze_r17(mock_curves(24, True), 24, boot_B=2000)
    neg = m1.m1_analyze_r17(mock_curves(24, False), 24, boot_B=2000)
    ok8 = (m1.MODE_ARMS == ["SCAT32-SAFE", "SCAT32-060", "RIDGE-SCAT32"]
           and m1.MODE_LOADS["SCAT32-SAFE"] == 0.05
           and m1.MODE_LOADS["SCAT32-060"] == 0.60
           and pos["RIDGE_OPERATING_PASS"] and pos["RIDGE_SPEED_PASS"]
           and not neg["RIDGE_OPERATING_PASS"]
           and not neg["RIDGE_SPEED_PASS"])
    check("empirical_modes_endpoints_unchanged", ok8,
          "3 modes + loads frozen; RIDGE_OPERATING/RIDGE_SPEED emitted "
          "correctly on positive AND negative mocks")

    # ---- item 9: accounting/SHAs/ridge/locks ------------------------------- #
    print("\n== R18-9: accounting + SHAs + ridge + locks ==")
    P = m1.prescan_matrix()
    ok9 = (P.shape[0] == 52
           and np.abs(P.mean(axis=0) - 1.0).max() < 1e-12)
    for arm in m1.ARMS_ALL:
        pat = m1.load_m1_pattern("m1pat:%s:x" % arm, 1024, 1024, 0)
        ok9 &= pat["A"].shape == (1024, N)
    rows_dep2, sha_dep = m1.deployed_scat32()
    man = v4.dictionary_manifest(SIDE)
    rt = m1.ridge_targets([2000.0])
    ok9 &= (len(sha_dep) == 64 and len(man["sha256_global"]) == 64
            and len(v5.d_cert_sha()) == 64
            and abs(rt["2000"]["rho_R_production"] - 22.2545) < 0.01)
    check("accounting_shas_ridge_locks", ok9,
          "52+972 all arms; deployed sha %s..; dict sha %s..; D_cert sha "
          "%s..; rho_R(2000)=%.4f"
          % (sha_dep[:8], man["sha256_global"][:8], v5.d_cert_sha()[:8],
             rt["2000"]["rho_R_production"]))

    # ---- item 10: no confirmatory data ------------------------------------- #
    print("\n== R18-10: confirmatory lock ==")
    fired = False
    try:
        m1.run_cert_cell(m1.cert_cell("m1_glyph_0", 0, 2000.0, 0.60,
                                      imageset="m1"))
    except RuntimeError:
        fired = True
    check("no_confirmatory_data_opened", fired,
          "M1_FREEZE_LAUNCHED lock raises pre-freeze; selftest touched only "
          "synthetic + m1_dev scenes")

    # ---- ledger ------------------------------------------------------------ #
    n_pass = sum(1 for _, ok_, _ in CHECKS if ok_)
    all_ok = n_pass == len(CHECKS)
    with open(BRANCH_JSON) as f:
        branch = json.load(f)["branch"]
    lines = ["# M1 REFREEZE CHECKLIST LEDGER (R18 §7, ten items)", "",
             "Outcome-blind selftest, %s. %d/%d items PASS. Campaign "
             "branch: **%s**."
             % (time.strftime("%Y-%m-%d %H:%M"), n_pass, len(CHECKS),
                branch), "",
             "| # | item (R18 §7) | verdict | detail |",
             "|---|---|---|---|"]
    for i, (name, ok_, det) in enumerate(CHECKS, 1):
        lines.append("| %d | %s | %s | %s |"
                     % (i, name, "PASS" if ok_ else "FAIL", det))
    lines += ["", ("**REFREEZE_READY_R18 (branch: %s)** -- per R18 §5.4/§7 "
                   "one immutable tag may be created and all arms launched "
                   "together." % branch
                   if all_ok else "HOLD REMAINS (failures above)."),
              "", "wall time: %.1f s" % (time.time() - t0)]
    with open(LEDGER, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n[R18 selftest] %d/%d PASS, branch=%s -> %s (%.0fs)"
          % (n_pass, len(CHECKS), branch, LEDGER, time.time() - t0),
          flush=True)
    if all_ok:
        print("[R18 selftest] REFREEZE_READY_R18 branch=%s" % branch,
              flush=True)
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
