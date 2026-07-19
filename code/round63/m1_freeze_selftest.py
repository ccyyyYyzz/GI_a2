"""M1 refreeze-authorization selftest — maps 1:1 onto R17 §7 / amendment §E.

Outcome-blind, DEV-only (synthetic scene + m1_dev 632900+ scenes; the
confirmatory guard is itself exercised as evidence). Writes
results/round63_m1/FREEZE_CHECKLIST_LEDGER.md with the nine §E items.

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

SIDE, N = 32, 1024
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


def mock_curves(n_images=24, positive=True):
    """DEV-only synthetic per-image curve set for the analyzer emission test
    (§E item 3). Family labels = 6 families x 4 instances."""
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
        rho_r = np.full(nus.size, 1.3)
        img = "m1_%s_%d" % (fams[u % 6], u // 6)
        curves[img] = {
            "family": fams[u % 6],
            "safe": (nus, np.full(nus.size, 0.05), safe_q),
            "ridge": (nus, rho_r, ridge_q),
            "q060_terminal": float(base + 6 * (1 - np.exp(-4.0))),
            "ridge_terminal": float(ridge_q[-1] + (0.6 if positive else -0.5)),
            "cert_cells": [positive] * 4}
    return curves


def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    x_true, xhat = synth_dev_scene()
    print("[refreeze selftest R17] outcome-blind, DEV-only.", flush=True)

    # ---- E1: retired semantics removed ------------------------------------ #
    print("\n== E1: retirements ==")
    prod_files = ["m1_runner.py", "m1_make_manifests.py", "campaign.py"]
    bad = []
    for fn in prod_files:
        src = open(os.path.join(HERE, fn), encoding="utf-8").read()
        for tok in ("METHOD_SPEED_PASS", "COMPLIANT_VIA_MIXTURE"):
            if tok in src:
                bad.append("%s:%s" % (fn, tok))
    ok = not bad
    ok &= not (set(m1.ARMS_ALL) & m1.RETIRED_ARMS)
    try:
        m1.load_m1_pattern("m1pat:OED-DT:x", 1024, 1024, 0)
        ok = False
        ret_msg = "OED-DT pattern NOT rejected"
    except ValueError as e:
        ret_msg = "retired-arm loader raises: %s" % str(e)[:40]
    pfa0 = v5.path_feasible_alpha({"m": 0})
    ok &= (pfa0["verdict"] == "ADAPTIVE_COLLAPSE_UNDER_GUARDS"
           and pfa0["is_pass"] is False and pfa0["dev_supplement_only"])
    check("oed_dt_and_m0_pass_semantics_removed", ok,
          "prod tokens absent%s; ARMS=%s; %s; m=0 -> %s (is_pass=%s)"
          % ("" if not bad else " EXCEPT %s" % bad, m1.ARMS_ALL, ret_msg,
             pfa0["verdict"], pfa0["is_pass"]))

    # ---- E2: deployed matrix + three modes frozen and hashed -------------- #
    print("\n== E2: deployed SCAT32 + operating modes ==")
    rows, sha = m1.deployed_scat32()
    dose = rows.sum(axis=0)
    dev = float(np.abs(dose / dose.mean() - 1.0).max())
    P = m1.prescan_matrix()
    import hashlib
    sha_p = hashlib.sha256(np.ascontiguousarray(P).tobytes()).hexdigest()
    ok = (rows.shape == (972, N) and dev <= 0.05
          and m1.MODE_LOADS["SCAT32-SAFE"] == 0.05
          and m1.MODE_LOADS["SCAT32-060"] == 0.60)
    check("deployed_scat32_and_modes_frozen_hashed", ok,
          "972x1024 sha=%s.. dose_dev=%.4f; modes SAFE=0.05 060=0.60 "
          "RIDGE=rho_R(nu) policy; prescan sha=%s.." % (sha[:12], dev,
                                                        sha_p[:12]))

    # ---- E3: verdict emission on positive AND negative mocks/toys --------- #
    print("\n== E3: analyzer + certificate verdict emission ==")
    pos = m1.m1_analyze_r17(mock_curves(24, True), 24, boot_B=2000)
    neg = m1.m1_analyze_r17(mock_curves(24, False), 24, boot_B=2000)
    names_v = ("RIDGE_OPERATING_PASS", "RIDGE_SPEED_PASS",
               "DOSE_SAFE_CERT_PASS")
    print("    positive: %s (medS=%.1f dQ=%.2f)"
          % ({k: pos[k] for k in names_v}, pos["median_S"],
             pos["median_dQ_dB"]))
    print("    negative: %s (all three emitted with FAIL values)"
          % {k: neg[k] for k in names_v})
    ok = all(pos[k] for k in names_v) and not any(neg[k] for k in names_v)
    tg = v5.toy_full_cert_check_gain()
    # negative certificate toy: the dual bound at a deliberately suboptimal
    # design must exceed the gate scale (verdict machinery emits FAIL)
    w_bad = np.array([1.0, 0.0, 0.0, 0.0])
    X = np.array([[1.0, 0.0], [0.0, 1.0],
                  [1 / math.sqrt(2), 1 / math.sqrt(2)],
                  [2 / math.sqrt(2), 2 / math.sqrt(2)]])
    Mm = np.outer(X[0], X[0]) + 1e-9 * np.eye(2)
    d_bad = np.array([x @ np.linalg.inv(Mm) @ x for x in X])
    G_bad = float(d_bad.max() - d_bad @ w_bad)      # theta=mu=0 dual point
    ok &= tg["agreement_ok"] and tg["dose_active"] and G_bad > 1.0
    print("    cert toys: optimum G=%.1e (PASS emission), suboptimal "
          "G>=%.1f (FAIL emission)" % (tg["G_at_opt"], G_bad))
    check("three_verdicts_emitted_pos_and_neg", ok,
          "RIDGE_OPERATING/RIDGE_SPEED/DOSE_SAFE_CERT on both cohorts; "
          "gain toy %.1e; negative toy G=%.1f" % (tg["G_at_opt"], G_bad))

    # ---- E4: D_cert complete + SHA-frozen --------------------------------- #
    print("\n== E4: D_cert = D_load U D_gain ==")
    t1 = time.time()
    ctxc = v5.setup_ctx_cert(xhat, 2000.0, 0.60, np.eye(N)[:, :200], 1e-6,
                             SIDE)
    sha_cert = v5.d_cert_sha(SIDE)
    n_load = int(ctxc["ALLOW"][:, :ctxc["L_load"]].sum())
    n_gain = int(ctxc["ALLOW"][:, ctxc["L_load"]:].sum())
    ok = (ctxc["L"] == ctxc["L_load"] + 5
          and ctxc["gammas"] == [0.2, 0.5, 1.0, 2.0, 5.0]
          and n_gain > 0 and n_load > 0)
    with open(os.path.join(OUT, "D_CERT_SHA.json"), "w") as f:
        json.dump({"sha256": sha_cert, "gammas": ctxc["gammas"],
                   "n_load_atoms_dev": n_load, "n_gain_atoms_dev": n_gain},
                  f, indent=1)
    check("d_cert_complete_and_sha_frozen", ok,
          "L=%d (load %d + gain 5); admissible: %d load + %d gain atoms; "
          "sha=%s.. (%.0fs)" % (ctxc["L"], ctxc["L_load"], n_load, n_gain,
                                sha_cert[:16], time.time() - t1))

    # ---- E5: expanded-class toy + scan tests ------------------------------ #
    print("\n== E5: expanded-class certificate toys ==")
    tb = v5.toy_full_cert_check()
    ok = (tb["agreement_ok"] and tg["agreement_ok"]
          and tb["dose_active"] and tg["dose_active"])
    check("expanded_class_cert_toys", ok,
          "base toy %.1e, gain-atom toy %.1e, dose band ACTIVE in both"
          % (tb["G_at_opt"], tg["G_at_opt"]))

    # ---- E6: DEV-only expanded-class feasibility + runtime ---------------- #
    print("\n== E6: DEV certificate cells (m1_dev scenes) ==")
    import m1_scenes
    dev_img = [nm for nm, _f, _s in m1_scenes._m1_dev_table()][0]
    results6 = []
    for (nu_c, b_c) in ((2000.0, 0.60), (200.0, 0.05)):
        cell = m1.cert_cell(dev_img, 0, nu_c, b_c, imageset="m1_dev")
        t2 = time.time()
        row = m1.run_cert_cell(cell)[0]
        row["wall_s"] = round(time.time() - t2, 1)
        results6.append(row)
        print("    (%s, s0, nu=%g, b=%g): status=%s G_full/r=%s "
              "primal=%s peak=%.1f wall=%.0fs"
              % (dev_img, nu_c, b_c, row["status"], row["G_full_per_r"],
                 row["primal_feasible"], row["deployed_peak"] or -1,
                 row["wall_s"]), flush=True)
    # HARD: the pipeline must terminate WITH a certified finite bound
    # (status OK, finite G_full/r, primal-feasible); an LP failure is NOT
    # "feasibility confirmed" (the freeze-8/9 lesson applies here too).
    ok = all(r["status"] == "OK" and r["G_full_per_r"] != ""
             and np.isfinite(float(r["G_full_per_r"]))
             and r["primal_feasible"] for r in results6)
    terminated = all(r["wall_s"] < 3600 for r in results6)
    gate_untouched = (v5.GFULL_TOL_PER_R == 1e-2)
    check("dev_expanded_class_feasibility_runtime", ok and terminated
          and gate_untouched,
          "certified finite bounds (walls: %s s); gate bar 1e-2 UNCHANGED; "
          "G_full/r: %s" % ([r["wall_s"] for r in results6],
                            [r["G_full_per_r"] for r in results6]))

    # ---- E7: 52+972 accounting + ridge disclosures ------------------------ #
    print("\n== E7: accounting + ridge disclosures ==")
    ok = P.shape[0] == 52
    for arm in m1.ARMS_ALL:
        pat = m1.load_m1_pattern("m1pat:%s:devimg" % arm, 1024, 1024, 0)
        ok &= pat["A"].shape == (1024, N)
        ok &= pat["meta"]["prescan_rows"] == 52
        ok &= pat["meta"]["main_rows"] == 972
    import campaign
    dev_imgs = campaign._images(SIDE, "all", imageset="m1_dev")
    cal = m1.ridge_scat32_calibration(dev_img, 0, imageset="m1_dev",
                                      x_true=dev_imgs[dev_img])
    print("    ridge calibration (%s, s0): " % dev_img)
    for i, nu in enumerate(cal["nu_grid"]):
        print("      nu=%-5g requested=%.4f achieved=%.4f clip=%d "
              "GUARD_CLIPPED=%d" % (nu, cal["requested"][i],
                                    cal["achieved"][i],
                                    cal["clip_applied"][i],
                                    cal["guard_clipped"][i]))
    disc = m1.ridge_disclosure(dev_img, 0, 2000.0, imageset="m1_dev")
    ok &= all(disc[c] != "" for c in ("achieved_mean_load",
                                     "incident_dose_ratio", "load_q95",
                                     "physical_peak", "deployed_sha"))
    check("accounting_52_plus_972_and_ridge_disclosures", ok,
          "5 arms x (52+972); 9-dwell calibration + full disclosure block "
          "(achieved=%.3f at nu=2000, ratio=%.2fx)"
          % (float(cal["achieved"][-1]), disc["incident_dose_ratio"] or -1))

    # ---- E8: manifests + expected cells regenerated ----------------------- #
    print("\n== E8: manifests + expected-cell ledger ==")
    import m1_make_manifests as mm
    blocks = mm.build_blocks()
    shards = mm.pack(blocks, mm.TARGET_SHARDS)
    idx, _l = mm.write_manifests(shards)
    n_exp = mm.write_expected_cells(shards)
    n_cert = sum(1 for s in shards if s["stage"] == "M1_CERT"
                 for _ in s["cells"])
    stale = []
    for s in shards:
        for e in mm.shard_frozen_inputs(s):
            if "designs/" in e["path"] and "scat32_deployed" not in e["path"]:
                stale.append(e["path"])
    ok = (n_cert == 480 and not stale and n_exp == idx["n_expected_rows"])
    check("manifests_and_expected_cells_r17", ok,
          "%d cells (%d cert of 480 required), %d shards; stale design "
          "frozen_inputs: %s" % (n_exp, n_cert, len(idx["default_shards"]),
                                 stale or "none"))

    # ---- E9: no confirmatory data opened ---------------------------------- #
    print("\n== E9: outcome-blindness / confirmatory lock ==")
    guard_fired = False
    try:
        m1.run_cert_cell(m1.cert_cell("m1_glyph_0", 0, 2000.0, 0.60,
                                      imageset="m1"))
    except RuntimeError:
        guard_fired = True
    check("no_confirmatory_data_opened", guard_fired,
          "confirmatory-imageset guard raises pre-freeze; selftest touched "
          "only the synthetic scene and m1_dev (632900+) instances")

    # ---- ledger ------------------------------------------------------------ #
    n_pass = sum(1 for _, ok_, _ in CHECKS if ok_)
    all_ok = n_pass == len(CHECKS)
    lines = ["# M1 REFREEZE CHECKLIST LEDGER (R17 §7 / amendment §E)", "",
             "Outcome-blind selftest, %s. %d/%d items PASS."
             % (time.strftime("%Y-%m-%d %H:%M"), n_pass, len(CHECKS)), "",
             "Architecture: R17 (issue #9). Deployed design = balanced "
             "exact 972-row SCAT32 multiset (sha %s..). Modes: SCAT32-SAFE "
             "0.05 / SCAT32-060 0.60 / RIDGE-SCAT32 rho_R(nu) policy. "
             "Context: SCAT16, LBLOB16 (descriptive). D_cert sha %s.."
             % (sha[:12], sha_cert[:12]), "",
             "Interpretations I1-I5 logged in m1_runner.py docstring.", "",
             "| # | item (R17 Sec 7) | verdict | detail |",
             "|---|---|---|---|"]
    for i, (name, ok_, det) in enumerate(CHECKS, 1):
        lines.append("| %d | %s | %s | %s |"
                     % (i, name, "PASS" if ok_ else "FAIL", det))
    lines += ["", ("**REFREEZE_READY** -- every Sec-7 item passes; one "
                   "immutable revised m1-freeze tag may be created and all "
                   "arms launched together (R17 Sec 5)."
                   if all_ok else "HOLD REMAINS (failures above)."),
              "", "wall time: %.1f s" % (time.time() - t0)]
    with open(LEDGER, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n[refreeze selftest] %d/%d PASS -> %s (%.1fs)"
          % (n_pass, len(CHECKS), LEDGER, time.time() - t0), flush=True)
    if all_ok:
        print("[refreeze selftest] REFREEZE_READY", flush=True)
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
