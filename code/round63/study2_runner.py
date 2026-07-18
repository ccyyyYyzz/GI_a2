"""ROUND63 STUDY-2 runner — the contrast-dead-time phase-diagram campaign
(frozen spec docs/ROUND63_GPT_ROUND8_RULING.md).

Study 2 tests, at FIXED mean detector load and FIXED incident-photon budget,
whether raising illumination-pattern contrast by lowering spatial occupancy
opens a photon-limited regime where dead-time-aware high-flux operation reduces
acquisition time. Geometry is frozen at 32x32, M=n=1024, tau=50 ns, sigma_b=0,
dark=0, nine dwell points nu in {5,10,20,50,100,200,500,1000,2000}, two operating
points rho_bar in {0.05, 0.60}, on the DETAIL-32 cohort (six families x four
confirmatory instances = 24). Every cell runs through campaign.run_cell — the
same production path as Study 1 — with the frozen analytic-lambda rule (C0 from
select_eta.frozen_C0()).

Modes (each writes a resume-safe CSV under results/round63_study2/):
  --primary     ONLY k=16 carries the positive gate (ruling §5). 24 images x 5
                seeds x both rho x 9 nu; arms RQL, POISSON-LIN, SAT-POISSON,
                PRECORRECT, GI. RQL descriptive audit ON (as Study 1).
  --robustness  k=32 replication arm, RQL only, same grid, NO gate (ruling §8).
  --controls    k=512 (dense) and k=1 (raster) demonstration controls, RQL only,
                same grid, NO gate (ruling §9).
  --fluxmap     descriptive flux map (ruling §12): one frozen representative per
                family (detail32_<fam>_0), rho in {0.05,0.1,0.3,0.6,1,2},
                nu in {20,200,2000}, 3 seeds, ALL FOUR k, RQL.
  --separation  mechanism-separation (ruling §13): 6 reps, nu {20,200,2000}, 3
                seeds, k {512,16}, both rho, TWO flux modes — normalized
                (A=(n/k)B, load fixed) vs unnormalized (A=B, detector rate drops
                with k), RQL.
  --scattering  dynamic-scattering secondary (ruling §14): k=16 geometry, 6 reps,
                rho {0.05,0.6}, nu {20,200,2000}, 3 seeds, AR(1) lognormal
                CV(alpha)=0.20 (rho_t=0.95) applied to lam BEFORE dead time, RQL.
  --analyze     reuse pilot_s1.endpoint_analysis (cohort='detail') on
                primary_rows.csv and write STUDY2_SUMMARY.md + study2_summary.json.

The RQL descriptive audit is OFF in every mode except --primary (ruling: audit on
as in Study 1). The 24 confirmatory instances are reconstructed ONLY inside the
frozen modes above; no dev instance may alter k, families, rho, nu, c, endpoints
or gates (ruling §15).

Run (from the repo root, with the py311 env):
  python code/round63/study2_runner.py --primary
  python code/round63/study2_runner.py --analyze
"""
import argparse
import csv
import json
import os
import sys
import time
import traceback
from types import SimpleNamespace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
ROOT = os.path.dirname(os.path.dirname(HERE))

from campaign import run_cell
from detail24 import _conf32_table, _dev32_table, FAMILIES
from select_eta import frozen_C0, C0_FILE
import pilot_s1

# ---- frozen Study-2 geometry (ruling §1) --------------------------------- #
SIDE = 32
M = SIDE * SIDE                                    # M = n = 1024 (square)
TAU = 50e-9
PATTERN = "sparsek"
NU_FULL = [5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0, 2000.0]
RHO_BOTH = [0.05, 0.60]
SEEDS5 = [0, 1, 2, 3, 4]
SEEDS3 = [0, 1, 2]
ARMS_PRIMARY = ["RQL", "POISSON-LIN", "SAT-POISSON", "PRECORRECT", "GI"]

# occupancy ladder (ruling §3): 50% / 3.1% / 1.6% / 0.1%
K_ALL = [512, 32, 16, 1]
K_PRIMARY = 16
K_ROBUST = 32
K_CONTROLS = [512, 1]
K_SEPARATION = [512, 16]

# descriptive flux-map / separation / scattering sub-grids (ruling §12-§14)
RHO_FLUX = [0.05, 0.1, 0.3, 0.6, 1.0, 2.0]
NU_MAP = [20.0, 200.0, 2000.0]
ALPHA_CV = 0.20
ALPHA_RHO_T = 0.95

# frozen image handles
CONF32 = [name for name, _f, _s in _conf32_table()]     # 24 confirmatory
REPS6 = ["detail32_%s_0" % fam for fam in FAMILIES]      # one rep per family

# ---- CSV schema: run_cell emission + runner-injected provenance ---------- #
MECH_COLS = ["cnr", "C_u", "Gamma", "S_det", "S_inc", "k_occupancy"]
PROV_COLS = ["k", "flux_mode", "alpha_cv", "study2_mode", "cell_id"]
STUDY2_FIELDNAMES = list(pilot_s1.FIELDNAMES) + MECH_COLS + PROV_COLS
# resume key: unique per (image, arm) row WITHIN one mode's CSV
RESUME_COLS = ["k", "flux_mode", "rho_bar", "nu", "M", "seed", "image", "arm"]

MODES = ("primary", "robustness", "controls", "fluxmap", "separation",
         "scattering")


# ---- cell construction ---------------------------------------------------- #
def _cell(mode, k, rho, nu, seed, arms, images, C0, flux_mode="normalized",
          alpha_cv=0.0, audit=False):
    cid = "%s_k%d_%s_r%.3g_nu%g_s%d" % (mode, k, flux_mode, rho, nu, seed)
    if alpha_cv > 0.0:
        cid += "_acv%.2f" % alpha_cv
    return dict(
        side=SIDE, pattern=PATTERN, k=int(k), rho_bar=float(rho), nu=float(nu),
        M=M, seed=int(seed), arms=list(arms), images=list(images),
        imageset="detail32", tau=TAU, sigma_b=0.0, fista_iters=200,
        select_iter=60, select_rule="discrepancy", use_lpips=False,
        flux_mode=flux_mode, alpha_cv=float(alpha_cv), alpha_rho_t=ALPHA_RHO_T,
        audit=bool(audit), C0=C0, cell_id=cid, study2_mode=mode)


def cells_primary(C0):
    """k=16 confirmatory arm (ruling §5). RQL audit ON as in Study 1."""
    return [_cell("primary", K_PRIMARY, rho, nu, seed, ARMS_PRIMARY, CONF32, C0,
                  audit=True)
            for rho in RHO_BOTH for nu in NU_FULL for seed in SEEDS5]


def cells_robustness(C0):
    """k=32 replication arm, RQL only, no gate (ruling §8)."""
    return [_cell("robustness", K_ROBUST, rho, nu, seed, ["RQL"], CONF32, C0)
            for rho in RHO_BOTH for nu in NU_FULL for seed in SEEDS5]


def cells_controls(C0):
    """k in {512,1} demonstration controls, RQL only, no gate (ruling §9)."""
    return [_cell("controls", k, rho, nu, seed, ["RQL"], CONF32, C0)
            for k in K_CONTROLS for rho in RHO_BOTH for nu in NU_FULL
            for seed in SEEDS5]


def cells_fluxmap(C0):
    """Descriptive flux map (ruling §12): 6 reps, all four k, RQL."""
    return [_cell("fluxmap", k, rho, nu, seed, ["RQL"], REPS6, C0)
            for k in K_ALL for rho in RHO_FLUX for nu in NU_MAP
            for seed in SEEDS3]


def cells_separation(C0):
    """Mechanism separation (ruling §13): k {512,16} x two flux modes x both
    rho x nu {20,200,2000} x 3 seeds, 6 reps, RQL. (rho unspecified by the
    ruling for this run -> both frozen operating points, documented.)"""
    return [_cell("separation", k, rho, nu, seed, ["RQL"], REPS6, C0,
                  flux_mode=fm)
            for k in K_SEPARATION
            for fm in ("normalized", "unnormalized")
            for rho in RHO_BOTH for nu in NU_MAP for seed in SEEDS3]


def cells_scattering(C0):
    """Dynamic-scattering secondary (ruling §14): k=16, 6 reps, both rho,
    nu {20,200,2000}, 3 seeds, AR(1) lognormal CV=0.20, RQL."""
    return [_cell("scattering", K_PRIMARY, rho, nu, seed, ["RQL"], REPS6, C0,
                  alpha_cv=ALPHA_CV)
            for rho in RHO_BOTH for nu in NU_MAP for seed in SEEDS3]


_CELL_BUILDERS = {
    "primary": cells_primary, "robustness": cells_robustness,
    "controls": cells_controls, "fluxmap": cells_fluxmap,
    "separation": cells_separation, "scattering": cells_scattering,
}


# ---- resumable sweep ------------------------------------------------------ #
def _resume_key_row(r):
    return tuple(str(r[c]) for c in RESUME_COLS)


def _expected_row_keys(cell):
    """Resume keys a cell will emit, mirroring run_cell's int()/float() coercions
    and the runner-injected k/flux_mode so the strings match _resume_key_row."""
    kk = str(int(cell["k"]))
    fm = str(cell["flux_mode"])
    rho = str(float(cell["rho_bar"]))
    nu = str(float(cell["nu"]))
    Mm = str(int(cell["M"]))
    sd = str(int(cell["seed"]))
    return [(kk, fm, rho, nu, Mm, sd, img, arm)
            for img in cell["images"] for arm in cell["arms"]]


def _inject(row, cell, mode):
    row["k"] = int(cell["k"])
    row["flux_mode"] = cell["flux_mode"]
    row["alpha_cv"] = float(cell["alpha_cv"])
    row["study2_mode"] = mode
    row["cell_id"] = cell["cell_id"]
    return row


def run_mode(mode, cells, out_dir):
    """Incremental, resumable sweep for one mode -> results/<mode>_rows.csv."""
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "%s_rows.csv" % mode)
    errlog = os.path.join(out_dir, "%s_errors.log" % mode)

    done = set()
    header_present = os.path.exists(path) and os.path.getsize(path) > 0
    if header_present:
        with open(path, newline="") as f:
            rd = csv.DictReader(f)
            existing = rd.fieldnames
            for r in rd:
                done.add(_resume_key_row(r))
        if existing is not None and list(existing) != STUDY2_FIELDNAMES:
            raise SystemExit(
                "[study2:%s] existing %s has an incompatible header\n  found:  "
                "%s\n  expect: %s\nremove/rename it to re-run."
                % (mode, path, existing, STUDY2_FIELDNAMES))

    n_cells = len(cells)
    t0 = time.time()
    with open(path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=STUDY2_FIELDNAMES,
                                extrasaction="ignore")
        if not header_present:
            writer.writeheader()
            f.flush()
        for ci, cell in enumerate(cells):
            keys = _expected_row_keys(cell)
            if all(kk in done for kk in keys):
                print("[study2:%s] SKIP %d/%d %s (%d rows present)"
                      % (mode, ci + 1, n_cells, cell["cell_id"], len(keys)),
                      flush=True)
                continue
            tc = time.time()
            try:
                rows = run_cell(cell)
            except Exception:
                tb = traceback.format_exc()
                with open(errlog, "a") as ef:
                    ef.write("[%s] %s\n%s\n"
                             % (time.strftime("%Y-%m-%d %H:%M:%S"),
                                cell["cell_id"], tb))
                print("[study2:%s] ERROR %d/%d %s -> logged to %s, continuing"
                      % (mode, ci + 1, n_cells, cell["cell_id"], errlog),
                      flush=True)
                continue
            n_new = 0
            for r in rows:
                _inject(r, cell, mode)
                kk = _resume_key_row(r)
                if kk in done:
                    continue
                writer.writerow(r)
                done.add(kk)
                n_new += 1
            f.flush()
            os.fsync(f.fileno())
            print("[study2:%s] cell %d/%d %s -> +%d rows (%.1fs, total %.0fs)"
                  % (mode, ci + 1, n_cells, cell["cell_id"], n_new,
                     time.time() - tc, time.time() - t0), flush=True)
    print("[study2:%s] SWEEP DONE wall=%.0fs rows=%s"
          % (mode, time.time() - t0, path), flush=True)
    return path


# ---- analysis (reuse pilot_s1 endpoint, cohort='detail') ------------------ #
def _analyze_cfg(boot_B):
    """A pilot_s1-compatible cfg for endpoint_analysis / summary formatting."""
    return SimpleNamespace(
        cohort="detail", boot_B=int(boot_B), rho_list=list(RHO_BOTH),
        images=list(CONF32), seeds=list(SEEDS5), arms=list(ARMS_PRIMARY),
        side=SIDE, pattern=PATTERN, M=M, nu_list=list(NU_FULL), tau=TAU,
        select_iter=60, fista_iters=200, C0=frozen_C0(), mode="study2_primary")


def analyze(out_dir, boot_B):
    path = os.path.join(out_dir, "primary_rows.csv")
    if not os.path.exists(path):
        raise SystemExit("[study2:analyze] no primary_rows.csv at %s (run "
                         "--primary first)" % path)
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SystemExit("[study2:analyze] %s is empty" % path)

    cfg = _analyze_cfg(boot_B)
    cen = pilot_s1.endpoint_analysis(rows, cfg)
    mft = pilot_s1.audit_table(rows)
    rtc = pilot_s1.runtime_calibration(rows)

    summary = {"config": {"cohort": "detail", "boot_B": int(boot_B),
                          "n_confirmatory_images": len(CONF32),
                          "k_primary": K_PRIMARY, "rho": RHO_BOTH,
                          "nu": NU_FULL, "seeds": SEEDS5},
               "n_rows": len(rows), "endpoint_analysis": cen,
               "descriptive_audit": mft, "runtime_calibration": rtc}
    with open(os.path.join(out_dir, "study2_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=pilot_s1._json_default)

    md = os.path.join(out_dir, "STUDY2_SUMMARY.md")
    _write_summary_md(cen, len(rows), md)

    pg = cen.get("primary_gate", {})
    fb = cen.get("fixed_budget_quality_gain", {})
    bt = cen.get("bootstrap", {})
    print("[study2:analyze] rows=%d  k=%d  images=%d  rho %g->%g"
          % (len(rows), K_PRIMARY, cen.get("n_images", 0),
             cen.get("rho_safe", float("nan")),
             cen.get("rho_fast", float("nan"))), flush=True)
    print("       primary Q90: median S_gate=%s  n(S_gate>1)=%d/%d  boot "
          "LB2.5=%s  DETAIL_SPEED_PASS=%s"
          % (pilot_s1._fmt_S(cen.get("S_median")), cen.get("n_Sgate_gt1", 0),
             cen.get("n_images", 0), pilot_s1._fmt_S(bt.get("S_gate_LB2.5")),
             pg.get("DETAIL_SPEED_PASS")), flush=True)
    print("       fixed-budget: median DeltaQ=%s dB  n(DeltaQ>0)=%d/%d  boot "
          "LB2.5=%s  pass=%s"
          % (pilot_s1._fmt_dB(fb.get("median_DeltaQ_dB")),
             fb.get("n_DeltaQ_pos", 0), cen.get("n_images", 0),
             pilot_s1._fmt_dB(bt.get("DeltaQ_LB2.5")), fb.get("pass")),
          flush=True)
    print("[study2:analyze] wrote %s , study2_summary.json" % md, flush=True)
    return summary


def _write_summary_md(cen, n_rows, path):
    pg = cen.get("primary_gate", {})
    fb = cen.get("fixed_budget_quality_gain", {})
    bt = cen.get("bootstrap", {})
    L = []
    L.append("# ROUND63 STUDY-2 SUMMARY — contrast-dead-time phase diagram")
    L.append("")
    L.append("Primary confirmatory arm k=%d on the DETAIL-32 cohort (six "
             "families x four instances = 24), reconstructed with RQL through "
             "the frozen production path; round-7 Q90 acquisition-speed "
             "endpoint, family-stratified nested bootstrap (ruling §5/§6)."
             % K_PRIMARY)
    L.append("")
    L.append("- rows: %d   images: %d   rho_safe=%g   rho_fast=%g   "
             "nu_budget=%s" % (n_rows, cen.get("n_images", 0),
                               cen.get("rho_safe", float("nan")),
                               cen.get("rho_fast", float("nan")),
                               cen.get("nu_budget")))
    L.append("")
    L.append("## Primary gate (ruling §6) — DETAIL_SPEED_PASS = %s"
             % pg.get("DETAIL_SPEED_PASS"))
    L.append("")
    L.append("- median S_gate = %s  (>=3: %s)"
             % (pilot_s1._fmt_S(cen.get("S_median")), pg.get("median>=3")))
    L.append("- bootstrap 2.5%% LB of median S_gate = %s  (>1: %s)"
             % (pilot_s1._fmt_S(bt.get("S_gate_LB2.5")), pg.get("LB>1")))
    L.append("- images with S_gate>1 = %d/%d  (>=18: %s)"
             % (cen.get("n_Sgate_gt1", 0), cen.get("n_images", 0),
                pg.get("n_gt1>=18")))
    L.append("")
    L.append("## Fixed-budget quality gain (ruling §7, secondary)")
    L.append("")
    L.append("- median DeltaQ = %s dB  (>=1.0: %s)"
             % (pilot_s1._fmt_dB(fb.get("median_DeltaQ_dB")),
                fb.get("median>=1.0")))
    L.append("- bootstrap 2.5%% LB of median DeltaQ = %s  (>0: %s)"
             % (pilot_s1._fmt_dB(bt.get("DeltaQ_LB2.5")), fb.get("LB>0")))
    L.append("- images with DeltaQ>0 = %d/%d  (>=18: %s)   pass = %s"
             % (fb.get("n_DeltaQ_pos", 0), cen.get("n_images", 0),
                fb.get("n_pos>=18"), fb.get("pass")))
    L.append("")
    sc = cen.get("status_counts", {})
    L.append("Censoring status counts: %s"
             % (", ".join("%s=%d" % (k, v) for k, v in sorted(sc.items()))
                or "none"))
    L.append("")
    if cen.get("per_image"):
        L.append("## Per-image endpoints (RQL, Q90)")
        L.append("")
        L.append("| image | family | status | S_gate | T_safe | T_fast | "
                 "R_j (dB) | Q90_j (dB) | DeltaQ_budget (dB) | seeds |")
        L.append("|---|---|---|---|---|---|---|---|---|---|")
        for d in cen["per_image"]:
            L.append("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %d |"
                     % (d["image"], d.get("family", "--"), d["status"],
                        pilot_s1._fmt_S(d["S_gate"]), pilot_s1._fmt_T(d["T_s"]),
                        pilot_s1._fmt_T(d["T_f"]), pilot_s1._fmt_dB(d.get("R_j")),
                        pilot_s1._fmt_dB(d.get("Q90_j")),
                        pilot_s1._fmt_dB(d.get("DeltaQ_budget")),
                        d.get("n_seeds", 0)))
        L.append("")
    if cen.get("warnings"):
        L.append("**Analysis-completeness warnings:**")
        for w in cen["warnings"]:
            L.append("- %s" % w)
        L.append("")
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")


# ---- driver --------------------------------------------------------------- #
def parse_args(argv):
    ap = argparse.ArgumentParser(
        description="ROUND63 Study-2 runner (frozen round-8 ruling).")
    for m in MODES:
        ap.add_argument("--%s" % m, action="store_true",
                        help="run the %s sweep" % m)
    ap.add_argument("--analyze", action="store_true",
                    help="reuse pilot_s1.endpoint_analysis (cohort=detail) on "
                         "primary_rows.csv; write STUDY2_SUMMARY.md")
    ap.add_argument("--out", type=str, default=None,
                    help="output dir (default results/round63_study2)")
    ap.add_argument("--boot-B", type=int, default=pilot_s1.BOOT_B_DEFAULT,
                    dest="boot_B", help="nested bootstrap replicates for "
                    "--analyze (default %d)" % pilot_s1.BOOT_B_DEFAULT)
    return ap.parse_args(argv)


def main(argv=None):
    a = parse_args(sys.argv[1:] if argv is None else argv)
    out = a.out or os.path.join(ROOT, "results", "round63_study2")
    selected = [m for m in MODES if getattr(a, m)]

    if not selected and not a.analyze:
        raise SystemExit("[study2] choose at least one of %s or --analyze"
                         % ", ".join("--%s" % m for m in MODES))

    if selected:
        C0 = frozen_C0()
        if C0 is None:
            raise SystemExit(
                "[study2] no frozen concentration threshold at %s — run "
                "`pilot_s1.py --pass-a --freeze-c0` first (ruling §5: the "
                "analytic lambda needs the development-calibrated C0)." % C0_FILE)
        print("[study2] out=%s  C0=%s  modes=%s"
              % (out, ("inf" if np.isinf(C0) else "%g" % C0),
                 ", ".join(selected)), flush=True)
        for mode in selected:
            cells = _CELL_BUILDERS[mode](C0)
            print("[study2:%s] %d cells" % (mode, len(cells)), flush=True)
            run_mode(mode, cells, out)

    if a.analyze:
        analyze(out, a.boot_B)
    print("[study2] DONE", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
