"""ROUND63 REAL end-to-end integration test (GPT round-3 digest §1.6).

Unlike the solver UNIT tests (test_solvers.py), this drives the ACTUAL production
path — patterns.make_patterns -> images63.build_image_set -> physics.simulate_counts
-> solvers.run_arm — through campaign.run_cell, for all three pattern kinds
{bern50, hadpair, gam4} and the arms {GI, POISSON-LIN, RQL}.

Discipline (digest §1.6): performance differences are NOT PASS gates — they are
printed as a report only. The PASS gates are STRUCTURAL:
  * finite main metrics (PSNR/SSIM) on every row;
  * finite radiometric NRMSE on every PHYSICAL arm row (RQL is a physical arm —
    this row also validates the campaign PHYSICAL_ARMS fix);
  * correct exposure accounting (hadpair -> 2*M physical optical rows / exposures,
    bern50 & gam4 -> M);
  * the REAL meta-key contract: patterns.make_patterns emits the pair table under
    'pair_indices', and hadamard_pair_combine consumes exactly that. The OLD unit
    test hand-built a {'pairs': ...} meta, which MASKED a KeyError bug on the real
    path (linear arm on real hadpair meta); this test uses the REAL meta so the
    bug cannot hide again.
  * complete row schema (exact key set) on every row.

~2 min: the default select_rule='discrepancy' now routes the iterative arms
(POISSON-LIN, RQL) through select_eta.select_eta_and_fit (spec D2 §4), so this
test also covers the production lam_TV selection path on all three pattern kinds
(including hadpair fold co-location). Run:
  D:/Anacondar/anaconda3/python.exe code/round63/test_integration.py
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

import patterns
import physics
import solvers
from campaign import run_cell, PHYSICAL_ARMS

SIDE = 32
N = SIDE * SIDE
KINDS = ("bern50", "hadpair", "gam4")
ARMS = ["GI", "POISSON-LIN", "RQL"]
M = 400
NU = 300.0
TAU = 50e-9
RHO = 0.5
SEED = 0
FISTA = 60

EXPECTED_KEYS = {
    "side", "pattern", "rho_bar", "nu", "M", "seed", "image", "arm",
    "PSNR", "SSIM", "LPIPS", "rad_nrmse", "flux_dev", "lam_tv",
    "mean_counts", "optical_time_s", "dark_frac", "tau_err", "runtime_s",
    # D2 §4 selection-routing columns (campaign.run_cell now emits these):
    "select_runtime_s", "MODEL_FAIL", "eta_star", "PSNR_rad",
    # F1 cell-level adequacy audit columns (RQL rows only; '' elsewhere):
    "gof_status", "gof_p", "leak_suspect",
}

PASS = []


def report(name, ok, detail=""):
    PASS.append(bool(ok))
    print("[%s] %s %s" % ("PASS" if ok else "FAIL", name, detail), flush=True)


# ----------------------------------------------------------------------
def check_meta_contract():
    """The REAL meta-key contract, exercised on the REAL meta (not a hand-built
    {'pairs': ...}). Proves 'pair_indices' from make_patterns flows through
    hadamard_pair_combine, and that the legacy 'pairs' key still round-trips."""
    n = 64
    pat = patterns.make_patterns("hadpair", 32, n, SEED)
    meta, A = pat["meta"], pat["A"]
    # the REAL key contract: make_patterns emits 'pair_indices', NOT 'pairs'
    key_ok = ("pair_indices" in meta) and ("pairs" not in meta)
    pairs = np.asarray(meta["pair_indices"], dtype=np.int64)
    xv = np.abs(np.random.default_rng(0).standard_normal(n)) + 0.1
    b = A @ xv
    # must NOT raise KeyError on the REAL meta (this was the masked bug)
    try:
        A_signed, b_signed = solvers.hadamard_pair_combine(A, b, meta)
        ran = True
    except KeyError:
        A_signed, b_signed, ran = None, None, False
    diff_ok = ran and A_signed.shape == (32, n) and b_signed.shape == (32,) \
        and np.allclose(A_signed, A[pairs[:, 0]] - A[pairs[:, 1]]) \
        and np.allclose(b_signed, b[pairs[:, 0]] - b[pairs[:, 1]])
    # backward compatibility: the legacy 'pairs' key path still works
    As2, bs2 = solvers.hadamard_pair_combine(A, b, {"pairs": pairs})
    legacy_ok = np.allclose(As2, A_signed) and np.allclose(bs2, b_signed)
    ok = key_ok and diff_ok and legacy_ok
    report("meta-key contract (pair_indices -> hadamard_pair_combine)", ok,
           "real_key=%s no_'pairs'=%s ran_on_real_meta=%s signed_ok=%s legacy_ok=%s"
           % ("pair_indices" in meta, "pairs" not in meta, ran, diff_ok, legacy_ok))


def check_kind(kind):
    """Drive campaign.run_cell for one pattern kind on one image, 3 arms."""
    cell = dict(side=SIDE, pattern=kind, rho_bar=RHO, nu=NU, M=M, seed=SEED,
                arms=ARMS, images=["text"], fista_iters=FISTA, tau=TAU)
    t0 = time.time()
    rows = run_cell(cell)                       # <- REAL pipeline (would KeyError
    dt = time.time() - t0                       #    on hadpair+GI if bug present)

    # --- schema + row count
    n_img = 1
    count_ok = len(rows) == n_img * len(ARMS)
    schema_ok = all(set(r.keys()) == EXPECTED_KEYS for r in rows)

    # --- finite metrics (PSNR/SSIM everywhere; rad_nrmse on physical arms)
    finite_ok = True
    phys_rad_ok = True
    for r in rows:
        if not (np.isfinite(r["PSNR"]) and np.isfinite(r["SSIM"])):
            finite_ok = False
        if r["arm"] in PHYSICAL_ARMS:
            v = r["rad_nrmse"]
            if v == "" or not np.isfinite(float(v)):
                phys_rad_ok = False

    # --- exposure accounting from the REAL meta
    T = NU * TAU
    pat = patterns.make_patterns(kind, M, N, SEED)
    meta = pat["meta"]
    if kind == "hadpair":
        exp_rows, exp_expo = 2 * M, 2
    else:
        exp_rows, exp_expo = M, 1
    exp_time = exp_rows * T
    meta_ok = (meta["n_physical_rows"] == exp_rows
               and meta["exposures_per_row"] == exp_expo
               and pat["A"].shape[0] == exp_rows)
    time_ok = all(abs(r["optical_time_s"] - exp_time) <= 1e-9 * exp_time
                  for r in rows)
    expo_ok = meta_ok and time_ok

    # --- RQL specifically must be a physical arm with a finite rad_nrmse
    rql = [r for r in rows if r["arm"] == "RQL"]
    rql_ok = (len(rql) == 1 and rql[0]["rad_nrmse"] != ""
              and np.isfinite(float(rql[0]["rad_nrmse"])))

    ok = (count_ok and schema_ok and finite_ok and phys_rad_ok and expo_ok
          and rql_ok)
    report("run_cell[%s] end-to-end" % kind, ok,
           "rows=%d schema=%s finite=%s phys_rad=%s expo(rows=%d/T=%.4g)=%s "
           "rql_phys=%s (%.1fs)"
           % (len(rows), schema_ok, finite_ok, phys_rad_ok, exp_rows, exp_time,
              expo_ok, rql_ok, dt))

    # --- performance: PRINTED REPORT ONLY (never a PASS gate; digest §1.6)
    for r in sorted(rows, key=lambda z: z["arm"]):
        rad = r["rad_nrmse"] if r["rad_nrmse"] != "" else "n/a"
        print("      %-11s PSNR=%6.2f SSIM=%.4f rad_nrmse=%s mean_counts=%.1f"
              % (r["arm"], r["PSNR"], r["SSIM"], rad, r["mean_counts"]),
              flush=True)


def main():
    t0 = time.time()
    np.seterr(over="ignore", invalid="ignore")
    check_meta_contract()
    for kind in KINDS:
        check_kind(kind)
    dt = time.time() - t0
    n_pass = sum(PASS)
    print("\n==== %d/%d integration checks PASS  (%.1fs) ===="
          % (n_pass, len(PASS), dt), flush=True)
    return 0 if n_pass == len(PASS) else 1


if __name__ == "__main__":
    sys.exit(main())
