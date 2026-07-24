#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
27-CELL ORTHOGONAL-ARRAY MC PLAN + COMPUTE-COST FORECAST (R41 sec 4.3 + 4.8).

Enumerates the precommitted 27-cell 3^(4-1) OA subset, chooses the per-(cell,eps) test T_eff from the
corrected analytic lam_mean, sizes the record counts (>=500 H0 + 500 H1 primary; >=250/class
specificity + mismatch), and forecasts total GPU-hours against the FROZEN 6-RTX-4060-GPU-hour ceiling
using MEASURED throughput (FIXED-COV Gaussian generator ~2e6 record-banks/s; FRESH-COV per-bank
~9 ms/bank). Exceeding the ceiling requires stopping and reporting -- never silently trimming.

Writes COST_FORECAST.md (+ .json). Read-only; no MC executed here (this is the pre-run plan).
"""
import json
import os
import time

import numpy as np
import torch

import sealed_common as sc

HERE = os.path.dirname(os.path.abspath(__file__))

# measured throughput (this RTX 4060; see benchmarks in COST_FORECAST) -- conservative
FIXED_RECORD_BANKS_PER_SEC = 8.0e5      # conservative: matches measured on-box rate incl. large-T slowdown
FRESH_MS_PER_BANK = 9.5                 # measured FRESH-COV-OPT per-bank (batched GPU, chunk=48)
GPU_HOUR_CEILING = 6.0
STORAGE_GB_CEILING = 64.0


def bench_throughput():
    """Quick on-box FIXED-COV throughput measurement (record-banks/s) at two T_eff points."""
    cell = sc.setup_cell(**sc.BEST_CELL)
    sc.gen_records(cell, 50, 500, "H0", None, rng=1)
    torch.cuda.synchronize() if sc.DEV == "cuda" else None
    rates = []
    for T in [800, 2000]:
        t0 = time.time()
        sc.gen_records(cell, 200, T, "H0", None, rng=2)
        if sc.DEV == "cuda":
            torch.cuda.synchronize()
        rates.append(200 * T / (time.time() - t0))
    return float(np.mean(rates))


def plan_primary():
    """Per-(cell,eps) test T_eff and record-bank counts for the 27-cell OA primary D1/D2 grid."""
    oa = sc.oa_grid_27()
    rows = []
    total_rb = 0
    for geo in oa:
        cell = sc.setup_cell(**geo, want_mc=False)
        lam = cell["lam_mean"]
        xn2 = cell["xnorm2"]
        for eps in sc.EPS_LEVELS:
            Tdet = sc.DP_STRONG2 / (lam * eps * eps * xn2)
            T_eff = int(min(max(round(Tdet), 32), sc.T_CAP))
            nrec = 500 + 500                                   # H0 + H1
            rb = nrec * T_eff
            total_rb += rb
            rows.append(dict(**geo, eps=eps, T_det_analytic=round(Tdet, 0), T_eff=T_eff,
                             n_records=nrec, record_banks=rb))
    return rows, total_rb


def plan_specificity():
    """6 event classes x >=250 records at a representative T_eff, on 3 base cells (best/mid/floor)."""
    cells = [sc.BEST_CELL, sc.MID_CELL, sc.FLOOR_CELL]
    total_rb = 0
    rows = []
    for geo in cells:
        cell = sc.setup_cell(**geo, want_mc=False)
        # specificity T_eff: enough for the beyond-band target d'~7 at eps=2% and to certify the
        # intended amplitude/lag scores (>=2048 banks certifies the lag class near d'=5)
        T_eff = int(min(2048, sc.T_CAP))
        n = 6 * 250                                            # six classes
        rb = n * T_eff
        total_rb += rb
        rows.append(dict(**geo, T_eff=T_eff, classes=6, records_per_class=250, record_banks=rb))
    return rows, total_rb


def plan_mismatch():
    """D5 axes x levels x >=250 records at 2 base cells (best/mid)."""
    axes_levels = sum(len(a["levels"]) for a in __import__("sealed_banks").MISMATCH_AXES)
    cells = [sc.BEST_CELL, sc.MID_CELL]
    total_rb = 0
    rows = []
    for geo in cells:
        cell = sc.setup_cell(**geo, want_mc=False)
        T_eff = 800
        n = axes_levels * (250 + 250)                         # H0 + H1 per axis-level
        rb = n * T_eff
        total_rb += rb
        rows.append(dict(**geo, T_eff=T_eff, axis_levels=axes_levels, records_per_level=250,
                         record_banks=rb))
    return rows, total_rb


def plan_fresh():
    """D4 FRESH-COV-OPT: analytic latency over the 27 OA cells (cheap) + MC validation at 6 points."""
    # analytic cost: 27 cells x 8 code draws x setup_cell(~0.27 s)
    analytic_sec = 27 * 8 * 0.27
    # MC validation: best/mid/floor x eps in {0.02, 0.05} x 40 rec x 2 hyp
    pts = []
    mc_bank_streams = 0
    for geo in [sc.BEST_CELL, sc.MID_CELL, sc.FLOOR_CELL]:
        cell = sc.setup_cell(**geo, want_mc=False)
        for eps in [0.02, 0.05]:
            Tdet = sc.DP_STRONG2 / (cell["lam_mean"] * eps * eps * cell["xnorm2"])
            T_eff = int(min(max(round(Tdet), 64), sc.T_CAP))
            streams = 40 * 2 * T_eff
            mc_bank_streams += streams
            pts.append(dict(**geo, eps=eps, T_eff=T_eff, records=40, bank_streams=streams))
    mc_sec = mc_bank_streams * FRESH_MS_PER_BANK / 1000.0
    return pts, analytic_sec, mc_sec, mc_bank_streams


def storage_forecast(total_records):
    """Records are reduced to scalar statistics on the fly (no per-bank frames stored). Only the
    per-record score scalars + summary JSON are persisted -> tiny. Report worst-case if M x M sample
    covariances were retained for audit (they are NOT; scalars only)."""
    scalar_gb = total_records * 8 / 1e9
    audit_cov_gb = total_records * (sc.M * sc.M) * 8 / 1e9      # if full covariances were kept (not)
    return dict(scalar_only_gb=round(scalar_gb, 4), if_full_cov_retained_gb=round(audit_cov_gb, 2))


def main():
    print("=== 27-cell OA MC plan + cost forecast ===", flush=True)
    fixed_rate = bench_throughput()
    print(f"measured FIXED-COV throughput: {fixed_rate:.2e} record-banks/s "
          f"(planning at {FIXED_RECORD_BANKS_PER_SEC:.1e}, conservative)", flush=True)

    prim_rows, prim_rb = plan_primary()
    spec_rows, spec_rb = plan_specificity()
    mis_rows, mis_rb = plan_mismatch()
    fresh_pts, fresh_an_sec, fresh_mc_sec, fresh_streams = plan_fresh()

    prim_sec = prim_rb / FIXED_RECORD_BANKS_PER_SEC
    spec_sec = spec_rb / FIXED_RECORD_BANKS_PER_SEC
    mis_sec = mis_rb / FIXED_RECORD_BANKS_PER_SEC
    d0_sec = 300.0        # two-engine Fisher + physical Poisson MC shot check (few cells)
    overhead_sec = 27 * 4 * 0.27 + 300     # per-cell setup builds + misc

    total_sec = prim_sec + spec_sec + mis_sec + fresh_an_sec + fresh_mc_sec + d0_sec + overhead_sec
    total_hr = total_sec / 3600.0

    total_records = (sum(r["n_records"] for r in prim_rows)
                     + sum(r["classes"] * r["records_per_class"] for r in spec_rows)
                     + sum(r["axis_levels"] * 2 * r["records_per_level"] for r in mis_rows)
                     + sum(p["records"] * 2 for p in fresh_pts))
    stor = storage_forecast(total_records)

    forecast = dict(
        gpu_hour_ceiling=GPU_HOUR_CEILING, storage_gb_ceiling=STORAGE_GB_CEILING,
        measured_fixed_rate_record_banks_per_sec=round(fixed_rate, 0),
        planning_fixed_rate=FIXED_RECORD_BANKS_PER_SEC, fresh_ms_per_bank=FRESH_MS_PER_BANK,
        blocks=dict(
            D0_mechanism=dict(sec=round(d0_sec, 0)),
            D1_D2_primary=dict(cells=27, eps=len(sc.EPS_LEVELS), record_banks=int(prim_rb),
                               sec=round(prim_sec, 0)),
            D3_specificity=dict(base_cells=len(spec_rows), record_banks=int(spec_rb),
                                sec=round(spec_sec, 0)),
            D5_mismatch=dict(base_cells=len(mis_rows), record_banks=int(mis_rb),
                             sec=round(mis_sec, 0)),
            D4_fresh=dict(analytic_sec=round(fresh_an_sec, 0), mc_sec=round(fresh_mc_sec, 0),
                          mc_bank_streams=int(fresh_streams), points=len(fresh_pts)),
            setup_overhead=dict(sec=round(overhead_sec, 0)),
        ),
        total_sec=round(total_sec, 0), total_gpu_hours=round(total_hr, 3),
        within_ceiling=bool(total_hr <= GPU_HOUR_CEILING),
        total_records=int(total_records), storage=stor,
        primary_grid=prim_rows, specificity_plan=spec_rows, mismatch_plan=mis_rows,
        fresh_validation_points=fresh_pts,
    )
    json.dump(forecast, open(os.path.join(HERE, "COST_FORECAST.json"), "w"), indent=2)
    print(f"\nTOTAL: {total_hr:.2f} GPU-hours  (ceiling {GPU_HOUR_CEILING})  -> "
          f"{'WITHIN' if forecast['within_ceiling'] else 'OVER'} ceiling", flush=True)
    for b, v in forecast["blocks"].items():
        print(f"  {b:18s} {v}", flush=True)
    print(f"records total={total_records}  storage(scalar)={stor['scalar_only_gb']} GB "
          f"(< {STORAGE_GB_CEILING} GB ceiling)", flush=True)
    return forecast


if __name__ == "__main__":
    main()
