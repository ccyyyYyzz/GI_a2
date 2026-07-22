#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
DLGI CONFIRMATORY CAMPAIGN -- calibration-bank launcher + shard plan (R34 sec.1-2).

The calibration bank constructs c_0.95(eta) ONLY (R34 sec.1).  Primary grid = 27
cells: t_c/Delta in {16,32,64} x CV in {0.05,0.15,0.40} x photon scale {0.5,1,2}
(R34 sec.2.1/2.2), >= 5000 full-pipeline replicates per cell across the sealed
CALIBRATION bank scenes.  Each replicate regenerates the gain path + Poisson noise
and RERUNS joint_dual_ledger (scene NOT held fixed); W(Y; eta_0=true) feeds the
conservative MC order-statistic c_0.95.

This file provides:
  * run_cell()   -- compute c_0.95 for one cell (the reusable unit of work).
  * run_shard()  -- a subset of cells (the Colab 5-route split).
  * merge()      -- combine partial tables into the frozen 27-cell critical-value
                    table + CritSurface, written as CALIBRATION_TABLE.json.
  * plan()       -- cost estimate from the smoke timing + shard manifests (this is
                    what tonight's scaffold emits; the RUN is a separate GO).

DISCIPLINE.  The calibration table is committed BEFORE the confirmatory bank is
opened; after that commit no tuning/regridding/estimator-switch/interval-repair is
permitted (R34 sec.2.1).  Run modes are explicit CLI flags; the default is --plan
(no compute).  CPU.
"""
import os
import sys
import json
import math
import time
import platform
from datetime import datetime, timezone

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dlgi_common as K
import dlgi_neyman as NY
import dl_common as C

OUT = os.path.join(K.CAMPAIGN, "calibration_plan")
TABLE = os.path.join(K.CAMPAIGN, "CALIBRATION_TABLE.json")

N_REP_PER_CELL = 5000        # R34 sec.2.1 minimum
N_ROUTES = 5                 # Colab 5-route split
SEC_PER_REP = 1.26          # measured end-to-end (sim + joint_dual_ledger + Whittle grid)


# ------------------------------------------------------------ unit of work
def _alloc(n, nsc):
    base = n // nsc
    return [base + (1 if k < n % nsc else 0) for k in range(nsc)]


def run_cell(tc, cv, ph, n_rep=N_REP_PER_CELL, nproc=None):
    """c_0.95 for one calibration cell (full-pipeline; scene re-fit each replicate)."""
    sig = C.sigma_l_of_cv(cv); ckey = K.cell_key(tc, cv, ph)
    scenes, _ = K.load_bank("calibration")
    ids = list(scenes)
    reps = {sid: list(range(a)) for sid, a in zip(ids, _alloc(n_rep, len(ids)))}
    t0 = time.time()
    res = NY._run_split("calibration", ckey, scenes, reps, tc, sig, ph, nproc)
    W = np.array([NY.W_at_true(d["W_grid"].astype(np.float64), tc, cv) for d in res])
    c95, n_used = NY.conservative_critical(W, NY.LEVEL)
    return dict(ckey=ckey, tc=tc, cv=cv, ph=ph, c_0_95=c95, n=n_used,
                W_median=float(np.median(W)), W_p50_p90_p99=[float(np.percentile(W, q))
                                                             for q in (50, 90, 99)],
                tc_hat_med=float(np.median([d["tc_hat"] for d in res])),
                cv_hat_med=float(np.median([d["cv_hat"] for d in res])),
                runtime_s=time.time() - t0)


def cells_for_route(route, cells):
    return [cells[i] for i in range(len(cells)) if i % N_ROUTES == route]


def run_shard(route, n_rep=N_REP_PER_CELL, nproc=None):
    """Run all cells assigned to one Colab route; write a partial table."""
    cells = K.primary_cells()
    mine = cells_for_route(route, cells)
    out = {}
    for (tc, cv, ph) in mine:
        r = run_cell(tc, cv, ph, n_rep, nproc)
        out[r["ckey"]] = r
        print(f"[shard {route}] {r['ckey']}: c_0.95={r['c_0_95']:.3f} "
              f"(n={r['n']}, {r['runtime_s']/60:.1f} min)", flush=True)
    fn = os.path.join(OUT, f"partial_route{route}.json")
    os.makedirs(OUT, exist_ok=True)
    json.dump(dict(route=route, cells=out), open(fn, "w"), indent=2)
    print(f"[shard {route}] wrote {fn}", flush=True)
    return out


def merge():
    """Combine partial route tables (or a single local table) into the frozen
    27-cell critical-value table + CritSurface."""
    parts = {}
    for route in range(N_ROUTES):
        fn = os.path.join(OUT, f"partial_route{route}.json")
        if os.path.exists(fn):
            parts.update(json.load(open(fn))["cells"])
    local_fn = os.path.join(OUT, "partial_local.json")
    if os.path.exists(local_fn):
        parts.update(json.load(open(local_fn))["cells"])
    cells = K.primary_cells()
    missing = [K.cell_key(*c) for c in cells if K.cell_key(*c) not in parts]
    if missing:
        raise RuntimeError("missing calibration cells: %s" % missing)
    # assemble c[tc,cv,ph]
    tcs, cvs, phs = K.TC_PRIMARY, K.CV_PRIMARY, K.PHOTON_PRIMARY
    cvals = np.empty((len(tcs), len(cvs), len(phs)))
    for a, tc in enumerate(tcs):
        for b, cv in enumerate(cvs):
            for d, ph in enumerate(phs):
                cvals[a, b, d] = parts[K.cell_key(tc, cv, ph)]["c_0_95"]
    surf = NY.CritSurface(tcs, cvs, phs, cvals)
    meta = dict(utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
                python=platform.python_version(), numpy=np.__version__,
                n_rep_per_cell=N_REP_PER_CELL, level=NY.LEVEL,
                estimator="joint_dual_ledger (frozen probe, verbatim)",
                whittle="dl_bar4_final (frozen probe, verbatim)",
                status="FROZEN calibration table; commit BEFORE opening confirmatory bank")
    blob = dict(meta=meta, cells=parts, surface=surf.to_dict())
    json.dump(blob, open(TABLE, "w"), indent=2)
    print(f"[merge] wrote frozen calibration table {TABLE} ({len(parts)} cells)", flush=True)
    return TABLE


# ------------------------------------------------------------ cost + plan
GB_PER_WORKER = 0.4          # measured: each worker imports numpy+scipy+frozen probe


def _mem_safe_cores(ncpu):
    """Windows commit-charge guard: each spawned worker re-imports numpy/scipy
    (~0.4 GB).  Cap workers so the pool does not exhaust available RAM+pagefile
    (the DLL-load / '页面文件太小' failure).  Falls back to a conservative 10."""
    try:
        import psutil
        avail = psutil.virtual_memory().available / 1e9
        return max(4, min(ncpu, int((avail - 1.0) / GB_PER_WORKER)))
    except Exception:
        return min(ncpu, 10)


def plan(sec_per_rep=SEC_PER_REP):
    os.makedirs(OUT, exist_ok=True)
    cells = K.primary_cells()
    n_cells = len(cells)
    total_rep = n_cells * N_REP_PER_CELL
    ncpu = os.cpu_count() or 8
    safe = _mem_safe_cores(ncpu)
    single_core_h = total_rep * sec_per_rep / 3600.0
    local_h = single_core_h / ncpu
    local_safe_h = single_core_h / safe

    # Colab 5-route split (assume ~2 usable CPU cores per Colab CPU runtime; the
    # pipeline is CPU-bound numpy, so an L4 GPU does not accelerate it)
    colab_cores = 2
    per_route_rep = max(len(cells_for_route(r, cells)) for r in range(N_ROUTES)) * N_REP_PER_CELL
    colab_route_h = per_route_rep * sec_per_rep / (3600.0 * colab_cores)

    # confirmatory run cost (for reference; separate GO): >=1000 records/cell x 4 arms
    conf_rep = n_cells * 1000
    conf_local_h = conf_rep * sec_per_rep / 3600.0 / ncpu

    routes = []
    for r in range(N_ROUTES):
        mine = cells_for_route(r, cells)
        routes.append(dict(route=r, n_cells=len(mine),
                           cells=[K.cell_key(*c) for c in mine],
                           replicates=len(mine) * N_REP_PER_CELL,
                           est_hours_colab_2core=len(mine) * N_REP_PER_CELL
                           * sec_per_rep / 3600.0 / colab_cores))
    plan_blob = dict(
        meta=dict(utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
                  python=platform.python_version(), numpy=np.__version__,
                  sec_per_replicate_measured=sec_per_rep, local_cores=ncpu),
        grid=dict(tc=K.TC_PRIMARY, cv=K.CV_PRIMARY, photon=K.PHOTON_PRIMARY,
                  n_cells=n_cells, n_rep_per_cell=N_REP_PER_CELL,
                  total_replicates=total_rep),
        cost=dict(single_core_hours=round(single_core_h, 1),
                  local_multiprocess_hours_all_cores=round(local_h, 2),
                  local_cores_physical=ncpu,
                  memory_safe_cores=safe, gb_per_worker=GB_PER_WORKER,
                  local_multiprocess_hours_memsafe=round(local_safe_h, 2),
                  under_12h_local=bool(local_safe_h < 12.0),
                  recommendation=("LOCAL multiprocess: ~%.1f h on the memory-safe %d "
                                  "workers (or ~%.1f h if all %d cores fit in RAM); "
                                  "both << 12 h -> Colab sharding unnecessary but "
                                  "specified below.  Each worker re-imports "
                                  "numpy/scipy (~%.1f GB); this 16.9 GB box overcommits "
                                  "at 32 workers (page-file failure) -- cap with --nproc."
                                  % (local_safe_h, safe, local_h, ncpu, GB_PER_WORKER))),
        colab_5route=dict(n_routes=N_ROUTES, assumed_cores_per_route=colab_cores,
                          max_route_hours=round(colab_route_h, 2), routes=routes),
        confirmatory_run_reference=dict(records_per_cell=1000, n_arms=4,
                                        approx_local_hours=round(conf_local_h, 2),
                                        note="separate GO; W-eval per record adds "
                                             "the same Whittle-grid cost"),
        discipline=("commit CALIBRATION_TABLE.json BEFORE opening the confirmatory "
                    "bank; no tuning/regridding/estimator-switch/interval-repair after"))
    fn = os.path.join(OUT, "COST_AND_SHARD_PLAN.json")
    json.dump(plan_blob, open(fn, "w"), indent=2)

    print("=== calibration cost + shard plan ===", flush=True)
    print(f"grid: {n_cells} cells x {N_REP_PER_CELL} reps = {total_rep} full-pipeline "
          f"replicates", flush=True)
    print(f"measured {sec_per_rep:.2f} s/replicate -> {single_core_h:.1f} core-hours",
          flush=True)
    print(f"LOCAL multiprocess: {local_safe_h:.2f} h on {safe} memory-safe workers "
          f"(or {local_h:.2f} h on all {ncpu} cores if RAM allows); "
          f"<12h: {local_safe_h < 12.0} -> RECOMMENDED", flush=True)
    print(f"COLAB 5-route split: <= {colab_route_h:.2f} h/route "
          f"(assuming {colab_cores} CPU cores/route)", flush=True)
    for r in routes:
        print(f"   route {r['route']}: {r['n_cells']} cells "
              f"({r['replicates']} reps, ~{r['est_hours_colab_2core']:.1f} h)", flush=True)
    print(f"wrote {fn}", flush=True)
    return plan_blob


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--plan", action="store_true", help="cost + shard manifests only (default)")
    ap.add_argument("--run-cell", default=None, help="tc,cv,ph  e.g. 64,0.40,1.0")
    ap.add_argument("--run-shard", type=int, default=None, help="Colab route id 0..4")
    ap.add_argument("--local-all", action="store_true", help="run ALL 27 cells locally (GO)")
    ap.add_argument("--merge", action="store_true", help="merge partials -> frozen table")
    ap.add_argument("--nrep", type=int, default=N_REP_PER_CELL)
    ap.add_argument("--nproc", type=int, default=None)
    a = ap.parse_args()
    if a.run_cell:
        tc, cv, ph = (float(v) for v in a.run_cell.split(","))
        print(json.dumps(run_cell(tc, cv, ph, a.nrep, a.nproc), indent=2))
    elif a.run_shard is not None:
        run_shard(a.run_shard, a.nrep, a.nproc)
    elif a.local_all:
        os.makedirs(OUT, exist_ok=True)
        out = {}
        for (tc, cv, ph) in K.primary_cells():
            r = run_cell(tc, cv, ph, a.nrep, a.nproc)
            out[r["ckey"]] = r
            print(f"  {r['ckey']}: c_0.95={r['c_0_95']:.3f} ({r['runtime_s']/60:.1f} min)",
                  flush=True)
        json.dump(dict(route="local", cells=out),
                  open(os.path.join(OUT, "partial_local.json"), "w"), indent=2)
        merge()
    elif a.merge:
        merge()
    else:
        plan()
