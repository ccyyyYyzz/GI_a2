#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
VM-AUTONOMOUS SEALED-PROBE SHARD RUNNER (R39 Phase-2 machinery).

Runs ONE shard of the beyond-band confirmatory on a Colab VM (or locally), fully autonomous so
token/network loss never kills computation, per COLAB_USAGE_GUIDE.md discipline:

  * task list lives in the shard MANIFEST on the VM (sealed_probe_planner.py output);
  * META-as-truth resumability: <shard>_meta.json records completed cell_ids; a re-invocation
    after a VM recycle SKIPS done cells and continues;
  * per-cell CHECKPOINT: each finished cell is appended to the shard CSV AND its cell_id added to
    meta.json atomically before the next cell starts;
  * HEARTBEAT JSON (epoch, shard, cells_done/total, rows) rewritten after every cell, so the
    watchdog (live_watch.sh) can detect stalls;
  * wall-budget guard: returns 2 (resumable) if the accumulated budget is exceeded.

CELL kinds (from the manifest): 'arms' (six arms + multi-start std), 'mismatch' (F5 degradation),
'aperture' (aperture-law in/out separation).  Consumes sealed banks under scene_banks/.

This is BUILD-ONLY infrastructure: the confirmatory is UNRUN until the R39 ruling freezes the
bars (bars64.THRESHOLDS).  A tiny n=16 selftest (--selftest) exercises the cell dispatch,
checkpointing, and resume path end-to-end.  Runnable under remote_lane.py or session_driver.sh.
"""
import argparse
import csv
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import arms64 as A                                       # noqa: E402
import fog_tracker as ft                                 # noqa: E402  (via arms64 sys.path)

BANKS_DIR = os.path.join(HERE, "scene_banks")


def load_scene(partition, scene_id):
    x = np.load(os.path.join(BANKS_DIR, partition, scene_id + ".npz"))["x"]
    return np.asarray(x, np.float64).ravel()


# --------------------------------------------------------------- cell runners
def run_arms_cell(cfg, cell, fr):
    """Six arms + covariance multi-start agreement std on one shared medium realization."""
    x = load_scene(cell["partition"], cell["scene_id"])
    med_rng, poi_rng = A.record_rngs(cell["partition"], cell["scene_id"], cell["replicate"])
    W = cfg.medium_field(fr["T"], med_rng)
    Pfix = cfg.bl_patterns(cfg.M, 10)                    # sealed fixed bank [R39]
    pat_base = 7_000_000 + A._u32(cell["scene_id"]) % 1_000_000
    row = dict(cell_id=cell["cell_id"], kind="arms", scene=cell["scene_id"],
               replicate=cell["replicate"])
    # oracle
    row["oracle"] = cfg.beyond_band_nmse(
        A.arm_oracle(cfg, x, W, np.random.default_rng(poi_rng.integers(1 << 31)), Pfix, fr["photons"]), x)
    # fixed moment (+ multi-start agreement std for B4)
    mu = np.einsum("mn,tn->tm", Pfix, W * x)
    B, R_diag = A._poisson_buckets(mu, np.random.default_rng(poi_rng.integers(1 << 31)), fr["photons"])
    rc = ft.cov_estimate(Pfix, cfg.Ub, B, cfg.coeff_sd, cfg.rho,
                         n_starts=fr["n_starts"], steps=fr["steps"], dev=A.DEV)
    per_start = [cfg.beyond_band_nmse(xs, x) for xs in rc["per_start_x"]]
    row["fixed_moment"] = float(per_start[int(np.argmin(rc["objs"]))])
    row["multistart_std"] = float(np.std(per_start))
    # fixed MLE (seeded by the best moment start)
    x_init = rc["per_start_x"][int(np.argmin(rc["objs"]))]
    xL = ft.kalman_em(Pfix, cfg.Ub, B, R_diag, cfg.rho, cfg.coeff_sd, 1e-4,
                      n_em=fr["n_em"], dev=A.DEV, x_init=x_init)
    row["fixed_mle"] = cfg.beyond_band_nmse(xL, x)
    # wall arms + fresh-cov-gls design comparator
    row["fresh_mean"] = cfg.beyond_band_nmse(
        A.arm_fresh_mean(cfg, x, W, np.random.default_rng(poi_rng.integers(1 << 31)),
                         fr["photons"], pat_base), x)
    row["classic_avg"] = cfg.beyond_band_nmse(
        A.arm_classic_avg(cfg, x, W, np.random.default_rng(poi_rng.integers(1 << 31)), Pfix,
                          fr["photons"]), x)
    row["fresh_cov_gls"] = cfg.beyond_band_nmse(
        A.arm_fresh_cov_gls(cfg, x, W, np.random.default_rng(poi_rng.integers(1 << 31)),
                            fr["photons"], pat_base), x)
    return [row]


def run_aperture_cell(cfg, cell, fr):
    """Aperture-law in/out coverage separation on a calibration scene."""
    x = load_scene(cell["partition"], cell["scene_id"])
    med_rng, poi_rng = A.record_rngs(cell["partition"], cell["scene_id"], 0)
    W = cfg.medium_field(fr["T"], med_rng)
    Pfix = cfg.bl_patterns(cfg.M, 10)
    mu = np.einsum("mn,tn->tm", Pfix, W * x)
    B, _ = A._poisson_buckets(mu, poi_rng, fr["photons"])
    rc = ft.cov_estimate(Pfix, cfg.Ub, B, cfg.coeff_sd, cfg.rho,
                         n_starts=fr["n_starts"], steps=fr["steps"], dev=A.DEV)
    xh = rc["per_start_x"][int(np.argmin(rc["objs"]))]
    s = np.dot(xh, x) / max(np.dot(xh, xh), 1e-12)
    E = cfg.D2(s * xh - x) ** 2; X0 = cfg.D2(x) ** 2
    inm = cfg.cov & (X0 > 1e-9); outm = (~cfg.cov) & (X0 > 1e-9)
    ine = float(E[inm].sum() / max(X0[inm].sum(), 1e-12))
    oute = float(E[outm].sum() / max(X0[outm].sum(), 1e-12))
    return [dict(cell_id=cell["cell_id"], kind="aperture", scene=cell["scene_id"],
                 in_coverage_err=ine, out_coverage_err=oute, separation=oute / max(ine, 1e-9))]


def run_mismatch_cell(cfg, cell, fr):
    """One F5 mismatch axis/condition: beyond-band NMSE degradation vs the matched estimator."""
    x = load_scene(cell["partition"], cell["scene_id"])
    med_rng, poi_rng = A.record_rngs(cell["partition"], cell["scene_id"], 0)
    Pfix = cfg.bl_patterns(cfg.M, 10)
    axis, cond = cell["axis"], cell["condition"]

    def moment_err(W, declared_Ub, declared_sd, declared_rho):
        mu = np.einsum("mn,tn->tm", Pfix, W * x)
        B, _ = A._poisson_buckets(mu, np.random.default_rng(poi_rng.integers(1 << 31)), fr["photons"])
        rc = ft.cov_estimate(Pfix, declared_Ub, B, declared_sd, declared_rho,
                             n_starts=3, steps=fr["steps"], dev=A.DEV)
        return cfg.beyond_band_nmse(rc["per_start_x"][int(np.argmin(rc["objs"]))], x)

    # matched baseline (true = frozen law, declared = frozen law)
    W0 = cfg.medium_field(fr["T"], med_rng)
    e_match = moment_err(W0, cfg.Ub, cfg.coeff_sd, cfg.rho)

    # perturbed condition
    mr2 = np.random.default_rng(med_rng.integers(1 << 31))
    if axis == "fine_rotation":
        eps = {"rot10": 0.10, "rot20": 0.20}[cond]
        e_pert = moment_err(cfg.medium_field(fr["T"], mr2), cfg.rotate_declared(eps),
                            cfg.coeff_sd, cfg.rho)
    elif axis == "band_width":
        hi = {"narrow": cfg.med_hi - 2, "wide": cfg.med_hi + 2}[cond]
        cfg2 = A.ProbeConfig(n=cfg.n, PB=cfg.PB, med_lo=cfg.med_lo, med_hi=hi,
                             sig_f=cfg.sig_f, tau=cfg.tau, M=cfg.M)
        e_pert = moment_err(cfg.medium_field(fr["T"], mr2), cfg2.Ub, cfg2.coeff_sd, cfg.rho)
    elif axis == "radial_profile":
        modes = [(i, j) for i in range(cfg.n) for j in range(cfg.n)
                 if cfg.med_lo <= i + j <= cfg.med_hi]
        prof = np.array([1.0 / (1.0 + np.hypot(i, j)) for (i, j) in modes])
        gen_sd = cfg.coeff_sd * prof / np.sqrt((prof ** 2).mean())
        Wp = cfg.medium_field_custom(fr["T"], mr2, gen_sd=gen_sd)
        e_pert = moment_err(Wp, cfg.Ub, cfg.coeff_sd, cfg.rho)     # declared flat
    elif axis == "geometry_alpha":
        from scipy.ndimage import gaussian_filter
        alpha = {"a0.1": 0.1, "a0.25": 0.25, "a0.5": 0.5}[cond]
        Wp = cfg.medium_field(fr["T"], mr2)
        PW = Pfix[None, :, :] * Wp[:, None, :]
        PWb = gaussian_filter(PW.reshape(fr["T"], cfg.M, cfg.n, cfg.n),
                              sigma=(0, 0, 0.8, 0.8), mode="reflect").reshape(fr["T"], cfg.M, -1)
        Aeff = (1 - alpha) * PW + alpha * PWb
        mu = np.einsum("tmn,n->tm", Aeff, x)
        B, _ = A._poisson_buckets(mu, np.random.default_rng(poi_rng.integers(1 << 31)), fr["photons"])
        rc = ft.cov_estimate(Pfix, cfg.Ub, B, cfg.coeff_sd, cfg.rho, n_starts=3,
                             steps=fr["steps"], dev=A.DEV)
        e_pert = cfg.beyond_band_nmse(rc["per_start_x"][int(np.argmin(rc["objs"]))], x)
    elif axis == "wrong_tau":
        tau_true = {"tau_half": cfg.tau / 2, "tau_double": cfg.tau * 2}[cond]
        Wp = cfg.medium_field_custom(fr["T"], mr2, rho=float(np.exp(-1.0 / tau_true)))
        e_pert = moment_err(Wp, cfg.Ub, cfg.coeff_sd, cfg.rho)     # declared frozen tau
    elif axis == "wrong_sigma_f":
        sf_true = {"sf_0.2": 0.2, "sf_0.4": 0.4}[cond]
        sd_true = sf_true * np.sqrt(cfg.N / cfg.db)
        Wp = cfg.medium_field_custom(fr["T"], mr2, coeff_sd=sd_true)
        e_pert = moment_err(Wp, cfg.Ub, cfg.coeff_sd, cfg.rho)     # declared frozen sigma_f
    else:
        raise ValueError(axis)

    degr = (e_pert - e_match) / max(e_match, 1e-9)
    return [dict(cell_id=cell["cell_id"], kind="mismatch", scene=cell["scene_id"],
                 axis=axis, condition=cond, err_matched=e_match, err_perturbed=e_pert,
                 degradation=degr)]


CELL_RUNNERS = {"arms": run_arms_cell, "aperture": run_aperture_cell, "mismatch": run_mismatch_cell}


# ------------------------------------------------- meta-as-truth + heartbeat
def _load_meta(meta_path):
    if os.path.exists(meta_path):
        return json.load(open(meta_path))
    return dict(done_cell_ids=[], rows=0)


def _atomic_write(path, obj):
    tmp = path + ".tmp"
    json.dump(obj, open(tmp, "w"))
    os.replace(tmp, path)


def _heartbeat(hb_path, shard, done, total, rows, rc=None):
    _atomic_write(hb_path, dict(epoch=int(time.time()), shard=shard, cells_done=done,
                                cells_total=total, rows=rows, returncode=rc))


def run_shard(manifest_path, wall_budget_s=21600.0, hb_dir=None):
    man = json.load(open(manifest_path))
    shard = man["shard_id"]; fr = man["frozen"]; cells = man["cells"]
    n = 64 if fr.get("scale", 64) == 64 else 16
    cfg = A.ProbeConfig(n=n, PB=8 if n == 64 else 3, med_lo=1, med_hi=16 if n == 64 else 6,
                        sig_f=0.30, tau=8.0, M=128 if n == 64 else 24)
    out_csv = os.path.join(HERE, "shards", f"{shard}.csv")
    meta_path = os.path.join(HERE, "shards", f"{shard}_meta.json")
    hb_dir = hb_dir or os.path.join(HERE, "shards")
    os.makedirs(os.path.join(HERE, "shards"), exist_ok=True)
    hb_path = os.path.join(hb_dir, f"{shard}.hb.json")
    meta = _load_meta(meta_path)
    done = set(meta["done_cell_ids"])
    fields = ["cell_id", "kind", "scene", "replicate", "oracle", "fixed_mle", "fixed_moment",
              "multistart_std", "fresh_mean", "classic_avg", "fresh_cov_gls", "in_coverage_err",
              "out_coverage_err", "separation", "axis", "condition", "err_matched",
              "err_perturbed", "degradation"]
    new_file = not os.path.exists(out_csv)
    t0 = time.time()
    with open(out_csv, "a", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        if new_file:
            wr.writeheader()
        for c in cells:
            if c["cell_id"] in done:
                continue
            if time.time() - t0 > wall_budget_s:
                _heartbeat(hb_path, shard, len(done), len(cells), meta["rows"], rc=2)
                print(f"[{shard}] wall-budget abort (resumable); {len(done)}/{len(cells)} done")
                return 2
            rows = CELL_RUNNERS[c["kind"]](cfg, c, fr)
            for r in rows:
                wr.writerow(r); meta["rows"] += 1
            fh.flush()
            done.add(c["cell_id"]); meta["done_cell_ids"] = sorted(done)
            _atomic_write(meta_path, meta)
            _heartbeat(hb_path, shard, len(done), len(cells), meta["rows"])
            print(f"[{shard}] cell {c['cell_id']} done ({len(done)}/{len(cells)})", flush=True)
    _heartbeat(hb_path, shard, len(done), len(cells), meta["rows"], rc=0)
    print(f"[{shard}] COMPLETE {len(done)}/{len(cells)} cells, {meta['rows']} rows -> {out_csv}")
    return 0


def _selftest():
    """Tiny n=16 end-to-end: plan a 3-cell shard, run it, verify resume skips done cells."""
    import tempfile
    fr = dict(T=512, photons=1e5, replicates=1, steps=1500, n_starts=3, n_em=10, scale=16)
    sid = "selftest_witness_0"                          # throwaway id -- never a sealed scene
    cells = [dict(kind="arms", partition="confirmatory", scene_id=sid,
                  replicate=0, cell_id="arms::w0::r0")]
    # need a 16x16 scene on disk: synthesize one into the confirmatory dir under a throwaway id
    cfg = A.ProbeConfig(n=16, PB=3, med_lo=1, med_hi=6, M=24)
    rs = np.random.default_rng(4)
    C = np.zeros((16, 16)); C[:4, :4] = rs.standard_normal((4, 4))
    for (i, j) in [(5, 2), (2, 5), (4, 4), (6, 1), (1, 6), (5, 5), (7, 0), (0, 7)]:
        C[i, j] = 2.0 * rs.choice([-1, 1])
    x = cfg.I2(C); x = x - x.min(); x = (x / x.max()).reshape(16, 16)
    d = os.path.join(BANKS_DIR, "confirmatory"); os.makedirs(d, exist_ok=True)
    np.savez(os.path.join(d, sid + ".npz"), x=x)
    man = dict(shard_id="selftest_shard", frozen=fr, cells=cells)
    mp = os.path.join(tempfile.gettempdir(), "selftest_manifest.json")
    json.dump(man, open(mp, "w"))
    rc1 = run_shard(mp)
    rc2 = run_shard(mp)                                  # resume: should skip the done cell
    meta = json.load(open(os.path.join(HERE, "shards", "selftest_shard_meta.json")))
    ok = rc1 == 0 and rc2 == 0 and len(meta["done_cell_ids"]) == 1
    # cleanup the throwaway 16x16 test scene + shard artifacts (sealed 64x64 banks untouched)
    for p in [os.path.join(d, sid + ".npz"),
              os.path.join(HERE, "shards", "selftest_shard.csv"),
              os.path.join(HERE, "shards", "selftest_shard_meta.json"),
              os.path.join(HERE, "shards", "selftest_shard.hb.json")]:
        if os.path.exists(p):
            os.remove(p)
    print(f"\nRUNNER SELFTEST: {'PASS' if ok else 'CHECK'} "
          f"(ran + resumed; done={meta['done_cell_ids']})")
    return ok


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="R39 beyond-band sealed-probe VM-autonomous shard runner")
    ap.add_argument("--manifest", help="shard manifest JSON (from sealed_probe_planner.py)")
    ap.add_argument("--wall-budget-s", type=float, default=21600.0)
    ap.add_argument("--heartbeat-dir", default=None)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if _selftest() else 1)
    if not a.manifest:
        ap.error("--manifest required (or --selftest)")
    sys.exit(run_shard(a.manifest, a.wall_budget_s, a.heartbeat_dir))
