#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
PHASE 2 DRY RUN -- CALIBRATION-SIDE D1 + attribution-simplex Gram (R41 sec 4.6 D1 + sec 9 Rank 2).

Runs ONLY on the calibration bank (confirmatory stays sealed). Two deliverables:
  D1 (calibration side): analytic vs empirical d' (FIXED-COV) across cells x anomaly sizes on the
                         calibration scenes; require the ratio in [0.80,1.20] (median ARE <= 10%).
  Simplex Gram         : the efficient-score canonical-correlation structure on a calibration cell,
                         + the five-class LDA balanced accuracy (bar D3 machinery).

Writes SIMPLEX_CALIBRATION.json (the .md is written from it). No confirmatory endpoint is touched.
"""
import json
import os
import time

import numpy as np

import sealed_common as sc
import sealed_banks as sb
import arms
import simplex as sx

HERE = os.path.dirname(os.path.abspath(__file__))
t0 = time.time()


def load_calibration_scenes():
    """Regenerate the calibration bank scenes deterministically from the committed manifest seeds."""
    man = sb.build_bank("calibration", write=False)
    scenes = {}
    for e in man["scenes"]:
        if e["family"] == "natural":
            x = sb.nat_scene(e["source"])
        else:
            x = sb.synth_scene(e["family"], e["seed"])
        # verify the committed hash
        assert sb.sha(x) == e["x_sha256"], f"hash mismatch {e['scene_id']}"
        scenes[e["scene_id"]] = np.asarray(x, float)
    return scenes


def d1_calibration(scenes):
    """FIXED-COV analytic vs empirical d' on calibration scenes across best/mid/floor x eps."""
    cell_geos = [("best", sc.BEST_CELL), ("mid", sc.MID_CELL), ("floor", sc.FLOOR_CELL)]
    # two calibration scenes: a witness (beyond-band structure) and a natural
    use = ["calibration_witness_0", "calibration_nat_camera"]
    ratios = []
    rows = []
    for scene_id in use:
        x = scenes[scene_id]
        for tag, geo in cell_geos:
            cell = sc.setup_cell(**geo, x_np=x)
            xn = np.sqrt(cell["xnorm2"])
            for eps in [0.01, 0.02, 0.05]:
                c = sc.energy_spread_delta(cell["db"], eps, xn, seed=int(1000 * eps) + 7)
                delta_px = cell["Phi_b"].cpu().numpy() @ c
                Tdet = sc.DP_STRONG2 / (cell["lam_mean"] * eps * eps * cell["xnorm2"])
                T_eff = int(min(max(round(Tdet), 64), sc.T_CAP))
                dpa = arms.dprime_analytic(cell, c, T_eff)
                mc = arms.fixed_cov_mc(cell, c, delta_px, T_eff, 300, rng0=int(1e5) + int(1e3 * eps))
                ratio = mc["dprime"] / max(dpa, 1e-9)
                ratios.append(ratio)
                rows.append(dict(scene=scene_id, cell=tag, eps=eps, T_eff=T_eff,
                                 dprime_analytic=round(dpa, 2), dprime_empirical=mc["dprime"],
                                 auc=mc["auc"], ratio=round(ratio, 3)))
    return rows, ratios


def main():
    import bars
    print("=== D1 (calibration side) + attribution simplex ===", flush=True)
    scenes = load_calibration_scenes()
    print(f"loaded + hash-verified {len(scenes)} calibration scenes", flush=True)

    d1_rows, ratios = d1_calibration(scenes)
    d1 = bars.eval_D1(ratios)
    print(f"\nD1 (calibration): {len(ratios)} points  median ratio={np.median(ratios):.3f}  "
          f"[{min(ratios):.3f},{max(ratios):.3f}]  median|ARE|={np.median(np.abs(np.array(ratios)-1)):.3f}",
          flush=True)
    for r in d1_rows:
        print(f"  {r['scene']:24s} {r['cell']:5s} eps={r['eps']:.2f} T={r['T_eff']:4d}  "
              f"d'_an={r['dprime_analytic']:.2f} d'_emp={r['dprime_empirical']:.2f} ratio={r['ratio']:.2f}",
              flush=True)
    print(f"  D1 checks: {d1['checks']}  -> {'PASS' if d1['passed'] else 'FAIL'}", flush=True)

    # simplex on a calibration cell (best geometry, calibration witness scene)
    print("\n=== attribution simplex (calibration witness, best-cell geometry) ===", flush=True)
    cell = sc.setup_cell(**sc.BEST_CELL, x_np=scenes["calibration_witness_0"])
    gram = sx.simplex_gram(cell)
    print(f"  joint dim={gram['joint_dim']}  max off-diag RAW cc={gram['max_offdiag_raw_cc']:.4f}  "
          f"max off-diag EFFICIENT cc={gram['max_offdiag_efficient_cc']:.2e}", flush=True)
    for pair, v in gram["raw_canonical_correlations"].items():
        print(f"    {pair:26s} max_cc={v['max_cc']:.4f}", flush=True)
    beyond_amp = gram["raw_canonical_correlations"]["beyond__amplitude"]["max_cc"]
    beyond_lag = gram["raw_canonical_correlations"]["beyond__lag"]["max_cc"]
    print(f"  KEY specificity orthogonality: beyond<->amplitude={beyond_amp:.4f}, "
          f"beyond<->lag={beyond_lag:.4f}  (target < {gram['target_max_offdiag']})", flush=True)

    # five-class classifier at a within-cap T_eff
    print("\n=== five-class attribution (LDA) ===", flush=True)
    ba_rows = {}
    for T in [1200, 2048]:
        pops = sx.five_class_populations(cell, T, 300, 0.02, rng0=2000 + T)
        ba = sx.balanced_accuracy(pops)
        ba_rows[T] = ba
        print(f"  T={T}: balanced_accuracy={ba['balanced_accuracy']:.3f}  recall={ba['per_class_recall']}",
              flush=True)

    out = dict(
        d1_calibration=dict(rows=d1_rows, eval=d1, n=len(ratios),
                            median_ratio=float(np.median(ratios))),
        simplex=gram,
        five_class=ba_rows,
        beyond_amplitude_cc=beyond_amp, beyond_lag_cc=beyond_lag,
        wall_s=round(time.time() - t0, 1))
    json.dump(out, open(os.path.join(HERE, "SIMPLEX_CALIBRATION.json"), "w"), indent=2)
    print(f"\n[D1-cal + simplex done {time.time()-t0:.0f}s]", flush=True)
    return out


if __name__ == "__main__":
    main()
