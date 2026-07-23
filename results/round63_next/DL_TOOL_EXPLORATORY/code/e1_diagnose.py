#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
EXPLORATORY / DEV-ONLY.  NOT preregistered.  No commits.  CPU.

E1 step (a) diagnosis.  Yesterday's head-to-head found the frozen DLGI loses the
fast-strong-drift corner (t_c=16, CV=0.40) on scene PSNR despite near-identical
gain-path CORRELATION.  That signature = a multiplicative gain-LEVEL/scale error,
not a path-shape error.  Here we CONFIRM that decomposition on the 4 test cells:

  * PSNR(x_D)                 raw DLGI reconstruction PSNR
  * PSNR after best scalar rescale (shape-limited PSNR); gap = pure LEVEL loss
  * mean(a_hat)/mean(a_true)  gain-level bias (should be >1 and grow with CV)
  * sum(x_D)/sum(x_true)      scene-brightness bias (inverse of the gain bias)
  * gain-path Pearson r       (should stay ~0.99 -> shape is fine)

If the level loss dominates and grows with CV, the analytic E[a]=1 gauge fix should
close most of the corner deficit WITHOUT any learning.
"""
import os, sys, json, time
from datetime import datetime, timezone
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dl_tool_common as T
import dl_common as DL

OUT = "D:/GI_another/results/round63_next/DL_TOOL_EXPLORATORY"
T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

TC_TEST = [16.0, 64.0]
CV_TEST = [0.15, 0.40]
SEEDS = [0, 1, 2]
SCENES = DL.SCENES


def main():
    cells = []
    for tc in TC_TEST:
        for cv in CV_TEST:
            sig = T.sigma_l_of_cv(cv)
            recs = []
            for sc in SCENES:
                x = DL.scene_x[sc]
                for seed in SEEDS:
                    Y, a_time, slot = DL.simulate_record(sc, seed, tc, sig, schedule="paired")
                    # frozen DLGI
                    jd = DL.joint_dual_ledger(Y, slot=slot)
                    xD = jd["x_hat"]; a_hat = jd["med"]["a_hat"]
                    # oracle
                    x_or = np.clip(DL.arm_A4_gaincorr(Y[0::2], Y[1::2],
                                   a_time[0::2], a_time[1::2], DL.LAM_TV), 0.0, None)
                    psnr_D = DL.psnr(xD, x)
                    psnr_D_scaled, kD = T.best_scale_psnr(xD, x)
                    psnr_or = DL.psnr(x_or, x)
                    recs.append(dict(
                        scene=sc, seed=seed,
                        psnr_D=psnr_D, psnr_D_scaled=psnr_D_scaled,
                        psnr_oracle=psnr_or,
                        level_gain_dB=psnr_D_scaled - psnr_D,   # pure scale recovery
                        shape_gap_dB=psnr_or - psnr_D_scaled,   # residual shape gap
                        gain_ratio=float(np.mean(a_hat) / np.mean(a_time)),
                        bright_ratio=float(np.sum(xD) / max(np.sum(x), 1e-9)),
                        kstar=kD,
                        corr=DL.gain_path_corr(a_hat, a_time),
                    ))
            def m(k): return float(np.mean([r[k] for r in recs]))
            c = dict(tc=tc, cv=cv, n=len(recs),
                     psnr_D=m("psnr_D"), psnr_D_scaled=m("psnr_D_scaled"),
                     psnr_oracle=m("psnr_oracle"),
                     level_gain_dB=m("level_gain_dB"), shape_gap_dB=m("shape_gap_dB"),
                     gain_ratio=m("gain_ratio"), bright_ratio=m("bright_ratio"),
                     kstar=m("kstar"), corr=m("corr"), records=recs)
            cells.append(c)
            log(f"tc={tc:>2.0f} cv={cv:.2f}  PSNR_D={c['psnr_D']:6.2f}  +scale={c['psnr_D_scaled']:6.2f} "
                f"(level +{c['level_gain_dB']:.2f}dB)  oracle={c['psnr_oracle']:6.2f}  "
                f"gain_ratio={c['gain_ratio']:.3f}  bright={c['bright_ratio']:.3f}  r={c['corr']:.4f}")

    meta = dict(exploratory=True, preregistered=False,
                utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
                note="E1 diagnosis: PSNR level-vs-shape decomposition of frozen DLGI",
                tc_test=TC_TEST, cv_test=CV_TEST, seeds=SEEDS, scenes=SCENES,
                runtime_s=time.time() - T0)
    json.dump(dict(meta=meta, cells=[{k: v for k, v in c.items() if k != "records"} for c in cells],
                   records={f"tc{int(c['tc'])}_cv{c['cv']}": c["records"] for c in cells}),
              open(os.path.join(OUT, "json", "e1_diagnose.json"), "w"), indent=2)
    log(f"saved e1_diagnose.json  ({time.time()-T0:.1f}s)")


if __name__ == "__main__":
    main()
