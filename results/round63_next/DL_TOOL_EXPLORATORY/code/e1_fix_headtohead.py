#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
EXPLORATORY / DEV-ONLY.  NOT preregistered.  No commits.

E1 steps (b): re-run yesterday's 4-cell head-to-head with the analytic-fix COPIES
of the DLGI estimator, to answer: does the +3.5 dB fast-strong-drift corner deficit
close WITHOUT any learning, using only the known-physics E[a]=1 gauge?

Arms (identical records, identical arm_A4+TV inversion; only the per-exposure gain
estimate differs):
  arm0        gain = 1  (no-correction floor)
  D0          frozen DLGI (renorm none)                            <- yesterday's Arm D
  D_arith     DLGI + E[a]=1 arithmetic renorm  (the analytic fix)
  D_geom      DLGI + geometric-mean(a)=1 renorm
  D_uni_arith DLGI + uniform gauge + E[a]=1 renorm
  D_noJ_arith DLGI + NO Jensen term + E[a]=1 renorm  (Jensen-contribution check)
  L_mix       SCGI-mix learned corrector (in-family reference; yesterday's best case)
  oracle      true per-exposure gain (ceiling)

Also reports the medium products (t_c, CV rel-err) to confirm the level fix does
NOT change the certified medium estimate (renorm is a pure brightness gauge).
"""
import os, sys, json, time, platform
from datetime import datetime, timezone
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "D:/GI_another/results/round63_next/DL_EXPLORATORY/code")  # scgi_common
import dl_tool_common as T
import dl_common as DL
import scgi_common as S

OUT = "D:/GI_another/results/round63_next/DL_TOOL_EXPLORATORY"
CKPT = "D:/GI_another/results/round63_next/DL_EXPLORATORY/checkpoints"
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

TC_TEST = [16.0, 64.0]
CV_TEST = [0.15, 0.40]
SEEDS = [0, 1, 2]
SCENES = DL.SCENES


def load_model(name):
    ck = torch.load(os.path.join(CKPT, f"{name}.pt"), map_location=DEVICE)
    net = S.UNet(c=ck["channels"]).to(DEVICE)
    net.load_state_dict(ck["state_dict"]); net.eval()
    return net, np.asarray(ck["baseline"], dtype=np.float64)


def relerr(hat, tru): return float(abs(hat - tru) / max(abs(tru), 1e-9))

# DLGI variant configs (all copies; frozen dl_common never mutated)
DVARIANTS = {
    "D0":          dict(gauge="counts",  jensen=True,  renorm="none"),
    "D_arith":     dict(gauge="counts",  jensen=True,  renorm="arith"),
    "D_geom":      dict(gauge="counts",  jensen=True,  renorm="geom"),
    "D_uni_arith": dict(gauge="uniform", jensen=True,  renorm="arith"),
    "D_noJ_arith": dict(gauge="counts",  jensen=False, renorm="arith"),
}


def main():
    net_mix, base_mix = load_model("scgi_mix")
    log(f"scgi_mix loaded  device={DEVICE}")
    cells = []
    for tc in TC_TEST:
        for cv in CV_TEST:
            sig = T.sigma_l_of_cv(cv)
            recs = []
            for sc in SCENES:
                x = DL.scene_x[sc]
                for seed in SEEDS:
                    Y, a_time, slot = DL.simulate_record(sc, seed, tc, sig, schedule="paired")
                    rec = dict(scene=sc, seed=seed, psnr={}, med={}, corr={})
                    # arm0
                    x0 = S.reconstruct_from_gain(Y, np.ones(S.N_EXP))
                    rec["psnr"]["arm0"] = DL.psnr(x0, x)
                    # DLGI variants
                    for name, cfg in DVARIANTS.items():
                        jd = T.joint_dual_ledger_cfg(Y, slot=slot, **cfg)
                        rec["psnr"][name] = DL.psnr(jd["x_hat"], x)
                        rec["med"][name] = dict(tc=jd["med"]["tc_hat"], cv=jd["med"]["cv_hat"])
                        rec["corr"][name] = DL.gain_path_corr(jd["med"]["a_hat"], a_time)
                    # SCGI-mix reference
                    g_mix = S.apply_correction(net_mix, base_mix, Y, DEVICE)
                    rec["psnr"]["L_mix"] = DL.psnr(S.reconstruct_from_gain(Y, g_mix), x)
                    medL = S.medium_from_gain(g_mix)
                    rec["med"]["L_mix"] = dict(tc=medL["tc_hat"], cv=medL["cv_hat"])
                    # oracle
                    x_or = np.clip(DL.arm_A4_gaincorr(Y[0::2], Y[1::2],
                                   a_time[0::2], a_time[1::2], DL.LAM_TV), 0.0, None)
                    rec["psnr"]["oracle"] = DL.psnr(x_or, x)
                    recs.append(rec)
            arms = ["arm0"] + list(DVARIANTS.keys()) + ["L_mix", "oracle"]
            def mp(a): return float(np.mean([r["psnr"][a] for r in recs]))
            def mre(a, key, tru): return float(np.median([relerr(r["med"][a][key], tru) for r in recs]))
            medarms = list(DVARIANTS.keys()) + ["L_mix"]
            c = dict(tc=tc, cv=cv, n=len(recs),
                     psnr={a: mp(a) for a in arms},
                     tc_relerr={a: mre(a, "tc", tc) for a in medarms},
                     cv_relerr={a: mre(a, "cv", cv) for a in medarms},
                     corr={a: float(np.mean([r["corr"][a] for r in recs])) for a in DVARIANTS},
                     records=recs)
            cells.append(c)
            log(f"tc={tc:>2.0f} cv={cv:.2f}  D0={c['psnr']['D0']:6.2f}  D_arith={c['psnr']['D_arith']:6.2f}  "
                f"D_geom={c['psnr']['D_geom']:6.2f}  L_mix={c['psnr']['L_mix']:6.2f}  oracle={c['psnr']['oracle']:6.2f}"
                f"  | D_arith-D0={c['psnr']['D_arith']-c['psnr']['D0']:+.2f}  D_arith-L_mix={c['psnr']['D_arith']-c['psnr']['L_mix']:+.2f}")

    meta = dict(exploratory=True, preregistered=False,
                utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
                python=platform.python_version(), torch=torch.__version__, device=DEVICE,
                note="E1(b) analytic-fix head-to-head; renorm=E[a]=1 physics gauge",
                variants=DVARIANTS, tc_test=TC_TEST, cv_test=CV_TEST, seeds=SEEDS, scenes=SCENES,
                runtime_s=time.time() - T0)
    json.dump(dict(meta=meta, cells=[{k: v for k, v in c.items() if k != "records"} for c in cells],
                   records={f"tc{int(c['tc'])}_cv{c['cv']}": c["records"] for c in cells}),
              open(os.path.join(OUT, "json", "e1_fix_headtohead.json"), "w"), indent=2)
    log("")
    log("=== SCENE PSNR (mean over 6 scenes x 3 seeds), dB ===")
    hdr = f"{'cell':14s}" + "".join(f"{a:>11s}" for a in ["arm0", "D0", "D_arith", "D_geom", "D_uni_arith", "D_noJ_arith", "L_mix", "oracle"])
    log(hdr)
    for c in cells:
        p = c["psnr"]
        log(f"tc={c['tc']:>2.0f} cv={c['cv']:.2f} " +
            "".join(f"{p[a]:>11.2f}" for a in ["arm0", "D0", "D_arith", "D_geom", "D_uni_arith", "D_noJ_arith", "L_mix", "oracle"]))
    log("")
    log("=== medium t_c / CV median rel-err (D_arith should ~= D0: renorm is pure brightness) ===")
    for c in cells:
        log(f"tc={c['tc']:>2.0f} cv={c['cv']:.2f}  tc: D0={c['tc_relerr']['D0']:.3f} D_arith={c['tc_relerr']['D_arith']:.3f} L_mix={c['tc_relerr']['L_mix']:.3f}  "
            f"| cv: D0={c['cv_relerr']['D0']:.3f} D_arith={c['cv_relerr']['D_arith']:.3f} L_mix={c['cv_relerr']['L_mix']:.3f}")
    log(f"saved e1_fix_headtohead.json  ({time.time()-T0:.1f}s)")


if __name__ == "__main__":
    main()
