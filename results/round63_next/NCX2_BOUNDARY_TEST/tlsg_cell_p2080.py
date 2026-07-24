#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""TEST 1(b) measurement: ONE new TLSG cell at p=2080 (M=64, M(M+1)/2=2080), single seed
20260724, CPU-forced (the GPU exam PID 19876 must not be contended).

The engine functions (sym_basis, make_cell, scores_ZT, power_curve, find_T50) are loaded
BYTE-IDENTICAL from the committed runner
  results/round63_next/TLSG/tlsg_A_three_ledger.py
by exec'ing its module prefix (everything before the 8-cell driver loop `results = {`).
Only device (cpu) and CHUNK (1, for the 3.3 GB RAM ceiling) differ from a GPU run; N_MC=1500
matches the committed bootstrap engine setting. The cell's per-cell torch generator is seeded
exactly as the runner: manual_seed(SEED + M*1000 + p). The numpy `rng` for the random tangent
subspace is re-seeded default_rng(SEED) so the (64,2080) subspace is the first draw (the blind
50%-power contour is realization-robust, so stream position is immaterial)."""
import os, sys, time, json
# Device policy: the GPU exam (PID 19876) has FINISHED (confirmed gone from nvidia-smi
# compute-apps), so the now-idle GPU may be used per the task's conditional. Set
# CUDA_VISIBLE_DEVICES=-1 in the environment to force CPU instead. CHUNK kept small so the
# BY tensor (CHUNK*Tmax*p*M float64) fits VRAM/RAM: CHUNK=2 -> ~3.4 GB.
os.environ["N_MC"] = os.environ.get("N_MC", "1500")
os.environ["CHUNK"] = os.environ.get("CHUNK", "2")
os.environ["SEED"] = os.environ.get("SEED", "20260724")
os.environ["LAM"] = os.environ.get("LAM", "0.15")

HERE = os.path.dirname(os.path.abspath(__file__))
RUNNER = os.path.abspath(os.path.join(HERE, "..", "TLSG", "tlsg_A_three_ledger.py"))
src = open(RUNNER, "r", encoding="utf-8").read()
prefix = src.split("\nresults = {")[0]           # engine functions + constants ONLY
ns = {}
exec(compile(prefix, RUNNER, "exec"), ns)
import numpy as np, torch

LAM = ns["LAM"]; SEED = ns["SEED"]; FPR = ns["FPR"]
M, p = 64, 2080
t0 = time.time()
def log(m):
    line = "[%.0fs] %s" % (time.time()-t0, m)
    print(line, flush=True)
    with open(os.path.join(HERE, "p2080.log"), "a", encoding="utf-8") as f:
        f.write(line + "\n")

log("start cell M=%d p=%d  DEV=%s N_MC=%s CHUNK=%s LAM=%s SEED=%s" %
    (M, p, ns["DEV"], ns["N_MC"], ns["CHUNK"], LAM, SEED))
ns["rng"] = np.random.default_rng(SEED)           # fresh -> subspace is first draw
B, trB, u, theta, chol, I = ns["make_cell"](M, p)
log("make_cell done B=%s" % (tuple(B.shape),))
Tmax = int(min(1600, max(60, 8 * 2.33 * np.sqrt(p) / LAM)))
Tgrid = sorted(set(int(round(x)) for x in np.unique(np.geomspace(4, Tmax, 16).astype(int)) if x >= 2))
log("Tmax=%d Tgrid=%s" % (Tmax, Tgrid))
gen = torch.Generator(device=ns["DEV"]).manual_seed(SEED + M*1000 + p)
ZT_null = ns["scores_ZT"](chol * 0 + I, B, trB, Tgrid, gen)
log("null scores done")
ZT_alt = ns["scores_ZT"](chol, B, trB, Tgrid, gen)
log("alt scores done")
matched_pw = []; blind_pw = []; rel_risk = []
for it, T in enumerate(Tgrid):
    zn = ZT_null[:, it]; za = ZT_alt[:, it]
    mp, _ = ns["power_curve"](zn @ u, za @ u)
    bp, _ = ns["power_curve"]((zn ** 2).sum(1), (za ** 2).sum(1))
    risk = float(((za / np.sqrt(T) - theta[None]) ** 2).sum(1).mean())
    matched_pw.append(mp); blind_pw.append(bp); rel_risk.append(risk / LAM)
T50m = ns["find_T50"](Tgrid, matched_pw); T50b = ns["find_T50"](Tgrid, blind_pw)
Tlam_50_matched = (T50m * LAM) if T50m else None
Tlam_50_blind = (T50b * LAM) if T50b else None
log("Tlam_50_blind=%s  Tlam_50_matched=%s" % (Tlam_50_blind, Tlam_50_matched))

# local-slope measurement vs committed p=528 endpoints (from TLSG_partA.json, seed 20260724)
committed_p528 = {"M64": 59.46055755490618, "M32": 57.75566906294253}
d50_exact_528 = 55.28833238  # scipy exact ncx2 boundary FPR=0.05 (this run, part a)
d50_exact_2080 = 107.91257382
pred_local_slope = (np.log(d50_exact_2080) - np.log(d50_exact_528)) / (np.log(2080) - np.log(528))
out = {
    "test": "TLSG_cell_p2080_partB_measurement",
    "engine": "tlsg_A_three_ledger.py prefix exec'd verbatim; DEV=%s; CHUNK=%s; N_MC=%s" %
              (ns["DEV"], ns["CHUNK"], ns["N_MC"]),
    "cell": dict(M=M, p_eff=p, LAM=LAM, FPR=FPR, seed=SEED, n_mc=int(ns["N_MC"]),
                 Tgrid=Tgrid, matched_power=matched_pw, blind_power=blind_pw),
    "Tlam_50_blind": Tlam_50_blind, "Tlam_50_matched": Tlam_50_matched,
    "committed_p528_Tlam_50_blind": committed_p528,
    "exact_boundary": dict(delta50_528=d50_exact_528, delta50_2080=d50_exact_2080,
                           predicted_local_slope=pred_local_slope),
}
if Tlam_50_blind:
    ls_M64 = (np.log(Tlam_50_blind) - np.log(committed_p528["M64"])) / (np.log(2080) - np.log(528))
    ls_M32 = (np.log(Tlam_50_blind) - np.log(committed_p528["M32"])) / (np.log(2080) - np.log(528))
    out["measured_local_slope_vs_M64_528"] = float(ls_M64)
    out["measured_local_slope_vs_M32_528"] = float(ls_M32)
    out["abs_err_vs_pred_M64"] = float(abs(ls_M64 - pred_local_slope))
    log("measured local slope vs (64,528)=%.4f  vs (32,528)=%.4f  pred=%.4f" %
        (ls_M64, ls_M32, pred_local_slope))
out["elapsed_s"] = round(time.time() - t0, 1)
with open(os.path.join(HERE, "p2080_result.json"), "w") as f:
    json.dump(out, f, indent=2)
log("DONE elapsed %.1fs -> p2080_result.json" % out["elapsed_s"])
