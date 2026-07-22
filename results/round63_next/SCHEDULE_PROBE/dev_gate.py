#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
R29 section 5.1 DEV gate (formal), Probe B.

Preregistered, outcome-blind decision gate. Reuses the IDENTICAL frozen setup of
probe_b.py (same 972-row SCAT32 bank doubled to N=1944 matched-pair frames, same
dose x4, same ridge lam=0.1*mean(diag(A^T A)), same OU AR(1) drift model, same
schedules, same seed scheme, same drift-BLIND ridge-LS reconstruction). No tuning,
no new arms, no estimator changes.

6-scene DEV cohort, hardest drift cell (sigma=0.40, tc=2.0):
  per scene Q_base = max(Q_ordered, Q_random) 5-seed mean, drift-BLIND recon only.
  dQ_j = Q_paired - Q_base.
Conditions (all must hold for DEV_GATE_PASS):
  (a) median_j dQ_j          >= 1.0 dB
  (b) #{j: dQ_j > 0}         >= 5 of 6
  (c) gain-free control (sigma=0): worst scene loss (Q_base - Q_paired) <= 0.5 dB
"""
import os, json, time
from datetime import datetime, timezone
import numpy as np
from scipy.linalg import cho_factor, cho_solve

T0=time.time()
REPO="D:/GI_another"
OUT=os.path.join(REPO,"results","round63_next","SCHEDULE_PROBE")

# ---------------- frozen inputs / config (identical to probe_b.py) ----------------
rows=np.load(os.path.join(REPO,"results","round63_bridge","library","L0.npz"))["rows"].astype(np.float64)
COHORT=["bridge_contour_1","bridge_microtex_1","bridge_twopop_1","bridge_control_2",
        "bridge_contour_3","bridge_microtex_3"]                       # 4 primary + 2 DEV extensions
scene_x={sc:np.load(os.path.join(REPO,"data","r63_bridge_scenes",sc+".npz"))["x"].astype(np.float64).ravel()
         for sc in COHORT}
patt=np.repeat(rows,2,axis=0); N=patt.shape[0]
tslots=np.linspace(-1.0,1.0,N)
DOSE=4.0; LAM_SCALE=0.1; TC=2.0
SIGMA_HARD=0.40; SIGMA_CTRL=0.0
SEEDS=list(range(5))

def paired_order(N):
    half=N//2; o=np.empty(N,dtype=int)
    for j in range(half): o[half-1-j]=2*j; o[half+j]=2*j+1
    return o
ORDERED=np.arange(N); PAIRED=paired_order(N)
def random_order(seed): return np.random.default_rng(9000+seed).permutation(N)
def ou_path(rng,sigma_l,tc,N):
    phi=np.exp(-1.0/tc); mu=-0.5*sigma_l**2
    eps=rng.standard_normal(N); l=np.empty(N)
    l[0]=mu+sigma_l*eps[0]; sd=sigma_l*np.sqrt(max(1-phi**2,0.0))
    for n in range(1,N): l[n]=mu+phi*(l[n-1]-mu)+sd*eps[n]
    return np.exp(l)

A=DOSE*patt; AtA=A.T@A
LAM=LAM_SCALE*float(np.diag(AtA).mean())
CF=cho_factor(AtA+LAM*np.eye(patt.shape[1]))
def recon_blind(order,Y): return cho_solve(CF,A[order].T@Y)
def psnr(xh,x): return 10*np.log10(1.0/np.mean((xh-x)**2))

def arm_mean_psnr(sc,sigma,order_spec):
    """5-seed mean drift-BLIND PSNR for one arm; order_spec None=random."""
    x=scene_x[sc]; vals=[]
    for seed in SEEDS:
        a_path=ou_path(np.random.default_rng(1000+seed),sigma,TC,N)   # shared-per-seed OU path
        oo=random_order(seed) if order_spec is None else order_spec
        s_slot=A[oo]@x
        Y=np.random.default_rng(5000+seed).poisson(a_path*s_slot).astype(np.float64)
        vals.append(psnr(recon_blind(oo,Y),x))
    return float(np.mean(vals))

def eval_cell(sigma):
    per={}
    for sc in COHORT:
        Qo=arm_mean_psnr(sc,sigma,ORDERED)
        Qr=arm_mean_psnr(sc,sigma,None)
        Qp=arm_mean_psnr(sc,sigma,PAIRED)
        Qb=max(Qo,Qr)
        per[sc]=dict(Q_ordered=Qo,Q_random=Qr,Q_paired=Qp,Q_base=Qb,dQ=Qp-Qb)
    return per

print(f"[{time.time()-T0:.1f}s] hardest cell sigma={SIGMA_HARD}, tc={TC} ...",flush=True)
hard=eval_cell(SIGMA_HARD)
print(f"[{time.time()-T0:.1f}s] gain-free control sigma={SIGMA_CTRL} ...",flush=True)
ctrl=eval_cell(SIGMA_CTRL)

per_scene_dQ={sc:hard[sc]["dQ"] for sc in COHORT}
dQ=np.array([per_scene_dQ[sc] for sc in COHORT])
median=float(np.median(dQ))
n_positive=int((dQ>0).sum())
# control loss = how far paired falls below the best baseline when there is no drift
ctrl_loss={sc:ctrl[sc]["Q_base"]-ctrl[sc]["Q_paired"] for sc in COHORT}
control_worst_loss=float(max(ctrl_loss.values()))

condA=dict(name="median dQ >= 1.0 dB", value=median, threshold=1.0, passed=bool(median>=1.0))
condB=dict(name="#scenes dQ>0 >= 5 of 6", value=n_positive, threshold=5, passed=bool(n_positive>=5))
condC=dict(name="control worst loss <= 0.5 dB (sigma=0)", value=control_worst_loss,
           threshold=0.5, passed=bool(control_worst_loss<=0.5))
verdict="DEV_GATE_PASS" if (condA["passed"] and condB["passed"] and condC["passed"]) else "DEV_GATE_FAIL"

# descriptive by-products (honest, not rescues): paired-ordered from this hard cell + primary pooled
paired_minus_ordered_hard=float(np.mean([hard[sc]["Q_paired"]-hard[sc]["Q_ordered"] for sc in COHORT]))
primary_pooled=None
pj=os.path.join(OUT,"probe_b_results.json")
if os.path.exists(pj):
    P=json.load(open(pj))
    ds=[]
    for sc in P["primary"]:
        for sg in P["primary"][sc]:
            o=np.mean(P["primary"][sc][sg]["ordered"]["blind"]["vals"])
            p=np.mean(P["primary"][sc][sg]["paired"]["blind"]["vals"])
            ds.append(p-o)
    primary_pooled=float(np.mean(ds))

RES=dict(
    generated_utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
    gate="R29 section 5.1 DEV gate (formal)",
    estimator="drift-BLIND ridge LS (identical to probe_b.py)",
    cell=dict(sigma=SIGMA_HARD,tc=TC,seeds=SEEDS,dose=DOSE,lam_scale=LAM_SCALE,N_frames=N),
    cohort=COHORT,
    per_scene=hard,
    per_scene_dQ=per_scene_dQ,
    median=median,
    n_positive=n_positive,
    control_worst_loss=control_worst_loss,
    control_sigma0=dict(per_scene={sc:dict(Q_base=ctrl[sc]["Q_base"],Q_paired=ctrl[sc]["Q_paired"],
                                           loss=ctrl_loss[sc]) for sc in COHORT},
                        control_worst_loss=control_worst_loss),
    conditions=dict(a=condA,b=condB,c=condC),
    verdict=verdict,
    descriptive_facts=dict(
        paired_minus_ordered_dB_primary_pooled=primary_pooled,
        paired_minus_ordered_dB_hard_cell_6scene=paired_minus_ordered_hard,
        note_anti_naive_ordering=("paired beats naive bank-order by ~1.3 dB pooled (up to ~2 dB at "
            "severe drift); this anti-naive-ordering effect is prior-art-adjacent (interleaving), "
            "and is the only image-level survivor of the Fisher-side ordering."),
        dominant_lever=("correlation time vs acquisition length: at the frozen tc=2 (<< N=%d) most OU "
            "power is uncancellable high-frequency residual r (paper model l=H beta + r), so the "
            "low-order moment cancellation the paired schedule performs has little image-level "
            "purchase." % N),
    ),
)
jp=os.path.join(OUT,"DEV_GATE_VERDICT.json")
with open(jp,"w",encoding="utf-8") as f: json.dump(RES,f,indent=2)
print(f"[{time.time()-T0:.1f}s] wrote {jp}",flush=True)
print("per-scene dQ:",{sc:round(v,3) for sc,v in per_scene_dQ.items()})
print(f"median dQ={median:.3f}  n_positive={n_positive}/6  control_worst_loss={control_worst_loss:.3f}")
print(f"conditions: a={condA['passed']} b={condB['passed']} c={condC['passed']}")
print(f"VERDICT: {verdict}")
