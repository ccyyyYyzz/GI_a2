#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
Probe B  --  ROUND63-NEXT software-method hunt.

Question: the Part II scattering theory (Thm O4-A) proves a pattern schedule is
drift-optimal iff the moment condition  M^T R^{-1} D_s H = 0  holds, and the
frozen Fisher-side verification (results/scatter_verify/, commit fb3feda) measured
the Schur-complement loss tr(I_xx - I_x|gain):  ordered 1770 >> random 63.6 >
PAIRED 0.0 -- a three-orders-of-magnitude Fisher-side advantage that is *pure
software* (pattern ORDER only).  Does that advantage survive to IMAGE level
(PSNR of actual reconstructions)?

This probe is outcome-blind in design and evaluated ONLY at image level.

Frozen drift model (lifted verbatim from code/scatter/verify_o1_o4.py):
  * multiplicative gain  a_n = exp(l_n),  l_n a stationary OU log-gain
  * frozen stationary sd  sigma_l  (swept: mild/moderate/severe = 0.05/0.15/0.40;
    the O1 verification used 0.30)
  * frozen mean            mu = -0.5 sigma_l^2   (so E[a] = 1 exactly)
  * frozen correlation     tc_main = 2.0 frames  (lag-1 corr exp(-1/tc)=0.607)
  * measurement            Y_n ~ Poisson(a_n * s_n),   s_n = m_n . x   (s_n>0)
The OU path is simulated as the exact AR(1) the Tauchen grid of verify_o1_o4.py
discretizes:  l_n = mu + phi (l_{n-1}-mu) + sigma_l sqrt(1-phi^2) eps_n.

Schedules (lifted verbatim from verify_o1_o4.py block2 (iii)):
  * matched-pair multiset: each distinct pattern duplicated -> patt[2k]=patt[2k+1].
  * ORDERED : patterns in bank order at increasing time (the naive schedule).
  * RANDOM  : frozen uniform permutation, one per seed.
  * PAIRED  : the mirror-chop that zeroed the moment in scatter_verify --
              pair j's two copies to time slots (half-1-j) and (half+j), i.e.
              symmetric +/- t about the centre, so odd temporal moments cancel.

Patterns: the deployed SCAT32 972-row bank rows (results/round63_bridge/library/
L0.npz, frozen); scenes: 4 frozen 32x32 bridge scenes.

Read-only on all inputs; writes only this folder.  CPU, minutes.
"""
import os, json, time, platform, hashlib
from datetime import datetime, timezone
import numpy as np
from scipy.linalg import cho_factor, cho_solve

T0=time.time()
REPO=os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","..",".."))
REPO="D:/GI_another"
OUT=os.path.join(REPO,"results","round63_next","SCHEDULE_PROBE")
os.makedirs(OUT,exist_ok=True)
def log(m): print(f"[{time.time()-T0:6.1f}s] {m}",flush=True)

# ---------------- frozen inputs ----------------
LIB=os.path.join(REPO,"results","round63_bridge","library","L0.npz")
rows=np.load(LIB)["rows"].astype(np.float64)          # (972,1024) nonneg patterns
def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest()
LIB_SHA=sha(LIB)
SCENES=["bridge_contour_1","bridge_microtex_1","bridge_twopop_1","bridge_control_2"]
scene_x={}; scene_sha={}
for sc in SCENES:
    p=os.path.join(REPO,"data","r63_bridge_scenes",sc+".npz")
    scene_x[sc]=np.load(p)["x"].astype(np.float64).ravel()
    scene_sha[sc]=sha(p)

# ---------------- frozen config ----------------
Md=rows.copy()                       # 972 distinct patterns
Npix=Md.shape[1]                     # 1024
patt=np.repeat(Md,2,axis=0)          # matched-pair multiset, N=1944 rows
N=patt.shape[0]
tslots=np.linspace(-1.0,1.0,N)       # fixed time slots (verify_o1_o4 block2 iii)
DOSE=4.0                             # photon-scale (mean count ~2200/frame; shot noise subdominant to drift by design)
LAM_SCALE=0.1                        # ridge = LAM_SCALE * mean(diag(A^T A)); same for ALL arms/scenes/drift
TC=2.0                               # frozen OU correlation time (frames)
SIGMAS=[0.05,0.15,0.40]              # mild / moderate / severe stationary log-gain sd
SEEDS=list(range(5))                 # 5 paired seeds
NITER_AWARE=3
H_AWARE_ORDERS=[1,2,3,4,5,6]         # declared smooth (low-frequency) drift modes for the aware estimator (DC excluded = gauge)
TAU=1e-6                             # weak prior precision on beta

# ---------------- schedule orders ----------------
def paired_order(N):
    half=N//2; o=np.empty(N,dtype=int)
    for j in range(half):
        o[half-1-j]=2*j; o[half+j]=2*j+1
    return o
ORDERED=np.arange(N)
PAIRED =paired_order(N)
def random_order(seed): return np.random.default_rng(9000+seed).permutation(N)

# ---------------- frozen OU gain path (exact AR(1)) ----------------
def ou_path(rng,sigma_l,tc,N):
    phi=np.exp(-1.0/tc); mu=-0.5*sigma_l**2
    eps=rng.standard_normal(N); l=np.empty(N)
    l[0]=mu+sigma_l*eps[0]; sd=sigma_l*np.sqrt(max(1-phi**2,0.0))
    for n in range(1,N): l[n]=mu+phi*(l[n-1]-mu)+sd*eps[n]
    return np.exp(l)

# ---------------- moment / Schur-loss metric (H=[t,t^3], R=I; reproduces scatter_verify) ----------------
Hmet=np.stack([tslots,tslots**3],axis=1)
def schur_metric(order,x):
    Mo=patt[order]; s=Mo@x; Rinv=np.ones(N)
    Ixx=Mo.T@(Rinv[:,None]*Mo)
    IxB=((Rinv*s)[:,None]*Mo).T@Hmet
    IBB=((s*Rinv*s)[:,None]*Hmet).T@Hmet+1e-3*np.eye(2)
    IxG=Ixx-IxB@np.linalg.inv(IBB)@IxB.T
    return float(np.trace(Ixx-IxG)), float(np.abs(IxB).max())

# ---------------- estimators (R=I; ridge fixed; one Cholesky reused everywhere) ----------------
A=DOSE*patt
AtA=A.T@A
LAM=LAM_SCALE*float(np.diag(AtA).mean())
CF=cho_factor(AtA+LAM*np.eye(Npix))
def legendre_modes(t,orders):
    cols=[t**k for k in range(max(orders)+1)]
    Q,_=np.linalg.qr(np.stack(cols,1))
    return Q[:,orders]
HAW=legendre_modes(tslots,H_AWARE_ORDERS)

def recon_blind(order,Y):
    return cho_solve(CF,A[order].T@Y)
def recon_aware(order,Y,niter=NITER_AWARE):
    Ao=A[order]; x=cho_solve(CF,Ao.T@Y)
    p=HAW.shape[1]
    for _ in range(niter):
        s=np.maximum(Ao@x,1e-6)
        U=s[:,None]*HAW                 # D_s H  (N x p)
        B=Ao.T@U                        # I_xbeta   (Npix x p)
        G=U.T@U+TAU*np.eye(p)
        rhsx=Ao.T@Y; rhsb=U.T@Y
        rr=rhsx-B@np.linalg.solve(G,rhsb)   # RHS of Schur system
        Pinv_rr=cho_solve(CF,rr); Pinv_B=cho_solve(CF,B)
        Sc=G-B.T@Pinv_B
        x=Pinv_rr+Pinv_B@np.linalg.solve(Sc,B.T@Pinv_rr)   # Woodbury on (P - B G^-1 B^T)
    return x
def psnr(xh,x):  # peak=1 (scenes in [0,1])
    return 10*np.log10(1.0/np.mean((xh-x)**2))

def noiseless_ceiling(x):
    return psnr(cho_solve(CF,A.T@(A@x)),x)

# ---------------- run grid ----------------
def run_cell(sc,sigma,tc):
    x=scene_x[sc]
    arms={"ordered":ORDERED,"random":None,"paired":PAIRED}
    out={a:{"blind":[],"aware":[]} for a in arms}
    for seed in SEEDS:
        a_path=ou_path(np.random.default_rng(1000+seed),sigma,tc,N)   # shared across arms (seed-paired)
        for aname,order in arms.items():
            oo=random_order(seed) if order is None else order
            s_slot=A[oo]@x
            Y=np.random.default_rng(5000+seed).poisson(a_path*s_slot).astype(np.float64)
            out[aname]["blind"].append(psnr(recon_blind(oo,Y),x))
            out[aname]["aware"].append(psnr(recon_aware(oo,Y),x))
    return out

def summ(v):
    v=np.array(v); return dict(mean=float(v.mean()),std=float(v.std(ddof=1)),vals=[float(z) for z in v])

log("primary grid: 4 scenes x 3 sigma x 3 arms x 5 seeds x {blind,aware}, tc=%.1f"%TC)
primary={}; metrics={}; ceilings={}
for sc in SCENES:
    ceilings[sc]=float(noiseless_ceiling(scene_x[sc]))
    metrics[sc]={}
    for aname,order in [("ordered",ORDERED),("random",random_order(0)),("paired",PAIRED)]:
        sl,mm=schur_metric(order,scene_x[sc])
        metrics[sc][aname]=dict(schur_loss=sl,moment_max=mm)
    primary[sc]={}
    for sigma in SIGMAS:
        cell=run_cell(sc,sigma,TC)
        primary[sc][f"{sigma}"]={a:{est:summ(cell[a][est]) for est in ("blind","aware")} for a in cell}
        log(f"  {sc:20s} sigma={sigma:.2f}: "
            +"  ".join(f"{a}=B{np.mean(cell[a]['blind']):.2f}/A{np.mean(cell[a]['aware']):.2f}" for a in cell))

# secondary: correlation-time sensitivity (representative scene, moderate drift)
log("secondary: correlation-time sensitivity (bridge_contour_1, sigma=0.15)")
TC_SWEEP=[2.0,10.0,50.0,200.0,1000.0]
secondary={}
sc0="bridge_contour_1"
for tc in TC_SWEEP:
    cell=run_cell(sc0,0.15,tc)
    secondary[f"{tc}"]={a:{est:summ(cell[a][est]) for est in ("blind","aware")} for a in cell}
    log(f"  tc={tc:7.1f}: "+"  ".join(f"{a}=B{np.mean(cell[a]['blind']):.2f}/A{np.mean(cell[a]['aware']):.2f}" for a in cell))

meta=dict(
    utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
    python=platform.python_version(), numpy=np.__version__, platform=platform.platform(),
    runtime_s=time.time()-T0,
    library=os.path.relpath(LIB,REPO).replace("\\","/"), library_sha256=LIB_SHA,
    scenes={sc:scene_sha[sc] for sc in SCENES},
    config=dict(N_frames=N, N_distinct=Md.shape[0], Npix=Npix, dose=DOSE,
                lam_scale=LAM_SCALE, lam_abs=LAM, tc_main=TC, sigmas=SIGMAS,
                mu_rule="-0.5*sigma_l^2 (E[a]=1)", seeds=SEEDS,
                H_metric="[t, t^3] (odd AC; R=I) -- reproduces scatter_verify block2(iii)",
                H_aware_orders=H_AWARE_ORDERS, niter_aware=NITER_AWARE, tau=TAU,
                peak_for_psnr=1.0),
    frozen_source="code/scatter/verify_o1_o4.py (commit fb3feda)",
    ceilings_noiseless_dB=ceilings,
)
RES=dict(meta=meta,moment_metrics=metrics,primary=primary,secondary=secondary,tc_sweep=TC_SWEEP)
jp=os.path.join(OUT,"probe_b_results.json")
with open(jp,"w",encoding="utf-8") as f: json.dump(RES,f,indent=2)
log(f"wrote {jp}")
log(f"TOTAL runtime {time.time()-T0:.1f}s")
