#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""COMSOL 6.3 (Wave Optics) micro-validation of the thin-phase-screen approximation.

2D TE (E_z out-of-plane), scattered-field formulation.  A dielectric slab (n=1.5, mean
thickness d0) with a GAUSSIAN-CORRELATED ROUGH TOP surface z=d0+h(x) is a physical phase
screen.  Normally-incident plane wave -> full-Maxwell transmitted field read on a line
above.  Compare speckle (contrast, grain) to the scalar thin-screen prediction with the
SAME h(x).  Graded-index single-rectangle build (robust headless); scattering BCs.

Env: N_REAL (default 6), LC (um, 5), SPHI (rad, 2pi developed), CORES (3).
"""
import os, time, json
import numpy as np

t0 = time.time()
LAMB = 0.532e-6; ND = 1.5
N_REAL = int(os.environ.get("N_REAL", "6"))
L_C_UM = float(os.environ.get("LC", "5.0"))
SIGMA_PHI = float(os.environ.get("SPHI", str(2*np.pi)))
D0_UM = 3.0; WIN_X_UM = 60.0; Z_TOP_UM = 12.0; Z_READ_UM = 8.0
MESH_LAMB_FRAC = 8.0
OUT = os.path.dirname(os.path.abspath(__file__))
sig_h_um = SIGMA_PHI*(LAMB*1e6)/(2*np.pi*(ND-1))

def make_surface(seed, nx=1024):
    rng = np.random.default_rng(seed)
    x = np.linspace(-WIN_X_UM/2, WIN_X_UM/2, nx); dx = x[1]-x[0]
    f = np.fft.fftfreq(nx, d=dx)
    Hk = np.exp(-2*np.pi**2*(L_C_UM**2)*f**2)
    h = np.fft.ifft(np.fft.fft(rng.standard_normal(nx))*np.sqrt(Hk)).real
    h = (h-h.mean())/(h.std()+1e-12)*sig_h_um
    return x, h

def thin_screen_pred(x_um, h_um, z_read_um, nx_out=512):
    nx=len(x_um); dx=(x_um[1]-x_um[0])*1e-6
    phi=(2*np.pi/LAMB)*(ND-1)*(h_um*1e-6)
    E=np.exp(1j*phi); f=np.fft.fftfreq(nx,d=dx)
    kz=2*np.pi*np.sqrt(np.maximum((1/LAMB)**2-f**2,0)).astype(complex)
    Eo=np.fft.ifft(np.fft.fft(E)*np.exp(1j*kz*z_read_um*1e-6))
    I=np.abs(Eo)**2
    xo=np.linspace(x_um[0],x_um[-1],nx_out)
    return xo, np.interp(xo, x_um, I)

def speckle_stats(x_um, I, central_frac=0.6):
    n=len(I); a=int(n*(1-central_frac)/2); b=n-a
    Ic=I[a:b]; xc=x_um[a:b]
    contrast=float(Ic.std()/(Ic.mean()+1e-30))
    dI=Ic-Ic.mean(); ac=np.correlate(dI,dI,'full'); ac=ac[ac.size//2:]; ac=ac/(ac[0]+1e-30)
    dx=xc[1]-xc[0]; below=np.where(ac<0.5)[0]
    grain=float(2*below[0]*dx) if len(below) else float('nan')
    return contrast, grain

import mph
client = mph.start(cores=int(os.environ.get("CORES","3")))
print("COMSOL", client.version, "started", flush=True)

def solve_one(seed):
    x_um,h_um = make_surface(seed)
    fpath=os.path.join(OUT,f"surf_{seed}.txt").replace("\\","/")
    np.savetxt(fpath, np.c_[x_um*1e-6, h_um*1e-6])
    model=client.create(f"m{seed}"); m=model.java
    m.param().set("lam",f"{LAMB}[m]"); m.param().set("nd",str(ND)); m.param().set("d0",f"{D0_UM}[um]")
    m.component().create("comp1",True)
    g=m.component("comp1").geom().create("geom1",2); g.lengthUnit("um")
    r=g.create("r1","Rectangle"); r.set("size",[WIN_X_UM, Z_TOP_UM+D0_UM+6.0]); r.set("pos",[-WIN_X_UM/2,-6.0]); g.run()
    fn=m.component("comp1").func().create("hint","Interpolation")
    fn.set("source","file"); fn.set("filename",fpath); fn.set("interp","linear"); fn.set("extrap","const")
    fn.setIndex("funcs","hx",0,0); fn.importData()
    # NOTE: COMSOL 2D in-plane coordinates are (x, y); y is the propagation/vertical axis.
    v=m.component("comp1").variable().create("var1")
    v.set("stepw","0.15[um]")
    v.set("nreal","1+(nd-1)*(0.5*(1+tanh(y/stepw)))*(0.5*(1+tanh((d0+hx(x)-y)/stepw)))")
    v.set("epsr","nreal^2")
    # material with spatially-varying relative permittivity (ewfd default constitutive relation)
    mat=m.component("comp1").material().create("mat1","Common"); mat.selection().all()
    mat.propertyGroup("def").set("relpermittivity",["epsr","0","0","0","epsr","0","0","0","epsr"])
    phys=m.component("comp1").physics().create("ewfd","ElectromagneticWavesFrequencyDomain","geom1")
    phys.prop("components").set("components","outofplane")   # Ez only (2D TE)
    phys.prop("BackgroundField").set("SolveFor","scatteredField")
    phys.prop("BackgroundField").set("Eb",["0","0","exp(-i*ewfd.k0*y)"])   # plane wave +y
    sbc=phys.create("sbc1","Scattering",1); sbc.selection().all()   # single rectangle: all 4 edges exterior
    mesh=m.component("comp1").mesh().create("mesh1")
    ftri=mesh.create("ftri","FreeTri")
    ftri.create("sz","Size").set("hmax",f"{LAMB*1e6/MESH_LAMB_FRAC}[um]")   # size subnode (unique tag)
    mesh.run()
    std=m.study().create("std1"); std.create("freq","Frequency")
    std.feature("freq").set("plist","c_const/lam"); std.run()
    # evaluate field + coords on the solution mesh nodes, extract a thin y-band at y=z_read
    Ival, xval, yval = model.evaluate(["abs(ewfd.Ez)^2","x","y"])
    Ival=np.asarray(Ival).ravel(); xval=np.asarray(xval).ravel(); yval=np.asarray(yval).ravel()
    scale = 1e6 if np.abs(yval).max() < 1e-2 else 1.0      # -> micrometres (auto units)
    xval=xval*scale; yval=yval*scale
    tol=0.4
    band = np.abs(yval - Z_READ_UM) < tol
    while band.sum() < 50 and tol < 3.0:
        tol*=1.5; band = np.abs(yval - Z_READ_UM) < tol
    xb, Ib = xval[band], Ival[band]
    order=np.argsort(xb); xb, Ib = xb[order], Ib[order]
    xs=np.linspace(-WIN_X_UM/2*0.9, WIN_X_UM/2*0.9, 512)   # central 90%
    I=np.interp(xs, xb, Ib)
    client.remove(f"m{seed}")
    return x_um,h_um,xs,I

res={"note":"COMSOL 6.3 Wave Optics 2D micro-validation of thin-phase-screen approx.",
     "lambda_nm":532,"n_d":ND,"l_c_um":L_C_UM,"sigma_phi":SIGMA_PHI,"sigma_h_um":sig_h_um,
     "d0_um":D0_UM,"win_x_um":WIN_X_UM,"z_read_um":Z_READ_UM,"mesh_lambda_frac":MESH_LAMB_FRAC,
     "n_real":N_REAL,"rows":[]}
fem_c=[]; fem_g=[]; ts_c=[]; ts_g=[]
for seed in range(N_REAL):
    try:
        x_um,h_um,xs,I = solve_one(seed)
        cF,gF = speckle_stats(xs, I)
        xo,Its = thin_screen_pred(x_um,h_um,Z_READ_UM); cT,gT = speckle_stats(xo,Its)
        fem_c.append(cF); fem_g.append(gF); ts_c.append(cT); ts_g.append(gT)
        res["rows"].append(dict(seed=seed, fem_contrast=round(cF,4), fem_grain_um=round(gF,3),
                                thinscreen_contrast=round(cT,4), thinscreen_grain_um=round(gT,3)))
        print(f"[seed {seed}] FEM contrast={cF:.3f} grain={gF:.2f}um | thin-screen contrast={cT:.3f} grain={gT:.2f}um  ({time.time()-t0:.0f}s)", flush=True)
    except Exception as e:
        print(f"[seed {seed}] ERROR {repr(e)[:300]}", flush=True)
        res["rows"].append(dict(seed=seed, error=repr(e)[:300]))

if fem_c:
    res["summary"]=dict(
        fem_contrast_mean=round(float(np.mean(fem_c)),4), thinscreen_contrast_mean=round(float(np.mean(ts_c)),4),
        fem_grain_mean_um=round(float(np.nanmean(fem_g)),3), thinscreen_grain_mean_um=round(float(np.nanmean(ts_g)),3),
        contrast_rel_diff=round(float(abs(np.mean(fem_c)-np.mean(ts_c))/(np.mean(ts_c)+1e-9)),4),
        grain_rel_diff=round(float(abs(np.nanmean(fem_g)-np.nanmean(ts_g))/(np.nanmean(ts_g)+1e-9)),4))
    print("SUMMARY", json.dumps(res["summary"]), flush=True)
json.dump(res, open(os.path.join(OUT,"COMSOL_MICRO_RESULTS.json"),"w"), indent=2)
print(f"saved COMSOL_MICRO_RESULTS.json  elapsed {time.time()-t0:.0f}s", flush=True)
