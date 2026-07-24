# GATE G3 -- SCALE + NATURAL IMAGES, shell-resolved band-extension DEPTH (R39 + coordinator spec).
#
# Ports the beyond-band cell to 64x64 (N=4096): pattern band fx,fy<=8 (M=128 band-limited nonneg
# patterns, M_eff~80), medium band annulus 1<=i+j<=16 (db~152, reach ~24).
#
# The beyond-band claim is a BAND-EXTENSION DEPTH, not an aggregate: recovery is SHELL-RESOLVED in
# Delta-k = max(i,j)-PB (Chebyshev distance beyond the pattern box). At the 16x16 toy (moment
# matching) the coordinator measured Dk=1->0.318, Dk=2->0.91, Dk=3->1.17, Dk=4->2.1 (deep shells
# actively corrupt), i.e. a viable depth ~1. G3 reports, at 64x64:
#   PART A  aperture-law in/out separation (does the law hold at scale).
#   PART B  per-shell beyond-band err(Delta-k=1..6) under FLAT and DECAYING (power~k^-2) media,
#           with moment matching (apples-to-apples vs the toy) and the two-stage MLE; the key
#           number is the VIABLE EXTENSION DEPTH (largest Dk with median err<=0.3) and whether it
#           GROWS vs the toy's depth of 1. Production form CLIPS beyond the viable depth.
# KILL: law separation < 3x, OR viable depth = 0 under BOTH spectra at generous budget (washout).
import json
import sys
import time

import numpy as np
import torch
from scipy.fftpack import dct, idct

sys.path.insert(0, "../FOG_DMD_PROBE")
import fog_tracker as ft  # noqa: E402

DEV = "cuda" if torch.cuda.is_available() else "cpu"
n = 64
N = n * n
PB = 8
MED_LO, MED_HI = 1, 16
SIG_F = 0.30
TAU = 8.0
PHOT = 1e5
RHO = float(np.exp(-1.0 / TAU))
M = 128
T_GEN = 4096
T_MLE = 1024                 # MLE cell: RTS smoother stores (T,d,d) -> memory-bounded at 64x64
N_STARTS = 3
STEPS = 4000
t0 = time.time()


def D2(a):
    return dct(dct(a.reshape(n, n), axis=0, norm="ortho"), axis=1, norm="ortho")


def I2(A):
    return idct(idct(A, axis=0, norm="ortho"), axis=1, norm="ortho").ravel()


def spike(i, j):
    g = np.zeros((n, n)); g[i, j] = 1.0
    return I2(g)


band = [(i, j) for i in range(n) for j in range(n) if MED_LO <= i + j <= MED_HI]
Ub = np.linalg.qr(np.array([spike(i, j) for (i, j) in band]).T)[0]
db = Ub.shape[1]
coeff_sd = SIG_F * np.sqrt(N / db)

# medium per-mode sd for FLAT and DECAYING (power ~ k^-2 => sd ~ k^-1) spectra, RMS-normalized
_r = np.array([np.hypot(i, j) for (i, j) in band])
_dec = 1.0 / np.maximum(_r, 1.0)
GEN_SD = {"flat": None,
          "decaying": coeff_sd * _dec / np.sqrt((_dec ** 2).mean())}

pat_box = np.zeros((n, n), bool); pat_box[:PB + 1, :PB + 1] = True
cov = np.zeros((n, n), bool)
for pi in range(PB + 1):
    for pj in range(PB + 1):
        for (ki, kj) in band:
            for si in (pi + ki, abs(pi - ki)):
                for sj in (pj + kj, abs(pj - kj)):
                    if si < n and sj < n:
                        cov[si, sj] = True
superres = cov & (~pat_box)

# Delta-k shell masks (Chebyshev distance beyond the pattern box), intersected with coverage
SHELLS = {}
for dk in range(1, 7):
    m = PB + dk
    mask = np.zeros((n, n), bool)
    for i in range(n):
        for j in range(n):
            if max(i, j) == m and cov[i, j]:
                mask[i, j] = True
    SHELLS[dk] = mask
print(f"G3 SCALE 64x64: N={N} PB={PB} medium {MED_LO}<=i+j<={MED_HI} db={db} M={M} T={T_GEN}", flush=True)
print(f"  coverage frac {cov.mean():.3f}  super-res-zone frac {superres.mean():.3f}", flush=True)
print(f"  shell sizes (freqs) Dk=1..6: {[int(SHELLS[k].sum()) for k in range(1,7)]}", flush=True)


def bl_patterns(M_, seed):
    rp = np.random.default_rng(seed); ps = []
    for _ in range(M_):
        Cp = np.zeros((n, n)); Cp[:PB + 1, :PB + 1] = rp.standard_normal((PB + 1, PB + 1)); Cp[0, 0] = 0
        f = I2(Cp); f /= np.abs(f).max(); ps.append(0.5 + 0.45 * f)
    return np.array(ps)


Pfix = bl_patterns(M, 10)


def medium_field(T, seed, gen_sd=None):
    r = np.random.default_rng(seed)
    sd = coeff_sd if gen_sd is None else np.asarray(gen_sd, float)
    z = sd * r.standard_normal(db)
    var_n = ((sd ** 2) * (Ub ** 2)).sum(axis=1) if gen_sd is not None \
        else (coeff_sd ** 2) * (Ub ** 2).sum(axis=1)
    W = np.zeros((T, N))
    for t in range(T):
        W[t] = np.exp(Ub @ z - 0.5 * var_n)
        z = RHO * z + np.sqrt(1 - RHO ** 2) * sd * r.standard_normal(db)
    return W


def buckets(P, W, x, seed, phot=PHOT):
    r = np.random.default_rng(seed)
    mu = np.einsum("mn,tn->tm", P, W * x)
    scp = phot / mu.mean()
    B = r.poisson(mu * scp) / scp
    R_diag = np.clip(mu, 1e-9, None) / scp
    return B, R_diag


def shell_err(xh, x, mask):
    s = np.dot(xh, x) / max(np.dot(xh, xh), 1e-12)
    E = D2(s * xh - x) ** 2; X0 = D2(x) ** 2
    d = X0[mask].sum()
    return float(E[mask].sum() / d) if d > 1e-12 else float("nan")


def moment_estimate(B, seed0=0):
    rc = ft.cov_estimate(Pfix, Ub, B, coeff_sd, RHO, n_starts=N_STARTS, steps=STEPS,
                         dev=DEV, seed0=seed0)
    return rc["per_start_x"][int(np.argmin(rc["objs"]))]


def viable_depth(shell_med, thr):
    """Largest contiguous Dk (from 1) with median err <= thr."""
    d = 0
    for dk in range(1, 7):
        if np.isfinite(shell_med[dk]) and shell_med[dk] <= thr:
            d = dk
        else:
            break
    return d


# witness with content at each shell (known per-shell energy) + in-band random
def scene_shell_witness(seed=4, per_shell=10, amp=2.5):
    C = np.zeros((n, n)); rs = np.random.default_rng(seed)
    C[:PB + 1, :PB + 1] = rs.standard_normal((PB + 1, PB + 1))
    for dk in range(1, 7):
        fi, fj = np.where(SHELLS[dk])
        if len(fi) == 0:
            continue
        pick = rs.choice(len(fi), size=min(per_shell, len(fi)), replace=False)
        for k in pick:
            C[fi[k], fj[k]] = amp * rs.choice([-1, 1])
    x = I2(C); x = x - x.min(); return x / x.max()


def scene_broadband(seed=7):
    rs = np.random.default_rng(seed)
    C = rs.standard_normal((n, n))
    fy, fx = np.mgrid[0:n, 0:n]
    C = C / (1.0 + np.hypot(fy, fx) / 8.0)
    x = I2(C); x = x - x.min(); return x / x.max()


def scene_natural(name):
    from skimage import data, transform
    src = {"cameraman": data.camera, "coins": data.coins, "moon": data.moon}[name]()
    img = src.astype(np.float64)
    img = transform.resize(img, (n, n), anti_aliasing=True, mode="reflect")
    return ((img - img.min()) / (np.ptp(img) + 1e-12)).ravel()


# =================== PART A: aperture-law in/out separation ===================
print("\n[A] aperture-law in/out separation (full-spectrum broadband scene)", flush=True)
xbb = scene_broadband()
ein, eout = [], []
for s in range(3):
    W = medium_field(T_GEN, 500 + s)
    B, _ = buckets(Pfix, W, xbb, 550 + s)
    xh = moment_estimate(B, seed0=10 * s)
    ein.append(shell_err(xh, xbb, cov & (D2(xbb) ** 2 > 1e-9)))
    eout.append(shell_err(xh, xbb, (~cov) & (D2(xbb) ** 2 > 1e-9)))
in_med, out_med = float(np.median(ein)), float(np.median(eout))
sep = out_med / max(in_med, 1e-9)
print(f"    in-coverage {in_med:.3f}  out-coverage {out_med:.3f}  SEPARATION {sep:.2f}x  (KILL <3x)", flush=True)

# =================== PART B: shell-resolved band-extension depth ===============
print("\n[B] shell-resolved beyond-band err(Delta-k) and viable extension depth", flush=True)
xw = scene_shell_witness()
sr_by_shell = {dk: float((D2(xw) ** 2)[SHELLS[dk]].sum()) for dk in range(1, 7)}
print(f"    witness per-shell energy Dk=1..6: {[round(sr_by_shell[k], 2) for k in range(1, 7)]}", flush=True)
partB = {}
for spec, gsd in GEN_SD.items():
    # moment matching (apples-to-apples vs the 16x16 toy)
    per_seed = {dk: [] for dk in range(1, 7)}
    for s in range(3):
        W = medium_field(T_GEN, 600 + 10 * s, gen_sd=gsd)
        B, _ = buckets(Pfix, W, xw, 650 + 10 * s)
        xh = moment_estimate(B, seed0=100 + 10 * s)
        for dk in range(1, 7):
            per_seed[dk].append(shell_err(xh, xw, SHELLS[dk]))
    shell_med = {dk: float(np.median(per_seed[dk])) for dk in range(1, 7)}
    depth03 = viable_depth(shell_med, 0.3); depth05 = viable_depth(shell_med, 0.5)
    partB[spec] = dict(estimator="moment", shell_median=shell_med,
                       shell_all={dk: per_seed[dk] for dk in range(1, 7)},
                       viable_depth_0p3=depth03, viable_depth_0p5=depth05)
    print(f"    [{spec:8s} moment] err Dk=1..6: "
          f"{[round(shell_med[k], 3) for k in range(1, 7)]}  depth(<=0.3)={depth03} depth(<=0.5)={depth05}",
          flush=True)

# two-stage MLE at T_MLE (flat spectrum) -- does the efficient recipe deepen the viable depth?
mle_block = None
try:
    W = medium_field(T_MLE, 700, gen_sd=None)
    B, R_diag = buckets(Pfix, W, xw, 750)
    tm = time.time()
    xm = moment_estimate(B, seed0=7)
    xL = ft.kalman_em(Pfix, Ub, B, R_diag, RHO, coeff_sd, 1e-4, n_em=15, dev=DEV, x_init=xm)
    shell_med_mle = {dk: shell_err(xL, xw, SHELLS[dk]) for dk in range(1, 7)}
    shell_med_mom = {dk: shell_err(xm, xw, SHELLS[dk]) for dk in range(1, 7)}
    mle_block = dict(T=T_MLE, mle_shell=shell_med_mle, moment_shell_sameT=shell_med_mom,
                     mle_depth_0p3=viable_depth(shell_med_mle, 0.3),
                     moment_depth_0p3=viable_depth(shell_med_mom, 0.3), wall_s=time.time() - tm)
    print(f"    [flat MLE T={T_MLE}] err Dk=1..6: {[round(shell_med_mle[k], 3) for k in range(1, 7)]}  "
          f"depth(<=0.3)={mle_block['mle_depth_0p3']} (moment sameT depth {mle_block['moment_depth_0p3']}) "
          f"[{mle_block['wall_s']:.0f}s]", flush=True)
except Exception as e:
    mle_block = dict(error=str(e))
    print(f"    [MLE] deferred at 64x64 (T-loop cost/memory): {e}", flush=True)

# =================== natural images (aggregate + shell-1/2 recovery) ==========
print("\n[B2] natural images: super-res-zone + shallow-shell recovery", flush=True)
nat = {}
for nm in ["cameraman", "coins", "moon"]:
    try:
        xs = scene_natural(nm)
    except Exception as e:
        print(f"    (skip {nm}: {e})", flush=True); continue
    W = medium_field(T_GEN, 800)
    B, _ = buckets(Pfix, W, xs, 850)
    xh = moment_estimate(B, seed0=200)
    sh = {dk: shell_err(xh, xs, SHELLS[dk]) for dk in (1, 2, 3)}
    agg = shell_err(xh, xs, superres)
    nat[nm] = dict(superres_nmse=agg, shell1=sh[1], shell2=sh[2], shell3=sh[3],
                   superres_energy_frac=float((D2(xs) ** 2)[superres].sum() / (D2(xs) ** 2).sum()))
    print(f"    {nm:11s} superres-NMSE {agg:.3f}  Dk1 {sh[1]:.3f}  Dk2 {sh[2]:.3f}  Dk3 {sh[3]:.3f}", flush=True)

# =================== verdict ==================================================
depth_flat = partB["flat"]["viable_depth_0p3"]
depth_dec = partB["decaying"]["viable_depth_0p3"]
kill_law = sep < 3.0
kill_depth = (depth_flat == 0 and depth_dec == 0)
verdict = "KILL" if (kill_law or kill_depth) else ("WATCH" if (depth_flat <= 1 and depth_dec == 0) else "PASS")
out = dict(gate="G3_scale_shell_depth", n=n, N=N, PB=PB, medium_band=[MED_LO, MED_HI], db=db, M=M,
           T=T_GEN, coverage_frac=float(cov.mean()), superres_frac=float(superres.mean()),
           shell_sizes={dk: int(SHELLS[dk].sum()) for dk in range(1, 7)},
           partA=dict(in_coverage_err=in_med, out_coverage_err=out_med, separation=sep,
                      in_all=ein, out_all=eout),
           partB=partB, mle=mle_block, naturals=nat,
           toy_depth_reference=1, viable_depth_flat=depth_flat, viable_depth_decaying=depth_dec,
           kill_rule="separation<3x OR viable depth=0 under BOTH spectra at generous budget",
           kill_law=kill_law, kill_depth=kill_depth, verdict=verdict, wall_s=time.time() - t0)
json.dump(out, open("G3_results.json", "w"), indent=2)
print(f"\nG3 VERDICT: {verdict}  (law {sep:.1f}x; viable depth flat={depth_flat} decaying={depth_dec}; "
      f"toy depth=1)  [{time.time()-t0:.0f}s]", flush=True)
