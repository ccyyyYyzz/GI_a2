# GATE G2 -- GEOMETRY MIXING (multiplicative thin-screen vs partial convolutive propagation).
#
# The multiplicative thin-screen model b=<P (.) w, x> is the DECLARED home (shower-glass regime,
# screen AT the object plane). At nonzero propagation distance the effective illumination is a
# CONVOLUTION (Fresnel) of the screen-modulated field, not a pure product. We contaminate the
# TRUE forward model with partial convolutive mixing and keep the PURE multiplicative estimator:
#
#   primary  (effective-illumination blend):  A_t = (1-a)*(P (.) w_t) + a*blur(P (.) w_t)
#   variant  (convolutive medium)          :  b   = <P (.) (w_t * h_a), x>   (medium pre-blurred)
#
# for a in {0, 0.1, 0.25, 0.5}, h = small Gaussian PSF (emulates propagation distance).
# Estimator keeps the declared fine band Ub and the pure multiplicative model. Report beyond-band
# hi_err vs alpha. KILL SIGNAL: collapse (>0.9) already at alpha=0.1.
import json
import time

import numpy as np
from scipy.ndimage import gaussian_filter

import bb_common as bb

T = 2048
N_SEED = 3
STEPS = 4000
N_STARTS = 3
BLUR_SIGMA = 0.8            # small propagation PSF (pixels)
t0 = time.time()

Ub, db, coeff_sd, Kg = bb.medium_subspace()
x, hi_mask = bb.make_scene()
Pfix = bb.bl_patterns(bb.M_PAT, 10)
n = bb.N_SIDE


def blur_rows(mat):
    """Spatial Gaussian blur applied to each (M,N) row reshaped to n x n (a small propagation PSF).
    Vectorized: blur only the two spatial axes of the (M,n,n) stack (sigma=0 on the pattern axis)."""
    stack = mat.reshape(mat.shape[0], n, n)
    out = gaussian_filter(stack, sigma=(0, BLUR_SIGMA, BLUR_SIGMA), mode="reflect")
    return out.reshape(mat.shape[0], -1)


def buckets_blend(W, alpha, seed, phot=bb.PHOT):
    """TRUE forward = effective-illumination blend; estimator will assume alpha=0.
    Vectorized over t: PW=(T,M,N), blur spatial axes only, A=blend, mu=<A,x>."""
    r = np.random.default_rng(seed)
    PW = Pfix[None, :, :] * W[:, None, :]               # (T,M,N)
    if alpha > 0:
        PWb = gaussian_filter(PW.reshape(T, Pfix.shape[0], n, n),
                              sigma=(0, 0, BLUR_SIGMA, BLUR_SIGMA),
                              mode="reflect").reshape(T, Pfix.shape[0], -1)
        A = (1 - alpha) * PW + alpha * PWb
    else:
        A = PW
    mu = np.einsum("tmn,n->tm", A, x)
    scp = phot / mu.mean()
    B = r.poisson(np.clip(mu, 0, None) * scp) / scp
    return B


def buckets_convmed(W, alpha, seed, phot=bb.PHOT):
    """TRUE forward = pre-blurred medium w_t*h_a blended in: w_eff=(1-a)w+a*blur(w)."""
    r = np.random.default_rng(seed)
    Wb = np.empty_like(W)
    for t in range(T):
        Wb[t] = gaussian_filter(W[t].reshape(n, n), BLUR_SIGMA, mode="reflect").ravel()
    Weff = (1 - alpha) * W + alpha * Wb
    mu = np.einsum("mn,tn->tm", Pfix, Weff * x)
    scp = phot / mu.mean()
    B = r.poisson(np.clip(mu, 0, None) * scp) / scp
    return B


def sweep(bucket_fn, tag):
    rows = []
    base = None
    for alpha in (0.0, 0.1, 0.25, 0.5):
        errs = []
        for s in range(N_SEED):
            W = bb.medium_field(T, 600 + s, Ub, coeff_sd)
            B = bucket_fn(W, alpha, 650 + s)
            xm, _ = bb.est_moment(Pfix, Ub, B, coeff_sd, n_starts=N_STARTS, steps=STEPS)
            errs.append(bb.hi_err(xm, x, hi_mask))
        m = float(np.median(errs))
        if alpha == 0.0:
            base = m
        rows.append(dict(alpha=alpha, err=m, all=errs,
                         degr=(m - base) / base if base else 0.0))
        print(f"  [{tag}] alpha={alpha:.2f}  err {m:.3f}  degr {((m-base)/base)*100:+.0f}%"
              f"  {[round(q,3) for q in errs]}", flush=True)
    return rows


print(f"G2 GEOMETRY MIXING  T={T} seeds={N_SEED} blur_sigma={BLUR_SIGMA}px", flush=True)
print("primary: effective-illumination blend  A=(1-a)(P.w)+a*blur(P.w)", flush=True)
rows_blend = sweep(buckets_blend, "blend")
print("variant: convolutive medium  w_eff=(1-a)w+a*blur(w)", flush=True)
rows_conv = sweep(buckets_convmed, "convmed")

# verdict: kill if collapse (>0.9) already at alpha=0.1 in the primary blend
a01 = [r for r in rows_blend if r["alpha"] == 0.1][0]["err"]
kill = a01 > 0.9
# graceful-degradation read across the sweep
worst = max(r["err"] for r in rows_blend)
verdict = "KILL" if kill else ("WATCH" if worst > 0.9 else "PASS")
out = dict(gate="G2_geometry_mixing", T=T, seeds=N_SEED, blur_sigma=BLUR_SIGMA,
           primary_blend=rows_blend, variant_convmed=rows_conv,
           err_at_alpha_0p1=a01, worst_err=worst,
           kill_rule="collapse (>0.9) already at alpha=0.1 (primary blend)",
           verdict=verdict, wall_s=time.time() - t0)
json.dump(out, open("G2_results.json", "w"), indent=2)
print(f"\nG2 VERDICT: {verdict}  (blend err@alpha=0.1 = {a01:.3f}; worst-across-sweep {worst:.3f})"
      f"  [{time.time()-t0:.0f}s]", flush=True)
