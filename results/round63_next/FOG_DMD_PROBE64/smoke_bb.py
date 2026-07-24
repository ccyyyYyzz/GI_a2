# Smoke test: reproduce the R39 beyond-band matched baseline (fixed+cov moment 0.47-0.79),
# the fresh+mean physics wall (1.000), and the oracle ceiling, to validate bb_common.
import time
import numpy as np
import bb_common as bb

t0 = time.time()
Ub, db, coeff_sd, Kg = bb.medium_subspace()
x, hi_mask = bb.make_scene()
Pfix = bb.bl_patterns(bb.M_PAT, 10)
print(f"cell: n={bb.N_SIDE} N={bb.N} PB={bb.PB} medium band {bb.MED_LO}<=i+j<={bb.MED_HI} "
      f"db={db} coeff_sd={coeff_sd:.3f} M={bb.M_PAT} hi-freqs={int(hi_mask.sum())}", flush=True)
print(f"beyond-band scene energy fraction (DCT power outside pattern band): "
      f"{(bb._D2(x)**2)[hi_mask].sum()/(bb._D2(x)**2).sum():.3f}", flush=True)

T = 2048
mom, orc, wall = [], [], []
for s in range(3):
    W = bb.medium_field(T, 600 + s, Ub, coeff_sd)
    B, Rd, scp = bb.buckets(Pfix, W, x, 650 + s)
    xm, _ = bb.est_moment(Pfix, Ub, B, coeff_sd, n_starts=3, steps=4000)
    mom.append(bb.hi_err(xm, x, hi_mask))
    xo = bb.est_oracle(Pfix, W, B)
    orc.append(bb.hi_err(xo, x, hi_mask))
xw = bb.fresh_wall(10, Ub, coeff_sd, x, T, 0)
wall.append(bb.hi_err(xw, x, hi_mask))

print(f"\nT={T}  (3 seeds)", flush=True)
print(f"  fixed+cov moment  beyond-band err: median {np.median(mom):.3f}  {[round(q,3) for q in mom]}", flush=True)
print(f"  oracle (medium known)            : median {np.median(orc):.3f}  {[round(q,3) for q in orc]}", flush=True)
print(f"  fresh+mean physics wall          : {wall[0]:.3f}  (expect ~1.000)", flush=True)
print(f"\nwall clock {time.time()-t0:.0f}s", flush=True)
