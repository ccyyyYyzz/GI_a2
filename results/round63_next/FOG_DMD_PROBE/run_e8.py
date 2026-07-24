# E8 -- BEYOND-MODULATOR-BAND (super-resolution) recovery, F4-immune. Pattern band-limited to
# the DMD pixel limit (fx,fy<=PB); medium band EXCEEDS it (fine speckle). Fresh band-limited
# patterns physically CANNOT recover beyond-band content (rel err = 1.000, the physics wall).
# The fixed-bank fluctuation route can. Reruns the coordinator's exact config with the E7a
# marginal-likelihood MLE (moment-matching init -> Kalman EM), + oracle ceiling.
import time, json, numpy as np, torch
from scipy.fftpack import idct, dct
import fog_tracker as ft
def log(*a): print(*a, flush=True)
dev = 'cuda' if torch.cuda.is_available() else 'cpu'
n = 16; N = n * n
D2 = lambda a: dct(dct(a.reshape(n, n), axis=0, norm='ortho'), axis=1, norm='ortho')
I2 = lambda A: idct(idct(A, axis=0, norm='ortho'), axis=1, norm='ortho').ravel()
SIG_F = 0.30; TAU = 8.0; PHOT = 1e5; rho = np.exp(-1 / TAU); PB = 3
# medium band 1<=i+j<=6 (beyond the pattern band PB=3)
band = [(i, j) for i in range(7) for j in range(7) if 1 <= i + j <= 6]
mm = []
for (i, j) in band:
    g = np.zeros((n, n)); g[i, j] = 1.0; mm.append(I2(g))
Ub = np.linalg.qr(np.array(mm).T)[0]
db = Ub.shape[1]
sc = SIG_F * np.sqrt(N / db)
Kg = (SIG_F ** 2 * N / db) * (Ub @ Ub.T)
# scene: in-band random + deliberate beyond-band content
C = np.zeros((n, n)); rs = np.random.default_rng(4)
C[:PB + 1, :PB + 1] = rs.standard_normal((PB + 1, PB + 1))
hi = [(5, 2), (2, 5), (4, 4), (6, 1), (1, 6), (5, 5), (7, 0), (0, 7)]
for (i, j) in hi: C[i, j] = 2.0 * rs.choice([-1, 1])
x = I2(C); x = x - x.min(); x = x / x.max()
hi_mask = np.zeros((n, n), bool)
for (i, j) in hi: hi_mask[i, j] = True

def bl_patterns(M_, seed):
    rp = np.random.default_rng(seed); ps = []
    for _ in range(M_):
        Cp = np.zeros((n, n)); Cp[:PB + 1, :PB + 1] = rp.standard_normal((PB + 1, PB + 1)); Cp[0, 0] = 0
        f = I2(Cp); f /= np.abs(f).max(); ps.append(0.5 + 0.45 * f)
    return np.array(ps)

M_ = 24
Pfix = bl_patterns(M_, 10)

def med(T, seed):
    r = np.random.default_rng(seed); z = sc * r.standard_normal(db); W = np.zeros((T, N))
    for t in range(T):
        W[t] = np.exp(Ub @ z - 0.5 * SIG_F ** 2); z = rho * z + np.sqrt(1 - rho ** 2) * sc * r.standard_normal(db)
    return W

def hi_err(xh):
    s = np.dot(xh, x) / max(np.dot(xh, xh), 1e-12); E = D2(s * xh - x) ** 2; X0 = D2(x) ** 2
    return float(E[hi_mask].sum() / X0[hi_mask].sum())

def fresh_wall(W, T, seed):
    AtA = np.zeros((N, N)); Atb = np.zeros(N); rb = np.random.default_rng(680 + seed)
    scp = PHOT / np.einsum('mn,tn->tm', Pfix, W * x).mean()
    for t in range(T):
        Pf = bl_patterns(M_, 7000 + 1000 * seed + t)
        mu2 = Pf @ (W[t] * x); B2 = rb.poisson(mu2 * scp) / scp
        AtA += Pf.T @ Pf; Atb += Pf.T @ B2
    lam = 1e-3 * np.trace(AtA) / N
    return hi_err(np.clip(np.linalg.solve(AtA + lam * np.eye(N), Atb), 0, None))

results = []; t0 = time.time()
log(f"E8 beyond-band super-resolution: pattern band<={PB}, medium band<=6, db={db}, hi-freqs={len(hi)}")
for T in [512, 1024, 2048]:
    mo, ml, orc = [], [], []
    for s in range(3):
        W = med(T, 600 + s); r = np.random.default_rng(650 + s)
        mu = np.einsum('mn,tn->tm', Pfix, W * x); scp = PHOT / mu.mean()
        B = r.poisson(mu * scp) / scp; Rd = np.clip(mu, 1e-9, None) / scp
        # oracle: medium known
        Wt = torch.tensor(W, device=dev, dtype=torch.float64)
        xo = ft.solve_oracle(Pfix, W, B, np.ones_like(B), 1e-6, dev); orc.append(hi_err(xo))
        # moment matching (E5a-style)
        rc = ft.cov_estimate(Pfix, Ub, B, sc, rho, n_starts=3, steps=4000, dev=dev)
        xm = rc['per_start_x'][int(np.argmin(rc['objs']))]; mo.append(hi_err(xm))
        # two-stage MLE (moment init -> Kalman EM)
        xL = ft.kalman_em(Pfix, Ub, B, Rd, rho, sc, 1e-4, n_em=25, dev=dev, x_init=xm); ml.append(hi_err(xL))
    # fresh wall (physics reference) at one seed
    W0 = med(T, 600); wall = fresh_wall(W0, T, 0)
    row = dict(T=T, fresh_wall=wall, moment=float(np.median(mo)), mle=float(np.median(ml)),
               oracle=float(np.median(orc)), moment_all=mo, mle_all=ml, oracle_all=orc)
    results.append(row); json.dump(dict(results=results, wall_s=time.time() - t0), open('E8_results.json', 'w'), indent=2)
    log(f"T={T:5d}: fresh-wall {wall:.3f} | moment {row['moment']:.3f} | MLE {row['mle']:.3f} | oracle {row['oracle']:.3f}")

best_mle = min(r['mle'] for r in results)
log(f"\n=== E8 VERDICT: best MLE beyond-band err = {best_mle:.3f} "
    f"({'STRONG super-resolution (<=0.3)' if best_mle <= 0.3 else 'qualitative only'}) ===")
log(f"total wall {time.time()-t0:.0f}s")
