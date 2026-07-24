# Oxygen test: does spatially-structured medium fluctuation reveal
# fixed-operator null-space content of a STATIC scene? (dev probe, clean linear)
import numpy as np

rng = np.random.default_rng(20260723)
n_side = 32; N = n_side * n_side
M = 64           # patterns per coherence epoch
T = 64           # epochs
SIG_W = 0.30     # medium fluctuation strength (30% CV)

# --- scene: smooth blobs + sharp edges (rich null content under any 64-row A)
yy, xx = np.mgrid[0:n_side, 0:n_side] / n_side
x_true = (np.exp(-((xx-.3)**2+(yy-.4)**2)/.02) + .8*np.exp(-((xx-.7)**2+(yy-.6)**2)/.01))
x_true[10:22, 15:18] += 1.0; x_true[5:8, 5:26] += .7   # bars (high-freq content)
x_true = (x_true / x_true.max()).ravel()

# --- fixed patterns (known), one set reused every epoch
P = rng.standard_normal((M, N)) / np.sqrt(N)

# smooth spatial basis for medium field: coarse grid, bilinear upsample, orthonormalized
def smooth_basis(c):
    B = []
    for i in range(c):
        for j in range(c):
            g = np.zeros((c, c)); g[i, j] = 1.0
            up = np.kron(g, np.ones((n_side//c, n_side//c)))
            # mild smoothing to avoid blocky artifacts
            up = up + .5*np.roll(up,1,0)+.5*np.roll(up,-1,0)+.5*np.roll(up,1,1)+.5*np.roll(up,-1,1)
            B.append(up.ravel())
    B = np.array(B).T
    Q, _ = np.linalg.qr(B)
    return Q  # N x c^2, orthonormal columns

def run_case(name, U):
    d_w = U.shape[1]
    Z = SIG_W * rng.standard_normal((T, d_w))          # true medium coeffs per epoch
    W = 1.0 + Z @ U.T                                   # T x N medium fields
    W = np.clip(W, 0.05, None)
    Bkt = np.einsum('mn,tn->tm', P, W * x_true)         # buckets: T x M

    # fixed-A projectors from the KNOWN pattern set
    Pi = np.linalg.pinv(P)
    x_range = lambda v: Pi @ (P @ v)                    # P_R
    null_true = x_true - x_range(x_true)
    nrm0 = np.linalg.norm(null_true)

    # --- joint alternating LS (init: medium flat -> classic fixed-A solution)
    z = np.zeros((T, d_w)); lam_z = 1e-3; lam_x = 1e-6
    for it in range(40):
        # solve x given z: stacked rows P_i * w_t
        rows = (W_est := 1.0 + z @ U.T, )[0]
        Aeff = (P[None, :, :] * rows[:, None, :]).reshape(T*M, N)
        b = Bkt.reshape(T*M)
        x_hat = np.linalg.solve(Aeff.T @ Aeff + lam_x*np.eye(N), Aeff.T @ b)
        # solve z_t given x: per-epoch small LS
        G = P @ (U * x_hat[:, None])                    # M x d_w
        GtG = G.T @ G + lam_z*np.eye(d_w)
        for t in range(T):
            r = Bkt[t] - P @ x_hat
            z[t] = np.linalg.solve(GtG, G.T @ r)
    # gauge: overall scale (E[w]=1 enforced by construction of U columns ~ zero-mean? not exactly)
    s = np.dot(x_hat, x_true) / np.dot(x_hat, x_hat)
    x_hat *= s
    null_err = np.linalg.norm((x_hat - x_range(x_hat)) - null_true) / nrm0
    tot_err = np.linalg.norm(x_hat - x_true) / np.linalg.norm(x_true)
    # baseline: classic fixed-A ridge (epoch-averaged buckets), zero null content by construction
    base_null = 1.0
    print(f"{name:28s} d_w={d_w:4d}  M={M}  null-err {null_err:.3f} (baseline {base_null:.1f})  total-err {tot_err:.3f}")
    return null_err

print("== can medium fluctuation reveal fixed-A null-space content? ==")
run_case("scalar gain (our old world)", np.ones((N,1))/np.sqrt(N))
run_case("spatial 4x4  (d_w=16 <M)",  smooth_basis(4))
run_case("spatial 8x8  (d_w=64 =M)",  smooth_basis(8))
run_case("spatial 16x16(d_w=256>M)",  smooth_basis(16))
