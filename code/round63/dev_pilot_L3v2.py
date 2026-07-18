"""L3 pilot v2: matched-intensity weighting on its PROPER target class.

Synthetic dev scene "L3target": low-amplitude COARSE structure (four
Gaussian blobs, sigma 3.5 px, peak +60% over background) — the class the
pre-scan CAN resolve (unlike bandpass microtexture). Compares scattered /
Lblob / Lblob x matched-weight (p=1,2) through the frozen path.
Dev-legal synthetic scene; no confirmatory data touched.
"""
import sys, os, csv, itertools
import numpy as np

os.chdir(r"D:\GI_another")
sys.path.insert(0, r"D:\GI_another\code\round63")
import campaign
import patterns as pat

SIDE, N, K = 32, 1024, 16
rng = np.random.RandomState(63)

yy, xx = np.mgrid[0:SIDE, 0:SIDE]
scene = np.ones((SIDE, SIDE))
for (cy, cx) in [(8, 9), (10, 23), (22, 7), (24, 21)]:
    scene += 0.6 * np.exp(-(((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * 3.5 ** 2)))
x = (scene / scene.sum()).ravel()
print(f"[L3v2] scene pixel contrast sd/mean = {x.std()/x.mean():.3f}", flush=True)

xi = x.reshape(SIDE, SIDE)
coarse = xi.reshape(SIDE // 4, 4, SIDE // 4, 4).mean(axis=(1, 3))
xhat = np.repeat(np.repeat(coarse, 4, axis=0), 4, axis=1).ravel()

def translate_family(offsets):
    mask = np.zeros((SIDE, SIDE))
    for (dy, dx) in offsets:
        mask[dy % SIDE, dx % SIDE] = 1.0
    B = np.empty((N, N)); r = 0
    for sy in range(SIDE):
        for sx in range(SIDE):
            B[r] = np.roll(np.roll(mask, sy, 0), sx, 1).ravel(); r += 1
    return B

LBLOB = [(0,0),(0,1),(0,2),(1,0),(1,1),(2,0),(2,1),(3,0),(3,1),(3,2),
         (4,1),(4,2),(4,3),(5,2),(5,3),(5,4)]
B = translate_family(LBLOB)

def matched(p):
    w = (xhat / xhat.mean()) ** p
    A = (N / K) * B * w[None, :]
    u_w = A @ xhat
    u_0 = ((N / K) * B) @ xhat
    return A * (u_0.mean() / u_w.mean())

CUSTOM = {
    "LblobT":    {"A": (N / K) * B},
    "Lmatch_p1": {"A": matched(1.0)},
    "Lmatch_p2": {"A": matched(2.0)},
}
for i, (kname, d) in enumerate(CUSTOM.items()):
    d.update({"exposures_per_row": 1,
              "meta": {"kind": kname, "M_signed": N, "n": N, "seed": 0,
                       "construction": "L3v2 pilot"}})
    pat.KIND_IDS[kname] = 110 + i

_orig_pat = campaign._patterns
campaign._patterns = lambda kind, M, n, seed, k=None: (
    CUSTOM[kind] if kind in CUSTOM else _orig_pat(kind, M, n, seed, k=k))
_orig_img = campaign._images
campaign._images = lambda side, spec, dev=False, imageset="conf": (
    {"L3target": x} if spec == ["L3target"]
    else _orig_img(side, spec, dev=dev, imageset=imageset))

# design-time C_u forecast per kind (the check that failed to exist last time)
for kname, d in CUSTOM.items():
    u = d["A"] @ x
    print(f"[L3v2] design C_u {kname}: {u.std()/u.mean():.3f}", flush=True)
u = (_orig_pat("sparsek", N, N, 0, k=K)["A"] @ x)
print(f"[L3v2] design C_u sparsek: {u.std()/u.mean():.3f}", flush=True)

KINDS = ["sparsek", "LblobT", "Lmatch_p1", "Lmatch_p2"]
NUS = [5, 20, 200, 1000, 2000]
OUT = r"D:\GI_another\results\round63_study2\dev_pilot_L3v2_rows.csv"
w_csv = csv.writer(open(OUT, "w", newline=""))
w_csv.writerow(["image", "patkind", "rho_bar", "nu", "seed", "PSNR_rad"])
res = {}
n_total = len(KINDS) * 2 * len(NUS) * 2
i_run = 0
for kind, rho, nu, seed in itertools.product(KINDS, [0.05, 0.6], NUS, [0, 1]):
    i_run += 1
    cell = {"side": SIDE, "nu": nu, "rho_bar": rho, "seed": seed, "M": N,
            "pattern": kind, "k": K, "arms": ["RQL"],
            "imageset": "detail32_dev", "images": ["L3target"],
            "cell_id": f"L3v2_{kind}_{rho}_{nu}_{seed}", "audit": False}
    row = campaign.run_cell(cell)[0]
    p_ = float(row["PSNR_rad"])
    res.setdefault((kind, rho, nu), []).append(p_)
    w_csv.writerow(["L3target", kind, rho, nu, seed, p_])
    print(f"[L3v2] {i_run}/{n_total} {kind} rho={rho} nu={nu} s{seed} PSNR_rad={p_:.2f}",
          flush=True)

def pava(y):
    y = list(y); n = len(y); wt = [1.0] * n; i = 0
    while i < n - 1:
        if y[i] > y[i + 1] + 1e-12:
            m = (y[i] * wt[i] + y[i + 1] * wt[i + 1]) / (wt[i] + wt[i + 1])
            y[i] = y[i + 1] = m; wt[i] = wt[i + 1] = wt[i] + wt[i + 1]
            i = max(i - 1, 0)
        else:
            i += 1
    return y

print("\n[L3v2] === endpoints ===", flush=True)
for kind in KINDS:
    curves = {}
    for rho in (0.05, 0.6):
        m = [np.mean(res[(kind, rho, nu)]) for nu in NUS]
        curves[rho] = (np.log(NUS), np.array(pava(m)))
    ls, qs = curves[0.05]; lf, qf = curves[0.6]
    R = qs[-1] - qs[0]
    dq = np.mean(res[(kind, 0.6, 2000)]) - np.mean(res[(kind, 0.05, 2000)])
    if R < 0.5:
        print(f"[L3v2] {kind:>10}: SAFE_RANGE_UNINFORMATIVE (R={R:.2f}) dQ={dq:+.2f}")
        continue
    Q90 = qs[0] + 0.9 * R
    Ts = np.exp(np.interp(Q90, qs, ls))
    if qf[-1] < Q90:
        print(f"[L3v2] {kind:>10}: FAST_RIGHT_CENSORED (R={R:.2f}) dQ={dq:+.2f}")
        continue
    Tf = np.exp(np.interp(Q90, qf, lf))
    print(f"[L3v2] {kind:>10}: R={R:.2f}  S_gate={Ts/Tf:.2f}  dQ(2000)={dq:+.2f}")
print("[L3v2] DONE", flush=True)
