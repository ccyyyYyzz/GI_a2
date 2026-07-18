"""DEV pilot add-on: L3 matched-intensity weighting on the hardest family.

Illumination A = (n/k) * B_Lblob * diag(w), w = coarse estimate of the
scene (4x4 box-blur of truth, simulating a pre-scan), normalized so the
MEAN DETECTOR LOAD IS UNCHANGED (mean_i u_i identical to unweighted).
Dev-legal display/exploration; frozen path verbatim; no gate touched.
"""
import sys, os, csv, itertools
import numpy as np
from PIL import Image

os.chdir(r"D:\GI_another")
sys.path.insert(0, r"D:\GI_another\code\round63")
import campaign
import patterns as pat

SIDE, N, K = 32, 1024, 16

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

IMAGE = "detail32_dev_microtexture"
img = np.asarray(Image.open(
    r"D:\GI_another\data\r63_images_detail32_dev\32\detail32_dev_microtexture.png"
).convert("L"), float)
x = (img / img.sum()).ravel()

# coarse pre-scan proxy: 4x4 box blur of truth (what a cheap coarse stage sees)
xi = x.reshape(SIDE, SIDE)
coarse = xi.reshape(SIDE//4, 4, SIDE//4, 4).mean(axis=(1, 3))
xhat = np.repeat(np.repeat(coarse, 4, axis=0), 4, axis=1).ravel()
w = xhat / xhat.mean()          # relative weighting, mean 1

A = (N / K) * B * w[None, :]
# equal mean detector load: rescale so mean_i (A x-generic) matches unweighted.
# Load normalization is scene-dependent; use the coarse estimate itself
# (deployment-legal: only pre-scan info) to set the global gain.
u_w = A @ xhat
u_0 = ((N / K) * B) @ xhat
A = A * (u_0.mean() / u_w.mean())

CUSTOM = {"Lblob_matched": {"A": A, "exposures_per_row": 1,
          "meta": {"kind": "Lblob_matched", "M_signed": N, "n": N, "seed": 0,
                   "construction": "Lblob translates x coarse-matched weights"}}}
_orig = campaign._patterns
campaign._patterns = lambda kind, M, n, seed, k=None: (
    CUSTOM[kind] if kind in CUSTOM else _orig(kind, M, n, seed, k=k))
pat.KIND_IDS["Lblob_matched"] = 103

NUS = [5, 20, 200, 1000, 2000]
OUT = r"D:\GI_another\results\round63_study2\dev_pilot_L3_rows.csv"
w_csv = csv.writer(open(OUT, "w", newline=""))
w_csv.writerow(["image", "patkind", "rho_bar", "nu", "seed", "PSNR_rad"])
res = {}
for rho, nu, seed in itertools.product([0.05, 0.6], NUS, [0, 1]):
    cell = {"side": SIDE, "nu": nu, "rho_bar": rho, "seed": seed, "M": N,
            "pattern": "Lblob_matched", "k": K, "arms": ["RQL"],
            "imageset": "detail32_dev", "images": [IMAGE],
            "cell_id": f"devL3_{rho}_{nu}_{seed}", "audit": False}
    row = campaign.run_cell(cell)[0]
    p = float(row["PSNR_rad"])
    res.setdefault((rho, nu), []).append(p)
    w_csv.writerow([IMAGE, "Lblob_matched", rho, nu, seed, p])
    print(f"[L3] rho={rho} nu={nu} s{seed} PSNR_rad={p:.2f}", flush=True)

def pava(y):
    y = list(y); n = len(y); wt = [1.0]*n; i = 0
    while i < n-1:
        if y[i] > y[i+1] + 1e-12:
            m = (y[i]*wt[i]+y[i+1]*wt[i+1])/(wt[i]+wt[i+1])
            y[i] = y[i+1] = m; wt[i] = wt[i+1] = wt[i]+wt[i+1]; i = max(i-1, 0)
        else:
            i += 1
    return y

curves = {}
for rho in (0.05, 0.6):
    m = [np.mean(res[(rho, nu)]) for nu in NUS]
    curves[rho] = (np.log(NUS), np.array(pava(m)))
ls, qs = curves[0.05]; lf, qf = curves[0.6]
R = qs[-1] - qs[0]
print(f"\n[L3] {IMAGE} Lblob_matched: R={R:.2f} dB", flush=True)
if R >= 0.5:
    Q90 = qs[0] + 0.9*R
    Ts = np.exp(np.interp(Q90, qs, ls))
    if qf[-1] >= Q90:
        Tf = np.exp(np.interp(Q90, qf, lf))
        print(f"[L3] S_gate={Ts/Tf:.2f}  dQ(2000)={np.mean(res[(0.6,2000)])-np.mean(res[(0.05,2000)]):+.2f} dB")
    else:
        print(f"[L3] FAST_RIGHT_CENSORED (fast max {qf[-1]:.2f} < {Q90:.2f})")
else:
    print("[L3] SAFE_RANGE_UNINFORMATIVE")
print("[L3] DONE")
