"""DEV-ONLY pilot: do clustered supports rescue the failing families?

Runs the FROZEN production path (campaign.run_cell verbatim, RQL arm,
frozen lambda_TV rule computing kappa_A/omega_A from the actual matrix)
on the three development-legal instances of the weak families, comparing
the campaign scattered k=16 pattern against clustered translate families.
No confirmatory image, seed, or gate is touched. Output: CSV of rows +
a Q90/S_gate descriptive endpoint per (image, pattern).
"""
import sys, os, csv, itertools
import numpy as np

os.chdir(r"D:\GI_another")
sys.path.insert(0, r"D:\GI_another\code\round63")
import campaign
import patterns as pat

SIDE, N, K = 32, 1024, 16

def translate_family(offsets):
    mask = np.zeros((SIDE, SIDE))
    for (dy, dx) in offsets:
        mask[dy % SIDE, dx % SIDE] = 1.0
    B = np.empty((N, N))
    r = 0
    for sy in range(SIDE):
        for sx in range(SIDE):
            B[r] = np.roll(np.roll(mask, sy, 0), sx, 1).ravel()
            r += 1
    return B

SHAPES = {
    "Lblob6x6": [(0,0),(0,1),(0,2),(1,0),(1,1),(2,0),(2,1),(3,0),(3,1),(3,2),
                 (4,1),(4,2),(4,3),(5,2),(5,3),(5,4)],
    "solid4x4": [(y, x) for y in range(4) for x in range(4)],
}
CUSTOM = {}
for name, offs in SHAPES.items():
    B = translate_family(offs)
    A = (N / K) * B                       # equal-load normalization
    CUSTOM[name] = {"A": A, "exposures_per_row": 1,
                    "meta": {"kind": name, "M_signed": N, "n": N, "seed": 0,
                             "construction": "all cyclic translates, dev pilot"}}

_orig_patterns = campaign._patterns
def routed(kind, M, n, seed, k=None):
    if kind in CUSTOM:
        return CUSTOM[kind]
    return _orig_patterns(kind, M, n, seed, k=k)
campaign._patterns = routed
# run_cell also derives its RNG stream from KIND_IDS[pattern]; register the
# pilot kinds with distinct ids far above the frozen ones (shared dict object,
# so campaign sees the update; process-local, nothing on disk changes).
pat.KIND_IDS["Lblob6x6"] = 101
pat.KIND_IDS["solid4x4"] = 102

IMAGES = ["detail32_dev_microtexture", "detail32_dev_contour", "detail32_dev_chirp"]
KINDS = {"detail32_dev_microtexture": ["sparsek", "Lblob6x6", "solid4x4"],
         "detail32_dev_contour": ["sparsek", "Lblob6x6"],
         "detail32_dev_chirp": ["sparsek", "Lblob6x6"]}
NUS = [5, 20, 200, 1000, 2000]
RHOS = [0.05, 0.6]
SEEDS = [0, 1]

OUT = r"D:\GI_another\results\round63_study2\dev_pilot_patch_rows.csv"
done = set()
if os.path.exists(OUT):
    for r in csv.DictReader(open(OUT, newline="")):
        done.add((r["image"], r["patkind"], r["rho_bar"], r["nu"], r["seed"]))
mode = "a" if done else "w"
fout = open(OUT, mode, newline="")
w = csv.writer(fout)
if mode == "w":
    w.writerow(["image", "patkind", "rho_bar", "nu", "seed", "PSNR_rad"])

total = sum(len(KINDS[i]) for i in IMAGES) * len(NUS) * len(RHOS) * len(SEEDS)
i_done = 0
for image in IMAGES:
    for kind in KINDS[image]:
        for rho, nu, seed in itertools.product(RHOS, NUS, SEEDS):
            key = (image, kind, str(float(rho)), str(float(nu)), str(seed))
            i_done += 1
            if key in done:
                continue
            cell = {"side": SIDE, "nu": nu, "rho_bar": rho, "seed": seed,
                    "M": N, "pattern": kind, "k": K, "arms": ["RQL"],
                    "imageset": "detail32_dev", "images": [image],
                    "cell_id": f"devpilot_{image}_{kind}_{rho}_{nu}_{seed}",
                    "audit": False}
            rows = campaign.run_cell(cell)
            hdr = campaign.CSV_COLUMNS if hasattr(campaign, "CSV_COLUMNS") else None
            row = rows[0]
            if isinstance(row, dict):
                psnr = row.get("PSNR_rad")
            else:
                psnr = dict(zip(hdr, row)).get("PSNR_rad") if hdr else row
            w.writerow([image, kind, rho, nu, seed, psnr])
            fout.flush()
            print(f"[devpilot] {i_done}/{total} {image} {kind} rho={rho} nu={nu} s{seed} PSNR_rad={psnr}",
                  flush=True)
fout.close()
print("[devpilot] SWEEP DONE")

# ---- descriptive Q90 endpoint per (image, patkind) ----
import collections
data = collections.defaultdict(dict)
for r in csv.DictReader(open(OUT, newline="")):
    ky = (r["image"], r["patkind"], float(r["rho_bar"]), float(r["nu"]))
    data[ky].setdefault("v", []).append(float(r["PSNR_rad"]))

def pava(y):
    y = list(y); n = len(y)
    w = [1.0] * n
    i = 0
    while i < n - 1:
        if y[i] > y[i + 1] + 1e-12:
            ym = (y[i] * w[i] + y[i + 1] * w[i + 1]) / (w[i] + w[i + 1])
            y[i] = y[i + 1] = ym; w[i] = w[i + 1] = w[i] + w[i + 1]
            i = max(i - 1, 0)
        else:
            i += 1
    return y

print("\n[devpilot] === Q90 endpoint (descriptive, dev-only) ===")
for image in IMAGES:
    for kind in KINDS[image]:
        curves = {}
        for rho in RHOS:
            pts = sorted((nu, np.mean(data[(image, kind, rho, float(nu))]["v"]))
                         for nu in NUS if (image, kind, rho, float(nu)) in data)
            nus = [p[0] for p in pts]
            iso = pava([p[1] for p in pts])
            curves[rho] = (np.log([n_ for n_ in nus]), np.array(iso))
        ls, qs = curves[0.05]
        lf, qf = curves[0.6]
        R = qs[-1] - qs[0]
        if R < 0.5:
            print(f"  {image:>28} {kind:>10}: SAFE_RANGE_UNINFORMATIVE (R={R:.2f} dB) S=1")
            continue
        Q90 = qs[0] + 0.9 * R
        Ts = np.exp(np.interp(Q90, qs, ls))
        if qf[-1] < Q90:
            print(f"  {image:>28} {kind:>10}: FAST_RIGHT_CENSORED (fast max {qf[-1]:.2f} < Q90 {Q90:.2f}) S=0")
            continue
        Tf = np.exp(np.interp(Q90, qf, lf))
        print(f"  {image:>28} {kind:>10}: R={R:.2f} dB  Q90={Q90:.2f}  S_gate={Ts/Tf:.2f}"
              f"  dQ(nu=2000)={np.mean(data[(image,kind,0.6,2000.0)]['v'])-np.mean(data[(image,kind,0.05,2000.0)]['v']):+.2f} dB")
