"""Ridge-servo dev pilot: does per-pattern gain control (load servo) beat
fixed-flux operation, and does information keep rising toward rho*(nu)?

Factorization-theorem prediction: with per-pattern attenuation alpha_i
placing every measurement's load at a common target rho_t, the optimal
rho_t is the ridge rho*(nu) (~22 at nu=2000, ~10 at nu=200), NOT the
rho=0.6 our campaigns stopped at. Servo uses only the coarse pre-scan
estimate (dev-legal). Frozen production path verbatim.

Modes per (scene, rho_t): "servo" (rows of A scaled by 1/u_hat_i so the
predicted per-pattern load is uniform) vs "fixed" (plain Lblob, load
spreads with u_i). Fixed-dwell comparison at nu in {200, 2000}.
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
    B = np.empty((N, N)); r = 0
    for sy in range(SIDE):
        for sx in range(SIDE):
            B[r] = np.roll(np.roll(mask, sy, 0), sx, 1).ravel(); r += 1
    return B

LBLOB = [(0,0),(0,1),(0,2),(1,0),(1,1),(2,0),(2,1),(3,0),(3,1),(3,2),
         (4,1),(4,2),(4,3),(5,2),(5,3),(5,4)]
B = translate_family(LBLOB)
A0 = (N / K) * B

from PIL import Image
def load_dev(name):
    x = np.asarray(Image.open(
        rf"D:\GI_another\data\r63_images_detail32_dev\32\{name}.png").convert("L"),
        float).ravel()
    return x / x.sum()

SCENES = {"detail32_dev_maze": load_dev("detail32_dev_maze"),
          "detail32_dev_glyph": load_dev("detail32_dev_glyph")}

CUSTOM = {}
for sname, x in SCENES.items():
    xi = x.reshape(SIDE, SIDE)
    coarse = xi.reshape(SIDE//4, 4, SIDE//4, 4).mean(axis=(1, 3))
    xhat = np.repeat(np.repeat(coarse, 4, axis=0), 4, axis=1).ravel()
    u_hat = A0 @ xhat                        # predicted per-pattern energy
    g = np.clip(u_hat.mean() / np.maximum(u_hat, 1e-12), 0.2, 5.0)
    As = A0 * g[:, None]  # servo with finite gain range (attenuator guard)
    kname = "servo_" + sname.split("_")[-1]
    CUSTOM[kname] = {"A": As, "exposures_per_row": 1,
                     "meta": {"kind": kname, "M_signed": N, "n": N, "seed": 0,
                              "construction": "Lblob translates, ridge servo"}}
CUSTOM["fixedLblob"] = {"A": A0, "exposures_per_row": 1,
                        "meta": {"kind": "fixedLblob", "M_signed": N, "n": N,
                                 "seed": 0, "construction": "Lblob translates"}}
for i, k_ in enumerate(CUSTOM):
    pat.KIND_IDS[k_] = 120 + i

_orig = campaign._patterns
campaign._patterns = lambda kind, M, n, seed, k=None: (
    CUSTOM[kind] if kind in CUSTOM else _orig(kind, M, n, seed, k=k))

RHOTS = [0.6, 2.0, 6.0, 22.0]
NUS = [200, 2000]
SEEDS = [0, 1]
OUT = r"D:\GI_another\results\round63_study2\dev_pilot_servo_rows.csv"
fout = open(OUT, "w", newline=""); w = csv.writer(fout)
w.writerow(["image", "mode", "rho_t", "nu", "seed", "PSNR_rad", "mean_counts"])
res = {}
jobs = [(s, mode) for s in SCENES for mode in ("servo", "fixed")]
tot = len(jobs) * len(RHOTS) * len(NUS) * len(SEEDS); i_run = 0
for (sname, mode), rho_t, nu, seed in itertools.product(jobs, RHOTS, NUS, SEEDS):
    i_run += 1
    kind = ("servo_" + sname.split("_")[-1]) if mode == "servo" else "fixedLblob"
    cell = {"side": SIDE, "nu": nu, "rho_bar": rho_t, "seed": seed, "M": N,
            "pattern": kind, "k": K, "arms": ["RQL"],
            "imageset": "detail32_dev", "images": [sname],
            "cell_id": f"servo_{sname}_{mode}_{rho_t}_{nu}_{seed}",
            "audit": False}
    row = campaign.run_cell(cell)[0]
    p = float(row["PSNR_rad"]); mc = row.get("mean_counts", "")
    res.setdefault((sname, mode, rho_t, nu), []).append(p)
    w.writerow([sname, mode, rho_t, nu, seed, p, mc]); fout.flush()
    print(f"[servo] {i_run}/{tot} {sname} {mode} rho_t={rho_t} nu={nu} s{seed} "
          f"PSNR_rad={p:.2f}", flush=True)
fout.close()

print("\n[servo] === fixed-dwell PSNR_rad vs target load ===", flush=True)
for sname in SCENES:
    for nu in NUS:
        line_s = "  ".join("%5.2f" % np.mean(res[(sname, "servo", r, nu)]) for r in RHOTS)
        line_f = "  ".join("%5.2f" % np.mean(res[(sname, "fixed", r, nu)]) for r in RHOTS)
        print(f"[servo] {sname} nu={nu}  rho_t={RHOTS}")
        print(f"[servo]    servo: {line_s}")
        print(f"[servo]    fixed: {line_f}")
print("[servo] DONE", flush=True)
