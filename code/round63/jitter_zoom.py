"""Zoom sweep: dense local grids around the R14-predicted peaks +
quadratic-interpolated continuous argmax. Resolution ~1%.

Targets (nu=2000): cv=0.02 pred 10.4; cv=0.05 exact-cubic 5.69;
cv=0.1 3.50; cv=0.2 2.30. Grid: 15 points, +-22% around prediction
(ratio ~1.03). n_mc=150k. Peak = quadratic fit through best point and
neighbors in (log rho, J).
"""
import sys, importlib.util, numpy as np
spec = importlib.util.spec_from_file_location("sfi",
    r"D:\tmp\claude\E--ns-mc-gan-gi-code-fcc-phase1\3da236b2-e20e-4d0c-b724-d25feb044b56\scratchpad\jitter_scalar_fi_v2.py")
# reuse counts_three/fi_at by importing the module body with patched globals
import types
sfi = types.ModuleType("sfi")
src = open(r"D:\tmp\claude\E--ns-mc-gan-gi-code-fcc-phase1\3da236b2-e20e-4d0c-b724-d25feb044b56\scratchpad\jitter_scalar_fi_v2.py", encoding="utf-8").read()
src = src.split('print(f"[sfi] nu=')[0]          # keep defs only
exec(compile(src, "sfi_defs", "exec"), sfi.__dict__)
sfi.N_MC = 150_000
sfi.CHUNK = 6_000

TARGETS = {0.02: 10.4, 0.05: 5.69, 0.1: 3.50, 0.2: 2.30}
print("[zoom] nu=2000 n_mc=150k local grids +-22%, quad-interp peaks",
      flush=True)
for cv, pred in TARGETS.items():
    rhos = pred * np.geomspace(0.78, 1.22, 15)
    js = []
    for i, rho in enumerate(rhos):
        js.append(sfi.fi_at(float(rho), cv, seed=int(rho * 10000) + 31 * i))
    js = np.array(js)
    k = int(np.argmax(js))
    k0 = min(max(k, 1), len(rhos) - 2)
    x = np.log(rhos[k0 - 1:k0 + 2]); y = js[k0 - 1:k0 + 2]
    denom = (x[0] - x[1]) * (x[0] - x[2]) * (x[1] - x[2])
    a = (x[2] * (y[1] - y[0]) + x[1] * (y[0] - y[2]) + x[0] * (y[2] - y[1])) / denom
    b = (x[2]**2 * (y[0] - y[1]) + x[1]**2 * (y[2] - y[0]) + x[0]**2 * (y[1] - y[2])) / denom
    rho_star = float(np.exp(-b / (2 * a))) if a < 0 else float(rhos[k])
    dev = (rho_star - pred) / pred * 100
    print(f"[zoom] cv={cv:<5g} pred={pred:6.2f}  measured rho*={rho_star:6.3f} "
          f"dev={dev:+5.1f}%  (grid argmax {rhos[k]:.3f}, Jmax={js[k]:.4f})",
          flush=True)
print("[zoom] DONE", flush=True)
