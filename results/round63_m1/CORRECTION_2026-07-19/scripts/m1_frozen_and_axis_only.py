"""Read-only independent reproduction of frozen M1 scoring and axis-only sensitivity.

This intentionally reimplements the currently frozen analyzer semantics, including
its historical PAVA and family-resampling choices, so their output can be compared
with the spec-correct implementation in m1_independent_recompute.py.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "independent", HERE / "m1_independent_recompute.py"
)
ind = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(ind)


def historical_pava(y):
    y = list(map(float, y))
    w = [1.0] * len(y)
    i = 0
    while i < len(y) - 1:
        if y[i] > y[i + 1] + 1e-12:
            mean = (y[i] * w[i] + y[i + 1] * w[i + 1]) / (w[i] + w[i + 1])
            y[i] = y[i + 1] = mean
            w[i] = w[i + 1] = w[i] + w[i + 1]
            i = max(i - 1, 0)
        else:
            i += 1
    return np.asarray(y)


def historical_crossing(nus, rhos, q, target, axis):
    qfit = historical_pava(q)
    x = np.asarray(nus) if axis == "nu" else np.asarray(nus) * np.asarray(rhos)
    if qfit[-1] < target:
        return None
    return float(np.exp(np.interp(target, qfit, np.log(x))))


def current_scalar_bootstrap(per_image, seed_tag):
    fam_map = {}
    for image, value in per_image.items():
        fam_map.setdefault(ind.family_from_image(image), []).append(value)
    families = sorted(fam_map)
    rng = ind.rng_from_frozen_tag(seed_tag)
    medians = np.empty(10_000)
    for b in range(10_000):
        selected_families = rng.integers(0, len(families), size=len(families))
        values = []
        for fi in selected_families:
            arr = fam_map[families[fi]]
            selected_images = rng.integers(0, len(arr), size=len(arr))
            values.extend(arr[k] for k in selected_images)
        medians[b] = np.median(values)
    return float(np.percentile(medians, 2.5))


def main():
    data = ind.load_data()
    images, ridge, ridge_rho, _ = ind.cube(data, "RIDGE-SCAT32")
    _, base060, _, _ = ind.cube(data, "SCAT32-060")
    _, safe, safe_rho, _ = ind.cube(data, "SCAT32-SAFE")

    delta = {
        image: float(np.mean(ridge[i, :, -1]) - np.mean(base060[i, :, -1]))
        for i, image in enumerate(images)
    }

    def speed(axis):
        out = {}
        for i, image in enumerate(images):
            qs = np.mean(safe[i], axis=0)
            qf = np.mean(ridge[i], axis=0)
            rs = np.mean(safe_rho[i], axis=0)
            rf = np.mean(ridge_rho[i], axis=0)
            target = float(qs[0] + 0.9 * (qs[-1] - qs[0]))
            ts = historical_crossing(ind.DWELLS, rs, qs, target, axis)
            tf = historical_crossing(ind.DWELLS, rf, qf, target, axis)
            out[image] = float(ts / tf) if ts and tf else 0.0
        return out

    frozen_speed = speed("nu_rho")
    axis_only_speed = speed("nu")

    def summary(values, seed_tag, pos_bar):
        vals = np.asarray(list(values.values()))
        return {
            "median": float(np.median(vals)),
            "lb2.5": current_scalar_bootstrap(values, seed_tag),
            "n_positive_bar": int(np.sum(vals > pos_bar)),
        }

    result = {
        "frozen_reproduction": {
            "RIDGE_OPERATING": summary(delta, 13, 0.0),
            "RIDGE_SPEED_nu_rho": summary(frozen_speed, 14, 1.0),
        },
        "axis_only_sensitivity_keep_all_other_frozen_code": {
            "RIDGE_SPEED_nu": summary(axis_only_speed, 14, 1.0)
        },
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
