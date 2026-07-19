from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(r"D:\GI_another")
SHARDS = ROOT / "results" / "round63_m1" / "shards"
DWELLS = np.array([5, 10, 20, 50, 100, 200, 500, 1000, 2000], dtype=float)
ARMS = ("LBLOB16", "RIDGE-SCAT32", "SCAT16", "SCAT32-060", "SCAT32-SAFE")
N_BOOT = 10_000
N_FAMILIES = 6


def arm_from_shard_id(shard_id: str) -> str:
    match = re.fullmatch(r"M1_(.+)_\d{2}", str(shard_id))
    if not match:
        raise ValueError(f"bad shard_id: {shard_id}")
    return match.group(1)


def family_from_image(image: str) -> str:
    match = re.fullmatch(r"m1_(.+)_([0-3])", str(image))
    if not match:
        raise ValueError(f"bad confirmatory image id: {image}")
    return match.group(1)


def pava_equal_weight(y: np.ndarray) -> np.ndarray:
    """Ordinary nondecreasing equal-weight PAVA, independently implemented."""
    y = np.asarray(y, dtype=float)
    means: list[float] = []
    weights: list[int] = []
    starts: list[int] = []
    for i, value in enumerate(y):
        means.append(float(value))
        weights.append(1)
        starts.append(i)
        while len(means) >= 2 and means[-2] > means[-1]:
            w = weights[-2] + weights[-1]
            m = (means[-2] * weights[-2] + means[-1] * weights[-1]) / w
            means[-2:] = [m]
            weights[-2:] = [w]
            starts[-2:] = [starts[-2]]
    out = np.empty_like(y, dtype=float)
    for mean, weight, start in zip(means, weights, starts):
        out[start : start + weight] = mean
    return out


def first_crossing(log_t: np.ndarray, curve: np.ndarray, target: float):
    if curve[0] >= target:
        return "LEFT", float(np.exp(log_t[0]))
    if curve[-1] < target:
        return "RIGHT", np.nan
    k = int(np.flatnonzero(curve >= target)[0])
    y0, y1 = float(curve[k - 1]), float(curve[k])
    x0, x1 = float(log_t[k - 1]), float(log_t[k])
    if y1 == y0:
        x = x1
    else:
        x = x0 + (target - y0) * (x1 - x0) / (y1 - y0)
    return "INTERIOR", float(np.exp(x))


def speed_endpoint(
    safe_curve_raw: np.ndarray,
    fast_curve_raw: np.ndarray,
    safe_time: np.ndarray,
    fast_time: np.ndarray,
):
    safe_curve = pava_equal_weight(safe_curve_raw)
    fast_curve = pava_equal_weight(fast_curve_raw)
    observed_range = float(safe_curve[-1] - safe_curve[0])
    if not np.isfinite(observed_range):
        return 0.0, "ANALYSIS_FAILURE", np.nan, safe_curve, fast_curve
    if observed_range < 0.50:
        return 1.0, "SAFE_RANGE_UNINFORMATIVE", np.nan, safe_curve, fast_curve
    target = float(safe_curve[0] + 0.90 * observed_range)
    safe_where, safe_cross = first_crossing(np.log(safe_time), safe_curve, target)
    fast_where, fast_cross = first_crossing(np.log(fast_time), fast_curve, target)
    if safe_where == "RIGHT":
        return 0.0, "ANALYSIS_FAILURE", target, safe_curve, fast_curve
    if fast_where == "RIGHT":
        return 0.0, "FAST_RIGHT_CENSORED", target, safe_curve, fast_curve
    if safe_where == "LEFT" and fast_where == "LEFT":
        return 1.0, "BOTH_LEFT_CENSORED", target, safe_curve, fast_curve
    if fast_where == "LEFT":
        return float(safe_cross / fast_time[0]), "FAST_LEFT_CENSORED", target, safe_curve, fast_curve
    if safe_where == "LEFT":
        return float(safe_time[0] / fast_cross), "SAFE_LEFT_FAST_INTERIOR", target, safe_curve, fast_curve
    return float(safe_cross / fast_cross), "NORMAL", target, safe_curve, fast_curve


def load_data() -> pd.DataFrame:
    pieces = []
    for path in sorted(SHARDS.glob("*.csv")):
        if not path.name.startswith("M1_CERT_"):
            frame = pd.read_csv(path)
            frame["source_file"] = path.name
            frame["arm_derived"] = frame["shard_id"].map(arm_from_shard_id)
            pieces.append(frame)
    data = pd.concat(pieces, ignore_index=True)
    data["family"] = data["image"].map(family_from_image)
    return data


def validate(data: pd.DataFrame) -> dict:
    grouped = data.groupby(["arm_derived", "image", "nu"], sort=True)
    counts = grouped.size()
    seed_sets = grouped["seed"].apply(lambda x: tuple(sorted(map(int, x))))
    expected_optical = data["M"].astype(float) * data["nu"].astype(float) * 50e-9
    return {
        "rows": int(len(data)),
        "arms": sorted(data["arm_derived"].unique().tolist()),
        "csv_arm_values": sorted(data["arm"].astype(str).unique().tolist()),
        "images": int(data["image"].nunique()),
        "families": data.groupby("family")["image"].nunique().sort_index().to_dict(),
        "dwells": sorted(data["nu"].astype(float).unique().tolist()),
        "all_group_counts_5": bool((counts == 5).all()),
        "group_count_minmax": [int(counts.min()), int(counts.max())],
        "all_seed_sets_0_to_4": bool(seed_sets.map(lambda x: x == (0, 1, 2, 3, 4)).all()),
        "duplicate_cell_ids": int(data["cell_id"].duplicated().sum()),
        "nonfinite_psnr_rad": int((~np.isfinite(data["PSNR_rad"].astype(float))).sum()),
        "bad_shard_to_filename": int(
            (data["source_file"].str.replace(".csv", "", regex=False) != data["shard_id"].astype(str)).sum()
        ),
        "optical_time_matches_M_nu_tau": bool(
            np.allclose(data["optical_time_s"].astype(float), expected_optical, rtol=0, atol=1e-15)
        ),
        "max_optical_time_abs_error": float(
            np.max(np.abs(data["optical_time_s"].astype(float) - expected_optical))
        ),
        "audit_status_values": sorted(data["audit_status"].dropna().astype(str).unique().tolist()),
    }


def cube(data: pd.DataFrame, arm: str):
    sub = data[data["arm_derived"] == arm].copy()
    images = sorted(sub["image"].unique())
    out = np.empty((len(images), 5, len(DWELLS)), dtype=float)
    rho = np.empty((len(images), 5, len(DWELLS)), dtype=float)
    optical_time = np.empty((len(images), 5, len(DWELLS)), dtype=float)
    for ii, image in enumerate(images):
        for ss, seed in enumerate(range(5)):
            block = sub[(sub["image"] == image) & (sub["seed"] == seed)].sort_values("nu")
            if not np.array_equal(block["nu"].to_numpy(float), DWELLS):
                raise AssertionError((arm, image, seed, block["nu"].tolist()))
            out[ii, ss] = block["PSNR_rad"].to_numpy(float)
            rho[ii, ss] = block["rho_bar"].to_numpy(float)
            optical_time[ii, ss] = block["optical_time_s"].to_numpy(float)
    return images, out, rho, optical_time


def family_index_sets(images: list[str]):
    families = sorted({family_from_image(image) for image in images})
    if len(families) != N_FAMILIES:
        raise AssertionError(families)
    by_family = []
    for family in families:
        idx = np.array([i for i, image in enumerate(images) if family_from_image(image) == family], dtype=int)
        if len(idx) != 4:
            raise AssertionError((family, idx))
        by_family.append(idx)
    return families, by_family


def rng_from_frozen_tag(seed_tag: int):
    # gi_core.config.SEED0=20260717 and rng_for(0, 63, seed_tag).
    return np.random.default_rng([20260717, 63, seed_tag])


def bootstrap_operating(images, ridge, base, seed_tag: int):
    rng = rng_from_frozen_tag(seed_tag)
    _, by_family = family_index_sets(images)
    medians = np.empty(N_BOOT, dtype=float)
    for b in range(N_BOOT):
        seed_draws = rng.integers(0, 5, size=(len(images), 5))
        delta = np.empty(len(images), dtype=float)
        for i in range(len(images)):
            seeds = seed_draws[i]
            delta[i] = float(np.mean(ridge[i, seeds, -1] - base[i, seeds, -1]))
            if not np.isfinite(delta[i]):
                delta[i] = 0.0
        sampled = np.concatenate([idx[rng.integers(0, 4, size=4)] for idx in by_family])
        medians[b] = float(np.median(delta[sampled]))
    return medians


def evaluate_speed_all(images, safe, ridge, safe_time, ridge_time, seed_draws=None):
    scores = np.empty(len(images), dtype=float)
    statuses = []
    targets = np.empty(len(images), dtype=float)
    for i in range(len(images)):
        seeds = np.arange(5) if seed_draws is None else seed_draws[i]
        safe_curve = np.mean(safe[i, seeds], axis=0)
        ridge_curve = np.mean(ridge[i, seeds], axis=0)
        score, status, target, _, _ = speed_endpoint(
            safe_curve,
            ridge_curve,
            np.mean(safe_time[i, seeds], axis=0),
            np.mean(ridge_time[i, seeds], axis=0),
        )
        scores[i] = score
        statuses.append(status)
        targets[i] = target
    return scores, statuses, targets


def bootstrap_speed(images, safe, ridge, safe_time, ridge_time, seed_tag: int):
    rng = rng_from_frozen_tag(seed_tag)
    _, by_family = family_index_sets(images)
    medians = np.empty(N_BOOT, dtype=float)
    for b in range(N_BOOT):
        seed_draws = rng.integers(0, 5, size=(len(images), 5))
        scores, _, _ = evaluate_speed_all(images, safe, ridge, safe_time, ridge_time, seed_draws)
        sampled = np.concatenate([idx[rng.integers(0, 4, size=4)] for idx in by_family])
        medians[b] = float(np.median(scores[sampled]))
    return medians


def summarize(data: pd.DataFrame):
    images, ridge, ridge_rho, ridge_t = cube(data, "RIDGE-SCAT32")
    images060, base060, _, _ = cube(data, "SCAT32-060")
    images_safe, safe, safe_rho, safe_t = cube(data, "SCAT32-SAFE")
    assert images == images060 == images_safe

    delta = np.mean(ridge[:, :, -1] - base060[:, :, -1], axis=1)

    # Independent interpretation of frozen optical time: use the CSV optical_time_s,
    # which equals M_total * nu * tau and is common to SAFE and RIDGE.
    speed, statuses, targets = evaluate_speed_all(images, safe, ridge, safe_t, ridge_t)

    # Diagnostic only: an alternative incident-exposure coordinate nu*rho. This is
    # not called optical time by the frozen spec, but is computed to catch axis mixups.
    safe_nurho = DWELLS[None, None, :] * safe_rho
    ridge_nurho = DWELLS[None, None, :] * ridge_rho
    speed_nurho, statuses_nurho, _ = evaluate_speed_all(
        images, safe, ridge, safe_nurho, ridge_nurho
    )

    operating_boot_by_seed = {
        "tag13": bootstrap_operating(images, ridge, base060, 13)
    }
    speed_boot_by_seed = {
        "tag14": bootstrap_speed(images, safe, ridge, safe_t, ridge_t, 14)
    }
    speed_nurho_boot_by_seed = {
        "tag14": bootstrap_speed(images, safe, ridge, safe_nurho, ridge_nurho, 14)
    }

    def bsum(values):
        return {
            "lb2.5": float(np.quantile(values, 0.025)),
            "median_of_bootstrap_medians": float(np.median(values)),
            "mean_of_bootstrap_medians": float(np.mean(values)),
            "min": float(np.min(values)),
            "max": float(np.max(values)),
            "n_unique": int(np.unique(values).size),
        }

    families, _ = family_index_sets(images)
    result = {
        "validation": validate(data),
        "family_order": families,
        "image_order": images,
        "RIDGE_OPERATING": {
            "per_image_delta": dict(zip(images, map(float, delta))),
            "median": float(np.median(delta)),
            "n_positive": int(np.sum(delta > 0)),
            "bootstrap": {k: bsum(v) for k, v in operating_boot_by_seed.items()},
        },
        "RIDGE_SPEED_optical_time": {
            "per_image_s_gate": dict(zip(images, map(float, speed))),
            "statuses": dict(zip(images, statuses)),
            "targets": dict(zip(images, map(float, targets))),
            "median": float(np.median(speed)),
            "n_gt_1": int(np.sum(speed > 1)),
            "bootstrap": {k: bsum(v) for k, v in speed_boot_by_seed.items()},
        },
        "RIDGE_SPEED_nu_times_rho_diagnostic": {
            "per_image_s_gate": dict(zip(images, map(float, speed_nurho))),
            "statuses": dict(zip(images, statuses_nurho)),
            "median": float(np.median(speed_nurho)),
            "n_gt_1": int(np.sum(speed_nurho > 1)),
            "bootstrap": {k: bsum(v) for k, v in speed_nurho_boot_by_seed.items()},
        },
    }
    print(json.dumps(result, indent=2, sort_keys=False, allow_nan=True))


if __name__ == "__main__":
    summarize(load_data())
