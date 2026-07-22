"""Q4 anomaly energy-split test for the GAN_FCC completion line.

Read-only reproduction of the lane0 rate05 / rate10 structured GI operators
(pure numpy, no torch), then measure how much of a small localized anomaly's
energy lands in the measured (range/row) subspace vs the null space.

Operator provenance (from config_rate05.yaml / config_rate10.yaml under
E:\\GAN_FCC_WORK\\artifacts\\gan_gi_journal_round59_recovery\\lane0):
  rate05: img_size=64, dct=128, hadamard=56, random=20, seed=772101, m=205
  rate10: img_size=64, dct=257, hadamard=112, random=40, seed=772011, m=410
Row builders reproduced verbatim from src/dc_balanced.py + build_structured_operator_rows.
Writes only under D:\\GI_another\\results\\round63_next\\GAN_LINE_AUDIT\\.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import hadamard

OUT = Path(r"D:/GI_another/results/round63_next/GAN_LINE_AUDIT")

# Expected operator SHA-256 from the frozen receipts (rate05 only recorded there).
EXPECTED_SHA = {"05": "62305c8e8309f772caa2b1b3601e1c45ded7224830adc667b7f6007bd97e79d6"}


# ---- row builders (verbatim logic from src/dc_balanced.py) ----
def normalize_row(row, *, zero_mean):
    r = np.asarray(row, dtype=np.float64).reshape(-1)
    if zero_mean:
        r = r - float(np.mean(r))
    norm = float(np.linalg.norm(r))
    return (r / norm).astype(np.float32)


def dc_row(dim):
    return np.full((int(dim),), 1.0 / math.sqrt(float(dim)), dtype=np.float32)


def random_zero_mean_rows(num_rows, dim, seed):
    rng = np.random.default_rng(int(seed))
    rows = []
    for _ in range(int(num_rows)):
        raw = rng.choice(np.array([-1.0, 1.0], dtype=np.float64), size=int(dim))
        rows.append(normalize_row(raw, zero_mean=True))
    return np.stack(rows, axis=0).astype(np.float32)


def dct2_basis_row(size, u, v):
    yy = np.arange(int(size), dtype=np.float64)
    xx = np.arange(int(size), dtype=np.float64)
    au = math.sqrt(1.0 / size) if int(u) == 0 else math.sqrt(2.0 / size)
    av = math.sqrt(1.0 / size) if int(v) == 0 else math.sqrt(2.0 / size)
    by = au * np.cos(np.pi * (2.0 * yy + 1.0) * int(u) / (2.0 * int(size)))
    bx = av * np.cos(np.pi * (2.0 * xx + 1.0) * int(v) / (2.0 * int(size)))
    return np.outer(by, bx).reshape(-1)


def dct_lowfreq_non_dc_rows(num_rows, img_size):
    size = int(img_size)
    coords = [(u, v) for u in range(size) for v in range(size) if not (u == 0 and v == 0)]
    coords.sort(key=lambda uv: (uv[0] * uv[0] + uv[1] * uv[1], uv[0] + uv[1], uv[0], uv[1]))
    rows = [normalize_row(dct2_basis_row(size, u, v), zero_mean=True) for u, v in coords[:int(num_rows)]]
    return np.stack(rows, axis=0).astype(np.float32)


def hadamard_lowsequency_non_dc_rows(num_rows, dim):
    h = hadamard(int(dim), dtype=np.int8)
    changes = np.sum(h[:, 1:] != h[:, :-1], axis=1)
    order = np.argsort(changes, kind="stable")
    rows = []
    for idx in order:
        row = np.asarray(h[int(idx)], dtype=np.float64)
        if abs(float(np.mean(row))) > 1e-12:
            continue
        rows.append(normalize_row(row, zero_mean=True))
        if len(rows) >= int(num_rows):
            break
    return np.stack(rows, axis=0).astype(np.float32)


def build_operator(img_size, dct_rows, hadamard_rows, random_rows, seed):
    dim = int(img_size) * int(img_size)
    rows = [dc_row(dim)[None, :]]
    rows.append(dct_lowfreq_non_dc_rows(dct_rows, img_size))
    rows.append(hadamard_lowsequency_non_dc_rows(hadamard_rows, dim))
    rows.append(random_zero_mean_rows(random_rows, dim, seed))
    A = np.concatenate(rows, axis=0).astype(np.float32)
    return A


def sha256_numpy(arr):
    return hashlib.sha256(np.ascontiguousarray(arr).tobytes()).hexdigest()


# ---- synthetic localized anomalies on 64x64 ----
def make_anomalies(img_size=64, seed=20260722):
    rng = np.random.default_rng(seed)
    S = int(img_size)
    anomalies = []
    specs = [
        ("gauss_sigma1.0", "bump", 1.0),
        ("gauss_sigma1.5", "bump", 1.5),
        ("gauss_sigma2.0", "bump", 2.0),
        ("gauss_sigma3.0", "bump", 3.0),
        ("square_3px", "square", 3),
        ("square_5px", "square", 5),
        ("stroke_h_5px", "stroke_h", 5),
        ("stroke_v_5px", "stroke_v", 5),
        ("point_2px", "square", 2),
        ("gauss_sigma1.2", "bump", 1.2),
    ]
    yy, xx = np.mgrid[0:S, 0:S]
    for name, kind, param in specs:
        img = np.zeros((S, S), dtype=np.float64)
        margin = 8
        cy = int(rng.integers(margin, S - margin))
        cx = int(rng.integers(margin, S - margin))
        if kind == "bump":
            sigma = float(param)
            img = np.exp(-((yy - cy) ** 2 + (xx - cx) ** 2) / (2.0 * sigma ** 2))
        elif kind == "square":
            h = int(param)
            img[cy:cy + h, cx:cx + h] = 1.0
        elif kind == "stroke_h":
            L = int(param)
            img[cy, cx:cx + L] = 1.0
        elif kind == "stroke_v":
            L = int(param)
            img[cy:cy + L, cx] = 1.0
        a = img.reshape(-1)
        anomalies.append((name, cy, cx, a))
    return anomalies


def range_null_split(A, anomalies):
    # Orthonormal row-space basis via SVD (matches GaugeGeometry construction).
    A64 = A.astype(np.float64)
    U, s, Vh = np.linalg.svd(A64, full_matrices=False)
    cutoff = 1e-12 * float(s.max())
    keep = s > cutoff
    rank = int(keep.sum())
    Q = Vh[keep]  # (rank, n), orthonormal rows spanning row(A) = range subspace
    # DC direction (row 0 normalized) to separate the trivially-measured mean.
    dc = A64[0] / np.linalg.norm(A64[0])
    results = []
    for name, cy, cx, a in anomalies:
        a = a.astype(np.float64)
        total = float(a @ a)
        coeff = Q @ a  # projections onto orthonormal row basis
        range_e = float(coeff @ coeff)
        dc_e = float((dc @ a) ** 2)
        results.append({
            "anomaly": name, "center": [cy, cx],
            "total_energy": total,
            "range_fraction": range_e / total,
            "null_fraction": 1.0 - range_e / total,
            "dc_fraction": dc_e / total,
            "range_fraction_excl_dc": (range_e - dc_e) / total,
        })
    return results, rank


def main():
    ops = {
        "05": dict(img_size=64, dct_rows=128, hadamard_rows=56, random_rows=20, seed=772101, m=205),
        "10": dict(img_size=64, dct_rows=257, hadamard_rows=112, random_rows=40, seed=772011, m=410),
    }
    anomalies = make_anomalies()
    out = {"n_pixels": 4096, "img_size": 64, "operators": {}}
    for rate, cfg in ops.items():
        A = build_operator(cfg["img_size"], cfg["dct_rows"], cfg["hadamard_rows"],
                           cfg["random_rows"], cfg["seed"])
        sha = sha256_numpy(A)
        sha_ok = EXPECTED_SHA.get(rate)
        results, rank = range_null_split(A, anomalies)
        rf = np.array([r["range_fraction"] for r in results])
        rf_nodc = np.array([r["range_fraction_excl_dc"] for r in results])
        out["operators"][f"rate{rate}"] = {
            "config": cfg,
            "operator_sha256": sha,
            "operator_sha256_expected": sha_ok,
            "operator_sha256_match": (sha == sha_ok) if sha_ok else None,
            "m_rows": int(A.shape[0]),
            "numerical_rank": rank,
            "sampling_ratio": A.shape[0] / A.shape[1],
            "per_anomaly": results,
            "range_fraction_mean": float(rf.mean()),
            "range_fraction_min": float(rf.min()),
            "range_fraction_max": float(rf.max()),
            "range_fraction_excl_dc_mean": float(rf_nodc.mean()),
        }
        print(f"\n=== rate{rate}  m={A.shape[0]} rank={rank} sha_match={out['operators'][f'rate{rate}']['operator_sha256_match']} ===")
        print(f"  mean range fraction = {rf.mean():.4f}  (excl DC = {rf_nodc.mean():.4f})   null = {1-rf.mean():.4f}")
        for r in results:
            print(f"  {r['anomaly']:16s} range={r['range_fraction']:.4f} "
                  f"(exclDC={r['range_fraction_excl_dc']:.4f}, DC={r['dc_fraction']:.4f}) null={r['null_fraction']:.4f}")
    (OUT / "anomaly_split.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nWROTE {OUT / 'anomaly_split.json'}")


if __name__ == "__main__":
    main()
