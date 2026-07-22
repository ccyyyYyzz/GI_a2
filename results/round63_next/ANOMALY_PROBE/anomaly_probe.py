"""ADVERSARIAL ANOMALY-FIDELITY PROBE (round63_next).

Decisive feasibility test for the "certified information-guided imaging" claim:
    Unconstrained DL priors overwrite/erase atypical content, while
    information-guided reconstruction (data subspace locked, prior confined to
    the null space) preserves resolvable anomalies at a modest background cost;
    and the range/null split flags sub-resolution anomalies as unwitnessed.

Chassis (all loaded READ-ONLY from E:):
  * Real lane0 rate05 structured GI operator A (205x4096), reproduced in pure
    numpy and SHA-256 verified against the frozen receipt (62305c8e...).
  * Real VQAE prior (E:/.../lane0/content/gan_operator_assets/priors/vqae.pt),
    the operator's own vector-quantized autoencoder manifold-projection prior.
  * Test images = `truth` from the frozen rate05 TEST cache
    (E:/.../lane0/caches/rate05_test_cache.pt), STL-10 test split, disjoint
    from the prior's train+unlabeled training pool.

Three arms share the SAME data estimate x_range = A^+ y and the SAME prior
output p = VQAE(clip(x_range)); they differ only in how range/null are filled:
  A1 DATA-CONSISTENT BASELINE : range = data,  null = 0 (min-norm), box-fiber projected. No learned prior.
  A2 UNCONSTRAINED DL         : range = prior, null = prior      (VQAE rewrites everything, no data lock).
  A3 INFORMATION-GUIDED       : range = data (locked), null = prior null component, box-fiber projected.
So A2 and A3 share the null component and differ ONLY in the measured (range)
directions: A2 lets the prior overwrite them, A3 locks them to the data. This
isolates the mechanism exactly.

Writes ONLY under D:/GI_another/results/round63_next/ANOMALY_PROBE/.
"""
from __future__ import annotations
import math, hashlib, json, time
from pathlib import Path
import numpy as np
from scipy.linalg import hadamard
import torch
import torch.nn as nn

OUT = Path(r"D:/GI_another/results/round63_next/ANOMALY_PROBE")
OUT.mkdir(parents=True, exist_ok=True)
VQAE_CKPT = r"E:/GAN_FCC_WORK/artifacts/gan_gi_journal_round59_recovery/lane0/content/gan_operator_assets/priors/vqae.pt"
CACHE = r"E:/GAN_FCC_WORK/artifacts/gan_gi_journal_round59_recovery/lane0/caches/rate05_test_cache.pt"
EXPECTED_SHA = "62305c8e8309f772caa2b1b3601e1c45ded7224830adc667b7f6007bd97e79d6"

# ----------------------------------------------------------------------------
# Operator (verbatim row builders from src/dc_balanced.py via anomaly_split_test.py)
# ----------------------------------------------------------------------------
def normalize_row(row, *, zero_mean):
    r = np.asarray(row, dtype=np.float64).reshape(-1)
    if zero_mean:
        r = r - float(np.mean(r))
    return (r / float(np.linalg.norm(r))).astype(np.float32)
def dc_row(dim):
    return np.full((dim,), 1.0 / math.sqrt(dim), dtype=np.float32)
def random_zero_mean_rows(n, dim, seed):
    rng = np.random.default_rng(int(seed)); rows = []
    for _ in range(n):
        raw = rng.choice(np.array([-1.0, 1.0]), size=dim)
        rows.append(normalize_row(raw, zero_mean=True))
    return np.stack(rows).astype(np.float32)
def dct2_basis_row(size, u, v):
    yy = np.arange(size, dtype=np.float64); xx = np.arange(size, dtype=np.float64)
    au = math.sqrt(1.0/size) if u == 0 else math.sqrt(2.0/size)
    av = math.sqrt(1.0/size) if v == 0 else math.sqrt(2.0/size)
    by = au*np.cos(np.pi*(2*yy+1)*u/(2*size)); bx = av*np.cos(np.pi*(2*xx+1)*v/(2*size))
    return np.outer(by, bx).reshape(-1)
def dct_lowfreq_non_dc_rows(n, S):
    coords = [(u, v) for u in range(S) for v in range(S) if not (u == 0 and v == 0)]
    coords.sort(key=lambda uv: (uv[0]**2+uv[1]**2, uv[0]+uv[1], uv[0], uv[1]))
    return np.stack([normalize_row(dct2_basis_row(S, u, v), zero_mean=True) for u, v in coords[:n]]).astype(np.float32)
def hadamard_lowseq_rows(n, dim):
    h = hadamard(dim, dtype=np.int8); changes = np.sum(h[:, 1:] != h[:, :-1], axis=1)
    order = np.argsort(changes, kind="stable"); rows = []
    for idx in order:
        row = np.asarray(h[int(idx)], dtype=np.float64)
        if abs(float(np.mean(row))) > 1e-12:
            continue
        rows.append(normalize_row(row, zero_mean=True))
        if len(rows) >= n:
            break
    return np.stack(rows).astype(np.float32)
def build_operator(S, dct, had, rnd, seed):
    dim = S*S
    return np.concatenate([dc_row(dim)[None, :], dct_lowfreq_non_dc_rows(dct, S),
        hadamard_lowseq_rows(had, dim), random_zero_mean_rows(rnd, dim, seed)], axis=0).astype(np.float32)

# ----------------------------------------------------------------------------
# VQAE prior (classes copied verbatim from measurement_conditioned_vqgan.py)
# ----------------------------------------------------------------------------
class VectorQuantizer(nn.Module):
    def __init__(self, codebook_size, embed_dim, beta=0.25):
        super().__init__(); self.codebook_size = int(codebook_size); self.embed_dim = int(embed_dim); self.beta = float(beta)
        self.embedding = nn.Embedding(self.codebook_size, self.embed_dim)
        self.embedding.weight.data.uniform_(-1.0/self.codebook_size, 1.0/self.codebook_size)
    def forward(self, z_e):
        z = z_e.permute(0, 2, 3, 1).contiguous(); flat = z.reshape(-1, self.embed_dim); emb = self.embedding.weight
        dist = flat.pow(2).sum(1, keepdim=True) + emb.pow(2).sum(1)[None, :] - 2.0*flat @ emb.t()
        indices = torch.argmin(dist, dim=1); z_q = emb[indices].view_as(z).permute(0, 3, 1, 2).contiguous()
        return z_e + (z_q - z_e).detach(), indices.view(z_e.shape[0], z_e.shape[2], z_e.shape[3]), torch.zeros(()), {}
    def lookup_indices(self, indices):
        return self.embedding(indices.long()).permute(0, 3, 1, 2).contiguous()
class VQAutoencoder(nn.Module):
    def __init__(self, *, codebook_size=128, z_dim=64, base=48, beta=0.25):
        super().__init__(); c = int(base); z = int(z_dim)
        self.encoder = nn.Sequential(
            nn.Conv2d(1, c, 3, padding=1), nn.GroupNorm(min(8, c), c), nn.SiLU(inplace=True),
            nn.Conv2d(c, c*2, 4, stride=2, padding=1), nn.GroupNorm(min(8, c*2), c*2), nn.SiLU(inplace=True),
            nn.Conv2d(c*2, c*4, 4, stride=2, padding=1), nn.GroupNorm(min(8, c*4), c*4), nn.SiLU(inplace=True),
            nn.Conv2d(c*4, z, 4, stride=2, padding=1), nn.GroupNorm(min(8, z), z), nn.SiLU(inplace=True),
            nn.Conv2d(z, z, 3, padding=1))
        self.quantizer = VectorQuantizer(codebook_size, z, beta=beta)
        self.decoder = nn.Sequential(
            nn.Conv2d(z, c*4, 3, padding=1), nn.GroupNorm(min(8, c*4), c*4), nn.SiLU(inplace=True), nn.Upsample(scale_factor=2, mode="nearest"),
            nn.Conv2d(c*4, c*2, 3, padding=1), nn.GroupNorm(min(8, c*2), c*2), nn.SiLU(inplace=True), nn.Upsample(scale_factor=2, mode="nearest"),
            nn.Conv2d(c*2, c, 3, padding=1), nn.GroupNorm(min(8, c), c), nn.SiLU(inplace=True), nn.Upsample(scale_factor=2, mode="nearest"),
            nn.Conv2d(c, c, 3, padding=1), nn.SiLU(inplace=True), nn.Conv2d(c, 1, 3, padding=1), nn.Sigmoid())
    def encode(self, x):
        return self.encoder(x)
    def decode_embeddings(self, z_q):
        return self.decoder(z_q)
    def forward(self, x):
        z_e = self.encode(x); z_q, idx, _, _ = self.quantizer(z_e); return self.decode_embeddings(z_q), idx

VQGAN_CKPT = r"E:/GAN_FCC_WORK/artifacts/gan_gi_journal_round59_recovery/lane0/content/gan_operator_assets/priors/vqgan.pt"
def load_prior(ckpt_path):
    ck = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    m = VQAutoencoder(codebook_size=256, z_dim=64, base=48, beta=0.25)
    m.load_state_dict(ck["model"], strict=True)
    m.eval()
    return m

@torch.no_grad()
def vqae_apply(model, imgs_flat):
    """imgs_flat: (B,4096) in [0,1] -> VQAE reconstruction (B,4096) in [0,1]."""
    x = torch.from_numpy(np.asarray(imgs_flat, dtype=np.float32)).reshape(-1, 1, 64, 64)
    rec, _ = model(x)
    return rec.reshape(-1, 4096).numpy().astype(np.float64)

# ----------------------------------------------------------------------------
# Anomaly templates
# ----------------------------------------------------------------------------
S = 64
DELTA_RES = 0.5   # additive contrast for resolvable anomalies
DELTA_SUB = 0.7   # additive contrast for sub-resolution anomalies (need high contrast to have any chance)

def disk_mask(cy, cx, radius):
    yy, xx = np.mgrid[0:S, 0:S]
    return ((yy - cy)**2 + (xx - cx)**2) <= radius**2

def make_delta(kind, cy, cx):
    """Return an additive field d (S,S). Bright = +, dark = -. Realized template
    (post-clip) is computed against the specific base image later."""
    d = np.zeros((S, S), dtype=np.float64)
    if kind == "blob5_bright":      # ~5px diameter bright disk
        d[disk_mask(cy, cx, 2.5)] = DELTA_RES
    elif kind == "blob6_dark":      # ~6px diameter dark disk
        d[disk_mask(cy, cx, 3.0)] = -DELTA_RES
    elif kind == "bar_3x8_bright":  # 3x8 bright bar
        d[cy:cy+3, cx:cx+8] = DELTA_RES
    elif kind == "ring_r5_bright":  # 5px-diameter ring (outer r2.5 minus inner r1.2)
        ring = disk_mask(cy, cx, 2.5) & (~disk_mask(cy, cx, 1.2))
        d[ring] = DELTA_RES
    elif kind == "point_1px":       # single pixel
        d[cy, cx] = DELTA_SUB
    elif kind == "point_2px":       # 2x2 point
        d[cy:cy+2, cx:cx+2] = DELTA_SUB
    else:
        raise ValueError(kind)
    return d

RES_TYPES = ["blob5_bright", "blob6_dark", "bar_3x8_bright", "ring_r5_bright"]
SUB_TYPES = ["point_1px", "point_2px"]
RES_CLASS = {t: "resolvable" for t in RES_TYPES}
SUB_CLASS = {t: "sub_resolution" for t in SUB_TYPES}
ALL_CLASS = {**RES_CLASS, **SUB_CLASS}

def bbox_of(delta, margin=2):
    ys, xs = np.nonzero(delta)
    y0 = max(0, ys.min()-margin); y1 = min(S, ys.max()+1+margin)
    x0 = max(0, xs.min()-margin); x1 = min(S, xs.max()+1+margin)
    return y0, y1, x0, x1

# ----------------------------------------------------------------------------
# Geometry / projections
# ----------------------------------------------------------------------------
class Geo:
    def __init__(self, A):
        self.A = A.astype(np.float64)
        U, s, Vh = np.linalg.svd(self.A, full_matrices=False)
        cut = 1e-12*s.max(); keep = s > cut
        self.rank = int(keep.sum())
        self.Q = Vh[keep]                 # (r, n) orthonormal row basis = range subspace
        self.Uk = U[:, keep]; self.sk = s[keep]
    def x_range(self, y):                 # A^+ y : min-norm range-consistent solution
        return self.Q.T @ ((self.Uk.T @ y) / self.sk)
    def P_R(self, x):
        return self.Q.T @ (self.Q @ x)
    def P_N(self, x):
        return x - self.P_R(x)
    def range_frac(self, a):
        return float((self.Q @ a) @ (self.Q @ a)) / float(a @ a)

def box_fiber_project(geo, x_init, xr, iters):
    """Dykstra alternating projection onto {x : P_R x = xr} (affine fiber,
    range locked to xr) and [0,1]^n. Returns box-feasible, range-locked x."""
    x = x_init.copy()
    p = np.zeros_like(x); q = np.zeros_like(x)
    xr_range = geo.P_R(xr)  # = xr since xr in range, but be safe
    for _ in range(iters):
        # project onto affine fiber: replace range component with xr's range
        y_ = x + p
        a = y_ - geo.P_R(y_) + xr_range     # null(y_) + xr_range
        p = y_ - a
        # project onto box
        y2 = a + q
        b = np.clip(y2, 0.0, 1.0)
        q = y2 - b
        x = b
    # final exact range lock (guarantee Ax=y): keep box-feasible null, reset range
    x = geo.P_N(x) + xr_range
    return x

def record_rel_err(geo, x, y):
    return float(np.linalg.norm(geo.A @ x - y) / (np.linalg.norm(y) + 1e-12))

# ----------------------------------------------------------------------------
# Differential-bucket Poisson measurement (honest counting channel)
# ----------------------------------------------------------------------------
def measure(A, Aplus, Aminus, x_true, N0, rng):
    lam_p = N0 * (Aplus @ x_true)
    lam_m = N0 * (Aminus @ x_true)
    Bp = rng.poisson(np.maximum(lam_p, 0.0))
    Bm = rng.poisson(np.maximum(lam_m, 0.0))
    return (Bp - Bm) / N0

# ----------------------------------------------------------------------------
# Metrics
# ----------------------------------------------------------------------------
def psnr_region(x, ref, mask):
    m = mask.reshape(-1).astype(bool)
    mse = float(np.mean((x.reshape(-1)[m] - ref.reshape(-1)[m])**2))
    return 10*math.log10(1.0/max(mse, 1e-12)), mse

def matched_recovery(a, x_planted, x_clean):
    """Recovered fraction of the anomaly in the matched direction:
    rho = <a, xhat_planted - xhat_clean> / <a,a>.  1=full recovery, 0=none."""
    d = x_planted - x_clean
    return float((a @ d) / (a @ a + 1e-12))

# ----------------------------------------------------------------------------
# Matched-filter detection helpers
# ----------------------------------------------------------------------------
def gaussian_blur(x2d, sigma=2.0):
    r = int(3*sigma); k = np.exp(-(np.arange(-r, r+1)**2)/(2*sigma**2)); k /= k.sum()
    tmp = np.apply_along_axis(lambda m: np.convolve(m, k, mode="same"), 0, x2d)
    return np.apply_along_axis(lambda m: np.convolve(m, k, mode="same"), 1, tmp)

def unit_patch(kind, S=64):
    """Centered unit-norm patch for a resolvable shape (used for invention scan)."""
    d = make_delta(kind, 32, 32)
    ys, xs = np.nonzero(d)
    patch = d[ys.min():ys.max()+1, xs.min():xs.max()+1]
    return patch / (np.linalg.norm(patch)+1e-12)

def main():
    import copy as _copy
    t0 = time.time()
    print("Building operator ...")
    A = build_operator(64, 128, 56, 20, 772101)
    sha = hashlib.sha256(np.ascontiguousarray(A).tobytes()).hexdigest()
    assert sha == EXPECTED_SHA, f"operator SHA mismatch {sha}"
    geo = Geo(A)
    A64 = A.astype(np.float64)
    Aplus = np.maximum(A64, 0.0); Aminus = np.maximum(-A64, 0.0)
    print(f"  operator OK sha={sha[:12]} rank={geo.rank}")
    print("Loading priors (VQAE + VQGAN) ...")
    vqae = load_prior(VQAE_CKPT)
    vqgan = load_prior(VQGAN_CKPT)
    print("Loading test cache ...")
    cache = torch.load(CACHE, map_location="cpu", weights_only=False)
    truth = cache["truth"].numpy().reshape(-1, 4096).astype(np.float64)
    src_index = cache["source_index"].numpy()
    N_all = truth.shape[0]

    # ---- image selection (deterministic) ----
    rng_sel = np.random.default_rng(20260722)
    idx_pool = np.sort(rng_sel.choice(N_all, size=20, replace=False))
    freeze_idx = idx_pool[:2].tolist()
    control_idx = idx_pool[2:8].tolist()
    test_idx = idx_pool[8:20].tolist()
    print(f"  freeze={freeze_idx}\n  control={control_idx}\n  test={test_idx}")

    # ---- anomaly plan (frozen) ----
    rng_pos = np.random.default_rng(90210)
    def plan_for(i):
        res_type = RES_TYPES[i % len(RES_TYPES)]; sub_type = SUB_TYPES[i % len(SUB_TYPES)]
        cy_r = int(rng_pos.integers(12, 30)); cx_r = int(rng_pos.integers(12, 52))
        cy_s = int(rng_pos.integers(38, 52)); cx_s = int(rng_pos.integers(12, 52))
        return [(res_type, cy_r, cx_r), (sub_type, cy_s, cx_s)]
    plans = {img: plan_for(k) for k, img in enumerate(test_idx)}

    SEEDS = [101, 202, 303]
    REC_THRESH_ABS = 0.5      # pre-registered absolute recovery bar
    WITNESS_FLAG = 0.5        # pre-registered range-fraction flag threshold
    DET_K = 5.0               # matched-filter detection at 5 sigma (frozen)

    def N0_for_peak(peak_target):
        peak = 0.0
        for img in freeze_idx:
            x = truth[img]; peak = max(peak, float((Aplus@x).max()), float((Aminus@x).max()))
        return peak_target/peak

    def pick_proj_iters(N0):
        for it in [50, 100, 200, 400]:
            ok = True
            for img in freeze_idx:
                yc = A64@truth[img]; xr = geo.x_range(yc)
                xhat = box_fiber_project(geo, np.clip(xr, 0, 1), xr, it)
                if record_rel_err(geo, xhat, yc) > 1e-4:
                    ok = False; break
            if ok:
                return it
        return 400

    def reconstruct(y, xr, proj_iters, prior):
        base = np.clip(xr, 0.0, 1.0)
        p = vqae_apply(prior, base[None, :])[0]
        a1 = box_fiber_project(geo, base, xr, proj_iters)
        a2 = p.copy()
        a3 = box_fiber_project(geo, xr + geo.P_N(p - base), xr, proj_iters)
        return {"A1": a1, "A2": a2, "A3": a3, "prior": p, "base": base}

    def measurement_snr(N0):
        vals = []
        for img in freeze_idx:
            yc = A64@truth[img]; yn = measure(A, Aplus, Aminus, truth[img], N0, np.random.default_rng(5))
            vals.append(float(np.linalg.norm(yc)/(np.linalg.norm(yn-yc)+1e-12)))
        return float(np.mean(vals))

    def run_regime(peak_target, tag, prior):
        N0 = N0_for_peak(peak_target)
        proj_iters = pick_proj_iters(N0)
        snr = measurement_snr(N0)
        print(f"\n[{tag}] N0={N0:.1f} peak~{peak_target:.0e} proj_iters={proj_iters} measSNR~{snr:.1f}")
        records = []; recon_store = {}
        for img in test_idx:
            clean = truth[img]
            deltas = []
            for (kind, cy, cx) in plans[img]:
                d = make_delta(kind, cy, cx)
                a = np.clip(clean + d.reshape(-1), 0, 1) - clean
                deltas.append((kind, cy, cx, d, a))
            x_true = clean.copy()
            for (_, _, _, d, _) in deltas:
                x_true = np.clip(x_true + d.reshape(-1), 0, 1)
            wf = {kind: geo.range_frac(a) for (kind, _, _, _, a) in deltas}
            for seed in SEEDS:
                y_pl = measure(A, Aplus, Aminus, x_true, N0, np.random.default_rng(seed*1000+img))
                y_cl = measure(A, Aplus, Aminus, clean, N0, np.random.default_rng(seed*1000+img))
                rec_pl = reconstruct(y_pl, geo.x_range(y_pl), proj_iters, prior)
                rec_cl = reconstruct(y_cl, geo.x_range(y_cl), proj_iters, prior)
                recon_store[(img, seed)] = {"x_true": x_true, "clean": clean, "deltas": deltas, "rec": rec_pl}
                for (kind, cy, cx, d, a) in deltas:
                    y0, y1, x0, x1 = bbox_of(d, margin=2)
                    bmask = np.zeros((S, S), dtype=bool); bmask[y0:y1, x0:x1] = True
                    cls = ALL_CLASS[kind]
                    anrm = float(np.sqrt(a@a))
                    for arm in ["A1", "A2", "A3"]:
                        psnr, mse = psnr_region(rec_pl[arm], x_true, bmask)
                        rho = matched_recovery(a, rec_pl[arm], rec_cl[arm])
                        R = float((a/(anrm+1e-12)) @ (rec_pl[arm]-rec_cl[arm]))  # matched-filter response
                        records.append({
                            "image": int(img), "src_index": int(src_index[img]), "seed": int(seed),
                            "anomaly": kind, "anomaly_class": cls, "center": [cy, cx],
                            "witnessed_fraction": wf[kind], "a_norm": anrm, "arm": arm,
                            "bbox_psnr_db": psnr, "bbox_mse": mse,
                            "recovered_fraction": rho, "mf_response": R,
                        })
        # ---- noise floor per arm & class from control images (independent-seed diff) ----
        sigma = {arm: {"resolvable": [], "sub_resolution": []} for arm in ["A1", "A2", "A3"]}
        rng_probe = np.random.default_rng(4242)
        for img in control_idx:
            clean = truth[img]
            # two independent noise realizations of the SAME clean image
            y_a = measure(A, Aplus, Aminus, clean, N0, np.random.default_rng(9000+img))
            y_b = measure(A, Aplus, Aminus, clean, N0, np.random.default_rng(9500+img))
            rec_a = reconstruct(y_a, geo.x_range(y_a), proj_iters, prior)
            rec_b = reconstruct(y_b, geo.x_range(y_b), proj_iters, prior)
            for cls, types in [("resolvable", RES_TYPES), ("sub_resolution", SUB_TYPES)]:
                for _ in range(6):
                    kind = types[int(rng_probe.integers(len(types)))]
                    cy = int(rng_probe.integers(12, 50)); cx = int(rng_probe.integers(12, 50))
                    d = make_delta(kind, cy, cx); a = d.reshape(-1); anrm = float(np.sqrt(a@a))
                    if anrm < 1e-9:
                        continue
                    ah = a/anrm
                    for arm in ["A1", "A2", "A3"]:
                        sigma[arm][cls].append(float(ah @ (rec_a[arm]-rec_b[arm])))
        sigma_val = {arm: {cls: (float(np.std(sigma[arm][cls])) if sigma[arm][cls] else 0.0)
                           for cls in ["resolvable", "sub_resolution"]} for arm in ["A1", "A2", "A3"]}
        # detection: survived if mf_response > DET_K * sigma_arm_cls
        for r in records:
            thr = DET_K * sigma_val[r["arm"]][r["anomaly_class"]]
            r["detect_threshold"] = thr
            r["survived"] = bool(r["mf_response"] > thr)
            r["survived_abs50"] = bool(r["recovered_fraction"] >= REC_THRESH_ABS)

        # ---- background PSNR ----
        bg_records = []
        for img in test_idx:
            pd = [make_delta(k, cy, cx) for (k, cy, cx) in plans[img]]
            outmask = np.ones((S, S), dtype=bool)
            for d in pd:
                y0, y1, x0, x1 = bbox_of(d, margin=2); outmask[y0:y1, x0:x1] = False
            clean = truth[img]
            for seed in SEEDS:
                rec = recon_store[(img, seed)]["rec"]
                for arm in ["A1", "A2", "A3"]:
                    psnr, _ = psnr_region(rec[arm], clean, outmask)
                    bg_records.append({"image": int(img), "seed": int(seed), "arm": arm, "bg_psnr_db": psnr})

        # ---- invention on anomaly-free controls (excess anomaly-strength structure) ----
        # A false anomaly = a localized bright/dark blob in the reconstruction that
        # exceeds the natural blob-structure level of real images. Threshold tau_inv is
        # calibrated from the TRUTH control images' own residual blob responses (the
        # natural rate), so counts are directly comparable across arms and regimes.
        templates = [unit_patch("blob5_bright"), unit_patch("blob6_dark")]
        def blob_local_maxima(x2d, tau):
            resid = x2d - gaussian_blur(x2d, 2.0); cnt = 0; resps = []
            for tpl in templates:
                th, tw = tpl.shape
                rmap = np.full((S, S), -1e9)
                for yy in range(6, S-6-th):
                    for xx in range(6, S-6-tw):
                        v = abs(float((resid[yy:yy+th, xx:xx+tw]*tpl).sum()))
                        rmap[yy, xx] = v; resps.append(v)
                if tau is not None:
                    for yy in range(7, S-7-th):
                        for xx in range(7, S-7-tw):
                            v = rmap[yy, xx]
                            if v > tau and v >= rmap[yy-1:yy+2, xx-1:xx+2].max():
                                cnt += 1
            return cnt, resps
        # calibrate tau_inv = 99.5th pct of TRUTH control-image residual blob responses
        truth_resps = []
        for img in control_idx:
            _, rs = blob_local_maxima(truth[img].reshape(S, S), None); truth_resps += rs
        tau_inv = float(np.percentile(truth_resps, 99.5))
        truth_hits = sum(blob_local_maxima(truth[img].reshape(S, S), tau_inv)[0] for img in control_idx)
        invention = {"A1": 0, "A2": 0, "A3": 0}; inv_detail = []
        for img in control_idx:
            clean = truth[img]
            yc = measure(A, Aplus, Aminus, clean, N0, np.random.default_rng(777000+img))
            rec = reconstruct(yc, geo.x_range(yc), proj_iters, prior)
            for arm in ["A1", "A2", "A3"]:
                h, _ = blob_local_maxima(rec[arm].reshape(S, S), tau_inv)
                invention[arm] += h; inv_detail.append({"image": int(img), "arm": arm, "hits": int(h)})

        # ---- aggregate ----
        def agg(cls, arm, field):
            return np.array([r[field] for r in records if r["anomaly_class"] == cls and r["arm"] == arm])
        summ = {}
        for cls in ["resolvable", "sub_resolution"]:
            summ[cls] = {}
            for arm in ["A1", "A2", "A3"]:
                bp = agg(cls, arm, "bbox_psnr_db"); rf = agg(cls, arm, "recovered_fraction")
                sv = np.array([1.0 if r["survived"] else 0.0 for r in records if r["anomaly_class"] == cls and r["arm"] == arm])
                wf = agg(cls, arm, "witnessed_fraction")
                eff = rf/np.maximum(wf, 1e-6)
                summ[cls][arm] = {"n": int(bp.size),
                    "bbox_psnr_mean": float(bp.mean()), "bbox_psnr_median": float(np.median(bp)),
                    "recovered_fraction_mean": float(rf.mean()), "recovered_fraction_median": float(np.median(rf)),
                    "recovery_efficiency_mean": float(eff.mean()),  # rho / witnessed_fraction
                    "survival_rate": float(sv.mean())}
        def paired(cls, field, aa, ab):
            kf = lambda r: (r["image"], r["seed"], r["anomaly"])
            da = {kf(r): r[field] for r in records if r["anomaly_class"] == cls and r["arm"] == aa}
            db = {kf(r): r[field] for r in records if r["anomaly_class"] == cls and r["arm"] == ab}
            return np.array([da[k]-db[k] for k in da if k in db])
        d_psnr = paired("resolvable", "bbox_psnr_db", "A3", "A2")
        d_rho = paired("resolvable", "recovered_fraction", "A3", "A2")
        headline = {"bbox_psnr_mean_db": float(d_psnr.mean()), "bbox_psnr_median_db": float(np.median(d_psnr)),
                    "bbox_psnr_win_rate": float((d_psnr > 0).mean()), "recovered_fraction_mean": float(d_rho.mean()),
                    "n_pairs": int(d_psnr.size)}
        bg = {}
        for arm in ["A1", "A2", "A3"]:
            v = np.array([r["bg_psnr_db"] for r in bg_records if r["arm"] == arm]); bg[arm] = float(v.mean())
        bg["A2_minus_A3_mean_db"] = float(np.mean([r["bg_psnr_db"] for r in bg_records if r["arm"] == "A2"]) -
                                          np.mean([r["bg_psnr_db"] for r in bg_records if r["arm"] == "A3"]))
        bg["A3_minus_A1_mean_db"] = float(np.mean([r["bg_psnr_db"] for r in bg_records if r["arm"] == "A3"]) -
                                          np.mean([r["bg_psnr_db"] for r in bg_records if r["arm"] == "A1"]))
        wf_res = np.array([r["witnessed_fraction"] for r in records if r["anomaly_class"] == "resolvable" and r["arm"] == "A1"])
        wf_sub = np.array([r["witnessed_fraction"] for r in records if r["anomaly_class"] == "sub_resolution" and r["arm"] == "A1"])
        flagging = {"resolvable_flagged": int((wf_res < WITNESS_FLAG).sum()), "resolvable_total": int(wf_res.size),
                    "sub_flagged": int((wf_sub < WITNESS_FLAG).sum()), "sub_total": int(wf_sub.size),
                    "threshold": WITNESS_FLAG, "wf_res_mean": float(wf_res.mean()), "wf_res_min": float(wf_res.min()),
                    "wf_res_max": float(wf_res.max()), "wf_sub_mean": float(wf_sub.mean()), "wf_sub_max": float(wf_sub.max()),
                    "separation_min_res_minus_max_sub": float(wf_res.min()-wf_sub.max())}
        result = {"tag": tag, "N0": N0, "peak_target": peak_target, "proj_iters": proj_iters,
                  "measurement_snr": snr, "sigma_noise_floor": sigma_val, "det_k": DET_K,
                  "by_class_arm": summ, "headline_resolvable_A3_minus_A2": headline,
                  "background": bg, "witness_flagging": flagging,
                  "invention_control": {"total_hits": invention, "tau_inv_natural_995pct": tau_inv,
                                        "truth_reference_hits": int(truth_hits),
                                        "n_control_images": len(control_idx), "detail": inv_detail},
                  "records": records, "bg_records": bg_records}
        return result, recon_store

    regimes = {}
    recon_for_fig = None
    for peak, tag, prior in [(1.0e6, "clean", vqae), (1.0e4, "noisy", vqae), (1.0e6, "clean_vqgan", vqgan)]:
        res, store = run_regime(peak, tag, prior)
        res["prior_name"] = "VQGAN" if tag.endswith("vqgan") else "VQAE"
        regimes[tag] = res
        if tag == "clean":
            recon_for_fig = store

    meta = {"operator_sha256": sha, "operator_sha_match": True, "rank": geo.rank,
            "n_measurements": int(A.shape[0]), "n_pixels": int(A.shape[1]),
            "sampling_ratio": float(A.shape[0]/A.shape[1]),
            "prior": "VQAE (real lane0 vqae.pt: codebook256 z64 base48), decode-encode manifold projection",
            "delta_resolvable": DELTA_RES, "delta_subresolution": DELTA_SUB,
            "freeze_idx": freeze_idx, "control_idx": control_idx, "test_idx": test_idx, "seeds": SEEDS,
            "survival_detector": f"matched-filter response > {DET_K}*sigma (sigma from control-image independent-seed diffs)",
            "recovery_metric": "recovered_fraction rho = <a, xhat_planted-xhat_clean>/<a,a>; efficiency = rho/witnessed_fraction",
            "n_anomaly_instances_per_regime": len(regimes["clean"]["records"])//3,
            "runtime_sec": round(time.time()-t0, 1)}
    out = {"meta": meta, "regimes": regimes}
    (OUT/"anomaly_probe_results.json").write_text(json.dumps(out, indent=2), encoding="utf-8")

    # stash clean-regime recon for figures
    flat = {}
    for (img, seed), v in recon_for_fig.items():
        for arm in ["A1", "A2", "A3", "prior", "base"]:
            flat[f"{img}_{seed}_{arm}"] = v["rec"][arm]
        flat[f"{img}_{seed}_xtrue"] = v["x_true"]; flat[f"{img}_{seed}_clean"] = v["clean"]
    np.savez_compressed(OUT/"_recon_cache.npz", **flat)
    (OUT/"_plans.json").write_text(json.dumps({"plans": {str(i): [(k, int(cy), int(cx)) for (k, cy, cx) in plans[i]] for i in test_idx}, "seeds": SEEDS}), encoding="utf-8")

    # ---- console ----
    for tag in ["clean", "noisy", "clean_vqgan"]:
        r = regimes[tag]
        print(f"\n=================== REGIME {tag} ({r['prior_name']} prior, peak~{r['peak_target']:.0e}, SNR~{r['measurement_snr']:.0f}) ===================")
        for cls in ["resolvable", "sub_resolution"]:
            print(f"  [{cls}]  (witnessed frac mean {r['witness_flagging']['wf_res_mean' if cls=='resolvable' else 'wf_sub_mean']:.3f})")
            for arm in ["A1", "A2", "A3"]:
                s = r["by_class_arm"][cls][arm]
                print(f"    {arm}: bboxPSNR {s['bbox_psnr_mean']:5.2f}dB | rho {s['recovered_fraction_mean']:.3f} | eff(rho/wf) {s['recovery_efficiency_mean']:.2f} | survive {s['survival_rate']*100:3.0f}%")
        h = r["headline_resolvable_A3_minus_A2"]
        print(f"  HEADLINE resolvable A3-A2: bboxPSNR {h['bbox_psnr_mean_db']:+.2f}dB (median {h['bbox_psnr_median_db']:+.2f}), win {h['bbox_psnr_win_rate']*100:.0f}%, drho {h['recovered_fraction_mean']:+.3f}")
        print(f"  Background: A1 {r['background']['A1']:.2f} A2 {r['background']['A2']:.2f} A3 {r['background']['A3']:.2f} dB | concession A2-A3 {r['background']['A2_minus_A3_mean_db']:+.2f} | A3-A1 {r['background']['A3_minus_A1_mean_db']:+.2f}")
        print(f"  Invention hits (control): {r['invention_control']['total_hits']} vs truth-ref {r['invention_control']['truth_reference_hits']} (tau={r['invention_control']['tau_inv_natural_995pct']:.4f})")
    print(f"\n  Witness flagging (clean): {regimes['clean']['witness_flagging']}")
    print("\nWROTE", OUT/"anomaly_probe_results.json", f"({meta['runtime_sec']}s)")

if __name__ == "__main__":
    main()
