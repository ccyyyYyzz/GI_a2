"""Metric protocol — spec §5. Two pipelines, never mixed:

MAIN (deployment-legal): xhat_+ = max(xhat, 0); flux matching
xhat_* = xhat_+ * (sum x / sum xhat_+). PSNR (data_range = truth max), SSIM,
LPIPS computed on xhat_*.

DIAGNOSTIC (DIRECTION_ONLY_DIAGNOSTIC): truth-LS scalar fit on xhat_+, plus
scale-invariant angle error and Pearson correlation. Reported side-by-side;
gates read MAIN only (UT7 checks non-contamination).
"""
import numpy as np

from .utils import pearson, psnr

_LPIPS_NET = None
_LPIPS_FAILED = False


def lpips_net():
    """Spec §1 degradation path: if lpips/torch is unavailable, PSNR/SSIM still
    run and LPIPS is recorded NaN (UNAVAILABLE); gates_b marks affected combos
    pending, never passed."""
    global _LPIPS_NET, _LPIPS_FAILED
    if _LPIPS_NET is None and not _LPIPS_FAILED:
        try:
            import lpips as lpips_mod
            import torch

            torch.set_num_threads(4)
            _LPIPS_NET = lpips_mod.LPIPS(net="alex", verbose=False)
            _LPIPS_NET.eval()
        except Exception as e:
            print("LPIPS UNAVAILABLE (%r) - recording NaN, gates pending" % (e,),
                  flush=True)
            _LPIPS_FAILED = True
    return _LPIPS_NET


def lpips_batch(recs, refs, data_range):
    """recs, refs: (B, side, side) arrays on the sum=1 scale. Mapped to [-1, 1]
    via division by data_range (truth max); recon clipped into range."""
    net = lpips_net()
    if net is None:
        return np.full(len(np.atleast_3d(recs)), np.nan)
    import torch
    r = np.clip(np.asarray(recs) / data_range, 0.0, 1.0) * 2.0 - 1.0
    g = np.clip(np.asarray(refs) / data_range, 0.0, 1.0) * 2.0 - 1.0
    with torch.no_grad():
        rt = torch.from_numpy(r).float()[:, None].repeat(1, 3, 1, 1)
        gt = torch.from_numpy(g).float()[:, None].repeat(1, 3, 1, 1)
        d = net(rt, gt)
    return d.squeeze().cpu().numpy().reshape(-1)


def flux_match(xhat, x_true):
    xp = np.maximum(xhat, 0.0)
    ssum = xp.sum()
    if ssum <= 0:
        return np.zeros_like(xp)  # total collapse (e.g. GAM1 centered) -> zero map
    return xp * (x_true.sum() / ssum)


def main_metrics(xhat, x_true, side, with_lpips=True):
    """Main-protocol metrics for one reconstruction (flattened, sum=1 scale)."""
    from skimage.metrics import structural_similarity

    xs = flux_match(xhat, x_true)
    dr = float(x_true.max())
    img_r = xs.reshape(side, side)
    img_t = x_true.reshape(side, side)
    out = {
        "PSNR": float(psnr(xs, x_true, dr)),
        "SSIM": float(structural_similarity(img_t, img_r, data_range=dr)),
    }
    # AlexNet-LPIPS is undefined below ~32px (conv stack collapses to 0 size);
    # Phase A (16x16) LPIPS is recorded as NaN/UNAVAILABLE — its gates are
    # PSNR-only. Phase B (64x64) computes LPIPS normally.
    if with_lpips and side >= 32:
        out["LPIPS"] = float(lpips_batch(img_r[None], img_t[None], dr)[0])
    else:
        out["LPIPS"] = float("nan")
    return out


def diagnostic_metrics(xhat, x_true, side):
    """Truth-LS scalar fit (on xhat_+) + angle error + Pearson."""
    from skimage.metrics import structural_similarity

    xp = np.maximum(xhat, 0.0)
    dr = float(x_true.max())
    denom = float(np.dot(xp, xp))
    alpha = float(np.dot(xp, x_true)) / denom if denom > 0 else 0.0
    xls = alpha * xp
    nx = np.linalg.norm(xp)
    nt = np.linalg.norm(x_true)
    if nx > 0 and nt > 0:
        cosv = float(np.clip(np.dot(xp, x_true) / (nx * nt), -1.0, 1.0))
        ang = float(np.degrees(np.arccos(cosv)))
    else:
        ang = 90.0
    return {
        "LS_PSNR": float(psnr(xls, x_true, dr)),
        "LS_SSIM": float(structural_similarity(
            x_true.reshape(side, side), xls.reshape(side, side), data_range=dr)),
        "ANGERR_DEG": ang,
        "PEARSON": pearson(xp, x_true),
    }
