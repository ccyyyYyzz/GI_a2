#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
EXPLORATORY / DEV-ONLY.  Train two SCGI-style correction U-Nets on OUR frozen
forward model, local GPU (RTX 4060).  No commits.

  SCGI-EXP : trained on THEIR assumption class only (exponential-decay gain).
  SCGI-MIX : trained on BOTH families (exponential + our OU class).

Loss = weighted MSE(delta_hat, delta_true) + lambda_prior * Gaussian random-walk
prior on delta_hat in acquisition/time order (the "Gaussian constraint" = slow drift).
Early stopping on held-out val loss.  Saves state_dict + frozen signal baseline.
"""
import os, sys, json, time, argparse
import numpy as np
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import scgi_common as S

OUT = os.path.join("D:/GI_another/results/round63_next/DL_EXPLORATORY")
CKPT = os.path.join(OUT, "checkpoints")
os.makedirs(CKPT, exist_ok=True)
T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)


def prior_rw(pred):
    """Gaussian random-walk prior: mean squared first-difference of the predicted
    log-gain in TIME order (rows of the (.,1,32,64) map are row-major time order)."""
    d = pred.reshape(pred.shape[0], S.N_EXP)
    dd = d[:, 1:] - d[:, :-1]
    return (dd * dd).mean()


def train_one(name, family, baseline, args, device):
    log(f"=== TRAIN {name}  family={family} ===")
    seed = {"exp": 11, "mix": 22, "ou": 33}[family]
    tr = S.build_dataset(args.n_train, family, seed=seed, baseline=baseline)
    if baseline is None:
        baseline = tr["baseline"]
        tr = dict(tr); tr["baseline"] = baseline
    va = S.build_dataset(args.n_val, family, seed=seed + 500, baseline=baseline)
    log(f"  data: train={tuple(tr['X'].shape)} val={tuple(va['X'].shape)}")

    Xtr, Dtr, Wtr = tr["X"].to(device), tr["D"].to(device), tr["W"].to(device)
    Xva, Dva, Wva = va["X"].to(device), va["D"].to(device), va["W"].to(device)

    net = S.UNet(c=args.channels).to(device)
    opt = torch.optim.Adam(net.parameters(), lr=args.lr)
    sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, factor=0.5, patience=6)

    def wmse(pred, tgt, w):
        return ((pred - tgt) ** 2 * w).sum() / w.sum().clamp_min(1.0)

    n = Xtr.shape[0]; bs = args.batch
    best = float("inf")
    best_state = None; bad = 0; hist = []
    for ep in range(args.epochs):
        net.train(); perm = torch.randperm(n, device=device)
        tl = 0.0; nb = 0
        for i in range(0, n, bs):
            idx = perm[i:i + bs]
            pred = net(Xtr[idx])
            loss = wmse(pred, Dtr[idx], Wtr[idx]) + args.prior * prior_rw(pred)
            opt.zero_grad(); loss.backward(); opt.step()
            tl += float(loss); nb += 1
        net.eval()
        with torch.no_grad():
            vp = net(Xva)
            vloss = float(wmse(vp, Dva, Wva) + args.prior * prior_rw(vp))
            vmse = float(wmse(vp, Dva, Wva))
        sched.step(vloss)
        hist.append(dict(ep=ep, train=tl / nb, val=vloss, val_mse=vmse))
        if vloss < best - 1e-6:
            best = vloss; best_state = {k: v.detach().cpu().clone() for k, v in net.state_dict().items()}
            bad = 0
        else:
            bad += 1
        if ep % 5 == 0 or bad == 0:
            log(f"  ep {ep:3d}  train={tl/nb:.4f}  val={vloss:.4f}  val_mse={vmse:.4f}  best={best:.4f}  bad={bad}")
        if bad >= args.patience:
            log(f"  early stop at ep {ep} (best val={best:.4f})"); break

    net.load_state_dict(best_state)
    ckpt = dict(state_dict=best_state, baseline=baseline, channels=args.channels,
                family=family, name=name, best_val=best,
                meta=dict(n_train=args.n_train, n_val=args.n_val, epochs_run=len(hist),
                          lr=args.lr, batch=args.batch, prior=args.prior,
                          ou_tc_train=S.OU_TC_TRAIN, ou_cv_train=S.OU_CV_TRAIN,
                          exp_lambda_range=[0.9995, 1.0]))
    fn = os.path.join(CKPT, f"{name}.pt")
    torch.save(ckpt, fn)
    json.dump(hist, open(os.path.join(CKPT, f"{name}_history.json"), "w"), indent=1)
    log(f"  saved {fn}  best_val={best:.4f}")
    return baseline


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n_train", type=int, default=4000)
    ap.add_argument("--n_val", type=int, default=800)
    ap.add_argument("--channels", type=int, default=24)
    ap.add_argument("--epochs", type=int, default=80)
    ap.add_argument("--patience", type=int, default=16)
    ap.add_argument("--batch", type=int, default=64)
    ap.add_argument("--lr", type=float, default=2e-3)
    ap.add_argument("--prior", type=float, default=0.1)
    args = ap.parse_args()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    torch.manual_seed(0); np.random.seed(0)
    log(f"device={device}  args={vars(args)}")

    # MIX first to fix the shared frozen baseline, then EXP reuses it (fair normalisation)
    baseline = train_one("scgi_mix", "mix", None, args, device)
    train_one("scgi_exp", "exp", baseline, args, device)
    log("DONE both models")


if __name__ == "__main__":
    main()
