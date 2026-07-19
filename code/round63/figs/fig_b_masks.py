"""Optics-Express figure, panel (b) — representative Study-2 sparse-k masks with
equal-dose bookkeeping.

Renders the FIRST pattern (row 0 of the binary occupancy matrix B^(k)) reshaped
to the 32x32 scene, for the occupancy ladder k = 512 / 32 / 16 / 1 (left->right),
as binary images. Above/below each panel: the occupancy label, the on-pixel
intensity multiplier n/k, and a verification line (row/col sums == k, full rank).
A shared caption records the equal-dose bookkeeping A = (n/k) B, Phi = rho_bar/tau.

PROVENANCE — the displayed masks are the ACTUAL campaign matrices. Study-2's
runner (code/round63/study2_runner.py) freezes SIDE=32, so M = n = 32^2 = 1024,
PATTERN = "sparsek". campaign.run_cell then builds the illumination via
    _patterns(pattern, M, n, seed, k=k) -> make_patterns(kind, M, n, seed, k=k)
i.e. the exact call reproduced here is

    patterns.make_patterns("sparsek", 1024, 1024, seed, k=k)

with the cell seed (SEEDS5 = [0,1,2,3,4]) passed straight through — no extra
derivation — and the internal rank-rebuild nonce advanced only on deficiency.
k=512 and k=1 are the --controls masks, k=32 the --robustness mask, k=16 the
--primary mask; every stage sweeps seeds [0..4], so seed=0 (the first seed, and
the generator's own smoke seed) is a genuine campaign realization for every k.
We therefore display the seed=0, row-0 pattern for each k.

Run (cwd = repo root D:\\GI_another):
    D:\\Anacondar\\anaconda3\\python.exe code\\round63\\figs\\fig_b_masks.py
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))       # code/round63/figs
R63 = os.path.dirname(HERE)                             # code/round63
if R63 not in sys.path:
    sys.path.insert(0, R63)
from patterns import make_patterns                      # noqa: E402

# ---- frozen Study-2 geometry (study2_runner.py: SIDE, M, PATTERN, seeds) ---- #
SIDE = 32
N = SIDE * SIDE                 # n = 1024
M = N                           # runner freezes M = n (square sparsek geometry)
PATTERN = "sparsek"
SEED = 0                        # representative: first campaign seed (SEEDS5[0])
K_LADDER = [512, 32, 16, 1]     # controls(512,1) / robustness(32) / primary(16)

OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(R63)), "paper", "figs")
OUT_PNG = os.path.join(OUT_DIR, "fig_b_masks.png")
OUT_PDF = os.path.join(OUT_DIR, "fig_b_masks.pdf")


def _occupancy_label(k):
    """'k=512 (50%)' style label; percent formatted like the ruling (50/3.1/1.6/0.1)."""
    pct = 100.0 * k / N
    if pct >= 10.0:
        pstr = "%.0f%%" % pct
    else:
        pstr = "%.1f%%" % pct
    return "k=%d (%s)" % (k, pstr)


def build_and_verify():
    """Reproduce the runner's make_patterns call for each k, recover B, verify."""
    results = []
    all_ok = True
    print("[fig_b] reproducing runner call: "
          "make_patterns(\"sparsek\", M=%d, n=%d, seed=%d, k=k) for k in %s"
          % (M, N, SEED, K_LADDER), flush=True)
    for k in K_LADDER:
        d = make_patterns(PATTERN, M, N, SEED, k=k)     # EXACT runner call
        A, meta = d["A"], d["meta"]
        # recover the binary occupancy matrix B from A = (n/k) B (ruling §4)
        B = np.rint(A * (float(k) / float(N))).astype(np.int64)
        bin_ok = bool(np.all((B == 0) | (B == 1)))
        row_ok = bool(np.all(B.sum(axis=1) == k))
        col_ok = bool(np.all(B.sum(axis=0) == k))
        rank = int(np.linalg.matrix_rank(B.astype(np.float64)))
        rank_ok = (rank == N)
        mult = N // k                                   # on-pixel multiplier n/k
        nonce = meta.get("nonce")
        cond = meta.get("cond")
        occ = meta.get("occupancy", float(k) / N)
        b_sha = meta.get("B_sha256", "")
        ok = bin_ok and row_ok and col_ok and rank_ok
        all_ok = all_ok and ok
        print("[fig_b] k=%-4d occ=%8.4f%%  n/k=x%-4d  binary=%s  "
              "row_sums==k=%s  col_sums==k=%s  rank=%d/%d(%s)  nonce=%s  "
              "cond=%.6g  B_sha256=%s  %s"
              % (k, 100.0 * occ, mult, bin_ok, row_ok, col_ok, rank, N, rank_ok,
                 nonce, (cond if cond is not None else float("nan")),
                 (b_sha[:16] + ".." if b_sha else "n/a"),
                 "PASS" if ok else "FAIL"), flush=True)
        results.append(dict(k=k, B=B, mult=mult, occ=occ, nonce=nonce, cond=cond,
                            row_ok=row_ok, col_ok=col_ok, rank=rank,
                            rank_ok=rank_ok))
    print("[fig_b] VERIFICATION %s (all k: row sums==k, col sums==k, rank==%d)"
          % (("ALL PASS" if all_ok else "FAILURES ABOVE"), N), flush=True)
    return results, all_ok


def make_figure(results):
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "mathtext.fontset": "dejavusans",
        "font.size": 8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    fig, axes = plt.subplots(1, len(results), figsize=(7.0, 2.2))
    # explicit margins: reserve a bottom band for the two per-panel annotation
    # lines + the shared caption (tight_layout does not account for transAxes
    # text drawn outside the axes, which caused a caption/label collision).
    fig.subplots_adjust(left=0.03, right=0.97, top=0.84, bottom=0.34,
                        wspace=0.18)
    for ax, r in zip(axes, results):
        k = r["k"]
        row0 = r["B"][0].reshape(SIDE, SIDE)            # FIRST pattern (row 0)
        ax.imshow(row0, cmap="gray_r", interpolation="nearest", vmin=0, vmax=1)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_linewidth(0.6)
        # line 1 (above): occupancy label
        ax.set_title(_occupancy_label(k), fontsize=7, pad=4)
        # lines 2 & 3 (below): on-pixel multiplier, then verification line
        verify = r"rows=cols=%d, rank %d" % (k, r["rank"])
        ax.text(0.5, -0.06, r"on-pixel $\times n/k=\times%d$" % r["mult"],
                transform=ax.transAxes, ha="center", va="top", fontsize=7)
        ax.text(0.5, -0.185, verify,
                transform=ax.transAxes, ha="center", va="top", fontsize=7)
    # shared equal-dose bookkeeping caption under the row
    fig.text(0.5, 0.04,
             r"equal mean load $\bar{\rho}$, equal incident dose per pixel "
             r"($A=(n/k)\,B$, $\Phi=\bar{\rho}/\tau$)",
             ha="center", va="bottom", fontsize=7)
    os.makedirs(OUT_DIR, exist_ok=True)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=200)
    plt.close(fig)
    print("[fig_b] wrote %s" % OUT_PDF, flush=True)
    print("[fig_b] wrote %s (200 dpi)" % OUT_PNG, flush=True)


def main():
    results, all_ok = build_and_verify()
    make_figure(results)
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
