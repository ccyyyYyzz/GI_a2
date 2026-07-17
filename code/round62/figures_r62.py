"""ROUND62 mandatory figures — spec §4 output contract.

Renders the required figures for the three ROUND62 parts over the committed
result CSV/JSON (never touches truth; this is a pure rendering layer):

  Part 1 (T1-G0, results/round62_g0/figures/):
    g0_permutation_collapse.png  — permutation/noise collapse bar chart
    g0_m_scaling.png             — MAVE-vs-WHITEN-LW gain vs M (smoothing-prior direction)
    g0_sensitivity.png           — bandwidth x ridge sensitivity heatmap (seed 0)
  Part 3 (P3, results/round62_p3/figures/):
    p3_psnr_vs_L.png             — PSNR-L crossover curve family, faceted CV x rho
  Part 2 (T1-G1, results/round62_g1/figures/, only if g1_scale_map.csv exists):
    g1_headroom_vs_n.png         — headroom vs n per (illum x link) at M=1e5
    g1_headroom_vs_Mn.png        — headroom vs M/n pooled scatter

Style follows code/figures.py: matplotlib Agg backend, 200 dpi, tight_layout,
each figure wrapped in try/except in main() so a missing/partial input skips
gracefully instead of aborting the batch.

Schema reference (read from the producing scripts):
  g0_metrics.csv : arm, config, seed, M, image, PSNR
  p3_crossover.csv : illum, cv, rho, seed, image, arm, readout, PSNR, SSIM
  g1_scale_map.csv : NOT yet frozen upstream — assumed columns
                     illum(/family), link, n, M, seed, image, arm, PSNR
                     with arm in {semi_oracle | opg | rmave, whiten_lw,
                     l_isotron, whiten_or}. See fig5 for the exact fallbacks.
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))  # D:\GI_another
RES = os.path.join(ROOT, "results")
G0 = os.path.join(RES, "round62_g0")
G1 = os.path.join(RES, "round62_g1")
P3 = os.path.join(RES, "round62_p3")

IMAGES = ["cat", "deer", "dog", "horse", "airplane", "car", "ship", "truck"]
PERM_GATE_DB = 12.0    # permuted / noise arms must sit at/below this
SEP_GATE_DB = 15.0     # required unperm - permuted separation

SENS_BW = [0.5, 1.0, 2.0]
SENS_RIDGE = [0.1, 1.0, 10.0]

P3_CVS = [0.02, 0.05, 0.1]
P3_RHOS = [0.99, 0.9]
P3_LS = [2, 4, 8, 16, 32, 64]

# ---------------------------------------------------------------- Part 1: G0


def fig1_g0_permutation_collapse():
    """Bar chart: 8-image 3-seed mean PSNR of the MAVE permutation-null arms
    (unperm / permuted-pool / noise) plus WHITEN-LW, L-ISOTRON, MLE-renewal for
    context. 12 dB gate line + annotated unperm-vs-permuted separation."""
    df = pd.read_csv(os.path.join(G0, "g0_metrics.csv"))
    df["PSNR"] = df["PSNR"].astype(float)
    full = df[df.M == df.M.max()]  # exclude M-scaling prefixes (M<max)

    def seed_means(mask):
        """Per-seed 8-image means -> (mean, std over seeds)."""
        s = full[mask]
        if s.empty:
            return np.nan, 0.0
        g = s.groupby("seed")["PSNR"].mean()
        return float(g.mean()), float(g.std(ddof=0))

    bars = [
        ("MAVE\nunperm", full.arm == "mave_unperm", "#2ca02c"),
        ("MAVE\npermuted", full.arm.str.startswith("mave_perm"), "#d62728"),
        ("MAVE\nnoise", full.arm == "mave_noise", "#7f7f7f"),
        ("WHITEN-LW", full.arm == "whiten_lw", "#1f77b4"),
        ("L-ISOTRON", full.arm == "l_isotron", "#17becf"),
        ("MLE-renewal", full.arm == "mle_renewal", "#9467bd"),
    ]
    labels, means, errs, colors = [], [], [], []
    for lab, mask, col in bars:
        m, e = seed_means(mask)
        labels.append(lab)
        means.append(m)
        errs.append(e)
        colors.append(col)

    fig, ax = plt.subplots(figsize=(9, 5))
    xs = np.arange(len(labels))
    ax.bar(xs, means, yerr=errs, color=colors, capsize=4, alpha=0.9,
           edgecolor="black", linewidth=0.5)
    for x, m in zip(xs, means):
        if np.isfinite(m):
            ax.text(x, m + 0.3, "%.2f" % m, ha="center", va="bottom", fontsize=8)

    ax.axhline(PERM_GATE_DB, color="red", lw=1.0, ls="--",
               label="permutation gate (<= %g dB)" % PERM_GATE_DB)

    # annotate unperm - permuted separation
    m_unperm, m_perm = means[0], means[1]
    if np.isfinite(m_unperm) and np.isfinite(m_perm):
        sep = m_unperm - m_perm
        y_hi = m_unperm + 0.6
        ax.annotate("", xy=(0, y_hi), xytext=(1, y_hi),
                    arrowprops=dict(arrowstyle="<->", color="black", lw=1.2))
        ax.text(0.5, y_hi + 0.2,
                "separation = %.1f dB\n(gate >= %g dB)" % (sep, SEP_GATE_DB),
                ha="center", va="bottom", fontsize=9,
                color=("black" if sep >= SEP_GATE_DB else "red"))

    ax.set_xticks(xs)
    ax.set_xticklabels(labels)
    ax.set_ylabel("8-image 3-seed mean PSNR (dB)")
    ax.set_title("G0 permutation-null collapse (16x16, GAM4xDT30, 3 seeds)\n"
                 "unperm bit-identical to the audited +6.9 dB anomaly")
    ax.legend(loc="upper right", fontsize=8)
    _save(fig, os.path.join(G0, "figures", "g0_permutation_collapse.png"))


def fig2_g0_m_scaling():
    """MAVE-unperm minus WHITEN-LW (8-image mean) vs M on log-x. Gain rising as
    M shrinks is the smoothing-prior signature (flagged suspicious in spec)."""
    df = pd.read_csv(os.path.join(G0, "g0_metrics.csv"))
    df["PSNR"] = df["PSNR"].astype(float)

    Ms, gains = [], []
    for M in sorted(df.M.unique()):
        sub = df[df.M == M]
        um = sub[sub.arm == "mave_unperm"].groupby("image")["PSNR"].mean()
        lw = sub[sub.arm == "whiten_lw"].groupby("image")["PSNR"].mean()
        if um.empty or lw.empty:
            continue
        common = um.index.intersection(lw.index)
        if len(common) == 0:
            continue
        Ms.append(int(M))
        gains.append(float((um[common] - lw[common]).mean()))

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.semilogx(Ms, gains, "o-", color="#2ca02c", lw=2, markersize=7)
    for m, g in zip(Ms, gains):
        ax.annotate("%.2f dB" % g, (m, g), textcoords="offset points",
                    xytext=(0, 8), ha="center", fontsize=8)

    suspicious = len(gains) >= 2 and all(
        gains[i] > gains[i + 1] for i in range(len(gains) - 1))
    if len(Ms) >= 2:
        ymid = (min(gains) + max(gains)) / 2.0
        ax.annotate("smaller M -> larger gain\n= smoothing-prior signature%s"
                    % ("\n(SUSPICIOUS)" if suspicious else ""),
                    xy=(Ms[0], gains[0]), xytext=(Ms[-1], ymid),
                    fontsize=9, color="#8B0000", ha="center",
                    arrowprops=dict(arrowstyle="->", color="#8B0000", lw=1.2))

    ax.set_xscale("log")
    ax.set_xticks(Ms)
    ax.set_xticklabels([str(m) for m in Ms])
    ax.set_xlabel("M (frames)")
    ax.set_ylabel("MAVE-unperm - WHITEN-LW, 8-image mean (dB)")
    ax.set_title("G0 M-scaling diagnostic (soft gate, report-only)")
    ax.grid(True, which="both", alpha=0.3)
    _save(fig, os.path.join(G0, "figures", "g0_m_scaling.png"))


def fig3_g0_sensitivity():
    """3x3 (bw_scale x ridge_scale) heatmap of seed-0 8-image mean PSNR from the
    mave_sens arm; center cell (bw=1,ridge=1) taken from mave_unperm seed 0.
    Worst (lowest-PSNR) config outlined."""
    df = pd.read_csv(os.path.join(G0, "g0_metrics.csv"))
    df["PSNR"] = df["PSNR"].astype(float)
    full = df[df.M == df.M.max()]

    Z = np.full((len(SENS_BW), len(SENS_RIDGE)), np.nan)
    for i, bw in enumerate(SENS_BW):
        for j, rs in enumerate(SENS_RIDGE):
            if bw == 1.0 and rs == 1.0:
                sub = full[(full.arm == "mave_unperm") & (full.seed == 0)]
            else:
                tag = "bw%g_r%g" % (bw, rs)
                sub = full[(full.arm == "mave_sens")
                           & (full.config == tag) & (full.seed == 0)]
            if not sub.empty:
                Z[i, j] = float(sub["PSNR"].mean())

    fig, ax = plt.subplots(figsize=(6.5, 5))
    im = ax.imshow(Z, cmap="viridis", aspect="auto")
    for i in range(len(SENS_BW)):
        for j in range(len(SENS_RIDGE)):
            if np.isfinite(Z[i, j]):
                ax.text(j, i, "%.2f" % Z[i, j], ha="center", va="center",
                        color="white", fontsize=10, fontweight="bold")

    if np.isfinite(Z).any():
        wi, wj = np.unravel_index(np.nanargmin(Z), Z.shape)
        ax.add_patch(plt.Rectangle((wj - 0.5, wi - 0.5), 1, 1, fill=False,
                                   edgecolor="red", lw=3))
        ax.text(wj, wi + 0.32, "worst", ha="center", va="center",
                color="red", fontsize=8, fontweight="bold")

    ax.set_xticks(range(len(SENS_RIDGE)))
    ax.set_xticklabels(["x%g" % r for r in SENS_RIDGE])
    ax.set_yticks(range(len(SENS_BW)))
    ax.set_yticklabels(["x%g" % b for b in SENS_BW])
    ax.set_xlabel("ridge scale")
    ax.set_ylabel("bandwidth scale")
    ax.set_title("G0 MAVE bandwidth x ridge sensitivity\n"
                 "seed-0 8-image mean PSNR (center = default mave_unperm)")
    fig.colorbar(im, ax=ax, shrink=0.85, label="PSNR (dB)")
    _save(fig, os.path.join(G0, "figures", "g0_sensitivity.png"))


# ---------------------------------------------------------------- Part 3: P3

_ILLUM_COLOR = {"GAM4": "#1f77b4", "MIX-LOGN": "#ff7f0e"}
_SOLID_COLOR = {
    ("GAM4", "plain"): "#1f77b4",
    ("GAM4", "whiten"): "#17becf",
    ("MIX-LOGN", "plain"): "#ff7f0e",
    ("MIX-LOGN", "whiten"): "#d62728",
}
_SOLID_MARKER = {"plain": "o", "whiten": "s"}
_FALLBACK = ["#8c564b", "#e377c2", "#bcbd22", "#7f7f7f"]


def _color_solid(illum, readout, k):
    return _SOLID_COLOR.get((illum, readout), _FALLBACK[k % len(_FALLBACK)])


def fig4_p3_psnr_vs_L():
    """PSNR-L crossover family, faceted rows=CV x cols=rho. Solid lines =
    (illum x {plain,whiten}) over numeric L; thin dashed = ratio readout;
    dotted horizontals = GLOBAL (per series) and RAW (per illum). argmax L
    marked per solid line."""
    df = pd.read_csv(os.path.join(P3, "p3_crossover.csv"))
    df["PSNR"] = df["PSNR"].astype(float)
    df["cv"] = df["cv"].astype(float)
    df["rho"] = df["rho"].astype(float)
    illums = [f for f in ("GAM4", "MIX-LOGN") if f in set(df.illum.unique())]
    if not illums:
        illums = sorted(df.illum.unique())

    def mean_of(sub, illum, arm, readout):
        s = sub[(sub.illum == illum) & (sub.arm == arm)
                & (sub.readout == readout)]
        return float(s["PSNR"].mean()) if not s.empty else np.nan

    fig, axes = plt.subplots(len(P3_CVS), len(P3_RHOS), figsize=(12, 12),
                             sharex=True, squeeze=False)
    for i, cv in enumerate(P3_CVS):
        for j, rho in enumerate(P3_RHOS):
            ax = axes[i][j]
            sub = df[(df.cv == cv) & (df.rho == rho)]
            for ki, illum in enumerate(illums):
                # solid: plain / whiten over L
                for readout in ("plain", "whiten"):
                    ys = np.array([mean_of(sub, illum, "L%d" % L, readout)
                                   for L in P3_LS])
                    if not np.isfinite(ys).any():
                        continue
                    col = _color_solid(illum, readout, ki)
                    ax.plot(P3_LS, ys, marker=_SOLID_MARKER.get(readout, "o"),
                            color=col, lw=1.8, markersize=5)
                    a = np.nanargmax(ys)
                    ax.plot(P3_LS[a], ys[a], marker="*", color=col,
                            markersize=15, markeredgecolor="black",
                            markeredgewidth=0.6, zorder=5)
                    ax.annotate("L*=%d" % P3_LS[a], (P3_LS[a], ys[a]),
                                textcoords="offset points", xytext=(4, 6),
                                fontsize=7, color=col)
                    # GLOBAL horizontal reference for this series
                    gv = mean_of(sub, illum, "GLOBAL", readout)
                    if np.isfinite(gv):
                        ax.axhline(gv, color=col, ls=":", lw=1.0, alpha=0.55)
                # ratio: thin dashed
                yr = np.array([mean_of(sub, illum, "L%d" % L, "ratio")
                               for L in P3_LS])
                if np.isfinite(yr).any():
                    ax.plot(P3_LS, yr, ls="--", lw=1.0,
                            color=_ILLUM_COLOR.get(illum, "#7f7f7f"), alpha=0.6)
                # RAW horizontal (plain readout, failure showcase)
                rv = mean_of(sub, illum, "RAW", "plain")
                if np.isfinite(rv):
                    ax.axhline(rv, color="#555555", ls=":", lw=0.9, alpha=0.4)

            ax.set_xscale("log", base=2)
            ax.set_xticks(P3_LS)
            ax.set_xticklabels(P3_LS)
            ax.set_title("CV=%g,  rho=%g" % (cv, rho), fontsize=10)
            ax.grid(True, which="both", alpha=0.25)
            if j == 0:
                ax.set_ylabel("PSNR (dB)")
            if i == len(P3_CVS) - 1:
                ax.set_xlabel("block length L")

    # figure-level legend (proxies)
    from matplotlib.lines import Line2D
    handles = []
    for illum in illums:
        for readout in ("plain", "whiten"):
            handles.append(Line2D([0], [0],
                           color=_color_solid(illum, readout, 0),
                           marker=_SOLID_MARKER.get(readout, "o"),
                           label="%s %s" % (illum, readout)))
    for illum in illums:
        handles.append(Line2D([0], [0], ls="--", lw=1.0,
                       color=_ILLUM_COLOR.get(illum, "#7f7f7f"),
                       label="%s ratio (report-only)" % illum))
    handles.append(Line2D([0], [0], ls=":", color="black",
                   label="GLOBAL (per series)"))
    handles.append(Line2D([0], [0], ls=":", color="#555555",
                   label="RAW (plain)"))
    handles.append(Line2D([0], [0], marker="*", color="black", lw=0,
                   markersize=12, label="argmax L"))
    fig.legend(handles=handles, loc="lower center",
               ncol=min(4, len(handles)), fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("P3 block-length crossover: PSNR vs L (64x64, LIN link, "
                 "8-image 3-seed mean)", fontsize=12)
    fig.tight_layout(rect=(0, 0.04, 1, 0.98))
    _save(fig, os.path.join(P3, "figures", "p3_psnr_vs_L.png"), tighten=False)


# ---------------------------------------------------------------- Part 2: G1

_LINK_STYLE = {"DT30": ("-", "o"), "SAT30": ("--", "s")}
_LINK_FALLBACK = [("-.", "^"), (":", "D"), ("-", "v")]


def _g1_headroom():
    """Load g1_scale_map.csv (schema not yet frozen upstream) and return a
    tidy headroom frame with columns illum, link, n, M, headroom, or None if
    the file/columns are unusable. Assumptions documented at module top."""
    path = os.path.join(G1, "g1_scale_map.csv")
    if not os.path.exists(path):
        print("figure skip: g1_scale_map.csv not found ->", path, flush=True)
        return None
    df = pd.read_csv(path)
    df.columns = [c.strip().lower() for c in df.columns]
    if "illum" not in df.columns and "family" in df.columns:
        df = df.rename(columns={"family": "illum"})
    needed = {"illum", "link", "n", "m", "arm", "psnr"}
    if not needed.issubset(df.columns):
        print("figure skip: g1_scale_map.csv missing columns; have %s need %s"
              % (sorted(df.columns), sorted(needed)), flush=True)
        return None
    df["psnr"] = df["psnr"].astype(float)
    df["arm"] = df["arm"].astype(str).str.strip().str.lower()

    grp = df.groupby(["illum", "link", "n", "m", "arm"])["psnr"].mean()
    piv = grp.unstack("arm")
    cols = set(piv.columns)

    def pick(names):
        for nm in names:
            if nm in cols:
                return piv[nm]
        return None

    semi = pick(["semi_oracle", "semi-oracle", "semioracle", "best_semi"])
    if semi is None:
        members = [c for c in ("opg", "rmave", "opg_rmave") if c in cols]
        semi = piv[members].max(axis=1) if members else None
    lw = pick(["whiten_lw", "whiten-lw", "whitenlw"])
    iso = pick(["l_isotron", "l-isotron", "lisotron", "isotron"])
    if semi is None or (lw is None and iso is None):
        print("figure skip: g1_scale_map.csv has no recognizable semi-oracle / "
              "baseline arms (cols=%s)" % sorted(cols), flush=True)
        return None

    base_parts = [s for s in (lw, iso) if s is not None]
    base = pd.concat(base_parts, axis=1).max(axis=1)
    hd = (semi - base).reset_index(name="headroom")
    hd["m"] = hd["m"].astype(float)
    hd["n"] = hd["n"].astype(float)
    return hd


def _link_style(link, k):
    return _LINK_STYLE.get(link, _LINK_FALLBACK[k % len(_LINK_FALLBACK)])


def fig5_g1_headroom():
    """Two G1 figures if g1_scale_map.csv exists: headroom vs n at M=1e5 (per
    illum x link line) and headroom vs M/n pooled scatter (color=illum,
    marker=link)."""
    hd = _g1_headroom()
    if hd is None:
        return
    illums = sorted(hd.illum.unique())
    links = sorted(hd.link.unique())
    illum_col = {il: _ILLUM_COLOR.get(il, plt.cm.tab10(k))
                 for k, il in enumerate(illums)}

    # --- headroom vs n at M = 1e5 (fall back to max M if 1e5 absent)
    target_M = 1e5
    sel = hd[np.isclose(hd.m, target_M)]
    if sel.empty:
        target_M = hd.m.max()
        sel = hd[np.isclose(hd.m, target_M)]
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for il in illums:
        for kl, lk in enumerate(links):
            g = sel[(sel.illum == il) & (sel.link == lk)].sort_values("n")
            if g.empty:
                continue
            ls, mk = _link_style(lk, kl)
            ax.plot(g.n, g.headroom, ls=ls, marker=mk, color=illum_col[il],
                    lw=1.8, markersize=6, label="%s x %s" % (il, lk))
    ax.axhline(0.0, color="black", lw=0.8, alpha=0.6)
    ax.axhline(1.0, color="red", lw=0.9, ls="--", alpha=0.6,
               label="G1 gate (+1.0 dB)")
    ax.set_xlabel("n (patterns / resolution cells)")
    ax.set_ylabel("headroom = SEMI-ORACLE - max(WHITEN-LW, L-ISOTRON) (dB)")
    ax.set_title("G1 head-space vs n at M=%g\n(each illum x link line; "
                 "CORR-LOGN expected negative)" % target_M)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    _save(fig, os.path.join(G1, "figures", "g1_headroom_vs_n.png"))

    # --- headroom vs M/n pooled scatter (color=illum, marker=link)
    hd = hd.copy()
    hd["Mn"] = hd.m / hd.n
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for il in illums:
        for kl, lk in enumerate(links):
            g = hd[(hd.illum == il) & (hd.link == lk)]
            if g.empty:
                continue
            _, mk = _link_style(lk, kl)
            ax.scatter(g.Mn, g.headroom, color=illum_col[il], marker=mk, s=55,
                       edgecolor="black", linewidth=0.4, alpha=0.85,
                       label="%s x %s" % (il, lk))
    ax.axhline(0.0, color="black", lw=0.8, alpha=0.6)
    ax.set_xscale("log")
    ax.set_xlabel("M / n (frames per resolution cell)")
    ax.set_ylabel("headroom (dB)")
    ax.set_title("G1 head-space vs M/n (pooled over all cells)")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    _save(fig, os.path.join(G1, "figures", "g1_headroom_vs_Mn.png"))


# ---------------------------------------------------------------- shared


def _save(fig, path, tighten=True):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if tighten:
        fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main():
    for fn in (fig1_g0_permutation_collapse, fig2_g0_m_scaling,
               fig3_g0_sensitivity, fig4_p3_psnr_vs_L, fig5_g1_headroom):
        try:
            fn()
            print("figure ok:", fn.__name__, flush=True)
        except FileNotFoundError as e:
            print("figure skip (input missing):", fn.__name__, repr(e),
                  flush=True)
        except Exception as e:
            print("figure FAILED:", fn.__name__, repr(e), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
