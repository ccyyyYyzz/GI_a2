"""ROUND63 main-figure panel (a): the (rho, nu) count-information map.

Spec: docs/ROUND63_SPEC_D2.md section 6(a) and the theory of section 5 /
paper/THEORY_SECTION_DRAFT.md sections X.2-X.4.

Renders panel (a) of the four-up main figure: filled contours of the exact
per-slot count information J_exact(rho, nu) about log-lambda, overlaid with

  * the principal information ridge rho*(nu) = argmax_rho J_exact (solid),
  * its closed form (6 nu)^(1/3) - 2/3 (dashed),
  * the 10% CLT-discrepancy boundary rho_0.9(nu) where J_exact/J_CLT = 0.9
    (dash-dot),
  * the RQL deployment band rho <= 1 (per-pattern rho_95 <= 1, shaded),
  * the deployment zoning of THEORY_SECTION_DRAFT section X.4, and
  * the conventional (rho = 0.05) and primary high-flux (rho = 0.6) operating
    points.

Convention (matched to code/round63/fisher_ridge.py EXACTLY):
  J_exact(rho, nu) = exact_fisher_analytic(lam=rho, T=nu, tau=1) * tau / T
                   = exact_fisher_analytic(rho, nu, 1) / nu
which equals fisher_ridge.py's `fisher_exact(lam,T,tau)*lam**2*tau/T` because
exact_fisher_analytic already returns the information about theta = log(lam)
per frame (i.e. it carries the lam**2 log-lambda Jacobian internally). The CLT
surrogate is J_CLT(rho) = rho/(1+rho), independent of nu. Both are cross-checked
against results/round63_theory/fisher_ridge.csv at load time (assertion).

This is a NEW, self-contained figure generator: it imports the frozen physics
kernel `exact_fisher_analytic` but modifies no existing file.

Usage:
  python fig_a_ridge_map.py            # use cached sweep if present
  python fig_a_ridge_map.py --force    # recompute the dense sweep
"""
import argparse
import csv
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from physics import exact_fisher_analytic  # frozen kernel; not modified

ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(ROOT, "results", "round63_theory")
NPZ = os.path.join(OUT, "fig_a_sweep.npz")
CSV = os.path.join(OUT, "fisher_ridge.csv")
JSON = os.path.join(OUT, "fisher_ridge.json")
PDF = os.path.join(OUT, "fig_a_ridge_map.pdf")
PNG = os.path.join(OUT, "fig_a_ridge_map.png")

# --- sweep grid (spec 6a: nu log 20..2000 >= 40 pts, rho log 0.01..64 >= 80)
TAU = 1.0                                   # scale-free: everything in tau units
ANCHORS = [20, 50, 100, 200, 500, 1000, 2000]   # fisher_ridge.json nu anchors
# nu: dense log grid, unioned with the json anchors so the assertions/table
# land on exact nu values.
NU_GRID = np.unique(np.concatenate([np.geomspace(20.0, 2000.0, 44),
                                    np.asarray(ANCHORS, float)]))
# rho: 481 log points -> adjacent spacing factor ~1.0184 (1.84% < 2%), so the
# raw argmax ridge has < 2% jitter (spec requirement 5).
RHO_GRID = np.geomspace(0.01, 64.0, 481)

# operating points and zoning
RHO_CONV = 0.05        # conventional photon-counting operating point
RHO_HIFLUX = 0.6       # campaign primary high-flux point (pre-registered)
RHO_DEPLOY = 1.0       # RQL deployment ceiling (per-pattern rho_95 <= 1)

# Wong colorblind-safe palette
C_RIDGE = "#000000"    # exact argmax ridge
C_CLOSED = "#D55E00"   # (6nu)^(1/3) - 2/3
C_CLT = "#0072B2"      # rho_0.9 CLT-10% boundary
C_CONV = "#009E73"     # rho = 0.05 line
C_HIFLUX = "#CC79A7"   # rho = 0.6 line


def j_exact(rho, nu, tau=TAU):
    """Per-slot exact count information about log-lambda (fisher_ridge convention)."""
    T = nu * tau
    return exact_fisher_analytic(rho, T, tau) * tau / T


def j_clt(rho):
    """CLT (Gaussian-renewal) per-slot information about log-lambda = rho/(1+rho)."""
    return rho / (1.0 + rho)


# ----------------------------------------------------------------------
# dense sweep (cached)
# ----------------------------------------------------------------------
def compute_sweep():
    t0 = time.time()
    nrho, nnu = RHO_GRID.size, NU_GRID.size
    J = np.empty((nrho, nnu), dtype=np.float64)
    for j, nu in enumerate(NU_GRID):
        T = nu * TAU
        for i, rho in enumerate(RHO_GRID):
            J[i, j] = exact_fisher_analytic(rho, T, TAU) * TAU / T
    dt = time.time() - t0
    Jclt = j_clt(RHO_GRID)                       # (nrho,) nu-independent
    ratio = J / Jclt[:, None]
    return {"nu_grid": NU_GRID, "rho_grid": RHO_GRID, "J": J,
            "J_clt": Jclt, "ratio": ratio, "sweep_seconds": np.array(dt)}


def load_or_compute(force=False):
    if (not force) and os.path.exists(NPZ):
        d = np.load(NPZ)
        if (d["rho_grid"].shape == RHO_GRID.shape
                and d["nu_grid"].shape == NU_GRID.shape
                and np.allclose(d["rho_grid"], RHO_GRID)
                and np.allclose(d["nu_grid"], NU_GRID)):
            print("Loaded cached sweep: %s" % NPZ)
            return {k: d[k] for k in d.files}
        print("Cached grid mismatch -> recomputing.")
    print("Computing dense sweep (%d rho x %d nu = %d cells)..."
          % (RHO_GRID.size, NU_GRID.size, RHO_GRID.size * NU_GRID.size))
    s = compute_sweep()
    os.makedirs(OUT, exist_ok=True)
    np.savez_compressed(NPZ, **s)
    print("Sweep cached to %s (%.2f s)" % (NPZ, float(s["sweep_seconds"])))
    return s


# ----------------------------------------------------------------------
# cross-check vs the frozen fisher_ridge.csv artifact (spec req 1)
# ----------------------------------------------------------------------
def spotcheck_against_csv(tol=1e-6):
    """Reproduce 3 interior csv rows of I_exact_log to <tol relative."""
    rows = []
    with open(CSV, newline="") as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            rows.append((float(row[0]), float(row[1]), float(row[2])))  # nu,rho,Ie
    # one interior row (rho closest to 5) per selected nu
    picks = []
    for nu_sel in (100.0, 500.0, 2000.0):
        cand = [rw for rw in rows if rw[0] == nu_sel]
        rw = min(cand, key=lambda t: abs(t[1] - 5.0))
        picks.append(rw)
    print("\nSpot-check vs fisher_ridge.csv (I_exact_log, log-lambda per slot):")
    print("   nu       rho     csv_value     recomputed     rel_err")
    worst = 0.0
    for nu, rho, ie in picks:
        rec = j_exact(rho, nu)
        rel = abs(rec - ie) / abs(ie)
        worst = max(worst, rel)
        print("%6.0f  %8.4f  %.9e  %.9e  %.2e" % (nu, rho, ie, rec, rel))
    assert worst < tol, "csv spot-check rel err %.2e exceeds %.1e" % (worst, tol)
    print("Spot-check OK: max rel err %.2e < %.1e" % (worst, tol))


# ----------------------------------------------------------------------
# ridge + boundary extraction
# ----------------------------------------------------------------------
def ridge_argmax(J):
    """Raw argmax over rho at each nu (spec: 'argmax over rho from the sweep')."""
    return RHO_GRID[np.argmax(J, axis=0)]


def closed_form_ridge(nu):
    return (6.0 * nu) ** (1.0 / 3.0) - 2.0 / 3.0


def rho_09(ratio):
    """First rho where J_exact/J_CLT crosses 0.9 (log-rho linear interpolation)."""
    lr = np.log(RHO_GRID)
    out = np.full(NU_GRID.size, np.nan)
    for j in range(NU_GRID.size):
        r = ratio[:, j]
        below = np.nonzero(r < 0.9)[0]
        if below.size == 0 or below[0] == 0:
            continue
        k = below[0]
        r0, r1 = r[k - 1], r[k]                 # r0 >= 0.9 > r1
        t = (r0 - 0.9) / (r0 - r1)
        out[j] = float(np.exp(lr[k - 1] + t * (lr[k] - lr[k - 1])))
    return out


# ----------------------------------------------------------------------
# figure
# ----------------------------------------------------------------------
def make_figure(s):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import patheffects as pe
    from matplotlib.ticker import LogLocator, LogFormatterMathtext

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "Times New Roman", "STIXGeneral"],
        "mathtext.fontset": "dejavuserif",
        "font.size": 8,
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 6.8,
        "axes.linewidth": 0.7,
        "lines.antialiased": True,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })

    J = s["J"]
    ratio = s["ratio"]
    ridge = ridge_argmax(J)
    rho09 = rho_09(ratio)
    cf = closed_form_ridge(NU_GRID)

    NU, RHO = np.meshgrid(NU_GRID, RHO_GRID)     # (nrho, nnu)

    fig, ax = plt.subplots(figsize=(3.35, 3.25))

    # filled contours of J_exact (ceiling = 1.0)
    levels = np.linspace(0.0, 1.0, 21)
    cf_fill = ax.contourf(NU, RHO, J, levels=levels, cmap="viridis",
                          vmin=0.0, vmax=1.0, extend="neither")
    # a few reference iso-information lines
    cl = ax.contour(NU, RHO, J, levels=[0.3, 0.5, 0.7, 0.9], colors="white",
                    linewidths=0.45, alpha=0.55)
    ax.clabel(cl, fmt="%.1f", fontsize=5.6, inline=True, inline_spacing=2)

    halo = [pe.withStroke(linewidth=2.6, foreground="white")]

    # principal ridge: exact argmax (solid) + closed form (dashed)
    ax.plot(NU_GRID, ridge, color=C_RIDGE, lw=2.0, solid_capstyle="round",
            path_effects=halo, label=r"$\rho^{*}(\nu)$ exact argmax", zorder=6)
    ax.plot(NU_GRID, cf, color=C_CLOSED, lw=1.7, ls=(0, (5, 2)),
            path_effects=halo, label=r"$(6\nu)^{1/3}-2/3$", zorder=6)
    # CLT 10% discrepancy boundary
    ax.plot(NU_GRID, rho09, color=C_CLT, lw=1.7, ls=(0, (4, 1.5, 1, 1.5)),
            path_effects=halo, label=r"$\rho_{0.9}$ (CLT $10\%$)", zorder=6)

    # RQL deployment band rho <= 1 (shaded)
    ax.axhspan(RHO_GRID.min(), RHO_DEPLOY, color="white", alpha=0.12, zorder=2)
    ax.axhline(RHO_DEPLOY, color="white", lw=0.8, alpha=0.7, zorder=3)

    # operating-point reference lines
    ax.axhline(RHO_CONV, color=C_CONV, lw=1.1, ls=(0, (1, 1)),
               path_effects=halo, zorder=5)
    ax.axhline(RHO_HIFLUX, color=C_HIFLUX, lw=1.3, ls=(0, (1, 1)),
               path_effects=halo, zorder=5)

    txt_halo = [pe.withStroke(linewidth=1.8, foreground="white")]

    def zone(x, y, t, **kw):
        ax.text(x, y, t, ha="center", va="center", fontsize=6.4,
                color="black", path_effects=txt_halo, zorder=7, **kw)

    # zoning annotations (THEORY_SECTION_DRAFT X.4)
    zone(250, 0.17, "RQL deployment\n" + r"(per-pattern $\rho_{95}\leq 1$)")
    zone(110, 2.6, "transition")
    zone(1250, 27.0, "information-decreasing\n" + r"($\rho \geq \rho^{*}$)")
    zone(31, 23, "exact-reference\n(short-window &\nextreme saturation)")

    # operating-point labels (left side, clear of the lower-right legend)
    xL = NU_GRID.min() * 1.06
    ax.text(xL, RHO_CONV * 1.16, r"$\rho=0.05$ conventional", ha="left",
            va="bottom", fontsize=6.2, color=C_CONV, path_effects=txt_halo,
            zorder=8)
    ax.text(xL, RHO_HIFLUX * 1.10, r"$\rho=0.6$ high-flux", ha="left",
            va="bottom", fontsize=6.2, color=C_HIFLUX, path_effects=txt_halo,
            zorder=8)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(NU_GRID.min(), NU_GRID.max())
    ax.set_ylim(RHO_GRID.min(), RHO_GRID.max())
    ax.set_xlabel(r"dead-time slots per frame  $\nu = T/\tau$")
    ax.set_ylabel(r"per-slot load  $\rho = \lambda\tau$")
    ax.xaxis.set_major_locator(LogLocator(base=10.0))
    ax.xaxis.set_major_formatter(LogFormatterMathtext(base=10.0))
    ax.tick_params(which="both", direction="out", length=2.6, width=0.6)

    cb = fig.colorbar(cf_fill, ax=ax, pad=0.02, fraction=0.046,
                      ticks=np.linspace(0, 1, 6))
    cb.set_label(r"$J_{\mathrm{exact}} = I_N/\nu$  (per-slot info about "
                 r"$\log\lambda$; ceiling $1$)", fontsize=8)
    cb.ax.tick_params(labelsize=7.5, length=2.2, width=0.5)
    cb.outline.set_linewidth(0.6)

    leg = ax.legend(loc="lower right", framealpha=0.92, handlelength=1.9,
                    borderpad=0.4, labelspacing=0.3, edgecolor="0.4")
    leg.get_frame().set_linewidth(0.5)

    fig.tight_layout(pad=0.4)
    fig.savefig(PDF)
    fig.savefig(PNG, dpi=300)
    plt.close(fig)
    return ridge, cf, rho09


# ----------------------------------------------------------------------
def report_and_assert(ridge):
    """Table of rho* sweep-vs-closed-form + assertion vs fisher_ridge.json."""
    jd = json.load(open(JSON))["ridge"]
    grid_factor = (RHO_GRID[-1] / RHO_GRID[0]) ** (1.0 / (RHO_GRID.size - 1))
    log_tol = np.log(grid_factor)
    print("\nRidge rho*(nu): sweep argmax vs closed form (6nu)^(1/3)-2/3")
    print("   nu   sweep_argmax   closed_form   json_rho_star   argmax/json")
    worst = 0.0
    for nu in (20, 100, 500, 2000):
        k = int(np.where(np.isclose(NU_GRID, nu))[0][0])
        sweep = ridge[k]
        cf = closed_form_ridge(float(nu))
        js = jd[str(nu)]["rho_star"]
        r = sweep / js
        worst = max(worst, abs(np.log(r)))
        print("%6d   %10.4f    %10.4f    %10.4f     %.4f"
              % (nu, sweep, cf, js, r))
    print("rho-grid spacing factor = %.5f (<= %.2f%% steps)"
          % (grid_factor, (grid_factor - 1) * 100))
    # assert every anchor's raw argmax matches json within one grid spacing
    for nu in ANCHORS:
        k = int(np.where(np.isclose(NU_GRID, nu))[0][0])
        js = jd[str(nu)]["rho_star"]
        dlog = abs(np.log(ridge[k] / js))
        assert dlog <= log_tol + 1e-12, (
            "nu=%d argmax %.4f vs json %.4f differs by %.4f > grid spacing %.4f"
            % (nu, ridge[k], js, dlog, log_tol))
    print("Assertion OK: every anchor argmax matches fisher_ridge.json "
          "within one rho-grid spacing (all 7 nu).")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true",
                    help="recompute the dense sweep even if the cache exists")
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    spotcheck_against_csv()
    s = load_or_compute(force=args.force)
    ridge, _, _ = make_figure(s)
    report_and_assert(ridge)
    print("\nWrote:\n  %s\n  %s" % (PDF, PNG))


if __name__ == "__main__":
    main()
