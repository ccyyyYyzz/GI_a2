"""Dead-time-aware optimal-experiment-design (OED) solver for single-pixel
imaging illumination patterns  (ROUND63, GI_another).

Design goal
-----------
Choose a set of k-sparse non-negative illumination rows (the SPI patterns) that
maximize the D-optimality criterion  log det M(xi)  for reconstructing an
n = side*side scene on the torus, under the campaign's **equal-detector-load**
rule.  Because every design atom is load-normalized to the SAME expected load
r.xhat = 1, the dead-time throttling kernel J(load) is a single constant factor
that multiplies every atom's Fisher contribution identically; it therefore
cancels out of every argmax and out of every log-det RATIO, so it is dropped
here (this is the "dead-time-aware" reduction -- equal load makes the design
problem dead-time-invariant).

Method
------
Frank-Wolfe / Fedorov-Wynn vertex-exchange for D-optimality with the
Kiefer-Wolfowitz (KW) equivalence-theorem certificate.  Atoms are k=16 sparse,
so the variance function d(a) = a^T M^{-1} a costs O(k^2) per atom (gather the
16x16 sub-block of M^{-1} and contract).  A full n x n inverse is refactored
every step (n = 1024 -> fast and numerically safe; correctness over speed).

Only numpy is imported.  Run as __main__ for the spec self-test.
"""

import numpy as np


# --------------------------------------------------------------------------- #
#  Default 16-pixel shapes (offset sets on the side x side torus)             #
# --------------------------------------------------------------------------- #
def _default_shapes(side=32):
    """Return the six default 16-pixel offset-set shapes."""
    shapes = {}

    # scatter16 : a fixed seeded 16-pixel scattered set (seed 63)
    rs = np.random.RandomState(63)
    idx = rs.choice(side * side, size=16, replace=False)
    shapes["scatter16"] = [(int(i // side), int(i % side)) for i in idx]

    # Lblob6x6 : the campaign's L-shaped blob
    shapes["Lblob6x6"] = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 0),
                          (2, 1), (3, 0), (3, 1), (3, 2), (4, 1), (4, 2),
                          (4, 3), (5, 2), (5, 3), (5, 4)]

    # solid4x4 : a filled 4x4 block
    shapes["solid4x4"] = [(i, j) for i in range(4) for j in range(4)]

    # rect2x8 : a 2-row x 8-col bar
    shapes["rect2x8"] = [(i, j) for i in range(2) for j in range(8)]

    # blocks2x2x4 : four 2x2 blocks at (0,0),(0,6),(6,0),(6,6)
    blk = []
    for (by, bx) in [(0, 0), (0, 6), (6, 0), (6, 6)]:
        for dy in range(2):
            for dx in range(2):
                blk.append((by + dy, bx + dx))
    shapes["blocks2x2x4"] = blk

    # ring5 : 16 pixels closest to a radius-2.5 ring inside a 6x6 window
    pts = [(i, j) for i in range(6) for j in range(6)]
    c = 2.5
    resid = np.array([abs(np.hypot(i - c, j - c) - 2.5) for (i, j) in pts])
    order = np.argsort(resid, kind="stable")[:16]          # deterministic
    shapes["ring5"] = [pts[int(o)] for o in order]

    for name, offs in shapes.items():
        assert len(offs) == 16, (name, len(offs))
        assert len(set(offs)) == 16, ("duplicate offsets in", name)
    return shapes


# --------------------------------------------------------------------------- #
#  Translation-support index builder                                          #
# --------------------------------------------------------------------------- #
def _shape_support(offsets, side):
    """All `side*side` cyclic translations of an offset set.

    Returns int array (n_trans, k) of flat pixel indices; translation t = sy*side+sx.
    """
    sy = np.repeat(np.arange(side), side)          # (n,)
    sx = np.tile(np.arange(side), side)            # (n,)
    oy = np.array([o[0] for o in offsets])[None, :]  # (1,k)
    ox = np.array([o[1] for o in offsets])[None, :]
    py = (oy + sy[:, None]) % side                 # (n,k)
    px = (ox + sx[:, None]) % side
    return (py * side + px).astype(np.int64)


# --------------------------------------------------------------------------- #
#  Dictionary of load-normalized design atoms                                 #
# --------------------------------------------------------------------------- #
def _build_dictionary(xhat, shapes, powers, side, cv_guard=1.0):
    """Assemble every (shape, translation, power) atom, load-normalized to load 1.

    Guard: drop atoms whose weight dispersion sd(w_on)/mean(w_on) > cv_guard.
    Returns IDX (A,k) int, VAL (A,k) float, and parallel meta arrays
    (shape names list, translation t int array, power float array).
    """
    n = side * side
    mean_x = xhat.mean()
    IDX_parts, VAL_parts, name_parts, t_parts, p_parts = [], [], [], [], []
    t_all = np.arange(n)

    for name, offs in shapes.items():
        sidx = _shape_support(offs, side)          # (n,k)
        xh_on = xhat[sidx]                          # (n,k)  intensity on support
        for p in powers:
            w = (xh_on / mean_x) ** float(p)        # (n,k)  illumination weights
            load = (w * xh_on).sum(axis=1)          # (n,)   raw r . xhat
            val = w / load[:, None]                 # scale so scaled_r . xhat = 1
            # weight-dispersion guard on the support
            wm = w.mean(axis=1)
            wsd = w.std(axis=1)
            cv = np.where(wm > 0, wsd / wm, 0.0)
            keep = cv <= cv_guard
            if not keep.any():
                continue
            IDX_parts.append(sidx[keep])
            VAL_parts.append(val[keep])
            name_parts.extend([name] * int(keep.sum()))
            t_parts.append(t_all[keep])
            p_parts.append(np.full(int(keep.sum()), float(p)))

    IDX = np.concatenate(IDX_parts, axis=0)
    VAL = np.concatenate(VAL_parts, axis=0)
    t_arr = np.concatenate(t_parts, axis=0)
    p_arr = np.concatenate(p_parts, axis=0)
    return IDX, VAL, name_parts, t_arr, p_arr


# --------------------------------------------------------------------------- #
#  Information matrix from a sparse dictionary + weight vector                 #
# --------------------------------------------------------------------------- #
def _info_matrix(IDX, VAL, xi, n):
    """M_design = sum_r xi_r a_r a_r^T  (eps*I added by caller). Scatter via bincount."""
    A, k = IDX.shape
    lin = (IDX[:, :, None] * n + IDX[:, None, :]).ravel()          # (A*k*k,)
    contrib = (xi[:, None, None] * VAL[:, :, None] * VAL[:, None, :]).ravel()
    return np.bincount(lin, weights=contrib, minlength=n * n).reshape(n, n)


# --------------------------------------------------------------------------- #
#  Main solver                                                                #
# --------------------------------------------------------------------------- #
def design(xhat, shapes=None, powers=(0, 0.5, 1.0), M_rows=1024,
           max_iter=400, tol=0.05):
    """Dead-time-aware D-optimal SPI illumination design.

    Parameters
    ----------
    xhat    : (n,) non-negative scene prior (n must be a perfect square).
    shapes  : dict {name: [(oy,ox), ...16]}   (default: six built-in shapes).
    powers  : intensity-matching weight powers p.
    M_rows  : number of discrete patterns produced by rounding.
    max_iter, tol : FW / KW-certificate budget and gap tolerance.

    Returns dict with keys A, atoms, gap, gap_traj, cert.
    """
    xhat = np.asarray(xhat, dtype=float).ravel()
    n = xhat.size
    side = int(round(np.sqrt(n)))
    if side * side != n:
        raise ValueError("xhat length must be a perfect square (torus side).")
    if xhat.min() < 0:
        raise ValueError("xhat must be non-negative.")
    if shapes is None:
        shapes = _default_shapes(side)

    # ---- dictionary ------------------------------------------------------- #
    IDX, VAL, names, t_arr, p_arr = _build_dictionary(xhat, shapes, powers, side)
    A = IDX.shape[0]

    # ---- initial (uniform) design + rank-guard eps ------------------------ #
    xi = np.full(A, 1.0 / A)
    Mdes = _info_matrix(IDX, VAL, xi, n)
    trace_scale = np.trace(Mdes) / n
    eps = 1e-9 * trace_scale
    I = np.eye(n)

    # ---- Frank-Wolfe / Fedorov-Wynn loop ---------------------------------- #
    gap_traj = []
    iters = 0
    for it in range(max_iter):
        iters = it + 1
        Minv = np.linalg.inv(Mdes + eps * I)
        # variance function d(a) = a^T M^{-1} a  for every atom (O(k^2) each)
        sub = Minv[IDX[:, :, None], IDX[:, None, :]]           # (A,k,k)
        d = np.einsum('ai,aij,aj->a', VAL, sub, VAL)           # (A,)
        astar = int(np.argmax(d))
        dstar = float(d[astar])
        gap = dstar / n - 1.0
        gap_traj.append(gap)
        if gap < tol:
            break
        # optimal D-optimal FW step, clipped to (0, 0.3]
        gamma = (dstar / n - 1.0) / (dstar - 1.0)
        gamma = min(max(gamma, 1e-6), 0.3)
        # M_design <- (1-g) M_design + g a* a*^T   (rank-1, 16-sparse)
        Mdes *= (1.0 - gamma)
        v = VAL[astar]
        ix = IDX[astar]
        Mdes[np.ix_(ix, ix)] += gamma * np.outer(v, v)
        # keep xi consistent with M_design
        xi *= (1.0 - gamma)
        xi[astar] += gamma

    final_gap = gap_traj[-1]

    # ---- exact M(xi*) and its log-det (continuous optimum) ---------------- #
    Mdes_exact = _info_matrix(IDX, VAL, xi, n)
    _, logdet_opt = np.linalg.slogdet(Mdes_exact + eps * I)

    # ---- rounding: largest-remainder apportionment over the support ------- #
    thr = xi.max() * 1e-6
    supp = np.where(xi > thr)[0]
    w_s = xi[supp]
    w_s = w_s / w_s.sum()
    quota = w_s * M_rows
    base = np.floor(quota).astype(np.int64)
    rem = int(M_rows - base.sum())
    frac = quota - base
    # deterministic tie-break: largest remainder first, then lowest atom index
    order = np.lexsort((supp, -frac))
    for i in range(rem):
        base[order[i]] += 1
    counts = base

    # ---- assemble discrete pattern matrix A and atom list ----------------- #
    Aout = np.zeros((M_rows, n), dtype=float)
    atoms = []
    row = 0
    for j, atom in enumerate(supp):
        c = int(counts[j])
        if c == 0:
            continue
        pattern = np.zeros(n, dtype=float)
        pattern[IDX[atom]] = VAL[atom]
        Aout[row:row + c] = pattern
        atoms.append((names[atom], int(t_arr[atom]), float(p_arr[atom]), c))
        row += c
    assert row == M_rows, (row, M_rows)

    cert = {"n": n,
            "final_gap": float(final_gap),
            "dict_size": int(A),
            "iters": int(iters),
            "eps": float(eps),
            "logdet": float(logdet_opt),
            "converged": bool(final_gap < tol)}

    return {"A": Aout,
            "atoms": atoms,
            "gap": float(final_gap),
            "gap_traj": [float(g) for g in gap_traj],
            "cert": cert}


# --------------------------------------------------------------------------- #
#  Self-test                                                                   #
# --------------------------------------------------------------------------- #
def _uniform_family_logdet(offsets, xhat, side, eps):
    """log det M for the pure, uniform-weight (p=0), load-normalized translate
    family of one shape -- used as a baseline in the self-test."""
    n = side * side
    sidx = _shape_support(offsets, side)               # (n,k)
    load = xhat[sidx].sum(axis=1)                       # p=0 -> w=1
    val = np.ones_like(sidx, dtype=float) / load[:, None]
    xi = np.full(n, 1.0 / n)
    Mdes = _info_matrix(sidx, val, xi, n)
    _, ld = np.linalg.slogdet(Mdes + eps * np.eye(n))
    return ld


if __name__ == "__main__":
    import time
    from collections import Counter

    SIDE, N = 32, 1024

    # ----- xhat : 4x4-block-blurred synthetic scene (matches dev_pilot_L3v2) #
    yy, xx = np.mgrid[0:SIDE, 0:SIDE]
    scene = np.ones((SIDE, SIDE))
    for (cy, cx) in [(8, 9), (10, 23), (22, 7), (24, 21)]:
        scene += 0.6 * np.exp(-(((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * 3.5 ** 2)))
    x = (scene / scene.sum()).ravel()
    xi_img = x.reshape(SIDE, SIDE)
    coarse = xi_img.reshape(SIDE // 4, 4, SIDE // 4, 4).mean(axis=(1, 3))
    xhat = np.repeat(np.repeat(coarse, 4, axis=0), 4, axis=1).ravel()
    print("[selftest] xhat: n=%d  min=%.3e  max=%.3e  sum=%.6f  cv=%.3f"
          % (xhat.size, xhat.min(), xhat.max(), xhat.sum(),
             xhat.std() / xhat.mean()), flush=True)

    # ----- (1) run the solver --------------------------------------------- #
    t0 = time.time()
    out = design(xhat)
    wall = time.time() - t0
    cert = out["cert"]
    print("\n[selftest] === (1) KW certificate ===", flush=True)
    print("  final KW gap        = %.6f  (tol 0.05, %s)"
          % (out["gap"], "REACHED" if cert["converged"] else "iter-capped"))
    print("  iterations          = %d" % cert["iters"])
    print("  dictionary size     = %d atoms" % cert["dict_size"])
    print("  gap_traj[:5]        = %s" % ["%.4f" % g for g in out["gap_traj"][:5]])
    print("  gap_traj[-5:]       = %s" % ["%.4f" % g for g in out["gap_traj"][-5:]])

    print("\n[selftest] atom composition histogram (shape, p): "
          "n_atoms / total_rows", flush=True)
    n_atoms = Counter()
    n_rows = Counter()
    for (name, t, p, c) in out["atoms"]:
        n_atoms[(name, p)] += 1
        n_rows[(name, p)] += c
    for key in sorted(n_atoms):
        print("  %-12s p=%.1f : %3d atoms / %4d rows"
              % (key[0], key[1], n_atoms[key], n_rows[key]))
    print("  TOTAL: %d distinct atoms / %d rows"
          % (len(out["atoms"]), sum(n_rows.values())))

    # ----- (2) verify every pattern row ----------------------------------- #
    Aout = out["A"]
    print("\n[selftest] === (2) pattern-row checks ===", flush=True)
    min_val = Aout.min()
    nnz = (Aout != 0).sum(axis=1)
    loads = Aout @ xhat
    load_rel = (loads.max() - loads.min()) / loads.mean()
    ok_nonneg = min_val >= 0.0
    ok_sparse = int(nnz.max()) <= 16
    ok_load = load_rel <= 1e-9
    print("  nonneg              : min entry = %.3e            -> %s"
          % (min_val, "PASS" if ok_nonneg else "FAIL"))
    print("  <=16 nonzeros/row   : max nnz   = %d               -> %s"
          % (int(nnz.max()), "PASS" if ok_sparse else "FAIL"))
    print("  equal load r.xhat   : load=%.6f  rel-spread=%.2e -> %s"
          % (loads.mean(), load_rel, "PASS" if ok_load else "FAIL"))

    # ----- (3) log-det sanity vs single-family baselines ------------------ #
    shp = _default_shapes(SIDE)
    ld_scatter = _uniform_family_logdet(shp["scatter16"], xhat, SIDE, cert["eps"])
    ld_lblob = _uniform_family_logdet(shp["Lblob6x6"], xhat, SIDE, cert["eps"])
    ld_opt = cert["logdet"]
    print("\n[selftest] === (3) log det M sanity (same eps=%.3e) ==="
          % cert["eps"], flush=True)
    print("  log det M(xi*)  optimized design      = %.4f" % ld_opt)
    print("  log det M       pure scattered (unif)  = %.4f" % ld_scatter)
    print("  log det M       pure Lblob     (unif)  = %.4f" % ld_lblob)
    ok_ld = (ld_opt >= ld_scatter - 1e-6) and (ld_opt >= ld_lblob - 1e-6)
    print("  optimized >= both baselines           -> %s"
          % ("PASS" if ok_ld else "FAIL"))

    # ----- (4) wall time -------------------------------------------------- #
    print("\n[selftest] === (4) wall time ===", flush=True)
    ok_time = wall < 300.0
    print("  design() wall time  = %.1f s  (< 300 s) -> %s"
          % (wall, "PASS" if ok_time else "FAIL"))

    all_ok = ok_nonneg and ok_sparse and ok_load and ok_ld and ok_time
    print("\n[selftest] RESULT: %s" % ("ALL CHECKS PASS" if all_ok
                                       else "*** FAILURE ***"), flush=True)
