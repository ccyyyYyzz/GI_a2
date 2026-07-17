"""ROUND63 illumination patterns — spec §3 (图样 {Bernoulli-50%, Hadamard 互补对,
GAM4 应力臂}).

All three constructions return a nonneg (M x n) float64 illumination matrix A whose
*per-pixel expectation is 1* (E[a_j] = 1 for every pixel j), so that for a
sum-normalized truth x (sum_j x_j = 1) the expected single-frame energy is
E[<a, x>] = sum_j E[a_j] x_j = 1.  This anchors the dead-time load axis
rho_bar = tau * E[lambda_i] = tau * Phi (spec §2): once the pattern carries unit
per-pixel mean, Phi alone sets the mean rate and the physical parameterisation is
decoupled from the pattern (the ROUND59 lambda = s*u coupling is NOT used here).

Kinds
-----
bern50   iid Bernoulli(0.5) in {0, 1}, scaled by 2.0 -> entries {0, 2}, E = 1.
         One physical exposure per row (exposures_per_row = 1).
hadpair  Random row/column-permuted Sylvester-Hadamard, returned as *complementary
         pairs*.  Each selected +/-1 Hadamard row h yields two physical 0/1 masks,
         (1+h)/2 and (1-h)/2, each scaled by 2.0 to entries {0, 2} (rows 1+h and
         1-h).  Both physical rows of a pair are returned in A (2*M_pairs rows,
         interleaved: rows 2k / 2k+1).  Estimator-side differencing of the pair
         b_+ - b_- (which recovers the signed +/-1 measurement) is downstream's
         job.  Per-pixel mean is *exactly* 1 for any row selection because each
         pair contributes (1+h_j) + (1-h_j) = 2 to pixel j across its two rows.
         A signed measurement costs two physical exposures (exposures_per_row = 2,
         spec §3 "按 2 次曝光计费").  n is padded up to the next power of two for
         scipy's Sylvester Hadamard; only n of the padded columns are used (random
         subset), which leaves the exact-unit-mean pair property intact.
gam4     iid Gamma(shape=4, scale=1/4), mean 4*(1/4) = 1 (spec GAM4 stress arm;
         same family used by ROUND59 Phase A/B).  exposures_per_row = 1.

RNG: rng_for(seed, 63, 1, kind_id) under the SEED0=20260717 system (spec §5,
utils.rng_for).  Stream tag 63 = ROUND63; sub-stream 1 = pattern layer;
kind_id disambiguates the three constructions.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from gi_core.utils import rng_for

R63 = 63              # ROUND63 stream tag (spec §5)
PATTERN_SUBSTREAM = 1  # pattern layer sub-stream (images63 uses sub-stream 2)

# stable kind ids for RNG derivation (frozen)
KIND_IDS = {"bern50": 1, "hadpair": 2, "gam4": 3}


def _next_pow2(n):
    """Smallest power of two >= n (scipy.linalg.hadamard needs a power-of-two
    order).  Powers of four (64^2 = 4096, 128^2 = 16384) are already powers of
    two, so the main-grid resolutions need no padding."""
    return 1 << int(np.ceil(np.log2(max(1, int(n)))))


def make_patterns(kind, M, n, seed):
    """Build an illumination matrix for one (kind, M, n, seed).

    Parameters
    ----------
    kind : {"bern50", "hadpair", "gam4"}
    M    : number of *signed* measurements (logical rows).  For hadpair this is
           the number of complementary pairs, so A has 2*M physical rows.
    n    : number of pixels (image side^2).
    seed : integer seed folded through rng_for(seed, 63, 1, kind_id).

    Returns
    -------
    dict with keys
      A                 : (R x n) float64, nonneg, per-pixel mean 1.  R = M for
                          bern50/gam4, R = 2*M for hadpair.
      exposures_per_row : physical exposures charged per signed measurement
                          (1 for bern50/gam4; 2 for hadpair pairs).
      meta              : dict describing the construction and RNG provenance.
    """
    if kind not in KIND_IDS:
        raise ValueError("unknown pattern kind %r (expected %s)"
                         % (kind, sorted(KIND_IDS)))
    M = int(M)
    n = int(n)
    kind_id = KIND_IDS[kind]
    rng = rng_for(seed, R63, PATTERN_SUBSTREAM, kind_id)
    meta = {"kind": kind, "M_signed": M, "n": n, "seed": int(seed),
            "rng_stream": [R63, PATTERN_SUBSTREAM, kind_id],
            "pixel_mean_target": 1.0, "nonneg": True}

    if kind == "bern50":
        A = 2.0 * (rng.random((M, n)) < 0.5).astype(np.float64)
        exposures_per_row = 1
        meta["construction"] = ("iid Bernoulli(0.5) in {0,1} scaled by 2.0 -> "
                                "{0,2}; per-pixel mean 1 in expectation")
        meta["n_physical_rows"] = M
        meta["total_exposures"] = M

    elif kind == "gam4":
        A = rng.gamma(shape=4.0, scale=0.25, size=(M, n)).astype(np.float64)
        exposures_per_row = 1
        meta["construction"] = ("iid Gamma(shape=4, scale=1/4); mean = 1, "
                                "Var = 1/4 (ROUND59 GAM4 stress family)")
        meta["n_physical_rows"] = M
        meta["total_exposures"] = M

    else:  # hadpair
        from scipy.linalg import hadamard

        n_had = _next_pow2(n)
        H = hadamard(n_had).astype(np.float64)          # {-1,+1}, row 0 all ones
        if M <= n_had:
            sel_rows = rng.permutation(n_had)[:M]        # distinct Hadamard rows
            row_sampling = "permutation-prefix (distinct rows)"
        else:  # more pairs than available rows -> sample with replacement
            sel_rows = rng.integers(0, n_had, size=M)
            row_sampling = "with-replacement (M > Hadamard order)"
        col_perm = rng.permutation(n_had)[:n]            # random n-column subset
        Hsub = H[np.ix_(sel_rows, col_perm)]             # (M, n) in {-1,+1}
        A = np.empty((2 * M, n), dtype=np.float64)
        A[0::2] = 1.0 + Hsub                             # 2*(1+h)/2, positive mask
        A[1::2] = 1.0 - Hsub                             # 2*(1-h)/2, complementary
        exposures_per_row = 2
        meta["construction"] = (
            "random row/col-permuted Sylvester-Hadamard; each signed row h -> "
            "interleaved physical pair (1+h, 1-h) in {0,2}; per-pixel mean "
            "EXACTLY 1 by pair construction")
        meta["n_had"] = n_had
        meta["padded"] = bool(n_had > n)
        meta["row_sampling"] = row_sampling
        meta["n_physical_rows"] = 2 * M
        meta["total_exposures"] = 2 * M
        meta["pair_layout"] = "interleaved: physical rows (2k, 2k+1) form pair k"
        meta["pair_indices"] = [[2 * k, 2 * k + 1] for k in range(M)]

    meta["exposures_per_row"] = exposures_per_row
    return {"A": A, "exposures_per_row": exposures_per_row, "meta": meta}


def _smoke():
    """Build each kind at n=4096, M=1024, check unit per-pixel mean, print shapes."""
    n, M = 4096, 1024
    print("[patterns smoke] n=%d M=%d seed=0" % (n, M), flush=True)
    all_ok = True
    for kind in ("bern50", "hadpair", "gam4"):
        d = make_patterns(kind, M, n, 0)
        A = d["A"]
        pm = A.mean(axis=0)                 # per-pixel (per-column) empirical mean
        agg = float(pm.mean())              # aggregate per-pixel mean == global mean
        # NOTE: the unit-mean spec is an ENSEMBLE property (E[a_j] = 1). For the iid
        # kinds individual columns carry sqrt sampling noise (~1/sqrt(M) ~ 3%/1.6%
        # for bern50/gam4 at M=1024), so the 2% tolerance is asserted on the
        # aggregate per-pixel mean; the per-column spread is printed for context.
        # hadpair is exact to machine precision.
        ok = abs(agg - 1.0) < 0.02 and bool((A >= 0.0).all())
        all_ok = all_ok and ok
        print("  %-8s A=%s exp/row=%d  pixmean agg=%.5f "
              "[min=%.4f max=%.4f std=%.4f] nonneg=%s  %s"
              % (kind, A.shape, d["exposures_per_row"], agg,
                 float(pm.min()), float(pm.max()), float(pm.std()),
                 bool((A >= 0.0).all()), "PASS" if ok else "FAIL"))
    print("[patterns smoke] %s" % ("ALL PASS" if all_ok else "FAILURES ABOVE"),
          flush=True)
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(_smoke())
