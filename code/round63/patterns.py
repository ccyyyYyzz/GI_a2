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
sparsek  STUDY-2 sparse illumination ladder (GPT round-8 ruling §3/§4). A binary
         occupancy-k measurement matrix B^(k) (every row and every column sum k,
         full numeric rank) is built from k distinct cyclic shifts of a seeded
         base permutation and then randomized by 20n accepted degree-preserving
         2-switches; the matrix handed to the campaign is the PHOTON-BUDGET-
         NORMALIZED A^(k) = (n/k) * B^(k), so that for a sum-normalized truth the
         mean detector load is identical across occupancy levels (E[<a,x>] = 1;
         ruling §4).  Requires the Study-2 square geometry M = n (both margins
         then equal k exactly).  exposures_per_row = 1.  k = 1 is a seeded
         permutation matrix (raster / point scanning).

RNG: rng_for(seed, 63, 1, kind_id) under the SEED0=20260717 system (spec §5,
utils.rng_for).  Stream tag 63 = ROUND63; sub-stream 1 = pattern layer;
kind_id disambiguates the constructions.  sparsek additionally folds (k, nonce)
into the stream so different occupancy levels and rank-rebuild nonces are
independent.
"""
import hashlib
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from gi_core.utils import rng_for

R63 = 63              # ROUND63 stream tag (spec §5)
PATTERN_SUBSTREAM = 1  # pattern layer sub-stream (images63 uses sub-stream 2)

# stable kind ids for RNG derivation (frozen)
KIND_IDS = {"bern50": 1, "hadpair": 2, "gam4": 3, "sparsek": 4}

# sparsek: number of accepted degree-preserving 2-switches per pixel (ruling §3
# "apply 20n = 20480 accepted degree-preserving 2-switches"); n = 1024 -> 20480.
SPARSEK_SWITCHES_PER_N = 20


def _next_pow2(n):
    """Smallest power of two >= n (scipy.linalg.hadamard needs a power-of-two
    order).  Powers of four (64^2 = 4096, 128^2 = 16384) are already powers of
    two, so the main-grid resolutions need no padding."""
    return 1 << int(np.ceil(np.log2(max(1, int(n)))))


def _sparsek_matrix(M, n, seed, k):
    """Build B^(k) for the sparse ladder (ruling §3), returning
    (B uint8, nonce, cond, n_switch_accepted, n_switch_attempts).

    Construction: a seeded base permutation pi is stacked with k DISTINCT cyclic
    column shifts {s_1..s_k}, i.e. B[i, (pi(i)+s_t) mod n] = 1.  With distinct
    shifts every row carries exactly k ones and — because pi is a bijection —
    every column also carries exactly k ones (one per shift).  The matrix is then
    randomized by 20n ACCEPTED degree-preserving 2-switches on its 1-edge list:
    two ones (i1,j1),(i2,j2) with i1!=i2, j1!=j2 and B[i1,j2]=B[i2,j1]=0 are
    swapped to (i1,j2),(i2,j1); this preserves every row AND column sum and keeps
    entries binary (edge-swap form of the checkerboard switch — high acceptance,
    unlike blind 4-tuple rejection).  numpy matrix_rank on float64 must be full
    (n); if deficient a preregistered nonce is advanced (folded into the RNG) and
    the matrix rebuilt SOLELY until full rank.  No condition-number threshold is
    imposed (ruling §3); the condition number is recorded only.

    Study-2 freezes M = n (square); both margins then equal k exactly."""
    if int(M) != int(n):
        raise ValueError("sparsek requires the Study-2 square geometry M == n; "
                         "got M=%d, n=%d" % (M, n))
    n = int(n)
    k = int(k)
    if not (1 <= k <= n):
        raise ValueError("sparsek occupancy k=%d out of range [1, n=%d]" % (k, n))
    kind_id = KIND_IDS["sparsek"]
    target = SPARSEK_SWITCHES_PER_N * n
    max_attempts = 200 * target + 100000        # generous anti-hang guard
    nonce = 0
    while True:
        rng = rng_for(seed, R63, PATTERN_SUBSTREAM, kind_id, k, nonce)
        perm = rng.permutation(n)                       # base permutation pi
        shifts = rng.permutation(n)[:k]                 # k DISTINCT cyclic shifts
        cols = (perm[:, None] + shifts[None, :]) % n    # (n, k) on-pixel columns
        B = np.zeros((n, n), dtype=np.uint8)
        B[np.arange(n)[:, None], cols] = 1              # row/col sums == k exactly

        edges = np.argwhere(B == 1)                     # (n*k, 2): [row, col]
        n_edges = edges.shape[0]
        accepted = attempts = 0
        while accepted < target:
            attempts += 1
            if attempts > max_attempts:
                raise RuntimeError(
                    "sparsek 2-switch did not reach %d accepted swaps in %d "
                    "attempts (k=%d n=%d)" % (target, attempts, k, n))
            e1 = int(rng.integers(0, n_edges))
            e2 = int(rng.integers(0, n_edges))
            if e1 == e2:
                continue
            i1, j1 = int(edges[e1, 0]), int(edges[e1, 1])
            i2, j2 = int(edges[e2, 0]), int(edges[e2, 1])
            if i1 == i2 or j1 == j2:
                continue
            if B[i1, j2] or B[i2, j1]:
                continue
            B[i1, j1] = 0
            B[i2, j2] = 0
            B[i1, j2] = 1
            B[i2, j1] = 1
            edges[e1, 1] = j2
            edges[e2, 1] = j1
            accepted += 1

        Bf = B.astype(np.float64)
        rank = int(np.linalg.matrix_rank(Bf))           # ruling §3: numpy rank
        if rank == n:
            cond = float(np.linalg.cond(Bf))
            return B, nonce, cond, accepted, attempts
        nonce += 1                                       # rebuild only until full


def make_patterns(kind, M, n, seed, k=None):
    """Build an illumination matrix for one (kind, M, n, seed[, k]).

    Parameters
    ----------
    kind : {"bern50", "hadpair", "gam4", "sparsek"}
    M    : number of *signed* measurements (logical rows).  For hadpair this is
           the number of complementary pairs, so A has 2*M physical rows.  For
           sparsek M must equal n (Study-2 square geometry).
    n    : number of pixels (image side^2).
    seed : integer seed folded through rng_for(seed, 63, 1, kind_id).
    k    : sparsek occupancy (on-pixels per row/column); required for sparsek,
           ignored otherwise.

    Returns
    -------
    dict with keys
      A                 : (R x n) float64, nonneg, per-pixel mean 1.  R = M for
                          bern50/gam4/sparsek, R = 2*M for hadpair.
      exposures_per_row : physical exposures charged per signed measurement
                          (1 for bern50/gam4/sparsek; 2 for hadpair pairs).
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

    elif kind == "sparsek":
        if k is None:
            raise ValueError("sparsek requires the occupancy argument k=...")
        k = int(k)
        B, nonce, cond, n_sw, n_att = _sparsek_matrix(M, n, seed, k)
        # PHOTON-BUDGET NORMALIZATION (ruling §4): A^(k) = (n/k) B^(k). With
        # every column sum exactly k and M = n, mean_i (A x)_i = (n/M) Sum_j x_j
        # = 1 for a sum-normalized x — identical mean detector load at every
        # occupancy level, on-pixel incident rate rising as n/k.
        A = (float(n) / float(k)) * B.astype(np.float64)
        exposures_per_row = 1
        occupancy = float(k) / float(n)
        b_sha = hashlib.sha256(
            np.ascontiguousarray(B, dtype=np.uint8).tobytes()).hexdigest()
        meta["construction"] = (
            "k=%d distinct cyclic shifts of a seeded base permutation + %d "
            "accepted degree-preserving 2-switches; A = (n/k) B, row/col sums "
            "== k, full numeric rank (ruling §3/§4)" % (k, n_sw))
        meta["rng_stream"] = [R63, PATTERN_SUBSTREAM, kind_id, k, nonce]
        meta["k"] = k
        meta["occupancy"] = occupancy
        meta["nonce"] = int(nonce)
        meta["cond"] = cond
        meta["row_sum"] = k
        meta["col_sum"] = k
        meta["normalization"] = "A = (n/k) * B"
        meta["B_sha256"] = b_sha
        meta["n_switch_accepted"] = int(n_sw)
        meta["n_switch_attempts"] = int(n_att)
        meta["switch_target"] = int(SPARSEK_SWITCHES_PER_N * n)
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
    all_ok = _smoke_sparsek() and all_ok
    print("[patterns smoke] %s" % ("ALL PASS" if all_ok else "FAILURES ABOVE"),
          flush=True)
    return 0 if all_ok else 1


def _smoke_sparsek(n=1024):
    """Study-2 sparse ladder checks: exact row/col sums k, full rank, exact
    unit mean load E[A x] = 1 for a sum-normalized x, and k=1 permutation."""
    print("[patterns sparsek] n=M=%d  ladder k in {512,32,16,1}" % n, flush=True)
    rng = np.random.default_rng(12345)
    x = rng.random(n)
    x = x / x.sum()                                    # sum-normalized truth
    ok = True
    for k in (512, 32, 16, 1):
        d = make_patterns("sparsek", n, n, 0, k=k)
        A, meta = d["A"], d["meta"]
        B = np.rint(A * (k / float(n))).astype(np.int64)   # recover B from A
        row_ok = bool(np.all(B.sum(axis=1) == k))
        col_ok = bool(np.all(B.sum(axis=0) == k))
        bin_ok = bool(np.all((B == 0) | (B == 1)))
        rank = int(np.linalg.matrix_rank(A.astype(np.float64)))
        rank_ok = rank == n
        emean = float(np.mean(A @ x))
        mean_ok = abs(emean - 1.0) < 1e-9
        perm_ok = True if k != 1 else (row_ok and col_ok
                                       and int(B.sum()) == n)
        this_ok = (row_ok and col_ok and bin_ok and rank_ok and mean_ok
                   and perm_ok and d["exposures_per_row"] == 1
                   and meta["k"] == k and abs(meta["occupancy"] - k / n) < 1e-15)
        ok = ok and this_ok
        print("  k=%-4d occ=%.4f rowsum=%s colsum=%s bin=%s rank=%d/%d "
              "E[Ax]=%.3e nonce=%d cond=%.3g swaps=%d/%d sha=%s.. %s"
              % (k, meta["occupancy"], row_ok, col_ok, bin_ok, rank, n,
                 emean, meta["nonce"], meta["cond"],
                 meta["n_switch_accepted"], meta["switch_target"],
                 meta["B_sha256"][:8], "PASS" if this_ok else "FAIL"),
              flush=True)
    print("[patterns sparsek] %s" % ("PASS" if ok else "FAIL"), flush=True)
    return ok


if __name__ == "__main__":
    sys.exit(_smoke())
