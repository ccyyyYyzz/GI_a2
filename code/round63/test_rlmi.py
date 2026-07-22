"""Unit tests for rlmi.py -- RLMI allocator + constrained materializer (T-C).

Covers the R25 acceptance items requested in the T-C brief:
  (a) two banks, bank 1 dominates -> w collapses to a corner;
  (b) symmetric banks -> tie-break selects the allocation closest to e_0;
  (c) KKT certificate residual < 1e-6 on a solved instance (R25.2);
  (d) materializer: exact 972 rows, guards hold, dose band verified (R25 SS9);
  (e) infeasible-guard case -> MIXTURE_MATERIALIZATION_FAIL + L0 fallback;
  (f) latency benchmark on K=8, S=16, r=200 (reports the number).
Plus:
  (g) R23 branch: loads beyond rho_c -> BRANCH_VIOLATION, U1 bound UNAVAILABLE;
  (h) kernel derivative + argmax checks.

All fixtures are SYNTHETIC and generated here (the T-C brief forbids reading the
bridge scenes).  Runnable directly (asserts + prints) or under pytest.
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rlmi


# ---------------------------------------------------------------- helpers ---- #
def _spd(rng, r, scale=1.0, ridge=1.0):
    M = rng.standard_normal((r, r))
    A = scale * (M @ M.T) + ridge * np.eye(r)
    return 0.5 * (A + A.T)


def _make_sm(H0, F, W=None, c=0.05, nu=2000.0):
    r = H0[0].shape[0]
    if W is None:
        W = np.eye(r)
    sm = rlmi.ScenarioMatrices(H0=H0, F=F, Rstar=None, vstar=None, W=W,
                               branch={}, kernel_c=c, nu=nu, r=r)
    return rlmi.compute_oracles(sm)


def _strip_bank(name, side, orient="row", amp=1.0, is_knob=False):
    """Bank of `side` full row/column strips (cyclic tiling): uniform weights,
    each pixel covered by exactly one atom -> uniform expected dose."""
    n = side * side
    rows = np.zeros((side, n))
    for i in range(side):
        for j in range(side):
            p = (i * side + j) if orient == "row" else (j * side + i)
            rows[i, p] = amp
    xhat = np.full(n, 0.5)
    return rlmi.synth_bank(name, rows, is_knob=is_knob, xhat=xhat)


# =============================================================== (h) kernel == #
def test_kernel_derivatives():
    for c in (0.0, 0.05, 0.1):
        for rho in (0.2, 1.0, 4.0, 9.0):
            h = 1e-6
            jp_fd = (rlmi.kernel_J(rho + h, c) - rlmi.kernel_J(rho - h, c)) / (2 * h)
            jpp_fd = (rlmi.kernel_Jp(rho + h, c) - rlmi.kernel_Jp(rho - h, c)) / (2 * h)
            assert abs(float(rlmi.kernel_Jp(rho, c)) - jp_fd) < 1e-5
            assert abs(float(rlmi.kernel_Jpp(rho, c)) - jpp_fd) < 1e-4
        rc = rlmi.kernel_rho_c(c)
        if np.isfinite(rc):
            assert abs(float(rlmi.kernel_Jp(rc, c))) < 1e-8
    print("[h] kernel derivative + argmax checks PASS")


# =============================================================== (a) corner === #
def test_corner_collapse():
    rng = np.random.default_rng(11)
    S, K, r = 5, 3, 6
    H0, F = [], []
    for s in range(S):
        H0.append(_spd(rng, r, scale=0.05, ridge=0.5))
        # bank 1 (library) strongly informative in all directions; 0 (knob), 2 weak
        Fk = [_spd(rng, r, scale=0.02, ridge=0.05),
              8.0 * np.eye(r) + _spd(rng, r, scale=0.1, ridge=0.0),
              _spd(rng, r, scale=0.02, ridge=0.05)]
        F.append([0.5 * (m + m.T) for m in Fk])
    sm = _make_sm(H0, F)
    mm = rlmi.standardized_maximin(sm)
    w = rlmi.tie_break_to_knob(sm, mm, knob_index=0)
    print("[a] corner: w*=%s t=%.5f kkt=%.2e" % (np.round(w, 4), mm.t_cont,
                                                 mm.kkt_residual))
    assert int(np.argmax(w)) == 1
    assert w[1] > 0.97, w
    assert mm.kkt_residual < 1e-6
    print("[a] dominant bank -> corner collapse PASS")


# =============================================================== (b) tiebreak = #
def test_symmetric_tiebreak():
    rng = np.random.default_rng(22)
    S, K, r = 4, 4, 6
    H0, F = [], []
    for s in range(S):
        H0.append(_spd(rng, r, scale=0.1, ridge=0.5))
        base = _spd(rng, r, scale=1.0, ridge=0.5)
        F.append([base.copy() for _ in range(K)])   # identical banks
    sm = _make_sm(H0, F)
    mm = rlmi.standardized_maximin(sm)
    w = rlmi.tie_break_to_knob(sm, mm, knob_index=0)
    print("[b] tiebreak: w*=%s t=%.5f (closest-to-e0)" % (np.round(w, 4), mm.t_cont))
    assert w[0] > 0.98, w                             # collapses to the knob corner
    assert np.linalg.norm(w - np.eye(K)[0]) < 2e-2
    print("[b] symmetric banks -> closest-to-e0 tie-break PASS")


# =============================================================== (c) KKT ====== #
def test_kkt_certificate():
    rng = np.random.default_rng(7)
    S, K, r = 6, 3, 6
    H0 = [_spd(rng, r, scale=1.0, ridge=2.0) for _ in range(S)]
    F = [[_spd(rng, r, scale=1.0, ridge=1.0),
          _spd(rng, r, scale=1.2, ridge=1.0),
          _spd(rng, r, scale=0.9, ridge=1.0)] for _ in range(S)]
    sm = _make_sm(H0, F)
    mm = rlmi.standardized_maximin(sm)
    # verify R25.2 directly from the returned (w*, alpha*)
    agk = (mm.alpha[:, None] * mm.g).sum(axis=0)      # sum_s alpha_s g_{s,k}
    supp = mm.support
    lam = agk[supp].mean()
    stat = np.max(np.abs(agk[supp] - lam))
    compl = np.max(mm.alpha * np.abs(mm.r_s - mm.t_cont))
    simplex = abs(mm.alpha.sum() - 1.0)
    print("[c] KKT: residual=%.2e (stat=%.1e compl=%.1e simplex=%.1e) t=%.5f"
          % (mm.kkt_residual, stat, compl, simplex, mm.t_cont))
    print("      w*=%s alpha active=%s" % (np.round(mm.w, 4),
                                           np.round(mm.alpha[mm.alpha > 1e-6], 3)))
    assert mm.kkt_residual < 1e-6
    assert stat < 1e-6 and compl < 1e-6 and simplex < 1e-9
    assert mm.t_cont >= 1.0 - 1e-9
    print("[c] KKT certificate residual < 1e-6 PASS")


# =============================================================== (d) material = #
def test_materializer_exact_972():
    side = 18                                         # 18 | 972
    knob = _strip_bank("L0", side, "row", amp=1.0, is_knob=True)
    lib = _strip_bank("L1", side, "col", amp=1.0)
    xhat = np.full(side * side, 0.5)
    w = np.array([0.4, 0.6])
    mat = rlmi.materialize([knob, lib], w, xhat, M_rows=972, dose_band=0.05,
                           incident_budget=None, peak_cap=1e9, side=side,
                           block_seed=650777)
    g = mat.guards
    print("[d] materialize: flag=%s count=%d dose_dev=%.4f (band %.2f) peak=%.2f"
          % (mat.flag, g["exact_count"], g["dose_dev"], g["dose_band"], g["peak"]))
    assert mat.flag == "OK"
    assert mat.rows.shape == (972, side * side)
    assert g["exact_count"] == 972 and g["count_ok"]
    assert g["dose_ok"] and g["dose_dev"] <= g["dose_band"]
    assert int(mat.counts.sum()) == 972
    # block order is a permutation of the exact multiset (frozen seed)
    assert sorted(mat.counts.tolist()) == sorted(np.bincount(
        np.repeat(np.arange(mat.counts.size), mat.counts),
        minlength=mat.counts.size).tolist())
    # realized bank fractions sum to 1
    assert abs(mat.what.sum() - 1.0) < 1e-9
    print("[d] exact 972 rows + guards + dose band PASS")


# ============================================================ (d2) materializer v2 #
def test_materializer_v2_power_repair():
    """R25 SS8 per-row source-power repair (rlmi._power_repair): on a design whose
    per-pixel dose sits ~7% over the band, the bounded gamma in [0.94,1.06] pulls
    it back inside.  Tested on a DECOUPLED design (each pixel covered by its own
    atom) where the repair is well-posed; asserts the trim fires, stays in bounds,
    and reaches the band.  (On the frozen T-B geometry banks the per-atom trims
    are coupled across the k pixels each atom covers, so the min-max dose cannot
    generally be driven below the ~1/k structured floor — documented in the T-D
    phase-1.5 report; see test_pure_banks_ek_real.)"""
    A = 24
    rng = np.random.default_rng(3)
    # one atom per pixel (n = A), amplitudes spread +-7% -> dose_dev ~0.07
    rows = np.zeros((A, A))
    amps = 1.0 + 0.07 * np.sign(rng.standard_normal(A))
    for a in range(A):
        rows[a, a] = amps[a]
    counts = np.full(A, 40, dtype=np.int64)                # 40 each -> 960 rows
    dose0 = (counts[:, None] * rows).sum(0)
    dev0 = float(np.abs(dose0 / dose0.mean() - 1).max())
    gamma, dev = rlmi._power_repair(counts, rows, M_rows=int(counts.sum()),
                                    dose_band=0.05, lo=0.94, hi=1.06)
    print("[d2] power-repair: dev %.4f -> %.4f  gamma[%.3f/%.3f/std%.4f]"
          % (dev0, dev, gamma.min(), gamma.max(), gamma.std()))
    assert dev0 > 0.05                                      # starts out of band
    assert dev <= 0.05 + 1e-9                               # repaired into band
    assert gamma.min() >= 0.94 - 1e-9 and gamma.max() <= 1.06 + 1e-9
    assert gamma.std() > 1e-6                               # the trim actually fired
    print("[d2] materializer-v2 bounded power repair PASS")


def test_pure_banks_ek_real():
    """DIAGNOSTIC (phase-1.5 ruling 1): materialize every pure library bank e_k,
    k=0..7, and REPORT genuine (non-fallback) vs L0-fallback status.  Reads the
    T-B library if present (bank measures, not scenes -> allowed).  Asserts the
    pipeline yields a valid exact-972 output for every bank and that L0 is genuine;
    per-bank genuine status is PRINTED, not asserted 8/8, because the current
    integer-rounding + bounded-trim materializer genuinely materializes only banks
    whose coverage granularity + coupling admit a +-5% 972-row realization (see
    the T-D phase-1.5 report: only L0, purpose-built by v5.fixed_dose_scat32, is
    genuine on the frozen T-B banks; L1-L7 fall back)."""
    libdir = os.path.join(os.path.dirname(HERE), "..", "results",
                          "round63_bridge", "library")
    libdir = os.path.normpath(os.path.join(HERE, "..", "..", "results",
                                           "round63_bridge", "library"))
    if not os.path.exists(os.path.join(libdir, "L0.npz")):
        print("[ek] library not present -> skipped (synthetic-only environment)")
        return
    banks = [rlmi.load_bank(os.path.join(libdir, "L%d.npz" % k),
                            name="L%d" % k, is_knob=(k == 0)) for k in range(8)]
    xhat = np.full(1024, 1.0 / 1024)
    genuine = []
    for k in range(8):
        mat = rlmi.materialize(banks, np.eye(8)[k], xhat, M_rows=972,
                               dose_band=0.05, incident_budget=None,
                               peak_cap=1536.0, side=32)
        gen = (mat.flag == "OK")
        genuine.append(gen)
        assert mat.rows.shape == (972, 1024)               # valid output always
        print("[ek] L%d: %-27s dose_dev=%.4f (pretrim=%.4f) genuine=%s"
              % (k, mat.flag, mat.guards.get("dose_dev", -1),
                 mat.guards.get("dose_dev_pretrim", -1), gen))
    assert genuine[0], "L0 (purpose-built dose-balanced) must materialize genuinely"
    print("[ek] pure-bank e_k report: %d/8 genuine (%s)"
          % (sum(genuine), "".join("L%d" % k for k in range(8) if genuine[k])))


# =============================================================== (e) infeasible #
def test_infeasible_guard_fallback():
    side = 12
    n = side * side
    # a bank whose atoms all pile onto the SAME few pixels -> dose cannot be
    # balanced to a tight band by any counts; knob L0 is a clean full-frame bank.
    conc = np.zeros((6, n))
    for i in range(6):
        conc[i, i] = 1.0                              # 6 single-pixel atoms
    xhat = np.full(n, 0.5)
    knob = _strip_bank("L0", side, "row", amp=1.0, is_knob=True)
    bad = rlmi.synth_bank("L1", conc, is_knob=False, xhat=xhat)
    w = np.array([0.0, 1.0])                          # force the infeasible bank
    mat = rlmi.materialize([knob, bad], w, xhat, M_rows=972, dose_band=0.05,
                           incident_budget=None, peak_cap=1e9, side=side)
    print("[e] infeasible: flag=%s rows=%s what=%s" % (mat.flag, mat.rows.shape,
                                                       np.round(mat.what, 3)))
    assert mat.flag == "MIXTURE_MATERIALIZATION_FAIL"
    assert mat.rows.shape == (972, n)                 # L0 realization returned
    assert mat.what[0] == 1.0                         # fell back to the knob corner
    print("[e] infeasible-guard -> MIXTURE_MATERIALIZATION_FAIL + L0 fallback PASS")


# =============================================================== (f) latency == #
def test_latency_benchmark():
    rng = np.random.default_rng(99)
    S, K, r = 16, 8, 200
    reps = 3
    times = []
    for _ in range(reps):
        H0 = [_spd(rng, r, scale=0.5, ridge=3.0) for _ in range(S)]
        F = [[_spd(rng, r, scale=0.3, ridge=0.4) for _ in range(K)]
             for _ in range(S)]
        t0 = time.perf_counter()
        sm = _make_sm(H0, F)
        mm = rlmi.standardized_maximin(sm)
        rlmi.tie_break_to_knob(sm, mm, knob_index=0)
        times.append(time.perf_counter() - t0)
    med = float(np.median(times))
    print("[f] LATENCY K=8 S=16 r=200: median=%.3fs  min=%.3fs  max=%.3fs "
          "(oracle+maximin+tiebreak; Gate D target <1s, ex-recon)"
          % (med, min(times), max(times)))
    # not a hard gate here -- report the number (Gate D is decided by T-D on the
    # reference CPU).  Sanity only:
    assert med < 30.0
    return med


# ==================================================== (f2) fast-linalg equiv == #
def test_fast_linalg_equivalence():
    """R27 §3.2 output-equivalence: the batched allocator (fast_linalg=True) must
    agree with the reference to |t|<=1e-8, ||w||_1<=1e-6, KKT within tol, and an
    IDENTICAL 972-row realization after the closest-to-e0 tie-break."""
    rng = np.random.default_rng(17)
    S, K, r = 16, 8, 40
    H0 = [_spd(rng, r, scale=0.4, ridge=2.5) for _ in range(S)]
    F = [[_spd(rng, r, scale=0.3, ridge=0.5) for _ in range(K)] for _ in range(S)]
    sm = _make_sm(H0, F)
    mm_r = rlmi.standardized_maximin(sm, fast_linalg=False)
    mm_f = rlmi.standardized_maximin(sm, fast_linalg=True)
    w_r = rlmi.tie_break_to_knob(sm, mm_r, knob_index=0)
    w_f = rlmi.tie_break_to_knob(sm, mm_f, knob_index=0)
    dt = abs(mm_f.t_cont - mm_r.t_cont)
    dw = float(np.abs(w_f - w_r).sum())
    # identical 972-row realization on a shared bank set built from w
    side = 18
    knob = _strip_bank("L0", side, "row", amp=1.0, is_knob=True)
    lib = _strip_bank("L1", side, "col", amp=1.0)
    xh = np.full(side * side, 0.5)
    # use a common 2-simplex projection of the (r-dim) weights for the realization
    w2r = np.array([w_r[0], 1 - w_r[0]]); w2r = np.clip(w2r, 0, None); w2r /= w2r.sum()
    w2f = np.array([w_f[0], 1 - w_f[0]]); w2f = np.clip(w2f, 0, None); w2f /= w2f.sum()
    m_r = rlmi.materialize([knob, lib], w2r, xh, M_rows=972, dose_band=0.05,
                           incident_budget=None, peak_cap=1e9, side=side)
    m_f = rlmi.materialize([knob, lib], w2f, xh, M_rows=972, dose_band=0.05,
                           incident_budget=None, peak_cap=1e9, side=side)
    ident = bool(np.array_equal(m_r.counts, m_f.counts))
    print("[f2] fast-linalg equiv: |dt|=%.2e ||dw||_1=%.2e kkt_r=%.1e kkt_f=%.1e realiz_ident=%s"
          % (dt, dw, mm_r.kkt_residual, mm_f.kkt_residual, ident))
    assert dt <= 1e-8 and dw <= 1e-6
    assert mm_f.kkt_residual <= 1e-6
    assert ident
    print("[f2] R27 §3.2 batched-allocator output-equivalence PASS")


def test_scenario_cache_identity():
    """R27 §3.2 cache: a second run_rlmi on the identical state (use_cache=True)
    reuses the cached scenario matrices and returns a BIT-IDENTICAL allocation +
    972-row realization."""
    rng = np.random.default_rng(8)
    side, n, r = 18, 18 * 18, 24
    xhat = np.clip(np.full(n, 0.5) + 0.02 * rng.standard_normal(n), 0.05, 1.0)
    Mr = rng.standard_normal((n, n)); B = np.linalg.eigh(Mr @ Mr.T)[1][:, -r:]
    Dx, Dy = rlmi._grad_ops(side)
    H0f = rlmi.tv_curvature(xhat, side, Dx, Dy) + 0.05 * np.eye(n)
    banks = [_strip_bank("L0", side, "row", amp=0.03, is_knob=True),
             _strip_bank("L1", side, "col", amp=0.03),
             _strip_bank("L2", side, "row", amp=0.05)]
    cfg = rlmi.RLMIConfig(n_scenarios=8, nu=200.0, c=0.05, eps=1e-4, side=side,
                          dose_band=0.08, peak_cap=1e9, jitter_scale=0.5,
                          use_cache=True)
    rlmi._SM_CACHE.clear()
    d1, m1_, *_ = rlmi.run_rlmi(xhat, banks, B, H0f, cfg)
    d2, m2_, *_ = rlmi.run_rlmi(xhat, banks, B, H0f, cfg)   # cache HIT
    assert np.array_equal(m1_.counts, m2_.counts)
    assert abs(d1["t_cont"] - d2["t_cont"]) == 0.0
    assert np.array_equal(np.asarray(d1["w_realized"]), np.asarray(d2["w_realized"]))
    print("[cache] scenario/oracle cache -> bit-identical re-run PASS")


# =============================================================== (g) branch === #
def test_branch_violation():
    side = 8
    n = side * side
    # a bank whose per-atom operating load a.xhat sits BEYOND rho_c(c) at high nu
    c = 0.1
    rho_c = rlmi.kernel_rho_c(c)                       # ~3.52
    xhat = np.full(n, 0.5)
    # single atom with load = 2*rho_c (well beyond the argmax)
    row = np.zeros((1, n)); row[0, :4] = (2.0 * rho_c) / (4 * 0.5)
    bank = rlmi.synth_bank("L0", row, is_knob=True, xhat=xhat)
    assert bank.rows[0] @ xhat > rho_c
    d = rlmi.r23_residual(bank, xhat, np.eye(n), nu=2000.0, c=c)
    print("[g] branch: max_load=%.3f rho_c=%.3f status=%s U1=%s"
          % (d["max_load"], d["rho_c"], d["branch_status"], d["U1_bound"]))
    assert d["branch_status"] == "BRANCH_VIOLATION"
    assert d["U1_bound"] is None                       # UNAVAILABLE, not vacuous
    # and an on-branch atom yields a finite bound
    row2 = np.zeros((1, n)); row2[0, :4] = (0.5 * rho_c) / (4 * 0.5)
    bank2 = rlmi.synth_bank("L0", row2, is_knob=True, xhat=xhat)
    d2 = rlmi.r23_residual(bank2, xhat, np.eye(n), nu=2000.0, c=c)
    assert d2["branch_status"] == "ON_BRANCH" and d2["U1_bound"] is not None
    print("[g] BRANCH_VIOLATION disclosure (bound UNAVAILABLE) PASS")


# =============================================================== end-to-end === #
def test_end_to_end_pipeline():
    """run_rlmi happy path on synthetic scenes: flag OK, exact 972, t_real>=t_cont
    (frozen-direction convention), certificate + disclosures populated."""
    rng = np.random.default_rng(5)
    side = 18
    n = side * side
    r = 24
    xhat = np.full(n, 0.5) + 0.02 * rng.standard_normal(n)
    xhat = np.clip(xhat, 0.05, 1.0)
    Mrand = rng.standard_normal((n, n))
    B = np.linalg.eigh(Mrand @ Mrand.T)[1][:, -r:]
    Dx, Dy = rlmi._grad_ops(side)
    H0_full = rlmi.tv_curvature(xhat, side, Dx, Dy) + 0.05 * np.eye(n)
    knob = _strip_bank("L0", side, "row", amp=0.03, is_knob=True)
    b1 = _strip_bank("L1", side, "col", amp=0.03)
    b2 = _strip_bank("L2", side, "row", amp=0.05)
    banks = [knob, b1, b2]
    cfg = rlmi.RLMIConfig(n_scenarios=8, nu=200.0, c=0.05, lambda_TV=1.0,
                          eps=1e-4, side=side, M_rows=972, dose_band=0.08,
                          peak_cap=1e9, jitter_scale=0.5, compute_loo=True)
    disc, mat, sm, scen, mm = rlmi.run_rlmi(xhat, banks, B, H0_full, cfg)
    print("[e2e] flag=%s w*=%s t_cont=%.4f t_real=%.4f regret=%.4f A=%.3f kkt=%.1e"
          % (disc["flag"], np.round(disc["w_star"], 3), disc["t_cont"],
             disc["t_real"], disc["materialization_regret"], disc["A_realized"],
             disc["kkt_residual_cont"]))
    assert disc["flag"] == "OK"
    assert mat.rows.shape == (972, n)
    assert disc["guards"]["exact_count"] == 972
    assert disc["t_cont"] >= 1.0 - 1e-6
    assert disc["t_real"] >= disc["t_cont"] - 1e-6      # frozen-dir consistency
    assert disc["kkt_residual_cont"] < 1e-6
    assert abs(sum(disc["w_star"]) - 1.0) < 1e-9
    assert disc["eps1_m_U1"]["status"] in ("ON_BRANCH", "BRANCH_VIOLATION")
    assert disc["delta_LOO"] is not None and disc["delta_LOO"] >= 0.0   # R25 SS3
    print("[e2e] end-to-end run_rlmi pipeline PASS (delta_LOO=%.3f)"
          % disc["delta_LOO"])


ALL = [test_kernel_derivatives, test_corner_collapse, test_symmetric_tiebreak,
       test_kkt_certificate, test_materializer_exact_972,
       test_materializer_v2_power_repair, test_pure_banks_ek_real,
       test_fast_linalg_equivalence, test_scenario_cache_identity,
       test_infeasible_guard_fallback, test_latency_benchmark,
       test_branch_violation, test_end_to_end_pipeline]


if __name__ == "__main__":
    t0 = time.perf_counter()
    med = None
    for fn in ALL:
        out = fn()
        if fn is test_latency_benchmark:
            med = out
    print("\nALL RLMI TESTS PASS  (%.1fs)  latency median K8S16r200 = %.3fs"
          % (time.perf_counter() - t0, med if med else -1))
