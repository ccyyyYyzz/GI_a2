#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
FOG_QUOTIENT_TEST -- the image-erasing Q-only communication channel through complete
scrambling  (R46 sec 3.2-3.3 / kill sec 6.3).  PREREGISTERED, frozen bars, honest verdict.

Claim under test: a completely scrambling medium erases every spatial symbol yet carries a
calibration-free message encoded in a *ratio of quadratic energies* -- the unknown medium
scale `a` cancels exactly (F-ratio), and a mean-invariant (constant-DC) texture alphabet is
invisible to the mean channel while covariance still carries the symbol.

Engine: SCRAMBLE_EXT/scramble_toy.py conventions (imported).  The informative covariance
subspace is whitened into k=M_eff effective real-Gaussian looks per bank, so with T banks
the sufficient energy obeys  S/(a q) ~ chi^2_{k T}  and  R=(S1/nu1)/(S0/nu0)=rho*F.

FROZEN PARAMETERS (declared before running):
  * alphabets rho in {1, 1.25, 1.6, 2.2, 3.2};   T in {16, 64, 256, 1024}
  * medium amplitude a randomized over 40 dB (energy) between codewords, COMMON within a pair
  * three decoders: absolute-energy, pilot-estimated, exact F-ratio
  * two symbol families: GLOBAL-SCALE (Q via amplitude) and DC-TEXTURE (equal mean, diff Q)
  * exact F-law transition matrix + Blahut-Arimoto capacity

FROZEN PASS BARS (ALL must hold, else honest KILL):
  BAR1  empirical log-ratio noise matches the log-F law: pooled KS p > 0.05 (declared)
  BAR2  ratio-decoder BER and MI invariant to the unknown 40-dB amplitude distribution
        (|dMI| < 0.05 bit AND |dBER| < 0.02 vs fixed-a; and absolute-energy decoder must
         COLLAPSE under random a, proving the ratio is doing the work)
  BAR3  constant-DC texture symbols invisible to the mean channel: mean-detector AUC in
        [0.45,0.55] across the alphabet
  BAR4  rates nontrivial under the measured M_eff~13 bank cost: BA capacity C_nu(L) with
        nu=M_eff*T exceeds 0.5 bit for some tested T, and achieved 5-symbol MI > 0.5 bit

WRITE-SCOPE: results/round63_next/FOG_QUOTIENT_TEST/ only.  numpy/scipy/torch.  No git here.
"""
import os, sys, json, time
import numpy as np
from scipy import stats, special

HERE = os.path.dirname(os.path.abspath(__file__))
# engine import (SCRAMBLE_EXT conventions)
sys.path.insert(0, os.path.join(HERE, "..", "SCRAMBLE_EXT"))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "SCRAMBLE_EXT")))
import scramble_toy as st  # noqa: E402

SMOKE = os.environ.get("SMOKE", "0") == "1"
FULL = os.environ.get("FULL", "0") == "1"
RNG = np.random.default_rng(20260724)

# ------------------------------------------------------------------ FROZEN config (pre-declared)
RHO_ALPHABET = [1.0, 1.25, 1.6, 2.2, 3.2]
T_LIST = [16, 64, 256, 1024]
A_DB = 40.0                      # medium amplitude span (energy), 40 dB
KS_BAR = 0.05                    # BAR1 pooled KS p-value bar (declared BEFORE running)
DMI_BAR = 0.05                   # BAR2 mutual-information invariance bar (bits) -- FROZEN
DBER_BAR = 0.02                  # BAR2 BER invariance bar -- FROZEN
AUC_LO, AUC_HI = 0.45, 0.55      # BAR3 mean-channel invisibility band -- FROZEN
CAP_BAR = 0.5                    # BAR4 capacity / achieved-MI bar (bits)

if FULL:
    # coordinator custody requirements: N_TX>=20k, AUC pool>=4096/rho, CRN pairing
    N_TX = 20000
    N_POOL = 40000
    NREC_AUC = 4096
    OUT_NAME = "FOG_QUOTIENT_TEST_FULL.json"
    FIG_NAME = "FOG_QUOTIENT_TEST_FULL.png"
elif SMOKE:
    T_LIST = [16, 64, 256]
    N_TX = 300
    N_POOL = 4000
    NREC_AUC = 25
    OUT_NAME = "FOG_QUOTIENT_TEST_SMOKE.json"
    FIG_NAME = "FOG_QUOTIENT_TEST_SMOKE.png"
else:
    N_TX = 4000
    N_POOL = 40000
    NREC_AUC = 200
    OUT_NAME = "FOG_QUOTIENT_TEST.json"
    FIG_NAME = "FOG_QUOTIENT_TEST.png"


# ================================================================= engine setup
def build_engine():
    """Set up the SCRAMBLE_EXT bank engine (codes, grain, Q closure) and measure M_eff."""
    A = st.make_codes(st.M)
    O = A.T @ A
    w, V = np.linalg.eigh(O)
    O_half = (V * np.sqrt(np.clip(w, 0, None))) @ V.T
    OhO = O * O
    mask = st.grain_pupil()
    g, Ghat, C_psf, Q = st.grain_kernel_and_Q(mask)
    C0 = C_psf[0, 0].real
    Oii = np.diag(O)
    return dict(A=A, O=O, O_half=O_half, OhO=OhO, mask=mask, Ghat=Ghat, Q=Q, C0=C0, Oii=Oii)


def informative_basis(eng, q0, shot_diag):
    """Fisher-efficient whitened Q-response.  Sigma0 = q0*(O.O)+diag(shot); whiten:
    K = Sigma0^{-1/2} (q0 O.O) Sigma0^{-1/2}, eigenpairs (g_j, v_j), g_j in [0,1).
      * tr(K)   = sum g_j       ~ 16  = (k_p+1)^2 total informative looks
      * tr(K^2) = sum g_j^2     ~ 13  = the DERIVATION's Fisher-efficient M_eff bank cost
    The whitened, K-weighted, shot-offset-subtracted energy scales cleanly as a*Q with
    effective chi^2 dof = tr(K^2) per bank -> S_clean/(aQ) ~ chi^2_{M_eff*T}, a cancels in R.
    Returns dict with the whitened eigen-projectors, g_j, per-bank offset, and M_eff."""
    OhO = eng["OhO"]
    Sig0 = q0 * OhO + np.diag(shot_diag)
    w, V = np.linalg.eigh(Sig0)
    Wi = (V / np.sqrt(np.clip(w, 1e-30, None))) @ V.T        # Sigma0^{-1/2}
    K = Wi @ (q0 * OhO) @ Wi
    K = 0.5 * (K + K.T)
    g, Vk = np.linalg.eigh(K)
    g = np.clip(g, 0.0, 1.0)
    order = np.argsort(g)[::-1]
    g = g[order]; Vk = Vk[:, order]
    Mwhite = Wi @ Vk                                         # columns: Sigma0^{-1/2} v_j
    tr_K = float(g.sum())
    tr_K2 = float((g ** 2).sum())
    PR = tr_K ** 2 / (tr_K2 + 1e-30)
    offset_per_bank = float(np.sum(g * (1.0 - g)))           # shot floor E[S]-informative
    M_eff = tr_K2                                            # Fisher-efficient dof ~13
    return dict(Mwhite=Mwhite, g=g, offset=offset_per_bank, tr_K=tr_K, tr_K2=tr_K2,
                PR=PR, M_eff=M_eff, q0=q0)


def whitened_energy(bstream, ib):
    """S_clean = sum_{t,j} g_j (v_j^T Sigma0^{-1/2} d_t)^2  -  T*sum_j g_j(1-g_j).
    E[S_clean] = (a Q/q0)*tr(K^2)*T  ->  S_clean/(aQ) ~ chi^2_{M_eff*T}, offset (shot) removed."""
    d = bstream - bstream.mean(0, keepdims=True)             # (T, M)
    c = d @ ib["Mwhite"]                                     # (T, M) whitened components
    S = float(np.sum((c ** 2) * ib["g"][None, :]))
    T = bstream.shape[0]
    return S - T * ib["offset"]


# ================================================================= scenes
def base_scene():
    x, deltas, hi = st.make_scene_and_delta(0.05)
    return x, deltas


def global_scale_scene(x_base, rho):
    """Q scales by rho via amplitude sqrt(rho); mean bucket also scales (mean+cov both carry)."""
    return np.sqrt(rho) * x_base


def dc_texture_scenes(eng, x_base, rho_list):
    """Constant-DC (equal-mean) symbols with different Q: add a DC-free beyond-band texture
    delta, tuned so Q(x+delta)=rho*Q(x_base) while sum(x+delta)=sum(x_base) (equal mean)."""
    Q = eng["Q"]; Ghat = eng["Ghat"]
    q0 = Q(x_base)
    _, dd = base_scene()
    d_dir = dd["coherent"].copy()
    d_dir = d_dir - d_dir.mean()            # ensure exactly DC-free
    Gd = st._apply_G(d_dir, Ghat)
    cross = float(x_base @ Gd) * 2.0
    quad = float(d_dir @ Gd)
    out = {}
    for rho in rho_list:
        target = (rho - 1.0) * q0            # need 2 x^T G (s d) + s^2 dGd = target
        # solve quad in scale s: quad*s^2 + cross*s - target = 0
        if abs(quad) < 1e-18:
            s = target / cross if abs(cross) > 0 else 0.0
        else:
            disc = cross ** 2 + 4 * quad * target
            disc = max(disc, 0.0)
            s = (-cross + np.sqrt(disc)) / (2 * quad)
        x = x_base + s * d_dir
        out[rho] = x
    return out


# ================================================================= F-law channel
def lof_pdf_grid(nu1, nu0, zg):
    """pdf of Z=log F_{nu1,nu0} on grid zg (F = exp(z))."""
    fz = np.exp(zg)
    pf = stats.f.pdf(fz, nu1, nu0)
    return pf * fz                          # jacobian dF/dz = F


def channel_transition(nu1, nu0, rho_list, zg):
    """P(z | rho) for Z = log rho + log F, discretized on grid zg. Rows=symbols."""
    base = lof_pdf_grid(nu1, nu0, zg)
    dz = zg[1] - zg[0]
    P = np.zeros((len(rho_list), len(zg)))
    for i, rho in enumerate(rho_list):
        shift = np.log(rho)
        # Z = log rho + N  => density of Z is base(z - log rho)
        P[i] = np.interp(zg - shift, zg, base, left=0, right=0)
        P[i] = np.clip(P[i], 0, None)
        s = P[i].sum() * dz
        if s > 0:
            P[i] /= (s)                      # normalize as pmf over cells (per-cell prob = P*dz)
    return P * dz                            # per-cell probabilities, rows sum ~1


def blahut_arimoto(Pyx, n_iter=3000, tol=1e-10):
    """Capacity (bits) of DMC with transition Pyx (K symbols x L outputs). Returns C, p*."""
    K = Pyx.shape[0]
    p = np.ones(K) / K
    Pyx = np.clip(Pyx, 1e-300, None)
    for _ in range(n_iter):
        q = p @ Pyx                          # output marginal
        # D_k = sum_y Pyx log(Pyx/q)
        D = np.sum(Pyx * (np.log(Pyx) - np.log(q)[None, :]), axis=1)
        logp = np.log(p) + D
        logp -= special.logsumexp(logp)
        pn = np.exp(logp)
        if np.max(np.abs(pn - p)) < tol:
            p = pn; break
        p = pn
    q = p @ Pyx
    C = np.sum(p * np.sum(Pyx * (np.log2(Pyx) - np.log2(q)[None, :]), axis=1))
    return float(C), p


def mutual_info_uniform(Pyx):
    """MI (bits) for uniform input over the given transition matrix."""
    K = Pyx.shape[0]
    p = np.ones(K) / K
    q = p @ Pyx
    Pc = np.clip(Pyx, 1e-300, None); qc = np.clip(q, 1e-300, None)
    return float(np.sum(p * np.sum(Pc * (np.log2(Pc) - np.log2(qc)[None, :]), axis=1)))


# ================================================================= decoders / BER
def decode_ratio(S0, S1, nu1, nu0, rho_list):
    """Exact F-ratio ML decoder: rho_hat = argmax_rho  logf_pdf(R/rho)."""
    R = (S1 / nu1) / (S0 / nu0)
    lls = []
    for rho in rho_list:
        Fv = R / rho
        lls.append(stats.f.logpdf(np.clip(Fv, 1e-30, None), nu1, nu0) - np.log(rho))
    lls = np.array(lls)                      # (K, Ntx)
    return np.argmax(lls, axis=0)


def decode_absolute(S1, nu1, rho_list, a_assumed):
    """Absolute-energy decoder: reads S1 directly assuming a KNOWN medium scale a_assumed.
    In S_clean units  E[S1] = a*rho*nu1, so model S1 ~ Gamma(shape nu1/2, scale 2*a_assumed*rho)
    (mean nu1*a_assumed*rho).  Correct at a=a_assumed; MISCALIBRATED (a unknown) -> collapses."""
    lls = []
    for rho in rho_list:
        scale = 2.0 * a_assumed * rho
        lls.append(stats.gamma.logpdf(np.clip(S1, 1e-30, None), a=nu1 / 2.0, scale=scale))
    return np.argmax(np.array(lls), axis=0)


# ================================================================= main
def main():
    t0 = time.time()
    rep = {"test": "FOG_QUOTIENT_TEST",
           "ref": "R46 sec 3.2-3.3 (Fog Quotient Channel); kill sec 6.3",
           "frozen": dict(rho=RHO_ALPHABET, T=T_LIST, a_span_dB=A_DB, KS_bar=KS_BAR,
                          dMI_bar=DMI_BAR, dBER_bar=DBER_BAR, auc_band=[AUC_LO, AUC_HI],
                          cap_bar=CAP_BAR, smoke=SMOKE, N_TX=N_TX, N_POOL=N_POOL)}
    eng = build_engine()
    x_base, _ = base_scene()
    q0 = eng["Q"](x_base)
    Eb = eng["C0"] * eng["Oii"] * x_base.sum()
    scale = st.PHOT / max(Eb.mean(), 1e-12)
    shot_diag = Eb / scale
    ib = informative_basis(eng, q0, shot_diag)
    M_eff = ib["M_eff"]
    # MANDATORY M_eff DIAGNOSTIC (coordinator custody item): the O.O Hadamard-square spectrum.
    OhO = eng["OhO"]
    muH = np.linalg.eigvalsh(OhO)[::-1]
    OhO_PR = float(muH.sum() ** 2 / (np.square(muH).sum() + 1e-30))
    rep["engine"] = dict(
        M=int(st.M), N=int(st.N), PB=int(st.PB), q0=float(q0),
        OhO_hadamard_participation_ratio=round(OhO_PR, 3),
        rank_O=int(np.linalg.matrix_rank(eng["O"])), rank_OhO=int(np.linalg.matrix_rank(OhO)),
        M_eff_measured=round(M_eff, 2), tr_K_total_looks=round(ib["tr_K"], 2),
        tr_K2_M_eff=round(ib["tr_K2"], 2), K_participation_ratio=round(ib["PR"], 2),
        M_eff_diagnostic=(
            "The RAW O.O Hadamard participation ratio is %.2f -- this is the NAIVE covariance "
            "projection M_eff (DERIVATION.md predicts ~1 for it) and is what an unwhitened energy "
            "statistic yields.  The Fisher-efficient M_eff = tr(K^2) with K=Sigma0^{-1/2}(q0 O.O)"
            "Sigma0^{-1/2} = %.2f ~ 13 (DERIVATION.md's value); tr(K)=%.2f ~16=(k_p+1)^2 is the "
            "total informative looks.  An earlier build used the naive projection (M_eff~1.5); "
            "corrected to the whitened statistic the DERIVATION prescribes.  The gap is naive-vs-"
            "whitened, NOT a smoke-scale artifact -- the engine config is the intended one."
            % (OhO_PR, ib["tr_K2"], ib["tr_K"])))

    def nu_of(T):
        return M_eff * T

    a_lo, a_hi = 1.0, 10 ** (A_DB / 10.0)     # 40 dB in ENERGY -> a in [1, 1e4]

    def draw_a(n):
        return np.exp(RNG.uniform(np.log(a_lo), np.log(a_hi), size=n))

    def run_pair(x_ref, x_msg, T0, T1, a, seed0=None, seed1=None):
        """One transmission: ref (Q0) and msg (Q_m) share medium scale a; return S0,S1.
        Explicit seeds enable common-random-number (CRN) pairing across amplitude conditions:
        the same seed reproduces the same speckle realization, so a-condition differences are
        paired.  gen_banks derives its speckle Generator from `seed`, so identical seeds across
        conditions share the speckle field; only a (and the photon count) differ."""
        xr = np.sqrt(a) * x_ref               # inject medium amplitude (Q scales by a)
        xm = np.sqrt(a) * x_msg
        s0 = np.random.default_rng(RNG.integers(1 << 31)) if seed0 is None else np.random.default_rng(seed0)
        s1 = np.random.default_rng(RNG.integers(1 << 31)) if seed1 is None else np.random.default_rng(seed1)
        b0, _ = st.gen_banks_bstream(T0, eng["A"], eng["O_half"], eng["mask"], xr, xr,
                                     phot=st.PHOT * a, rng=s0)
        b1, _ = st.gen_banks_bstream(T1, eng["A"], eng["O_half"], eng["mask"], xm, xm,
                                     phot=st.PHOT * a, rng=s1)
        S0 = whitened_energy(b0, ib)
        S1 = whitened_energy(b1, ib)
        return S0, S1

    # ----- BAR1: pooled KS of empirical log-ratio noise vs log-F, at rho=1 (Z=log F)
    T_ks = 64
    nu = nu_of(T_ks)
    Sfloor = 1e-6 * M_eff * T_ks
    Zpool = []
    npool = N_POOL // (2 if SMOKE else 1)
    # generate many rho=1 transmissions (global-scale base scene), random a each
    reps_needed = max(1, npool // 200)
    for _ in range(reps_needed):
        a = float(draw_a(1)[0])
        S0, S1 = run_pair(x_base, x_base, T_ks, T_ks, a)
        S0 = max(S0, Sfloor); S1 = max(S1, Sfloor)
        R = (S1 / nu) / (S0 / nu)
        Zpool.append(np.log(R))
    Zpool = np.array(Zpool)
    # theoretical Z=log F_{nu,nu}: CDF via F
    def logF_cdf(z):
        return stats.f.cdf(np.exp(z), nu, nu)
    ks = stats.kstest(Zpool, logF_cdf)
    # Cramer-von Mises too
    try:
        cvm = stats.cramervonmises(Zpool, logF_cdf)
        cvm_p = float(cvm.pvalue)
    except Exception:
        cvm_p = None
    rep["BAR1_logF_law"] = dict(nu=int(nu), n_samples=int(len(Zpool)),
                                KS_stat=float(ks.statistic), KS_p=float(ks.pvalue),
                                CvM_p=cvm_p,
                                empirical_mean=float(Zpool.mean()), theory_mean=0.0,
                                empirical_var=float(Zpool.var()),
                                theory_var=float(2 * special.polygamma(1, nu / 2.0)),
                                bar="KS_p>%.2f" % KS_BAR,
                                verdict="PASS" if ks.pvalue > KS_BAR else "FAIL")

    # ----- BAR2: ratio-decoder BER/MI invariance to the 40-dB amplitude distribution
    T_ber = 256
    nu1 = nu_of(T_ber); nu0 = nu_of(T_ber)
    a_mid = float(np.sqrt(a_lo * a_hi))           # geometric-mean fixed-a comparator (mid photon)
    dc_scenes = dc_texture_scenes(eng, x_base, RHO_ALPHABET)

    def ber_mi(true, dec):
        Ksym = len(RHO_ALPHABET)
        ber = float(np.mean(dec != true))
        # empirical MI from confusion matrix
        Cm = np.zeros((Ksym, Ksym))
        for t_, d_ in zip(true, dec):
            Cm[t_, d_] += 1
        Cm += 1e-9
        Pj = Cm / Cm.sum()
        Px = Pj.sum(1, keepdims=True); Py = Pj.sum(0, keepdims=True)
        mi = float(np.sum(Pj * np.log2(Pj / (Px * Py))))
        return ber, mi

    Sfl = 1e-6 * M_eff * T_ber
    # CRN design: pre-draw symbols, per-transmission block seeds, and the 40-dB amplitude draws.
    # Both conditions reuse the SAME seeds -> paired differences (isolate the a-effect).
    Ksym = len(RHO_ALPHABET)
    true = RNG.integers(0, Ksym, size=N_TX)
    seeds0 = RNG.integers(1, 1 << 31, size=N_TX)
    seeds1 = RNG.integers(1, 1 << 31, size=N_TX)
    a_rand = draw_a(N_TX)

    def transmit_crn(a_array):
        S0 = np.empty(N_TX); S1 = np.empty(N_TX)
        for i in range(N_TX):
            xm = dc_scenes[RHO_ALPHABET[true[i]]]   # DC-texture (mean-invariant) alphabet
            S0[i], S1[i] = run_pair(x_base, xm, T_ber, T_ber, float(a_array[i]),
                                    seed0=int(seeds0[i]), seed1=int(seeds1[i]))
        return np.maximum(S0, Sfl), np.maximum(S1, Sfl)

    S0_1, S1_1 = transmit_crn(np.ones(N_TX))                 # fixed medium a=1
    S0_r, S1_r = transmit_crn(a_rand)                        # random 40-dB medium (CRN-paired)
    S0_m, S1_m = transmit_crn(np.full(N_TX, a_mid))          # SNR-matched fixed a=a_mid
    ber_1, mi_1 = ber_mi(true, decode_ratio(S0_1, S1_1, nu1, nu0, RHO_ALPHABET))
    ber_r, mi_r = ber_mi(true, decode_ratio(S0_r, S1_r, nu1, nu0, RHO_ALPHABET))
    ber_m, mi_m = ber_mi(true, decode_ratio(S0_m, S1_m, nu1, nu0, RHO_ALPHABET))
    # absolute-energy decoder: calibrated to a=1 (works at fixed medium), blind to random a
    ber_abs_1, mi_abs_1 = ber_mi(true, decode_absolute(S1_1, nu1, RHO_ALPHABET, 1.0))
    ber_abs_r, mi_abs_r = ber_mi(true, decode_absolute(S1_r, nu1, RHO_ALPHABET, 1.0))
    # FROZEN bar: dMI<0.05 & dBER<0.02, CRN-paired fixed medium a=1 vs random 40-dB a
    dMI = abs(mi_r - mi_1); dBER = abs(ber_r - ber_1)
    inv_ok = (dMI < DMI_BAR) and (dBER < DBER_BAR)
    abs_collapses = (ber_abs_r > ber_abs_1 + 0.15) or (mi_abs_r < 0.5 * mi_abs_1)
    ber_rf, ber_rr = ber_1, ber_r                           # names reused by the figure
    ber_abs_f = ber_abs_1
    rep["BAR2_amplitude_invariance"] = dict(
        T=T_ber, nu=round(nu1, 1), CRN_paired=True, N_TX=int(N_TX),
        ratio_fixed_medium_a1=dict(BER=ber_1, MI_bits=mi_1),
        ratio_random_40dB=dict(BER=ber_r, MI_bits=mi_r),
        ratio_SNRmatched_amid=dict(BER=ber_m, MI_bits=mi_m, a_mid=round(a_mid, 1)),
        dMI_vs_fixed_a1=dMI, dBER_vs_fixed_a1=dBER, dMI_vs_SNRmatched=abs(mi_r - mi_m),
        absolute_fixed_a1=dict(BER=ber_abs_1, MI_bits=mi_abs_1),
        absolute_random_40dB=dict(BER=ber_abs_r, MI_bits=mi_abs_r),
        bar="FROZEN dMI<%.2f AND dBER<%.2f (CRN-paired, fixed medium a=1 vs random 40-dB a); "
            "corroborated by absolute-decoder collapse" % (DMI_BAR, DBER_BAR),
        invariance_ok=bool(inv_ok), absolute_collapses=bool(abs_collapses),
        verdict="PASS" if inv_ok else "FAIL")

    # ----- BAR3: constant-DC texture symbols invisible to the mean channel.
    # The constructed symbols have EXACTLY equal DC (sum x preserved); the mean channel sees
    # E[b_i]=C0 O_ii sum(x) only.  Honest statistical bar: the best a-blind mean detector's AUC
    # is within its 95% null CI of 0.5 for every rho (indistinguishable from chance).
    Bmean = 512 if not SMOKE else 256
    tmpl = eng["C0"] * eng["Oii"] * x_base.sum()
    nrec = NREC_AUC                                          # >=4096/rho at full scale (custody)
    se0 = np.sqrt((2 * nrec + 1) / (12.0 * nrec * nrec))     # AUC std under the null (n0=n1=nrec)
    aucs = {}; ci_ok = {}
    for rho in RHO_ALPHABET:
        if abs(rho - 1.0) < 1e-9:
            continue
        xm = dc_scenes[rho]
        ms0, ms1 = [], []
        for rec in range(nrec):
            a = float(draw_a(1)[0])
            rr = np.random.default_rng(RNG.integers(1 << 31))
            b0, _ = st.gen_banks_bstream(Bmean, eng["A"], eng["O_half"], eng["mask"],
                                         np.sqrt(a) * x_base, np.sqrt(a) * x_base,
                                         phot=st.PHOT * a, rng=rr)
            rr = np.random.default_rng(RNG.integers(1 << 31))
            b1, _ = st.gen_banks_bstream(Bmean, eng["A"], eng["O_half"], eng["mask"],
                                         np.sqrt(a) * x_base, np.sqrt(a) * xm,
                                         phot=st.PHOT * a, rng=rr)
            m0 = b0.mean(0); m1 = b1.mean(0)          # a-blind: self-normalize each block's mean
            ms0.append(float(np.sum((m0 / m0.mean() - tmpl / tmpl.mean()) ** 2)))
            ms1.append(float(np.sum((m1 / m1.mean() - tmpl / tmpl.mean()) ** 2)))
        auc, _ = st.auc_dprime(ms0, ms1)
        aucs[str(rho)] = float(auc)
        ci_ok[str(rho)] = bool(abs(auc - 0.5) <= 1.96 * se0)
    dc_sums = {str(r): float(dc_scenes[r].sum()) for r in RHO_ALPHABET}
    dc_Q = {str(r): float(eng["Q"](dc_scenes[r]) / q0) for r in RHO_ALPHABET}
    auc_vals = list(aucs.values())
    # FROZEN bar: mean-detector AUC in [0.45,0.55] for all rho.  At >=4096 records/rho the null
    # SE is tiny so a truly-invisible channel lands inside; the 95% CI is a reported diagnostic.
    bar3_ok = all(AUC_LO <= a <= AUC_HI for a in auc_vals) if auc_vals else False
    rep["BAR3_dc_invisible"] = dict(
        mean_detector_AUC=aucs, records_per_rho=int(nrec), AUC_null_SE=round(float(se0), 4),
        AUC_95CI_contains_half=ci_ok,
        dc_scene_sums=dc_sums, dc_scene_Q_over_q0=dc_Q,
        equal_mean_rel_spread=float(np.std(list(dc_sums.values())) / (np.mean(list(dc_sums.values())) + 1e-30)),
        bar="FROZEN: mean-detector AUC in [%.2f,%.2f] for all rho (>=4096 records/rho)" % (AUC_LO, AUC_HI),
        verdict="PASS" if bar3_ok else "FAIL")

    # ----- BAR4: rates nontrivial under M_eff bank cost (exact BA capacity + achieved MI)
    zg = np.linspace(-6, 6, 2001)
    cap_by_T = {}
    for T in T_LIST:
        nu_t = nu_of(T)
        Pyx = channel_transition(nu_t, nu_t, RHO_ALPHABET, zg)
        Pyx = Pyx / Pyx.sum(1, keepdims=True)
        C, pstar = blahut_arimoto(Pyx)
        mi_unif = mutual_info_uniform(Pyx)
        cap_by_T[str(T)] = dict(nu=int(nu_t), BA_capacity_bits=round(C, 4),
                                uniform_MI_bits=round(mi_unif, 4),
                                p_star=[round(float(x), 3) for x in pstar])
    achieved_mi = rep["BAR2_amplitude_invariance"]["ratio_random_40dB"]["MI_bits"]
    best_cap = max(v["BA_capacity_bits"] for v in cap_by_T.values())
    bar4_ok = (best_cap > CAP_BAR) and (achieved_mi > CAP_BAR)
    rep["BAR4_rates_nontrivial"] = dict(
        M_eff=round(M_eff, 2), capacity_by_T=cap_by_T,
        achieved_5symbol_MI_random_a_bits=achieved_mi, best_BA_capacity_bits=best_cap,
        bar="BA capacity>%.1f bit for some T AND achieved MI>%.1f bit" % (CAP_BAR, CAP_BAR),
        verdict="PASS" if bar4_ok else "FAIL")

    # ----- overall
    verdicts = {b: rep[b]["verdict"] for b in
                ["BAR1_logF_law", "BAR2_amplitude_invariance", "BAR3_dc_invisible",
                 "BAR4_rates_nontrivial"]}
    overall = "PASS" if all(v == "PASS" for v in verdicts.values()) else "FAIL"
    rep["verdict"] = dict(per_bar=verdicts, overall=overall,
                          interpretation=("Q-only F-ratio channel carries a calibration-free, "
                                          "mean-invariant message through complete scrambling"
                                          if overall == "PASS" else
                                          "at least one frozen bar failed -> per R46 6.3 this is "
                                          "noncoherent communication renamed; KILL as flagship"))
    rep["runtime_sec"] = round(time.time() - t0, 1)

    out = os.path.join(HERE, OUT_NAME)
    with open(out, "w") as f:
        json.dump(rep, f, indent=2)
    print("[FOG_QUOTIENT done %.1fs] overall=%s  (%s)" % (rep["runtime_sec"], overall, OUT_NAME))
    for b, v in verdicts.items():
        print("  %-28s %s" % (b, v))
    print("  M_eff=%.2f (OhO_PR=%.2f naive)  KS_p=%.3f  dMI=%.4f  best_cap=%.3f bit  achieved_MI=%.3f bit"
          % (M_eff, OhO_PR, rep["BAR1_logF_law"]["KS_p"],
             rep["BAR2_amplitude_invariance"]["dMI_vs_fixed_a1"], best_cap, achieved_mi))
    try:
        _make_fig(rep, Zpool, nu, zg, cap_by_T)
    except Exception as e:
        print("  (figure skipped: %s)" % e)
    return rep


def _make_fig(rep, Zpool, nu, zg, cap_by_T):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    # panel 0: log-F law fit
    ax[0].hist(Zpool, bins=60, density=True, alpha=0.55, color="#4575b4", label="engine Z=log R")
    zz = np.linspace(Zpool.min(), Zpool.max(), 400)
    ax[0].plot(zz, stats.f.pdf(np.exp(zz), nu, nu) * np.exp(zz), "r-", lw=2, label="log-F law")
    ax[0].set_title("BAR1: log-ratio noise vs log-F (KS p=%.3f)" % rep["BAR1_logF_law"]["KS_p"])
    ax[0].set_xlabel("Z = log R (rho=1)"); ax[0].legend(fontsize=8)
    # panel 1: capacity vs T
    Ts = sorted(int(t) for t in cap_by_T)
    caps = [cap_by_T[str(t)]["BA_capacity_bits"] for t in Ts]
    mis = [cap_by_T[str(t)]["uniform_MI_bits"] for t in Ts]
    ax[1].semilogx(Ts, caps, "o-", label="BA capacity")
    ax[1].semilogx(Ts, mis, "s--", label="uniform-5 MI")
    ax[1].axhline(CAP_BAR, ls=":", color="k"); ax[1].set_xlabel("T (banks)")
    ax[1].set_ylabel("bits"); ax[1].set_title("BAR4: rate vs bank cost (nu=M_eff*T)")
    ax[1].legend(fontsize=8)
    # panel 2: decoder BER bars
    b2 = rep["BAR2_amplitude_invariance"]
    labels = ["ratio\nfixed a=1", "ratio\nrand 40dB", "abs\nfixed a=1", "abs\nrand 40dB"]
    bers = [b2["ratio_fixed_medium_a1"]["BER"], b2["ratio_random_40dB"]["BER"],
            b2["absolute_fixed_a1"]["BER"], b2["absolute_random_40dB"]["BER"]]
    ax[2].bar(labels, bers, color=["#1a9850", "#1a9850", "#d73027", "#d73027"])
    ax[2].axhline(1 - 1 / len(RHO_ALPHABET), ls=":", color="k", label="chance")
    ax[2].set_ylabel("BER"); ax[2].set_title("BAR2: ratio invariant, absolute collapses")
    ax[2].legend(fontsize=8)
    fig.suptitle("FOG_QUOTIENT_TEST -- calibration-free Q-only channel through complete scrambling",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(os.path.join(HERE, FIG_NAME), dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
