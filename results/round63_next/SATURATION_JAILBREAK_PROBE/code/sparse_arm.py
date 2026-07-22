"""
SPARSE-MASK ARM  (missing piece 1) + DUAL RESOURCE FRAMING (missing piece 2).

The committed T2/P2 runs used 50%-density masks (n_eff 214-316) and failed the
matched-photon bars by x1.6e5 - x1.1e6 (T2) / x2.84e5 (P2 null-pair d'=3).  The
Pro dose law (appendix D2) is
        B_3sigma,knownp1 >= 25.33 * n_eff^2 / f^2   (primary opportunities/mask),
so masks with small n_eff are ~(n_eff_dense/n_eff_sparse)^2 cheaper.  Here we
rerun the SAME T2 Fisher/noise pipeline and the SAME P2 null-pair gate with
k=32-support BINARY masks (uniform-random 32-pixel supports, K=51).  Because a
32-support masked scene has at most 32 active pixels, n_eff <= 32 by construction.

We REUSE the committed framework verbatim: saturation_core (Q_exp, three_level_
design, annealed_fim, default_a_levels) and gi_operator (scenes).  crb_se_p2 and
dprime_pair are re-imported from the committed t2/p2p3 modules where possible;
the CRB helper is small and identical, reproduced here to avoid a name clash on
the module-level constants of t2_fisher_noise (which runs its report on import).

DUAL RESOURCE FRAMING (piece 2): from each sparse mask's photon budget for 10%
p2 we convert to TIME -- gates G = B/(C t*) at the per-photon ridge t*=2.821439,
C=3600, and wall time at gate rates {100 kHz, 1 MHz}.  Note C*t* = 10156 ~ B_LIN
= 1e4, so ONE linear-campaign photon budget ~ ONE gate: the photon multiple in
units of B_LIN is (to a few %) the number of gates, hence directly a microsecond
count at MHz rates.
"""
import json
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import saturation_core as sc
import gi_operator as op

OUT = os.path.dirname(HERE)
DATA = "data/r63_bridge_scenes"

C = 3600
L = 12
ORDER = 4
TARGET = 0.10
B_LIN = 1.0e4
K = op.K                                  # 51 masks (same as linear system)
SUPPORT = 32                              # k=32-support sparse binary masks
T_STAR = 2.821439372122079               # per-photon information ridge (p2)
GATE_RATES = {"100kHz": 1.0e5, "1MHz": 1.0e6}

SCENES = ["bridge_contour_0", "bridge_twopop_1", "bridge_microtex_2",
          "bridge_control_0", "bridge_contour_2", "bridge_twopop_3"]


# --------------------------------------------------------------------------- #
# sparse operator: K masks each with EXACTLY `support` random pixels = 1
# --------------------------------------------------------------------------- #
def build_sparse_operator(K=K, P=op.P, support=SUPPORT, seed=651_032):
    rng = np.random.default_rng(seed)
    M = np.zeros((K, P))
    for i in range(K):
        idx = rng.choice(P, size=support, replace=False)
        M[i, idx] = 1.0
    return M


def sparse_mask_vec(P=op.P, support=SUPPORT, seed=0):
    rng = np.random.default_rng(seed)
    m = np.zeros(P)
    m[rng.choice(P, size=support, replace=False)] = 1.0
    return m


def crb_se_p2(a_lv, v, C, M_per_level, order=ORDER):
    """CRB se(p2) with flat M per level (identical to committed t2_fisher_noise)."""
    v = np.asarray(v, float); v = v[v > 0]
    Q = sc.Q_exp(a_lv, v)
    Q2 = sc.Q_exp(2 * a_lv, v)
    sig2 = (Q - Q2) / C
    dQ = np.stack([Q * ((-1.0) ** k) * a_lv ** k / k
                   for k in range(1, order + 1)], axis=1)
    I = np.zeros((order, order))
    Ml = np.atleast_1d(M_per_level) * np.ones(len(a_lv))
    for l in range(len(a_lv)):
        I += Ml[l] * np.outer(dQ[l], dQ[l]) / sig2[l]
    cov = np.linalg.inv(I)
    return float(np.sqrt(cov[1, 1])), Q, dict(cond=float(np.linalg.cond(I)))


# =========================================================================== #
# PART A -- sparse-mask T2 (Fisher effect size) + n_eff^2 dense/sparse law
# =========================================================================== #
report = {"C": C, "L": L, "order": ORDER, "target_rel_se": TARGET, "B_LIN": B_LIN,
          "support": SUPPORT, "K": K, "t_star": T_STAR, "scenes": {}}
print("SPARSE-MASK ARM  (32-support binary masks)")
print("scene                n_eff   p1      p2     relse@M1   B_sat(10%)   x_linear")
mults = []
for si, name in enumerate(SCENES):
    x = op.load_scene(os.path.join(DATA, name + ".npz"))
    m = sparse_mask_vec(seed=651_600 + si)          # one 32-support mask/scene
    v = m * x
    vp = v[v > 0]
    p = sc.power_sums(vp, 4); p1, p2 = float(p[1]), float(p[2])
    n_eff = p1 ** 2 / p2
    a_lv = sc.default_a_levels(vp, L=L)
    se1, Q, dico = crb_se_p2(a_lv, vp, C, 1.0)
    rel1 = se1 / p2
    M_need = (rel1 / TARGET) ** 2
    sum_occ = float(np.sum(a_lv * p1))
    B_sat = M_need * C * sum_occ
    mult = B_sat / B_LIN
    mults.append(mult)
    report["scenes"][name] = dict(support=SUPPORT, n_active=int(vp.size),
        n_eff=float(n_eff), p1=p1, p2=p2, relse_at_M1=float(rel1),
        M_need=float(M_need), B_sat=float(B_sat), budget_multiple=float(mult),
        fisher_cond=dico["cond"])
    print("%-20s %6.1f %7.1f %7.1f  %8.2e  %9.2e  %8.2e"
          % (name, n_eff, p1, p2, rel1, B_sat, mult))

# --- dense-vs-sparse n_eff^2 law: same scenes, 50% vs 32-support -------------
print("\nn_eff^2 SCALING LAW  (dense 50%% vs sparse 32-support, same scene):")
law = {}
for si, name in enumerate(SCENES):
    x = op.load_scene(os.path.join(DATA, name + ".npz"))
    row = {}
    for tag, m in (("dense", (np.random.default_rng(651_400).random(op.P) < 0.5).astype(float)),
                   ("sparse", sparse_mask_vec(seed=651_600 + si))):
        v = (m * x); vp = v[v > 0]
        p = sc.power_sums(vp, 4); p1, p2 = float(p[1]), float(p[2])
        n_eff = p1 ** 2 / p2
        a_lv = sc.default_a_levels(vp, L=L)
        se1, _, _ = crb_se_p2(a_lv, vp, C, 1.0)
        B_sat = (se1 / p2 / TARGET) ** 2 * C * float(np.sum(a_lv * p1))
        row[tag] = dict(n_eff=float(n_eff), B_sat=float(B_sat))
    r_neff = row["dense"]["n_eff"] / row["sparse"]["n_eff"]
    r_B = row["dense"]["B_sat"] / row["sparse"]["B_sat"]
    row["neff_ratio"] = float(r_neff)
    row["B_ratio"] = float(r_B)
    row["B_ratio_over_neff_ratio_sq"] = float(r_B / r_neff ** 2)
    law[name] = row
    print("  %-20s n_eff %6.1f->%5.1f (x%.2f)  B_sat ratio %.1f  vs (n_eff ratio)^2 %.1f"
          "  -> ratio-of-ratios %.3f"
          % (name, row["dense"]["n_eff"], row["sparse"]["n_eff"], r_neff,
             r_B, r_neff ** 2, row["B_ratio_over_neff_ratio_sq"]))
report["neff2_dense_vs_sparse"] = law
rr = [law[n]["B_ratio_over_neff_ratio_sq"] for n in law]
report["neff2_law_ratioofratios_range"] = [float(min(rr)), float(max(rr))]
print("  B_ratio/(n_eff ratio)^2 in [%.3f, %.3f]  (==1 confirms B_sat ~ n_eff^2)"
      % (min(rr), max(rr)))

report["summary_T2"] = dict(
    min_budget_multiple=float(np.min(mults)),
    max_budget_multiple=float(np.max(mults)),
    median_budget_multiple=float(np.median(mults)),
    BAR="rel se(p2) <= 10% at per-mask budget <= 4x linear arm",
    PASS=bool(np.min(mults) <= 4.0),
    note="sparse 32-support cuts n_eff ~10x -> ~100x cheaper than dense, still "
         ">>4x on matched PHOTONS (see dual-resource framing for the TIME view)")
print("\nsparse-mask T2 budget multiple: min=%.2e  median=%.2e  max=%.2e"
      % (np.min(mults), np.median(mults), np.max(mults)))


# =========================================================================== #
# PART B -- sparse-mask P2 null-pair gate (the decisive economics test)
# =========================================================================== #
print("\n=== GATE P2 (sparse 32-support masks): null-pair bank d' ===")
Msp = build_sparse_operator()
rankMsp = int(np.linalg.matrix_rank(Msp))
U, S, Vt = np.linalg.svd(Msp, full_matrices=True)
kerM = Vt[rankMsp:].T
print("  rank(M_sparse)=%d  dim ker=%d" % (rankMsp, kerM.shape[1]))

LO, HI = 0.08, 0.92
xbase = np.clip(op.load_scene(os.path.join(DATA, "bridge_microtex_0.npz")), LO, HI)


def valid_pair(h_dir):
    h = h_dir / np.abs(h_dir).max()
    pos = h > 1e-12; neg = h < -1e-12
    amax = np.inf
    if pos.any(): amax = min(amax, np.min((1 - xbase[pos]) / h[pos]))
    if neg.any(): amax = min(amax, np.min((0.0 - xbase[neg]) / h[neg]))
    return xbase, xbase + 0.9 * amax * h, 0.9 * amax * h


def dprime_pair(x, xp, budget_per_mask=B_LIN):
    d2 = 0.0
    neffs = []
    for k in range(K):
        m = Msp[k]
        vx = (m * x); vxa = vx[vx > 0]
        p1 = float(vxa.sum())
        if p1 <= 0:
            continue
        p2 = float((vxa ** 2).sum())
        neffs.append(p1 ** 2 / p2)
        dp2 = float(np.sum(m * (xp ** 2 - x ** 2)))
        a_lv, G = sc.three_level_design(p1, budget_incident=budget_per_mask, C=C)
        I = sc.annealed_fim(a_lv, G, vxa, C, order=2)
        I_p2 = I[1, 1] - I[1, 0] * I[0, 0] ** -1 * I[0, 1]
        d2 += dp2 ** 2 * I_p2
    return float(np.sqrt(max(d2, 0.0))), float(np.mean(neffs)), float(np.max(neffs))

rng = np.random.default_rng(651_919)
NBANK = 40
dprimes = []; neff_means = []; neff_maxes = []
for b in range(NBANK):
    h_dir = kerM @ rng.standard_normal(kerM.shape[1])
    x, xp, h = valid_pair(h_dir)
    dpr, nm, nmax = dprime_pair(x, xp)
    dprimes.append(dpr); neff_means.append(nm); neff_maxes.append(nmax)
dprimes = np.array(dprimes)
median_d = float(np.median(dprimes))
frac2 = float(np.mean(dprimes >= 2)); frac3 = float(np.mean(dprimes >= 3))
budget_mult_for_d3 = (3.0 / median_d) ** 2
print("  bank %d  median d'=%.4f  frac>=2=%.0f%%  frac>=3=%.0f%%  (range %.4f..%.4f)"
      % (NBANK, median_d, 100 * frac2, 100 * frac3, dprimes.min(), dprimes.max()))
print("  per-mask n_eff: mean-of-means=%.1f  max=%.1f  (<=%d by construction)"
      % (np.mean(neff_means), np.max(neff_maxes), SUPPORT))
print("  budget multiple for median d'=3: %.1fx the matched linear campaign"
      % budget_mult_for_d3)

# d' table: representative percentiles + gamma-k degradation
dtab = {"p10": float(np.percentile(dprimes, 10)),
        "median": median_d,
        "p90": float(np.percentile(dprimes, 90)),
        "max": float(dprimes.max())}
gk = {f"k{k}": dict(median_dprime=float(median_d / k),
                    budget_mult_for_d3=float((3.0 / (median_d / k)) ** 2))
      for k in (3, 10)}
report["P2_sparse"] = dict(support=SUPPORT, K=K, rank_M=rankMsp,
    linear_campaign=float(K * B_LIN), bank_size=NBANK,
    median_dprime=median_d, frac_ge2=frac2, frac_ge3=frac3,
    dprime_min=float(dprimes.min()), dprime_max=float(dprimes.max()),
    dprime_table=dtab, per_mask_neff_mean=float(np.mean(neff_means)),
    per_mask_neff_max=float(np.max(neff_maxes)),
    budget_mult_for_median_d3=float(budget_mult_for_d3),
    gamma_k_degradation=gk,
    BAR="median d'>=3 AND >=80% pairs >=2 at matched budget",
    PASS=bool(median_d >= 3 and frac2 >= 0.8))

# reference dense P2 numbers (from committed p2p3_gates.json) for the ratio story
DENSE_P2_MEDIAN_D = 0.005629944340455816
DENSE_P2_BUDMULT = 283945.053144013
report["P2_sparse"]["dense_reference"] = dict(
    dense_median_dprime=DENSE_P2_MEDIAN_D, dense_budget_mult_for_d3=DENSE_P2_BUDMULT,
    sparse_over_dense_dprime=float(median_d / DENSE_P2_MEDIAN_D),
    dense_over_sparse_budmult=float(DENSE_P2_BUDMULT / budget_mult_for_d3))
print("  vs committed DENSE P2: median d' %.4f->%.4f (x%.1f);  budget-mult %.2e->%.2e (x%.1f cheaper)"
      % (DENSE_P2_MEDIAN_D, median_d, median_d / DENSE_P2_MEDIAN_D,
         DENSE_P2_BUDMULT, budget_mult_for_d3, DENSE_P2_BUDMULT / budget_mult_for_d3))


# =========================================================================== #
# PART C -- DUAL RESOURCE FRAMING (piece 2): photons vs TIME, per sparse mask
# =========================================================================== #
print("\n=== DUAL RESOURCE FRAMING (photons vs time), sparse masks ===")
print("  C*t* = %.0f  (~B_LIN=%.0f -> photon-multiple ~ #gates)" % (C * T_STAR, B_LIN))
dual = {"C_times_tstar": float(C * T_STAR), "t_star": T_STAR, "B_LIN": B_LIN,
        "gate_rates_Hz": GATE_RATES, "per_scene": {}}
print("  %-20s  photonx   B_10        gates      t@100kHz   t@1MHz"
      % "scene")
for name in SCENES:
    B10 = report["scenes"][name]["B_sat"]
    pmult = B10 / B_LIN
    G = B10 / (C * T_STAR)
    walls = {r: G / f for r, f in GATE_RATES.items()}
    dual["per_scene"][name] = dict(B_10=float(B10), photon_multiple=float(pmult),
        gates_equiv=float(G),
        wall_s_100kHz=float(walls["100kHz"]), wall_s_1MHz=float(walls["1MHz"]))
    print("  %-20s  %.1e   %.2e   %.2e   %.2e s  %.2e s"
          % (name, pmult, B10, G, walls["100kHz"], walls["1MHz"]))

# explicit x30 reference point the strategy ruling asked for
for xmult in (30.0, 100.0):
    B = xmult * B_LIN
    G = B / (C * T_STAR)
    dual[f"reference_x{int(xmult)}"] = dict(photon_multiple=xmult, gates_equiv=float(G),
        wall_s_100kHz=float(G / 1e5), wall_s_1MHz=float(G / 1e6))
    print("  REFERENCE  x%3d photons -> %.1f gates -> %.1f us @1MHz / %.1f us @100kHz"
          % (xmult, G, 1e6 * G / 1e6, 1e6 * G / 1e5))

# also the ACTUAL decisive (P2 null-pair d'=3) budget in the time frame
Bd3 = budget_mult_for_d3 * (K * B_LIN)            # total over the K-mask campaign
Bd3_permask = budget_mult_for_d3 * B_LIN          # per mask
G_d3 = Bd3_permask / (C * T_STAR)
dual["decisive_P2_d3"] = dict(photon_multiple_vs_linear_campaign=float(budget_mult_for_d3),
    per_mask_B=float(Bd3_permask), per_mask_gates=float(G_d3),
    per_mask_wall_s_100kHz=float(G_d3 / 1e5), per_mask_wall_s_1MHz=float(G_d3 / 1e6))
print("  DECISIVE P2 d'=3 (matched-photon x%.0f): %.1f gates/mask -> %.2e s @1MHz, %.2e s @100kHz"
      % (budget_mult_for_d3, G_d3, G_d3 / 1e6, G_d3 / 1e5))
dual["dose_caveat"] = (
    "The two frames diverge by design. On the PHOTON axis MOLT's second channel "
    "costs ~n_eff^2 more primary opportunities than the linear bucket, so for "
    "DOSE-LIMITED scenes (photodamage, single-molecule, quantum-light budgets) the "
    "photon frame governs and MOLT LOSES. On the TIME axis those same photons are "
    "delivered in G=B/(C t*) gates at 0.1-1 MHz, so even the decisive matched-photon "
    "multiple lands at milliseconds/mask and a x30 multiple at microseconds/mask; for "
    "LIGHT-CHEAP, robust scenes the time frame governs and MOLT is essentially free.")
report["dual_resource"] = dual
print("  DOSE CAVEAT:", dual["dose_caveat"][:96], "...")

with open(os.path.join(OUT, "sparse_arm.json"), "w") as fh:
    json.dump(report, fh, indent=2, default=float)
print("\nwrote sparse_arm.json")
