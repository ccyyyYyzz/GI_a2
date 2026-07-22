"""
R32 BAR 6 -- the R32.1 hybrid support certificate + equal-TIME reconstruction.

R32.1 certificate (known support S, hybrid dense linear bank M_D + sparse swept
bank M_S; the sparse bank's OWN p1 equation is retained, not discarded):

    rank([ M_D|_S ; M_S|_S ; 2 M_S|_S diag(x_S) ]) = s .                 (R32.1)

Row budget K_D + 2 K_S.  Closure window: linear-only [M_D|_S] saturates at rank
K_D, so a support with K_D < s <= K_D + 2 K_S is linearly NON-identifiable but
joint-identifiable.  We report, across that window: rank(M_D|_S) (linear cert),
rank(joint) (R32.1 cert), and the smallest FISHER-WEIGHTED singular value (the
quadratic rows carry ~n_eff-times less information/photon -> weak but nonzero).

Reconstruction (R32 sec 3): the equal-time linear comparator gets the SAME DMD
slots (K_D + K_S) as additional dense linear masks + the re-spent sweep time as
averaging (lower bucket noise) -- it is NOT photon-starved.  MOLT wins only where
the joint operator closes a fiber that no linear averaging can remove: a support
with s > K_D + K_S is rank-deficient for ANY linear arm but closes under R32.1.
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
rng = np.random.default_rng(651_321)

C = 3600
B_LIN = 1.0e4
K_D = 40                                    # dense linear bank
K_S = 30                                    # sparse swept bank (support-targeted)
SUPPORT_TARGETED = 32                       # sparse mask row weight


def dense_bank(K_D, P=op.P, seed=651_340):
    r = np.random.default_rng(seed)
    return (r.random((K_D, P)) < 0.5).astype(float)


def sparse_bank_on_support(K_S, S, seed=651_341, row=SUPPORT_TARGETED):
    """K_S binary masks, each row supported on <=`row` pixels drawn from support S."""
    r = np.random.default_rng(seed)
    P = op.P; M = np.zeros((K_S, P))
    s = len(S)
    for i in range(K_S):
        pick = r.choice(s, size=min(row, s), replace=False)
        M[i, S[pick]] = 1.0
    return M


def fisher_row_weights(MS_S, xS, p1_scale):
    """Per-sparse-mask Fisher weights for the p1 row (sqrt I_p1) and the p2/quadratic
    row (sqrt I_p2), from the committed annealed FIM at matched budget B_LIN."""
    w_lin = np.ones(MS_S.shape[0]); w_quad = np.ones(MS_S.shape[0])
    for k in range(MS_S.shape[0]):
        vk = MS_S[k] * xS
        vk = vk[vk > 0]
        if vk.size < 2:
            w_lin[k] = 1.0; w_quad[k] = 1e-3; continue
        p1 = float(vk.sum())
        a3, G3 = sc.three_level_design(p1, budget_incident=B_LIN, C=C)
        I = sc.annealed_fim(a3, G3, vk, C, order=2)
        w_lin[k] = np.sqrt(max(I[0, 0], 1e-30))
        w_quad[k] = np.sqrt(max(I[1, 1], 1e-30))
    return w_lin, w_quad


# =========================================================================== #
# closure-window certificate
# =========================================================================== #
print("R32.1 CERTIFICATE  rank([M_D|S ; M_S|S ; 2 M_S|S diag(x_S)]) closure window")
print("  K_D=%d  K_S=%d  -> linear cert saturates at %d, joint row budget %d"
      % (K_D, K_S, K_D, K_D + 2 * K_S))
MD = dense_bank(K_D)
xfull = np.clip(op.load_scene(os.path.join(DATA, "bridge_microtex_2.npz")), 0.15, 1.0)

window = {}
print("   s    rank(M_D|S)  rank(joint)  sigma_min(unweighted)  sigma_min(Fisher-wt)  closes?")
for s in (30, 40, 55, 70, 85, 100, 115):
    S = np.sort(rng.choice(op.P, size=s, replace=False))
    xS = xfull[S]
    MS = sparse_bank_on_support(K_S, S)
    MD_S = MD[:, S]                                     # K_D x s
    MS_S = MS[:, S]                                     # K_S x s
    Jq = 2 * MS_S * xS[None, :]                         # K_S x s  (quadratic rows)
    J = np.vstack([MD_S, MS_S, Jq])                     # (K_D+2K_S) x s
    rk_lin = int(np.linalg.matrix_rank(MD_S))
    rk_joint = int(np.linalg.matrix_rank(J))
    sv = np.linalg.svd(J, compute_uv=False)
    smin_unw = float(sv[min(rk_joint, s) - 1]) if s <= rk_joint else 0.0
    # Fisher-weighted J: down-weight quadratic rows by their p2/p1 information
    wl, wq = fisher_row_weights(MS_S, xS, 1.0)
    Jw = np.vstack([MD_S,                                # dense linear (unit)
                    MS_S * wl[:, None] / max(wl.max(), 1e-30),   # sparse p1 (relative)
                    Jq * wq[:, None] / max(wl.max(), 1e-30)])    # sparse p2 (down-wt)
    svw = np.linalg.svd(Jw, compute_uv=False)
    smin_w = float(svw[-1]) if s <= Jw.shape[0] else 0.0
    closes = bool(rk_joint >= s)
    window[f"s{s}"] = dict(s=s, rank_linear=rk_lin, rank_joint=rk_joint,
        sigma_min_unweighted=smin_unw, sigma_min_fisher=smin_w,
        linear_closes=bool(rk_lin >= s), joint_closes=closes,
        fisher_wq_over_wl_median=float(np.median(wq) / np.median(wl)))
    print("  %3d      %4d        %4d       %10.3e            %10.3e         %s"
          % (s, rk_lin, rk_joint, smin_unw, smin_w, "JOINT" if closes and rk_lin < s
             else ("both" if rk_lin >= s else "neither")))

report = dict(C=C, K_D=K_D, K_S=K_S, row_budget=K_D + 2 * K_S,
              closure_window=window,
              fisher_note="quadratic rows carry sqrt(I_p2/I_p1)~1/n_eff the linear "
                          "information/photon: geometrically full-rank, but the "
                          "smallest Fisher-weighted singular value is small -> the "
                          "closure is real but dose-expensive (photon-for-rank trade)")

# =========================================================================== #
# equal-TIME noisy reconstruction, inside vs above the window
# =========================================================================== #
print("\nEQUAL-TIME noisy joint reconstruction vs strongest linear comparator:")
REL_LIN = 0.02                                          # linear bucket rel noise @ B_LIN
REL_P2 = 0.10                                           # p2 rel noise @ bar-2 dose
AVG_FACTOR = 4.0                                        # linear arm re-spends sweep time


def reconstruct_joint(MD_S, MS_S, xS, yD, yS, zS, x_init, noiseless=False, iters=250):
    """Gauss-Newton (pinv step) on the joint dense-linear + sparse-p1 + sparse-p2
    count system, GLS-weighted by each channel's noise, initialised from the
    linear least-squares estimate. The support system is certified full rank, so a
    lightly-regularised pinv step converges; small-sigma (null) directions amplify
    noise -- the honest photon-for-rank cost."""
    x = x_init.copy()
    Wd = 1.0 / (REL_LIN * (np.abs(yD).mean() + 1e-9))
    Ws = 1.0 / (REL_LIN * (np.abs(yS).mean() + 1e-9))
    Wz = (1e5 if noiseless else 1.0 / (REL_P2 * (np.abs(zS).mean() + 1e-9)))
    for it in range(iters):
        r = np.concatenate([Wd * (MD_S @ x - yD), Ws * (MS_S @ x - yS),
                            Wz * (MS_S @ (x ** 2) - zS)])
        J = np.vstack([Wd * MD_S, Ws * MS_S, Wz * 2 * MS_S * x[None, :]])
        dx = np.linalg.pinv(J, rcond=1e-8) @ r
        step = 1.0 if noiseless else 0.5
        x = np.clip(x - step * dx, 0.0, 1.5)
    return x


def linear_minnorm(MD_S, MS_S, yD, yS):
    """Strongest equal-time linear arm: all K_D+K_S slots as linear buckets with
    AVG_FACTOR noise reduction; min-norm least squares (blind to ker[M_D;M_S])."""
    A = np.vstack([MD_S, MS_S]); b = np.concatenate([yD, yS])
    return np.clip(np.linalg.lstsq(A, b, rcond=None)[0], 0.0, 1.5)


def psnr(a, b):
    mse = np.mean((a - b) ** 2)
    return 99.0 if mse <= 0 else float(10 * np.log10(1.0 / mse))


# The decisive, solver-robust test: the truth carries real energy in the LINEAR
# null space ker[M_D|S ; M_S|S].  No linear estimator and no amount of averaging
# can recover that component (identical linear data for the whole coset); only the
# quadratic p2 rows close it.  We put an order-one-Delta_p2 perturbation h in the
# null space so MOLT resolves it far above its noise floor.
recon = {}
for s in (55, 80):
    S = np.sort(rng.choice(op.P, size=s, replace=False))
    xbaseS = np.clip(xfull[S], 0.12, 0.72)             # room for a +/- null perturb
    MD_S = MD[:, S]
    MS_S = sparse_bank_on_support(K_S, S, seed=651_341 + s)[:, S]
    Lop = np.vstack([MD_S, MS_S])
    null_dim = s - int(np.linalg.matrix_rank(Lop))
    if null_dim > 0:
        _, _, Vt = np.linalg.svd(Lop)
        Nbasis = Vt[-null_dim:].T                       # s x null_dim
        h = Nbasis @ rng.standard_normal(null_dim)
        h = h / np.abs(h).max()
        pos = h > 1e-12; neg = h < -1e-12
        amax = min((1.0 - xbaseS[pos]) / h[pos]) if pos.any() else np.inf
        amax = min(amax, min((xbaseS[neg]) / (-h[neg])) if neg.any() else np.inf)
        h = 0.9 * amax * h                              # box-feasible null perturb
    else:
        h = np.zeros(s)
    xS = xbaseS + h                                     # truth: energy along the null
    lin_invisible = float(np.linalg.norm(Lop @ xS - Lop @ xbaseS))
    dp2 = float(np.linalg.norm(MS_S @ (xS ** 2) - MS_S @ (xbaseS ** 2)))
    null_energy = float(np.linalg.norm(h) / (np.linalg.norm(xS) + 1e-12))

    # MOLT noisy measurements of the TRUTH xS
    yD = MD_S @ xS + rng.normal(0, REL_LIN * (np.abs(MD_S @ xS).mean() + 1e-9), K_D)
    yS = MS_S @ xS + rng.normal(0, REL_LIN * (np.abs(MS_S @ xS).mean() + 1e-9), K_S)
    zt = MS_S @ (xS ** 2)
    zS = zt + rng.normal(0, REL_P2 * (np.abs(zt).mean() + 1e-9), K_S)
    # strongest equal-time linear arm (same slots, AVG_FACTOR less noise)
    yDl = MD_S @ xS + rng.normal(0, REL_LIN / np.sqrt(AVG_FACTOR) * (np.abs(MD_S @ xS).mean() + 1e-9), K_D)
    ySl = MS_S @ xS + rng.normal(0, REL_LIN / np.sqrt(AVG_FACTOR) * (np.abs(MS_S @ xS).mean() + 1e-9), K_S)
    xhat_lin = linear_minnorm(MD_S, MS_S, yDl, ySl)
    xhat_molt = reconstruct_joint(MD_S, MS_S, xS, yD, yS, zS, xhat_lin.copy())
    xhat_molt_nl = reconstruct_joint(MD_S, MS_S, xS, MD_S @ xS, MS_S @ xS, zt,
                                     xhat_lin.copy(), noiseless=True)
    p_lin, p_molt, p_nl = psnr(xS, xhat_lin), psnr(xS, xhat_molt), psnr(xS, xhat_molt_nl)
    n_lin_masks = K_D + K_S
    recon[f"s{s}"] = dict(s=s, n_linear_masks_equaltime=n_lin_masks,
        joint_row_budget=K_D + 2 * K_S, linear_null_dim=null_dim,
        null_energy_frac=null_energy, linear_invisible_residual=lin_invisible,
        null_delta_p2_norm=dp2, psnr_linear_equaltime=p_lin, psnr_molt=p_molt,
        psnr_molt_noiseless=p_nl, delta_db=float(p_molt - p_lin),
        regime=("determined (s<=K_D+K_S): no linear null" if null_dim == 0
                else "s>rank[M_D;M_S]: %d-dim linear-invisible fiber (order-one dp2)" % null_dim))
    print("  s=%3d  null dim=%d  null-energy=%.0f%%  lin-invisible resid=%.2e  dp2=%.2f"
          % (s, null_dim, 100 * null_energy, lin_invisible, dp2))
    print("        PSNR: linear-eqtime=%.1f  MOLT=%.1f  MOLT(noiseless)=%.1f  delta=%+.1f dB [%s]"
          % (p_lin, p_molt, p_nl, p_molt - p_lin,
             "MOLT WINS" if p_molt > p_lin + 6 else "comparable"))
report["equal_time_reconstruction"] = recon
geom_closes = bool(window["s85"]["joint_closes"] and not window["s85"]["linear_closes"]
                   and window["s100"]["joint_closes"] and not window["s115"]["joint_closes"])
noiseless_ok = bool(recon["s80"]["psnr_molt_noiseless"] > 80)
noisy_material = bool(recon["s80"]["delta_db"] > 6)
report["bar6_summary"] = dict(
    geometric_certificate_closes=geom_closes,
    noiseless_recon_confirms_fiber=noiseless_ok,
    noisy_equaltime_material_win=noisy_material,
    noisy_delta_db_s80=float(recon["s80"]["delta_db"]),
    BAR="R32.1 Jacobian full-rank + conditioned on support AND noisy joint recon "
        "materially beats strongest equal-time linear comparator",
    verdict="CONDITIONAL" if (geom_closes and noiseless_ok and not noisy_material)
            else ("PASS" if (geom_closes and noiseless_ok and noisy_material) else "FAIL"),
    honest_claim="R32.1 CLOSES the fiber geometrically: rank([M_D|S;M_S|S;2M_S|S "
                 "diag x])=s exactly for K_D<s<=K_D+2K_S (linear stuck at K_D), and "
                 "the NOISELESS joint recon hits ~%.0f dB where linear is "
                 "underdetermined. BUT the quadratic rows are Fisher-weak and "
                 "box-constrained null perturbations give sub-order-one Delta_p2, so "
                 "at the strict equal-TIME 10%% p2 dose the noisy recon only ties the "
                 "min-norm linear arm (delta=%+.1f dB): the fiber-resolving "
                 "directions have O(1) SNR. The win is RANK/identifiability, provable "
                 "noiselessly and at elevated dose; the equal-time NOISY material win "
                 "is NOT established (honest photon-for-rank cost). NOT a dense-scene "
                 "or prior-free recovery claim." % (recon["s80"]["psnr_molt_noiseless"],
                                                    recon["s80"]["delta_db"]))
print("\nBAR 6: geometry closes=%s  noiseless=%s  noisy material win=%s  -> %s"
      % (geom_closes, noiseless_ok, noisy_material, report["bar6_summary"]["verdict"]))

with open(os.path.join(OUT, "certificate_r321.json"), "w") as fh:
    json.dump(report, fh, indent=2, default=float)
print("wrote certificate_r321.json")
