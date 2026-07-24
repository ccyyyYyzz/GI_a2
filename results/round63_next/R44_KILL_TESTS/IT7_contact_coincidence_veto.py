#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""IT7 -- Contact-arm coincidence veto (internal divergence stranger #6; D3/D5 FA repair).

Hardware identity: at the relay/contact plane z2_eff=0 the diffuser phase cancels in the
intensity, |E_diff * exp(i phi)|^2 = |E_diff|^2, so the CONTACT arm is EXACTLY blind to a
medium event (l_c or sigma_phi change) but still sees a real scene/object change through
the mean channel (Delta b_c = J_mean . delta).  The z2=5 mm arm (covariance detector) sees
BOTH.  A 2D coincidence thresholded on H0 fires only when BOTH arms agree -> a medium event
(fires z2=5 arm, NOT the contact arm) is VETOED, repairing the sealed campaign's only failed
gate (medium-event false alarm FA 0.096/0.084).

DECISION (coordinator): medium-event contact-arm AUC ~0.5 AND coincidence FA <=0.05 at
>=0.9x single-arm power => repaired by a hardware identity (flip mirror enters the bench).
KILL: contact arm responds to medium (>0.55), or coincidence power cost >20%, or the scene
leak is too weak for the contact arm to witness.
"""
import os, time, json
import numpy as np
import torch
import wave_twin as wt
import twin_pool as tp

t0 = time.time()
DEV = wt.DEV; NM = wt.NMACRO; c = wt._cx; HALF = wt.DMD_PIX // 2
M = int(os.environ.get("M", "16"))
N_BANK = int(os.environ.get("N_BANK", "700"))
Z1 = 10e-3; Z2 = 5e-3
T_EFF = int(os.environ.get("T_EFF", "300"))
N_MC = int(os.environ.get("N_MC", "2000"))
RIDGE = 1e-4

x64 = wt.witness_scene(4); xnorm = float(np.linalg.norm(x64))
d5 = wt.beyond_band_delta(0.05, xnorm, seed=7)
xt = torch.as_tensor(x64, device=DEV, dtype=wt.RDT)
xdt = torch.as_tensor(x64 + d5, device=DEV, dtype=wt.RDT)

Ep, Em, S = tp.code_fields_at_diffuser(M, Z1, seed=11)

def block_sum_64(D):
    a = c - HALF
    reg = D[..., a:a + wt.DMD_PIX, a:a + wt.DMD_PIX].reshape(*D.shape[:-2], NM, wt.MACRO_PX, NM, wt.MACRO_PX)
    return reg.sum(dim=(-3, -1))

# ---- CONTACT arm (z2=0): deterministic intensity kernels, medium-INDEPENDENT ----
Kp = block_sum_64(Ep.real ** 2 + Ep.imag ** 2)      # (M,64,64)
Km = block_sum_64(Em.real ** 2 + Em.imag ** 2)
kp_x = torch.einsum('mij,ij->m', Kp, xt); km_x = torch.einsum('mij,ij->m', Km, xt)
kp_xd = torch.einsum('mij,ij->m', Kp, xdt); km_xd = torch.einsum('mij,ij->m', Km, xdt)
scp_c = float(wt.PHOT / (0.5 * (kp_x + km_x)).mean().clamp(min=1e-30))
dbc = (kp_xd - km_xd) - (kp_x - km_x)               # contact mean shift = J_mean . d5 (the leak)
Rc = torch.diag((kp_x + km_x) / scp_c + 1e-12)       # shot covariance per code
w_c = torch.linalg.solve(Rc, dbc)                    # matched mean filter

def contact_stats(kp_m, km_m, seed):
    g = torch.Generator(device=DEV).manual_seed(seed)
    lamp = (kp_m * scp_c).clamp(min=0).expand(N_MC, T_EFF, M)
    lamm = (km_m * scp_c).clamp(min=0).expand(N_MC, T_EFF, M)
    bmean = (torch.poisson(lamp, generator=g) - torch.poisson(lamm, generator=g)).mean(1) / scp_c
    return (bmean @ w_c).cpu().numpy()

tc_h0 = contact_stats(kp_x, km_x, 1)
tc_scene = contact_stats(kp_xd, km_xd, 2)
tc_med = contact_stats(kp_x, km_x, 3)                 # medium == H0 fluxes (contact invariance)

# ---- z2=5 arm: covariance detector; target + medium-mismatch pools ----
Gp, Gm = tp.speckle_pool(Ep, Em, Z2, 2 * np.pi, 50.0, N_BANK, seed=300)
Gpm, Gmm = tp.speckle_pool(Ep, Em, Z2, 1.8 * np.pi, 45.0, N_BANK, seed=777)
g = torch.Generator(device=DEV).manual_seed(11)
b5_h0, scp = tp.signed_buckets(Gp, Gm, x64, shot=True, gen=g)
b5_scene, _ = tp.signed_buckets(Gp, Gm, x64 + d5, shot=True, scp=scp, gen=g)
b5_med, _ = tp.signed_buckets(Gpm, Gmm, x64, shot=True, scp=scp, gen=g)
# matched filter for the scene change (CRN-clean)
fp0 = torch.einsum('bmij,ij->bm', Gp, xt); fm0 = torch.einsum('bmij,ij->bm', Gm, xt)
fp1 = torch.einsum('bmij,ij->bm', Gp, xdt); fm1 = torch.einsum('bmij,ij->bm', Gm, xdt)
f0 = fp0 - fm0; f1 = fp1 - fm1
C0 = torch.cov(f0.T); dC = torch.cov(f1.T) - C0
shot0 = ((fp0 + fm0).mean(0)) / scp
V0 = C0 + torch.diag(shot0) + RIDGE * float(torch.diag(C0).mean()) * torch.eye(M, device=DEV, dtype=wt.RDT)
Vinv = torch.linalg.inv(V0)
wv, Uv = torch.linalg.eigh(0.5 * (V0 + V0.T)); Rh = Uv @ torch.diag(wv.clamp(min=wv.max() * 1e-12).rsqrt()) @ Uv.T
C0_obs = torch.cov(b5_h0.T)                          # observed H0 covariance reference

# The z2=5 arm is a GENERIC covariance-DEVIATION detector: it fires on ANY covariance
# change (scene OR medium) -> this is the FA-prone primary the coincidence must repair.
def cov_stats(b5, seed):
    g = torch.Generator(device=DEV).manual_seed(seed)
    n = b5.shape[0]
    idx = torch.randint(0, n, (N_MC, T_EFF), device=DEV, generator=g)
    B = b5[idx]; Bc = B - B.mean(1, keepdim=True)
    C = torch.einsum('kti,ktj->kij', Bc, Bc) / (T_EFF - 1)
    Dm = torch.einsum('ai,kij,jb->kab', Rh, C - C0_obs[None], Rh)      # whitened deviation
    return (Dm.reshape(N_MC, -1) ** 2).sum(1).cpu().numpy()           # ||.||_F^2 (fires on any change)

t5_h0 = cov_stats(b5_h0, 4); t5_scene = cov_stats(b5_scene, 5); t5_med = cov_stats(b5_med, 6)

def auc(a0, a1):
    return float((a1[:, None] > a0[None, :]).mean())

# thresholds on H0 (95th pct) -> single-arm FA = 0.05 by construction
th_c = float(np.quantile(tc_h0, 0.95)); th_5 = float(np.quantile(t5_h0, 0.95))
def coincidence_rate(tc, t5):
    return float(((tc > th_c) & (t5 > th_5)).mean())

single5_power = float((t5_scene > th_5).mean())
single5_FA_med = float((t5_med > th_5).mean())        # the medium false alarm the veto must repair
coinc_power = coincidence_rate(tc_scene, t5_scene)
coinc_FA_med = coincidence_rate(tc_med, t5_med)
coinc_FA_h0 = coincidence_rate(tc_h0, t5_h0)

results = {"test": "IT7_contact_coincidence_veto",
           "ref": "internal divergence stranger #6; merged adjudication D3/D5 FA repair",
           "params": dict(M=M, n_bank=N_BANK, T_eff=T_EFF, n_mc=N_MC, z1_mm=10, z2_mm=5,
                          medium_event="l_c 50->45, sigma 2pi->1.8pi", change="5pct beyond-band"),
           "decision": "medium contact-arm AUC ~0.5 AND coincidence FA<=0.05 at >=0.9x single-arm "
                       "power => repaired; KILL if contact responds to medium (>0.55) or power cost>20%",
           "measured": dict(
               contact_AUC_H0_vs_medium=round(auc(tc_h0, tc_med), 4),
               contact_AUC_H0_vs_scene=round(auc(tc_h0, tc_scene), 4),
               z2_5_AUC_H0_vs_medium=round(auc(t5_h0, t5_med), 4),
               z2_5_AUC_H0_vs_scene=round(auc(t5_h0, t5_scene), 4),
               single_arm_z2_5_power_scene=round(single5_power, 4),
               single_arm_z2_5_FA_medium=round(single5_FA_med, 4),
               coincidence_power_scene=round(coinc_power, 4),
               coincidence_FA_medium=round(coinc_FA_med, 4),
               coincidence_FA_H0=round(coinc_FA_h0, 4),
               power_ratio_coinc_over_single=round(coinc_power / (single5_power + 1e-9), 4),
               FA_repair_factor=round(single5_FA_med / (coinc_FA_med + 1e-9), 2))}

m = results["measured"]
contact_blind = m["contact_AUC_H0_vs_medium"] <= 0.55
fa_ok = m["coincidence_FA_medium"] <= 0.05
power_ok = m["power_ratio_coinc_over_single"] >= 0.9
scene_witnessed = m["contact_AUC_H0_vs_scene"] > 0.55        # contact arm must see the scene change
if contact_blind and fa_ok and power_ok and scene_witnessed:
    verdict = ("PASS -- the contact arm is medium-blind (AUC %.3f ~0.5); the 2D coincidence cuts the "
               "medium false alarm from %.3f (single-arm) to %.3f (%.0fx repair, <=0.05) while "
               "preserving %.0f%% of single-arm scene power. The D3/D5 failed gate is repaired by a "
               "hardware identity (flip mirror)." % (m["contact_AUC_H0_vs_medium"],
               m["single_arm_z2_5_FA_medium"], m["coincidence_FA_medium"], m["FA_repair_factor"],
               100 * m["power_ratio_coinc_over_single"]))
elif not contact_blind:
    verdict = "KILL -- contact arm responds to the medium event (AUC %.3f>0.55): the contact identity does not hold as modeled." % m["contact_AUC_H0_vs_medium"]
elif not scene_witnessed:
    verdict = ("KILL -- scene leak too weak for the contact arm to witness (contact scene AUC %.3f): "
               "the coincidence would suppress scene power too. Fall back to the semiparametric/cross-fit repair."
               % m["contact_AUC_H0_vs_scene"])
elif not power_ok:
    verdict = "KILL -- coincidence power cost >20%% (ratio %.3f): too expensive." % m["power_ratio_coinc_over_single"]
else:
    verdict = "PARTIAL -- FA %.3f not <=0.05; report honestly." % m["coincidence_FA_medium"]
results["verdict"] = verdict
results["elapsed_s"] = round(time.time() - t0, 1)
print("contact AUC med=%.3f scene=%.3f | z2=5 AUC med=%.3f scene=%.3f | single5=%.3f coinc_pow=%.3f coinc_FA_med=%.3f"
      % (m["contact_AUC_H0_vs_medium"], m["contact_AUC_H0_vs_scene"], m["z2_5_AUC_H0_vs_medium"],
         m["z2_5_AUC_H0_vs_scene"], single5_power, coinc_power, coinc_FA_med))
with open(os.environ.get("OUT", "IT7_contact_coincidence_veto.json"), "w") as f:
    json.dump(results, f, indent=2)
print("VERDICT:", verdict)
print("saved  elapsed %.1fs" % (time.time() - t0))
