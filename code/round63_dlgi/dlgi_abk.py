#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
DLGI CONFIRMATORY CAMPAIGN -- preregistered (A,B,K) predictions (R34 sec.2.4).

BEFORE any confirmatory data exist, compute for every schedule and cell the exact
marginal-hyperparameter Fisher blocks and the canonical-confusion geometry
(R33 Thm 2), and FREEZE the predicted schedule ordering.  This file is committed
before the confirmatory bank is opened; the empirical ordering may NOT be selected
after seeing results (R34 bar C6).

    A = I_xx (scene oracle ledger),  B = I_thetatheta (medium oracle ledger),
    C = I_xtheta (confusion),        K = A^{-1/2} C B^{-1/2},
    J_x     = A^{1/2}(I - K K^T)A^{1/2},   J_theta = B^{1/2}(I - K^T K)B^{1/2}.

Reciprocity (bar C6):  det(J_x)/det(A) = det(J_theta)/det(B) = prod_i (1 - kappa_i^2).

Machinery reused VERBATIM from the feasibility probe's dl_t3_reciprocity
(marginal_fisher, ledger_geometry, compact_instance, schedule_slots,
verify_theorem2_generic).  The exact 1024-pixel marginal Fisher is intractable
(2048x2048 Sigma with per-pixel derivatives); the reciprocity theorem is exact on
any instance, so it is evaluated on the frozen compact GI witness (K=16 pixels,
n_pair=24), identical to the probe.  The photon-SNR axis scales Sigma and leaves the
reciprocity structure invariant, so the frozen predictions are stated per (t_c, CV)
schedule and hold across the SNR grid.

PREDICTED ORDERING (frozen).  The full-record Fisher information (A,B,K and the two
efficiencies) is ~schedule-invariant: a schedule is a permutation of the identical
pattern multiset, so it preserves the joint information (small ||K|| spread).  The
paired schedule is the PRIMARY acquisition and is predicted best for BOTH ledgers at
the ESTIMATOR level, because the complementary-pair differential is well-conditioned
only when a pair's members are temporally adjacent; random/ordered scramble that and
degrade the differential SCENE reconstruction (a data-processing loss, R33 sec.3),
while the medium estimate stays schedule-robust.  There is NO forced scene-vs-medium
trade.

Writes results/round63_dlgi_campaign/predictions/ABK_PREDICTIONS.json.  CPU.
"""
import os
import sys
import json
import platform
from datetime import datetime, timezone

import numpy as np
from numpy.linalg import inv, slogdet

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dlgi_common as K
# frozen reciprocity machinery (verbatim import from the probe):
sys.path.insert(0, os.path.join(K.REPO, "results", "round63_next", "DUAL_LEDGER_PROBE", "code"))
import dl_t3_reciprocity as T3

OUT = os.path.join(K.CAMPAIGN, "predictions")
COMPACT = dict(Kpix=16, n_pair=24)                      # frozen compact witness (== probe)


def compute():
    M, x = T3.compact_instance(K=COMPACT["Kpix"], n_pair=COMPACT["n_pair"])
    n_pair = COMPACT["n_pair"]

    generic = T3.verify_theorem2_generic()              # frozen identity check

    per_cell = {}
    reciprocity_worst = 0.0
    for tc in K.TC_PRIMARY:
        for cv in K.CV_PRIMARY:
            tag = K.cell_key(tc, cv, 1.0).rsplit("_ph", 1)[0]    # tc..cv.. (SNR-invariant)
            geo_by_sched = {}
            for sch in ["paired", "random", "split"]:            # probe schedule names
                A, B, Cc = T3.marginal_fisher(M, x, T3.schedule_slots(sch, n_pair), tc, cv)
                geo = T3.ledger_geometry(A, B, Cc)
                # reciprocity residual on THIS GI Fisher (bar C6)
                effx = float(np.exp(slogdet(A - Cc @ inv(B) @ Cc.T)[1] - slogdet(A)[1]))
                efft = float(np.exp(slogdet(B - Cc.T @ inv(A) @ Cc)[1] - slogdet(B)[1]))
                geo["reciprocity_abserr"] = abs(effx - efft)
                reciprocity_worst = max(reciprocity_worst, geo["reciprocity_abserr"])
                geo_by_sched[sch] = geo
            # map probe names -> R34 names for the frozen record
            rename = {"paired": "paired", "random": "random", "split": "ordered"}
            geo_named = {rename[s]: geo_by_sched[s] for s in geo_by_sched}
            kfrob = {rename[s]: geo_by_sched[s]["K_frob"] for s in geo_by_sched}
            kspread = ((max(kfrob.values()) - min(kfrob.values()))
                       / (np.mean(list(kfrob.values())) + 1e-30))
            per_cell[tag] = dict(
                geometry=geo_named,
                K_frob=kfrob,
                K_frob_spread_frac=float(kspread),
                # FROZEN predictions:
                predicted_best_scene_schedule="paired",
                predicted_best_medium_schedule="paired",
                predicted_medium_schedule_robust=True,
                predicted_no_forced_trade=True,
                predicted_fundamental_info_schedule_invariant=bool(kspread < 0.05))
    return dict(generic_theorem2=generic,
                reciprocity_worst_abserr=float(reciprocity_worst),
                per_cell=per_cell)


def main():
    os.makedirs(OUT, exist_ok=True)
    res = compute()
    meta = dict(utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"),
                python=platform.python_version(), numpy=np.__version__,
                machinery="dl_t3_reciprocity (frozen probe, verbatim import)",
                compact_witness=COMPACT,
                snr_axis_note="reciprocity structure is SNR-invariant; predictions "
                              "stated per (t_c,CV) hold across the photon-SNR grid",
                schedule_names="R34 {paired, random, ordered}; ordered == probe 'split'",
                predicted_ordering=("full-record Fisher info ~schedule-invariant; "
                                    "paired predicted best for BOTH ledgers at the "
                                    "estimator level; no forced scene-vs-medium trade"),
                status="FROZEN before confirmatory data (R34 sec.2.4)")
    blob = dict(meta=meta, **res)
    fn = os.path.join(OUT, "ABK_PREDICTIONS.json")
    json.dump(blob, open(fn, "w"), indent=2)
    g = res["generic_theorem2"]
    print("=== (A,B,K) preregistered predictions ===", flush=True)
    print(f"generic Thm2 identities: relerr Jx={g['worst_relerr_Jx']:.1e} "
          f"Jtheta={g['worst_relerr_Jtheta']:.1e} det-recip abserr="
          f"{g['worst_abserr_det_reciprocity']:.1e}", flush=True)
    print(f"GI-Fisher reciprocity worst abserr over cells/schedules: "
          f"{res['reciprocity_worst_abserr']:.2e}", flush=True)
    for tag, d in res["per_cell"].items():
        print(f"  {tag}: ||K|| paired={d['K_frob']['paired']:.3e} "
              f"random={d['K_frob']['random']:.3e} ordered={d['K_frob']['ordered']:.3e} "
              f"spread={100*d['K_frob_spread_frac']:.1f}%  "
              f"pred_best=paired invariant={d['predicted_fundamental_info_schedule_invariant']}",
              flush=True)
    print(f"wrote {fn}", flush=True)


if __name__ == "__main__":
    main()
