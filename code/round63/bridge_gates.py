"""ROUND63 DEV bridge — kill-gate evaluation (T-D deliverable d).

Decision rules ONLY (R24 §3.4 verbatim; Gate C replaced by R25 §12; Gate D
renamed + reclassified by R27 §3; the four-gap decomposition R24 §3.5).  No
constant is tunable: every threshold below is quoted from the ruling.

R27 §3 decision-tree split:
  * Gates A, B, C are the METHOD-SCIENCE go/no-go for M2 — if any fails, no M2.
  * The former Gate D is renamed `PLUGIN_LATENCY_PASS` and is a CLAIM-
    CLASSIFICATION gate only: it NEVER kills M2.  PASS (median<1.0 s & p95<2.0 s)
    -> the method may be described as subsecond/plugin-like; FAIL -> M2 may still
    launch if A-C pass, but the method is described as BATCH/TWO-SHOT adaptive
    illumination with the full latency distribution disclosed (median/p95/max/
    CPU/threads/cache state).  No substitute threshold.

These functions consume a flat list of per-cell result records and return
verdicts; they are unit-tested on SYNTHETIC outcomes (test_bridge_gates.py) and
are NOT run on real bridge data before the grid.

Record schema (one per scene x arm x nu x c x seed)
  {"scene": str, "group": "contour|twopop|microtex|control",
   "arm": "SCAT32-060|RIDGE-SCAT32|TRUE-X-FW|XHAT-FW|RLMI|ORACLE-LIB",
   "nu": float, "c": float, "seed": int, "PSNR_rad": float,
   # RLMI only:
   "A_realized": float,            # R25.4  A_j = 1 - w_hat_{j,0}
   "w0_realized": float,           # w_hat_{j,0}
   "route_latency_s": float,       # Gate D (rlmi alloc_time_s, ex-recon)
   # ORACLE-LIB only:
   "oracle_best_psnr": float,      # best bank PSNR_rad (the arm's reported value)
   "oracle_knob_psnr": float}      # the L0 (knob) bank sub-acquisition PSNR_rad

The hard bridge decision cell is (c=0.05, nu=2000) (R24 §3.2); other cells are
diagnostic.  Per-scene values aggregate the 5 paired seeds by MEAN (paired
radiometric PSNR); Q_base,j = max{Q_SCAT32-060,j, Q_RIDGE-SCAT32,j} (R24 §3.3).
"""
from __future__ import annotations

import numpy as np

STRESS_GROUPS = ("contour", "twopop", "microtex")   # 12 misalignment-stress scenes
CONTROL_GROUP = "control"                            # 4 aligned controls
DECISION_NU = 2000.0
DECISION_C = 0.05
BOOT_CONF = 0.90                                     # one-sided lower-bound level
BOOT_B = 2000


# --------------------------------------------------------------------------- #
#  aggregation
# --------------------------------------------------------------------------- #
def _sel(records, nu, c):
    return [r for r in records if r["nu"] == nu and r["c"] == c]


def scene_means(records, arm, nu, c, field="PSNR_rad"):
    """{scene: mean over the paired seeds of `field` for `arm`} at (nu,c)."""
    acc = {}
    for r in _sel(records, nu, c):
        if r["arm"] != arm:
            continue
        acc.setdefault(r["scene"], []).append(float(r[field]))
    return {sc: float(np.mean(v)) for sc, v in acc.items()}


def scene_groups_map(records):
    return {r["scene"]: r["group"] for r in records}


def q_base(records, nu, c):
    """Q_base,j = max{Q_SCAT32-060,j, Q_RIDGE-SCAT32,j} per scene (R24 §3.3)."""
    a = scene_means(records, "SCAT32-060", nu, c)
    b = scene_means(records, "RIDGE-SCAT32", nu, c)
    return {sc: max(a[sc], b[sc]) for sc in a if sc in b}


def paired_gains(records, arm, nu, c, groups=STRESS_GROUPS):
    """Per-scene gain Q_arm,j - Q_base,j restricted to `groups`.  Returns
    {scene: gain} (scene means; the pairing is the shared Q_base per scene)."""
    gmap = scene_groups_map(records)
    qa = scene_means(records, arm, nu, c)
    qb = q_base(records, nu, c)
    return {sc: qa[sc] - qb[sc] for sc in qa
            if sc in qb and gmap.get(sc) in groups}


# --------------------------------------------------------------------------- #
#  scene bootstrap (lower bound on the median)
# --------------------------------------------------------------------------- #
def scene_bootstrap_median_lb(values, conf=BOOT_CONF, B=BOOT_B, seed=650000):
    """One-sided `conf` lower bound on the median by resampling the SCENES with
    replacement (R24 §3.4 "90% scene-bootstrap lower bound on the median").
    LB = the (1-conf) quantile of the bootstrap median distribution."""
    v = np.asarray(values, dtype=np.float64)
    if v.size == 0:
        return float("nan")
    rng = np.random.default_rng(seed)
    meds = np.empty(B)
    for b in range(B):
        idx = rng.integers(0, v.size, size=v.size)
        meds[b] = np.median(v[idx])
    return float(np.quantile(meds, 1.0 - conf))


# --------------------------------------------------------------------------- #
#  LIBRARY_REACHABILITY_PASS — finite-library image-reachability gate (R28 §1)
# --------------------------------------------------------------------------- #
#  R28 replaces R24's Gate A: the science go/no-go now tests whether the DECLARED
#  finite library contains a robust image-rescuing design (the "finite-library
#  image oracle" / "declared-library reachability reference" — NEVER "class
#  ceiling").  It reads ONLY the genuine, exact-972, postchecked ORACLE-LIB bank
#  materializations, never an FW arm.  Per scene the oracle bank is chosen by the
#  FIVE-SEED MEAN (one bank per scene, min-index tie-break); no bank is chosen
#  per noise seed (R28 §1).
FW_ARMS = ("TRUE-X-FW", "XHAT-FW")                    # descriptive only; never gated
LIBRARY_BANKS = tuple("L%d" % k for k in range(8))
KNOB_BANK = "L0"


def oracle_bank_scene_means(records, nu, c):
    """{scene: {bank: five-seed-mean PSNR_rad}} from the ORACLE-LIB per_bank
    disclosures (each ORACLE-LIB record is one (scene, seed) with per_bank =
    [{bank, PSNR_rad}, ...] over the 8 genuinely-materialized library banks)."""
    acc = {}
    for r in _sel(records, nu, c):
        if r["arm"] != "ORACLE-LIB":
            continue
        for e in (r.get("per_bank") or []):
            (acc.setdefault(r["scene"], {})
                .setdefault(e["bank"], []).append(float(e["PSNR_rad"])))
    return {sc: {b: float(np.mean(v)) for b, v in banks.items()}
            for sc, banks in acc.items()}


def oracle_lib_reference(records, nu, c):
    """Per scene (R28 §1): k*_j = min arg max_k mean-seed Q_{j,k}; the oracle value
    Q_ORACLE-LIB,j = mean-seed Q_{j,k*_j}; the reachability gain
    dQ^A_j = Q_ORACLE-LIB,j - Q_base,j.  Returns {scene: {kstar,q_oracle,q_base,
    gain}}."""
    bm = oracle_bank_scene_means(records, nu, c)
    qb = q_base(records, nu, c)
    out = {}
    for sc, banks in bm.items():
        if sc not in qb:
            continue
        kstar, best = None, -np.inf
        for k in range(8):                            # L0..L7 -> min-index tie-break
            name = "L%d" % k
            if name in banks and banks[name] > best + 1e-12:
                best, kstar = banks[name], name
        if kstar is None:
            continue
        out[sc] = {"kstar": kstar, "q_oracle": best, "q_base": qb[sc],
                   "gain": best - qb[sc]}
    return out


def _oracle_stress_gains(records, nu, c):
    ref = oracle_lib_reference(records, nu, c)
    gmap = scene_groups_map(records)
    return {sc: ref[sc]["gain"] for sc in ref if gmap.get(sc) in STRESS_GROUPS}, ref


def gate_library_reachability(records, nu=DECISION_NU, c=DECISION_C,
                              boot_seed=650000):
    """LIBRARY_REACHABILITY_PASS (R28 §1) — the finite-library image oracle gains
    over the fixed composite baseline on the 12 stress scenes.  PASS iff median
    >= 1.0 dB AND >= 9/12 positive AND the frozen 90% scene-bootstrap LB > 0.50 dB.
    Fail -> the declared bank has no robust image-rescuing design; RLMI/M2 stop."""
    stress, ref = _oracle_stress_gains(records, nu, c)
    scenes = sorted(stress)
    vals = np.array([stress[s] for s in scenes])
    med = float(np.median(vals)) if vals.size else float("nan")
    n_pos = int((vals > 0).sum())
    lb = scene_bootstrap_median_lb(vals, seed=boot_seed)
    c1, c2, c3 = med >= 1.0, n_pos >= 9, lb > 0.50
    return {"gate": "LIBRARY_REACHABILITY_PASS", "nu": nu, "c": c,
            "n_scenes": len(scenes), "median_gain_dB": med, "n_positive": n_pos,
            "boot_median_lb_dB": lb, "c1_median>=1.0": bool(c1),
            "c2_>=9of12_positive": bool(c2), "c3_bootlb>0.50": bool(c3),
            "PASS": bool(c1 and c2 and c3), "per_scene": stress,
            "kstar_per_scene": {s: ref[s]["kstar"] for s in scenes},
            "reference_terminology": "finite-library image oracle "
            "(declared-library reachability reference) -- NOT a class ceiling"}


# back-compat alias: the M2 go/no-go 'Gate A' is now the library-reachability gate
def gate_A(records, nu=DECISION_NU, c=DECISION_C, boot_seed=650000):
    return gate_library_reachability(records, nu, c, boot_seed)


# --------------------------------------------------------------------------- #
#  Gate B — the deployable router retains the gain (R24 §3.4)
# --------------------------------------------------------------------------- #
def gate_B(records, nu=DECISION_NU, c=DECISION_C, boot_seed=650001):
    gr = paired_gains(records, "RLMI", nu, c, STRESS_GROUPS)
    scenes = sorted(gr)
    vals = np.array([gr[s] for s in scenes])
    med = float(np.median(vals)) if vals.size else float("nan")
    n_pos = int((vals > 0).sum())
    lb = scene_bootstrap_median_lb(vals, seed=boot_seed)
    # R28 §2: the 60% capture denominator is the ORACLE-LIB reachability gain, NOT
    # TRUE-X-FW (an FW arm never enters a gate).
    ostress, _ = _oracle_stress_gains(records, nu, c)
    med_oracle = (float(np.median([ostress[s] for s in sorted(ostress)]))
                  if ostress else float("nan"))
    # controls
    gc = paired_gains(records, "RLMI", nu, c, (CONTROL_GROUP,))
    cvals = np.array([gc[s] for s in sorted(gc)])
    med_ctrl = float(np.median(cvals)) if cvals.size else float("nan")
    n_below1 = int((cvals < -1.0).sum())
    c1 = med >= 1.0
    c2 = n_pos >= 9
    c3 = lb > 0.0
    c4 = med >= 0.60 * med_oracle if np.isfinite(med_oracle) else False
    c5 = (med_ctrl >= -0.25) and (n_below1 <= 1)
    return {"gate": "B", "nu": nu, "c": c, "n_stress": len(scenes),
            "median_gain_dB": med, "n_positive": n_pos, "boot_median_lb_dB": lb,
            "median_oracle_dB": med_oracle, "60pct_oracle_dB": 0.60 * med_oracle,
            "control_median_dB": med_ctrl, "n_controls_below_-1dB": n_below1,
            "c1_median>=1.0": bool(c1), "c2_>=9of12_positive": bool(c2),
            "c3_bootlb>0": bool(c3), "c4_>=60pct_oracle": bool(c4),
            "c5_noharm_controls": bool(c5),
            "PASS": bool(c1 and c2 and c3 and c4 and c5),
            "per_scene_stress": gr, "per_scene_control": gc}


# --------------------------------------------------------------------------- #
#  Gate C — allocation informativeness (R25 §12; A_j from materialized weights)
# --------------------------------------------------------------------------- #
def _oracle_knob_advantage(records, nu, c):
    """Per-scene ORACLE-LIB - knob(L0) advantage (DEV-only label; R25 §12), from
    the per-bank five-seed means: best = max_k mean-seed Q_{j,k}, knob = mean-seed
    Q_{j,L0}.  Consistent with the R28 §1 oracle reference (never reads an FW arm)."""
    bm = oracle_bank_scene_means(records, nu, c)
    out = {}
    for sc, banks in bm.items():
        if KNOB_BANK in banks and banks:
            out[sc] = max(banks.values()) - banks[KNOB_BANK]
    return out


def gate_C(records, nu=DECISION_NU, c=DECISION_C):
    gmap = scene_groups_map(records)
    adv = _oracle_knob_advantage(records, nu, c)
    A = scene_means(records, "RLMI", nu, c, field="A_realized")
    w0 = scene_means(records, "RLMI", nu, c, field="w0_realized")
    # C1: rescue-needed scenes (ORACLE-LIB beats knob by > 1 dB)
    resc = [sc for sc in adv if adv[sc] > 1.0]
    A_resc = [A[sc] for sc in resc if sc in A]
    meanA = float(np.mean(A_resc)) if A_resc else float("nan")
    c1 = (meanA >= 0.80) if A_resc else False
    # C2: aligned controls where ORACLE-LIB advantage < 0.25 dB
    ctrl = [sc for sc in adv if gmap.get(sc) == CONTROL_GROUP and adv[sc] < 0.25]
    w0_ctrl = [w0[sc] for sc in ctrl if sc in w0]
    meanW0 = float(np.mean(w0_ctrl)) if w0_ctrl else float("nan")
    c2 = (meanW0 >= 0.75) if w0_ctrl else False
    return {"gate": "C", "nu": nu, "c": c,
            "n_rescue_needed": len(resc), "mean_A_j": meanA,
            "n_aligned_lowadv": len(ctrl), "mean_w0_j": meanW0,
            "c1_meanA>=0.80": bool(c1), "c2_meanw0>=0.75": bool(c2),
            "PASS": bool(c1 and c2),
            "rescue_scenes": resc, "aligned_control_scenes": ctrl}


# --------------------------------------------------------------------------- #
#  PLUGIN_LATENCY_PASS — route-latency CLAIM-CLASSIFICATION gate (R27 §3)
# --------------------------------------------------------------------------- #
#  R24's Gate D is renamed and reclassified by R27: it is NOT a science go/no-go
#  and NEVER kills M2.  It only classifies the deployment claim:
#    PASS  (median < 1.0 s AND p95 < 2.0 s) -> the method may be described as
#          subsecond / plugin-like / real-time.
#    FAIL  -> M2 may still launch if A-C pass, but the method must be described as
#          a BATCH / TWO-SHOT adaptive illumination method (never real-time /
#          subsecond).  The full latency distribution is disclosed.
#  There is no substitute threshold (R27 §3.1): the measured latency is reported
#  as-is, never converted into a new pass bar.
def gate_plugin_latency(records, field="route_latency_s", cpu=None, threads=None,
                        cache_state=None):
    lat = [float(r[field]) for r in records
           if r["arm"] == "RLMI" and r.get(field) is not None]
    if not lat:
        return {"gate": "PLUGIN_LATENCY_PASS", "PASS": False,
                "reason": "no RLMI latencies", "gates_M2": False,
                "median_s": float("nan"), "p95_s": float("nan"),
                "max_s": float("nan"), "claim_label": "UNKNOWN"}
    med = float(np.median(lat))
    p95 = float(np.percentile(lat, 95))
    mx = float(np.max(lat))
    plugin = (med < 1.0) and (p95 < 2.0)
    return {"gate": "PLUGIN_LATENCY_PASS", "n": len(lat),
            "median_s": med, "p95_s": p95, "max_s": mx,
            "c1_median<1.0s": bool(med < 1.0), "c2_p95<2.0s": bool(p95 < 2.0),
            "PASS": bool(plugin), "gates_M2": False,      # R27 §3: never kills M2
            "claim_label": ("subsecond / plugin-like / real-time" if plugin
                            else "batch / two-shot adaptive illumination"),
            # full distribution disclosure (R27 §3.1)
            "cpu": cpu, "threads": threads, "cache_state": cache_state}


# back-compat alias (old name) — same object, claim-classification semantics
def gate_D(records, **kw):
    return gate_plugin_latency(records, **kw)


# --------------------------------------------------------------------------- #
#  four-gap decomposition (R24 §3.5) — reported, not gated
# --------------------------------------------------------------------------- #
def _median_paired_diff(records, arm_hi, arm_lo, nu, c, groups=STRESS_GROUPS):
    a = scene_means(records, arm_hi, nu, c)
    b = scene_means(records, arm_lo, nu, c)
    gmap = scene_groups_map(records)
    diffs = [a[sc] - b[sc] for sc in a
             if sc in b and gmap.get(sc) in groups]
    return float(np.median(diffs)) if diffs else float("nan")


def four_gap(records, nu=DECISION_NU, c=DECISION_C):
    """R28 §4 four-quantity decomposition (medians over the 12 stress scenes);
    the FW-based quantities 3-4 are DESCRIPTIVE (dose-relaxed, nondeployable) and
    NEVER gate anything:
      1 in-class finite-library reachability = ORACLE-LIB - fixed baseline (Gate A);
      2 deployable allocation loss           = RLMI - ORACLE-LIB (Gate B/capture);
      3 plug-in FW loss                      = XHAT-FW - TRUE-X-FW (descriptive);
      4 dose-relaxation/library contrast     = TRUE-X-FW - ORACLE-LIB (descriptive;
        conflates dose-band removal + FW surrogate + finite bank -- NOT a pure
        finite-bank or Fisher-to-image gap).
    Optional descriptive: TRUE-X-FW - baseline (nondeployable surrogate-to-image).
    Quantities 3,4,opt read the NONDEPLOYABLE_DOSE_RELAXED_FW_DIAGNOSTIC arms."""
    ostress, _ = _oracle_stress_gains(records, nu, c)
    q1 = float(np.median(list(ostress.values()))) if ostress else float("nan")
    # RLMI - ORACLE-LIB per scene = (RLMI - base) - (ORACLE - base)
    rlmi_g = paired_gains(records, "RLMI", nu, c, STRESS_GROUPS)
    q2 = (float(np.median([rlmi_g[s] - ostress[s] for s in ostress if s in rlmi_g]))
          if ostress else float("nan"))
    q3 = _median_paired_diff(records, "XHAT-FW", "TRUE-X-FW", nu, c)   # descriptive
    q4 = _median_paired_diff(records, "TRUE-X-FW", "ORACLE-LIB", nu, c)  # descriptive
    txb = float(np.median(list(paired_gains(
        records, "TRUE-X-FW", nu, c, STRESS_GROUPS).values() or [np.nan])))
    return {"nu": nu, "c": c,
            "q1_inclass_library_reachability_dB": q1,      # ORACLE-LIB - base (Gate A)
            "q2_deployable_allocation_loss_dB": q2,        # RLMI - ORACLE-LIB
            "q3_plugin_fw_loss_dB_descriptive": q3,        # XHAT-FW - TRUE-X-FW
            "q4_dose_relax_library_contrast_dB_descriptive": q4,  # TRUE-X-FW - ORACLE-LIB
            "opt_truex_minus_baseline_dB_nondeployable": txb,
            "note": "q3,q4,opt are NONDEPLOYABLE_DOSE_RELAXED_FW diagnostics; not gated"}


# --------------------------------------------------------------------------- #
#  R28 §5 smoke-exposure provenance + leave-out sensitivity (no-gate)
# --------------------------------------------------------------------------- #
SMOKE_EXPOSED_PREGRID = ("bridge_twopop_0", "bridge_control_0")


def leave_smoke_out_sensitivity(records, nu=DECISION_NU, c=DECISION_C):
    """No-gate (R28 §5): recompute the library-reachability + Gate B verdicts on
    the cohort MINUS the pre-grid smoke-exposed scenes, purely as a robustness
    disclosure.  The OFFICIAL gates always use the complete frozen cohort."""
    filt = [r for r in records if r["scene"] not in SMOKE_EXPOSED_PREGRID]
    return {"excluded_smoke_exposed": list(SMOKE_EXPOSED_PREGRID),
            "gate_A_leaveout": gate_library_reachability(filt, nu, c),
            "gate_B_leaveout": gate_B(filt, nu, c),
            "disclosure_only": True}


# --------------------------------------------------------------------------- #
#  R28 §3/§7 invariant: no Gate A-C statistic may read an FW arm
# --------------------------------------------------------------------------- #
def gates_read_no_fw(records, nu=DECISION_NU, c=DECISION_C):
    """Return True iff perturbing the FW arms (TRUE-X-FW / XHAT-FW) to extreme
    values leaves gate_library_reachability, gate_B, and gate_C bit-identical
    (R28 §3: no science gate reads a dose-relaxed FW diagnostic)."""
    def key(rec):
        return (gate_library_reachability(rec, nu, c)["PASS"],
                gate_library_reachability(rec, nu, c)["median_gain_dB"],
                gate_B(rec, nu, c)["PASS"], gate_B(rec, nu, c)["median_oracle_dB"],
                gate_C(rec, nu, c)["PASS"], gate_C(rec, nu, c).get("mean_A_j"))
    base = key(records)
    perturbed = []
    for r in records:
        r2 = dict(r)
        if r2.get("arm") in FW_ARMS:
            r2["PSNR_rad"] = -999.0
        perturbed.append(r2)
    return _key_eq(base, key(perturbed))


def _key_eq(a, b):
    for x, y in zip(a, b):
        if isinstance(x, float) and isinstance(y, float):
            if not (np.isnan(x) and np.isnan(y)) and abs(x - y) > 1e-12:
                return False
        elif x != y:
            return False
    return True


def evaluate_all(records, nu=DECISION_NU, c=DECISION_C, cpu=None, threads=None,
                 cache_state=None):
    """Full gate report at the hard decision cell.

    R27 §3 decision tree: Gates A, B, C are the method-science go/no-go for M2;
    PLUGIN_LATENCY_PASS (the former Gate D) NEVER gates M2 — it only classifies
    the deployment claim (subsecond/plugin vs batch/two-shot).  The four-gap
    decomposition (R24 §3.5) is reported, never gated."""
    A = gate_library_reachability(records, nu, c)     # R28 §1 (former Gate A)
    B = gate_B(records, nu, c)
    C = gate_C(records, nu, c)
    L = gate_plugin_latency(records, cpu=cpu, threads=threads,
                            cache_state=cache_state)
    gaps = four_gap(records, nu, c)
    # M2 go/no-go depends ONLY on A(library-reachability), B, C (R27 §3.1 / R28).
    # No FW arm enters this decision (R28 §3) -- verified by gates_read_no_fw.
    if not A["PASS"]:
        m2 = ("LIBRARY_REACHABILITY_FAIL -> declared bank has no robust image-"
              "rescuing design; RLMI/M2 stop (R28 §1)")
    elif not B["PASS"]:
        m2 = "GATE_B_FAIL -> do not launch M2 (no rescue redesign; R24 hard stop)"
    elif not C["PASS"]:
        m2 = ("GATE_C_FAIL -> A,B pass but allocation not informative: may be a "
              "static composite design; MUST NOT lead M2 as adaptive (R25 §12)")
    else:
        m2 = "GATES A-C PASS -> M2 preregistration; RLMI leads paper 2"
    # The latency result only sets the deployment claim label (R27 §3).
    claim = ("PLUGIN_LATENCY_PASS -> method may be described as subsecond/plugin"
             if L["PASS"] else
             "PLUGIN_LATENCY_FAIL -> describe as BATCH/TWO-SHOT adaptive "
             "illumination (not real-time/subsecond); disclose full latency dist")
    return {"decision_cell": {"nu": nu, "c": c},
            "gate_library_reachability": A, "gate_A": A,  # alias
            "gate_B": B, "gate_C": C, "plugin_latency": L, "four_gap": gaps,
            "no_gate_reads_fw": bool(gates_read_no_fw(records, nu, c)),
            "leave_smoke_out_sensitivity": leave_smoke_out_sensitivity(records, nu, c),
            "m2_verdict": m2, "claim_label": claim,
            "verdict": "%s | %s" % (m2, claim)}
