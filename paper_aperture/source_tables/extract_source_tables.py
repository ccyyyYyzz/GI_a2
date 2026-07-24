#!/usr/bin/env python3
"""
FROZEN SOURCE-TABLE EXTRACTION  (R40 flagship: "The information aperture of a bucket record")

RULE ZERO (R40 ruling, writing order step 3): every plotted number must have a
machine-readable source table; no prose values typed by hand into a figure script.

This script reads the ORIGINAL repo evidence files and emits one table per figure
panel into paper_aperture/source_tables/.  Each table carries a PROVENANCE field
naming the exact repo file (and section, where the value lives in prose/tex) it
came from.  Figures read ONLY from source_tables/, never from the raw evidence.

Two provenance classes are distinguished honestly:
  (A) COPIED   - transcribed programmatically from a machine-readable JSON/CSV in
                 this repo (ridge curves, aperture, P1, 81-cell map, pocket).
  (B) FROZEN   - transcribed from a frozen prose/tex document in this repo because
                 the canonical machine-readable source lives outside this repo
                 (cohort R19 statistics, five-intervention closure).  Every such
                 value cites file + section and is a frozen-verbatim number.

Run:  python paper_aperture/source_tables/extract_source_tables.py
      (cwd = repo root  D:\\GI_another)
"""
import csv, json, os, sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT  = os.path.join(REPO, "paper_aperture", "source_tables")
os.makedirs(OUT, exist_ok=True)

def rp(*p):  # repo path
    return os.path.join(REPO, *p)

def dump_json(name, obj):
    with open(os.path.join(OUT, name), "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)
    print("  wrote", name)

# ---------------------------------------------------------------------------
# T_fig2_ridge  (curves + peaks + operating point)  -- CLASS A (+ frozen law)
# ---------------------------------------------------------------------------
def t_fig2_ridge():
    src_csv = rp("results", "round63_theory", "fisher_ridge.csv")
    src_jsn = rp("results", "round63_theory", "fisher_ridge.json")
    # --- curves: copy fisher_ridge.csv verbatim into the panel table ---
    rows = list(csv.DictReader(open(src_csv)))
    out_csv = os.path.join(OUT, "T_fig2_ridge_curves.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["# PROVENANCE: results/round63_theory/fisher_ridge.csv (COPIED verbatim)"])
        w.writerow(["# I_exact = exact count-Fisher information about log-lambda per unit time (tau units);"])
        w.writerow(["# it RISES to a peak at rho_star then DECLINES (dead-time count saturation)."])
        w.writerow(["# I_clt   = CLT proxy rho/(1+rho) (no-dead-time reference, monotone to 1)."])
        w.writerow(["nu", "rho", "I_exact", "I_clt", "ratio"])
        for r in rows:
            w.writerow([r["nu"], r["rho"], r["I_exact_log"], r["I_clt_log"], r["ratio"]])
    print("  wrote T_fig2_ridge_curves.csv")

    # --- peaks (ridge locus) + operating point ---
    fj = json.load(open(src_jsn))
    ridge = fj["ridge"]
    peaks = []
    for nu, d in ridge.items():
        peaks.append({
            "nu": int(nu),
            "rho_star": d["rho_star"],          # deterministic ridge peak location
            "I_at_peak": d["I_at_peak"],         # peak per-unit-time info
            "I_asymptote_CLT": d["I_asymptote_CLT"],
            "qmle_trust_rho_ratio0.9": d["qmle_trust_rho_ratio0.9"],
        })
    peaks.sort(key=lambda x: x["nu"])
    obj = {
        "panel": "Fig.2a / Fig.1 mean-row : hidden-state information ridge",
        "PROVENANCE": {
            "curves": "results/round63_theory/fisher_ridge.csv (COPIED)",
            "peaks":  "results/round63_theory/fisher_ridge.json -> ridge{} (COPIED)",
            "ridge_law": "paper2/main_m1.tex sec4 + docs/FLAGSHIP_MATERIALS_MAP.md sec1a (FROZEN)",
            "jitter_cap": "paper2/main_m1.tex sec3 Thm3 + supplement_m1.tex S4 (FROZEN)",
            "operating_point": "paper2/main_m1.tex sec4/sec6, R19-frozen (FROZEN)",
        },
        "definition": fj.get("definition", ""),
        "ridge_peaks": peaks,
        "ridge_law": {
            "rho_star_of_nu": "(6*nu)**(1/3) - 2/3 + O(nu**-1/3)",
            "J_at_peak_of_nu": "1 - 0.8255*nu**(-1/3)",
            "class": "FROZEN (deterministic dead-time ridge law, paper2 sec4)",
        },
        "jitter_cap": {
            "note": "finite illumination jitter c_v caps the reachable load below the deterministic ridge",
            "optimum_condition": "c_v**2 * rho_c**2 * (1 + 2*rho_c) = 1",
            "small_jitter_expansion": "rho_c = 2**(-1/3) * c_v**(-2/3) - 1/6 + ...",
            "pooled_MC_slope": -0.658,
            "pooled_MC_slope_target": -2/3,
            "example_cv0p05_nu2000": {"det_rho_star": 22.3, "capped_rho_c": 5.7,
                                      "info_forfeited_frac": 0.55},
            "MC_peaks_by_nu_desc": [22.297, 10.739, 5.173, 3.862, 2.153, 1.607],
            "class": "FROZEN (paper2 supplement_m1.tex S4, jitter_sfi_v2 100k frames)",
        },
        "selected_operating_point": {
            # the M1 power-control policy operates on the deterministic ridge at nu=2000
            "nu": 2000,
            "rho_operating": 22.25,          # rho_R(2000) policy load (paper2 sec4)
            "rho_star_deterministic": 22.15503546089851,
            "I_at_peak": 0.9353615833281874,
            "label": "certified mean-channel operating point (living optimum)",
            "certified_result_ref": "T_fig2_cohorts.json (+1.87 dB / 19.13x)",
            "class": "FROZEN",
        },
    }
    dump_json("T_fig2_ridge_points.json", obj)

# ---------------------------------------------------------------------------
# T_fig2_cohorts  (two certified statistics)  -- CLASS B (FROZEN R19)
# ---------------------------------------------------------------------------
def t_fig2_cohorts():
    obj = {
        "panel": "Fig.2b : two certified mean-channel cohort statistics",
        "PROVENANCE": ("paper2/main_m1.tex sec6 + docs/FLAGSHIP_MATERIALS_MAP.md sec9a "
                       "(R19-frozen verbatim). Canonical machine source "
                       "results/round63_m1/CORRECTION_2026-07-19/M1_VERDICTS_SPEC_CORRECTED_R19.json "
                       "lives in the E:/GAN_FCC_WORK program repo, not GI_a2."),
        "class": "FROZEN",
        "materiality_bar_dB": 1.0,
        "cohorts": [
            {
                "key": "fixed_dwell_quality",
                "label": "fixed-dwell radiometric quality gain",
                "metric": "median dQ over matched 0.60 comparator (dB)",
                "value": 1.87,
                "lower_bound": 1.41,      # family-stratified LB
                "lb_kind": "family-stratified",
                "n_positive": 19, "n_total": 24,
                "unit": "dB",
            },
            {
                "key": "elapsed_speed",
                "label": "elapsed-time speed ratio",
                "metric": "median conservative speed ratio S_gate (x)",
                "value": 19.13,
                "value_exact": 19.127043091646133,
                "lower_bound": 18.33,
                "lower_bound_exact": 18.328492357080282,
                "lb_kind": "family-stratified nested-bootstrap 2.5th pct",
                "n_positive": 21, "n_total": 24,
                "unit": "x",
            },
        ],
        "resource_frame": {
            "note": "benefits bought with power, not photon efficiency",
            "incident_dose_x": 37.1, "detected_counts_x": 2.6,
        },
    }
    dump_json("T_fig2_cohorts.json", obj)

# ---------------------------------------------------------------------------
# T_fig2_closure  (five preregistered interventions)  -- CLASS B (FROZEN)
# ---------------------------------------------------------------------------
def t_fig2_closure():
    # values frozen in docs/SOFTWARE_SATURATION_VERDICT.md (table) and
    # docs/FLAGSHIP_MATERIALS_MAP.md sec3a.
    rows = [
        # key, label, verdict, governing_number, unit, extra
        ["bridge",    "scene-adaptive geometry",     "FAIL",      0.68,  "dB",
         "median dQ +0.680 < 1.0; boot LB +0.543; 11/12 pos"],
        ["allocator", "robust maximin routing",      "HARM",     -8.3,   "dB",
         "allocation loss -8.33 dB (harmful arm)"],
        ["dops",      "moment-orthogonal schedule",  "FAIL",      0.039, "dB",
         "median dQ +0.039 < 1.0; 5/6 pos"],
        ["probeA",    "estimator vs Fisher ceiling", "SATURATED", 0.20,  "dB",
         "~99% var-efficient; <=0.2 dB residual headroom"],
        ["cpl",       "conditional-likelihood gain", "KILL",      0.412, "dB",
         "median dQ +0.412 < 1.0; 5/6 pos"],
    ]
    out_csv = os.path.join(OUT, "T_fig2_closure.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["# PROVENANCE: docs/SOFTWARE_SATURATION_VERDICT.md (table) + "
                    "docs/FLAGSHIP_MATERIALS_MAP.md sec3a  (FROZEN)"])
        w.writerow(["# materiality_bar_dB = 1.0 ; probeA governing_number is a headroom CEILING (<=0.2)"])
        w.writerow(["key", "mechanism", "verdict", "governing_number", "unit",
                    "materiality_bar_dB", "detail"])
        for r in rows:
            w.writerow([r[0], r[1], r[2], r[3], r[4], 1.0, r[5]])
    print("  wrote T_fig2_closure.csv")

# ---------------------------------------------------------------------------
# T_fig3_aperture  (aperture-law validation)  -- CLASS A (+ E8 duel frozen)
# ---------------------------------------------------------------------------
def t_fig3_aperture():
    g3 = json.load(open(rp("results", "round63_next", "FOG_DMD_PROBE64", "G3_results.json")))
    pa = g3["partA"]
    obj = {
        "panel": "Fig.3b : empirical aperture-law validation (in/out coverage separation)",
        "PROVENANCE": {
            "separation_64x64": "results/round63_next/FOG_DMD_PROBE64/G3_results.json -> partA (COPIED)",
            "mean_wall_duel":   "results/round63_next/FOG_DMD_PROBE/FOG_DMD_PROBE_REPORT.md sec14 (E8) (FROZEN)",
            "16x16_pair":       "NOT FOUND in repo (see missing_inputs)",
        },
        "geometry_64x64": {"n": g3["n"], "N": g3["N"], "PB": g3["PB"],
                            "medium_band": g3["medium_band"], "M": g3["M"], "T": g3["T"]},
        "separation_64x64": {
            "in_coverage_err":  pa["in_coverage_err"],
            "out_coverage_err": pa["out_coverage_err"],
            "separation":       pa["separation"],
            "in_all":  pa["in_all"],
            "out_all": pa["out_all"],
            "class": "COPIED",
        },
        # E8 physics-wall in-band duel = the EXACT first-moment wall the covariance
        # aperture is measured against.  A band-limited operator cannot leave its
        # band, so fresh band-limited patterns are pinned at rel-err 1.000 on
        # beyond-band content while the fixed-bank covariance route recovers it.
        "mean_wall_duel_16x16": {
            "note": "16x16, band-limited patterns M=24, beyond-band rel-err (1.000 = recovers nothing)",
            "T_eff":                 [512, 1024, 2048],
            "fresh_bandlimited_wall":[1.000, 1.000, 1.000],
            "fixed_bank_moment":     [0.650, 0.547, 0.480],
            "fixed_bank_mle":        [0.505, 0.439, 0.481],
            "oracle_medium_known":   [0.389, 0.221, 0.091],
            "class": "FROZEN",
        },
        "missing_inputs": {
            "requested_16x16_in_out": {"in": 0.073, "out": 1.213},
            "status": "UNRECOVERABLE",
            "note": ("The R40 task listed a 16x16 in-coverage 0.073 vs out 1.213 pair. "
                     "No committed result in D:/GI_another produces this pair; the digits "
                     "0.073 and 1.213 appear only coincidentally in unrelated scatter/probe "
                     "JSONs (paper3_scatter T_s and study2 select_max_s). The canonical "
                     "aperture-law validation retained here is the G3 64x64 18.76x separation "
                     "plus the E8 exact-mean-wall duel."),
        },
    }
    dump_json("T_fig3_aperture.json", obj)

# ---------------------------------------------------------------------------
# T_fig3_p1  (P1 prognosis + T_eff^-1/2 scaling)  -- CLASS A
# ---------------------------------------------------------------------------
def t_fig3_p1():
    p1 = json.load(open(rp("results", "round63_next", "FOG_DMD_PROBE64", "P1_results.json")))
    prim = {d["tag"]: d for d in p1["primary_point"]}
    sweep = sorted(p1["T_eff_sweep"], key=lambda d: d["T_eff"])
    obj = {
        "panel": "Fig.3c : P1 profiled-Fisher prognosis (prelaunch kill) + T_eff^-1/2 scaling",
        "PROVENANCE": "results/round63_next/FOG_DMD_PROBE64/P1_results.json (COPIED)",
        "class": "COPIED",
        "frozen_geometry": p1["frozen_geometry"],
        "bar": {"f_rec_snr3_required": 0.70,
                "nmse_crb_synth_max": 0.25, "nmse_crb_nat_median_max": 0.35,
                "eta_prof_p10_min": 0.10, "eta_prof_median_min": 0.25},
        "primary_point": {
            "witness":       {"f_rec_snr3": prim["witness"]["f_rec_snr3"],
                              "nmse_crb": prim["witness"]["nmse_crb"],
                              "eta_prof_p10": prim["witness"]["eta_prof_p10"],
                              "eta_prof_median": prim["witness"]["eta_prof_median"],
                              "min_eig_Jbar": prim["witness"]["min_eig_Jbar"],
                              "exact_null": prim["witness"]["exact_null"]},
            "nat_cameraman": {"f_rec_snr3": prim["nat_cameraman"]["f_rec_snr3"],
                              "nmse_crb": prim["nat_cameraman"]["nmse_crb"],
                              "eta_prof_p10": prim["nat_cameraman"]["eta_prof_p10"]},
            "nat_coins":     {"f_rec_snr3": prim["nat_coins"]["f_rec_snr3"],
                              "nmse_crb": prim["nat_coins"]["nmse_crb"]},
            "nat_moon":      {"f_rec_snr3": prim["nat_moon"]["f_rec_snr3"],
                              "nmse_crb": prim["nat_moon"]["nmse_crb"]},
        },
        "natural_nrmse_crb_median": p1["natural_nrmse_crb_median"],
        "T_eff_sweep_witness": {
            "T_eff":   [d["T_eff"] for d in sweep],
            "nmse_crb":[d["nmse_crb"] for d in sweep],
            "f_rec_snr3":[d["f_rec_snr3"] for d in sweep],
            "note": "nmse_crb ~ T_eff^-1/2 (no floor); f_rec grows far too slowly to reach 0.70",
        },
        "checks": p1["checks"],
        "P1_pass": p1["P1_pass"],
        "verdict": p1["verdict"],
    }
    dump_json("T_fig3_p1.json", obj)

# ---------------------------------------------------------------------------
# T_fig3_map  (full 81-cell living-region grid)  -- CLASS A
# ---------------------------------------------------------------------------
def t_fig3_map():
    lm = json.load(open(rp("results", "round63_next", "FOG_DMD_PROBE64", "LIVING_REGION_MAP.json")))
    cols = ["sigma_f", "sigma_f_eff", "shape", "k_w_over_kp", "claim", "db",
            "f_rec_witness", "f_rec_nat_med", "nmse_crb_witness", "nmse_crb_nat_med",
            "eta_p10", "eta_med", "t_req_0p30", "exact_null", "P1_pass"]
    out_csv = os.path.join(OUT, "T_fig3_map.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["# PROVENANCE: results/round63_next/FOG_DMD_PROBE64/LIVING_REGION_MAP.json -> grid[] (COPIED)"])
        w.writerow([f"# n_cells={lm['n_cells']} n_pass={lm['n_pass']} verdict={lm['verdict']} "
                    f"(P1 bar: witness f_rec>=0.70, met by NONE)"])
        w.writerow(cols)
        for c in lm["grid"]:
            w.writerow([c[k] for k in cols])
    print("  wrote T_fig3_map.csv (%d cells)" % len(lm["grid"]))
    # sidecar json with grid meta + pocket pointer
    meta = {"PROVENANCE": "results/round63_next/FOG_DMD_PROBE64/LIVING_REGION_MAP.json (COPIED)",
            "class": "COPIED",
            "n_cells": lm["n_cells"], "n_pass": lm["n_pass"], "verdict": lm["verdict"],
            "fixed": lm["fixed"], "sigma_f_eff": lm["sigma_f_eff"],
            "axes": {"sigma_f": [0.3, 0.6, 1.0], "shape": ["flat", "k^-1", "k^-2"],
                     "k_w_over_kp": [1, 2, 4], "claim": [1.25, 1.5, 1.8]},
            "mc_cross_check": lm["mc_cross_check"]}
    dump_json("T_fig3_map_meta.json", meta)

# ---------------------------------------------------------------------------
# T_fig3_pocket  (pocket cell + PENDING branch)  -- CLASS A (+ md sec4)
# ---------------------------------------------------------------------------
def t_fig3_pocket():
    lm = json.load(open(rp("results", "round63_next", "FOG_DMD_PROBE64", "LIVING_REGION_MAP.json")))
    # pocket = sigma_f=0.3, k^-2, k_w/k_p=1, claim 1.25 (LIVING_REGION_MAP.md sec4 "most defensible")
    cell = next(c for c in lm["grid"] if c["sigma_f"] == 0.3 and c["shape"] == "k^-2"
                and c["k_w_over_kp"] == 1 and c["claim"] == 1.25)
    obj = {
        "panel": "Fig.3d : natural-scene 1.25x pocket (branch inset)",
        "PROVENANCE": ("results/round63_next/FOG_DMD_PROBE64/LIVING_REGION_MAP.json -> grid (COPIED) "
                       "+ LIVING_REGION_MAP.md sec4 'most defensible' (FROZEN wording)"),
        "class": "COPIED",
        "cell": {"sigma_f": cell["sigma_f"], "sigma_f_eff": cell["sigma_f_eff"],
                 "shape": cell["shape"], "k_w_over_kp": cell["k_w_over_kp"],
                 "claim_x_kp": cell["claim"], "db": cell["db"]},
        "f_rec_nat": cell["f_rec_nat_med"],       # 0.792  (~0.79)
        "f_rec_witness": cell["f_rec_witness"],   # 0.458  (~0.46, still < 0.70 bar)
        "nmse_crb_witness": cell["nmse_crb_witness"],  # 0.066 (~0.07)
        "nmse_crb_nat_med": cell["nmse_crb_nat_med"],  # 0.097 (~0.10)
        "t_req_0p30": cell["t_req_0p30"],         # 197.0 epochs
        "P1_pass": cell["P1_pass"],               # False (sub-P1, modest, prior-art-adjacent)
        "BRANCH": "PENDING",
        "branch_note": ("Blind lifted-GLS demo verdict arrives separately. "
                        "POCKET_ACHIEVED -> filled amber marker + one reconstruction inset; "
                        "ESTIMATOR_GAP_PERSISTS -> hollow amber marker "
                        "'Fisher-feasible, algorithmically unrealized'."),
    }
    dump_json("T_fig3_pocket.json", obj)

# ---------------------------------------------------------------------------
# T_fig1_atlas  (derived master-atlas summary)  -- DERIVED from the above
# ---------------------------------------------------------------------------
def t_fig1_atlas():
    # Units: spatial-frequency / task mode measured in multiples of the modulator
    # band edge k_p (so k_p == 1.0).  formal_support_extent and fisher_usable_extent
    # are the outer radii (in k_p units) of the pale and dark fills respectively.
    obj = {
        "panel": "Fig.1 : master information atlas (two rows x two apertures)",
        "PROVENANCE": ("DERIVED from T_fig2_ridge_points.json, T_fig3_aperture.json, "
                       "T_fig3_p1.json, T_fig3_map(_meta).json, T_fig3_pocket.json. "
                       "Each field's origin is in derivation_notes."),
        "class": "DERIVED",
        "axis": {"name": "spatial-frequency / task mode",
                 "unit": "multiples of modulator band edge k_p (k_p = 1.0)"},
        "rows": {
            "mean": {
                "label": "mean channel",
                "formal_support_extent": 1.0,   # modulator band B_p = [0, k_p]
                "formal_support_band": "B_p (modulator band)",
                "fisher_usable_extent": 1.0,    # full band usable; headroom closed at ridge
                "fisher_usable_note": "full band Fisher-usable; ridge optimum reached; residual headroom closed",
                "exact_wall": 1.0,              # first moment exactly blind outside B_p
                "wall_kind": "exact first-moment wall",
            },
            "covariance": {
                "label": "covariance channel",
                "formal_support_extent": 2.0,   # B_p (+) B_w = 2 k_p
                "formal_support_band": "B_p (+) B_w (pattern (+) medium band)",
                "fisher_usable_extent": 1.1,    # thin Fisher-usable rim ~1.1 k_p (P1)
                "fisher_usable_note": "thin Fisher-usable rim ~1.1 k_p; REGION_EMPTY at publication scale",
                "pocket_extent": 1.25,          # single natural-scene pocket
            },
        },
        "markers": [
            {"id": "M1_operating", "row": "mean", "x_kp": 1.0, "color": "green",
             "style": "dot",
             "label": "certified mean operating point",
             "source": "T_fig2_ridge_points.json.selected_operating_point + T_fig2_cohorts.json"},
            {"id": "killed_1p8", "row": "covariance", "x_kp": 1.8, "color": "red",
             "style": "cross",
             "label": "killed 1.8x covariance claim (P1 f_rec=0.05)",
             "source": "T_fig3_p1.json (witness f_rec_snr3=0.05 at claim 1.8)"},
            {"id": "pocket_1p25", "row": "covariance", "x_kp": 1.25, "color": "amber",
             "style": "filled_or_hollow",
             "label": "natural-scene 1.25x pocket (branch: PENDING)",
             "source": "T_fig3_pocket.json (f_rec_nat=0.79; branch PENDING -> A/B variant)"},
        ],
        "legend_line": "support is not supply",
        "derivation_notes": {
            "mean.formal_support_extent=1.0":
                "modulator band B_p; first moment exactly blind outside B_p (exact wall).",
            "mean.fisher_usable_extent=1.0":
                "M1 living optimum reached on the ridge (T_fig2_ridge_points.selected_operating_point); "
                "five preregistered interventions found no material residual headroom (T_fig2_closure).",
            "cov.formal_support_extent=2.0":
                "aperture law B_p (+) B_w = 2 k_p; validated 18.76x in/out separation at 64x64 "
                "(T_fig3_aperture.separation_64x64).",
            "cov.fisher_usable_extent=1.1":
                "P1 profiled-Fisher: only 5% of the 1.0-1.8 k_p annulus clears SNR>=3; usable aperture "
                "collapses to ~1.1 k_p (T_fig3_p1 + RISK_GATES_REPORT synthesis).",
            "cov.pocket_extent=1.25":
                "sole sub-P1 residue: natural-scene 1.25x pocket under k^-2 medium "
                "(T_fig3_pocket, f_rec_nat=0.79, NRMSE 0.10, T_req 197).",
            "marker.killed_1p8":
                "claim 1.8 k_p witness f_rec_snr3 = 0.05 << 0.70 bar (T_fig3_map / T_fig3_p1).",
        },
    }
    dump_json("T_fig1_atlas.json", obj)

# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("Extracting frozen source tables into paper_aperture/source_tables/ ...")
    t_fig2_ridge()
    t_fig2_cohorts()
    t_fig2_closure()
    t_fig3_aperture()
    t_fig3_p1()
    t_fig3_map()
    t_fig3_pocket()
    t_fig1_atlas()
    print("DONE.")
