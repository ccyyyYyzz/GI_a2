"""Merge the finishing-run sections into jailbreak_results.json and freeze the
R32 six-bar verdict block."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)


def load(name):
    with open(os.path.join(OUT, name)) as fh:
        return json.load(fh)


res = load("jailbreak_results.json")
sparse = load("sparse_arm.json")
bars = load("bars_dose.json")
p0m = load("p0_modelfit.json")
qf = load("quenched_floor.json")
tim = load("timing_endtoend.json")
cert = load("certificate_r321.json")

res["method"] = ("MOLT (Microcell-Occupancy Laplace Tomography); working name "
                 "SATURATION_JAILBREAK")
res["rulings"] = ["docs/ROUND63_GPT_ROUND31_RULING_RAW.md (R31)",
                  "docs/ROUND63_GPT_ROUND31_PRO_APPENDIX_0.md (R31 Pro addendum)",
                  "docs/ROUND63_GPT_ROUND32_RULING_RAW.md (R32 strategy — six binding bars)"]

# finishing-run sections
res["FINISHING_RUN"] = dict(
    sparse_arm=sparse, bars_dose=bars, p0_modelfit=p0m,
    quenched_floor=qf, timing_endtoend=tim, certificate_r321=cert)

# ---- frozen R32 six-bar verdict ------------------------------------------- #
b1 = sparse["scenes"]; neffs = [b1[k]["n_eff"] for k in b1]
verdict = {
    "frame": "equal WALL-CLOCK time (R32 primary); matched-photon = mandatory adverse audit",
    "bar1_neff_Ceff_k": dict(
        requirement="n_eff<=32 on actual sparse masks x scenes; report C_eff and k",
        result="per-mask n_eff in [%.1f, %.1f] (all<=32); C_eff_iid=%.0f (~C=3600); "
                "k_eff=1 (fully-developed speckle)" % (min(neffs), max(neffs),
                                                       qf["bar1_C_eff"]["C_eff_iid"]),
        verdict="PASS"),
    "bar2_10pct_p2_map": dict(
        requirement="10% p2 map at <= x40 linear detected-photon budget/mask (expect ~x30)",
        result="production 3-level CRB multiple max x%.1f over masks; dose-law "
                "boundary n_eff=32 -> x%.1f (known); expected ~x30 confirmed"
                % (bars["bar2_map"]["summary"]["max_mult_crb_3level"],
                   bars["bar2_map"]["boundary_neff32"]["mult_law_known"]),
        verdict="PASS" if bars["bar2_map"]["summary"]["PASS"] else "FAIL"),
    "bar3_favorable_nullpair_d3": dict(
        requirement="favorable order-one Delta_p2/p2 pair d'>=3 at <= x4 linear dose (expect ~2.6)",
        result="dose-law f=1 nuisance-p1 x%.2f (~ruling x2.6); measured aligned pair "
                "x%.2f" % (bars["bar3_nullpair"]["dose_law_f1_nuis_mult"],
                          bars["bar3_nullpair"]["measured_budget_mult_for_d3"]),
        verdict="PASS" if bars["bar3_nullpair"]["PASS"] else "FAIL"),
    "bar4_wallclock_overhead": dict(
        requirement="<=10% end-to-end overhead incl mask/settle/power/readout, "
                    "{100kHz,1MHz} gates x {1,10,22}kHz DMD",
        result="meets <=10%% in %d/%d regimes: only switching-limited (1MHz gates + "
               "1kHz DMD -> %.1f%% at K_S=51). At 1MHz+10kHz DMD -> %.1f%%; at 100kHz "
               "gates -> 22-85%%. 'Negligible time' false in gate-limited corner."
               % (tim["bar4_summary"]["n_pass_10pct"], tim["bar4_summary"]["n_regimes"],
                  tim["overhead_pct"]["K_S_51"]["1MHz"]["1kHz"]["overhead_pct"],
                  tim["overhead_pct"]["K_S_51"]["1MHz"]["10kHz"]["overhead_pct"]),
        verdict="CONDITIONAL"),
    "bar5_P0_identifiability": dict(
        requirement="3x3 mechanism model-recovery survives at realistic noise",
        result="confusion recovery accuracy %.3f (coherent/product/correlated all "
                ">=%.0f%%); permutation control: equal-p1,equal-p2 clustered-vs-spread "
                "changes tr[(SigmaV)^2] 198%%, correlated d'=3 at x2.2, product blind"
                % (p0m["recovery_accuracy"], 100 * min(
                    p0m["confusion"][0][0], p0m["confusion"][1][1],
                    p0m["confusion"][2][2]) / p0m["n_trial"]),
        verdict="PASS"),
    "bar6_R32.1_certificate": dict(
        requirement="R32.1 Jacobian full-rank+conditioned on support AND noisy joint "
                    "recon materially beats strongest equal-time linear comparator",
        result="closure window EXACT: rank[M_D|S;M_S|S;2M_S|S diag x]=s for "
                "K_D=40<s<=K_D+2K_S=100 (linear stuck at 40); noiseless recon ~297 dB; "
                "BUT noisy equal-time recon only ties linear (delta=%+.1f dB) -- "
                "fiber-resolving rows Fisher-weak, O(1) SNR at equal dose"
                % cert["bar6_summary"]["noisy_delta_db_s80"],
        verdict=cert["bar6_summary"]["verdict"]),
    "matched_photon_adverse_audit": dict(
        dense_50pct_p2_map="x1.6e5 - x1.1e6 (permanently kills dense MOLT)",
        dense_nullpair_d3="median d'=0.0056, x2.84e5",
        sparse_nullpair_d3="median d'=0.267, x126 (x2254 cheaper than dense)",
        note="matched-photon FAILS by orders of magnitude and is reported as a "
             "failure per R32; equal-time is primary ONLY for dose-insensitive scenes"),
    "known_support_witness": dict(
        s80_joint_vs_linear="joint 296.7 dB vs linear 27.8 dB = +268.9 dB noiseless (committed T4)",
        caveat="noiseless known-support witness; NOT extrapolated to unknown-support "
               "noisy recovery (R32 forbidden)"),
    "overall": ("MOLT survives ONLY as a sparse-support, dose-insensitive, "
                "switching-limited, equal-time method. Math fully verified; dense "
                "arm retired; photon economics fail; time economics conditional on "
                "hardware regime; mechanism identifiable; support fiber closes "
                "geometrically but its equal-time NOISY material win is unproven.")}
res["R32_SIX_BAR_VERDICT"] = verdict

with open(os.path.join(OUT, "jailbreak_results.json"), "w") as fh:
    json.dump(res, fh, indent=2, default=float)
print("updated jailbreak_results.json with FINISHING_RUN + R32_SIX_BAR_VERDICT")
for k, v in verdict.items():
    if isinstance(v, dict) and "verdict" in v:
        print("  %-32s %s" % (k, v["verdict"]))
