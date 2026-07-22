"""Post-grid gate evaluation runner (frozen record of the 2026-07-22 verdict).

Reads all 1920 cell records, applies the schema-flatten adapter (the harness
writes A_realized / route_latency_s / w_realized nested under arm_info; the
frozen gates read them at top level -- pure field promotion, same values,
no science change), runs bridge_gates.evaluate_all, writes GATE_VERDICT.json.
"""
import json, glob, sys, os

sys.path.insert(0, os.path.dirname(__file__))
import bridge_gates as bg


def load_flattened(cell_glob="results/round63_bridge/cells/*.json"):
    records = []
    for f in glob.glob(cell_glob):
        r = json.load(open(f))
        ai = r.get("arm_info") or {}
        rd = r.get("rlmi_disclosures") or {}
        if "A_realized" not in r and "A_realized" in ai:
            r["A_realized"] = ai["A_realized"]
        if "route_latency_s" not in r and "route_latency_s" in ai:
            r["route_latency_s"] = ai["route_latency_s"]
        if "w0_realized" not in r:
            w = ai.get("w_realized") or rd.get("w_realized")
            if w is not None:
                r["w0_realized"] = float(w[0])
        records.append(r)
    return records


if __name__ == "__main__":
    records = load_flattened()
    out = bg.evaluate_all(records)
    with open("results/round63_bridge/GATE_VERDICT.json", "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(out["verdict"])
