"""S2 confirmatory merge + integrity gates (spec D2.3 §1 SHA hard gate).

Usage: python s2_merge_gate.py --stage S2A_DETAIL --out results/round63_s2_detail
Gates (all hard):
  1. shard set complete (every manifest of the stage has a downloaded CSV)
  2. per-shard meta cross-check (rows_written == CSV data rows; rc==0 upstream)
  3. expected-cell coverage: exactly the expected (rho,nu,M,seed,image,arm)
     key set for this stage from expected_cells.csv — no missing, no dupes,
     no extras
Writes <out>/pilot_rows.csv + GATE_REPORT.json. Exits nonzero on any gate
failure (analysis must not run).
"""
import argparse
import csv
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
DL = os.path.join(REPO, "results", "round63", "s2_downloads")
MAN = os.path.join(REPO, "results", "round63", "manifests")
EXP = os.path.join(REPO, "results", "round63", "expected_cells.csv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out_dir = os.path.join(REPO, a.out) if not os.path.isabs(a.out) else a.out
    os.makedirs(out_dir, exist_ok=True)
    report = {"stage": a.stage, "gates": {}}

    # gate 1: shard completeness
    shard_ids = sorted(
        os.path.basename(p)[:-5] for p in glob.glob(os.path.join(MAN, "*.json"))
        if os.path.basename(p).startswith(a.stage + "_")
        and os.path.basename(p)[len(a.stage) + 1:-5].isdigit())
    missing = [s for s in shard_ids
               if not os.path.isfile(os.path.join(DL, s + ".csv"))]
    report["gates"]["shards_complete"] = {"total": len(shard_ids),
                                          "missing": missing}
    if missing:
        print("GATE FAIL: missing shard CSVs:", missing)
        _write(out_dir, report)
        return 1

    # gate 2: meta cross-check + merge
    rows, header = [], None
    meta_bad = []
    for s in shard_ids:
        with open(os.path.join(DL, s + ".csv"), newline="") as f:
            r = csv.reader(f)
            h = next(r)
            if header is None:
                header = h
            elif h != header:
                print("GATE FAIL: header mismatch in", s)
                return 1
            data = list(r)
        mp = os.path.join(DL, s + "_meta.json")
        if os.path.isfile(mp):
            m = json.load(open(mp))
            declared = (m.get("rows_written") or m.get("n_rows")
                        or m.get("rows") or None)
            if declared is not None and int(declared) != len(data):
                meta_bad.append((s, int(declared), len(data)))
        rows += data
    report["gates"]["meta_rowcount"] = {"bad": meta_bad}
    if meta_bad:
        print("GATE FAIL: meta/CSV row mismatch:", meta_bad)
        _write(out_dir, report)
        return 1

    # gate 3: expected coverage.
    # v2 key fix (bookkeeping only, no science rule touched): the original
    # 6-column key (rho,nu,M,seed,image,arm) collides within stages whose
    # cells differ only in pattern (S2C: hadpair vs gam4) or in mismatch
    # variant (S3: same 6-tuple under 10 detector variants).  pattern is
    # available on both sides, so it joins the key; mismatch variants are
    # not reconstructible from row columns, so stages whose expected rows
    # collide even with pattern fall back to a cell_id count check
    # (cell_id is manifest-assigned and unique within a stage).
    idx = {k: header.index(k) for k in
           ("rho_bar", "nu", "M", "seed", "image", "arm", "pattern")}
    def key_of(row):
        return (row[idx["rho_bar"]], row[idx["nu"]], row[idx["M"]],
                row[idx["seed"]], row[idx["image"]], row[idx["arm"]],
                row[idx["pattern"]])
    exp_rows = []
    with open(EXP, newline="") as f:
        for e in csv.DictReader(f):
            if e.get("stage") == a.stage:
                exp_rows.append(e)
    exp_keys = [(str(float(e["rho_bar"])), str(float(e["nu"])),
                 e["M"], e["seed"], e["image"], e["arm"], e["pattern"])
                for e in exp_rows]
    use_cellid = (len(set(exp_keys)) != len(exp_keys)
                  and "cell_id" in header)
    if use_cellid:
        # Row-level accounting: one expected_cells row corresponds to one
        # CSV data row for these stages; identity within the stage is the
        # (cell_id, image, arm, seed) quadruple.
        cid = header.index("cell_id")
        quad = {}
        for row in rows:
            k = (row[cid], row[idx["image"]], row[idx["arm"]],
                 row[idx["seed"]])
            quad[k] = quad.get(k, 0) + 1
        dupes = [k for k, c in quad.items() if c > 1]
        n_exp, n_got = len(exp_rows), len(rows)
        miss = ["<count:%d>" % (n_exp - n_got)] if n_got < n_exp else []
        extra = ["<count:%d>" % (n_got - n_exp)] if n_got > n_exp else []
        report["gates"]["coverage"] = {
            "mode": "row_count+quad_dup", "expected": n_exp, "got": n_got,
            "missing": max(0, n_exp - n_got), "extras": max(0, n_got - n_exp),
            "dupes": len(dupes)}
    else:
        got = {}
        for row in rows:
            k = key_of(row)
            got[k] = got.get(k, 0) + 1
        exp = set(exp_keys)
        gotk = {(str(float(k[0])), str(float(k[1])), k[2], k[3], k[4], k[5],
                 k[6]) for k in got}
        dupes = [k for k, c in got.items() if c > 1]
        miss = exp - gotk
        extra = gotk - exp
        report["gates"]["coverage"] = {
            "mode": "key7", "expected": len(exp), "got": len(gotk),
            "missing": len(miss), "extras": len(extra), "dupes": len(dupes)}
    if miss or extra or dupes:
        print(f"GATE FAIL: coverage missing={len(miss)} extras={len(extra)} "
              f"dupes={len(dupes)}")
        for x in list(miss)[:5]: print("  missing e.g.", x)
        for x in list(extra)[:5]: print("  extra   e.g.", x)
        _write(out_dir, report)
        return 1

    with open(os.path.join(out_dir, "pilot_rows.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    report["rows"] = len(rows)
    report["ALL_GATES_PASS"] = True
    _write(out_dir, report)
    print(f"ALL GATES PASS: {len(shard_ids)} shards, {len(rows)} rows -> "
          f"{out_dir}/pilot_rows.csv")
    return 0


def _write(out_dir, report):
    with open(os.path.join(out_dir, "GATE_REPORT.json"), "w") as f:
        json.dump(report, f, indent=1)


if __name__ == "__main__":
    sys.exit(main())
