#!/usr/bin/env python3
"""ROUND63 S2 lane planner — split the default campaign shards into five Colab
session plans (LOCAL, deterministic, read-only).

Reads ``results/round63/manifests/MANIFEST_INDEX.json`` and partitions the 82
``default_shards`` into five session plans:

    pro2 x3 sessions  ->  odd-parity shards   (account_hint == "pro2")
    pro1 x2 sessions  ->  even-parity shards   (account_hint == "pro1")

The frozen shard<->account parity (even index -> pro1, odd -> pro2) is already
baked into each shard's ``account_hint`` by make_manifests, so this planner
trusts ``account_hint`` as the parity assignment rather than re-deriving it.

Balancing within each account is greedy **LPT** (Longest Processing Time):
shards are sorted by descending ``est_hours`` and each is placed on the session
that currently has the least accumulated ``est_hours`` (ties broken by session
name, then shard_id). LPT is a 4/3-approximation to the min-makespan multiway
partition and is fully deterministic.

Outputs (under ``code/round63/colab/plans/`` by default):

    session_pro1a.txt  session_pro1b.txt
    session_pro2a.txt  session_pro2b.txt  session_pro2c.txt
    PLAN_SUMMARY.md

Each ``session_*.txt`` lists one shard_id per line (``#`` comments ignored by
session_driver.sh). The planner writes no state anywhere else and never touches
the manifests. Deterministic: same index in -> byte-identical plans out.
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
# code/round63/colab -> repo root is three levels up
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
DEFAULT_INDEX = os.path.join(ROOT, "results", "round63", "manifests",
                             "MANIFEST_INDEX.json")
DEFAULT_OUTDIR = os.path.join(HERE, "plans")

# Fixed session roster. Order is load-bearing for deterministic tie-breaks and
# for the summary table.
ACCOUNT_SESSIONS = {
    "pro1": ["pro1a", "pro1b"],       # even-parity shards, 2 sessions
    "pro2": ["pro2a", "pro2b", "pro2c"],  # odd-parity shards, 3 sessions
}
SESSION_ORDER = ["pro1a", "pro1b", "pro2a", "pro2b", "pro2c"]


def load_default_shards(index_path):
    with open(index_path) as f:
        idx = json.load(f)
    shards = idx.get("default_shards", [])
    if not shards:
        raise SystemExit("no default_shards in %s" % index_path)
    return idx, shards


def lpt_partition(shards, session_names):
    """Greedy LPT: assign shards (each a dict with est_hours, shard_id) to the
    given session_names, minimising the maximum per-session est_hours sum.

    Returns dict session_name -> list of shard dicts (in assignment order,
    i.e. descending est_hours)."""
    plan = {s: [] for s in session_names}
    load = {s: 0.0 for s in session_names}
    # descending est_hours; deterministic tie-break by shard_id
    ordered = sorted(shards, key=lambda c: (-float(c["est_hours"]),
                                            c["shard_id"]))
    for shard in ordered:
        # least-loaded session; tie-break by session_names order (stable min)
        target = min(session_names, key=lambda s: (load[s], session_names.index(s)))
        plan[target].append(shard)
        load[target] += float(shard["est_hours"])
    return plan


def build_plans(shards):
    """Partition all default shards into the five sessions. Returns
    (plan, unknown) where plan maps session_name -> [shard dict] and unknown is
    the list of shards whose account_hint is not pro1/pro2."""
    by_acct = {"pro1": [], "pro2": []}
    unknown = []
    for s in shards:
        acct = s.get("account_hint")
        if acct in by_acct:
            by_acct[acct].append(s)
        else:
            unknown.append(s)
    plan = {}
    for acct, sessions in ACCOUNT_SESSIONS.items():
        plan.update(lpt_partition(by_acct[acct], sessions))
    return plan, unknown, by_acct


def write_plan_files(plan, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    written = []
    for name in SESSION_ORDER:
        shards = plan[name]
        total = sum(float(s["est_hours"]) for s in shards)
        acct = "pro1" if name.startswith("pro1") else "pro2"
        path = os.path.join(out_dir, "session_%s.txt" % name)
        lines = [
            "# ROUND63 S2 lane plan -- session %s (account %s)" % (name, acct),
            "# %d shards, %.2f est_hours total (greedy-LPT balanced)" % (
                len(shards), total),
            "# one shard_id per line; blank lines and '#' comments ignored.",
        ]
        # shard ids sorted ascending for a stable, readable file
        for s in sorted(shards, key=lambda c: c["shard_id"]):
            lines.append(s["shard_id"])
        with open(path, "w", newline="\n") as f:
            f.write("\n".join(lines) + "\n")
        written.append(path)
    return written


def summary_rows(plan):
    rows = []
    for name in SESSION_ORDER:
        shards = plan[name]
        ests = [float(s["est_hours"]) for s in shards]
        acct = "pro1" if name.startswith("pro1") else "pro2"
        rows.append({
            "session": name,
            "account": acct,
            "n_shards": len(shards),
            "sum_hours": sum(ests),
            "max_shard": max(ests) if ests else 0.0,
            "min_shard": min(ests) if ests else 0.0,
        })
    return rows


def render_summary(rows, by_acct, n_total, out_dir, index_path):
    lines = []
    lines.append("# ROUND63 S2 lane plan summary")
    lines.append("")
    lines.append("Source index: `%s`" % os.path.relpath(index_path, ROOT))
    lines.append("")
    lines.append("Parity: pro2 x3 sessions = odd-parity shards "
                 "(account_hint=pro2); pro1 x2 = even (account_hint=pro1).")
    lines.append("Balance: greedy LPT on est_hours within each account.")
    lines.append("")
    lines.append("| session | account | shards | est_hours | max shard | min shard |")
    lines.append("|---------|---------|-------:|----------:|----------:|----------:|")
    for r in rows:
        lines.append("| %s | %s | %d | %.2f | %.3f | %.3f |" % (
            r["session"], r["account"], r["n_shards"], r["sum_hours"],
            r["max_shard"], r["min_shard"]))
    lines.append("")
    # per-account balance spread
    for acct in ("pro1", "pro2"):
        sums = [r["sum_hours"] for r in rows if r["account"] == acct]
        n = sum(r["n_shards"] for r in rows if r["account"] == acct)
        spread = (max(sums) - min(sums)) if sums else 0.0
        lines.append("- **%s**: %d shards across %d sessions, "
                     "%.2f est_hours total, per-session spread %.2f h "
                     "(max-min)." % (acct, n, len(sums), sum(sums), spread))
    lines.append("- **grand total**: %d shards, %.2f est_hours." % (
        n_total, sum(r["sum_hours"] for r in rows)))
    lines.append("")
    body = "\n".join(lines) + "\n"
    with open(os.path.join(out_dir, "PLAN_SUMMARY.md"), "w", newline="\n") as f:
        f.write(body)
    return body


def print_table(rows, by_acct, n_total):
    hdr = "%-8s %-7s %7s %11s %10s %10s" % (
        "session", "account", "shards", "est_hours", "max", "min")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print("%-8s %-7s %7d %11.2f %10.3f %10.3f" % (
            r["session"], r["account"], r["n_shards"], r["sum_hours"],
            r["max_shard"], r["min_shard"]))
    print("-" * len(hdr))
    tot_shards = sum(r["n_shards"] for r in rows)
    tot_hours = sum(r["sum_hours"] for r in rows)
    print("%-8s %-7s %7d %11.2f" % ("TOTAL", "", tot_shards, tot_hours))
    print()
    print("parity split: pro1(even)=%d shards, pro2(odd)=%d shards, total=%d"
          % (len(by_acct["pro1"]), len(by_acct["pro2"]), n_total))


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Split ROUND63 default shards into 5 Colab session plans "
                    "(deterministic, greedy-LPT balanced).")
    ap.add_argument("--index", default=DEFAULT_INDEX,
                    help="path to MANIFEST_INDEX.json (default: %(default)s)")
    ap.add_argument("--out-dir", default=DEFAULT_OUTDIR,
                    help="directory for session_*.txt + PLAN_SUMMARY.md "
                         "(default: %(default)s)")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the split table only; write no files")
    args = ap.parse_args(argv)

    idx, shards = load_default_shards(args.index)
    plan, unknown, by_acct = build_plans(shards)
    n_total = len(shards)
    if unknown:
        print("WARNING: %d shard(s) with unrecognised account_hint (skipped): %s"
              % (len(unknown), ", ".join(s["shard_id"] for s in unknown)),
              file=sys.stderr)

    rows = summary_rows(plan)
    print_table(rows, by_acct, n_total)

    assigned = sum(r["n_shards"] for r in rows)
    if assigned != n_total - len(unknown):
        raise SystemExit("BUG: assigned %d != %d expected"
                         % (assigned, n_total - len(unknown)))

    if args.dry_run:
        print("\n[dry-run] no files written.")
        return 0

    written = write_plan_files(plan, args.out_dir)
    render_summary(rows, by_acct, n_total, args.out_dir, args.index)
    print("\nwrote %d session plans + PLAN_SUMMARY.md to %s"
          % (len(written), args.out_dir))
    for p in written:
        print("  %s" % p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
