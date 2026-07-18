"""F1 freeze ledger builder (spec D2.3 §1 SHA hard gate; round-7 §7 locks).

Collects into results/round63/F1_LEDGER.json:
  - git commit + dirty flag
  - python / numpy / scipy versions + full pip freeze (hashed, stored aside)
  - sha256 of every frozen input surface: image caches (confirmatory natural,
    detail24 confirmatory + dev, S1 dev), C0_FROZEN.json, detail24 params
    manifests, spec + ruling docs, and the frozen figure preregistration
  - manifests directory aggregate hash (if generated)
Rerunnable: output is deterministic given the same tree (no wall clock in the
hashed payload; generated_utc kept OUTSIDE the self-hash).
"""
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT_DIR = os.path.join(REPO, "results", "round63")


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_tree(root, exts=None):
    """Aggregate sha256 over a directory tree (sorted relative paths)."""
    h = hashlib.sha256()
    for dirpath, dirnames, filenames in sorted(os.walk(root)):
        dirnames.sort()
        for fn in sorted(filenames):
            if exts and not fn.lower().endswith(exts):
                continue
            p = os.path.join(dirpath, fn)
            rel = os.path.relpath(p, root).replace("\\", "/")
            h.update(rel.encode())
            h.update(sha256_file(p).encode())
    return h.hexdigest()


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    git = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO,
                         capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain"], cwd=REPO,
                           capture_output=True, text=True).stdout.strip()
    import numpy, scipy
    freeze = subprocess.run([sys.executable, "-m", "pip", "freeze"],
                            capture_output=True, text=True).stdout
    freeze_path = os.path.join(OUT_DIR, "env_lock.txt")
    with open(freeze_path, "w") as f:
        f.write(freeze)

    inputs = {}
    surfaces = {
        "images_confirmatory_natural": os.path.join(REPO, "data",
                                                    "r63_images", "64"),
        "images_detail24_confirmatory": os.path.join(
            REPO, "data", "r63_images_detail24", "64"),
        "images_detail24_dev": os.path.join(
            REPO, "data", "r63_images_detail24_dev", "64"),
        "images_s1_dev": os.path.join(REPO, "data", "r63_images_dev", "64"),
    }
    for k, root in surfaces.items():
        inputs[k] = sha256_tree(root) if os.path.isdir(root) else None
    for k, p in {
        "C0_FROZEN": os.path.join(HERE, "C0_FROZEN.json"),
        "detail24_params_manifest": os.path.join(
            REPO, "data", "r63_images_detail24", "64", "params_manifest.json"),
        "spec_D2": os.path.join(REPO, "docs", "ROUND63_SPEC_D2.md"),
        "ruling_r6": os.path.join(REPO, "docs", "ROUND63_GPT_ROUND6_RULING.md"),
        "ruling_r7": os.path.join(REPO, "docs", "ROUND63_GPT_ROUND7_RULING.md"),
        "figure_prereg": os.path.join(REPO, "docs",
                                      "ROUND63_FIGURE_PREREG.md"),
        "env_lock": freeze_path,
    }.items():
        inputs[k] = sha256_file(p) if os.path.exists(p) else None
    man_dir = os.path.join(OUT_DIR, "manifests")
    inputs["manifests_tree"] = (sha256_tree(man_dir)
                                if os.path.isdir(man_dir) else None)
    exp = os.path.join(OUT_DIR, "expected_cells.csv")
    inputs["expected_cells"] = sha256_file(exp) if os.path.exists(exp) else None

    payload = {
        "git_commit": git,
        "git_dirty": bool(dirty),
        "python": sys.version.split()[0],
        "numpy": numpy.__version__,
        "scipy": scipy.__version__,
        "inputs_sha256": inputs,
    }
    self_hash = hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode()).hexdigest()
    ledger = {"payload": payload, "payload_sha256": self_hash,
              "generated_utc": datetime.now(timezone.utc).isoformat()}
    out = os.path.join(OUT_DIR, "F1_LEDGER.json")
    with open(out, "w") as f:
        json.dump(ledger, f, indent=1)
    print("F1 ledger ->", out)
    print("payload_sha256:", self_hash)
    print("git:", git, "(dirty)" if dirty else "(clean)")
    missing = [k for k, v in inputs.items() if v is None]
    print("missing surfaces:", missing or "none")
    return 0


if __name__ == "__main__":
    sys.exit(main())
