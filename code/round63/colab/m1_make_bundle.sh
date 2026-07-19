#!/usr/bin/env bash
# ROUND63 **M1** bundle assembler -- builds ONE self-contained tarball that
# carries everything a Colab VM needs to run any M1 shard. Modeled on
# make_bundle.sh (S2), retargeted at the M1 code/data/manifests.
#
# The staged tree reproduces the repo-root layout EXACTLY, so shard_runner.py's
# ROOT = dirname(dirname(HERE)) resolves output_csv / meta_json / frozen_inputs;
# detail24.DATA = ROOT/data resolves the m1 image caches; and oed_design_v2's
# J_NPZ/J_CSV = ROOT/results/round63_theory/... resolve unchanged in the bundle.
# It unpacks into a single "round63_m1_bundle/" dir -> launch untars at /content
# so the bundle root is /content/round63_m1_bundle.
#
# Deterministic: fixed order (tar --sort=name), zeroed mtimes/owners, gzip -n.
# A MANIFEST_SHA256 of every bundled file is written. >50 MB is split into 45 MB
# chunks (reassemble: cat ...part_* > ....tar.gz).
#
# Contents (task list + load-bearing deps flagged [dep]):
#   code/round63/*.py                    (all M1 + shared modules)
#   code/round63/C0_FROZEN.json          [dep] select_eta (imported by campaign.run_cell)
#   code/round63/colab/session_driver.sh + remote_lane.py   (on-VM drivers)
#   code/gi_core/*.py                    [dep] campaign/m1_runner: `from gi_core ...`
#   data/r63_images_m1/     + data/r63_images_m1_dev/        (m1 scene caches, side 32)
#   results/round63_theory/fig_a_sweep.npz + fisher_ridge.csv (J tables for oed_design_v2+)
#   results/round63_m1/manifests/                            (40 shard manifests + index)
#   results/round63_m1/designs/          [dep] OED/MATCH1/RIDGE frozen_inputs (*.npz)
#                                        NOTE: incomplete until the m1-freeze fix
#                                        builds the OED/RIDGE design caches.
#
# NOT included (verified unnecessary for the "m1" imageset): detail24/detail32
# caches. m1_scenes builds its own scenes via detail24._build32 into the m1 cache
# roots; FAMILIES/PARAMS/PARAMS32 are code-level constants (no detail32 param
# manifest is read), and campaign.load_detail32_rois fires only for detail32*
# imagesets -- M1 rows carry blank CNR (m1_runner A10).
#
# Usage:
#   bash m1_make_bundle.sh [--repo DIR] [--out DIR] [--name NAME] [--dry-run]
#     --repo DIR    repo root        (default: /mnt/d/GI_another)
#     --out  DIR    build output dir (default: /mnt/d/tmp/round63_m1_build)
#     --name NAME   bundle basename  (default: round63_m1_bundle)
#     --dry-run     list every file that WOULD be packed (no staging, no tar)
set -euo pipefail

REPO="/mnt/d/GI_another"
OUT="/mnt/d/tmp/round63_m1_build"
NAME="round63_m1_bundle"
DRYRUN=0
SPLIT_THRESHOLD=$((50 * 1024 * 1024))   # 50 MB
SPLIT_CHUNK="45M"

while [ $# -gt 0 ]; do
  case "$1" in
    --repo) REPO="$2"; shift 2 ;;
    --out)  OUT="$2";  shift 2 ;;
    --name) NAME="$2"; shift 2 ;;
    --dry-run) DRYRUN=1; shift ;;
    -h|--help) sed -n '2,45p' "$0"; exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done

command -v tar >/dev/null       || { echo "need tar" >&2; exit 1; }
command -v gzip >/dev/null      || { echo "need gzip" >&2; exit 1; }
command -v sha256sum >/dev/null || { echo "need sha256sum" >&2; exit 1; }
command -v split >/dev/null     || { echo "need split" >&2; exit 1; }
[ -d "$REPO" ] || { echo "repo not found: $REPO" >&2; exit 1; }
[ -d "$REPO/code/round63" ] || { echo "code/round63 not found under $REPO" >&2; exit 1; }

STAGE="$OUT/$NAME"
NPLAN=0
BYTES=0
MISSING=0

# ---- packing primitives: mirror repo-relative paths exactly ----------------- #
pack_one() {   # $1 = abs src, $2 = repo-relative dest path
  local src="$1" rel="$2" sz
  sz="$(stat -c%s "$src")"
  if [ "$DRYRUN" = 1 ]; then
    printf '  PACK  %-58s %10s B\n' "$rel" "$sz"
  else
    mkdir -p "$STAGE/$(dirname "$rel")"
    cp "$src" "$STAGE/$rel"
  fi
  NPLAN=$((NPLAN + 1)); BYTES=$((BYTES + sz))
}

pack_file() {  # $1 = repo-relative file (must exist)
  local rel="$1" src="$REPO/$1"
  if [ -e "$src" ]; then pack_one "$src" "$rel"
  else echo "  MISSING (file): $rel" >&2; MISSING=$((MISSING + 1)); fi
}

pack_glob() {  # $1 = repo-relative shell glob (top-level, non-recursive)
  local g="$1" f rel n=0
  for f in $REPO/$g; do
    [ -e "$f" ] || continue
    rel="${f#"$REPO"/}"; pack_one "$f" "$rel"; n=$((n + 1))
  done
  [ "$n" = 0 ] && { echo "  MISSING (glob): $g" >&2; MISSING=$((MISSING + 1)); }
  return 0
}

pack_tree() {  # $1 = repo-relative dir (recursive; all files, deterministic order)
  local rel="$1" src="$REPO/$1" f frel n=0
  [ -d "$src" ] || { echo "  MISSING (dir): $rel" >&2; MISSING=$((MISSING + 1)); return 0; }
  while IFS= read -r -d '' f; do
    frel="${f#"$REPO"/}"; pack_one "$f" "$frel"; n=$((n + 1))
  done < <(find "$src" -type f ! -path '*/__pycache__/*' -print0 | LC_ALL=C sort -z)
  [ "$n" = 0 ] && echo "  NOTE: $rel is empty" >&2
  return 0
}

echo "[m1_make_bundle] repo=$REPO  name=$NAME  dry_run=$DRYRUN" >&2
[ "$DRYRUN" = 1 ] || { rm -rf "$STAGE"; mkdir -p "$STAGE"; }
echo "[m1_make_bundle] --- packing plan ---" >&2

# (1) round63 python modules (all: campaign, shard_runner, m1_*, oed_*, select_eta, ...)
pack_glob "code/round63/*.py"
# (2) frozen C0 [dep]: select_eta.C0_FILE = HERE/C0_FROZEN.json
pack_file "code/round63/C0_FROZEN.json"
# (3) on-VM drivers
pack_file "code/round63/colab/session_driver.sh"
pack_file "code/round63/colab/remote_lane.py"
# (4) gi_core package [dep]: campaign/m1_runner import gi_core.{metrics,utils}
pack_glob "code/gi_core/*.py"
# (5) m1 scene caches (side 32): frozen_inputs of the scatter/lblob shards + m1 imageset
pack_tree "data/r63_images_m1"
pack_tree "data/r63_images_m1_dev"
# (6) theory J tables: oed_design_v2 J_NPZ / J_CSV
pack_file "results/round63_theory/fig_a_sweep.npz"
pack_file "results/round63_theory/fisher_ridge.csv"
# (7) shard manifests + index (40 shards)
pack_glob "results/round63_m1/manifests/*.json"
pack_glob "results/round63_m1/manifests/*.md"
# (8) design caches [dep]: OED/MATCH1/RIDGE frozen_inputs (*.npz) -- incomplete pre-freeze
pack_tree "results/round63_m1/designs"

echo "[m1_make_bundle] planned files=$NPLAN  bytes=$BYTES  missing=$MISSING" >&2

if [ "$DRYRUN" = 1 ]; then
  MB="$(awk "BEGIN{printf \"%.2f\", $BYTES/1048576}")"
  echo ""
  echo "============= M1 BUNDLE DRY-RUN ============="
  echo " would pack : $NPLAN files"
  echo " total size : $BYTES bytes ($MB MB)"
  echo " missing    : $MISSING"
  [ "$MISSING" -gt 0 ] && echo " (missing entries above are expected pre-freeze, e.g. OED/RIDGE design caches)"
  echo " (dry-run: nothing staged, no tarball built)"
  echo "============================================"
  exit 0
fi

# ---- normalise perms for reproducibility ------------------------------------ #
find "$STAGE" -type d -exec chmod 755 {} +
find "$STAGE" -type f -exec chmod 644 {} +
[ -f "$STAGE/code/round63/colab/session_driver.sh" ] && \
  chmod 755 "$STAGE/code/round63/colab/session_driver.sh"

# ---- MANIFEST_SHA256 of every bundled file ---------------------------------- #
( cd "$STAGE" && find . -type f ! -name MANIFEST_SHA256 -print0 \
    | LC_ALL=C sort -z | xargs -0 sha256sum | sed 's#  \./#  #' > MANIFEST_SHA256 )
NFILES="$(wc -l < "$STAGE/MANIFEST_SHA256")"

# ---- deterministic tar + gzip -n -------------------------------------------- #
TAR="$OUT/$NAME.tar"; GZ="$TAR.gz"
rm -f "$TAR" "$GZ" "$GZ".part_* 2>/dev/null || true
tar --sort=name --format=gnu --mtime='@0' --owner=0 --group=0 --numeric-owner \
    -C "$OUT" -cf "$TAR" "$NAME"
gzip -n -9 "$TAR"      # -> $GZ

SZ="$(stat -c%s "$GZ")"
SHA="$(sha256sum "$GZ" | awk '{print $1}')"
SZ_MB="$(awk "BEGIN{printf \"%.2f\", $SZ/1048576}")"
echo ""
echo "=============== ROUND63 M1 BUNDLE ==============="
echo " bundle : $GZ"
echo " files  : $NFILES"
echo " size   : $SZ bytes ($SZ_MB MB)"
echo " sha256 : $SHA"

if [ "$SZ" -gt "$SPLIT_THRESHOLD" ]; then
  split -b "$SPLIT_CHUNK" -d -a 3 "$GZ" "$GZ.part_"
  echo " split  : >50 MB -> $SPLIT_CHUNK chunks:"
  for p in "$GZ".part_*; do
    echo "            $(basename "$p")  $(stat -c%s "$p") bytes  $(sha256sum "$p" | awk '{print $1}')"
  done
  echo " reassemble on the VM:  cat $NAME.tar.gz.part_* > $NAME.tar.gz"
else
  echo " split  : not needed (<=50 MB) -- upload $NAME.tar.gz whole"
fi
echo "================================================"
echo "verify:  echo \"$SHA  $NAME.tar.gz\" | sha256sum -c -"
echo "unpack:  tar -xzf $NAME.tar.gz   # -> $NAME/ (repo-root layout)"
