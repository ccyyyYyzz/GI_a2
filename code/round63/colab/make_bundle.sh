#!/usr/bin/env bash
# ROUND63 S2 bundle assembler -- builds ONE self-contained tarball that carries
# everything a Colab VM needs to run any default shard: the round63 code, the
# gi_core package, the frozen image caches, the frozen C0, the shard manifests,
# and the on-VM drivers (remote_lane.py / session_driver.sh) + lane plans.
#
# The staged tree reproduces the repo-root layout exactly, so shard_runner.py's
# ROOT = dirname(dirname(HERE)) resolves output_csv / frozen_inputs / data paths
# and select_eta's HERE/C0_FROZEN.json all resolve unchanged inside the bundle.
#
# Deterministic: fixed file order (tar --sort=name), zeroed mtimes/owners, and
# gzip -n (no name/timestamp in the gzip header) -> byte-reproducible tarball.
# A MANIFEST_SHA256 of every bundled file is included. If the gzip exceeds
# 50 MB it is split into <=45 MB chunks (reassemble: cat ...part_* > ....tar.gz).
#
# Usage:
#   bash make_bundle.sh [--repo DIR] [--out DIR] [--name NAME]
#     --repo DIR   repo root       (default: /mnt/d/GI_another)
#     --out  DIR   build output dir (default: /mnt/d/tmp/round63_s2_build)
#     --name NAME  bundle basename  (default: round63_s2_bundle)
set -euo pipefail

REPO="/mnt/d/GI_another"
OUT="/mnt/d/tmp/round63_s2_build"
NAME="round63_s2_bundle"
SPLIT_THRESHOLD=$((50 * 1024 * 1024))   # 50 MB
SPLIT_CHUNK="45M"

while [ $# -gt 0 ]; do
  case "$1" in
    --repo) REPO="$2"; shift 2 ;;
    --out)  OUT="$2";  shift 2 ;;
    --name) NAME="$2"; shift 2 ;;
    -h|--help)
      sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done

command -v tar >/dev/null       || { echo "need tar" >&2; exit 1; }
command -v gzip >/dev/null      || { echo "need gzip" >&2; exit 1; }
command -v sha256sum >/dev/null || { echo "need sha256sum" >&2; exit 1; }
command -v split >/dev/null     || { echo "need split" >&2; exit 1; }

[ -d "$REPO" ] || { echo "repo not found: $REPO" >&2; exit 1; }
CODE="$REPO/code/round63"
[ -d "$CODE" ] || { echo "code/round63 not found under $REPO" >&2; exit 1; }

STAGE="$OUT/$NAME"
rm -rf "$STAGE"
mkdir -p "$STAGE/code/round63/colab" "$STAGE/code/gi_core" \
         "$STAGE/data" "$STAGE/results/round63/manifests"

echo "[make_bundle] staging from $REPO -> $STAGE" >&2

# ---- (1) round63 top-level python modules (incl shard_runner.py, campaign.py)  #
for f in "$CODE"/*.py; do
  cp "$f" "$STAGE/code/round63/"
done

# ---- (2) frozen C0 (select_eta expects it at code/round63/C0_FROZEN.json) ----- #
[ -f "$CODE/C0_FROZEN.json" ] || { echo "missing C0_FROZEN.json" >&2; exit 1; }
cp "$CODE/C0_FROZEN.json" "$STAGE/code/round63/C0_FROZEN.json"

# ---- (3) on-VM drivers + lane plans (bundle mirrors code/round63/colab) ------- #
cp "$CODE/colab/remote_lane.py"    "$STAGE/code/round63/colab/"
cp "$CODE/colab/session_driver.sh" "$STAGE/code/round63/colab/"
[ -f "$CODE/colab/LAUNCH_RUNBOOK.md" ] && \
  cp "$CODE/colab/LAUNCH_RUNBOOK.md" "$STAGE/code/round63/colab/"
if [ -d "$CODE/colab/plans" ]; then
  mkdir -p "$STAGE/code/round63/colab/plans"
  # only the plan text/markdown, no stray files
  for p in "$CODE"/colab/plans/*.txt "$CODE"/colab/plans/*.md; do
    [ -e "$p" ] && cp "$p" "$STAGE/code/round63/colab/plans/"
  done
else
  echo "[make_bundle] NOTE: no colab/plans/ yet -- run lane_plan.py first to "\
       "embed session plans (bundle is still valid; pass a plan path at launch)." >&2
fi

# ---- (4) gi_core package (python only) ---------------------------------------- #
for f in "$REPO"/code/gi_core/*.py; do
  cp "$f" "$STAGE/code/gi_core/"
done

# ---- (5) frozen image caches (three side-64 sets the default shards consume) -- #
copy_cache() {   # $1 = repo-relative cache dir
  local src="$REPO/$1"
  [ -d "$src" ] || { echo "missing data cache: $1" >&2; exit 1; }
  mkdir -p "$STAGE/$1"
  # PNGs + sha256.txt + params_manifest.json; skip nothing, no subdirs expected
  find "$src" -maxdepth 1 -type f -print0 | while IFS= read -r -d '' f; do
    cp "$f" "$STAGE/$1/"
  done
}
copy_cache "data/r63_images/64"
copy_cache "data/r63_images_detail24/64"
copy_cache "data/r63_images_detail24_dev/64"

# ---- (6) shard manifests + index ---------------------------------------------- #
for f in "$REPO"/results/round63/manifests/*.json "$REPO"/results/round63/manifests/*.md; do
  [ -e "$f" ] && cp "$f" "$STAGE/results/round63/manifests/"
done

# ---- normalise perms for reproducibility -------------------------------------- #
find "$STAGE" -type d -exec chmod 755 {} +
find "$STAGE" -type f -exec chmod 644 {} +
chmod 755 "$STAGE/code/round63/colab/session_driver.sh"

# ---- MANIFEST_SHA256 of every bundled file (sorted, bundle-relative) ---------- #
( cd "$STAGE" && find . -type f ! -name MANIFEST_SHA256 -print0 \
    | LC_ALL=C sort -z \
    | xargs -0 sha256sum \
    | sed 's#  \./#  #' > MANIFEST_SHA256 )
NFILES="$(wc -l < "$STAGE/MANIFEST_SHA256")"

# ---- deterministic tar + gzip -n ---------------------------------------------- #
TAR="$OUT/$NAME.tar"
GZ="$TAR.gz"
rm -f "$TAR" "$GZ" "$GZ".part_* 2>/dev/null || true
# archive the whole staged dir so it unpacks into one self-contained
# "$NAME/" directory (BDIR); deterministic ordering + zeroed metadata.
tar \
    --sort=name --format=gnu --mtime='@0' \
    --owner=0 --group=0 --numeric-owner \
    -C "$OUT" -cf "$TAR" "$NAME"
gzip -n -9 "$TAR"       # -> $GZ

SZ="$(stat -c%s "$GZ")"
SHA="$(sha256sum "$GZ" | awk '{print $1}')"
SZ_MB="$(awk "BEGIN{printf \"%.2f\", $SZ/1048576}")"

echo ""
echo "================ ROUND63 S2 BUNDLE ================"
echo " bundle : $GZ"
echo " files  : $NFILES"
echo " size   : $SZ bytes ($SZ_MB MB)"
echo " sha256 : $SHA"

# ---- split if over threshold -------------------------------------------------- #
if [ "$SZ" -gt "$SPLIT_THRESHOLD" ]; then
  split -b "$SPLIT_CHUNK" -d -a 3 "$GZ" "$GZ.part_"
  echo " split  : >50 MB -> $SPLIT_CHUNK chunks:"
  for p in "$GZ".part_*; do
    echo "            $(basename "$p")  $(stat -c%s "$p") bytes  $(sha256sum "$p" | awk '{print $1}')"
  done
  echo " reassemble on the VM:"
  echo "            cat $NAME.tar.gz.part_* > $NAME.tar.gz"
  echo "            sha256sum -c <(echo \"$SHA  $NAME.tar.gz\")"
else
  echo " split  : not needed (<=50 MB) -- upload $NAME.tar.gz whole"
fi
echo "=================================================="
echo ""
echo "verify after transfer:  echo \"$SHA  $NAME.tar.gz\" | sha256sum -c -"
echo "unpack:                 tar -xzf $NAME.tar.gz    # -> $NAME/ (repo-root layout)"
echo "self-check bundled files: ( cd $NAME && sha256sum -c MANIFEST_SHA256 --quiet )"
