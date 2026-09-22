#!/bin/bash
# Verify + rebuild per-part answer notebooks.
# Usage:  ./_check.sh                 # all parts whose notebook carries no evaluated outputs
#         ./_check.sh Part-05.md      # one or more named parts
#         ./_check.sh --force Part-11.md   # rebuild even if the notebook carries evaluated outputs
# For each part: runs every wl cell in a fresh kernel (flags errors/messages),
# then rebuilds the sibling .nb via md2nb without outputs. A part whose existing
# notebook holds Output cells is skipped unless --force is given, because the
# rebuild would drop those outputs. Exits 1 if any part flagged a cell or failed to build.
DIR="$(cd "$(dirname "$0")" && pwd)"
FORCE=0
args=()
for a in "$@"; do
  case "$a" in --force) FORCE=1 ;; *) args+=("$a") ;; esac
done
set -- "${args[@]}"
if [ "$#" -eq 0 ]; then set -- "$DIR"/Part-*.md; fi
rc=0
for f in "$@"; do
  case "$f" in /*) p="$f" ;; *) p="$DIR/$f" ;; esac
  nb="${p%.md}.nb"
  if [ "$FORCE" -eq 0 ] && [ -f "$nb" ] && grep -q '"Output"' "$nb"; then
    echo "$(basename "$p"): skipped, $(basename "$nb") carries evaluated outputs (use --force to rebuild without them)" >&2
    continue
  fi
  wolframscript -file "$DIR/_check-one.wls" "$p" || rc=1
done
exit "$rc"
