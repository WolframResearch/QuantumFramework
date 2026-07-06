#!/bin/bash
# Verify + rebuild per-part answer notebooks.
# Usage:  ./_check.sh                 # all parts
#         ./_check.sh Part-05.md      # one or more named parts
# For each part: runs every wl cell in a fresh kernel (flags errors/messages),
# then rebuilds the sibling .nb via md2nb.
DIR="$(cd "$(dirname "$0")" && pwd)"
if [ "$#" -eq 0 ]; then set -- "$DIR"/Part-*.md; fi
for f in "$@"; do
  case "$f" in /*) p="$f" ;; *) p="$DIR/$f" ;; esac
  wolframscript -file "$DIR/_check-one.wls" "$p"
done
