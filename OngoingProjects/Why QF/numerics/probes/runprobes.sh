#!/bin/zsh
# run each probe in its own kernel under an OS timeout; results go straight to the log
S="$(cd "$(dirname "$0")" && pwd)"
LOG=$S/results.log
for f in "$@"; do
  echo "### $f  $(date +%H:%M:%S)" >> $LOG
  timeout -k 15 ${PROBE_TIMEOUT:-420} wolframscript -f $S/$f >> $LOG 2>&1
  echo "### exit $? for $f  $(date +%H:%M:%S)" >> $LOG
done
