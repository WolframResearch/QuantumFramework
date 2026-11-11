#!/bin/bash
# Runs every row of before-after.wls in three environments, three times each, alternating the
# environments within each repeat, and appends the raw results to before-after.tsv:
#   repeat, environment, TN/Arrays versions, 1-minute load average, row, timings, check.
# Summarize with: python3 -I before-after-table.py
#
# Environments:
#   before  QF 0bd75989 with TensorNetworks 1.0.11 and Arrays 1.4.1, from a separate user base
#   middle  QF 0bd75989 with the installed TensorNetworks 1.1.0 and Arrays 1.4.2
#   after   QF abe5a3d7 with the installed TensorNetworks 1.1.0 and Arrays 1.4.2
#
# Setup, once, in $WORK:
#   git archive 0bd75989 QuantumFramework | tar -x -C $WORK/qf-0bd75989   (likewise for abe5a3d7)
#   mkdir $WORK/userbase-before && ln -s ~/Library/Wolfram/Licensing $WORK/userbase-before/Licensing
#   WOLFRAM_USERBASE=$WORK/userbase-before wolframscript -code \
#     'PacletInstall[{"Wolfram/TensorNetworks", "1.0.11"}]; PacletInstall[{"Wolfram/Arrays", "1.4.1"}]'
# GUARD kills a run whose memory passes 8 GB (memguard.sh from the qf-networks folder).
PROBES=$(cd "$(dirname "$0")" && pwd)
WORK=${WORK:-$HOME/src/wolfram/qf-networks/runs}
GUARD=${GUARD:-$HOME/src/wolfram/qf-networks/memguard.sh}
OUT=$PROBES/before-after.tsv
rows="build-qft14 build-oracle10 controlled-gate apply-qft14-zero apply-qft14-random apply-oracle10-zero apply-oracle10-uniform quench-12 quench-16 quench-20 fourier-basis-16"
echo "# started $(date '+%Y-%m-%d %H:%M:%S') on $(sysctl -n machdep.cpu.brand_string), $(sysctl -n hw.memsize | awk '{print $1 / 2^30}') GB" >> "$OUT"
for rep in 1 2 3; do
    for row in $rows; do
        for env in before middle after; do
            case $env in
                before) ub=$WORK/userbase-before; qf=$WORK/qf-0bd75989/QuantumFramework ;;
                middle) ub=""; qf=$WORK/qf-0bd75989/QuantumFramework ;;
                after) ub=""; qf=$WORK/qf-abe5a3d7/QuantumFramework ;;
            esac
            load=$(uptime | sed 's/.*load averages: //' | cut -d' ' -f1)
            if [ -n "$ub" ]; then
                res=$(WOLFRAM_USERBASE=$ub "$GUARD" 8 timeout -k 10 900 wolframscript -file "$PROBES/before-after.wls" "$qf" $row 2>&1)
            else
                res=$("$GUARD" 8 timeout -k 10 900 wolframscript -file "$PROBES/before-after.wls" "$qf" $row 2>&1)
            fi
            versions=$(echo "$res" | grep "^ENV" | cut -f2-3 | tr '\t' '/')
            line=$(echo "$res" | grep "^RESULT" | cut -f2-)
            [ -z "$line" ] && line="$row	no result	$(echo "$res" | grep -E 'MEMGUARD|::' | head -1)"
            printf "%s\t%s\t%s\t%s\t%s\n" "$rep" "$env" "$versions" "$load" "$line" >> "$OUT"
        done
    done
done
echo "# finished $(date '+%Y-%m-%d %H:%M:%S')" >> "$OUT"
