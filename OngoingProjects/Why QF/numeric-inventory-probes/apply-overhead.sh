#!/bin/bash
# QF before (d682b341, pushed) and after (2180faba) the apply-overhead commits, both with the
# installed TensorNetworks and Arrays, alternating, three repeats, each row in a fresh kernel.
# Output lines: rep <tab> env <tab> row <tab> fields...
R=/Users/swish/src/wolfram/qf-networks/runs; M=/Users/swish/src/wolfram/qf-networks/memguard.sh
cd /Users/swish/src/wolfram/QuantumFramework
for rep in 1 2 3; do
  for env in before after; do
    case $env in before) qf=$R/qf-d682b341/QuantumFramework ;; after) qf=$R/qf-2180faba/QuantumFramework ;; esac
    for row in first-state oracle10-new-states; do
      r=$($M 6 timeout -k 10 400 wolframscript -file $R/apply-overhead.wls $qf $row 2>&1 | grep "^RESULT" | cut -f3-)
      echo -e "$rep\t$env\t$row\t${r:-killed}"
    done
    for row in apply-oracle10-zero apply-oracle10-uniform apply-qft14-zero apply-qft14-random quench-12 quench-16 quench-20; do
      r=$($M 6 timeout -k 10 400 wolframscript -file $R/before-after.wls $qf $row 2>&1 | grep "^RESULT" | cut -f3-)
      echo -e "$rep\t$env\t$row\t${r:-killed}"
    done
    for mode in greedy default; do
      out=$($M 7 timeout -k 10 400 wolframscript -file $R/grover.wls $qf $mode 2>&1)
      echo -e "$rep\t$env\tgrover8-$mode\t$(echo "$out" | grep "apply s" | awk '{print $NF}')\t$(echo "$out" | grep "max probability" | awk '{print $NF}')"
    done
  done
done
