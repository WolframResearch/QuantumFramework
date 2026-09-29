#!/bin/bash
# Every run behind the report, one kernel at a time, on paclet copies beside this
# script: base (main at a1799b61), proto (sections 5.1-5.3) and protos (proto plus
# section 5.4), each exported with git archive and edited by apply_proto.py, with
# QuantumOperatorDiagonalExp.wlt copied into its Tests directory.
cd "$(dirname "$0")"
for t in base proto protos; do wolframscript -file $t/Tests/RunTests.wls DiagonalExp > file-$t.out 2>&1; done
wolframscript -file proto/Tests/RunTests.wls > suite-proto.out 2>&1
wolframscript -file protos/Tests/RunTests.wls > suite-protos.out 2>&1
for t in base proto; do wolframscript -file correct.wls $t > correct-$t.out 2>&1; done
for t in base proto; do wolframscript -file accuracy.wls $t > accuracy-$t.out 2>&1; done
wolframscript -file followon.wls base > followon-base.out 2>&1
for s in resonant edge nested; do for t in base proto; do wolframscript -file $s.wls $t; done > $s.out 2>&1; done
for s in reach liouvillian-matrix; do for t in proto protos; do wolframscript -file $s.wls $t; done > $s.out 2>&1; done
wolframscript -file tolerance.wls > tolerance.out 2>&1
wolframscript -file structured.wls > structured.out 2>&1
for t in proto protos base; do
  for g in detect ising liouvillian qaoa controls; do
    if [ "$t" = "protos" ] && [ "$g" != "liouvillian" ]; then continue; fi
    wolframscript -file bench.wls $t $g
  done
done > bench.out 2>&1
date > run-all.done
