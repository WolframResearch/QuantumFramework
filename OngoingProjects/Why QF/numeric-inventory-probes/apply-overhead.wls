(* Rows of the apply-overhead table that before-after.wls does not have, in a fresh kernel.
   Usage: apply-overhead.wls <QF paclet dir> first-state|oracle10-new-states
   first-state: a kernel's first QuantumState from a digit string, in ms.
   oracle10-new-states: the 10-variable phase oracle applied to five new random states after
   two applies to another, the median in s, and a digest of the last result. *)
{qfDir, row} = $ScriptCommandLine[[2 ;; 3]];
PacletDirectoryLoad[qfDir];
Needs["Wolfram`QuantumFramework`"];
out[fields___] := WriteString["stdout", StringRiffle[ToString /@ {"RESULT", row, fields}, "\t"], "\n"];
Switch[row,
  "first-state", out[Round[1000 First[AbsoluteTiming[QuantumState["0000000000"]]], 0.1]],
  "oracle10-new-states", Module[{qc, states, times, r},
    SeedRandom[3]; qc = QuantumCircuitOperator["PhaseOracle"[BooleanFunction[RandomInteger[2^(2^10) - 1], 10]]];
    SeedRandom[7]; states = Table[QuantumState[Normalize[RandomComplex[{-1 - I, 1 + I}, 2^10]], ConstantArray[2, 10]], 6];
    qc[states[[1]]]; qc[states[[1]]];
    times = Table[First[AbsoluteTiming[r = MemoryConstrained[qc[states[[i]]], 4 10^9]]], {i, 2, 6}];
    out[Round[Median[times], 0.001], StringTake[IntegerString[Hash[Round[N[Normal[r["StateVector"]]], 10^-10], "SHA256"], 16], 12]]
  ]
];
