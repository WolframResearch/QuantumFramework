(* One row of the before/after table, in a fresh kernel: QuantumFramework from the given directory,
   with the TensorNetworks and Arrays the user base provides. Only the workload is timed.
   Usage: before-after.wls <QF paclet dir> <row>
   Prints ENV <tab> TN version <tab> Arrays version, then
   RESULT <tab> row <tab> seconds (or milliseconds per gate) <tab> a check of the result. *)
{qfDir, row} = $ScriptCommandLine[[2 ;; 3]];
PacletDirectoryLoad[qfDir];
Needs["Wolfram`QuantumFramework`"];
WriteString["stdout", "ENV\t", First[PacletFind["Wolfram/TensorNetworks"]]["Version"], "\t", First[PacletFind["Wolfram/Arrays"]]["Version"], "\n"];

SetAttributes[sec, HoldFirst];
sec[e_, limit_ : 240] := MemoryConstrained[
    TimeConstrained[With[{t = AbsoluteTiming[e]}, {ToString[N[Round[First[t], 1 / 1000]], InputForm], Last[t]}], limit, {"timeout", None}],
    4 10^9, {"memory", None}];
digest[qs_QuantumState] := StringTake[IntegerString[Hash[Round[N[Normal[qs["StateVector"]]], 10^-10], "SHA256"], 16], 12];
digest[other_] := StringTake[ToString[other, InputForm], UpTo[40]];
out[fields___] := WriteString["stdout", StringRiffle[ToString /@ {"RESULT", row, fields}, "\t"], "\n"];
zero[n_] := QuantumState[StringJoin[ConstantArray["0", n]]];
oracle[] := (SeedRandom[3]; QuantumCircuitOperator["PhaseOracle"[BooleanFunction[RandomInteger[2^(2^10) - 1], 10]]]);
quench[n_] := QuantumCircuitOperator @ Flatten @ Table[Join[("R"[0.4, "ZZ"] -> # &) /@ Partition[Range[n], 2, 1], ("RX"[0.6] -> # &) /@ Range[n]], {4}];
twice[qc_, qs_] := Module[{a1, r1, a2, r2},
    {a1, r1} = sec[qc[qs]];
    {a2, r2} = sec[qc[qs]];
    out[a1, a2, digest[r1]]
];

Switch[row,
    "build-qft14", With[{t = sec[QuantumCircuitOperator["Fourier"[14]]]}, out[First[t], Length[Last[t]["Elements"]]]],
    "build-oracle10", With[{t = sec[oracle[]]}, out[First[t], Length[Last[t]["Elements"]]]],
    "controlled-gate", Module[{gates, t},
        gates[m_] := Catenate @ Table[{# + 1, {# + i + m}, {i + m}} & /@ Range[14 - i], {i, 14}];
        Scan[QuantumOperator["Controlled"["PhaseShift"[#[[1]]], #[[2]]], #[[3]]] &, gates[20]];
        t = First[AbsoluteTiming[Scan[QuantumOperator["Controlled"["PhaseShift"[#[[1]]], #[[2]]], #[[3]]] &, gates[40]]]];
        out[ToString[N[Round[1000 t / 91, 1 / 100]], InputForm], "ms per gate"]
    ],
    "apply-qft14-zero", twice[QuantumCircuitOperator["Fourier"[14]], zero[14]],
    "apply-qft14-random", With[{qs = (SeedRandom[1]; QuantumState[Normalize[RandomComplex[{-1 - I, 1 + I}, 2^14]]])}, twice[QuantumCircuitOperator["Fourier"[14]], qs]],
    "apply-oracle10-zero", twice[oracle[], zero[10]],
    "apply-oracle10-uniform", twice[oracle[], QuantumState[ConstantArray[1 / Sqrt[2.^10], 2^10], ConstantArray[2, 10]]],
    "quench-12" | "quench-16" | "quench-20", With[{n = ToExpression[StringDrop[row, 7]]},
        Module[{t, psi},
            {t, psi} = sec[quench[n][zero[n]]];
            out[t, If[MatchQ[psi, _QuantumState], ToString[Chop[(psi["Dagger"] @ QuantumOperator["Z", {Ceiling[n / 2]}] @ psi)["Scalar"]], InputForm], "none"]]
        ]
    ],
    "fourier-basis-16", With[{s = (SeedRandom[2]; QuantumState[Normalize[RandomComplex[{-1 - I, 1 + I}, 16]]]), b = QuantumBasis["Fourier"[16]]},
        With[{t = sec[QuantumState[s, b], 60]}, out[First[t], digest[Last[t]]]]
    ],
    _, out["unknown row"]
];
