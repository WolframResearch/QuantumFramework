(* ::Package:: *)

(* harness.wl : correctness oracle (QF dense) + exact differential testing +     *)
(* randomized-circuit generation for the path-sum simulator.                     *)
(* Assumes Wolfram`QuantumFramework` and pathsum.wl are loaded.                   *)

(* ---- QF dense oracle: exact <y|C|0...0> via Method->"Schrodinger" ---- *)
PathSum`qfGate = {
    {"H", q_} :> "H" -> q, {"T", q_} :> "T" -> q, {SuperDagger["T"], q_} :> SuperDagger["T"] -> q,
    {"S", q_} :> "S" -> q, {"Z", q_} :> "Z" -> q, {"CZ", a_, b_} :> "CZ" -> {a, b},
    {"CCZ", a_, b_, c_} :> QuantumOperator[DiagonalMatrix[{1, 1, 1, 1, 1, 1, 1, -1}], {a, b, c}],
    {"CNOT", a_, b_} :> "CNOT" -> {a, b}
};

PathSum`qfAmp[n_Integer, gates_List, y_List] := With[
    {qco = QuantumCircuitOperator[gates /. PathSum`qfGate]},
    (Normal[qco[QuantumState[StringRepeat["0", n]], Method -> "Schrodinger"]["StateVector"]])[[1 + FromDigits[y, 2]]]
];

(* ---- exact differential checks (RootReduce equality, no floating point) ---- *)
PathSum`exactEqualQ[a_, b_] := TrueQ[RootReduce[a - b] === 0] || PossibleZeroQ[RootReduce[a - b]];

(* builder (un-reduced brute-force amplitude) vs QF, on all 2^n outputs *)
PathSum`builderMatchesQ[n_Integer, gates_List] := With[{ps = PathSum`build[n, gates]},
    AllTrue[Tuples[{0, 1}, n], PathSum`exactEqualQ[PathSum`amp[ps, n, #], PathSum`qfAmp[n, gates, #]] &]
];

(* reduced (eliminated) amplitude vs QF, on all 2^n outputs *)
PathSum`reducedMatchesQ[n_Integer, gates_List] := With[{r = PathSum`reduce[PathSum`build[n, gates]]},
    AllTrue[Tuples[{0, 1}, n], PathSum`exactEqualQ[PathSum`amp[r, n, #], PathSum`qfAmp[n, gates, #]] &]
];

(* reduced vs un-reduced (elimination preserves the amplitude) *)
PathSum`reductionSoundQ[n_Integer, gates_List] := With[
    {ps = PathSum`build[n, gates], r = PathSum`reduce[PathSum`build[n, gates]]},
    AllTrue[Tuples[{0, 1}, n], PathSum`exactEqualQ[PathSum`amp[r, n, #], PathSum`amp[ps, n, #]] &]
];

(* ---- randomized circuit generator over {H,T,S,Z,CZ,CCZ,CNOT} ---- *)
PathSum`randomGate[n_Integer] := Module[{kinds, k},
    kinds = Join[{"H", "T", "S", "Z"}, If[n >= 2, {"CZ", "CNOT"}, {}], If[n >= 3, {"CCZ"}, {}]];
    k = RandomChoice[kinds];
    Switch[k,
        "H" | "T" | "S" | "Z", {k, RandomInteger[{1, n}]},
        "CZ" | "CNOT", {k, Sequence @@ RandomSample[Range[n], 2]},
        "CCZ", {k, Sequence @@ RandomSample[Range[n], 3]}
    ]
];
PathSum`randomCircuit[n_Integer, depth_Integer] := Table[PathSum`randomGate[n], {depth}];
