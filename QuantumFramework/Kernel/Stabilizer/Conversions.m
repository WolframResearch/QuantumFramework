Package["Wolfram`QuantumFramework`"]

PackageScope[PauliRow]
PackageScope[QuantumOperatorTableau]



(* ============================================================================ *)
(* Local helpers used by PauliRow.                                              *)
(* ============================================================================ *)

firstIndex[row_] := Replace[FirstPosition[Normal[row], x_ /; Abs[x] == 1, {1}, Heads -> False], {i_} :> i - 1]

bitvector[x_, n_, d_ : 2] := Reverse @ IntegerDigits[x, d, n]


(* ============================================================================ *)
(* PauliRow: decode a Pauli matrix back to (xbits, zbits, phase).               *)
(* Used by QuantumOperatorTableau to extract a tableau row from a conjugated    *)
(* Pauli matrix.                                                                *)
(* ============================================================================ *)

PauliRow[mat_ ? MatrixQ, n_Integer ? Positive, d : _Integer ? Positive : 2] := Enclose @ Block[
    {xint, xbits, zbits, entries, positivePhase, phase, coef},
    xint = Confirm @ firstIndex[mat[[1]]];
    xbits = bitvector[xint, n, d];
    entries = MapIndexed[
        With[{i = #2[[1]] - 1, index = Confirm @ firstIndex[#1]},
            ConfirmAssert[bitvector[index, n, d] == Mod[xbits + Threaded[bitvector[i, n, d]], d]
                && MemberQ[Table[Exp[I a 2 Pi], {a, 0, 1 - 1 / 2 / d, 1 / 2 / d}], #1[[index + 1]]]];
            #1[[index + 1]]
        ] &,
        mat
    ];

    zbits = Confirm @ Replace[entries[[d ^ # + 1]] / entries[[1]], {1 -> 0, -1 -> 1, _ -> Missing[]}] & /@ Range[0, n - 1];
    positivePhase = (-I) ^ Inner[Mod[#1 #2, d] &, xbits, zbits];
    phase = Confirm @ Which[
        positivePhase == entries[[1]], 0,
        positivePhase == - entries[[1]], 1,
        True, Missing[]
    ];
    coef = ((-1) ^ phase) positivePhase;
    ConfirmAssert[entries == Flatten[coef kroneckerProduct @@ Replace[Reverse[zbits], {1 -> Diagonal[pauliMatrix[3, d]], 0 -> ConstantArray[1, d]}, {1}]]];

    (* Reverse bit convention *)
    Join[Reverse[xbits], Reverse[zbits], {phase}]
]


QuantumOperatorTableau[qo_QuantumOperator] /; qo["InputQudits"] == qo["OutputQudits"] && Equal @@ qo["Dimensions"] := Enclose @ With[
    {n = qo["InputQudits"], d = First @ qo["Dimensions"]},
    Join[
        Confirm @ PauliRow[(qo @ QuantumOperator[{"X", d} -> #] @ qo["Dagger"])["Matrix"], n, d] & /@ Sort[qo["InputOrder"]],
        Confirm @ PauliRow[(qo @ QuantumOperator[{"Z", d} -> #] @ qo["Dagger"])["Matrix"], n, d] & /@ Sort[qo["InputOrder"]]
    ]
]


(* ============================================================================ *)
(* PauliStabilizer -> QuantumState.                                             *)
(* Phase 1 cleanup: dropped "QuantumSttate" typo alias                          *)
(* (OngoingProjects/Stabilizer/paulistabilizer-source-audit.md \[Section]13.7).  *)
(* Phase 5c: optional "GlobalPhase" association key (set by                     *)
(* PauliStabilizer[qs_QuantumState] / PauliStabilizer[qo_QuantumOperator])      *)
(* recovers the overall phase that the stabilizer tableau drops, so that         *)
(* PauliStabilizer[qs]["State"] === qs exactly. Default 1 when absent.           *)
(* Cost: O(4^n) -- materializes 2^n x 2^n matrices. Practical limit n <= 10.    *)
(* ============================================================================ *)

ps_PauliStabilizer["State" | "QuantumState"] := With[{
    canonical = Normalize @ Total @ NullSpace[
        Total @ MapThread[(IdentityMatrix[Length[#2]] - #1 #2) &, {
            ps["StabilizerSigns"],
            kroneckerProduct @@@
                Map[PauliMatrix, Replace[Transpose[ps["StabilizerTableau"], {3, 2, 1}], {{0, 0} -> 0, {1, 0} -> 1, {1, 1} -> 2, {0, 1} -> 3}, {2}], {2}]
        }]
    ],
    phase = Lookup[First[ps], "GlobalPhase", 1]
},
    QuantumState[phase canonical]
]


(* ============================================================================ *)
(* PauliStabilizer -> QuantumCircuitOperator (AG canonical decomposition).      *)
(* Greedy: for each qubit, set destab-X to a single 1, zero destab-Z, zero     *)
(* stab-X. Returns a Clifford circuit whose dagger equals `ps`.                 *)
(* Performance: O(n^3) gates. Circuit length is NOT minimized -- Reid24,        *)
(* Winderl23 propose better (Phase 5).                                          *)
(* ============================================================================ *)

ps_PauliStabilizer["Circuit" | "QuantumCircuit" | "QuantumCircuitOperator"] := Block[{
    clifford = ps,
    n = ps["Qudits"],
    gates = {},
    destabX, destabZ, stabX, stabZ, destabP, stabP,
    append, setQubitX1, setRowX0, setRowZ0
},
    destabX[q_] := Thread[clifford["DestabilizerX"][[All, q]] == 1];
    destabZ[q_] := Thread[clifford["DestabilizerZ"][[All, q]] == 1];
    stabX[q_]   := Thread[clifford["StabilizerX"][[All, q]] == 1];
    stabZ[q_]   := Thread[clifford["StabilizerZ"][[All, q]] == 1];
    destabP[q_] := clifford["Phase"][[q]] == 1;
    stabP[q_]   := clifford["Phase"][[n + q]] == 1;
    append[gate_] := (clifford = clifford[gate]; AppendTo[gates, gate]);
    setQubitX1[q_] := Catch @ With[{x := destabX[q], z := destabZ[q]},
        If[x[[q]], Throw[0]];
        Do[If[x[[i]], append["SWAP" -> {i, q}]; Throw[0]], {i, q + 1, n}];
        Do[If[z[[i]], append["H" -> i]; If[i != q, append["SWAP" -> {i, q}]]; Throw[0]], {i, q, n}]
    ];
    setRowX0[q_] := With[{x := destabX[q], z := destabZ[q]},
        Do[If[x[[i]], append["CNOT" -> {q, i}]], {i, q + 1, n}];
        If[Or @@ z[[q ;;]], If[! z[[q]], append["S" -> q]]; Do[If[z[[i]], append["CNOT" -> {i, q}]], {i, q + 1, n}]; append["S" -> q]]
    ];
    setRowZ0[q_] := With[{x := stabX[q], z := stabZ[q]},
        If[Or @@ z[[q + 1 ;;]], Do[If[z[[i]], append["CNOT" -> {i, q}]], {i, q + 1, n}]];
        If[Or @@ x[[q ;;]], append["H" -> q]; Do[If[x[[i]], append["CNOT" -> {q, i}]], {i, q + 1, n}]; If[z[[q]], append["S" -> q]]; append["H" -> q]]
    ];
    Do[setQubitX1[q]; setRowX0[q]; setRowZ0[q], {q, n}];
    Do[If[destabP[q], append["Z" -> q]]; If[stabP[q], append["X" -> q]], {q, n}];
    gates = FixedPoint[SequenceReplace[{g_, g_} -> Nothing], gates];
    QuantumCircuitOperator[gates]["Dagger"]
]


ps_PauliStabilizer["Operator" | "QuantumOperator"] := With[{
    canonical = QuantumOperator["I", Range[ps["Qudits"]]] @ ps["Circuit"]["QuantumOperator"],
    phase = Lookup[First[ps], "GlobalPhase", 1]
},
    phase canonical
]


(* ============================================================================ *)
(* UpValues: PauliStabilizer integration with the rest of QuantumFramework.     *)
(* ============================================================================ *)

qo_QuantumOperator[ps_PauliStabilizer] ^:= PauliStabilizerApply[QuantumCircuitOperator[qo], ps]
QuantumState[ps_PauliStabilizer] ^:= ps["State"]
QuantumCircuitOperator[ps_PauliStabilizer] := ps["Circuit"]
QuantumOperator[ps_PauliStabilizer] := ps["Circuit"]["QuantumOperator"]

(* Note: qmo_QuantumMeasurementOperator[ps_PauliStabilizer] used to be defined  *)
(* here as PauliStabilizerApply[QuantumCircuitOperator[qmo], ps]. Phase 7.1     *)
(* (2026-05-06) replaces that generic path with a basis-aware dispatch in       *)
(* Stabilizer/HybridInterop.m so Pauli-string measurements stay in the tableau  *)
(* (cost O(n^2)) instead of routing through circuit conversion.                 *)
