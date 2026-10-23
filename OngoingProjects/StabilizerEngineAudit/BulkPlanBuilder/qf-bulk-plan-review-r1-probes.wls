PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
$ProgressReporting = False;
Echo[{$Version, FindFile["Wolfram`QuantumFramework`"], PacletFind["Wolfram/Arrays"][[1]]["Version"], PacletFind["Wolfram/TensorNetworks"][[1]]["Version"]}, "ENVIRONMENT"];

rvCircuit[ops_List] := QuantumCircuitOperator[<|"Elements" -> (FromCircuitOperatorShorthand /@ ops), "Label" -> None|>];
rvEqual[a_PauliStabilizer, b_PauliStabilizer] := a["Tableau"] === b["Tableau"] && a["Signs"] === b["Signs"];
rvEqual[_, _] := False;
rvPauli[s_String] := With[{negative = StringStartsQ[s, "-"]},
    (-1)^Boole[negative] (KroneckerProduct @@ (Lookup[<|"I" -> PauliMatrix[0], "X" -> PauliMatrix[1], "Y" -> PauliMatrix[2], "Z" -> PauliMatrix[3]|>, Characters[If[negative, StringDrop[s, 1], s]]]))
];
rvGenerators[n_] := Join[Table[StringReplacePart[StringRepeat["I", n], "X", {j, j}], {j, n}], Table[StringReplacePart[StringRepeat["I", n], "Z", {j, j}], {j, n}]];
(* A direct basis-index embedding, independently of circuit compilation and of
   all stabilizer classification/remapping helpers. QF matrices index sorted wires. *)
rvEmbed[u_, wires_, n_] := With[{bits = IntegerDigits[Range[0, 2^n - 1], 2, n], rest = Complement[Range[n], wires]},
    Table[If[bits[[r, rest]] === bits[[c, rest]],
        u[[FromDigits[bits[[r, wires]], 2] + 1, FromDigits[bits[[c, wires]], 2] + 1]], 0], {r, 2^n}, {c, 2^n}]
];
rvDenseU[ops_, n_] := Fold[#2 . #1 &, IdentityMatrix[2^n], rvEmbed[Normal[# ["Sort"]["Matrix"]], Sort[# ["InputOrder"]], n] & /@ ops];
rvDenseCheck[ops_, n_, numeric_:False] := With[{u = rvDenseU[ops, n], ps = rvCircuit[ops][PauliStabilizer[n]]},
    Head[ps] === PauliStabilizer && With[{actual = rvPauli /@ Join[ps["Destabilizers"], ps["Stabilizers"]]},
        If[numeric,
            Max[Abs[N[actual - (u . rvPauli[#] . ConjugateTranspose[u] & /@ rvGenerators[n])]]] < 10^-8,
            Simplify[actual - (u . rvPauli[#] . ConjugateTranspose[u] & /@ rvGenerators[n])] === ConstantArray[0, {2 n, 2^n, 2^n}]
        ]
    ]
];
rvLegacySpecs[qc_] := Catenate[rvLegacyElement /@ qc["Operators"]];
rvLegacyElement[qc_QuantumCircuitOperator] := rvLegacySpecs[qc];
rvLegacyElement[op_QuantumOperator] := stabilizerOperatorSpecs[op];
rvLegacyElement[other_] := QuantumShortcut[other];
rvEquivalent[qc_] := With[{current = qc[PauliStabilizer[Max[1, qc["Max"]]]]},
    With[{legacy = Block[{Wolfram`QuantumFramework`GateMatrix`PackagePrivate`$stabilizerCircuitMemo = <||>},
        Internal`InheritedBlock[{stabilizerCircuitSpecs}, Clear[stabilizerCircuitSpecs]; stabilizerCircuitSpecs[x_QuantumCircuitOperator] := rvLegacySpecs[x]; qc[PauliStabilizer[Max[1, qc["Max"]]]]]]},
        Which[Head[current] === PauliStabilizer && Head[legacy] === PauliStabilizer, rvEqual[current, legacy],
            Head[current] === StabilizerFrame && Head[legacy] === StabilizerFrame, Simplify[Normal[current["State"]["StateVector"] - legacy["State"]["StateVector"]]] === ConstantArray[0, 2^Max[1, qc["Max"]]],
            True, current === legacy]
    ]
];

rvH = {{1, 1}, {1, -1}}/Sqrt[2];
rvCX = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}};
rvCZ = DiagonalMatrix[{1, 1, 1, -1}];
rvCustom = rvCX . KroneckerProduct[rvH, DiagonalMatrix[{1, I}]] . rvCZ;
rvCases = <|
    "AsymmetricInputOutputPermutation" -> {QuantumOperator[rvCX, {{4, 2}, {2, 4}}, "Label" -> "H"], QuantumOperator[rvCX, {{1, 3}, {3, 1}}, "Label" -> "I"], QuantumOperator[rvCX, {{4, 1}, {1, 4}}, "Label" -> "T"]},
    "RepeatedCustomEntangler" -> {QuantumOperator[rvCustom, {4, 1}, "Label" -> "I"], QuantumOperator[rvCustom, {1, 3}, "Label" -> "CNOT"], QuantumOperator[rvCustom, {4, 2}, "Label" -> "T"]},
    "AlgebraicGlobalPhase" -> {QuantumOperator[Exp[I Pi/11] rvCustom, {4, 1}], QuantumOperator[-I rvH, {3}], QuantumOperator[-rvCustom, {2, 4}]},
    "MixedStoredRepresentations" -> {QuantumOperator[rvCX, {4, 2}], QuantumOperator[SparseArray[rvCX], {1, 3}], QuantumOperator[N[rvH], {2}], QuantumOperator[N[rvCustom, 50], {4, 1}], QuantumOperator[rvCustom, {2, 3}]}
|>;
rvDenseResults = KeyValueMap[Function[{name, ops}, name -> rvDenseCheck[ops, 4, name === "MixedStoredRepresentations"]], rvCases];
Echo[rvDenseResults, "INDEPENDENT_DENSE"];

SeedRandom[104729];
rvSpecs[n_, count_] := Table[With[{w = RandomSample[Range[n], 2], g = RandomChoice[{"H", "S", "X", "Y", "Z", "V", SuperDagger["V"], "CNOT", "CZ", "SWAP"}]}, g -> If[MemberQ[{"CNOT", "CZ", "SWAP"}, g], w, First[w]]], count];
rvRandomResults = Table[With[{ops = FromCircuitOperatorShorthand /@ rvSpecs[4, 32]}, rvDenseCheck[ops, 4]], 8];
Echo[rvRandomResults, "RANDOM_DENSE_EXACT"];

rvLargeResults = Table[With[{specs = Join[{"H" -> n, "CNOT" -> {n, 1}, "S" -> 1, "CZ" -> {1, n}}, rvSpecs[n, 96]]}, With[{qc = rvCircuit[specs], initial = PauliStabilizer[n]["ApplyCircuit", rvSpecs[n, 12]]},
    rvEqual[qc[initial], Fold[#1[#2] &, PauliStabilizer[<|"Signs" -> initial["Signs"], "Tableau" -> initial["Tableau"]|>], specs]]]], {n, {31, 32, 61, 62, 63, 124, 125, 257, 1024}}];
Echo[rvLargeResults, "PACKED_VS_SCALAR"];

rvSymbolicInitial = PauliStabilizer[<|"Signs" -> Array[(-1)^rvBit[#] &, 10], "Tableau" -> PauliStabilizer[5]["Tableau"]|>];
rvSymbolicSpecs = {"H" -> 5, "S" -> 1, "Y" -> 3, "CNOT" -> {5, 2}, "CZ" -> {3, 1}, "SWAP" -> {4, 2}};
Echo[With[{a = rvCircuit[rvSymbolicSpecs][rvSymbolicInitial], b = Fold[#1[#2] &, rvSymbolicInitial, rvSymbolicSpecs]}, a["Tableau"] === b["Tableau"] && And @@ (With[{rules = Thread[Array[rvBit, 10] -> #]}, (a["Signs"] /. rules) === (b["Signs"] /. rules)] & /@ Tuples[{0, 1}, 10])], "SYMBOLIC_SIGNS"];
Echo[rvEquivalent[rvCircuit[{"H" -> 1, "P"[rvAngle] -> 1, "CNOT" -> {1, 2}, "H" -> 2, "P"[-rvAngle] -> 2}]], "SYMBOLIC_PHASE_LEGACY"];
Echo[With[{n = 257, specs = Catenate[ConstantArray[{"H" -> 257, "CNOT" -> {257, 1}, "S" -> 1, "CZ" -> {1, 257}}, 250]]}, With[{qc = rvCircuit[specs]}, rvEqual[Nest[qc, PauliStabilizer[n], 5], Fold[#1[#2] &, PauliStabilizer[<|"Signs" -> PauliStabilizer[n]["Signs"], "Tableau" -> PauliStabilizer[n]["Tableau"]|>], Catenate[ConstantArray[specs, 5]]]]]], "LONG_5000_VS_SCALAR"];

rvReport = TestReport[{
    VerificationTest[rvCircuit[{}][PauliStabilizer[0]], PauliStabilizer[0], TestID -> "ZeroRegisterEmpty"],
    VerificationTest[rvEquivalent[rvCircuit[{rvCircuit[{}], "Barrier", rvCircuit[{"H" -> 2, "CNOT" -> {2, 1}}], rvCircuit[{}], "S" -> 2}]], True, TestID -> "NestedEmptyMixed"],
    VerificationTest[rvCircuit[{"H" -> 1, QuantumMeasurementOperator["Computational", {1}], "S" -> 1}][PauliStabilizer[1]], $Failed, {PauliStabilizer::nonclifford}, TestID -> "MixedMeasurementRefusal"],
    VerificationTest[rvCircuit[{QuantumOperator[rvH, {{1}, {2}}], QuantumOperator[rvH, {{2}, {1}}]}][PauliStabilizer[2]], $Failed, {PauliStabilizer::nonclifford}, TestID -> "MovedWiresSamePattern"],
    VerificationTest[rvCircuit[{QuantumOperator[rvCX, {{1, 3}, {1, 3}}], QuantumOperator[rvCX, {{2, 3}, {1, 2}}]}][PauliStabilizer[3]], $Failed, {PauliStabilizer::nonclifford}, TestID -> "MovedWiresTwoQubit"],
    VerificationTest[rvCircuit[{QuantumOperator[rvH, {0}], QuantumOperator[rvH, {2}]}][PauliStabilizer[2]], $Failed, {PauliStabilizer::wires}, TestID -> "WireZero"],
    VerificationTest[rvCircuit[{QuantumOperator[rvH, {3}]}][PauliStabilizer[2]]["Stabilizers"], {"ZII", "IZI", "IIX"}, TestID -> "ExistingPublicCircuitRegisterExtension"],
    VerificationTest[rvCircuit[{QuantumOperator[rvH, {1}], QuantumOperator[QuantumOperator[rvH, {2}], QuantumBasis["PauliX"]]}][PauliStabilizer[2]], $Failed, {PauliStabilizer::nonclifford}, TestID -> "AtomicAlternateBasis"],
    VerificationTest[rvCircuit[{QuantumOperator["X"[4]], QuantumOperator[rvCX]}][PauliStabilizer[2]], $Failed, {PauliStabilizer::nonclifford}, TestID -> "FourLevelIsNotTwoQubits"],
    VerificationTest[rvCircuit[{QuantumOperator[{{1, 0}, {0, 2}}], "H" -> 1}][PauliStabilizer[1]], $Failed, {PauliStabilizer::nonclifford}, TestID -> "NonunitaryDiagonal"],
    VerificationTest[rvCircuit[{QuantumOperator["P"[ArcCos[rvTheta]]], QuantumOperator[QuantumOperator["P"[ArcCos[rvTheta]]], {2}, "Label" -> None]}][PauliStabilizer[2]], $Failed, {PauliStabilizer::nonclifford}, TestID -> "SymbolicMissingLabel"],
    VerificationTest[rvCircuit[{QuantumOperator["T"] @ QuantumOperator["H"], QuantumOperator[QuantumOperator["T"] @ QuantumOperator["H"], "Label" -> "I"]}][PauliStabilizer[1]], $Failed, {PauliStabilizer::nonclifford}, TestID -> "CompositeDishonestIdentityLabel"]
}];
Echo[{rvReport["TestsSucceededCount"], rvReport["TestsFailedCount"], ({#["TestID"], #["Outcome"], #["ActualMessages"]} & /@ Values[rvReport["TestResults"]])}, "EDGE_REPORT"];
Echo[AllTrue[Join[Last /@ rvDenseResults, rvRandomResults, rvLargeResults], TrueQ], "FINAL_DENSE_PACKED_PASS"];
rvSymbolicVector = Normal[rvCircuit[{"H" -> 1, "P"[rvAngle] -> 1, "CNOT" -> {1, 2}, "H" -> 2, "P"[-rvAngle] -> 2}][PauliStabilizer[2]]["State"]["StateVector"]];
rvSymbolicReference = {1, Exp[-I rvAngle], Exp[I rvAngle], -1}/2;
Echo[FullSimplify[Outer[Times, rvSymbolicVector, Conjugate[rvSymbolicVector]] - Outer[Times, rvSymbolicReference, Conjugate[rvSymbolicReference]], Element[rvAngle, Reals]] === ConstantArray[0, {4, 4}], "SYMBOLIC_PHASE_DENSE_DENSITY"];
Exit[If[rvReport["TestsFailedCount"] === 0 && AllTrue[Join[Last /@ rvDenseResults, rvRandomResults, rvLargeResults], TrueQ], 0, 1]];
