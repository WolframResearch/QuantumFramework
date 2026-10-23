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
            Head[current] === StabilizerFrame && Head[legacy] === StabilizerFrame, Simplify[current["State"]["StateVector"] - legacy["State"]["StateVector"]] === ConstantArray[0, 2^Max[1, qc["Max"]]],
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

rvSymbolicInitial = PauliStabilizer[<|"Signs" -> Array[(-1)^rvBit[#] &, 10], "Tableau" -> PauliStabilizer[5]["Tableau"]|>];
rvSymbolicSpecs = {"H" -> 5, "S" -> 1, "Y" -> 3, "CNOT" -> {5, 2}, "CZ" -> {3, 1}, "SWAP" -> {4, 2}};
rvA = rvCircuit[rvSymbolicSpecs][rvSymbolicInitial]; rvB = Fold[#1[#2] &, rvSymbolicInitial, rvSymbolicSpecs];
Echo[{PauliStabilizerQ[rvSymbolicInitial], rvA["Tableau"] === rvB["Tableau"], rvA["Signs"], rvB["Signs"], FullSimplify[rvA["Signs"] - rvB["Signs"], Element[Array[rvBit, 10], Integers]]}, "SIGNS_DETAIL"];
rvQ = rvCircuit[{"H" -> 1, "P"[rvAngle] -> 1, "CNOT" -> {1, 2}, "H" -> 2, "P"[-rvAngle] -> 2}];
rvNew = rvQ[PauliStabilizer[2]];
rvOld = Block[{Wolfram`QuantumFramework`GateMatrix`PackagePrivate`$stabilizerCircuitMemo = <||>}, Internal`InheritedBlock[{stabilizerCircuitSpecs}, Clear[stabilizerCircuitSpecs]; stabilizerCircuitSpecs[x_QuantumCircuitOperator] := rvLegacySpecs[x]; rvQ[PauliStabilizer[2]]]];
Echo[{Head[rvNew], Head[rvOld], rvNew["State"]["StateVector"], rvOld["State"]["StateVector"], FullSimplify[rvNew["State"]["StateVector"] - rvOld["State"]["StateVector"], Element[rvAngle, Reals]]}, "PHASE_DETAIL"];
rvSmall = rvCircuit[{QuantumOperator[rvH, {3}]}];
Echo[{rvSmall[PauliStabilizer[2]]["Stabilizers"], Block[{Wolfram`QuantumFramework`GateMatrix`PackagePrivate`$stabilizerCircuitMemo = <||>}, Internal`InheritedBlock[{stabilizerCircuitSpecs}, Clear[stabilizerCircuitSpecs]; stabilizerCircuitSpecs[x_QuantumCircuitOperator] := rvLegacySpecs[x]; rvSmall[PauliStabilizer[2]]["Stabilizers"]]]}, "EXTEND_DETAIL"];
Exit[];
