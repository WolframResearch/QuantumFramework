Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

BeginTestSection["BulkPlanBuilder"];

bpbCircuit[ops_List] := QuantumCircuitOperator[<|"Elements" -> (FromCircuitOperatorShorthand /@ ops), "Label" -> None|>];
bpbImages[ps_PauliStabilizer] := Join[ps["Destabilizers"], ps["Stabilizers"]];
bpbImages[_] := $Failed;
bpbEqual[a_PauliStabilizer, b_PauliStabilizer] := a["Tableau"] === b["Tableau"] && a["Signs"] === b["Signs"];
bpbEqual[_, _] := False;
bpbPauli[s_String] := KroneckerProduct @@ Replace[Characters[s],
    {"I" -> PauliMatrix[0], "X" -> PauliMatrix[1], "Y" -> PauliMatrix[2], "Z" -> PauliMatrix[3]}, {1}];
bpbDenseImages[qc_, n_] := With[{u = Normal @ bpbCircuit[
    Prepend[qc["Operators"], QuantumOperator[IdentityMatrix[2^n]]]]["QuantumOperator"]["Matrix"]},
    stabilizerPauliFromMatrix[Simplify[u . bpbPauli[#] . ConjugateTranspose[u]], n] & /@
        Join[Table[StringReplacePart[StringRepeat["I", n], "X", {j, j}], {j, n}],
            Table[StringReplacePart[StringRepeat["I", n], "Z", {j, j}], {j, n}]]
];
bpbRandomSpecs[n_, m_] := Table[With[{
    gate = RandomChoice[{"H", "S", SuperDagger["S"], "X", "Y", "Z", "V", SuperDagger["V"], "CNOT", "CZ", "SWAP"}],
    wires = RandomSample[Range[n], 2]
}, gate -> If[MemberQ[{"CNOT", "CZ", "SWAP"}, gate], wires, First[wires]]], {m}];

(* These are execution contracts: physical fallback alone would not detect
   accidental failure to extract an atomic QF object's local data. *)
VerificationTest[With[{
    local = Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorData[QuantumOperator["H"]],
    placed = Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorData[QuantumOperator["H", {3}, "Label" -> "T"]]
}, ListQ[local] && local === placed], True, {}, TestID -> "BPB-LocalDataExcludesLabelsAndAbsoluteWires"];

VerificationTest[With[{
    local = Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorData[QuantumOperator["CNOT", {1, 2}]],
    reversed = Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorData[QuantumOperator["CNOT", {2, 1}]]
}, ListQ[local] && ListQ[reversed] && local =!= reversed], True, {}, TestID -> "BPB-LocalDataRetainsWireOrderPattern"];

VerificationTest[Module[{calls = 0, specs},
    stabilizerClearCaches[];
    specs = Block[{
        Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorClassification =
            Function[Null, calls++; Missing["NotClifford"]]
    }, stabilizerCircuitSpecs[bpbCircuit[{
        QuantumOperator["H", {3}, "Label" -> "T"], QuantumOperator["S", {2}, "Label" -> "X"],
        QuantumOperator["CNOT", {1, 3}, "Label" -> "Y"], QuantumOperator["CZ", {1, 2}, "Label" -> None],
        QuantumOperator["V", {2}, "Label" -> "H"]
    }]]];
    calls === 0 && FreeQ[specs, _Missing] &&
        bpbEqual[PauliStabilizer[3]["ApplyCircuit", specs],
            PauliStabilizer[3]["ApplyCircuit", {"H" -> 3, "S" -> 2, "CNOT" -> {1, 3}, "CZ" -> {1, 2}, "V" -> 2}]]
], True, {}, TestID -> "BPB-VerifiedPrimitivesBypassGeneralClassifier"];

(* Check primitive actions against dense U P U-dagger, including signed images.
   Relabelled and reverse-wire operators must retain their actual action. *)
Scan[Function[spec,
    VerificationTest[With[{qc = bpbCircuit[{spec}]}, bpbImages[qc[PauliStabilizer[2]]]],
        bpbDenseImages[bpbCircuit[{spec}], 2], {},
        TestID -> "BPB-PrimitiveDense-" <> ToString[spec, InputForm]]
], Join[Thread[{"I", "H", "S", SuperDagger["S"], "X", "Y", "Z", "V", SuperDagger["V"]} -> ConstantArray[2, 9]],
    Thread[{"CNOT", "CZ", "SWAP"} -> ConstantArray[{1, 2}, 3]],
    Thread[{"CNOT", "CZ", "SWAP"} -> ConstantArray[{2, 1}, 3]]]];

SeedRandom[20261022];
Scan[Function[n,
    Scan[Function[trial,
        VerificationTest[With[{qc = bpbCircuit[bpbRandomSpecs[n, 24]]},
            bpbImages[qc[PauliStabilizer[n]]] === bpbDenseImages[qc, n]], True, {},
            TestID -> "BPB-RandomDense-n" <> ToString[n] <> "-" <> ToString[trial]]
    ], Range[3]]
], {2, 3, 4}];

VerificationTest[With[{ops = {
    QuantumOperator["H", {2}, "Label" -> "T"], QuantumOperator["CNOT", {2, 3}, "Label" -> "Y"],
    QuantumOperator["S", {1}, "Label" -> "H"], QuantumOperator["H", {3}, "Label" -> None],
    QuantumOperator["CZ", {3, 1}, "Label" -> "I"], QuantumOperator["S", {2}, "Label" -> "X"]
}}, With[{qc = bpbCircuit[ops]}, bpbImages[qc[PauliStabilizer[3]]] === bpbDenseImages[qc, 3]]],
    True, {}, TestID -> "BPB-ContentRecognitionIgnoresMisleadingLabels"];

VerificationTest[With[{specs = {"I" -> 1, "H" -> 2, "I" -> 3, "S" -> 1, "I" -> 2,
    "CNOT" -> {2, 3}, "CZ" -> {3, 1}, "H" -> 2, "I" -> 1}},
    bpbEqual[bpbCircuit[specs][PauliStabilizer[3]], PauliStabilizer[3]["ApplyCircuit", DeleteCases[specs, "I" -> _]]]],
    True, {}, TestID -> "BPB-InterleavedIdentityBlocksPreserveGateOrder"];

VerificationTest[bpbEqual[bpbCircuit[{"I" -> 1, "I" -> 3, "I" -> 2}][PauliStabilizer[3]], PauliStabilizer[3]],
    True, {}, TestID -> "BPB-AllIdentityGroups"];
VerificationTest[bpbEqual[bpbCircuit[{}][PauliStabilizer[3]], PauliStabilizer[3]],
    True, {}, TestID -> "BPB-EmptyCircuit"];

VerificationTest[With[{a = bpbCircuit[{"H" -> 1, "CNOT" -> {1, 3}}], b = bpbCircuit[{"S" -> 2, "CZ" -> {3, 2}}]},
    bpbEqual[bpbCircuit[{a, "Barrier", b, a}][PauliStabilizer[3]],
        PauliStabilizer[3]["ApplyCircuit", {"H" -> 1, "CNOT" -> {1, 3}, "S" -> 2, "CZ" -> {3, 2}, "H" -> 1, "CNOT" -> {1, 3}}]]],
    True, {}, TestID -> "BPB-NestedCircuitsAndBarriers"];

Scan[Function[n,
    VerificationTest[With[{specs = Join[{"H" -> n, "S" -> n, "CNOT" -> {n, 1}, "CZ" -> {1, n}}, bpbRandomSpecs[n, 60]]},
        With[{initial = PauliStabilizer[n]["ApplyCircuit", {"H" -> 1, "Y" -> n}]},
            bpbEqual[bpbCircuit[specs][initial], initial["ApplyCircuit", specs]]]],
        True, {}, TestID -> "BPB-PackedBoundary-" <> ToString[n]]
], {31, 32, 61, 62, 63, 124, 125, 128, 1024}];

VerificationTest[With[{specs = bpbRandomSpecs[128, 1000]},
    With[{qc = bpbCircuit[specs], initial = PauliStabilizer[128]},
        bpbEqual[Nest[qc, initial, 10], initial["ApplyCircuit", Catenate[ConstantArray[specs, 10]]]]]],
    True, {}, TestID -> "BPB-TenThousandGates-ReusedPlan"];

Scan[Function[precision,
    VerificationTest[With[{h = N[{{1, 1}, {1, -1}} / Sqrt[2], precision]},
        With[{qc = bpbCircuit[{QuantumOperator[h, {3}, "Label" -> "T"], "CNOT" -> {3, 1}, QuantumOperator[h, {2}]}]},
            bpbImages[qc[PauliStabilizer[3]]] === bpbDenseImages[bpbCircuit[{"H" -> 3, "CNOT" -> {3, 1}, "H" -> 2}], 3]]],
        True, {}, TestID -> "BPB-ApproximateCustomGates-" <> ToString[precision]]
], {MachinePrecision, 40}];

VerificationTest[With[{ops = {
    QuantumOperator[I {{1, 1}, {1, -1}} / Sqrt[2], {1}], QuantumOperator[-PauliMatrix[1], {3}],
    QuantumOperator[Exp[I Pi/7] DiagonalMatrix[{1, I}], {2}]
}}, With[{qc = bpbCircuit[ops]}, bpbImages[qc[PauliStabilizer[3]]] === bpbDenseImages[qc, 3]]],
    True, {}, TestID -> "BPB-ExactGlobalPhasesUseGeneralClassifier"];

VerificationTest[With[{specs = {"H" -> 1, "S" -> 2, "CNOT" -> {2, 3}, "V" -> 1, "CZ" -> {3, 1}, "Y" -> 2},
    initial = PauliStabilizer[<|"Tableau" -> PauliStabilizer[3]["Tableau"],
        "Signs" -> 1 - 2 {bpbA, bpbB, bpbC, bpbD, bpbE, bpbF}|>]},
    With[{actual = bpbCircuit[specs][initial], expected = initial["ApplyCircuit", specs]},
        actual["Tableau"] === expected["Tableau"] &&
            And @@ (With[{rules = Thread[{bpbA, bpbB, bpbC, bpbD, bpbE, bpbF} -> #]},
                (actual["Signs"] /. rules) === (expected["Signs"] /. rules)] & /@ Tuples[{0, 1}, 6])]],
    True, {}, TestID -> "BPB-SymbolicSigns-AllBinaryAssignments"];

VerificationTest[With[{specs = {"H" -> 1, "CNOT" -> {1, 2}, "T" -> 2, "V" -> 1, "S" -> 2}},
    With[{actual = bpbCircuit[specs][PauliStabilizer[2]], expected = PauliStabilizer[2]["ApplyCircuit", specs]},
        actual["State"]["StateVector"] === expected["State"]["StateVector"]]],
    True, {}, TestID -> "BPB-NonCliffordFrameFallback"];

VerificationTest[With[{p = QuantumOperator["P"[ArcCos[bpbTheta]]]},
    With[{qc = bpbCircuit[{p, QuantumOperator[p, {2}], "H" -> 1}]},
        Length[stabilizerCircuitSpecs[qc]] == 3 && FreeQ[stabilizerCircuitSpecs[qc], _Missing]]],
    True, {}, TestID -> "BPB-SymbolicTokensRemappedPerOccurrence"];

VerificationTest[With[{p = QuantumOperator["P"[ArcCos[bpbTheta]]]},
    bpbCircuit[{p, QuantumOperator[p, {2}, "Label" -> "H"]}][PauliStabilizer[2]]],
    $Failed, {PauliStabilizer::nonclifford}, TestID -> "BPB-SymbolicTokensDoNotBorrowAnotherLabel"];

VerificationTest[With[{op = QuantumOperator["T"] @ QuantumOperator["H"]},
    bpbCircuit[{op, QuantumOperator[op, "Label" -> "H"]}][PauliStabilizer[1]]],
    $Failed, {PauliStabilizer::nonclifford}, TestID -> "BPB-CompositeTokensDoNotBorrowAnotherLabel"];

VerificationTest[bpbCircuit[{QuantumOperator[PauliMatrix[1], {{1}, {1}}],
    QuantumOperator[PauliMatrix[1], {{2}, {1}}]}][PauliStabilizer[2]],
    $Failed, {PauliStabilizer::nonclifford}, TestID -> "BPB-SharedDataCannotHideMovedWires"];

VerificationTest[bpbCircuit[{QuantumOperator["H"], QuantumOperator[QuantumOperator["H"], QuantumBasis["PauliX"]]}][PauliStabilizer[1]],
    $Failed, {PauliStabilizer::nonclifford}, TestID -> "BPB-ComputationalAndNoncomputationalBasesRemainDistinct"];

VerificationTest[bpbCircuit[{QuantumOperator["X"[3]], QuantumOperator["X"]}][PauliStabilizer[1]],
    $Failed, {PauliStabilizer::nonclifford}, TestID -> "BPB-QuditDimensionsRemainDistinct"];

VerificationTest[bpbCircuit[{"H" -> 1, QuantumOperator[2 PauliMatrix[1]], "S" -> 1}][PauliStabilizer[1]],
    $Failed, {PauliStabilizer::nonclifford}, TestID -> "BPB-NonunitaryGateRefused"];

VerificationTest[With[{qc = bpbCircuit[{"H" -> 3, "CNOT" -> {3, 1}}]},
    bpbEqual[qc[Method -> "Stabilizer"], qc[PauliStabilizer[3]]]],
    True, {}, TestID -> "BPB-DefaultRegisterUsesHighestWire"];

VerificationTest[With[{specs = {"H" -> 1, "S" -> 2, "CNOT" -> {1, 3}, "CZ" -> {3, 2}}},
    With[{ps = PauliStabilizer[3]["ApplyCircuit", specs], qc = bpbCircuit[specs]},
        And @@ (bpbEqual[qc[#], #["ApplyCircuit", specs]] & /@ {
            PauliStabilizer[<|"Signs" -> ps["Signs"], "Tableau" -> ps["Tableau"]|>],
            psFromPacked[psGetPacked[ps]], ps})]],
    True, {}, TestID -> "BPB-AllStorageLayouts"];

EndTestSection[];
