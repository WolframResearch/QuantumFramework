BeginTestSection["ExecutionContracts"];

ecLayouts[ps_] := {
    PauliStabilizer[<|"Signs" -> ps["Signs"], "Tableau" -> ps["Tableau"]|>],
    psFromPacked[psGetPacked[ps]], ps
};

Scan[Function[layout,
    Scan[Function[gate,
        VerificationTest[layout[[2]][gate, 0], $Failed,
            {PauliStabilizer::invalidorder},
            TestID -> "EC-ScalarLowerBound-" <> layout[[1]] <> "-" <> ToString[gate, InputForm]]
    ], {"H", "S", SuperDagger["S"], "X", "Y", "Z", "V", SuperDagger["V"],
        "T", SuperDagger["T"], "P"[Pi/8]}];
    Scan[Function[gate,
        VerificationTest[layout[[2]][gate, 2, 2], $Failed,
            {PauliStabilizer::invalidorder},
            TestID -> "EC-ScalarDistinctWires-" <> layout[[1]] <> "-" <> gate]
    ], {"CNOT", "CX", "CZ", "SWAP"}];
    Scan[Function[wire,
        VerificationTest[layout[[2]]["H", wire], $Failed,
            {PauliStabilizer::invalidorder},
            TestID -> "EC-ScalarUpperAndExtreme-" <> layout[[1]] <> "-" <> ToString[wire]]
    ], {4, -2^63, 2^63}]
], Transpose[{{"Canonical", "Chunk", "Generator"}, ecLayouts[PauliStabilizer[3]]}]];

VerificationTest[Module[{calls = 0, invalid, outputs},
    invalid = {
        PauliStabilizer[<|"Signs" -> {1}, "Tableau" -> ConstantArray[0, {2, 1, 2}]|>],
        PauliStabilizer[<|"Signs" -> {1, 1}, "Tableau" -> ConstantArray[0, {1, 1, 2}]|>],
        PauliStabilizer[Append[First[PauliStabilizer[3]], "Signs" -> {}]],
        PauliStabilizer[Append[First[PauliStabilizer[3]], "GeneratorZ" -> {{2^62}}]],
        PauliStabilizer[Append[First[psFromPacked[psGetPacked[PauliStabilizer[3]]]], "Signs" -> {1}]]
    };
    outputs = Block[{compiledCliffordFold = Function[Null, calls++; $Failed]},
        {PauliStabilizerQ[#], psConcreteFastQ[#], PauliStabilizerQ[#["ApplyCircuit", {"H" -> 1}]]} & /@ invalid];
    calls === 0 && outputs === ConstantArray[{False, False, False}, Length[invalid]]],
    True, TestID -> "EC-InvalidRepresentationsCannotEnterNativeFold"];

VerificationTest[With[{tableau = PauliStabilizer[2]["Tableau"]},
    {PauliStabilizerQ[#], psConcreteFastQ[#]} & /@ {
        PauliStabilizer[<|"Signs" -> {1, 1, 1 - 2 ecOutcome, 1}, "Tableau" -> tableau|>],
        PauliStabilizer[<|"Signs" -> {1, 1, 1, 1}, "Tableau" -> -tableau|>]}],
    {{True, False}, {True, False}}, TestID -> "EC-SymbolicAndSignedTableauxKeepCanonicalFallback"];

VerificationTest[Module[{specs = Table[{"H" -> wire}, {wire, 20}], encoded},
    stabilizerClearCaches[];
    encoded = encodeStabilizerGates /@ specs;
    0 < Length[$stabilizerEncodeMemo] <= 16 && AllTrue[encoded, Developer`PackedArrayQ] &&
        encodeStabilizerGates[First[specs]] === {{0, 0, 0}}],
    True, TestID -> "EC-EncodingCacheBoundAndEvictedKeyRecompute"];

VerificationTest[Module[{circuits, plans},
    circuits = QuantumCircuitOperator[{"H" -> #}] & /@ Range[20];
    stabilizerClearCaches[];
    plans = stabilizerExecutionPlan /@ circuits;
    0 < Length[Wolfram`QuantumFramework`GateMatrix`PackagePrivate`$stabilizerCircuitMemo] <= 16 &&
        AllTrue[plans, AssociationQ] &&
        stabilizerExecutionPlan[First[circuits]] === First[plans]],
    True, TestID -> "EC-ExecutionPlanCacheBoundAndEvictedKeyRecompute"];

VerificationTest[Module[{h, x, plan},
    h = QuantumCircuitOperator[{"H" -> 1}]; x = QuantumCircuitOperator[{"X" -> 1}];
    stabilizerClearCaches[]; plan = stabilizerExecutionPlan[h];
    Wolfram`QuantumFramework`GateMatrix`PackagePrivate`$stabilizerCircuitMemo[Hash[x]] = {h, plan};
    x[PauliStabilizer[1]]["Stabilizers"]], {"-Z"},
    TestID -> "EC-ExecutionPlanHashCollisionCannotChangePhysicalGate"];

Scan[Function[wire,
    VerificationTest[PauliStabilizer[3]["V", wire]["Stabilizers"],
        ReplacePart[{"ZII", "IZI", "IIZ"}, wire -> "-" <>
            StringReplacePart["III", "Y", {wire, wire}]],
        TestID -> "EC-ScalarVGeneratorImages-wire" <> ToString[wire]];
    VerificationTest[PauliStabilizer[3][SuperDagger["V"], wire]["Stabilizers"],
        ReplacePart[{"ZII", "IZI", "IIZ"}, wire ->
            StringReplacePart["III", "Y", {wire, wire}]],
        TestID -> "EC-ScalarVAdjointGeneratorImages-wire" <> ToString[wire]]
], Range[3]];

EndTestSection[];
