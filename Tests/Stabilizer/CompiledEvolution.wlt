BeginTestSection["CompiledEvolution"];
SeedRandom[20261022];
ceSame[a_, b_] := a["Signs"] === b["Signs"] && a["Tableau"] === b["Tableau"];
ceCanonical[ps_] := PauliStabilizer[<|"Signs" -> ps["Signs"], "Tableau" -> ps["Tableau"]|>];
ceSpecs[1, count_] := Table[RandomChoice[{"H", "S", SuperDagger["S"], "X", "Y", "Z"}] -> 1, {count}];
ceSpecs[n_, count_] := Table[With[{w = RandomSample[Range[n], 2],
    g = RandomChoice[{"H", "S", SuperDagger["S"], "X", "Y", "Z", "CNOT", "CZ", "SWAP"}]},
    g -> If[MemberQ[{"CNOT", "CZ", "SWAP"}, g], w, First[w]]], {count}];
ceOmega[n_] := ArrayFlatten[{{ConstantArray[0, {n, n}], IdentityMatrix[n]},
    {IdentityMatrix[n], ConstantArray[0, {n, n}]}}];

Scan[Function[n,
    VerificationTest[With[{ps = PauliStabilizer[n], specs = ceSpecs[n, 160]},
        With[{bulk = ps["ApplyCircuit", specs], scalar = Fold[#1[#2] &, ceCanonical[ps], specs]},
            ceSame[bulk, scalar] && PauliStabilizerQ[bulk] &&
            KeyExistsQ[First[bulk], "GeneratorX"] && !KeyExistsQ[First[bulk], "Tableau"]]],
        True, TestID -> "CE-ExactScalarAndResident-n" <> ToString[n]];
    VerificationTest[With[{ps = PauliStabilizer[n], specs = ceSpecs[n, 120]},
        ceSame[ps["ApplyCircuit", specs], Fold[#1["ApplyCircuit", #2] &, ps, Partition[specs, 20]]]],
        True, TestID -> "CE-SplitCircuitExact-n" <> ToString[n]];
    VerificationTest[With[{ps = PauliStabilizer[n]["ApplyCircuit", ceSpecs[n, 100]], j = ceOmega[n]},
        Mod[ps["Matrix"] . j . Transpose[ps["Matrix"]], 2] === j],
        True, TestID -> "CE-Symplectic-n" <> ToString[n]]
], {1, 2, 5, 31, 32, 61, 62, 63, 124, 125}];

VerificationTest[With[{ps = PauliStabilizer[1024], specs = ceSpecs[1024, 1000]},
    ceSame[ps["ApplyCircuit", specs], Fold[#1[#2] &, ceCanonical[ps], specs]]],
    True, TestID -> "CE-LargeExactScalar-n1024"];

VerificationTest[Module[{ps = ceCanonical[PauliStabilizer[32]], before, after, bad},
    before = With[{value = ps}, System`Private`HoldValidQ[value]];
    PauliStabilizerQ[ps]; after = With[{value = ps}, System`Private`HoldValidQ[value]];
    bad = ps; bad[[1, "Signs"]] = {};
    {before, after, PauliStabilizerQ[bad]}],
    {True, True, False}, TestID -> "CE-ConstructionValidityMarkerAndMutation"];

VerificationTest[Module[{ps = PauliStabilizer[32], bad},
    bad = PauliStabilizer[Append[First[ps], "Signs" -> {}]];
    {PauliStabilizerQ[ps], PauliStabilizerQ[bad]}], {True, False},
    TestID -> "CE-NewAssociationDoesNotInheritValidity"];

VerificationTest[With[{a = First[PauliStabilizer[32]]},
    And @@ (Not[PauliStabilizerQ[PauliStabilizer[#]]] & /@ {
        KeyDrop[a, "GeneratorZ"], Append[a, "GeneratorX" -> {{1}}],
        Append[a, "GeneratorX" -> ConstantArray[-1, {32, 2}]],
        Append[a, "GeneratorZ" -> ConstantArray[2^62, {32, 2}]],
        Append[a, "GeneratorZ" -> ConstantArray[4, {32, 2}]]})],
    True, TestID -> "CE-MalformedGeneratorPackedRejected"];

VerificationTest[With[{ps = psFromPacked[psGetPacked[PauliStabilizer[63]]]},
    PauliStabilizerQ[PauliStabilizer[Append[First[ps], "PackedZ" -> {{1}}]]]],
    False, TestID -> "CE-MalformedChunkPackedRejected"];

Scan[Function[g, Scan[Function[q,
    VerificationTest[With[{ps = PauliStabilizer[4]["ApplyCircuit", ceSpecs[4, 40]]},
        ceSame[ps["ApplyCircuit", {g -> q}], Fold[#1[#2] &, ceCanonical[ps],
            {"H" -> q, If[g === "V", "S", SuperDagger["S"]] -> q, "H" -> q}]]],
        True, TestID -> "CE-V-CliffordDecomposition-" <> ToString[{g, q}, InputForm]]
], Range[4]]], {"V", SuperDagger["V"]}];

Scan[Function[spec,
    VerificationTest[PauliStabilizer[3]["ApplyCircuit", {spec}], $Failed,
        {PauliStabilizer::invalidorder}, TestID -> "CE-BadWire-" <> ToString[spec, InputForm]]
], {"H" -> 0, "H" -> 4, "S" -> -1, "CNOT" -> {1, 1}, "SWAP" -> {2, 2}, "CZ" -> {2, 5}}];

VerificationTest[With[{ps = PauliStabilizer[5]}, ps["ApplyCircuit", {}] === ps],
    True, TestID -> "CE-EmptyIdentityPreservesObject"];

VerificationTest[Module[{a = {"H" -> 1}, b = {"S" -> 2}, expected},
    stabilizerClearCaches[]; expected = encodeStabilizerGates[b];
    $stabilizerEncodeMemo[Hash[b]] = {a, encodeStabilizerGates[a]};
    encodeStabilizerGates[b] === expected], True,
    TestID -> "CE-EncodingCacheCollisionChecksOriginalKey"];

VerificationTest[Module[{qc, first, second, count1, count2},
    qc = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 3}, "S" -> 2}];
    stabilizerClearCaches[];
    first = stabilizerExecutionPlan[qc]; count1 = Length[$stabilizerEncodeMemo];
    second = stabilizerExecutionPlan[qc]; count2 = Length[$stabilizerEncodeMemo];
    first === second && count1 === count2 && ceSame[qc[PauliStabilizer[3]],
        PauliStabilizer[3]["ApplyCircuit", first["Specs"]]]], True,
    TestID -> "CE-ValidatedCircuitPlanReuse"];

VerificationTest[Module[{h, x, a, b},
    (* Use normalized H and an X with the same display label. *)
    h = QuantumCircuitOperator[{QuantumOperator[(PauliMatrix[1] + PauliMatrix[3])/Sqrt[2], "Label" -> "X"]}];
    x = QuantumCircuitOperator[{QuantumOperator[PauliMatrix[1], "Label" -> "X"]}];
    stabilizerClearCaches[]; a = h[PauliStabilizer[1]]; b = x[PauliStabilizer[1]];
    {a["Stabilizers"], b["Stabilizers"]}], {{"X"}, {"-Z"}},
    TestID -> "CE-PlanCacheKeepsMatrixHonesty"];

VerificationTest[Module[{q1, q3},
    q1 = QuantumCircuitOperator[{"H" -> 1}]; q3 = QuantumCircuitOperator[{"H" -> 3}];
    stabilizerClearCaches[]; q1[PauliStabilizer[3]];
    q3[PauliStabilizer[3]]["Stabilizers"]], {"ZII", "IZI", "IIX"},
    TestID -> "CE-PlanCacheKeepsWirePlacement"];

VerificationTest[Module[{qc, s, t},
    qc = QuantumCircuitOperator[{"P"[ceAngle] -> 1}, "Parameters" -> {ceAngle}];
    s = qc[<|ceAngle -> Pi/2|>]; t = qc[<|ceAngle -> Pi/4|>];
    stabilizerClearCaches[];
    {Head[s[PauliStabilizer[1]]], Head[t[PauliStabilizer[1]]]}],
    {PauliStabilizer, StabilizerFrame}, TestID -> "CE-PlanCacheKeepsParameterValues"];

VerificationTest[Module[{ps, sym, evolved, concrete},
    ps = PauliStabilizer[2]["ApplyCircuit", {"H" -> 1, "CNOT" -> {1, 2}}];
    sym = ps["SymbolicMeasure", 1]; evolved = sym["ApplyCircuit", {"X" -> 2, "H" -> 1}];
    concrete = evolved["SampleOutcomes", 10];
    AllTrue[concrete, PauliStabilizerQ]], True,
    TestID -> "CE-SymbolicSignsRetainFallback"];

VerificationTest[With[{ps = PauliStabilizer[3]["ApplyCircuit", ceSpecs[3, 50]]},
    With[{canonical = ceCanonical[ps]},
        KeyExistsQ[First[ps["Dagger"]], "GeneratorX"] &&
        ceSame[ps["Dagger"], canonical["Dagger"]] &&
        ceSame[ps["Dagger"]["Dagger"], ps] && ps["Dagger"]["Dagger"] === ps]], True,
    TestID -> "CE-InversePreservesResidentRepresentation"];

VerificationTest[Module[{ps = PauliStabilizer[3]["ApplyCircuit", {"H" -> 1, "CNOT" -> {1, 3}, "S" -> 3}], canonical},
    canonical = ceCanonical[ps];
    ps["Stabilizers"] === canonical["Stabilizers"] && ps["TableauPhase"] === canonical["TableauPhase"] &&
    ps["Entropy", {1}] === canonical["Entropy", {1}] &&
    ps["Expectation", "XIZ"] === canonical["Expectation", "XIZ"] &&
    Keys[ps["M", 1]] === Keys[canonical["M", 1]]], True,
    TestID -> "CE-ResidentPropertiesAndMeasurement"];

VerificationTest[With[{ps = PauliStabilizer[3]["ApplyCircuit", ceSpecs[3, 30]]},
    With[{a = Normal[ps["State"]["StateVector"]], b = Normal[ceCanonical[ps]["State"]["StateVector"]]},
        Simplify[(Conjugate[a] . b) (Conjugate[b] . a)] === 1]], True,
    TestID -> "CE-ResidentStateMaterialization"];
VerificationTest[Module[{specs, expected, first},
    specs = {"H" -> 1, "CNOT" -> {2, 1}, "V" -> 3,
        SuperDagger["V"] -> {1}, "CZ" -> {1, 2}, "SWAP" -> {3, 1}};
    expected = {{0, 0, 0}, {3, 1, 0}, {0, 2, 0}, {1, 2, 0}, {0, 2, 0},
        {0, 0, 0}, {2, 0, 0}, {0, 0, 0}, {0, 1, 0}, {3, 0, 1}, {0, 1, 0}, {4, 2, 0}};
    stabilizerClearCaches[]; first = encodeStabilizerGates[specs];
    first === expected && Developer`PackedArrayQ[first] && encodeStabilizerGates[specs] === first],
    True, TestID -> "CE-VectorizedEncodingZeroBasedRows"];

VerificationTest[encodeStabilizerGates[{{"V", 2}, "C"["NOT" -> {3}, {1}], "C"["Z" -> {1}, {3}]}],
    {{0, 1, 0}, {1, 1, 0}, {0, 1, 0}, {3, 0, 2}, {0, 0, 0}, {3, 2, 0}, {0, 0, 0}},
    TestID -> "CE-VectorizedEncodingLegacyAndControlledForms"];

VerificationTest[{encodeStabilizerGates[{}], encodeStabilizerGates[{"H" -> 1, "T" -> 2}]},
    {{}, $Failed}, TestID -> "CE-VectorizedEncodingEmptyAndFallback"];

VerificationTest[encodeStabilizerGates[{"H" -> 2^63, "CNOT" -> {-2^63, 2^63}}],
    {{0, 2^63 - 1, 0}, {3, -2^63 - 1, 2^63 - 1}},
    TestID -> "CE-VectorizedEncodingExactIntegerOverflowBoundary"];

Scan[Function[wire,
    VerificationTest[PauliStabilizer[3]["ApplyCircuit", {"H" -> wire}], $Failed,
        {PauliStabilizer::invalidorder}, TestID -> "CE-VectorizedEncodingRejectsExtremeWire-" <> ToString[wire]]
], {-2^63, 2^63}];

EndTestSection[];
