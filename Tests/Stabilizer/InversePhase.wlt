BeginTestSection["InversePhase"];
SeedRandom[20261023];
ivSame[a_, b_] := a["Signs"] === b["Signs"] && Normal[a["Tableau"]] === Normal[b["Tableau"]];
ivKey[p_] := {Normal[p["Tableau"]], p["Signs"]};
ivCanonical[p_] := PauliStabilizer[<|"Signs" -> p["Signs"], "Tableau" -> Normal[p["Tableau"]]|>];
ivLayouts[p_] := <|"Canonical" -> ivCanonical[p], "Chunk" -> psFromPacked[psGetPacked[p]],
    "Resident" -> psToGeneratorPacked[p]|>;
ivStorage[p_] := Which[psGeneratorPackedQ[p], "Resident", psPackedQ[p], "Chunk", True, "Canonical"];
ivIdentity[n_] := PauliStabilizer[n];
ivLaws[p_] := With[{q = p["Dagger"], id = ivIdentity[p["Qubits"]]},
    ivSame[p[q], id] && ivSame[q[p], id] && ivSame[q["Dagger"], p] &&
    q["Dagger"] === p && p["Inverse"] === q && ivStorage[q] === ivStorage[p]];
ivAdjoint[g_ -> w_] := Replace[g, {"S" -> SuperDagger["S"], SuperDagger["S"] -> "S",
    "V" -> SuperDagger["V"], SuperDagger["V"] -> "V"}] -> w;
ivSpecs[1, count_] := Table[RandomChoice[{"H", "S", SuperDagger["S"], "X", "Y", "Z", "V", SuperDagger["V"]}] -> 1, {count}];
ivSpecs[n_, count_] := Table[With[{w = RandomSample[Range[n], 2],
    g = RandomChoice[{"H", "S", SuperDagger["S"], "X", "Y", "Z", "V", SuperDagger["V"], "CNOT", "CZ", "SWAP"}]},
    g -> If[MemberQ[{"CNOT", "CZ", "SWAP"}, g], w, First[w]]], {count}];

(* Exhaust the complete one-qubit Clifford action, including every sign pattern. *)
ivOneCliffords = FixedPoint[DeleteDuplicatesBy[Join[#, Catenate[({#["H", 1], #["S", 1]} &) /@ #]], ivKey] &,
    {PauliStabilizer[1]}];
VerificationTest[Length[ivOneCliffords], 24, TestID -> "Inverse-All24SingleQubitCliffordsEnumerated"];
MapIndexed[Function[{p, index}, KeyValueMap[Function[{layout, value},
    VerificationTest[ivLaws[value], True,
        TestID -> "Inverse-All24-" <> ToString[First[index]] <> "-" <> layout]
], ivLayouts[p]]], ivOneCliffords];

Scan[Function[g, VerificationTest[
    ivSame[PauliStabilizer[1]["ApplyCircuit", {g -> 1}]["Dagger"], PauliStabilizer[1]["ApplyCircuit", {ivAdjoint[g -> 1]}]],
    True, TestID -> "Inverse-PrimitiveAdjoint-" <> ToString[g, InputForm]]],
    {"H", "S", SuperDagger["S"], "X", "Y", "Z", "V", SuperDagger["V"]}];

ivRegression = PauliStabilizer[<|"Signs" -> {1, 1, 1, -1, 1, 1},
    "Tableau" -> {{{1, 1, 0, 1, 1, 1}, {1, 1, 0, 1, 1, 0}, {1, 0, 0, 1, 0, 1}},
        {{1, 1, 1, 0, 0, 0}, {0, 0, 1, 1, 0, 1}, {0, 0, 0, 1, 1, 0}}}|>];
KeyValueMap[Function[{layout, p}, VerificationTest[ivLaws[p], True,
    TestID -> "Inverse-ReportedPhaseRegression-" <> layout]], ivLayouts[ivRegression]];

(* Test every possible input sign vector for the reported three-qubit matrix. *)
Scan[Function[bits, KeyValueMap[Function[{layout, p}, VerificationTest[ivLaws[p], True,
    TestID -> "Inverse-All64SignVectors-" <> StringJoin[ToString /@ bits] <> "-" <> layout]],
    ivLayouts[PauliStabilizer[<|"Matrix" -> ivRegression["Matrix"], "Signs" -> 1 - 2 bits|>]]]], Tuples[{0, 1}, 6]];

Scan[Function[n, Scan[Function[trial,
    VerificationTest[With[{specs = ivSpecs[n, 60]},
        With[{p = PauliStabilizer[n]["ApplyCircuit", specs],
            expected = PauliStabilizer[n]["ApplyCircuit", ivAdjoint /@ Reverse[specs]]},
            ivSame[p["Dagger"], expected] && ivLaws[p]]], True,
        TestID -> "Inverse-ReversedRandomCircuit-n" <> ToString[n] <> "-" <> ToString[trial]]
], Range[8]]], {2, 3, 4, 8}];

(* Independent System-matrix reference: no tableau synthesis or composition. *)
ivTensor[matrices_] := Fold[KroneckerProduct, {{1}}, matrices];
ivEmbed[a_, q_, n_] := ivTensor[ReplacePart[ConstantArray[IdentityMatrix[2], n], q -> a]];
ivGate[g_ -> w_, n_] := With[{wires = Flatten[{w}]},
    If[MemberQ[{"CNOT", "CZ", "SWAP"}, g],
        With[{q = wires[[1]], r = wires[[2]], bits = IntegerDigits[#, 2, n] & /@ Range[0, 2^n - 1]},
            If[g === "CZ", DiagonalMatrix[(-1)^(#[[q]] #[[r]]) & /@ bits],
                With[{outputs = If[g === "CNOT",
                    (ReplacePart[#, r -> Mod[#[[r]] + #[[q]], 2]] &) /@ bits,
                    (ReplacePart[#, {q -> #[[r]], r -> #[[q]]}] &) /@ bits]},
                    Normal[SparseArray[Thread[Transpose[{(FromDigits[#, 2] + 1 &) /@ outputs, Range[2^n]}] -> 1], {2^n, 2^n}]]]
            ]
        ],
        ivEmbed[Replace[g, {"H" -> {{1, 1}, {1, -1}}/Sqrt[2], "S" -> {{1, 0}, {0, I}},
            SuperDagger["S"] -> {{1, 0}, {0, -I}}, "X" -> PauliMatrix[1], "Y" -> PauliMatrix[2], "Z" -> PauliMatrix[3],
            "V" -> {{1 + I, 1 - I}, {1 - I, 1 + I}}/2,
            SuperDagger["V"] -> {{1 - I, 1 + I}, {1 + I, 1 - I}}/2}], First[wires], n]
    ]
];
ivRow[p_, j_] := With[{t = p["Tableau"]}, p["Signs"][[j]] ivTensor[
    Replace[Transpose[{t[[1, All, j]], t[[2, All, j]]}],
        {{0, 0} -> PauliMatrix[0], {1, 0} -> PauliMatrix[1], {1, 1} -> PauliMatrix[2], {0, 1} -> PauliMatrix[3]}, {1}]]];
Scan[Function[n, Scan[Function[trial, VerificationTest[
    With[{specs = ivSpecs[n, 25]}, With[{p = PauliStabilizer[n]["ApplyCircuit", specs],
        u = Fold[ivGate[#2, n] . #1 &, IdentityMatrix[2^n], specs]},
        With[{inverse = p["Dagger"]}, Simplify[ConjugateTranspose[u] . u] === IdentityMatrix[2^n] &&
            And @@ Map[Function[j, With[{generator = ivEmbed[PauliMatrix[If[j <= n, 1, 3]], Mod[j - 1, n] + 1, n]},
                Simplify[ivRow[inverse, j] - ConjugateTranspose[u] . generator . u] === ConstantArray[0, {2^n, 2^n}]]], Range[2 n]]]]],
    True, TestID -> "Inverse-DenseAllGeneratorConjugation-n" <> ToString[n] <> "-" <> ToString[trial]]
], Range[3]]], {1, 2, 3, 4}];

Scan[Function[n, VerificationTest[With[{p = PauliStabilizer[n]["ApplyCircuit", ivSpecs[n, 100]]},
    With[{q = p["Dagger"], id = IdentityMatrix[2 n]},
        ivLaws[p] && Mod[p["Matrix"] . q["Matrix"], 2] === id &&
        KeyExistsQ[First[q], "GeneratorX"]]], True,
    TestID -> "Inverse-PackedBoundary-n" <> ToString[n]]], {31, 32, 61, 62, 63, 124, 125}];

VerificationTest[Module[{base, symbolic, inverse, variables = Array[ivBit, 4]},
    base = PauliStabilizer[2]["ApplyCircuit", {"H" -> 1, "S" -> 1, "CNOT" -> {1, 2}, "V" -> 2}];
    symbolic = PauliStabilizer[<|"Matrix" -> base["Matrix"], "Signs" -> 1 - 2 variables|>];
    inverse = symbolic["Dagger"];
    And @@ Map[Function[bits, With[{rules = Thread[variables -> bits]},
        (* Association replacement can retain held arithmetic; evaluate properties
           when constructing the concrete comparison tableau. *)
        With[{concrete = ivCanonical[symbolic /. rules]},
            ivSame[inverse /. rules, concrete["Dagger"]] && ivLaws[concrete]]]], Tuples[{0, 1}, 4]]],
    True, TestID -> "Inverse-SymbolicPhaseSubstitution-All16Assignments"];

Scan[Function[alias, VerificationTest[PauliStabilizer[{"XX", "XX"}][alias], $Failed,
    {PauliStabilizer::singular}, TestID -> "Inverse-SingularInput-" <> alias]], {"Dagger", "Inverse"}];
EndTestSection[];
