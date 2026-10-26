(* ==========================================================================
   Tests/Stabilizer/MidCircuitMeasurement.wlt -- Method -> "Stabilizer" on
   circuits with computational-basis measurements, exact mode.

   qc[..., Method -> "Stabilizer"] returns <|record -> PauliStabilizer|>: every
   record the circuit can produce with its normalized post-measurement tableau,
   each with probability 1 / Length. Records follow the default method's
   outcome order: sorted by measured wire, time order within one wire.

   Every result is compared as <|record -> {probability, density matrix}|>,
   the density matrix normalized and read in the computational basis on wires
   1..n. Closed-form cases are checked against hand-derived kets, and the
   default method and Method -> "Schrodinger" are checked against the same
   kets. Random circuits are checked against the default method.
   ========================================================================== *)

Get[FileNameJoin[{ParentDirectory[PacletObject["Wolfram/QuantumFramework"]["Location"]], "Tests", "Stabilizer", "fixtures", "MeasurementHelpers.wl"}]]

stabMidM = QuantumMeasurementOperator;

stabMidDensity[v_List] := With[{u = v / Norm[v]}, Outer[Times, u, Conjugate[u]]]

(* Records of a stabilizer result: probability 1 / Length, the tableau's state. *)
stabMidFromTableaux[result_Association] := Map[{1 / Length[result], stabMeasTableauDensity[#]} &, result]

stabMidBits[QuditName[x_, ___]] := Flatten[{x}]

(* Records of a dense result with nonzero probability. A branch state is       *)
(* unnormalized: its density matrix has trace p, the record's probability.     *)
(* ("ProbabilitiesList" is not read: with 7 or more record bits it calls       *)
(* Eigensystem on a sparse matrix and emits Eigensystem::arhm.)                 *)
stabMidFromDense[qm_QuantumMeasurement] := Association @ MapThread[
    Function[{name, st},
        With[{rho = Normal[st["Computational"]["DensityMatrix"]]},
            With[{p = Tr[rho]}, If[PossibleZeroQ[p], Nothing, stabMidBits[name] -> {p, rho / p}]]
        ]
    ],
    {qm["Outcomes"], qm["States"]}
]

stabMidExpected[rules_List] := Association @ Map[First[#] -> {#[[2]], stabMidDensity[Last[#]]} &, rules]

(* Same records, exactly equal probabilities and density matrices. *)
stabMidSameQ[a_Association, b_Association] :=
    Sort[Keys[a]] === Sort[Keys[b]] && AllTrue[Keys[a], stabMeasExactSameQ[a[#], b[#]] &]

(* The three routes against expected kets: stabilizer, default, Schrodinger. *)
stabMidClosedForm[ops_List, input_QuantumState, expected_List] := With[{qc = QuantumCircuitOperator[ops], e = stabMidExpected[expected]}, {
    stabMidSameQ[stabMidFromTableaux[qc[input, Method -> "Stabilizer"]], e],
    stabMidSameQ[stabMidFromDense[qc[input]], e],
    stabMidSameQ[stabMidFromDense[qc[input, Method -> "Schrodinger"]], e]
}]

stabMidFailureReason[f_Failure] := If[f[[1]] === "StabilizerMeasurement", f["Reason"], f]
stabMidFailureReason[other_] := other


(* ---------- closed forms ---------- *)

VerificationTest[
    stabMidClosedForm[{stabMidM[{1}]}, QuantumState["1"], {{{1}, 1, {0, 1}}}],
    {True, True, True},
    TestID -> "MidMeasurement-Z-On-One"
]

VerificationTest[
    stabMidClosedForm[{"H" -> 1, "CNOT" -> {1, 2}, stabMidM[{1}], stabMidM[{2}]}, QuantumState["00"],
        {{{0, 0}, 1/2, {1, 0, 0, 0}}, {{1, 1}, 1/2, {0, 0, 0, 1}}}],
    {True, True, True},
    TestID -> "MidMeasurement-Bell"
]

VerificationTest[
    stabMidClosedForm[{"H" -> 1, stabMidM[{1}], stabMidM[{1}]}, QuantumState["0"],
        {{{0, 0}, 1/2, {1, 0}}, {{1, 1}, 1/2, {0, 1}}}],
    {True, True, True},
    TestID -> "MidMeasurement-RepeatedZ"
]

(* the gate after the first measurement must act: records disagree *)
VerificationTest[
    stabMidClosedForm[{"H" -> 1, stabMidM[{1}], "X" -> 1, stabMidM[{1}]}, QuantumState["0"],
        {{{0, 1}, 1/2, {0, 1}}, {{1, 0}, 1/2, {1, 0}}}],
    {True, True, True},
    TestID -> "MidMeasurement-GateAfterMeasurement"
]

(* four records, two distinct final states, all four kept *)
VerificationTest[
    stabMidClosedForm[{"H" -> 1, stabMidM[{1}], "H" -> 1, stabMidM[{1}]}, QuantumState["0"],
        {{{0, 0}, 1/4, {1, 0}}, {{0, 1}, 1/4, {0, 1}}, {{1, 0}, 1/4, {1, 0}}, {{1, 1}, 1/4, {0, 1}}}],
    {True, True, True},
    TestID -> "MidMeasurement-HMHM"
]

VerificationTest[
    With[{r = QuantumCircuitOperator[{"H" -> 1, stabMidM[{1}], "H" -> 1, stabMidM[{1}]}][Method -> "Stabilizer"]},
        {Length[r], Length[DeleteDuplicates[{#["Tableau"], #["Signs"]} & /@ Values[r]]]}
    ],
    {4, 2},
    TestID -> "MidMeasurement-IdenticalTableauxStaySeparate"
]

(* one instrument on wires {3, 1}: bits sorted by wire. Method -> "Schrodinger" *)
(* is left out: it does not apply an instrument on several wires (it returns    *)
(* the circuit unevaluated or a Failure).                                       *)
VerificationTest[
    Most @ stabMidClosedForm[{stabMidM[{3, 1}]}, QuantumState["100"], {{{1, 0}, 1, {0, 0, 0, 0, 1, 0, 0, 0}}}],
    {True, True},
    TestID -> "MidMeasurement-MultiWireOrder"
]

(* wires 3, 1, X on 3, 3, 2, 4 on |0110>: by wire, time order within wire 3 *)
VerificationTest[
    stabMidClosedForm[{stabMidM[{3}], stabMidM[{1}], "X" -> 3, stabMidM[{3}], stabMidM[{2}], stabMidM[{4}]}, QuantumState["0110"],
        {{{0, 1, 1, 0, 0}, 1, UnitVector[16, 5]}}],
    {True, True, True},
    TestID -> "MidMeasurement-RecordOrder"
]

VerificationTest[
    stabMidClosedForm[{"H" -> 1, "CNOT" -> {1, 2}, "CNOT" -> {2, 3}, stabMidM[{2}]}, QuantumState["000"],
        {{{0}, 1/2, UnitVector[8, 1]}, {{1}, 1/2, UnitVector[8, 8]}}],
    {True, True, True},
    TestID -> "MidMeasurement-GHZ3-MiddleWire"
]

(* a relative phase that survives the measurement *)
VerificationTest[
    stabMidClosedForm[{"H" -> 1, "H" -> 2, "CZ" -> {1, 2}, stabMidM[{1}], "S" -> 2}, QuantumState["00"],
        {{{0}, 1/2, {1, I, 0, 0} / Sqrt[2]}, {{1}, 1/2, {0, 0, 1, - I} / Sqrt[2]}}],
    {True, True, True},
    TestID -> "MidMeasurement-RelativePhase"
]

(* record key {wire-1 bit, wire-2 bit} though wire 2 is measured first *)
VerificationTest[
    stabMidClosedForm[{"H" -> 1, "CNOT" -> {1, 2}, stabMidM[{2}], "H" -> 1, stabMidM[{1}]}, QuantumState["00"],
        {{{0, 0}, 1/4, {1, 0, 0, 0}}, {{0, 1}, 1/4, {0, 1, 0, 0}}, {{1, 0}, 1/4, {0, 0, 1, 0}}, {{1, 1}, 1/4, {0, 0, 0, 1}}}],
    {True, True, True},
    TestID -> "MidMeasurement-WireOrderAcrossMeasurements"
]


(* ---------- random circuits against the default method ---------- *)

stabMidRandomMeasurement[n_] := If[RandomReal[] < 0.3, stabMidM[RandomSample[Range[n], 2]], stabMidM[{RandomInteger[{1, n}]}]]

(* {n, preparation gates, body with 1 to 4 measurements} *)
stabMidRandomCase[seed_] := BlockRandom[
    SeedRandom[seed];
    With[{n = RandomInteger[{3, 6}]},
        {
            n,
            Table[stabMeasRandomGate[n], 12],
            Catenate @ Table[Append[Table[stabMeasRandomGate[n], RandomInteger[{0, 4}]], stabMidRandomMeasurement[n]], RandomInteger[{1, 4}]]
        }
    ]
]

(* For one case: the three input forms agree with each other and with the default *)
(* method; all branches share one set of X/Z bits; the count is a power of two.   *)
stabMidRandomCheck[seed_] := Replace[stabMidRandomCase[seed], {n_, prep_, ops_} :>
    With[{ps = PauliStabilizer[n]["ApplyCircuit", prep], body = QuantumCircuitOperator[ops]},
        With[{
            fromDefault = QuantumCircuitOperator[Join[prep, ops]][Method -> "Stabilizer"],
            fromTableau = body[ps, Method -> "Stabilizer"],
            fromState = body[ps["State"], Method -> "Stabilizer"],
            dense = stabMidFromDense[body[ps["State"]]]
        },
            AllTrue[{fromDefault, fromTableau, fromState}, AssociationQ] &&
            AllTrue[{fromDefault, fromTableau, fromState}, stabMidSameQ[stabMidFromTableaux[#], dense] &] &&
            Length[DeleteDuplicates[#["Tableau"] & /@ Values[fromTableau]]] == 1
        ]
    ]
]

VerificationTest[
    Select[Range[50], ! stabMidRandomCheck[#] &],
    {},
    TestID -> "MidMeasurement-RandomCircuits"
]


(* ---------- contract ---------- *)

(* a measurement beyond the input register pads it with |0> *)
VerificationTest[
    With[{r = QuantumCircuitOperator[{"H" -> 3, stabMidM[{3}]}][PauliStabilizer[1], Method -> "Stabilizer"]},
        {Keys[r], #["Qubits"] & /@ Values[r]}
    ],
    {{{0}, {1}}, {3, 3}},
    TestID -> "MidMeasurement-Padding"
]

VerificationTest[
    Head[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}][Method -> "Stabilizer"]],
    PauliStabilizer,
    TestID -> "MidMeasurement-MeasurementFreeUnchanged"
]

(* branch count follows the random measurements only *)
(* for record b, the state is fixed by Z_i Z_(i+1) = +1 for all i and Z_1 = (-1)^b *)
stabMidGHZStateQ[ps_PauliStabilizer, b_] :=
    AllTrue[Range[49], Keys[ps["M", stabMeasZ[50, {#, # + 1}]]] === {0} &] && Keys[ps["M", stabMeasZ[50, {1}]]] === {b}

VerificationTest[
    {
        With[{r = QuantumCircuitOperator[Join[stabMeasGHZ[50], stabMidM[{#}] & /@ Range[50]]][Method -> {"Stabilizer", "MaxBranches" -> 2}]},
            {Keys[r], stabMidGHZStateQ[r[ConstantArray[0, 50]], 0], stabMidGHZStateQ[r[ConstantArray[1, 50]], 1]}],
        Length[QuantumCircuitOperator[Join[stabMeasGHZ[50], stabMidM[{#}] & /@ Range[50], stabMidM[{#}] & /@ Range[50]]][Method -> {"Stabilizer", "MaxBranches" -> 2}]]
    },
    {{{ConstantArray[0, 50], ConstantArray[1, 50]}, True, True}, 2},
    TestID -> "MidMeasurement-CertainMeasurementsDoNotCount"
]

(* three random measurements with certain ones between them *)
stabMidRandomCount[j_] := Catenate @ Table[{"H" -> i, stabMidM[{i}], stabMidM[{i}], stabMidM[{1}]}, {i, j}]

VerificationTest[
    {
        Length[QuantumCircuitOperator[stabMidRandomCount[3]][Method -> {"Stabilizer", "MaxBranches" -> 8}]],
        With[{f = QuantumCircuitOperator[stabMidRandomCount[4]][Method -> {"Stabilizer", "MaxBranches" -> 8}]},
            {stabMidFailureReason[f], f["Count"]}]
    },
    {8, {"MaxBranches", 16}},
    TestID -> "MidMeasurement-MaxBranchesBoundary"
]


(* ---------- refusals ---------- *)

VerificationTest[
    stabMidFailureReason /@ {
        QuantumCircuitOperator[{"T" -> 1, stabMidM[{1}]}][Method -> "Stabilizer"],
        QuantumCircuitOperator[{"H" -> 1, stabMidM[{1}], "T" -> 1}][Method -> "Stabilizer"],
        QuantumCircuitOperator[{"H" -> 1, QuantumChannel["BitFlip"[1/3]], stabMidM[{1}]}][Method -> "Stabilizer"],
        QuantumCircuitOperator[{"H" -> 1, stabMidM["X", {1}]}][Method -> "Stabilizer"],
        QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, stabMidM["Lueders"[QuantumOperator["ZZ"]], {1, 2}]}][Method -> "Stabilizer"],
        QuantumCircuitOperator[{"H" -> 1, stabMidM[{1}], QuantumOperator["C"["X" -> {2}, {0}]]}][Method -> "Stabilizer"],
        QuantumCircuitOperator[{stabMidM[{1}]}][QuantumState[{1, Exp[I Pi / 4]} / Sqrt[2]], Method -> "Stabilizer"],
        QuantumCircuitOperator[{stabMidM[{1}]}][QuantumState[IdentityMatrix[2] / 2], Method -> "Stabilizer"],
        QuantumCircuitOperator[{stabMidM[{1}]}][PauliStabilizer[1]["H" -> 1]["SymbolicMeasure", 1], Method -> "Stabilizer"],
        QuantumCircuitOperator[{"H" -> 1, stabMidM[{1}]}][Method -> {"Stabilizer", "MaxBranches" -> 0}],
        QuantumCircuitOperator[{"H" -> 1, stabMidM[{1}]}][Method -> {"Stabilizer", "Foo" -> 1}],
        QuantumCircuitOperator[{"H" -> 1, stabMidM[{1}]}][Method -> {"Stabilizer", "Compress" -> "PhasePolynomial"}],
        QuantumCircuitOperator[{stabMidM[QuantumBasis[3], {1}]}][Method -> "Stabilizer"],
        QuantumCircuitOperator[{stabMidM[QuantumBasis[3], {1}]}][QuantumState[{1, 0, 0}], Method -> "Stabilizer"]
    },
    {
        "NotClifford", "NotClifford", "NotClifford",
        "NotComputationalMeasurement", "NotComputationalMeasurement",
        "RecordWire",
        "NotStabilizerState", "NotStabilizerState", "SymbolicSigns",
        "InvalidOption", "UnknownOption", "InvalidOption",
        "NotQubit", "NotQubit"
    },
    TestID -> "MidMeasurement-Refusals"
]

(* an instrument on more wires than the matrix cap is refused before its operators *)
(* are read, whether or not the same instrument was classified before, and is      *)
(* accepted again once the cap allows it                                            *)
VerificationTest[
    With[{qc = QuantumCircuitOperator[{stabMidM[{5, 7}]}]},
        {
            AssociationQ[qc[Method -> "Stabilizer"]],
            Block[{Wolfram`QuantumFramework`PackageScope`$stabilizerMatrixMaxQubits = 1}, stabMidFailureReason[qc[Method -> "Stabilizer"]]],
            Keys[qc[Method -> "Stabilizer"]]
        }
    ],
    {True, "TooLarge", {{0, 0}}},
    TestID -> "MidMeasurement-TooLarge"
]


(* ---------- general case: the number of records ---------- *)

(* Measuring every wire of a stabilizer state gives the records b0 + span of   *)
(* the generators' X parts over GF(2): 2^rank records, and every difference    *)
(* lies in that span (adding the differences does not raise the rank). States  *)
(* of small X rank make the span a small part of GF(2)^n: H on wires 1 to 3,    *)
(* S on 2, then random CNOTs (which keep the X rank) give rank 3 at any n.     *)
stabMidAffineLawQ[ps_PauliStabilizer] := With[{
    records = Keys[QuantumCircuitOperator[stabMidM[{#}] & /@ Range[ps["Qubits"]]][ps, Method -> "Stabilizer"]],
    xs = Transpose[ps["StabilizerX"]]
},
    With[{rank = MatrixRank[xs, Modulus -> 2]},
        {rank, Length[records] == 2 ^ rank, MatrixRank[Join[Mod[# - First[records] & /@ records, 2], xs], Modulus -> 2] == rank}
    ]
]

stabMidLowRank[n_] := PauliStabilizer[n]["ApplyCircuit",
    Join[{"H" -> 1, "H" -> 2, "H" -> 3, "S" -> 2}, Table["CNOT" -> RandomSample[Range[n], 2], 3 n]]]

VerificationTest[
    BlockRandom[SeedRandom[61];
        {stabMidAffineLawQ[stabMidLowRank[12]], stabMidAffineLawQ[stabMidLowRank[30]], Rest[stabMidAffineLawQ[PauliStabilizer["Random"[7]]]]}
    ],
    {{3, True, True}, {3, True, True}, {True, True}},
    TestID -> "MidMeasurement-RecordsAreAnAffineSpan"
]

(* Mid-circuit measurements at n = 24, gates continuing between them: the k     *)
(* random measurements give 2^k records whose differences span k dimensions.    *)
VerificationTest[
    Map[
        Function[seed, BlockRandom[SeedRandom[seed];
            With[{records = Keys[QuantumCircuitOperator[Catenate @ Table[Append[Table[stabMeasRandomGate[24], 6], stabMidM[{RandomInteger[{1, 24}]}]], 14]][
                    stabMidLowRank[24], Method -> "Stabilizer"]]},
                Log2[Length[records]] == MatrixRank[Mod[# - First[records] & /@ records, 2], Modulus -> 2]
            ]
        ]],
        {1, 2, 3}
    ],
    {True, True, True},
    TestID -> "MidMeasurement-LargeMidCircuitRecordLaw"
]


(* ---------- input forms and instrument labels ---------- *)

(* a pure state given as a density matrix, and an unnormalized ket *)
VerificationTest[
    {
        stabMidSameQ[
            stabMidFromTableaux[QuantumCircuitOperator[{"CNOT" -> {1, 2}, stabMidM[{2}]}][QuantumState[{{1, 0, 1, 0}, {0, 0, 0, 0}, {1, 0, 1, 0}, {0, 0, 0, 0}} / 2], Method -> "Stabilizer"]],
            stabMidExpected[{{{0}, 1/2, {1, 0, 0, 0}}, {{1}, 1/2, {0, 0, 0, 1}}}]
        ],
        Tr[stabMeasTableauDensity[QuantumCircuitOperator[{stabMidM[{1}]}][QuantumState[{2, 0}], Method -> "Stabilizer"][{0}]]]
    },
    {True, 1},
    TestID -> "MidMeasurement-DensityMatrixAndUnnormalizedInput"
]

(* a named basis labels its outcomes +/-: bit 0 is the first outcome, + *)
VerificationTest[
    With[{qc = QuantumCircuitOperator[{stabMidM["Z", {1}]}]},
        {
            Keys[qc[QuantumState["0"], Method -> "Stabilizer"]], qc[QuantumState["0"]]["ProbabilitiesList"],
            Keys[qc[QuantumState["1"], Method -> "Stabilizer"]], qc[QuantumState["1"]]["ProbabilitiesList"]
        }
    ],
    {{{0}}, {1, 0}, {{1}}, {0, 1}},
    TestID -> "MidMeasurement-NamedBasisOutcomeOrder"
]

(* wires listed in an order that is its own inverse are read correctly; a 3-cycle *)
(* order is refused, built-in or listed by hand                                  *)
VerificationTest[
    {
        stabMidSameQ[
            stabMidFromTableaux[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 3}, stabMidM[{3, 2, 1}]}][QuantumState["010"], Method -> "Stabilizer"]],
            stabMidFromDense[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 3}, stabMidM[{3, 2, 1}]}][QuantumState["010"]]]
        ],
        stabMidFailureReason @ QuantumCircuitOperator[{stabMidM[{2, 3, 1}]}][Method -> "Stabilizer"],
        stabMidFailureReason @ QuantumCircuitOperator[{stabMidM[QuantumOperator[DiagonalMatrix[UnitVector[8, #]], {2, 3, 1}] & /@ Range[8], {2, 3, 1}]}][Method -> "Stabilizer"]
    },
    {True, "WireOrder", "WireOrder"},
    TestID -> "MidMeasurement-WireOrderRefused"
]

(* projectors listed out of computational order would relabel the outcomes *)
VerificationTest[
    stabMidFailureReason @ QuantumCircuitOperator[{stabMidM[{QuantumOperator[{{0, 0}, {0, 1}}], QuantumOperator[{{1, 0}, {0, 0}}]}, {1}]}][Method -> "Stabilizer"],
    "NotComputationalMeasurement",
    TestID -> "MidMeasurement-OutOfOrderProjectorsRefused"
]
