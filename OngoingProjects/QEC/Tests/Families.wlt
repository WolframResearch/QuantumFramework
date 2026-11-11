(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Families.wlt

   Named codes and parametric families, all reached through the QECCode object
   itself: QECCode["SteaneCode"], QECCode["Repetition", 5], QECCode["Hamming", 4].
   The scalable families of roadmap item (f) -- the surface and toric codes --
   are tested at the end of the file.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

(* Transitional: the rebuilt QEC core still lives under OngoingProjects/QEC/.
   Once it moves into Kernel/QEC/ and PacletInfo.wl lists its context, this Get
   disappears and the Needs above is enough.  The repo root is found from the
   loaded paclet, which RunTests.wls points at this checkout through
   PacletDirectoryLoad, so the tests always run against the source tree. *)
Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

(* Derived properties are memoised on the code's data.  A kernel that already
   held results from an earlier load would let those outrank the definitions just
   read, so the tests would silently check the old package. *)
QECClearCache[];

(* The row functions and helpers that the public QECPauli and QECCode["CSS", ...]
   stand on, in PackageScope since step 5 of the API redesign. *)
qecHammingMatrix = Symbol["Wolfram`QuantumFramework`QEC`PackageScope`classicalHammingMatrix"];


(* ============================================================================
   Named codes
   ============================================================================ *)

VerificationTest[QECCode["5QubitCode"]["Parameters"], {5, 1, 3}, TestID -> "QEC-Families-five-parameters"]

VerificationTest[
    QECCode["5QubitCode"]["Generators"],
    {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"},
    TestID -> "QEC-Families-five-generators"
]

VerificationTest[QECCode["SteaneCode"]["Parameters"], {7, 1, 3}, TestID -> "QEC-Families-steane-parameters"]

VerificationTest[QECCode["ShorCode"]["Parameters"], {9, 1, 3}, TestID -> "QEC-Families-shor-parameters"]

VerificationTest[QECCode["BitFlipCode"]["Generators"], {"ZZI", "IZZ"}, TestID -> "QEC-Families-bitflip-generators"]

VerificationTest[QECCode["PhaseFlipCode"]["Generators"], {"XXI", "IXX"}, TestID -> "QEC-Families-phaseflip-generators"]

VerificationTest[Quiet[QECCode["NoSuchCode"]], $Failed, TestID -> "QEC-Families-unknown-name"]

VerificationTest[
    Quiet[Check[QECCode["NoSuchCode"], "fired", QECCode::unknown]],
    "fired",
    TestID -> "QEC-Families-unknown-name-message"
]


(* ============================================================================
   Repetition codes
   ============================================================================ *)

VerificationTest[
    QECCode["Repetition", 3]["Generators"],
    QECCode["BitFlipCode"]["Generators"],
    TestID -> "QEC-Families-repetition-three-is-bitflip"
]

VerificationTest[
    QECCode["PhaseRepetition", 3]["Generators"],
    QECCode["PhaseFlipCode"]["Generators"],
    TestID -> "QEC-Families-phase-repetition-three-is-phaseflip"
]

VerificationTest[QECCode["Repetition", 5]["Parameters"], {5, 1, 1}, TestID -> "QEC-Families-repetition-five"]

VerificationTest[Quiet[QECCode["Repetition", 1]], $Failed, TestID -> "QEC-Families-repetition-too-small"]

(* The repetition code detects bit flips and is blind to phase flips, which is why
   its quantum distance is 1 whatever n is. *)
VerificationTest[
    Table[QECCode["Repetition", n]["Distance"], {n, 2, 6}],
    {1, 1, 1, 1, 1},
    TestID -> "QEC-Families-repetition-distance-is-one"
]

VerificationTest[
    Table[QECCode["Repetition", n]["Syndrome", StringJoin[ReplacePart[ConstantArray["I", n], 1 -> "X"]]], {n, 3, 5}],
    {{1, 0}, {1, 0, 0}, {1, 0, 0, 0}},
    TestID -> "QEC-Families-repetition-detects-bit-flips"
]


(* ============================================================================
   The [[n, n-2, 2]] family
   ============================================================================ *)

VerificationTest[QECCode["DistanceTwo", 4]["Parameters"], {4, 2, 2}, TestID -> "QEC-Families-distance-two-4"]

VerificationTest[QECCode["DistanceTwo", 6]["Parameters"], {6, 4, 2}, TestID -> "QEC-Families-distance-two-6"]

VerificationTest[Quiet[QECCode["DistanceTwo", 5]], $Failed, TestID -> "QEC-Families-distance-two-odd"]

VerificationTest[Quiet[QECCode["DistanceTwo", 2]], $Failed, TestID -> "QEC-Families-distance-two-too-small"]

(* It saturates the quantum Singleton bound, n - k >= 2(d - 1). *)
VerificationTest[
    AllTrue[{4, 6, 8},
        With[{p = QECCode["DistanceTwo", #]["Parameters"]},
            p[[1]] - p[[2]] === 2 (p[[3]] - 1)
        ] &
    ],
    True,
    TestID -> "QEC-Families-distance-two-saturates-singleton"
]


(* ============================================================================
   The Hamming CSS family
   ============================================================================ *)

VerificationTest[Dimensions[qecHammingMatrix[3]], {3, 7}, TestID -> "QEC-Families-hamming-matrix-shape"]

VerificationTest[
    qecHammingMatrix[3],
    {{0, 0, 0, 1, 1, 1, 1}, {0, 1, 1, 0, 0, 1, 1}, {1, 0, 1, 0, 1, 0, 1}},
    TestID -> "QEC-Families-hamming-matrix"
]

VerificationTest[
    Sort[QECCode["Hamming", 3]["Generators"]] === Sort[QECCode["SteaneCode"]["Generators"]],
    True,
    TestID -> "QEC-Families-hamming-three-is-steane"
]

VerificationTest[QECCode["Hamming", 4]["Parameters"], {15, 7, 3}, TestID -> "QEC-Families-hamming-four"]

VerificationTest[Quiet[QECCode["Hamming", 2]], $Failed, TestID -> "QEC-Families-hamming-too-small"]

(* [[2^r - 1, 2^r - 1 - 2r, 3]]. *)
VerificationTest[
    Table[Take[QECCode["Hamming", r]["Parameters"], 2], {r, 3, 5}],
    Table[{2^r - 1, 2^r - 1 - 2 r}, {r, 3, 5}],
    TestID -> "QEC-Families-hamming-parameters-formula"
]


(* ============================================================================
   The catalog
   ============================================================================ *)

VerificationTest[Head[QECCode["Catalog"]], Dataset, TestID -> "QEC-Families-catalog-head"]

VerificationTest[
    KeyExistsQ[Normal[QECCode["Catalog"]], "SteaneCode"],
    True,
    TestID -> "QEC-Families-catalog-contains-steane"
]

VerificationTest[
    AllTrue[Values[Normal[QECCode["Catalog"]]], Keys[#] === {"n", "k", "d"} &],
    True,
    TestID -> "QEC-Families-catalog-columns"
]

VerificationTest[
    Normal[QECCode["Catalog"]]["ShorCode"],
    <|"n" -> 9, "k" -> 1, "d" -> 3|>,
    TestID -> "QEC-Families-catalog-shor-row"
]


(* ============================================================================
   The scalable families of roadmap item (f): the rotated surface code
   QECCode["Surface", d] = [[d^2, 1, d]] and the toric code
   QECCode["Toric", L] = [[2 L^2, 2, L]] (Got26 ch. 17).

   The parameters come out of the layer's own machinery -- rank for k, the
   exact distance search for d -- so the layout is tested by what it implies,
   not by re-reading the construction.  The geometry is checked as locality:
   every check of the rotated code acts on one 2 x 2 block of the grid.  And
   the circuits are the ones the StimCrossCheck run of November 2026 sampled
   with Stim (surface d = 3 at circuit level, toric L = 3 at phenomenological
   level, two million shots: every detector within 2.5 sigma, PyMatching
   decoding both unchanged); the detector counts pinned here are the ones
   Stim reported.
   ============================================================================ *)

surface3 = QECCode["Surface", 3];
toric3 = QECCode["Toric", 3];

(* The (row, column) cells of the d x d grid a check touches. *)
checkCells[c_QECCode, d_Integer] := With[{n = c["Qubits"]},
    Function[row, QuotientRemainder[# - 1, d] & /@ Flatten[Position[BitOr[row[[1 ;; n]], row[[n + 1 ;; 2 n]]], 1]]] /@
        c["GeneratorVectors"]
]


(* ============================================================================
   The rotated surface code
   ============================================================================ *)

VerificationTest[
    Table[QECCode["Surface", d]["Parameters"], {d, 2, 4}],
    {{4, 1, 2}, {9, 1, 3}, {16, 1, 4}},
    TestID -> "QEC-Topological-surface-parameters"
]

(* Exact distance search on 25 qubits: the slowest test in the file, and the one
   that shows the layout is right beyond the smallest cases. *)
VerificationTest[
    QECCode["Surface", 5]["Parameters"],
    {25, 1, 5},
    TestID -> "QEC-Topological-surface-distance-five"
]

VerificationTest[
    surface3["Generators"],
    {"IXXIIIIII", "XXIXXIIII", "IIIIXXIXX", "IIIIIIXXI", "ZIIZIIIII", "IZZIZZIII", "IIIZZIZZI", "IIIIIZIIZ"},
    TestID -> "QEC-Topological-surface-three-generators"
]

(* (d - 1)^2 weight-four checks in the bulk, 2 (d - 1) weight-two on the boundary. *)
VerificationTest[
    Table[Sort[Counts[QECPauli[#]["Weight"] & /@ QECCode["Surface", d]["Generators"]]], {d, 3, 5}],
    {<|2 -> 4, 4 -> 4|>, <|2 -> 6, 4 -> 9|>, <|2 -> 8, 4 -> 16|>},
    TestID -> "QEC-Topological-surface-check-weights"
]

VerificationTest[
    Table[
        AllTrue[checkCells[QECCode["Surface", d], d], Max[#[[All, 1]]] - Min[#[[All, 1]]] <= 1 && Max[#[[All, 2]]] - Min[#[[All, 2]]] <= 1 &],
        {d, 3, 5}
    ],
    {True, True, True},
    TestID -> "QEC-Topological-surface-checks-are-local"
]

VerificationTest[
    {surface3["CSSQ"], QECPauli[First[surface3["LogicalX"]]]["Weight"], QECPauli[First[surface3["LogicalZ"]]]["Weight"]},
    {True, 3, 3},
    TestID -> "QEC-Topological-surface-is-CSS-with-weight-d-logicals"
]

(* The CNOT order: X checks row by row, Z checks column by column, so the pair a
   mid-check ancilla fault leaves is horizontal for X and vertical for Z --
   perpendicular to the logical of the same letter.  The Xbar consequence is
   checked in Matching.wlt; a code without one gets its supports in qubit order.  ("Schedule" is
   the time layering of circuitSchedule, a different thing.) *)
VerificationTest[
    {surface3["CheckOrder"], QECCode["SteaneCode"]["CheckOrder"][[1]]},
    {{{2, 3}, {1, 2, 4, 5}, {5, 6, 8, 9}, {7, 8}, {1, 4}, {2, 5, 3, 6}, {4, 7, 5, 8}, {6, 9}}, {4, 5, 6, 7}},
    TestID -> "QEC-Topological-surface-check-order-keeps-hooks-off-the-logicals"
]

VerificationTest[
    Cases[QECSyndromeCircuit[surface3]["Instructions"], {"CNOT", q_, 15} :> q],
    {2, 5, 3, 6},
    TestID -> "QEC-Topological-surface-circuit-follows-the-check-order"
]

(* The exact rate is a polynomial whose leading order is p^2, as a distance-3 code
   owes; and the exact route agrees with the sampled one. *)
VerificationTest[
    With[{rate = QECLogicalErrorRate[surface3, QECNoiseModel["BitFlip", qecP]]},
        {Coefficient[Expand[rate], qecP, 1], Coefficient[Expand[rate], qecP, 2]}
    ],
    {0, 18},
    TestID -> "QEC-Topological-surface-rate-starts-at-p-squared"
]

VerificationTest[
    BlockRandom[SeedRandom[7];
        With[{noise = QECNoiseModel["Depolarizing", 1/20]},
            Abs[N[QECLogicalErrorRate[surface3, noise]] - QECLogicalErrorRate[surface3, noise, 200000]]
        ]
    ] < 0.002,
    True,
    TestID -> "QEC-Topological-surface-exact-and-sampled-rates-agree"
]

VerificationTest[
    0.097 < N[QECPseudothreshold[surface3, QECNoiseModel["Depolarizing", qecP]]["Value"]] < 0.098,
    True,
    TestID -> "QEC-Topological-surface-pseudothreshold"
]

VerificationTest[
    QECDetectorModel[surface3, QECNoiseModel["Circuit", 1/500], 3]["Detectors"],
    32,
    TestID -> "QEC-Topological-surface-detector-count-as-Stim-reports"
]

VerificationTest[
    StringQ[QECStim[surface3, QECNoiseModel["Circuit", 1/500], 3]],
    True,
    TestID -> "QEC-Topological-surface-exports-to-Stim"
]


(* ============================================================================
   The toric code
   ============================================================================ *)

VerificationTest[
    Table[QECCode["Toric", l]["Parameters"], {l, 2, 4}],
    {{8, 2, 2}, {18, 2, 3}, {32, 2, 4}},
    TestID -> "QEC-Topological-toric-parameters"
]

(* Every check acts on the four edges of a face or a star, and every edge sits on two
   faces and two stars; one of each is left out as dependent, so an edge is covered by
   the remaining generators two to four times, and never fewer. *)
VerificationTest[
    With[{gens = toric3["GeneratorVectors"], n = 18},
        {Length[gens], Min[Total[BitOr[#[[1 ;; n]], #[[n + 1 ;; 2 n]]] & /@ gens]], Union[Total /@ (BitOr[#[[1 ;; n]], #[[n + 1 ;; 2 n]]] & /@ gens)]}
    ],
    {16, 2, {4}},
    TestID -> "QEC-Topological-toric-checks-and-coverage"
]

VerificationTest[
    {toric3["CSSQ"], Length[toric3["LogicalX"]], QECPauli[toric3["MinimumWeightLogical"]]["Weight"]},
    {True, 2, 3},
    TestID -> "QEC-Topological-toric-two-logical-qubits-of-weight-L"
]

VerificationTest[
    QECDetectorModel[toric3, QECNoiseModel["Depolarizing", 1/100, "MeasurementError" -> 1/100], 2]["Detectors"],
    48,
    TestID -> "QEC-Topological-toric-detector-count-as-Stim-reports"
]


(* ============================================================================
   The scalable families: refusals and the catalog
   ============================================================================ *)

VerificationTest[
    QECCode["Surface", 1],
    $Failed,
    {QECCode::surfacesize},
    TestID -> "QEC-Topological-surface-refuses-d-one"
]

VerificationTest[
    QECCode["Toric", 1],
    $Failed,
    {QECCode::toricsize},
    TestID -> "QEC-Topological-toric-refuses-L-one"
]

VerificationTest[
    With[{cat = Normal[QECCode["Catalog"]]},
        {cat["QECCode[\"Surface\", 3]"], cat["QECCode[\"Toric\", 3]"]}
    ],
    {<|"n" -> 9, "k" -> 1, "d" -> 3|>, <|"n" -> 18, "k" -> 2, "d" -> 3|>},
    TestID -> "QEC-Topological-families-are-in-the-catalog"
]
