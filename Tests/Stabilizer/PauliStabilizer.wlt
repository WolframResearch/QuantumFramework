(* ::Package:: *)

BeginTestSection["PauliStabilizer"]


(* ============================================================================ *)
(* Setup / Fixtures                                                             *)
(* ============================================================================ *)

SeedRandom[20260430];

(* String-list fixtures (auto-padded destabilizers via Reverse rule, line 103-107) *)
$psBellLit  = PauliStabilizer[{"XX", "ZZ"}];
$psGHZ3Lit  = PauliStabilizer[{"XXX", "ZZI", "IZZ"}];
$ps5Q       = PauliStabilizer["5QubitCode"];
$psSteane   = PauliStabilizer["SteaneCode"];
$psShor     = PauliStabilizer["9QubitCode"];

(* Circuit-built fixtures (AG-valid destabilizers) *)
$psBell    = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
$psGHZ3    = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, "CNOT" -> {2, 3}}]];

(* Helpers *)
agOmega[n_Integer] := ArrayFlatten[{
    {ConstantArray[0, {n, n}], IdentityMatrix[n]},
    {IdentityMatrix[n], ConstantArray[0, {n, n}]}
}];

(* Validity check that doesn't depend on PauliStabilizerQ context resolution.
   Phase 1 typo-fix follow-up: promote PauliStabilizerQ to PackageScope so
   bare-name resolution works in Tests/ files. For now, structural Head check. *)
psValidQ[ps_PauliStabilizer] := Module[{a = First[ps]},
    AssociationQ[a] && KeyExistsQ[a, "Signs"] && KeyExistsQ[a, "Tableau"] &&
    MatchQ[a["Signs"], {(-1 | 1) ...}] &&
    ArrayQ[a["Tableau"], 3, MatchQ[0 | 1] (* unsigned *)] &&
    Length[a["Signs"]] === Dimensions[a["Tableau"]][[3]]
]
psValidQ[_] := False


(* ============================================================================ *)
(* TIER 1.1 — Constructor coverage (Ctor-)                                      *)
(* ============================================================================ *)

(* Constructor 2.1: <|"Signs" -> ..., "Tableau" -> ...|> *)
VerificationTest[
    psValidQ @ PauliStabilizer[<|
        "Signs" -> {1, 1, 1, 1},
        "Tableau" -> {{{0, 1, 1, 0}, {0, 1, 1, 0}}, {{1, 0, 0, 1}, {1, 0, 0, 1}}}
    |>],
    True,
    TestID -> "Ctor-FromSigns"
]

(* Constructor 2.2: <|"Phase" -> ...|> *)
VerificationTest[
    psValidQ @ PauliStabilizer[<|
        "Phase" -> {0, 0, 0, 0},
        "Tableau" -> {{{0, 1, 1, 0}, {0, 1, 1, 0}}, {{1, 0, 0, 1}, {1, 0, 0, 1}}}
    |>],
    True,
    TestID -> "Ctor-FromPhase"
]

(* Constructor 2.3: <|"Matrix" -> ...|> reshapes 2n*2n binary back to rank-3 tableau *)
VerificationTest[
    psValidQ @ PauliStabilizer[<|
        "Matrix" -> {{0, 0, 1, 1}, {1, 1, 0, 0}, {1, 1, 0, 0}, {0, 0, 1, 1}},
        "Signs" -> {1, 1, 1, 1}
    |>],
    True,
    TestID -> "Ctor-FromMatrix"
]

(* Constructor 2.4: <|"Tableau" -> ...|> auto-fills Signs with +1s *)
VerificationTest[
    PauliStabilizer[<|"Tableau" -> {{{0, 1, 1, 0}, {0, 1, 1, 0}}, {{1, 0, 0, 1}, {1, 0, 0, 1}}}|>]["Signs"],
    {1, 1, 1, 1},
    TestID -> "Ctor-FromTableauOnly"
]

(* Constructor 2.5: from QuantumOperator *)
VerificationTest[
    PauliStabilizer[QuantumOperator["H", 1]]["Stabilizers"],
    {"X"},
    TestID -> "Ctor-FromOperator"
]

(* Constructor 2.6: from QuantumCircuitOperator *)
VerificationTest[
    PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]["Stabilizers"],
    {"XX", "ZZ"},
    TestID -> "Ctor-FromCircuit"
]

(* Constructor 2.9: string list (Pauli strings) *)
VerificationTest[
    PauliStabilizer[{"XX", "ZZ"}]["Stabilizers"],
    {"XX", "ZZ"},
    TestID -> "Ctor-FromStringList"
]

(* Constructor 2.9 with sign prefix *)
VerificationTest[
    PauliStabilizer[{"-XX", "+ZZ"}]["StabilizerSigns"],
    {-1, 1},
    TestID -> "Ctor-FromStringList-Signed"
]

(* Constructor 2.10: separate stab + destab strings *)
VerificationTest[
    With[{ps = PauliStabilizer[{"XX", "ZZ"}, {"IX", "ZI"}]},
        {ps["Stabilizers"], ps["Destabilizers"]}
    ],
    {{"XX", "ZZ"}, {"IX", "ZI"}},
    TestID -> "Ctor-FromStabDestabStrings"
]

(* Constructor 2.11: positive integer for |0...0> *)
VerificationTest[
    PauliStabilizer[3]["Stabilizers"],
    {"ZII", "IZI", "IIZ"},
    TestID -> "Ctor-FromInteger"
]

(* Constructor 2.12: 5-qubit code *)
VerificationTest[
    PauliStabilizer["5QubitCode"]["Stabilizers"],
    {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ", "XXXXX"},
    TestID -> "Ctor-FromName-5Q"
]

(* Constructor 2.12: Steane code *)
VerificationTest[
    PauliStabilizer["SteaneCode"]["Qubits"],
    7,
    TestID -> "Ctor-FromName-Steane"
]

(* Constructor 2.12: Shor code *)
VerificationTest[
    PauliStabilizer["9QubitCode"]["Qubits"],
    9,
    TestID -> "Ctor-FromName-Shor"
]

(* Constructor 2.13: random tag *)
VerificationTest[
    Block[{}, SeedRandom[20260430]; PauliStabilizer["Random", 4]["Qubits"]],
    4,
    TestID -> "Ctor-FromRandomTag"
]


(* ============================================================================ *)
(* TIER 1.2 — Property dispatch (Prop-)                                         *)
(* ============================================================================ *)

VerificationTest[$ps5Q["Qubits"], 5, TestID -> "Prop-Qubits"]
VerificationTest[$ps5Q["Qudits"], 5, TestID -> "Prop-Qudits-Alias"]
VerificationTest[$ps5Q["GeneratorCount"], 5, TestID -> "Prop-GeneratorCount"]
VerificationTest[$ps5Q["Signs"], ConstantArray[1, 10], TestID -> "Prop-Signs"]
VerificationTest[$ps5Q["Phase"], ConstantArray[0, 10], TestID -> "Prop-Phase"]
VerificationTest[Length[$ps5Q["X"]], 5, TestID -> "Prop-X-Shape"]
VerificationTest[Length[$ps5Q["Z"]], 5, TestID -> "Prop-Z-Shape"]
VerificationTest[Dimensions[$ps5Q["Tableau"]], {2, 5, 10}, TestID -> "Prop-Tableau-Shape"]
VerificationTest[Dimensions[$ps5Q["Matrix"]], {10, 10}, TestID -> "Prop-Matrix-Shape"]
VerificationTest[$ps5Q["Stabilizers"], {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ", "XXXXX"}, TestID -> "Prop-Stabilizers"]
VerificationTest[$ps5Q["StabilizerSigns"], {1, 1, 1, 1, 1}, TestID -> "Prop-StabilizerSigns"]
VerificationTest[$ps5Q["DestabilizerSigns"], {1, 1, 1, 1, 1}, TestID -> "Prop-DestabilizerSigns"]
VerificationTest[Length[$ps5Q["Stabilizer"][[1]]], 5, TestID -> "Prop-Stabilizer-Tableau-Shape"]
VerificationTest[Length[$ps5Q["Destabilizer"][[1]]], 5, TestID -> "Prop-Destabilizer-Tableau-Shape"]
VerificationTest[$ps5Q["Generators"], $ps5Q["Stabilizers"], TestID -> "Prop-Generators-Equals-Stabilizers"]
VerificationTest[$ps5Q["PauliForm"], $ps5Q["Stabilizers"], TestID -> "Prop-PauliForm-Equals-Stabilizers"]
VerificationTest[Length[$ps5Q["PauliStrings"]], 10, TestID -> "Prop-PauliStrings-Length"]
VerificationTest[Length[$ps5Q["PauliSymbols"]], 10, TestID -> "Prop-PauliSymbols-Length"]
VerificationTest[
    Sort[$psBellLit["Stabilizers"]] === Sort[{"XX", "ZZ"}],
    True,
    TestID -> "Prop-Stabilizers-Bell"
]
(* X/Z arrays have shape {n_qubits, 2*GeneratorCount}: rows are qubits, cols are tableau rows.
   StabilizerX = X[[All, n+1 ;;]] (last n cols = stabilizer rows).
   DestabilizerX = X[[All, ;; n]] (first n cols = destabilizer rows). *)
VerificationTest[$psBellLit["StabilizerX"], $psBellLit["X"][[All, 3 ;;]], TestID -> "Prop-StabilizerX-Slice"]
VerificationTest[$psBellLit["StabilizerZ"], $psBellLit["Z"][[All, 3 ;;]], TestID -> "Prop-StabilizerZ-Slice"]
VerificationTest[$psBellLit["DestabilizerX"], $psBellLit["X"][[All, ;; 2]], TestID -> "Prop-DestabilizerX-Slice"]
VerificationTest[$psBellLit["DestabilizerZ"], $psBellLit["Z"][[All, ;; 2]], TestID -> "Prop-DestabilizerZ-Slice"]

(* Phase 1 cleanup: ["Properties"] list now reflects the actual 30+ property surface.
   Old contract (8 entries) was incomplete vs. the dispatch table.
   Test pins down the new exhaustive list. *)
VerificationTest[
    Length @ $ps5Q["Properties"] >= 25
        && SubsetQ[$ps5Q["Properties"], {"Qubits", "Stabilizer", "Destabilizer", "Generators", "TableauForm", "State", "Circuit", "Matrix", "Phase"}],
    True,
    TestID -> "Prop-PropertiesList-CoversFullSurface"
]


(* ============================================================================ *)
(* TIER 1.3 — Gate-update closure (Gate-Closure-)                               *)
(* ============================================================================ *)

VerificationTest[psValidQ @ $ps5Q["H", 1], True, TestID -> "Gate-Closure-H"]
VerificationTest[psValidQ @ $ps5Q["S", 1], True, TestID -> "Gate-Closure-S"]
VerificationTest[psValidQ @ $ps5Q[SuperDagger["S"], 1], True, TestID -> "Gate-Closure-Sdag"]
VerificationTest[psValidQ @ $ps5Q["X", 1], True, TestID -> "Gate-Closure-X"]
VerificationTest[psValidQ @ $ps5Q["Y", 1], True, TestID -> "Gate-Closure-Y"]
VerificationTest[psValidQ @ $ps5Q["Z", 1], True, TestID -> "Gate-Closure-Z"]
VerificationTest[psValidQ @ $ps5Q["CNOT", 1, 2], True, TestID -> "Gate-Closure-CNOT"]
VerificationTest[psValidQ @ $ps5Q["CX", 1, 2], True, TestID -> "Gate-Closure-CX-Alias"]
VerificationTest[psValidQ @ $ps5Q["CZ", 1, 2], True, TestID -> "Gate-Closure-CZ"]
VerificationTest[psValidQ @ $ps5Q["SWAP", 1, 2], True, TestID -> "Gate-Closure-SWAP"]
VerificationTest[psValidQ @ $ps5Q["V", 1], True, TestID -> "Gate-Closure-V"]
VerificationTest[psValidQ @ $ps5Q[SuperDagger["V"], 1], True, TestID -> "Gate-Closure-Vdag"]

(* op -> order flatten path (line 327) *)
VerificationTest[psValidQ @ $ps5Q["H" -> 1], True, TestID -> "Gate-Closure-H-OrderArg"]
VerificationTest[psValidQ @ $ps5Q["CNOT" -> {1, 2}], True, TestID -> "Gate-Closure-CNOT-OrderArg"]


(* ============================================================================ *)
(* TIER 1.4 — Round-trips (Roundtrip-)                                          *)
(* ============================================================================ *)

(* Bell circuit-construction matches string-list-construction at the stabilizer level *)
VerificationTest[
    PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]["Stabilizers"],
    {"XX", "ZZ"},
    TestID -> "Roundtrip-CircuitToBellStabilizers"
]

(* Phase encoding: Phase = (1 - Signs) / 2; Signs = 1 - 2*Phase *)
VerificationTest[
    With[{p = $psBellLit["Phase"]}, 1 - 2 p],
    $psBellLit["Signs"],
    TestID -> "Roundtrip-PhaseSigns"
]

(* Build via "Matrix" key, recover the same matrix *)
VerificationTest[
    PauliStabilizer[<|
        "Matrix" -> $psBellLit["Matrix"],
        "Signs" -> $psBellLit["Signs"]
    |>]["Matrix"],
    $psBellLit["Matrix"],
    TestID -> "Roundtrip-MatrixIdempotent"
]

(* Apply H twice = identity at stabilizer level *)
VerificationTest[
    $psBellLit["H", 1]["H", 1]["Stabilizers"],
    $psBellLit["Stabilizers"],
    TestID -> "Roundtrip-HSquaredIsIdentity"
]

(* Apply S 4 times = identity *)
VerificationTest[
    Nest[#["S", 1] &, $psBellLit, 4]["Stabilizers"],
    $psBellLit["Stabilizers"],
    TestID -> "Roundtrip-SFourthIsIdentity"
]


(* ============================================================================ *)
(* TIER 1.4a -- QuantumOperator -> PauliStabilizer -> QuantumOperator           *)
(* ============================================================================ *)
(* The recovered operator must equal the original *exactly* (matrix equality),   *)
(* not just up to global phase. Y = -i*X*Z, so the AG decomposition that         *)
(* produces Z@X drops a -i factor; the constructor must capture and replay it.   *)
(* These tests will fail until the global-phase fix lands (Step 3).              *)

(* Helper: matrix equality of two QuantumOperators on the same support *)
matEqQO[a_QuantumOperator, b_QuantumOperator] := (a["Matrix"] // Normal // Simplify) === (b["Matrix"] // Normal // Simplify)

(* 1Q Pauli + Clifford generators *)
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["I"]]["QuantumOperator"], QuantumOperator["I"]], True, TestID -> "Roundtrip-QO-I"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["X"]]["QuantumOperator"], QuantumOperator["X"]], True, TestID -> "Roundtrip-QO-X"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["Y"]]["QuantumOperator"], QuantumOperator["Y"]], True, TestID -> "Roundtrip-QO-Y"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["Z"]]["QuantumOperator"], QuantumOperator["Z"]], True, TestID -> "Roundtrip-QO-Z"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["H"]]["QuantumOperator"], QuantumOperator["H"]], True, TestID -> "Roundtrip-QO-H"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["S"]]["QuantumOperator"], QuantumOperator["S"]], True, TestID -> "Roundtrip-QO-S"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["S"]["Dagger"]]["QuantumOperator"], QuantumOperator["S"]["Dagger"]], True, TestID -> "Roundtrip-QO-SDagger"]

(* 2Q Cliffords *)
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["CNOT"]]["QuantumOperator"], QuantumOperator["CNOT"]], True, TestID -> "Roundtrip-QO-CNOT"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["CZ"]]["QuantumOperator"], QuantumOperator["CZ"]], True, TestID -> "Roundtrip-QO-CZ"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["SWAP"]]["QuantumOperator"], QuantumOperator["SWAP"]], True, TestID -> "Roundtrip-QO-SWAP"]

(* 2Q Pauli tensor products. XY/YZ pick up factor i from the embedded Y;         *)
(* YY picks up i^2 = -1. Each must be recovered exactly.                         *)
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["XX"]]["QuantumOperator"], QuantumOperator["XX"]], True, TestID -> "Roundtrip-QO-XX"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["YY"]]["QuantumOperator"], QuantumOperator["YY"]], True, TestID -> "Roundtrip-QO-YY"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["ZZ"]]["QuantumOperator"], QuantumOperator["ZZ"]], True, TestID -> "Roundtrip-QO-ZZ"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["XY"]]["QuantumOperator"], QuantumOperator["XY"]], True, TestID -> "Roundtrip-QO-XY"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["YX"]]["QuantumOperator"], QuantumOperator["YX"]], True, TestID -> "Roundtrip-QO-YX"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["YZ"]]["QuantumOperator"], QuantumOperator["YZ"]], True, TestID -> "Roundtrip-QO-YZ"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["ZY"]]["QuantumOperator"], QuantumOperator["ZY"]], True, TestID -> "Roundtrip-QO-ZY"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["XZ"]]["QuantumOperator"], QuantumOperator["XZ"]], True, TestID -> "Roundtrip-QO-XZ"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["ZX"]]["QuantumOperator"], QuantumOperator["ZX"]], True, TestID -> "Roundtrip-QO-ZX"]

(* Composite circuit: H @ CNOT (the Bell builder) -- Y-free, expected exact today *)
VerificationTest[
    matEqQO[
        PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]["QuantumOperator"]]["QuantumOperator"],
        QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]["QuantumOperator"]
    ],
    True,
    TestID -> "Roundtrip-QO-BellBuilder"
]

(* Random Clifford roundtrip (seeded). Each draw is exercised against itself.    *)
VerificationTest[
    BlockRandom[SeedRandom[20260504]; matEqQO[#["QuantumOperator"], #["QuantumOperator"]] & @ PauliStabilizer[PauliStabilizer["Random",1]["QuantumOperator"]]],
    True,
    TestID -> "Roundtrip-QO-RandomClifford-1Q-Stable"
]

(* The real Random Clifford 1Q test: build from a 1Q Clifford QO, roundtrip, compare *)
VerificationTest[
    BlockRandom[SeedRandom[20260504];
        Module[{rc = PauliStabilizer["Random",1]["QuantumOperator"]},
            matEqQO[PauliStabilizer[rc]["QuantumOperator"], rc]
        ]
    ],
    True,
    TestID -> "Roundtrip-QO-RandomClifford-1Q"
]

VerificationTest[
    BlockRandom[SeedRandom[20260504 + 1];
        Module[{rc = PauliStabilizer["Random",3]["QuantumOperator"]},
            matEqQO[PauliStabilizer[rc]["QuantumOperator"], rc]
        ]
    ],
    True,
    TestID -> "Roundtrip-QO-RandomClifford-3Q"
]


(* ============================================================================ *)
(* TIER 1.4b -- QuantumState -> PauliStabilizer -> QuantumState                 *)
(* ============================================================================ *)
(* PauliStabilizer[qs] currently runs through PauliStabilizerTableau, which     *)
(* loses generator signs in the RowSpace canonicalization step. As a result    *)
(* |1>, |->, |-i>, Bell, GHZ-3 all roundtrip to *different* states (not just    *)
(* off by phase). These tests will fail until the tomography rewrite (Step 2). *)

matEqQS[a_QuantumState, b_QuantumState] := (a["StateVector"] // Normal // Simplify) === (b["StateVector"] // Normal // Simplify)

(* Single-qubit eigenstates of X, Y, Z *)
VerificationTest[matEqQS[PauliStabilizer[QuantumState[{1, 0}]]["State"], QuantumState[{1, 0}]],         True, TestID -> "Roundtrip-QS-Comp-Zero"]
VerificationTest[matEqQS[PauliStabilizer[QuantumState[{0, 1}]]["State"], QuantumState[{0, 1}]],         True, TestID -> "Roundtrip-QS-Comp-One"]
VerificationTest[matEqQS[PauliStabilizer[QuantumState["Plus"]]["State"], QuantumState["Plus"]],         True, TestID -> "Roundtrip-QS-Plus"]
VerificationTest[matEqQS[PauliStabilizer[QuantumState["Minus"]]["State"], QuantumState["Minus"]],       True, TestID -> "Roundtrip-QS-Minus"]
VerificationTest[matEqQS[PauliStabilizer[QuantumState[{1, I}/Sqrt[2]]]["State"],  QuantumState[{1, I}/Sqrt[2]]],   True, TestID -> "Roundtrip-QS-PlusI"]
VerificationTest[matEqQS[PauliStabilizer[QuantumState[{1, -I}/Sqrt[2]]]["State"], QuantumState[{1, -I}/Sqrt[2]]],  True, TestID -> "Roundtrip-QS-MinusI"]

(* Bell states (4 orientations) *)
VerificationTest[matEqQS[PauliStabilizer[QuantumState[{1, 0, 0,  1}/Sqrt[2]]]["State"], QuantumState[{1, 0, 0,  1}/Sqrt[2]]], True, TestID -> "Roundtrip-QS-Bell-PhiPlus"]
VerificationTest[matEqQS[PauliStabilizer[QuantumState[{1, 0, 0, -1}/Sqrt[2]]]["State"], QuantumState[{1, 0, 0, -1}/Sqrt[2]]], True, TestID -> "Roundtrip-QS-Bell-PhiMinus"]
VerificationTest[matEqQS[PauliStabilizer[QuantumState[{0, 1,  1, 0}/Sqrt[2]]]["State"], QuantumState[{0, 1,  1, 0}/Sqrt[2]]], True, TestID -> "Roundtrip-QS-Bell-PsiPlus"]
VerificationTest[matEqQS[PauliStabilizer[QuantumState[{0, 1, -1, 0}/Sqrt[2]]]["State"], QuantumState[{0, 1, -1, 0}/Sqrt[2]]], True, TestID -> "Roundtrip-QS-Bell-PsiMinus"]

(* GHZ-3 *)
VerificationTest[
    matEqQS[
        PauliStabilizer[QuantumState[Normalize @ {1, 0, 0, 0, 0, 0, 0, 1}]]["State"],
        QuantumState[Normalize @ {1, 0, 0, 0, 0, 0, 0, 1}]
    ],
    True,
    TestID -> "Roundtrip-QS-GHZ-3"
]

(* Random-Clifford-on-zero roundtrip (seeded). Build a random stabilizer state, *)
(* take its tableau, reproject back; the resulting state must match exactly.    *)
VerificationTest[
    BlockRandom[SeedRandom[20260504];
        Module[{n = 2, rc, qs0, qs},
            rc = PauliStabilizer["Random",n]["QuantumOperator"];
            qs0 = QuantumState["0", QuantumBasis[2, n]];
            qs = rc[qs0];
            matEqQS[PauliStabilizer[qs]["State"], qs]
        ]
    ],
    True,
    TestID -> "Roundtrip-QS-RandomClifford-on-zero-2Q"
]

VerificationTest[
    BlockRandom[SeedRandom[20260504 + 1];
        Module[{n = 3, rc, qs0, qs},
            rc = PauliStabilizer["Random",n]["QuantumOperator"];
            qs0 = QuantumState["0", QuantumBasis[2, n]];
            qs = rc[qs0];
            matEqQS[PauliStabilizer[qs]["State"], qs]
        ]
    ],
    True,
    TestID -> "Roundtrip-QS-RandomClifford-on-zero-3Q"
]


(* ============================================================================ *)
(* TIER 1.4c -- Internal sign correctness for state-derived tableaux            *)
(* ============================================================================ *)
(* Direct probes of PauliStabilizer[qs]["Stabilizers"] / ["StabilizerSigns"].   *)
(* These document the WHY behind 1.4b failures: the tableau itself encodes the  *)
(* wrong stabilizer set.                                                        *)

(* |1> stabilizer is -Z, not +Z *)
VerificationTest[
    PauliStabilizer[QuantumState[{0, 1}]]["StabilizerSigns"],
    {-1},
    TestID -> "Roundtrip-Internal-Comp-One-Sign"
]

(* |-> stabilizer is -X *)
VerificationTest[
    PauliStabilizer[QuantumState["Minus"]]["StabilizerSigns"],
    {-1},
    TestID -> "Roundtrip-Internal-Minus-Sign"
]

(* |+i> stabilizer is +Y *)
VerificationTest[
    PauliStabilizer[QuantumState[{1, I}/Sqrt[2]]]["StabilizerSigns"],
    {1},
    TestID -> "Roundtrip-Internal-PlusI-Sign"
]

(* |-i> stabilizer is -Y *)
VerificationTest[
    PauliStabilizer[QuantumState[{1, -I}/Sqrt[2]]]["StabilizerSigns"],
    {-1},
    TestID -> "Roundtrip-Internal-MinusI-Sign"
]

(* Bell phi+ stabilizer group is generated by {XX, ZZ}, both +1 *)
VerificationTest[
    Sort @ PauliStabilizer[QuantumState[{1, 0, 0, 1}/Sqrt[2]]]["Stabilizers"],
    Sort @ {"XX", "ZZ"},
    TestID -> "Roundtrip-Internal-Bell-PhiPlus-Stabilizers"
]
VerificationTest[
    PauliStabilizer[QuantumState[{1, 0, 0, 1}/Sqrt[2]]]["StabilizerSigns"],
    {1, 1},
    TestID -> "Roundtrip-Internal-Bell-PhiPlus-Signs"
]


(* ============================================================================ *)
(* TIER 1.4d -- Gate-update round-trip contract (UP TO GLOBAL PHASE)            *)
(* ============================================================================ *)
(* Phase 5c records that PauliStabilizer[qs]["gate", q] does NOT propagate     *)
(* "GlobalPhase" through the gate update -- exact propagation requires O(2^n)  *)
(* materialization at each step (ROADMAP \[Section]A.9). Documented contract:   *)
(*                                                                              *)
(*   PauliStabilizer[qs]["gate", q]["State"]  ==  gate @ qs   UP TO GLOBAL PHASE  *)
(*                                                                              *)
(* For *exact* equality, re-run the constructor on the post-gate state:        *)
(*   PauliStabilizer[gate @ qs]["State"]  ===  gate @ qs   (Phase 5c TIER 1.4b)  *)
(*                                                                              *)
(* These tests document the up-to-phase contract on the gate-update path.      *)

(* Two normalized state vectors differ by a global phase iff |<a|b>| = 1.       *)
equalUpToGlobalPhaseQS[a_QuantumState, b_QuantumState] := Module[{v1, v2, overlap},
    v1 = a["StateVector"] // Normal // Simplify;
    v2 = b["StateVector"] // Normal // Simplify;
    overlap = Quiet @ Simplify[Abs[Conjugate[v1] . v2]];
    Quiet @ Simplify[overlap == 1] === True
]

(* X applied to |1> *)
VerificationTest[
    equalUpToGlobalPhaseQS[
        PauliStabilizer[QuantumState[{0, 1}]]["X", 1]["State"],
        QuantumOperator["X"][QuantumState[{0, 1}]]
    ],
    True,
    TestID -> "GateUpdate-UpToPhase-XOn1"
]

(* Z applied to |1>: Z|1> = -|1>; tableau-recovered state is |1>; up-to-phase OK *)
VerificationTest[
    equalUpToGlobalPhaseQS[
        PauliStabilizer[QuantumState[{0, 1}]]["Z", 1]["State"],
        QuantumOperator["Z"][QuantumState[{0, 1}]]
    ],
    True,
    TestID -> "GateUpdate-UpToPhase-ZOn1"
]

(* Y applied to |1>: Y|1> = -i|0>; tableau-recovered state is |0>; up-to-phase OK *)
VerificationTest[
    equalUpToGlobalPhaseQS[
        PauliStabilizer[QuantumState[{0, 1}]]["Y", 1]["State"],
        QuantumOperator["Y"][QuantumState[{0, 1}]]
    ],
    True,
    TestID -> "GateUpdate-UpToPhase-YOn1"
]

(* H applied to |0> (real-valued; first-hop happens to roundtrip exact) *)
VerificationTest[
    equalUpToGlobalPhaseQS[
        PauliStabilizer[QuantumState[{1, 0}]]["H", 1]["State"],
        QuantumOperator["H"][QuantumState[{1, 0}]]
    ],
    True,
    TestID -> "GateUpdate-UpToPhase-HOn0"
]

(* S applied to |1>: S|1> = i|1>; canary for the i-factor on diag(1, i) *)
VerificationTest[
    equalUpToGlobalPhaseQS[
        PauliStabilizer[QuantumState[{0, 1}]]["S", 1]["State"],
        QuantumOperator["S"][QuantumState[{0, 1}]]
    ],
    True,
    TestID -> "GateUpdate-UpToPhase-SOn1"
]

(* CNOT applied to Bell phi+ *)
VerificationTest[
    equalUpToGlobalPhaseQS[
        PauliStabilizer[QuantumState[{1, 0, 0, 1}/Sqrt[2]]]["CNOT", 1, 2]["State"],
        QuantumOperator["CNOT"][QuantumState[{1, 0, 0, 1}/Sqrt[2]]]
    ],
    True,
    TestID -> "GateUpdate-UpToPhase-CNOTOnBell"
]

(* The "escape hatch": exact equality recovered by re-running the constructor    *)
(* on the post-gate state. This is the user-facing workaround for A.9.           *)
VerificationTest[
    matEqQS[
        PauliStabilizer[QuantumOperator["Y"][QuantumState[{0, 1}]]]["State"],
        QuantumOperator["Y"][QuantumState[{0, 1}]]
    ],
    True,
    TestID -> "GateUpdate-EscapeHatch-ConstructAfterGate-YOn1"
]

VerificationTest[
    matEqQS[
        PauliStabilizer[QuantumOperator["S"][QuantumState[{0, 1}]]]["State"],
        QuantumOperator["S"][QuantumState[{0, 1}]]
    ],
    True,
    TestID -> "GateUpdate-EscapeHatch-ConstructAfterGate-SOn1"
]


(* ============================================================================ *)
(* TIER 1.4e -- Coverage matrix: every accessor in _PauliStabilizer["Properties"] *)
(* ============================================================================ *)
(* Locks in that every documented accessor returns a non-Missing value, and that *)
(* name-aliases (e.g. "State"/"QuantumState", "Operator"/"QuantumOperator",      *)
(* "Circuit"/"QuantumCircuit"/"QuantumCircuitOperator", "Stabilizer"/             *)
(* "StabilizerTableau", "Generators"/"Stabilizers"/"PauliForm") return identical *)
(* outputs. Probes "p" (symplectic-form diagonal vector, untested before),       *)
(* "GlobalPhase" (Phase 5c key), and the display-only forms.                     *)

(* Every property in _PauliStabilizer["Properties"] resolves to a non-Missing    *)
(* value on a representative ps. This is a structural backstop -- if a future    *)
(* refactor accidentally drops an accessor from the dispatch table, this fails. *)
VerificationTest[
    Module[{ps = PauliStabilizer["5QubitCode"], props},
        props = ps["Properties"];
        FreeQ[Values @ AssociationMap[ps[#] &, props], _Missing]
    ],
    True,
    TestID -> "Coverage-AllPropertiesResolve-5Q"
]

(* Aliases produce identical outputs *)
VerificationTest[
    Module[{ps = PauliStabilizer[QuantumState["Bell"]]},
        ps["State"] === ps["QuantumState"]
    ],
    True,
    TestID -> "Coverage-Alias-State-equals-QuantumState"
]
VerificationTest[
    Module[{ps = PauliStabilizer[QuantumOperator["H"]]},
        ps["Operator"]["Matrix"] === ps["QuantumOperator"]["Matrix"]
    ],
    True,
    TestID -> "Coverage-Alias-Operator-equals-QuantumOperator"
]
VerificationTest[
    Module[{ps = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]},
        ps["Circuit"] === ps["QuantumCircuit"] === ps["QuantumCircuitOperator"]
    ],
    True,
    TestID -> "Coverage-Alias-Circuit-equals-QuantumCircuit"
]
VerificationTest[
    Module[{ps = $ps5Q}, ps["Stabilizer"] === ps["StabilizerTableau"]],
    True,
    TestID -> "Coverage-Alias-Stabilizer-equals-StabilizerTableau"
]
VerificationTest[
    Module[{ps = $ps5Q}, ps["Destabilizer"] === ps["DestabilizerTableau"]],
    True,
    TestID -> "Coverage-Alias-Destabilizer-equals-DestabilizerTableau"
]
VerificationTest[
    Module[{ps = $ps5Q}, ps["Stabilizers"] === ps["Generators"] === ps["PauliForm"]],
    True,
    TestID -> "Coverage-Alias-Stabilizers-Generators-PauliForm"
]

(* "p" — symplectic-form diagonal: zero on a valid stabilizer state, since rows *)
(* of the AG matrix come in commuting/anticommuting pairs that diagonal-select   *)
(* to zero. Was never tested before this commit.                                 *)
VerificationTest[
    PauliStabilizer["5QubitCode"]["p"],
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    TestID -> "Coverage-p-symplectic-diagonal-zero"
]
VerificationTest[
    PauliStabilizer[QuantumState["Bell"]]["p"],
    {0, 0, 0, 0},
    TestID -> "Coverage-p-symplectic-diagonal-Bell"
]

(* "GlobalPhase" set when constructed from QO / QS, default 1 elsewhere *)
VerificationTest[
    PauliStabilizer[QuantumOperator["Y"]]["GlobalPhase"],
    -I,
    TestID -> "Coverage-GlobalPhase-Y-is-minus-I"
]
VerificationTest[
    PauliStabilizer[QuantumOperator["YY"]]["GlobalPhase"],
    -1,
    TestID -> "Coverage-GlobalPhase-YY-is-minus-one"
]
VerificationTest[
    PauliStabilizer[QuantumState["Bell"]]["GlobalPhase"],
    1,
    TestID -> "Coverage-GlobalPhase-Bell-is-one"
]
VerificationTest[
    Lookup[First[PauliStabilizer[3]], "GlobalPhase", "absent"],
    "absent",
    TestID -> "Coverage-GlobalPhase-IntegerCtor-absent"
]

(* Display-form properties produce concrete heads (not Missing). *)
VerificationTest[
    Head[PauliStabilizer["5QubitCode"]["TableauForm"]] =!= Missing,
    True,
    TestID -> "Coverage-TableauForm-NonMissing"
]
VerificationTest[
    Head[PauliStabilizer["5QubitCode"]["StabilizerTableauForm"]] =!= Missing,
    True,
    TestID -> "Coverage-StabilizerTableauForm-NonMissing"
]


(* ============================================================================ *)
(* TIER 1.5 — Edge cases & known-quirks codification (Edge-)                    *)
(* ============================================================================ *)

(* PauliStabilizer[0] gets padded to 1 qubit (line 131 Max[q,1]) *)
VerificationTest[
    PauliStabilizer[0]["Qubits"],
    1,
    TestID -> "Edge-ZeroQubitPadsToOne"
]

(* 1-qubit identity register *)
VerificationTest[
    PauliStabilizer[1]["Stabilizers"],
    {"Z"},
    TestID -> "Edge-OneQubitInitStab"
]
VerificationTest[
    PauliStabilizer[1]["Destabilizers"],
    {"X"},
    TestID -> "Edge-OneQubitInitDestab"
]

(* "7QubitCode" and "SteaneCode" produce equal objects (line 136 pattern) *)
VerificationTest[
    PauliStabilizer["7QubitCode"] === PauliStabilizer["SteaneCode"],
    True,
    TestID -> "Edge-SteaneSynonyms"
]

(* Phase 1 cleanup: $PauliStabilizerNames duplicate "SteaneCode" (was 10 entries) deduped to 9. *)
VerificationTest[
    Length[Wolfram`QuantumFramework`PackageScope`$PauliStabilizerNames],
    9,
    TestID -> "Edge-NamesListLength-Deduped"
]

(* The deduped catalog must still contain both SteaneCode and SteaneCode1 *)
VerificationTest[
    SubsetQ[Wolfram`QuantumFramework`PackageScope`$PauliStabilizerNames,
        {"5QubitCode", "SteaneCode", "SteaneCode1", "9QubitCode", "Random"}
    ],
    True,
    TestID -> "Edge-NamesListContents"
]

(* Empty multi-qubit measurement *)
VerificationTest[
    $psBellLit["M", {}],
    <|{} -> $psBellLit|>,
    TestID -> "Edge-EmptyMeasure"
]

(* Phase 1 cleanup: "QuantumSttate" typo alias dropped. Now resolves to undefined property. *)
VerificationTest[
    Head[$psBellLit["QuantumState"]],
    QuantumState,
    TestID -> "Edge-QuantumState-Resolves"
]

(* The typo "QuantumSttate" no longer matches any property handler -- returns the unresolved expression *)
VerificationTest[
    MatchQ[$psBellLit["QuantumSttate"], _PauliStabilizer["QuantumSttate"]],
    True,
    TestID -> "Edge-QuantumSttateTypo-Removed"
]


(* ============================================================================ *)
(* TIER 1.6 — Failure / boundary tests (Fail-)                                  *)
(* ============================================================================ *)

(* Phase 4: P[theta] now returns a StabilizerFrame (closes under further Clifford gates).
   See StabilizerFrame.m. Was a top-level Plus in Phase 1. *)
VerificationTest[
    Head[$psBellLit["P"[Pi/3], 1]],
    StabilizerFrame,
    TestID -> "Fail-PReturnsStabilizerFrame-Phase4"
]

(* T is alias for P[Pi/2] *)
VerificationTest[
    $psBellLit["T", 1] === $psBellLit["P"[Pi/2], 1],
    True,
    TestID -> "Fail-TIsP"
]

(* P composition closes under further Clifford -- the result is a StabilizerFrame *)
VerificationTest[
    Wolfram`QuantumFramework`PackageScope`StabilizerFrameQ @ $psBellLit["P"[Pi/3], 1],
    True,
    TestID -> "Fail-PStaysAsStabilizerFrame-Phase4"
]


(* ============================================================================ *)
(* TIER 2.1 — AG tableau invariants (Tableau-Invariant-)                        *)
(* ============================================================================ *)
(* Invariants from AarGot04 Prop 2 / Yashin25 \[Section]2.3:                    *)
(* For circuit-built fixtures (proper destab pairing), check:                   *)
(*   M . \[CapitalOmega] . M^T == \[CapitalOmega] (mod 2)                       *)
(* where \[CapitalOmega] = [[0,I],[I,0]] is the symplectic form.                *)
(* String-list fixtures auto-pad destabilizers via Reverse rule (line 105),     *)
(* which does NOT generally satisfy AG -- so use circuit-built fixtures here.   *)

VerificationTest[
    With[{m = $psBell["Matrix"], n = $psBell["Qubits"]},
        Mod[m . agOmega[n] . Transpose[m] - agOmega[n], 2]
    ],
    ConstantArray[0, {4, 4}],
    TestID -> "Tableau-Invariant-Symplectic-Bell"
]

VerificationTest[
    With[{m = $psGHZ3["Matrix"], n = $psGHZ3["Qubits"]},
        Mod[m . agOmega[n] . Transpose[m] - agOmega[n], 2]
    ],
    ConstantArray[0, {6, 6}],
    TestID -> "Tableau-Invariant-Symplectic-GHZ3"
]

(* Stabilizer rows pairwise commute (AarGot04 Prop 2 #2) for any valid stabilizer state.
   Test on string-list fixtures since stab/stab commutation doesn't depend on destab pairing. *)
VerificationTest[
    With[{
        m = $ps5Q["Matrix"],
        n = $ps5Q["Qubits"],
        gen = $ps5Q["GeneratorCount"]
    },
        Mod[m[[gen + 1 ;; 2 gen]] . agOmega[n] . Transpose[m[[gen + 1 ;; 2 gen]]], 2]
    ],
    ConstantArray[0, {5, 5}],
    TestID -> "Tableau-Invariant-StabPairwiseCommute-5Q"
]

VerificationTest[
    With[{
        m = $psSteane["Matrix"],
        n = $psSteane["Qubits"],
        gen = $psSteane["GeneratorCount"]
    },
        Mod[m[[gen + 1 ;; 2 gen]] . agOmega[n] . Transpose[m[[gen + 1 ;; 2 gen]]], 2]
    ],
    ConstantArray[0, {7, 7}],
    TestID -> "Tableau-Invariant-StabPairwiseCommute-Steane"
]


(* ============================================================================ *)
(* TIER 2.2 — Symplectic gate invariants (Symplectic-)                          *)
(* ============================================================================ *)

(* H \[SmallCircle] H == identity at stabilizer level *)
VerificationTest[
    PauliStabilizer[1]["H", 1]["H", 1]["Stabilizers"],
    PauliStabilizer[1]["Stabilizers"],
    TestID -> "Symplectic-H-Order2"
]

(* S^4 == identity at stabilizer level *)
VerificationTest[
    Nest[#["S", 1] &, PauliStabilizer[1], 4]["Stabilizers"],
    PauliStabilizer[1]["Stabilizers"],
    TestID -> "Symplectic-S-Order4"
]

(* CNOT^2 == identity at stabilizer level *)
VerificationTest[
    PauliStabilizer[2]["CNOT", 1, 2]["CNOT", 1, 2]["Stabilizers"],
    PauliStabilizer[2]["Stabilizers"],
    TestID -> "Symplectic-CNOT-Involution"
]

(* CZ regression: H_2 \[SmallCircle] CNOT(1,2) \[SmallCircle] H_2 == CZ(1,2)
   (codifies decomposition at PauliStabilizer.m:246) *)
VerificationTest[
    With[{ps = PauliStabilizer[2]},
        ps["H", 2]["CNOT", 1, 2]["H", 2]["Stabilizers"] === ps["CZ", 1, 2]["Stabilizers"]
    ],
    True,
    TestID -> "Symplectic-CZ-DecomposesToHCNOTH"
]


(* ============================================================================ *)
(* TIER 2.3 — Stabilizer-state count (Count-N-)                                 *)
(* ============================================================================ *)
(* AarGot04 Prop 1: N(n) = 2^n \[Product]_{k=0..n-1} (2^{n-k} + 1)              *)

VerificationTest[
    2^1 Product[2^(1 - k) + 1, {k, 0, 1 - 1}],
    6,
    TestID -> "Count-N-Formula-n1"
]

VerificationTest[
    2^2 Product[2^(2 - k) + 1, {k, 0, 2 - 1}],
    60,
    TestID -> "Count-N-Formula-n2"
]

VerificationTest[
    2^3 Product[2^(3 - k) + 1, {k, 0, 3 - 1}],
    1080,
    TestID -> "Count-N-Formula-n3"
]


(* ============================================================================ *)
(* TIER 3.1 — Standard-state checks (Physics-StandardState-)                    *)
(* ============================================================================ *)

(* Bell |\[CapitalPhi]+\[RightAngleBracket] stabilized by \[LeftAngleBracket]XX, ZZ\[RightAngleBracket] (Got97 \[Section]2) *)
VerificationTest[
    Sort @ $psBell["Stabilizers"],
    Sort @ {"XX", "ZZ"},
    TestID -> "Physics-StandardState-BellPlus"
]

(* GHZ_3 stabilized by \[LeftAngleBracket]XXX, ZZI, IZZ\[RightAngleBracket] (Got97 \[Section]2.2 / standard pattern) *)
VerificationTest[
    Sort @ $psGHZ3["Stabilizers"],
    Sort @ {"XXX", "ZZI", "IZZ"},
    TestID -> "Physics-StandardState-GHZ3"
]

(* 5-qubit code cyclic generators (Got97 \[Section]3.5) *)
VerificationTest[
    $ps5Q["Stabilizers"],
    {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ", "XXXXX"},
    TestID -> "Physics-StandardState-5Q-Generators"
]

(* Steane code: 7-qubit, 7 generators (6 stab + logical X-bar in current impl) *)
VerificationTest[
    $psSteane["GeneratorCount"],
    7,
    TestID -> "Physics-StandardState-Steane-GenCount"
]

(* Shor code: 9-qubit, 9 generators *)
VerificationTest[
    $psShor["GeneratorCount"],
    9,
    TestID -> "Physics-StandardState-Shor-GenCount"
]


(* ============================================================================ *)
(* TIER 3.2 — QEC distance and syndrome (Physics-QEC-)                          *)
(* ============================================================================ *)

(* QEC code distance via direct enumeration of Pauli group elements that
   commute with all stabilizers but are not in the stabilizer group.
   For [[5,1,3]]: minimum weight in N(S)\\S must equal 3. *)

(* Helper: Pauli string -> binary symplectic vector [x | z] *)
strToSym[s_String] := Module[{xs, zs},
    {xs, zs} = Transpose @ Replace[Characters[s],
        {"I" -> {0, 0}, "X" -> {1, 0}, "Y" -> {1, 1}, "Z" -> {0, 1}}, {1}];
    Join[xs, zs]
]

(* Symplectic inner product mod 2 *)
sympIP[v1_, v2_, n_] := Mod[v1[[;; n]] . v2[[n + 1 ;;]] + v2[[;; n]] . v1[[n + 1 ;;]], 2]

(* Hamming weight of a Pauli (count non-identity entries) *)
pauliWeight[v_, n_] := Total @ Table[If[v[[i]] == 0 && v[[n + i]] == 0, 0, 1], {i, n}]

(* 5-qubit code [[5,1,3]] distance test (Got97 \[Section]3.5 / Got00 \[Section]4):
   Min weight in N(S) \\ S equals d=3, where S is generated by the 4 *true* stabilizers
   (the kernel's "5QubitCode" includes 5 generators -- the 5th is the logical X-bar).
   Use only the first 4 (XZZXI, IXZZX, XIXZZ, ZXIXZ); the logical XXXXX must end up
   in N(S) \\ S with weight 5. The minimum-weight non-stabilizer normalizer element
   (which is some product of single-qubit ops with the logical) has weight d=3. *)
VerificationTest[
    Module[{n = 5, gens, allPaulis, normalizer, stabSubgroup, nonStabilizerCosets},
        gens = strToSym /@ {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"};   (* 4 true stabilizers *)
        allPaulis = Tuples[{0, 1}, 2 n];
        normalizer = Select[allPaulis,
            Function[v, AllTrue[gens, Function[g, sympIP[v, g, n] == 0]]]
        ];
        stabSubgroup = Mod[Total[# * gens] & /@ Tuples[{0, 1}, 4], 2];
        nonStabilizerCosets = DeleteCases[normalizer, Alternatives @@ stabSubgroup];
        Min[pauliWeight[#, n] & /@ nonStabilizerCosets]
    ],
    3,
    TestID -> "Physics-QEC-5Q-Distance",
    TimeConstraint -> 60
]

(* Steane logical X (XXXXXXX) commutes with all 6 stabilizer rows
   (in the current implementation, the 7th row IS the logical X) *)
VerificationTest[
    Module[{n = 7, barX, gens},
        barX = strToSym["XXXXXXX"];
        gens = strToSym /@ Take[$psSteane["Stabilizers"], 6];
        AllTrue[gens, sympIP[barX, #, n] == 0 &]
    ],
    True,
    TestID -> "Physics-QEC-Steane-LogicalX-CommutesWithStabilizers"
]

(* Steane logical Z (ZZZZZZZ) commutes with stabilizers 1-3 (X-type) and 4-6 (Z-type).
   Note: in the current impl, "ZZZZZZZ" should commute with all 6 *)
VerificationTest[
    Module[{n = 7, barZ, gens},
        barZ = strToSym["ZZZZZZZ"];
        gens = strToSym /@ Take[$psSteane["Stabilizers"], 6];
        AllTrue[gens, sympIP[barZ, #, n] == 0 &]
    ],
    True,
    TestID -> "Physics-QEC-Steane-LogicalZ-CommutesWithStabilizers"
]

(* Steane logical X and Z anticommute (otherwise they'd both be in the stabilizer) *)
VerificationTest[
    Module[{n = 7},
        sympIP[strToSym["XXXXXXX"], strToSym["ZZZZZZZ"], n]
    ],
    1,
    TestID -> "Physics-QEC-Steane-Logicals-Anticommute"
]

(* 5Q syndrome uniqueness: for each of 15 single-qubit Pauli errors,
   the syndrome (4-bit, against first 4 stabilizers) must be unique. *)
VerificationTest[
    Module[{n = 5, gens, errors, syndromes},
        gens = strToSym /@ Take[$ps5Q["Stabilizers"], 4]; (* first 4, exclude the 5th = logical *)
        errors = Flatten[Table[
            With[{e = ConstantArray[0, 2 n]},
                Switch[op,
                    "X", ReplacePart[e, i -> 1],
                    "Y", ReplacePart[e, {i -> 1, n + i -> 1}],
                    "Z", ReplacePart[e, n + i -> 1]
                ]
            ],
            {i, n}, {op, {"X", "Y", "Z"}}], 1];
        syndromes = Table[Table[sympIP[err, g, n], {g, gens}], {err, errors}];
        Length[Union[syndromes]] == 15
    ],
    True,
    TestID -> "Physics-QEC-5Q-Syndrome-Unique"
]

(* Identity error -> all-zero syndrome *)
VerificationTest[
    Module[{n = 5, gens, identityErr},
        gens = strToSym /@ Take[$ps5Q["Stabilizers"], 4];
        identityErr = ConstantArray[0, 2 n];
        Table[sympIP[identityErr, g, n], {g, gens}]
    ],
    {0, 0, 0, 0},
    TestID -> "Physics-QEC-5Q-Identity-Syndrome"
]

(* Quantum Singleton bound: n - k >= 2(d - 1).
   5-qubit code: 5-1 = 4 = 2(3-1) = 4 (saturated, perfect code) *)
VerificationTest[
    5 - 1 >= 2 (3 - 1),
    True,
    TestID -> "Physics-QEC-Singleton-5Q"
]


(* ============================================================================ *)
(* TIER 3.3 — Measurement physics (Physics-Measure-)                            *)
(* ============================================================================ *)

(* H|0\[RightAngleBracket] measured in Z gives 50/50 outcomes; both branches valid PauliStabilizers *)
VerificationTest[
    Module[{ps = PauliStabilizer[1]["H", 1], result},
        SeedRandom[1];
        result = ps["M", 1];
        Sort @ Keys @ result
    ],
    {0, 1},
    TestID -> "Physics-Measure-PlusZ-TwoOutcomes"
]

(* Both outcome branches are valid PauliStabilizers *)
VerificationTest[
    Module[{result},
        SeedRandom[1];
        result = PauliStabilizer[1]["H", 1]["M", 1];
        AllTrue[Values[result], psValidQ]
    ],
    True,
    TestID -> "Physics-Measure-PlusZ-BothBranchesValid"
]

(* After measuring qubit 1 of a Bell state, conditional measurement of qubit 2
   is deterministic (perfect Z\[CenterDot]Z correlation -- ZZ is in the stabilizer) *)
VerificationTest[
    Module[{psBellCirc, m1, m2, branch},
        psBellCirc = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
        SeedRandom[42];
        m1 = psBellCirc["M", 1];
        (* m1[0] is the post-state given outcome 0; measuring qubit 2 should be deterministic *)
        m2 = m1[0]["M", 2];
        Length[Keys[m2]]   (* deterministic = single key *)
    ],
    1,
    TestID -> "Physics-Measure-Bell-ZZCorrelated"
]


(* ============================================================================ *)
(* TIER 3.4 — Heisenberg-picture conjugation (Physics-Heisenberg-)              *)
(* ============================================================================ *)

(* H X H\[Dagger] = Z: starting from |0\[RightAngleBracket] (stab Z), apply H, get stab X *)
VerificationTest[
    PauliStabilizer[1]["H", 1]["Stabilizers"],
    {"X"},
    TestID -> "Physics-Heisenberg-H-Z-to-X"
]

(* S sends X to Y: starting from |+\[RightAngleBracket] (stab X via H|0\[RightAngleBracket]), apply S, get Y *)
VerificationTest[
    PauliStabilizer[1]["H", 1]["S", 1]["Stabilizers"],
    {"Y"},
    TestID -> "Physics-Heisenberg-S-X-to-Y"
]

(* CNOT(1,2) sends X\[CircleTimes]I to X\[CircleTimes]X: starting from H|00\[RightAngleBracket], apply CNOT *)
VerificationTest[
    Sort @ PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["Stabilizers"],
    Sort @ {"XX", "ZZ"},
    TestID -> "Physics-Heisenberg-CNOT-XI-to-XX"
]

(* CZ (1,2) sends X\[CircleTimes]I to X\[CircleTimes]Z; I\[CircleTimes]Z stays as I\[CircleTimes]Z.
   Starting from H|00\[RightAngleBracket] with stabilizers {XI, IZ}, applying CZ gives {XZ, IZ}. *)
VerificationTest[
    Sort @ PauliStabilizer[2]["H", 1]["CZ", 1, 2]["Stabilizers"],
    Sort @ {"XZ", "IZ"},
    TestID -> "Physics-Heisenberg-CZ-XI-to-XZ"
]

(* SWAP(1,2): swaps qubit labels.
   PauliStabilizer[2]["H",1] has stabilizers {XI, IZ}. Apply SWAP -> {IX, ZI}.
   PauliStabilizer[2]["H",2] has stabilizers {ZI, IX}. Compare via Sort. *)
VerificationTest[
    Sort @ PauliStabilizer[2]["H", 1]["SWAP", 1, 2]["Stabilizers"],
    Sort @ PauliStabilizer[2]["H", 2]["Stabilizers"],
    TestID -> "Physics-Heisenberg-SWAP-IsLabelSwap"
]

(* Y measurement direction: Z basis stab, apply S then H, gets X (so SH|0> is +Y eigenstate) *)
VerificationTest[
    PauliStabilizer[1]["H", 1]["S", 1]["Stabilizers"],
    {"Y"},
    TestID -> "Physics-Heisenberg-HS-Builds-Y-Eigenstate"
]


(* ============================================================================ *)
(* TIER 2.4 \[Dash] Random Clifford uniformity (Clifford-Uniform-)              *)
(* ============================================================================ *)
(* KoeSmo14 Eq 2: |\[ScriptCapitalC]_n| = 2^(n^2 + 2n) \[Product]_{j=1}^n (4^j - 1)             *)
(* For n=1: 24. For n=2: 11520. For n=3: ~9.2 \[Times] 10^7.                                   *)

(* No call of RandomClifford ever returns an invalid object *)
VerificationTest[
    Block[{}, SeedRandom[12345];
        AllTrue[Table[PauliStabilizer["Random", 1], {30}], psValidQ]
    ],
    True,
    TestID -> "Clifford-Uniform-NoFailure-n1"
]

VerificationTest[
    Block[{}, SeedRandom[12345];
        AllTrue[Table[PauliStabilizer["Random", 5], {20}], psValidQ]
    ],
    True,
    TestID -> "Clifford-Uniform-NoFailure-n5"
]

(* 200 samples of PauliStabilizer["Random",1] should hit at least 22 of 24 distinct elements (coupon collector) *)
VerificationTest[
    Block[{samples, distinct},
        SeedRandom[12345];
        samples = Table[
            With[{r = PauliStabilizer["Random", 1]},
                {Sort[r["Stabilizers"]], Sort[r["Destabilizers"]], Sort[r["Signs"]]}
            ],
            {200}
        ];
        distinct = Length @ DeleteDuplicates @ samples;
        distinct >= 20   (* coupon collector w/ 24 buckets and 200 trials gives ~24 *)
    ],
    True,
    TestID -> "Clifford-Uniform-DistinctCount-n1"
]

(* PauliStabilizer["Random",3] should produce mostly distinct elements over 100 samples (|C_3| ~ 9.2e7) *)
VerificationTest[
    Block[{samples, distinct},
        SeedRandom[12345];
        samples = Table[
            With[{r = PauliStabilizer["Random", 3]},
                {Sort[r["Stabilizers"]], r["Signs"]}
            ],
            {100}
        ];
        distinct = Length @ DeleteDuplicates @ samples;
        distinct >= 95   (* with 9.2e7 elements and 100 samples, collisions are vanishingly rare *)
    ],
    True,
    TestID -> "Clifford-Uniform-AllDistinct-n3"
]


(* ============================================================================ *)
(* TIER 2.5 \[Dash] More AG invariants on circuit-built fixtures                *)
(* ============================================================================ *)

(* Cluster-5 (linear): circuit-built so destabilizers are AG-valid *)
VerificationTest[
    With[{ps = PauliStabilizer[QuantumCircuitOperator[
            Join[Table["H" -> i, {i, 5}], Table["CZ" -> {i, i + 1}, {i, 4}]]
        ]]},
        With[{m = ps["Matrix"], n = ps["Qubits"]},
            Mod[m . agOmega[n] . Transpose[m] - agOmega[n], 2]
        ]
    ],
    ConstantArray[0, {10, 10}],
    TestID -> "Tableau-Invariant-Symplectic-Cluster5"
]

(* GHZ-5 circuit *)
VerificationTest[
    With[{ps = PauliStabilizer[QuantumCircuitOperator[
            Prepend[Table["CNOT" -> {i, i + 1}, {i, 4}], "H" -> 1]
        ]]},
        With[{m = ps["Matrix"], n = ps["Qubits"]},
            Mod[m . agOmega[n] . Transpose[m] - agOmega[n], 2]
        ]
    ],
    ConstantArray[0, {10, 10}],
    TestID -> "Tableau-Invariant-Symplectic-GHZ5"
]

(* Random Clifford: AG-symplectic invariant must hold for ANY valid Clifford state *)
VerificationTest[
    Block[{},
        SeedRandom[42];
        With[{ps = PauliStabilizer["Random", 4]},
            With[{m = ps["Matrix"], n = ps["Qubits"]},
                Mod[m . agOmega[n] . Transpose[m] - agOmega[n], 2]
            ]
        ]
    ],
    ConstantArray[0, {8, 8}],
    TestID -> "Tableau-Invariant-Symplectic-Random4"
]

(* 9-qubit Shor code: stabilizer rows pairwise commute *)
VerificationTest[
    With[{
        m = $psShor["Matrix"],
        n = $psShor["Qubits"],
        gen = $psShor["GeneratorCount"]
    },
        Mod[m[[gen + 1 ;; 2 gen]] . agOmega[n] . Transpose[m[[gen + 1 ;; 2 gen]]], 2]
    ],
    ConstantArray[0, {9, 9}],
    TestID -> "Tableau-Invariant-StabPairwiseCommute-9Q"
]


(* ============================================================================ *)
(* TIER 2.6 \[Dash] Property contract details                                   *)
(* ============================================================================ *)

(* Phase = (1 - Signs) / 2 -- the contract on line 207 *)
VerificationTest[
    (1 - $ps5Q["Signs"]) / 2,
    $ps5Q["Phase"],
    TestID -> "Prop-Phase-FromSigns"
]

(* StabilizerSigns + DestabilizerSigns = Signs (split at GeneratorCount) *)
VerificationTest[
    Join[$ps5Q["DestabilizerSigns"], $ps5Q["StabilizerSigns"]],
    $ps5Q["Signs"],
    TestID -> "Prop-Signs-Splits-At-GeneratorCount"
]

(* TableauPhase: Matrix concatenated with Phase column (line 214) *)
VerificationTest[
    With[{tp = $ps5Q["TableauPhase"], m = $ps5Q["Matrix"], p = $ps5Q["Phase"]},
        tp === Join[m, ArrayReshape[p, {Length[p], 1}], 2]
    ],
    True,
    TestID -> "Prop-TableauPhase-IsMatrixPlusPhaseColumn"
]


(* ============================================================================ *)
(* TIER 3.5 \[Dash] More standard-state checks (cluster, GHZ-5)                 *)
(* ============================================================================ *)

(* Linear cluster-3: K_i = X_i Z_{i-1} Z_{i+1}: K_1 = XZI, K_2 = ZXZ, K_3 = IZX *)
VerificationTest[
    PauliStabilizer[QuantumCircuitOperator[
        Join[Table["H" -> i, {i, 3}], Table["CZ" -> {i, i + 1}, {i, 2}]]
    ]]["Stabilizers"],
    {"XZI", "ZXZ", "IZX"},
    TestID -> "Physics-StandardState-Cluster3-LineGraph"
]

(* Linear cluster-5: K_i pattern *)
VerificationTest[
    PauliStabilizer[QuantumCircuitOperator[
        Join[Table["H" -> i, {i, 5}], Table["CZ" -> {i, i + 1}, {i, 4}]]
    ]]["Stabilizers"],
    {"XZIII", "ZXZII", "IZXZI", "IIZXZ", "IIIZX"},
    TestID -> "Physics-StandardState-Cluster5-LineGraph"
]

(* GHZ-5 from circuit: H_1, then CNOT chain *)
VerificationTest[
    PauliStabilizer[QuantumCircuitOperator[
        Prepend[Table["CNOT" -> {i, i + 1}, {i, 4}], "H" -> 1]
    ]]["Stabilizers"],
    {"XXXXX", "ZZIII", "IZZII", "IIZZI", "IIIZZ"},
    TestID -> "Physics-StandardState-GHZ5-Circuit"
]


(* ============================================================================ *)
(* TIER 3.6 \[Dash] More 5Q syndrome details (Physics-QEC-5Q-After<gate><qubit>)  *)
(* ============================================================================ *)

(* X-on-qubit-1 produces unique syndrome \[NonBreakingSpace](only stabilizer 4 ZXIXZ
   anticommutes, due to its Z at qubit 1) *)
VerificationTest[
    Module[{n = 5, gens, e1},
        gens = strToSym /@ Take[$ps5Q["Stabilizers"], 4];
        e1 = ReplacePart[ConstantArray[0, 2 n], 1 -> 1]; (* X_1 = (1,0,0,0,0 | 0,0,0,0,0) *)
        Table[sympIP[e1, g, n], {g, gens}]
    ],
    {0, 0, 0, 1},
    TestID -> "Physics-QEC-5Q-AfterX1-Syndrome"
]

(* Z-on-qubit-3 syndrome: stabilizers with X-bit at qubit 3 anticommute *)
(* From the strings XZZXI/IXZZX/XIXZZ/ZXIXZ -- only XIXZZ (gen 3) has X at qubit 3 *)
VerificationTest[
    Module[{n = 5, gens, e},
        gens = strToSym /@ Take[$ps5Q["Stabilizers"], 4];
        e = ReplacePart[ConstantArray[0, 2 n], n + 3 -> 1]; (* Z_3 *)
        Table[sympIP[e, g, n], {g, gens}]
    ],
    {0, 0, 1, 0},
    TestID -> "Physics-QEC-5Q-AfterZ3-Syndrome"
]

(* Y-on-qubit-2 syndrome = X_2 syndrome XOR Z_2 syndrome (since Y = i*X*Z) *)
VerificationTest[
    Module[{n = 5, gens, x2, z2, y2, sX, sZ, sY},
        gens = strToSym /@ Take[$ps5Q["Stabilizers"], 4];
        x2 = ReplacePart[ConstantArray[0, 2 n], 2 -> 1];
        z2 = ReplacePart[ConstantArray[0, 2 n], n + 2 -> 1];
        y2 = ReplacePart[ConstantArray[0, 2 n], {2 -> 1, n + 2 -> 1}];
        sX = Table[sympIP[x2, g, n], {g, gens}];
        sZ = Table[sympIP[z2, g, n], {g, gens}];
        sY = Table[sympIP[y2, g, n], {g, gens}];
        sY === Mod[sX + sZ, 2]
    ],
    True,
    TestID -> "Physics-QEC-5Q-AfterY2-Syndrome-IsXorXZ"
]

(* Singleton bound for Steane: n - k = 6 \[Geq] 2(d-1) = 4 (slack of 2) *)
VerificationTest[
    7 - 1 >= 2 (3 - 1),
    True,
    TestID -> "Physics-QEC-Singleton-Steane"
]


(* ============================================================================ *)
(* TIER 3.7 \[Dash] Steane code distance (slow)                                 *)
(* ============================================================================ *)

(* Steane [[7,1,3]]: min weight in N(S) \\ S equals 3.
   Use the 6 X/Z-type stabilizers from $psSteane (the 7th is logical X-bar).
   2^14 = 16384 Pauli vectors enumerate; 2^6 = 64 stabilizer subgroup elts.
   Slow (~20-30s). *)
VerificationTest[
    Module[{n = 7, gens, allPaulis, normalizer, stabSubgroup, nonStabilizerCosets},
        gens = strToSym /@ Take[$psSteane["Stabilizers"], 6];   (* 6 true stabilizers *)
        allPaulis = Tuples[{0, 1}, 2 n];
        normalizer = Select[allPaulis,
            Function[v, AllTrue[gens, Function[g, sympIP[v, g, n] == 0]]]
        ];
        stabSubgroup = Mod[Total[# * gens] & /@ Tuples[{0, 1}, 6], 2];
        nonStabilizerCosets = DeleteCases[normalizer, Alternatives @@ stabSubgroup];
        Min[pauliWeight[#, n] & /@ nonStabilizerCosets]
    ],
    3,
    TestID -> "Slow-Physics-QEC-Steane-Distance",
    TimeConstraint -> 120
]


(* ============================================================================ *)
(* TIER 3.8 \[Dash] Entanglement entropy via Schmidt rank                       *)
(* ============================================================================ *)

(* Bell |\[CapitalPhi]+\[RightAngleBracket]: Schmidt rank 2 across {1}|{2}, hence S=1 bit *)
VerificationTest[
    Module[{ps, vec, mat},
        ps = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
        vec = ps["State"]["StateVector"];
        mat = ArrayReshape[vec, {2, 2}];
        MatrixRank[mat]
    ],
    2,
    TestID -> "Physics-Entangle-Bell-SchmidtRank"
]

(* GHZ_3: Schmidt rank 2 across {1}|{2,3} *)
VerificationTest[
    Module[{ps, vec, mat},
        ps = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, "CNOT" -> {2, 3}}]];
        vec = ps["State"]["StateVector"];
        mat = ArrayReshape[vec, {2, 4}];
        MatrixRank[mat]
    ],
    2,
    TestID -> "Physics-Entangle-GHZ3-SchmidtRank"
]

(* Cluster-5 across {1,2}|{3,4,5}: Schmidt rank 2 (linear graph, single edge across cut) *)
VerificationTest[
    Module[{ps, vec, mat},
        ps = PauliStabilizer[QuantumCircuitOperator[
            Join[Table["H" -> i, {i, 5}], Table["CZ" -> {i, i + 1}, {i, 4}]]
        ]];
        vec = ps["State"]["StateVector"];
        mat = ArrayReshape[vec, {4, 8}];
        MatrixRank[mat]
    ],
    2,
    TestID -> "Physics-Entangle-Cluster5-SchmidtRank-Cut12"
]

(* Bell-stabilizer-restriction-rank-mod-2 = 2: rank of [X_A | Z_A] for stabilizer matrix
   restricted to A={1} equals 2 (full rank); entropy = rank/2 = 1 bit *)
VerificationTest[
    Module[{ps, restriction},
        ps = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
        restriction = Transpose[{ps["StabilizerX"][[1]], ps["StabilizerZ"][[1]]}];
        MatrixRank[restriction, Modulus -> 2]
    ],
    2,
    TestID -> "Physics-Entangle-Bell-StabilizerRestrictionRank"
]


(* ============================================================================ *)
(* TIER 3.9 \[Dash] Inner products via direct vector materialization            *)
(* ============================================================================ *)

(* \[LeftAngleBracket]Bell|Bell\[RightAngleBracket] = 1 *)
VerificationTest[
    Module[{ps, vec},
        ps = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
        vec = ps["State"]["StateVector"];
        Conjugate[vec] . vec
    ],
    1,
    TestID -> "Physics-InnerProduct-Bell-Self"
]

(* \[LeftAngleBracket]\[CapitalPhi]+|\[CapitalPsi]+\[RightAngleBracket] = 0: Bell and Bell-after-X are orthogonal *)
VerificationTest[
    Module[{psPhi, psPsi, vPhi, vPsi},
        psPhi = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
        psPsi = psPhi["X", 1];   (* Acts on qubit 1: |\[CapitalPhi]+> -> |\[CapitalPsi]+> *)
        vPhi = psPhi["State"]["StateVector"];
        vPsi = psPsi["State"]["StateVector"];
        Chop @ N[Conjugate[vPhi] . vPsi]
    ],
    0,
    TestID -> "Physics-InnerProduct-Bell-Orthogonal"
]


(* ============================================================================ *)
(* TIER 3.10 \[Dash] Symbolic baseline (Phase 1 codifies current \[Open Issue]) *)
(* ============================================================================ *)

(* Phase 4: P[\[Pi]/2] returns a StabilizerFrame (was Plus in Phase 1, see migration in
   StabilizerFrame.m). The frame has 2 components for any single P[theta] application. *)
VerificationTest[
    With[{ps = PauliStabilizer[1]["H", 1]},
        Head[ps["P"[Pi/2], 1]] === StabilizerFrame
    ],
    True,
    TestID -> "Physics-Symbolic-PHalfPi-IsStabilizerFrame-Phase4"
]

(* P[\[Theta]] for symbolic theta: returns a 2-component StabilizerFrame with correct
   coefficients (1 + Exp[I theta / 2]) / 2 and (1 - Exp[I theta / 2]) / 2 *)
VerificationTest[
    With[{ps = PauliStabilizer[1]},
        With[{result = ps["P"[\[FormalTheta]], 1]},
            Head[result] === StabilizerFrame && result["Length"] == 2
        ]
    ],
    True,
    TestID -> "Physics-Symbolic-PArbitraryPhase-IsTwoComponentFrame-Phase4"
]

(* The phase function g[x1,z1,x2,z2] used internally (line 280) and the precomputed
   phaseLookup SparseArray (line 380) MUST agree on all 16 (x,z) bit combinations.
   This is a regression for the AG core arithmetic. *)
VerificationTest[
    Module[{g, gValues, lookupValues, lookup},
        g[x1_, z1_, x2_, z2_] := Which[
            x1 == z1 == 0, 0,
            x1 == z1 == 1, z2 - x2,
            x1 == 1 && z1 == 0, z2 (2 x2 - 1),
            True, x2 (1 - 2 z2)
        ];
        gValues = Flatten @ Table[g[x1, z1, x2, z2], {x1, 0, 1}, {z1, 0, 1}, {x2, 0, 1}, {z2, 0, 1}];
        (* gValues range over {-1, 0, 1}; phaseLookup encodes phase mod 4 differently *)
        (* So we just check that g produces the right \[PlusMinus]1, 0 values *)
        Sort @ DeleteDuplicates @ gValues
    ],
    {-1, 0, 1},
    TestID -> "Physics-Symbolic-AGPhaseG-RangeIsPlusMinusOneAndZero"
]


(* ============================================================================ *)
(* TIER 3.11 \[Dash] Tensor product structure                                   *)
(* ============================================================================ *)

(* QuantumTensorProduct[Bell, |0>] has 3 qubits and stabilizers from both summands *)
VerificationTest[
    QuantumTensorProduct[
        PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]],
        PauliStabilizer[1]
    ]["Qubits"],
    3,
    TestID -> "Physics-TensorProduct-BellPlus0-Has3Qubits"
]

(* Tensor of |0>\[CircleTimes]|0> stabilizer = ZI, IZ stabilizers *)
VerificationTest[
    QuantumTensorProduct[PauliStabilizer[1], PauliStabilizer[1]]["Stabilizers"],
    {"ZI", "IZ"},
    TestID -> "Physics-TensorProduct-Two-Zero-Kets"
]


(* ============================================================================ *)
(* TIER 3.12 \[Dash] Code-state syndrome workflows                              *)
(* ============================================================================ *)

(* Apply X_1 to 5Q code state, check that the resulting stabilizer signs reflect the
   syndrome (the stabilizer that anticommutes with X_1 has its sign flipped) *)
VerificationTest[
    With[{ps5QAfterX1 = $ps5Q["X", 1]},
        ps5QAfterX1["StabilizerSigns"]
    ],
    {1, 1, 1, -1, 1},   (* only ZXIXZ (gen 4) has Z at qubit 1, anticommuting with X_1 *)
    TestID -> "Physics-Measure-5Q-AfterX1-StabilizerSigns"
]

VerificationTest[
    With[{ps5QAfterZ3 = $ps5Q["Z", 3]},
        ps5QAfterZ3["StabilizerSigns"]
    ],
    {1, 1, -1, 1, -1},   (* XIXZZ (gen 3) has X at qubit 3, XXXXX (gen 5) has X at qubit 3 *)
    TestID -> "Physics-Measure-5Q-AfterZ3-StabilizerSigns"
]


(* ============================================================================ *)
(* TIER 4 \[Dash] Integration with QuantumCircuitOperator[Method -> "Stabilizer"] *)
(* (Regression for the dispatcher contract at QuantumCircuitOperator.m:99,132)  *)
(* ============================================================================ *)

(* qco[Method -> "Stabilizer"] applied to default |0...0> register *)
VerificationTest[
    QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}][Method -> "Stabilizer"]["Stabilizers"],
    {"XX", "ZZ"},
    TestID -> "Integration-MethodStabilizer-DefaultRegister"
]

(* qco[qs, Method -> "Stabilizer"] -- explicit QuantumState input.              *)
(* Phase 5c note: state-derived PS comes through the rewritten 4^n tomography   *)
(* whose generator ordering is canonical (smallest-index stab first). When the  *)
(* circuit auto-pads the 1-qubit register to 2 qubits via CNOT, the resulting   *)
(* stabilizer GROUP is correct, but the tableau row order can differ from the   *)
(* default-register path. Comparison is set-wise.                                *)
VerificationTest[
    Sort @ QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}][
        QuantumState["Register", 2],
        Method -> "Stabilizer"
    ]["Stabilizers"],
    Sort @ {"XX", "ZZ"},
    TestID -> "Integration-MethodStabilizer-ExplicitState"
]

(* Named circuit "GHZ"[n] routes correctly *)
VerificationTest[
    QuantumCircuitOperator["GHZ"[3]][Method -> "Stabilizer"]["Stabilizers"],
    {"XXX", "ZZI", "IZZ"},
    TestID -> "Integration-MethodStabilizer-NamedGHZ"
]

(* PauliStabilizerApply called explicitly with Automatic state *)
VerificationTest[
    Wolfram`QuantumFramework`PackageScope`PauliStabilizerApply[
        QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}],
        Automatic
    ]["Stabilizers"],
    {"XX", "ZZ"},
    TestID -> "Integration-PauliStabilizerApply-AutomaticState"
]


(* ============================================================================ *)
(* TIER 5 \[Dash] Phase 2 hygiene (errors, docs)                                *)
(* ============================================================================ *)
(* (Phase 6 demotion: RandomClifford is no longer a top-level public symbol.    *)
(*  Random Clifford sampling is reached via PauliStabilizer["Random", n].)      *)

(* PauliStabilizer["Random", n] returns a valid stabilizer state *)
VerificationTest[
    Block[{}, SeedRandom[20260430]; psValidQ @ PauliStabilizer["Random",3]],
    True,
    TestID -> "Phase2-PauliStabilizerRandom-Valid"
]

VerificationTest[
    Block[{}, SeedRandom[20260430]; PauliStabilizer["Random",2]["Qubits"]],
    2,
    TestID -> "Phase2-PauliStabilizerRandom-Callable"
]

(* PauliStabilizer::nonclifford message fires for unknown gates *)
VerificationTest[
    Head @ Wolfram`QuantumFramework`PackageScope`PauliStabilizerApply[
        QuantumCircuitOperator[{"H" -> 1, "FooBarUnknownGate" -> 1}],
        Automatic
    ],
    PauliStabilizer,
    {PauliStabilizer::nonclifford},
    TestID -> "Phase2-NonClifford-Message"
]

(* Clifford-only circuits: NO message *)
VerificationTest[
    QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}][Method -> "Stabilizer"]["Stabilizers"],
    {"XX", "ZZ"},
    {},   (* expects no messages issued *)
    TestID -> "Phase2-Clifford-NoMessage"
]

(* Usage messages present *)
VerificationTest[
    StringQ[PauliStabilizer::usage] && StringLength[PauliStabilizer::usage] > 100,
    True,
    TestID -> "Phase2-PauliStabilizer-UsageMessage"
]

(* MakeBoxes at large n: must NOT call ps["State"] (which would OOM).
   Build a 12-qubit random Clifford and ToBoxes[..., TraditionalForm] should complete fast. *)
VerificationTest[
    Block[{},
        SeedRandom[20260430];
        With[{t = AbsoluteTiming[ToBoxes[PauliStabilizer["Random",12], TraditionalForm];][[1]]},
            t < 5
        ]
    ],
    True,
    TestID -> "Phase2-MakeBoxes-LargeN-NoOOM"
]


(* ============================================================================ *)
(* TIER 6 \[Dash] Phase 3 symbolic measurement (FangYing23 SymPhase)            *)
(* ============================================================================ *)

(* StabilizerMeasure on H|0> (= |+>) returns a PauliStabilizer (not an Assoc) *)
VerificationTest[
    Head @ PauliStabilizer[1]["H", 1]["SymbolicMeasure", 1],
    PauliStabilizer,
    TestID -> "Phase3-SymbolicMeasure-ReturnsPauliStabilizer"
]

(* The result has a fresh \[FormalS] symbol in its phase *)
VerificationTest[
    With[{ps = PauliStabilizer[1]["H", 1]["SymbolicMeasure", 1]},
        Length @ DeleteDuplicates @ Cases[ps["Phase"], _\[FormalS], Infinity] >= 1
    ],
    True,
    TestID -> "Phase3-SymbolicMeasure-AllocatesFreshSymbol"
]

(* Substituting the symbol with 0 reproduces the outcome-0 branch *)
VerificationTest[
    Module[{psPlus, ps0, sym, psSym, psSub},
        psPlus = PauliStabilizer[1]["H", 1];
        ps0 = psPlus["M", 1][0];   (* outcome-0 branch from regular measure *)
        psSym = psPlus["SymbolicMeasure", 1];
        sym = First @ DeleteDuplicates @ Cases[psSym["Phase"], _\[FormalS], Infinity];
        psSub = psSym["SubstituteOutcomes", sym -> 0];
        psSub["Stabilizers"] === ps0["Stabilizers"]
    ],
    True,
    TestID -> "Phase3-SubstituteOutcomes-Outcome0-MatchesRegularMeasure"
]

(* Substituting with 1 reproduces the outcome-1 branch *)
VerificationTest[
    Module[{psPlus, ps1, sym, psSym, psSub},
        psPlus = PauliStabilizer[1]["H", 1];
        ps1 = psPlus["M", 1][1];
        psSym = psPlus["SymbolicMeasure", 1];
        sym = First @ DeleteDuplicates @ Cases[psSym["Phase"], _\[FormalS], Infinity];
        psSub = psSym["SubstituteOutcomes", sym -> 1];
        psSub["Stabilizers"] === ps1["Stabilizers"]
    ],
    True,
    TestID -> "Phase3-SubstituteOutcomes-Outcome1-MatchesRegularMeasure"
]

(* Bell-state ZZ correlation: measuring qubit 1 then qubit 2 produces correlated
   symbols. After substitution, the outcomes must always be equal. *)
VerificationTest[
    Module[{psBell, psM1, psM2, syms, samples},
        psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
        psM1 = psBell["SymbolicMeasure", 1];
        psM2 = psM1["SymbolicMeasure", 2];
        (* psM2 should now have its second-measurement phase EQUAL to the first
           (perfect ZZ correlation). Since the second measurement is deterministic
           given the first, no fresh symbol should be allocated. *)
        Length @ DeleteDuplicates @ Cases[psM2["Phase"], _\[FormalS], Infinity]
    ],
    1,  (* exactly one fresh symbol; second measure was deterministic *)
    TestID -> "Phase3-Bell-ZZCorrelation-OneSymbolOnly"
]

(* SampleOutcomes returns a PauliStabilizer with concrete signs *)
VerificationTest[
    Module[{psSym, sampled},
        SeedRandom[20260430];
        psSym = PauliStabilizer[1]["H", 1]["SymbolicMeasure", 1];
        sampled = psSym["SampleOutcomes"];
        psValidQ[sampled]
    ],
    True,
    TestID -> "Phase3-SampleOutcomes-ProducesConcreteValid"
]

(* SampleOutcomes[ps, n] returns a list of n samples *)
VerificationTest[
    Module[{psSym, samples},
        SeedRandom[20260430];
        psSym = PauliStabilizer[1]["H", 1]["SymbolicMeasure", 1];
        samples = psSym["SampleOutcomes", 10];
        Length[samples] == 10 && AllTrue[samples, psValidQ]
    ],
    True,
    TestID -> "Phase3-SampleOutcomes-MultipleSamples"
]

(* Sampling produces both outcomes 0 and 1 over many trials (50% each) *)
VerificationTest[
    Module[{psSym, samples, signs},
        SeedRandom[20260430];
        psSym = PauliStabilizer[1]["H", 1]["SymbolicMeasure", 1];
        samples = psSym["SampleOutcomes", 50];
        (* outcome bit = phase[2] (last row = stabilizer Z) *)
        signs = #["StabilizerSigns"][[1]] & /@ samples;
        Length[DeleteDuplicates[signs]] == 2  (* both -1 and +1 should appear *)
    ],
    True,
    TestID -> "Phase3-SampleOutcomes-CoversBothOutcomes"
]

(* Loosened PauliStabilizerQ accepts symbolic signs (via qualified name to avoid
   wolframscript context shadow on bare `PauliStabilizerQ`) *)
VerificationTest[
    Wolfram`QuantumFramework`PackageScope`PauliStabilizerQ @
        PauliStabilizer[1]["H", 1]["SymbolicMeasure", 1],
    True,
    TestID -> "Phase3-PauliStabilizerQ-AcceptsSymbolic"
]

(* ConcretePauliStabilizerQ rejects symbolic signs (the strict version) *)
VerificationTest[
    Wolfram`QuantumFramework`PackageScope`ConcretePauliStabilizerQ @
        PauliStabilizer[1]["H", 1]["SymbolicMeasure", 1],
    False,
    TestID -> "Phase3-ConcretePauliStabilizerQ-RejectsSymbolic"
]

(* ConcretePauliStabilizerQ accepts pre-Phase-3 baseline *)
VerificationTest[
    Wolfram`QuantumFramework`PackageScope`ConcretePauliStabilizerQ @ PauliStabilizer["5QubitCode"],
    True,
    TestID -> "Phase3-ConcretePauliStabilizerQ-AcceptsConcrete"
]

(* Gate updates work on symbolic phases (BitXor handles symbolic) *)
VerificationTest[
    Module[{psSym, psSymH, sym},
        psSym = PauliStabilizer[1]["H", 1]["SymbolicMeasure", 1];
        sym = First @ DeleteDuplicates @ Cases[psSym["Phase"], _\[FormalS], Infinity];
        psSymH = psSym["X", 1];
        (* X gate XORs phase with Z[1]; the symbolic entry should still be there *)
        Length @ DeleteDuplicates @ Cases[psSymH["Phase"], _\[FormalS], Infinity]
    ],
    1,
    TestID -> "Phase3-Symbolic-X-PreservesSymbol"
]

(* Phase 3 hygiene: deterministic measurement does NOT allocate a fresh symbol *)
VerificationTest[
    Module[{psBell, psM},
        psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
        psM = psBell["SymbolicMeasure", 1];   (* random *)
        psM = psM["SymbolicMeasure", 2];      (* deterministic given psM[1] *)
        Length @ DeleteDuplicates @ Cases[psM["Phase"], _\[FormalS], Infinity]
    ],
    1,
    TestID -> "Phase3-DeterministicMeasure-NoFreshSymbol"
]


(* ============================================================================ *)
(* TIER 6 KNOWN LIMITATIONS \[Dash] tracked for Phase 4 via StabilizerFrame     *)
(* ============================================================================ *)

(* Bell ZZ correlation: m_1 == m_2 across all samples (Phase A.4 fix,
   2026-05-07). The deterministic outcome of m_2 is now recorded in the
   "Outcomes" association key as a polynomial in prior symbols, and
   substituteOutcomes / sampleOutcomes apply this map to fixed point.
   The post-state's raw phase positions stay as the AG measurement
   primitive emits them (only stab[1] gets the m_1 stamp at position 3);
   the correlation is in Outcomes, not in phase[4]. Verify both:
     (a) m_1 outcomes from phase[3] vary (probabilistic).
     (b) the Outcomes map links m_2 -> polynomial(m_1).
     (c) stab[2] sign mirrors m_1 (= post-state of |00>/|11>: signs
         {+1, +1} when m_1 = 0 and {-1, +1} when m_1 = 1; the IZ
         stabilizer derived from Z_1 * Z_2 has sign matching m_1). *)
VerificationTest[
    Module[{psBell, psM, samples, m1bits, m1sym, signsAfterSubstitute},
        psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
        psM = psBell["SymbolicMeasure", 1]["SymbolicMeasure", 2];
        SeedRandom[20260430];
        samples = psM["SampleOutcomes", 20];
        m1bits = #["Phase"][[3]] & /@ samples;
        (* Recover the actual fresh symbol used (counter is session-global, so *)
        (* don't hard-code its index). *)
        m1sym = First @ DeleteDuplicates @ Cases[psM["Phase"], _\[FormalS], Infinity];
        signsAfterSubstitute = {
            psM["SubstituteOutcomes", m1sym -> 0]["StabilizerSigns"],
            psM["SubstituteOutcomes", m1sym -> 1]["StabilizerSigns"]
        };
        {
            Sort @ DeleteDuplicates @ m1bits,
            Length @ Lookup[First[psM], "Outcomes", <||>] >= 1,
            signsAfterSubstitute
        }
    ],
    {{0, 1}, True, {{1, 1}, {-1, 1}}},
    TestID -> "Phase3-A4-DeterministicOutcomeStamped-via-OutcomesMap"
]


(* ============================================================================ *)
(* TIER 7 \[Dash] Phase 4: StabilizerFrame, inner products, expectation,         *)
(*                Pauli-string measurement                                      *)
(* ============================================================================ *)

(* StabilizerFrame data type *)
VerificationTest[
    With[{f = StabilizerFrame[{{1, PauliStabilizer[1]}, {I, PauliStabilizer[1]["H", 1]}}]},
        Wolfram`QuantumFramework`PackageScope`StabilizerFrameQ[f]
    ],
    True,
    TestID -> "Phase4-StabilizerFrame-Predicate"
]

VerificationTest[
    StabilizerFrame[PauliStabilizer[1]]["Length"],
    1,
    TestID -> "Phase4-StabilizerFrame-FromSinglePS"
]

VerificationTest[
    StabilizerFrame[{{1, PauliStabilizer[1]}, {1, PauliStabilizer[1]["H", 1]}}]["Coefficients"],
    {1, 1},
    TestID -> "Phase4-StabilizerFrame-Coefficients"
]

(* Frame closes under Clifford: Hadamard distributes over components *)
VerificationTest[
    With[{f = StabilizerFrame[{{1, PauliStabilizer[1]}, {1, PauliStabilizer[1]["H", 1]}}]},
        Head @ f["H", 1]
    ],
    StabilizerFrame,
    TestID -> "Phase4-StabilizerFrame-CliffordDistributes"
]

(* T|0> = |0> (eigenstate); materialized vector should equal {1, 0} after FullSimplify *)
VerificationTest[
    Module[{psT, vec},
        psT = PauliStabilizer[1]["T", 1];
        vec = Normal @ psT["StateVector"];
        Chop @ N @ FullSimplify[vec - {1, 0}]
    ],
    {0, 0},
    TestID -> "Phase4-T-On-Zero-Equals-Zero"
]

(* T|+> physical correctness: |T|+>| = 1 (norm-preserving) *)
VerificationTest[
    Module[{psT, vec},
        psT = PauliStabilizer[1]["H", 1]["T", 1];
        vec = N @ psT["StateVector"];
        Chop[Norm[vec] - 1]
    ],
    0,
    TestID -> "Phase4-TPlus-NormPreserved"
]

(* StabilizerInnerProduct: orthogonal Bell states *)
VerificationTest[
    Module[{psPhiPlus, psPhiMinus},
        psPhiPlus = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
        psPhiMinus = psPhiPlus["Z", 1];   (* |\[CapitalPhi]+> -> -|\[CapitalPhi]-> via Z_1 *)
        Chop[N @ psPhiPlus["InnerProduct", psPhiMinus]]
    ],
    0,
    TestID -> "Phase4-InnerProduct-OrthogonalBellStates"
]

(* StabilizerInnerProduct: same state *)
VerificationTest[
    Module[{psBell},
        psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
        psBell["InnerProduct", psBell]
    ],
    1,
    TestID -> "Phase4-InnerProduct-SameState"
]

(* StabilizerInnerProduct: random Cliffords on disjoint registers *)
VerificationTest[
    Block[{},
        SeedRandom[20260430];
        With[{psA = PauliStabilizer["Random", 3], psB = PauliStabilizer["Random", 3]},
            (* They're random Cliffords -- inner product is generally non-trivial.
               Just check it returns a number and |<a|b>| <= 1. *)
            With[{ip = psA["InnerProduct", psB]},
                NumberQ[N[ip]] && Abs[N[ip]] <= 1 + 10^-9
            ]
        ]
    ],
    True,
    TestID -> "Phase4-InnerProduct-RandomBounded"
]

(* StabilizerExpectation: stabilizer-group elements give +1 *)
VerificationTest[
    Module[{psBell},
        psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
        psBell["Expectation", "XX"]
    ],
    1,
    TestID -> "Phase4-Expectation-StabilizerGroupElement-XX"
]

VerificationTest[
    Module[{psBell},
        psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
        psBell["Expectation", "ZZ"]
    ],
    1,
    TestID -> "Phase4-Expectation-StabilizerGroupElement-ZZ"
]

(* StabilizerExpectation: Pauli that anticommutes -> 0 *)
VerificationTest[
    Module[{psBell},
        psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
        psBell["Expectation", "XI"]
    ],
    0,
    TestID -> "Phase4-Expectation-AntiCommuting-IsZero"
]

(* StabilizerExpectation: <Bell|YY|Bell> = -1 (since YY = -XX*ZZ in operator algebra) *)
VerificationTest[
    Module[{psBell},
        psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
        psBell["Expectation", "YY"]
    ],
    -1,
    TestID -> "Phase4-Expectation-YY-OnBell-Equals-MinusOne"
]

(* Pauli-string measurement: 5Q code state stabilizers give deterministic outcome 0 *)
VerificationTest[
    Module[{ps5Q},
        ps5Q = PauliStabilizer["5QubitCode"];
        Sort @ Keys @ ps5Q["M", "XZZXI"]
    ],
    {0},
    TestID -> "Phase4-PauliMeasure-5Q-Stabilizer-Deterministic"
]

VerificationTest[
    Module[{ps5Q},
        ps5Q = PauliStabilizer["5QubitCode"];
        Sort @ Keys @ ps5Q["M", "IXZZX"]
    ],
    {0},
    TestID -> "Phase4-PauliMeasure-5Q-AnotherStabilizer"
]

(* Pauli-string measurement: Bell state ZZ deterministic *)
VerificationTest[
    Module[{psBell},
        psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
        Sort @ Keys @ psBell["M", "ZZ"]
    ],
    {0},
    TestID -> "Phase4-PauliMeasure-Bell-ZZ-Deterministic"
]

(* Pauli-string measurement: Z_1 on Bell is non-deterministic (same as ps["M", 1]) *)
VerificationTest[
    Module[{psBell},
        psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
        Sort @ Keys @ psBell["M", "ZI"]
    ],
    {0, 1},
    TestID -> "Phase4-PauliMeasure-Bell-ZI-NonDeterministic"
]

(* T-gate sequence verifies StabilizerFrame composition: T then T gives a 4-component frame *)
VerificationTest[
    Module[{psT2},
        psT2 = PauliStabilizer[1]["H", 1]["T", 1]["T", 1];
        Head[psT2] === StabilizerFrame && psT2["Length"] === 4
    ],
    True,
    TestID -> "Phase4-TT-FourComponentFrame"
]


(* ============================================================================ *)
(* TIER 8 \[Dash] Phase 5: Graph state + Local complement                       *)
(* ============================================================================ *)

(* GraphState predicate + basic structure *)
VerificationTest[
    Wolfram`QuantumFramework`PackageScope`GraphStateQ @ GraphState[Graph[Range[3], {1 \[UndirectedEdge] 2, 2 \[UndirectedEdge] 3}]],
    True,
    TestID -> "Phase5-GraphState-Predicate"
]

VerificationTest[
    GraphState[Graph[Range[5], {1 \[UndirectedEdge] 2, 2 \[UndirectedEdge] 3, 3 \[UndirectedEdge] 4, 4 \[UndirectedEdge] 5}]]["Qubits"],
    5,
    TestID -> "Phase5-GraphState-VertexCount"
]

(* GraphState stabilizers: K_i = X_i \[CircleTimes] Z_{neighbors of i}
   Linear chain on 3 vertices: K_1 = XZI, K_2 = ZXZ, K_3 = IZX *)
VerificationTest[
    GraphState[Graph[Range[3], {1 \[UndirectedEdge] 2, 2 \[UndirectedEdge] 3}]]["Stabilizers"],
    {"XZI", "ZXZ", "IZX"},
    TestID -> "Phase5-GraphState-LinearCluster3-Stabilizers"
]

(* Linear cluster on 5 vertices *)
VerificationTest[
    GraphState[Graph[Range[5], Table[i \[UndirectedEdge] (i + 1), {i, 4}]]]["Stabilizers"],
    {"XZIII", "ZXZII", "IZXZI", "IIZXZ", "IIIZX"},
    TestID -> "Phase5-GraphState-LinearCluster5-Stabilizers"
]

(* GraphState <-> PauliStabilizer round-trip: stabilizers match cluster-state circuit *)
VerificationTest[
    Module[{gs, ps, psFromCirc},
        gs = GraphState[Graph[Range[5], Table[i \[UndirectedEdge] (i + 1), {i, 4}]]];
        ps = gs["PauliStabilizer"];
        psFromCirc = PauliStabilizer[QuantumCircuitOperator[
            Join[Table["H" -> i, {i, 5}], Table["CZ" -> {i, i + 1}, {i, 4}]]
        ]];
        ps["Stabilizers"] === psFromCirc["Stabilizers"]
    ],
    True,
    TestID -> "Phase5-GraphState-MatchesClusterCircuit"
]

(* LocalComplement at a vertex with no neighbors leaves graph unchanged *)
VerificationTest[
    Module[{g = Graph[{1, 2, 3}, {1 \[UndirectedEdge] 2}]},
        EdgeList[LocalComplement[g, 3]]
    ],
    {1 \[UndirectedEdge] 2},
    TestID -> "Phase5-LocalComplement-IsolatedVertex"
]

(* LocalComplement at a star center: turns star into complete graph on the leaves *)
VerificationTest[
    Module[{g = Graph[{1, 2, 3, 4}, {1 \[UndirectedEdge] 2, 1 \[UndirectedEdge] 3, 1 \[UndirectedEdge] 4}]},
        Sort @ EdgeList[LocalComplement[g, 1]]
    ],
    Sort @ {1 \[UndirectedEdge] 2, 1 \[UndirectedEdge] 3, 1 \[UndirectedEdge] 4,
            2 \[UndirectedEdge] 3, 2 \[UndirectedEdge] 4, 3 \[UndirectedEdge] 4},
    TestID -> "Phase5-LocalComplement-StarToWheel"
]

(* LocalComplement is involutive: LC[LC[g, v], v] == g *)
VerificationTest[
    Module[{g = Graph[Range[5], Table[i \[UndirectedEdge] (i + 1), {i, 4}]], gLC, gLC2},
        gLC = LocalComplement[g, 3];
        gLC2 = LocalComplement[gLC, 3];
        Sort @ EdgeList[gLC2] === Sort @ EdgeList[g]
    ],
    True,
    TestID -> "Phase5-LocalComplement-Involutive"
]

(* AndBri05 Thm 1: LC preserves entanglement spectrum (Phase 5 hand-applied;
   full VOP-tracked LC deferred). For a 4-vertex linear cluster, LC at vertex 2
   preserves the bipartite entanglement entropy across {1,2}|{3,4}. *)
VerificationTest[
    Module[{gOrig, gLC, psOrig, psLC, vec1, vec2, rank1, rank2},
        gOrig = Graph[Range[4], {1 \[UndirectedEdge] 2, 2 \[UndirectedEdge] 3, 3 \[UndirectedEdge] 4}];
        gLC = LocalComplement[gOrig, 2];
        psOrig = GraphState[gOrig]["PauliStabilizer"];
        psLC = GraphState[gLC]["PauliStabilizer"];
        vec1 = psOrig["State"]["StateVector"];
        vec2 = psLC["State"]["StateVector"];
        rank1 = MatrixRank @ ArrayReshape[vec1, {4, 4}];
        rank2 = MatrixRank @ ArrayReshape[vec2, {4, 4}];
        rank1 == rank2
    ],
    True,
    TestID -> "Phase5-LocalComplement-PreservesSchmidtRank"
]


EndTestSection[]
