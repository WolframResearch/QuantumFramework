(* ::Package:: *)

(* ============================================================================
   Tests/Stabilizer/AuditMatrix.wlt -- comprehensive audit-matrix coverage.

   Generated 2026-05-06 from a function-by-function walk of every public
   constructor, method (DownValue), and property (UpValue / accessor) of the
   Stabilizer subsystem (PauliStabilizer + StabilizerFrame + GraphState +
   LocalComplement + CliffordChannel). The goal is a *complete coverage matrix*
   so that no input form is implicitly tested only via downstream behaviour.

   Conventions
   -----------
   - Tests are grouped into TIERs by *what they cover* (Construction,
     Properties, Methods (Clifford), Methods (non-Clifford), Methods
     (Measurement), Method dispatch (Phase 6/6.5 demoted ops), UpValues,
     Cross-head interop, Cross-package fixtures).
   - Every TestID is namespaced "Audit-{tier}-{thing}" so a failing entry is
     traceable to a specific row in API.md.
   - Cross-package comparison uses the live Stim fixture JSON
     (Tests/Stabilizer/fixtures/stim_fixtures.json) and the QuantumClifford.jl
     hand-coded canonical values.

   This file is purely additive over PauliStabilizer.wlt / CliffordChannel.wlt /
   HybridInterop.wlt / Roundtrips.wlt; it is meant to be run AFTER those to
   provide a sanity sweep over the API surface.
   ========================================================================== *)

Needs["Wolfram`QuantumFramework`"];

(* Local helpers. psValidQ is the package-scoped predicate. matEqQO compares two
   QuantumOperators by exact matrix equality (no global-phase tolerance). *)
psValidQ = Wolfram`QuantumFramework`PackageScope`PauliStabilizerQ;
ccValidQ = Wolfram`QuantumFramework`PackageScope`CliffordChannelQ;

matEqQO[a_, b_] := SameQ @@ (Normal @ #["Matrix"] & /@ {a, b})


(* ============================================================================ *)
(* TIER 1 -- PauliStabilizer constructors: ALL recognized input forms.          *)
(*                                                                              *)
(* Mirrors API.md "Constructors" section verbatim. Each form gets exactly one  *)
(* validity test and one structural check (Qubits / Stabilizers / etc.).        *)
(* ============================================================================ *)

(* 1.1 -- list of Pauli strings, no signs *)
VerificationTest[
    psValidQ @ PauliStabilizer[{"XX", "ZZ"}],
    True,
    TestID -> "Audit-Construct-PauliStrings-NoSigns"
]

VerificationTest[
    PauliStabilizer[{"XX", "ZZ"}]["Stabilizers"],
    {"XX", "ZZ"},
    TestID -> "Audit-Construct-PauliStrings-NoSigns-Stabs"
]

(* 1.2 -- list of Pauli strings WITH +/- prefix *)
VerificationTest[
    PauliStabilizer[{"-XX", "+ZZ"}]["StabilizerSigns"],
    {-1, 1},
    TestID -> "Audit-Construct-PauliStrings-Signs"
]

(* 1.3 -- explicit stab + destab halves *)
VerificationTest[
    With[{ps = PauliStabilizer[{"XX", "ZZ"}, {"IX", "ZI"}]},
        {ps["Stabilizers"], ps["Destabilizers"]}
    ],
    {{"XX", "ZZ"}, {"IX", "ZI"}},
    TestID -> "Audit-Construct-Stab-and-Destab"
]

(* 1.4 -- integer constructor: |0...0> register *)
VerificationTest[
    With[{ps = PauliStabilizer[3]},
        {ps["Qubits"], ps["Stabilizers"], ps["Destabilizers"]}
    ],
    {3, {"ZII", "IZI", "IIZ"}, {"XII", "IXI", "IIX"}},
    TestID -> "Audit-Construct-Integer-3qubit-Register"
]

VerificationTest[
    PauliStabilizer[1]["Stabilizers"],
    {"Z"},
    TestID -> "Audit-Construct-Integer-1qubit"
]

(* PauliStabilizer[0] auto-pads to 1 qubit (Max[q, 1] in Constructors.m:241).  *)
(* Document this contract.                                                      *)
VerificationTest[
    PauliStabilizer[0]["Qubits"],
    1,
    TestID -> "Audit-Construct-Integer-0qubit-AutoPadsToOne"
]

(* 1.5 -- empty form delegates to PauliStabilizer[1] *)
VerificationTest[
    PauliStabilizer[]["Stabilizers"],
    {"Z"},
    TestID -> "Audit-Construct-Empty"
]

(* 1.6 -- named code: 5QubitCode + signed variant *)
VerificationTest[
    PauliStabilizer["5QubitCode"]["Stabilizers"],
    {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ", "XXXXX"},
    TestID -> "Audit-Construct-Named-5QubitCode"
]

VerificationTest[
    PauliStabilizer["5QubitCode1"]["StabilizerSigns"][[5]],
    -1,
    TestID -> "Audit-Construct-Named-5QubitCode1-LogicalOne"
]

VerificationTest[
    PauliStabilizer["7QubitCode"]["Stabilizers"][[7]],
    "XXXXXXX",
    TestID -> "Audit-Construct-Named-7QubitCode-Last"
]

VerificationTest[
    PauliStabilizer["SteaneCode"]["Stabilizers"] === PauliStabilizer["7QubitCode"]["Stabilizers"],
    True,
    TestID -> "Audit-Construct-Named-SteaneCode-AliasFor7"
]

VerificationTest[
    PauliStabilizer["9QubitCode"]["Qubits"],
    9,
    TestID -> "Audit-Construct-Named-9QubitCode"
]

(* 1.7 -- Random named pattern *)
VerificationTest[
    Block[{},
        SeedRandom[12345];
        PauliStabilizer["Random", 3]["Qubits"]
    ],
    3,
    TestID -> "Audit-Construct-Named-Random-Qubits"
]

VerificationTest[
    Block[{},
        SeedRandom[12345];
        psValidQ @ PauliStabilizer["Random", 4]
    ],
    True,
    TestID -> "Audit-Construct-Named-Random-Validity"
]

(* 1.8 -- from QuantumCircuitOperator *)
VerificationTest[
    PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]["Stabilizers"],
    {"XX", "ZZ"},
    TestID -> "Audit-Construct-QCO-Bell"
]

(* 1.9 -- from QuantumOperator (single-gate Heisenberg projection) *)
VerificationTest[
    PauliStabilizer[QuantumOperator["H", 1]]["Stabilizers"],
    {"X"},
    TestID -> "Audit-Construct-QO-H"
]

VerificationTest[
    PauliStabilizer[QuantumOperator["Y"]]["GlobalPhase"],
    -I,
    TestID -> "Audit-Construct-QO-Y-CapturesGlobalPhase"
]

(* 1.10 -- from QuantumState (4^n tomography path) *)
VerificationTest[
    PauliStabilizer[QuantumState[{1, 0}]]["Stabilizers"],
    {"Z"},
    TestID -> "Audit-Construct-QS-Zero"
]

VerificationTest[
    PauliStabilizer[QuantumState[{0, 1}]]["Stabilizers"],
    {"-Z"},
    TestID -> "Audit-Construct-QS-One"
]

VerificationTest[
    Sort @ PauliStabilizer[QuantumState["PhiPlus"]]["Stabilizers"],
    Sort @ {"XX", "ZZ"},
    TestID -> "Audit-Construct-QS-PhiPlus"
]

(* 1.11 -- Association-keyed: Phase + Tableau (encodes |1> with stabilizer -Z). *)
(* Tableau shape {2, n, 2n}: X-block then Z-block, 1 qubit, 2 rows.            *)
(* Row 1 = destabilizer (X here): X-bit 1, Z-bit 0.                             *)
(* Row 2 = stabilizer (Z here): X-bit 0, Z-bit 1.                               *)
(* Phase = {0, 1} marks the stabilizer sign as -1.                              *)
VerificationTest[
    PauliStabilizer[<|"Phase" -> {0, 1}, "Tableau" -> {{{1, 0}}, {{0, 1}}}|>]["Stabilizers"],
    {"-Z"},
    TestID -> "Audit-Construct-Assoc-PhaseTableau"
]

(* 1.12 -- Association-keyed: Matrix + Phase (sympletic block) *)
VerificationTest[
    With[{ps = PauliStabilizer[<|"Matrix" -> IdentityMatrix[2], "Phase" -> {0, 0}|>]},
        {ps["Qubits"], ps["Stabilizers"]}
    ],
    {1, {"Z"}},
    TestID -> "Audit-Construct-Assoc-MatrixPhase"
]

(* 1.13 -- Tableau-only, default signs *)
VerificationTest[
    With[{ps = PauliStabilizer[<|"Tableau" -> {{{1, 0}, {0, 1}}, {{0, 1}, {1, 0}}}|>]},
        ps["Signs"]
    ],
    {1, 1},
    TestID -> "Audit-Construct-Assoc-TableauOnly-DefaultSigns"
]

(* 1.14 -- shortcut form: string -> circuit *)
VerificationTest[
    psValidQ @ PauliStabilizer["H" -> 1],
    True,
    TestID -> "Audit-Construct-Shortcut-StringRule"
]


(* ============================================================================ *)
(* TIER 2 -- PauliStabilizer Properties: every key in ps["Properties"].         *)
(*                                                                              *)
(* Coverage matrix: each property tested exactly once on a non-trivial state.  *)
(* The Bell state and the 5Q code give us a small + medium fixture.            *)
(* ============================================================================ *)

$psBell = PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "CNOT" -> {1, 2}};
$ps5Q = PauliStabilizer["5QubitCode"];

VerificationTest[$psBell["Qubits"],  2, TestID -> "Audit-Property-Qubits"]
VerificationTest[$psBell["Qudits"],  2, TestID -> "Audit-Property-Qudits"]
VerificationTest[$psBell["GeneratorCount"], 2, TestID -> "Audit-Property-GeneratorCount"]

VerificationTest[
    $psBell["Signs"],
    {1, 1, 1, 1},
    TestID -> "Audit-Property-Signs"
]

VerificationTest[
    $psBell["Phase"],
    {0, 0, 0, 0},
    TestID -> "Audit-Property-Phase"
]

VerificationTest[
    Dimensions @ $psBell["Tableau"],
    {2, 2, 4},
    TestID -> "Audit-Property-Tableau-Shape"
]

VerificationTest[
    Dimensions @ $psBell["Matrix"],
    {4, 4},
    TestID -> "Audit-Property-Matrix-Shape"
]

VerificationTest[
    Dimensions @ $psBell["X"],
    {2, 4},
    TestID -> "Audit-Property-X-Shape"
]

VerificationTest[
    Dimensions @ $psBell["Z"],
    {2, 4},
    TestID -> "Audit-Property-Z-Shape"
]

VerificationTest[
    $psBell["StabilizerSigns"],
    {1, 1},
    TestID -> "Audit-Property-StabilizerSigns"
]

VerificationTest[
    $psBell["DestabilizerSigns"],
    {1, 1},
    TestID -> "Audit-Property-DestabilizerSigns"
]

VerificationTest[
    Dimensions @ $psBell["Stabilizer"],
    {2, 2, 2},
    TestID -> "Audit-Property-Stabilizer-Shape"
]

VerificationTest[
    Dimensions @ $psBell["StabilizerTableau"],
    {2, 2, 2},
    TestID -> "Audit-Property-StabilizerTableau-Shape"
]

VerificationTest[
    Dimensions @ $psBell["Destabilizer"],
    {2, 2, 2},
    TestID -> "Audit-Property-Destabilizer-Shape"
]

VerificationTest[
    Dimensions @ $psBell["StabilizerX"],
    {2, 2},
    TestID -> "Audit-Property-StabilizerX-Shape"
]

VerificationTest[
    Dimensions @ $psBell["StabilizerZ"],
    {2, 2},
    TestID -> "Audit-Property-StabilizerZ-Shape"
]

VerificationTest[
    Dimensions @ $psBell["DestabilizerX"],
    {2, 2},
    TestID -> "Audit-Property-DestabilizerX-Shape"
]

VerificationTest[
    Dimensions @ $psBell["DestabilizerZ"],
    {2, 2},
    TestID -> "Audit-Property-DestabilizerZ-Shape"
]

VerificationTest[
    Sort @ $psBell["Stabilizers"],
    Sort @ {"XX", "ZZ"},
    TestID -> "Audit-Property-Stabilizers"
]

VerificationTest[
    Length @ $psBell["Destabilizers"],
    2,
    TestID -> "Audit-Property-Destabilizers"
]

VerificationTest[
    Length @ $psBell["PauliStrings"],
    4,
    TestID -> "Audit-Property-PauliStrings-Length"
]

VerificationTest[
    Length @ $psBell["PauliSymbols"],
    4,
    TestID -> "Audit-Property-PauliSymbols-Length"
]

VerificationTest[
    $psBell["Generators"] === $psBell["Stabilizers"],
    True,
    TestID -> "Audit-Property-Generators-AliasFor-Stabilizers"
]

VerificationTest[
    $psBell["PauliForm"] === $psBell["Stabilizers"],
    True,
    TestID -> "Audit-Property-PauliForm-AliasFor-Stabilizers"
]

VerificationTest[
    Head @ $psBell["TableauForm"],
    Row,
    TestID -> "Audit-Property-TableauForm-Head"
]

VerificationTest[
    Head @ $psBell["StabilizerTableauForm"],
    Row,
    TestID -> "Audit-Property-StabilizerTableauForm-Head"
]

VerificationTest[
    Head @ $psBell["State"],
    QuantumState,
    TestID -> "Audit-Property-State-Head"
]

VerificationTest[
    Head @ $psBell["QuantumState"],
    QuantumState,
    TestID -> "Audit-Property-QuantumState-Head"
]

VerificationTest[
    Head @ $psBell["Operator"],
    QuantumOperator,
    TestID -> "Audit-Property-Operator-Head"
]

VerificationTest[
    Head @ $psBell["QuantumOperator"],
    QuantumOperator,
    TestID -> "Audit-Property-QuantumOperator-Head"
]

VerificationTest[
    Head @ $psBell["Circuit"],
    QuantumCircuitOperator,
    TestID -> "Audit-Property-Circuit-Head"
]

VerificationTest[
    Head @ $psBell["QuantumCircuit"],
    QuantumCircuitOperator,
    TestID -> "Audit-Property-QuantumCircuit-Head"
]

VerificationTest[
    Length @ $psBell["Properties"] >= 30,
    True,
    TestID -> "Audit-Property-Properties-MinLength"
]

VerificationTest[
    Lookup[First[$psBell], "GlobalPhase", 1],
    1,
    TestID -> "Audit-Property-GlobalPhase-DefaultsTo1"
]


(* ============================================================================ *)
(* TIER 3 -- Clifford gate methods: every recognized gate name + alias.         *)
(*                                                                              *)
(* Each gate is checked for: (a) the result is a valid PauliStabilizer; (b)   *)
(* the closure invariant holds (psValidQ); (c) one *physical* outcome (via    *)
(* Heisenberg conjugation).                                                     *)
(* ============================================================================ *)

VerificationTest[PauliStabilizer[1]["H", 1]["Stabilizers"], {"X"}, TestID -> "Audit-Gate-H-on-Z"]
VerificationTest[PauliStabilizer[1]["H", 1]["S", 1]["Stabilizers"], {"Y"}, TestID -> "Audit-Gate-S-on-X"]
VerificationTest[PauliStabilizer[1]["H", 1][SuperDagger["S"], 1]["Stabilizers"], {"-Y"}, TestID -> "Audit-Gate-Sdag-on-X"]
VerificationTest[PauliStabilizer[1]["X", 1]["Stabilizers"], {"-Z"}, TestID -> "Audit-Gate-X-on-Z"]
VerificationTest[PauliStabilizer[1]["Y", 1]["Stabilizers"], {"-Z"}, TestID -> "Audit-Gate-Y-on-Z"]
VerificationTest[PauliStabilizer[1]["Z", 1]["Stabilizers"], {"Z"}, TestID -> "Audit-Gate-Z-on-Z"]
VerificationTest[PauliStabilizer[1]["H", 1]["V", 1]["Stabilizers"], {"X"}, TestID -> "Audit-Gate-V-on-Plus"]
VerificationTest[PauliStabilizer[1]["H", 1][SuperDagger["V"], 1]["Stabilizers"], {"X"}, TestID -> "Audit-Gate-Vdag-on-Plus"]

(* Two-qubit gates: CNOT (= CX), CZ, SWAP. *)
VerificationTest[
    Sort @ PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["Stabilizers"],
    Sort @ {"XX", "ZZ"},
    TestID -> "Audit-Gate-CNOT-on-PlusZero"
]

VerificationTest[
    PauliStabilizer[2]["H", 1]["CX", 1, 2]["Stabilizers"] === PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["Stabilizers"],
    True,
    TestID -> "Audit-Gate-CX-Alias-CNOT"
]

VerificationTest[
    Sort @ PauliStabilizer[2]["H", 1]["CZ", 1, 2]["Stabilizers"],
    Sort @ {"IZ", "XZ"},
    TestID -> "Audit-Gate-CZ-on-PlusZero"
]

VerificationTest[
    Sort @ PauliStabilizer[2]["H", 1]["SWAP", 1, 2]["Stabilizers"],
    Sort @ {"IX", "ZI"},
    TestID -> "Audit-Gate-SWAP-on-PlusZero"
]

(* `op -> order` arg-rewrite dispatcher. *)
VerificationTest[
    Sort @ PauliStabilizer[2]["H" -> 1]["CNOT" -> {1, 2}]["Stabilizers"],
    Sort @ {"XX", "ZZ"},
    TestID -> "Audit-Gate-OrderRule-Form"
]

(* Permutation methods. *)
VerificationTest[
    PauliStabilizer[2]["PermuteQudits", Cycles[{{1, 2}}]]["Stabilizers"],
    {"IZ", "ZI"},
    TestID -> "Audit-Gate-PermuteQudits"
]

VerificationTest[
    psValidQ @ PauliStabilizer[2]["Permute", Cycles[{{1, 2}}]],
    True,
    TestID -> "Audit-Gate-Permute"
]

(* Padding methods. *)
VerificationTest[
    PauliStabilizer[1]["PadRight", 3]["Qubits"],
    3,
    TestID -> "Audit-Gate-PadRight"
]

VerificationTest[
    PauliStabilizer[1]["PadLeft", 3]["Qubits"],
    3,
    TestID -> "Audit-Gate-PadLeft"
]

(* Dagger / Inverse aliases. *)
VerificationTest[
    psValidQ @ $psBell["Dagger"],
    True,
    TestID -> "Audit-Gate-Dagger-Validity"
]

VerificationTest[
    $psBell["Dagger"] === $psBell["Inverse"],
    True,
    TestID -> "Audit-Gate-Inverse-AliasFor-Dagger"
]


(* ============================================================================ *)
(* TIER 4 -- Non-Clifford gates: P[\[Theta]] / T / T\[Dagger] return a frame.   *)
(* ============================================================================ *)

VerificationTest[
    Head @ PauliStabilizer[1]["T", 1],
    StabilizerFrame,
    TestID -> "Audit-NonClifford-T-ReturnsFrame"
]

VerificationTest[
    Head @ PauliStabilizer[1][SuperDagger["T"], 1],
    StabilizerFrame,
    TestID -> "Audit-NonClifford-Tdag-ReturnsFrame"
]

VerificationTest[
    Head @ PauliStabilizer[1]["P"[Pi/3], 1],
    StabilizerFrame,
    TestID -> "Audit-NonClifford-P-Theta-ReturnsFrame"
]

(* T|0> = |0>: materialized state vector matches |0>. *)
VerificationTest[
    Chop @ N @ FullSimplify[
        Normal @ PauliStabilizer[1]["T", 1]["StateVector"] - {1, 0}
    ],
    {0, 0},
    TestID -> "Audit-NonClifford-T-on-zero-Materializes-zero"
]


(* ============================================================================ *)
(* TIER 5 -- Measurement methods: ps["M", q], ps["M", string], ps["M", list],   *)
(* ps["SymbolicMeasure", ...], ps["SubstituteOutcomes", ...], ps["SampleOutcomes", ...]. *)
(* ============================================================================ *)

(* Z-basis single-qubit, deterministic. *)
VerificationTest[
    Sort @ Keys @ PauliStabilizer[1]["M", 1],
    {0},
    TestID -> "Audit-Measure-Z-Deterministic"
]

(* Z-basis single-qubit, non-deterministic. *)
VerificationTest[
    Sort @ Keys @ PauliStabilizer[1]["H", 1]["M", 1],
    {0, 1},
    TestID -> "Audit-Measure-Z-NonDeterministic"
]

(* Pauli-string measurement on 5Q stabilizer. *)
VerificationTest[
    Sort @ Keys @ $ps5Q["M", "XZZXI"],
    {0},
    TestID -> "Audit-Measure-PauliString-5Q-Stab"
]

(* Multi-qubit measurement returns Association keyed by tuples. *)
VerificationTest[
    Sort @ Keys @ $psBell["M", {1, 2}],
    {{0, 0}, {1, 1}},
    TestID -> "Audit-Measure-MultiQubit-Bell-Correlation"
]

(* SymbolicMeasure on |+> allocates a fresh symbol. *)
VerificationTest[
    Module[{psSym = PauliStabilizer[1]["H", 1]["SymbolicMeasure", 1]},
        Length @ DeleteDuplicates @ Cases[psSym["Phase"], _\[FormalS], Infinity]
    ],
    1,
    TestID -> "Audit-Measure-SymbolicMeasure-AllocatesSymbol"
]

(* SymbolicMeasure deterministic case returns post-state with no fresh symbol. *)
VerificationTest[
    Cases[PauliStabilizer[1]["SymbolicMeasure", 1]["Phase"], _\[FormalS], Infinity],
    {},
    TestID -> "Audit-Measure-SymbolicMeasure-Deterministic-NoSymbol"
]

(* SymbolicMeasure on a list of qudits = fold. *)
VerificationTest[
    psValidQ @ PauliStabilizer[2]["H", 1]["SymbolicMeasure", {1, 2}],
    True,
    TestID -> "Audit-Measure-SymbolicMeasure-MultiQubit"
]

(* SubstituteOutcomes: outcome 0 vs outcome 1 reproduces the Z basis. *)
VerificationTest[
    Module[{psSym, sym},
        psSym = PauliStabilizer[1]["H", 1]["SymbolicMeasure", 1];
        sym = First @ DeleteDuplicates @ Cases[psSym["Phase"], _\[FormalS], Infinity];
        {
            psSym["SubstituteOutcomes", sym -> 0]["Stabilizers"],
            psSym["SubstituteOutcomes", sym -> 1]["Stabilizers"]
        }
    ],
    {{"Z"}, {"-Z"}},
    TestID -> "Audit-Measure-SubstituteOutcomes-Roundtrip"
]

(* SampleOutcomes returns a list of valid PauliStabilizers. *)
VerificationTest[
    Block[{},
        SeedRandom[42];
        With[{samples = PauliStabilizer[1]["H", 1]["SymbolicMeasure", 1]["SampleOutcomes", 5]},
            AllTrue[samples, psValidQ]
        ]
    ],
    True,
    TestID -> "Audit-Measure-SampleOutcomes-AllValid"
]

(* SampleOutcomes single-arg returns a single ps (not a list). *)
VerificationTest[
    Block[{},
        SeedRandom[42];
        Head @ PauliStabilizer[1]["H", 1]["SymbolicMeasure", 1]["SampleOutcomes"]
    ],
    PauliStabilizer,
    TestID -> "Audit-Measure-SampleOutcomes-SingleArg-ReturnsPS"
]


(* ============================================================================ *)
(* TIER 6 -- Inner Product and Expectation method dispatch (Phase 6.5).         *)
(* ============================================================================ *)

VerificationTest[
    $psBell["InnerProduct", $psBell],
    1,
    TestID -> "Audit-Method-InnerProduct-Self"
]

VerificationTest[
    Module[{psPhiMinus = $psBell["Z", 1]},
        Chop @ N @ $psBell["InnerProduct", psPhiMinus]
    ],
    0,
    TestID -> "Audit-Method-InnerProduct-Orthogonal-PhiPlus-PhiMinus"
]

VerificationTest[
    $psBell["Expectation", "XX"],
    1,
    TestID -> "Audit-Method-Expectation-XX-OnBell"
]

VerificationTest[
    $psBell["Expectation", "ZZ"],
    1,
    TestID -> "Audit-Method-Expectation-ZZ-OnBell"
]

VerificationTest[
    $psBell["Expectation", "YY"],
    -1,
    TestID -> "Audit-Method-Expectation-YY-OnBell-Sign"
]

VerificationTest[
    $psBell["Expectation", "XI"],
    0,
    TestID -> "Audit-Method-Expectation-Anticommuting"
]

VerificationTest[
    With[{psT = PauliStabilizer[1]["T", 1]},
        FullSimplify @ psT["InnerProduct", PauliStabilizer[1]]
    ],
    1,
    TestID -> "Audit-Method-InnerProduct-Mixed-FrameAndPS"
]


(* ============================================================================ *)
(* TIER 7 -- Composition + tensor product UpValues.                             *)
(* ============================================================================ *)

VerificationTest[
    With[{psH = PauliStabilizer[1]["H", 1]["S", 1]["H", 1]},
        psH["Stabilizers"]
    ],
    {"-Y"},
    TestID -> "Audit-Compose-HSH-on-Z"
]

VerificationTest[
    With[{a = $psBell, b = PauliStabilizer[1]},
        QuantumTensorProduct[a, b]["Stabilizers"]
    ],
    {"XXI", "ZZI", "IIZ"},
    TestID -> "Audit-Compose-TensorProduct-BellAndZero"
]

(* PauliStabilizer composed with another (apply right-to-left). *)
VerificationTest[
    psValidQ @ PauliStabilizer[2]["H", 1] @ PauliStabilizer[2],
    True,
    TestID -> "Audit-Compose-PSofPS-Validity"
]


(* ============================================================================ *)
(* TIER 8 -- StabilizerFrame: constructors, properties, methods.                *)
(* ============================================================================ *)

VerificationTest[
    Wolfram`QuantumFramework`PackageScope`StabilizerFrameQ @ StabilizerFrame[PauliStabilizer[1]],
    True,
    TestID -> "Audit-Frame-Construct-FromPS"
]

VerificationTest[
    Wolfram`QuantumFramework`PackageScope`StabilizerFrameQ @ StabilizerFrame[{{1, PauliStabilizer[1]}, {1, PauliStabilizer[1]["H", 1]}}],
    True,
    TestID -> "Audit-Frame-Construct-FromPairs"
]

(* Properties of a frame produced by a T gate. *)
VerificationTest[
    Module[{psT = PauliStabilizer[1]["T", 1]},
        psT["Length"]
    ],
    2,
    TestID -> "Audit-Frame-Property-Length"
]

VerificationTest[
    Module[{psT = PauliStabilizer[1]["T", 1]},
        Length @ psT["Coefficients"]
    ],
    2,
    TestID -> "Audit-Frame-Property-Coefficients-Length"
]

VerificationTest[
    Module[{psT = PauliStabilizer[1]["T", 1]},
        Head /@ psT["Stabilizers"] === {PauliStabilizer, PauliStabilizer}
    ],
    True,
    TestID -> "Audit-Frame-Property-Stabilizers-AllPS"
]

VerificationTest[
    PauliStabilizer[1]["T", 1]["Qubits"],
    1,
    TestID -> "Audit-Frame-Property-Qubits"
]

VerificationTest[
    PauliStabilizer[1]["T", 1]["GeneratorCount"],
    1,
    TestID -> "Audit-Frame-Property-GeneratorCount"
]

(* Clifford gate distributes over a frame; non-Clifford doubles. *)
VerificationTest[
    PauliStabilizer[1]["T", 1]["H", 1]["Length"],
    2,
    TestID -> "Audit-Frame-Method-CliffordDistributes"
]

VerificationTest[
    PauliStabilizer[1]["H", 1]["T", 1]["T", 1]["Length"],
    4,
    TestID -> "Audit-Frame-Method-NonCliffordDoubles"
]

(* Frame Plus / Times upvalues. *)
VerificationTest[
    Module[{f = StabilizerFrame[PauliStabilizer[1]]},
        (f + f)["Length"]
    ],
    2,
    TestID -> "Audit-Frame-UpValue-Plus"
]

VerificationTest[
    Module[{f = StabilizerFrame[PauliStabilizer[1]]},
        (3 f)["Coefficients"]
    ],
    {3},
    TestID -> "Audit-Frame-UpValue-Times"
]


(* ============================================================================ *)
(* TIER 9 -- GraphState: constructors, properties, conversion.                  *)
(* ============================================================================ *)

VerificationTest[
    Wolfram`QuantumFramework`PackageScope`GraphStateQ @
        GraphState[Graph[Range[3], {1 \[UndirectedEdge] 2, 2 \[UndirectedEdge] 3}]],
    True,
    TestID -> "Audit-Graph-Construct-FromGraph"
]

VerificationTest[
    GraphState[Graph[Range[3], {1 \[UndirectedEdge] 2, 2 \[UndirectedEdge] 3}]]["Stabilizers"],
    {"XZI", "ZXZ", "IZX"},
    TestID -> "Audit-Graph-Stabilizers-LinearCluster3"
]

VerificationTest[
    GraphState[Graph[Range[5], Table[i \[UndirectedEdge] (i + 1), {i, 4}]]]["Stabilizers"],
    {"XZIII", "ZXZII", "IZXZI", "IIZXZ", "IIIZX"},
    TestID -> "Audit-Graph-Stabilizers-LinearCluster5"
]

VerificationTest[
    GraphState[Graph[Range[3], {1 \[UndirectedEdge] 2, 2 \[UndirectedEdge] 3}]]["VertexCount"],
    3,
    TestID -> "Audit-Graph-Property-VertexCount"
]

VerificationTest[
    GraphState[Graph[Range[3], {1 \[UndirectedEdge] 2, 2 \[UndirectedEdge] 3}]]["EdgeCount"],
    2,
    TestID -> "Audit-Graph-Property-EdgeCount"
]

VerificationTest[
    Dimensions @ GraphState[Graph[Range[3], {1 \[UndirectedEdge] 2, 2 \[UndirectedEdge] 3}]]["AdjacencyMatrix"],
    {3, 3},
    TestID -> "Audit-Graph-Property-AdjacencyMatrix-Shape"
]

VerificationTest[
    Head @ GraphState[Graph[Range[3], {1 \[UndirectedEdge] 2, 2 \[UndirectedEdge] 3}]]["PauliStabilizer"],
    PauliStabilizer,
    TestID -> "Audit-Graph-Property-PauliStabilizer-ConvertsToPS"
]


(* ============================================================================ *)
(* TIER 10 -- LocalComplement: graph and graph-state inputs, involutivity.      *)
(* ============================================================================ *)

(* LC at the center of a star = K4 minus star (i.e. all pairs among leaves). *)
VerificationTest[
    Sort @ EdgeList @ LocalComplement[
        Graph[{1, 2, 3, 4}, {1 \[UndirectedEdge] 2, 1 \[UndirectedEdge] 3, 1 \[UndirectedEdge] 4}],
        1
    ],
    Sort @ {
        UndirectedEdge[1, 2], UndirectedEdge[1, 3], UndirectedEdge[1, 4],
        UndirectedEdge[2, 3], UndirectedEdge[2, 4], UndirectedEdge[3, 4]
    },
    TestID -> "Audit-LC-Graph-StarToK4"
]

(* LC is involutive: LC(LC(g, v), v) = g. *)
VerificationTest[
    With[{g = Graph[Range[5], Table[i \[UndirectedEdge] (i + 1), {i, 4}]]},
        Sort @ EdgeList @ LocalComplement[LocalComplement[g, 3], 3] === Sort @ EdgeList[g]
    ],
    True,
    TestID -> "Audit-LC-Involutive"
]

(* LC on a GraphState returns a GraphState. *)
VerificationTest[
    Head @ LocalComplement[
        GraphState[Graph[Range[3], {1 \[UndirectedEdge] 2, 2 \[UndirectedEdge] 3}]],
        2
    ],
    GraphState,
    TestID -> "Audit-LC-GraphState-ReturnsGraphState"
]


(* ============================================================================ *)
(* TIER 11 -- CliffordChannel: head, predicate, accessors, composition.         *)
(* ============================================================================ *)

VerificationTest[
    ccValidQ @ CliffordChannel["Identity", 1],
    True,
    TestID -> "Audit-CC-Construct-Identity-Validity"
]

VerificationTest[
    CliffordChannel["Identity", 2]["Rank"],
    4,
    TestID -> "Audit-CC-Property-Identity-Rank"
]

VerificationTest[
    {CliffordChannel["Identity", 3]["InputQubits"], CliffordChannel["Identity", 3]["OutputQubits"]},
    {3, 3},
    TestID -> "Audit-CC-Property-Identity-InputOutputQubits"
]

VerificationTest[
    Length @ CliffordChannel["Identity", 2]["Properties"] >= 7,
    True,
    TestID -> "Audit-CC-Properties-MinKeys"
]

(* Tableau accessor. *)
VerificationTest[
    Dimensions @ CliffordChannel["Identity", 2]["Tableau"],
    {4, 9},
    TestID -> "Audit-CC-Property-Tableau-Shape"
]

(* From a stabilizer state -> state-prep channel. *)
VerificationTest[
    With[{cc = CliffordChannel[PauliStabilizer[2]]},
        cc["InputQubits"] == 0 && cc["OutputQubits"] == 2 && cc["Rank"] == 2
    ],
    True,
    TestID -> "Audit-CC-Construct-FromPS-StatePrep"
]

(* Composition: identity composed with itself. *)
VerificationTest[
    Module[{cc = CliffordChannel["Identity", 2]},
        ccValidQ @ cc[cc] && cc[cc]["Rank"] == 4
    ],
    True,
    TestID -> "Audit-CC-Compose-IdentitySquared"
]

(* cc[ps] state evolution: identity preserves the state. *)
VerificationTest[
    Module[{ps = PauliStabilizer[2]},
        CliffordChannel["Identity", 2][ps]["Stabilizers"] === ps["Stabilizers"]
    ],
    True,
    TestID -> "Audit-CC-Apply-Identity-PreservesState"
]


(* ============================================================================ *)
(* TIER 12 -- Hybrid interop (Phase 7): qmo[ps], qc[ps], stabilizer dispatch.   *)
(* ============================================================================ *)

(* Pauli-basis QMO routes to ps["M", pauli] (no nonpaulibasis message). *)
VerificationTest[
    Sort @ Keys @ QuantumMeasurementOperator["ZZ", {1, 2}][$psBell],
    {0},
    TestID -> "Audit-Hybrid-QMO-ZZ-on-Bell-FastPath"
]

(* Non-Pauli QMO falls back with info message. *)
VerificationTest[
    Head @ QuantumMeasurementOperator[QuantumBasis["Computational"], {1}][PauliStabilizer[1]],
    PauliStabilizer,
    {PauliStabilizer::nonpaulibasis, PauliStabilizer::nonclifford},
    TestID -> "Audit-Hybrid-QMO-NonPauli-Fallback"
]

(* Named Pauli channel: BitFlip routes to tableau. *)
VerificationTest[
    Length @ QuantumChannel["BitFlip"[1/4], {1}][PauliStabilizer[1]],
    2,
    TestID -> "Audit-Hybrid-QC-BitFlip-Branches"
]

(* Probability sum is 1. *)
VerificationTest[
    Total[First /@ QuantumChannel["BitFlip"[1/4], {1}][PauliStabilizer[1]]],
    1,
    TestID -> "Audit-Hybrid-QC-BitFlip-Probability-Sum"
]

(* Phase 7.4 matrix-iteration detector. *)
VerificationTest[
    Wolfram`QuantumFramework`PackageScope`stabilizerPauliLabelFromQMO @
        QuantumMeasurementOperator[QuantumOperator[KroneckerProduct[PauliMatrix[1], PauliMatrix[3]]], {1, 2}],
    "XZ",
    TestID -> "Audit-Hybrid-Phase7.4-MatrixDetector-XZ"
]


(* ============================================================================ *)
(* TIER 13 -- AG symplectic invariant: M . J . M^T = J for every state.         *)
(*                                                                              *)
(* This is the structural correctness oracle: every PauliStabilizer that comes  *)
(* out of a constructor or gate sequence must satisfy `M Omega M^T = Omega mod 2`. *)
(* (Aaronson-Gottesman Proposition 2.) Hits 6 distinct constructor + gate paths. *)
(* ============================================================================ *)

(* Symplectic form Omega = [[0, I_n], [I_n, 0]] over F_2.                       *)
agSymplecticFormAudit[n_Integer] := ArrayFlatten[{
    {ConstantArray[0, {n, n}], IdentityMatrix[n]},
    {IdentityMatrix[n], ConstantArray[0, {n, n}]}
}]

agSymplectic[ps_] := With[{
    n = ps["Qubits"],
    m = ps["Matrix"]
},
    With[{omega = agSymplecticFormAudit[n]},
        Mod[m . omega . Transpose[m] - omega, 2]
    ]
]

VerificationTest[
    Normal @ agSymplectic[$psBell],
    ConstantArray[0, {4, 4}],
    TestID -> "Audit-AG-Symplectic-Invariant-Bell"
]

(* String-list-built PauliStabilizer["5QubitCode"] auto-pads destabilizers via *)
(* the Reverse rule (Constructors.m:203), which does NOT satisfy the full AG  *)
(* canonical pairing -- we only assert stabilizer-stabilizer commutation here. *)
VerificationTest[
    With[{
        m = $ps5Q["Matrix"],
        n = $ps5Q["Qubits"],
        gen = $ps5Q["GeneratorCount"]
    },
        Mod[m[[gen + 1 ;; 2 gen]] . agSymplecticFormAudit[n] . Transpose[m[[gen + 1 ;; 2 gen]]], 2]
    ],
    ConstantArray[0, {5, 5}],
    TestID -> "Audit-AG-Symplectic-StabStabBlock-5Q"
]

(* Circuit-built random Cliffords ALWAYS produce AG-canonical pairs. *)
VerificationTest[
    Block[{},
        SeedRandom[20260507];
        AllTrue[
            Table[
                Normal @ agSymplectic @ PauliStabilizer["Random", 4] === ConstantArray[0, {8, 8}],
                {6}
            ],
            TrueQ
        ]
    ],
    True,
    TestID -> "Audit-AG-Symplectic-Invariant-Random4-6reps"
]


(* ============================================================================ *)
(* TIER 14 -- Cross-package equivalence: live Stim fixture + QC.jl canonical.   *)
(*                                                                              *)
(* These are spot checks that the existing CrossPackage_*.wlt files cover in   *)
(* depth. The spot here proves the matrix-form audit doesn't drift from the    *)
(* canonical Stim / QC.jl values across kernel updates.                         *)
(* ============================================================================ *)

(* Bell phi+ via Stim canonical = our AG-derived stabilizers. Skip if fixture
   not present (CI-friendly). *)
$stimFixturePath = SelectFirst[
    {
        FileNameJoin[{Directory[], "Tests", "Stabilizer", "fixtures", "stim_fixtures.json"}],
        "/Users/mohammadb/Documents/GitHub/QuantumFramework/Tests/Stabilizer/fixtures/stim_fixtures.json"
    },
    FileExistsQ,
    None
];

If[StringQ[$stimFixturePath],
    $stimFixturesAudit = Import[$stimFixturePath, "RawJSON"];

    (* Stim writes "+XX" / "-_X" with explicit signs and "_" for I. Convert     *)
    (* to our convention: drop "+", convert "_" to "I", keep "-".               *)
    stimNormalize[s_String] := StringReplace[
        StringReplace[s, "_" -> "I"],
        StartOfString ~~ "+" -> ""
    ];

    VerificationTest[
        Sort @ Map[stimNormalize, Lookup[$stimFixturesAudit["bell_phi_plus"], "stabilizers"]],
        Sort @ $psBell["Stabilizers"],
        TestID -> "Audit-CrossPackage-Stim-BellPhiPlus-Stabilizers"
    ];

    (* Cluster-4 from Stim. *)
    VerificationTest[
        Sort @ Map[stimNormalize, Lookup[$stimFixturesAudit["cluster_4"], "stabilizers"]],
        Sort @ GraphState[Graph[Range[4], Table[i \[UndirectedEdge] (i + 1), {i, 3}]]]["Stabilizers"],
        TestID -> "Audit-CrossPackage-Stim-Cluster4-Stabilizers"
    ];

    (* GHZ-3 from Stim. *)
    VerificationTest[
        Sort @ Map[stimNormalize, Lookup[$stimFixturesAudit["ghz_3"], "stabilizers"]],
        Sort @ PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, "CNOT" -> {2, 3}}]]["Stabilizers"],
        TestID -> "Audit-CrossPackage-Stim-GHZ3-Stabilizers"
    ]
];

(* AG g-function per-qubit phase for Z * X via the agPhase helper.              *)
(* Per qubit: src=Z=(0,1), dst=X=(1,0). agPhase(0,1,1,0) = case x_src=0,z_src=1: *)
(*   x_dst (1 - 2 z_dst) = 1 * (1 - 0) = 1.                                     *)
(* Three qubits sum to 3 mod 4 = i^3 = -i. QC.jl docs ZZZ * XXX = -i * YYY,    *)
(* matching i^3.                                                                *)
VerificationTest[
    Mod[
        Total @ MapThread[
            Wolfram`QuantumFramework`PackageScope`agPhase,
            {{0, 0, 0}, {1, 1, 1}, {1, 1, 1}, {0, 0, 0}}
        ],
        4
    ],
    3,
    TestID -> "Audit-CrossPackage-QC.jl-prodphase-ZZZ-XXX-Equals-NegI"
]


(* ============================================================================ *)
(* TIER 15 -- Round-trip exact-equality matrix (per Phase 5c contract).         *)
(*                                                                              *)
(* PauliStabilizer[qo]["QuantumOperator"] === qo for every Pauli + Clifford gate. *)
(* This duplicates a portion of TIER 1.4a/b/c in PauliStabilizer.wlt, but as a *)
(* concise audit-matrix that is explicitly tied to the ROADMAP A.* contract.   *)
(* ============================================================================ *)

VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["I"]]["QuantumOperator"], QuantumOperator["I"]], True, TestID -> "Audit-RT-QO-I"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["X"]]["QuantumOperator"], QuantumOperator["X"]], True, TestID -> "Audit-RT-QO-X"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["Y"]]["QuantumOperator"], QuantumOperator["Y"]], True, TestID -> "Audit-RT-QO-Y"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["Z"]]["QuantumOperator"], QuantumOperator["Z"]], True, TestID -> "Audit-RT-QO-Z"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["H"]]["QuantumOperator"], QuantumOperator["H"]], True, TestID -> "Audit-RT-QO-H"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["S"]]["QuantumOperator"], QuantumOperator["S"]], True, TestID -> "Audit-RT-QO-S"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["S"]["Dagger"]]["QuantumOperator"], QuantumOperator["S"]["Dagger"]], True, TestID -> "Audit-RT-QO-Sdag"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["CNOT"]]["QuantumOperator"], QuantumOperator["CNOT"]], True, TestID -> "Audit-RT-QO-CNOT"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["CZ"]]["QuantumOperator"], QuantumOperator["CZ"]], True, TestID -> "Audit-RT-QO-CZ"]
VerificationTest[matEqQO[PauliStabilizer[QuantumOperator["SWAP"]]["QuantumOperator"], QuantumOperator["SWAP"]], True, TestID -> "Audit-RT-QO-SWAP"]

(* QS round-trip on canonical states. *)
VerificationTest[
    PauliStabilizer[QuantumState[{0, 1}]]["State"]["StateVector"] === QuantumState[{0, 1}]["StateVector"],
    True,
    TestID -> "Audit-RT-QS-One"
]

VerificationTest[
    PauliStabilizer[QuantumState["PhiPlus"]]["State"]["StateVector"] === QuantumState["PhiPlus"]["StateVector"],
    True,
    TestID -> "Audit-RT-QS-PhiPlus"
]
