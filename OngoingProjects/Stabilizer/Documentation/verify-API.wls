#!/usr/bin/env wolframscript

(* Verifier for OngoingProjects/Stabilizer/Documentation/API.md.
   Each block-id heading marks a code example referenced in the API doc; the
   verifier re-runs every block and prints labeled output. The MD's output
   blocks should match this script's output verbatim.

   Run from repo root:
     wolframscript -f OngoingProjects/Stabilizer/Documentation/verify-API.wls
*)

(* Resolve paclet directory from this script's location: ../../QuantumFramework
   relative to OngoingProjects/Stabilizer/Documentation/. Works for both the
   canonical repo and any worktree. *)
$pacletDir = FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "..", "QuantumFramework"}];
Quiet[PacletUninstall /@ PacletFind["Wolfram/QuantumFramework"]];
PacletDirectoryLoad[$pacletDir];
Needs["Wolfram`QuantumFramework`"];

block[id_, expr_] := Module[{result},
    result = expr;
    Print["=== ", id, " ==="];
    Print[ToString[result, InputForm]];
    Print[]
];


(* ============================================================================ *)
(* PauliStabilizer constructors                                                 *)
(* ============================================================================ *)

block["PS-ctor-string-list",
    PauliStabilizer[{"XX", "ZZ"}]["Stabilizers"]
];

block["PS-ctor-string-list-signed",
    PauliStabilizer[{"-XX", "+ZZ"}]["StabilizerSigns"]
];

block["PS-ctor-stab-destab",
    With[{ps = PauliStabilizer[{"XX", "ZZ"}, {"IX", "ZI"}]},
        {ps["Stabilizers"], ps["Destabilizers"]}
    ]
];

block["PS-ctor-integer",
    With[{ps = PauliStabilizer[3]},
        {ps["Qubits"], ps["Stabilizers"], ps["Destabilizers"]}
    ]
];

block["PS-ctor-named-5Q",
    With[{ps = PauliStabilizer["5QubitCode"]},
        {ps["Qubits"], ps["GeneratorCount"], ps["Stabilizers"]}
    ]
];

block["PS-ctor-named-Steane",
    With[{ps = PauliStabilizer["SteaneCode"]},
        {ps["Qubits"], ps["GeneratorCount"], ps["Stabilizers"][[7]]}
    ]
];

block["PS-ctor-from-circuit",
    PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]["Stabilizers"]
];

block["PS-ctor-from-operator",
    PauliStabilizer[QuantumOperator["H", 1]]["Stabilizers"]
];

block["PS-ctor-random",
    Block[{},
        SeedRandom[20260430];
        With[{ps = PauliStabilizer["Random", 3]},
            {ps["Qubits"], ps["Stabilizers"]}
        ]
    ]
];


(* ============================================================================ *)
(* PauliStabilizer properties                                                   *)
(* ============================================================================ *)

block["PS-prop-shape",
    With[{ps = PauliStabilizer["5QubitCode"]},
        <|
            "Qubits"           -> ps["Qubits"],
            "GeneratorCount"   -> ps["GeneratorCount"],
            "Tableau-shape"    -> Dimensions[ps["Tableau"]],
            "Matrix-shape"     -> Dimensions[ps["Matrix"]],
            "Signs-length"     -> Length[ps["Signs"]]
        |>
    ]
];

block["PS-prop-stabilizer-vs-destabilizer",
    With[{ps = PauliStabilizer["5QubitCode"]},
        <|
            "Stabilizers"        -> ps["Stabilizers"],
            "Destabilizers"      -> ps["Destabilizers"],
            "StabilizerSigns"    -> ps["StabilizerSigns"],
            "DestabilizerSigns"  -> ps["DestabilizerSigns"]
        |>
    ]
];

block["PS-prop-XZ-bits",
    With[{ps = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]},
        <|
            "X" -> ps["X"],
            "Z" -> ps["Z"],
            "Phase" -> ps["Phase"]
        |>
    ]
];

block["PS-prop-pauli-strings",
    With[{ps = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]},
        ps["PauliStrings"]
    ]
];


(* ============================================================================ *)
(* PauliStabilizer Clifford gate methods                                        *)
(* ============================================================================ *)

block["PS-gate-H",
    PauliStabilizer[1]["H", 1]["Stabilizers"]
];

block["PS-gate-S",
    (* H|0> has stab X; apply S -> stab Y *)
    PauliStabilizer[1]["H", 1]["S", 1]["Stabilizers"]
];

block["PS-gate-CNOT",
    Sort @ PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["Stabilizers"]
];

block["PS-gate-CZ",
    Sort @ PauliStabilizer[2]["H", 1]["CZ", 1, 2]["Stabilizers"]
];

block["PS-gate-SWAP",
    Sort @ PauliStabilizer[2]["H", 1]["SWAP", 1, 2]["Stabilizers"]
];

block["PS-gate-X",
    (* X gate flips phase row where X-bit = 1 *)
    With[{ps5 = PauliStabilizer["5QubitCode"]},
        ps5["X", 1]["StabilizerSigns"]
    ]
];

block["PS-gate-Z",
    With[{ps5 = PauliStabilizer["5QubitCode"]},
        ps5["Z", 3]["StabilizerSigns"]
    ]
];

block["PS-gate-Sdag-V-Vdag",
    With[{ps = PauliStabilizer[1]["H", 1]},
        <|
            "Sdag-on-X" -> ps[SuperDagger["S"], 1]["Stabilizers"],
            "V-on-X"    -> ps["V", 1]["Stabilizers"],
            "Vdag-on-X" -> ps[SuperDagger["V"], 1]["Stabilizers"]
        |>
    ]
];

block["PS-gate-arrow-form",
    (* op -> order convenience syntax *)
    With[{ps = PauliStabilizer[2]["H" -> 1]["CNOT" -> {1, 2}]},
        Sort @ ps["Stabilizers"]
    ]
];

block["PS-gate-T-returns-frame",
    (* T is non-Clifford; result is StabilizerFrame *)
    Head @ PauliStabilizer[1]["T", 1]
];


(* ============================================================================ *)
(* PauliStabilizer measurement                                                  *)
(* ============================================================================ *)

block["PS-measure-Z-deterministic",
    (* |0> measured in Z is deterministic, outcome 0 *)
    PauliStabilizer[1]["M", 1]
];

block["PS-measure-Z-random",
    (* H|0> = |+> measured in Z is random, both outcomes *)
    Sort @ Keys @ PauliStabilizer[1]["H", 1]["M", 1]
];

block["PS-measure-Pauli-string",
    (* 5Q code state stabilizer measurement gives 0 deterministically *)
    With[{ps5 = PauliStabilizer["5QubitCode"]},
        Keys @ ps5["M", "XZZXI"]
    ]
];

block["PS-measure-multi-qubit",
    (* Multi-qubit measurement returns Association keyed by tuples *)
    With[{psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]},
        Sort @ Keys @ psBell["M", {1, 2}]
    ]
];


(* ============================================================================ *)
(* PauliStabilizer composition + tensor                                         *)
(* ============================================================================ *)

block["PS-compose",
    (* Apply ps2 after ps1: ps1[ps2] (note: this is composition of Cliffords) *)
    Module[{psH = PauliStabilizer[1]["H", 1]["S", 1]["H", 1], expected},
        (* HSH on |0>: H -> |+>, S -> |+i>, H -> ... *)
        psH["Stabilizers"]
    ]
];

block["PS-tensor-product",
    With[{a = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]],
          b = PauliStabilizer[1]},
        QuantumTensorProduct[a, b]["Stabilizers"]
    ]
];


(* ============================================================================ *)
(* PauliStabilizer conversions                                                  *)
(* ============================================================================ *)

block["PS-conv-State",
    (* Materialize Bell state *)
    With[{psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]},
        Normal @ psBell["State"]["StateVector"]
    ]
];

block["PS-conv-Circuit-equivalence",
    (* ps['Circuit'] returns a Clifford circuit whose dagger equals ps.
       Verification: applying the circuit to |0...0> reproduces the original
       state vector (up to global phase). String equality of stabilizer
       generators does NOT hold (different but equivalent generating sets). *)
    Module[{psBell, circ, fromCircuit, vec1, vec2, ratio},
        psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
        circ = psBell["Circuit"];
        fromCircuit = PauliStabilizer[circ];
        vec1 = N @ Normal @ psBell["State"]["StateVector"];
        vec2 = N @ Normal @ fromCircuit["State"]["StateVector"];
        (* Up to global phase: the inner product magnitude should be 1 *)
        Chop @ Abs[Conjugate[vec1] . vec2] - 1
    ]
];


(* ============================================================================ *)
(* Round-trip contract (Phase 5c) -- block order matches API.md doc order       *)
(* ============================================================================ *)

block["Roundtrip-QO-Y-GlobalPhase",
    (* First-hop exact: Y is the canary (AG decomposition recovers Z@X = i*Y;
       constructor captures the missing -I as "GlobalPhase"). *)
    With[{ps = PauliStabilizer[QuantumOperator["Y"]]},
        <|"GlobalPhase" -> ps["GlobalPhase"], "matches" -> ps["QuantumOperator"]["Matrix"] === QuantumOperator["Y"]["Matrix"]|>
    ]
];

block["Roundtrip-QS-One-GlobalPhase",
    (* First-hop exact: |1> stabilizer-tomography returns the right state. *)
    With[{ps = PauliStabilizer[QuantumState[{0, 1}]]},
        <|"GlobalPhase" -> ps["GlobalPhase"], "matches" -> ps["State"]["StateVector"] === QuantumState[{0, 1}]["StateVector"]|>
    ]
];

block["Roundtrip-GateUpdate-UpToPhase-YOn1",
    (* Gate-update path is up-to-phase only (ROADMAP A.9): GP is dropped. *)
    With[{
        actual   = PauliStabilizer[QuantumState[{0, 1}]]["Y", 1]["State"]["StateVector"] // Normal,
        expected = (QuantumOperator["Y"][QuantumState[{0, 1}]])["StateVector"] // Normal
    },
        <|"actual" -> actual, "expected" -> expected, "abs-overlap" -> Quiet @ Simplify[Abs[Conjugate[actual] . expected]]|>
    ]
];

block["Roundtrip-EscapeHatch-Y-on-1",
    (* Escape hatch: re-construct from post-gate state to recover GP. *)
    With[{qs = QuantumState[{0, 1}], gate = QuantumOperator["Y"]},
        PauliStabilizer[gate[qs]]["State"]["StateVector"] // Normal // Simplify
    ]
];


(* ============================================================================ *)
(* Integration: Method -> "Stabilizer" (matches API.md doc order)               *)
(* ============================================================================ *)

block["Integration-Method-Stabilizer",
    QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}][Method -> "Stabilizer"]["Stabilizers"]
];

block["Integration-named-circuits",
    QuantumCircuitOperator["GHZ"[3]][Method -> "Stabilizer"]["Stabilizers"]
];

block["Integration-UpValues",
    (* QuantumOperator['H', 1] @ ps applies via PauliStabilizerApply *)
    With[{ps = PauliStabilizer[1]},
        (QuantumOperator["H", 1] @ ps)["Stabilizers"]
    ]
];


(* ============================================================================ *)
(* PauliStabilizer["Random", n]  (Phase 6: was RandomClifford[n])               *)
(* ============================================================================ *)

block["RandomClifford-basic",
    Block[{},
        SeedRandom[42];
        With[{r = PauliStabilizer["Random", 3]},
            <|"Qubits" -> r["Qubits"], "Stabilizers" -> r["Stabilizers"]|>
        ]
    ]
];

block["RandomClifford-uniformity-n1",
    Block[{},
        SeedRandom[42];
        Length @ DeleteDuplicates @ Table[
            With[{r = PauliStabilizer["Random", 1]}, {r["Stabilizers"], r["Signs"]}],
            {200}
        ]
    ]
];


(* ============================================================================ *)
(* Symbolic-measurement methods on PauliStabilizer  (Phase 6: were              *)
(*   StabilizerMeasure / SubstituteOutcomes / SampleOutcomes top-level symbols).*)
(* ============================================================================ *)

block["StabilizerMeasure-deterministic",
    (* |0> Z-measurement is deterministic, returns post-state directly *)
    PauliStabilizer[1]["SymbolicMeasure", 1]["Stabilizers"]
];

block["StabilizerMeasure-symbolic",
    (* H|0> Z-measurement is non-deterministic, allocates a fresh symbol *)
    Module[{psSym = PauliStabilizer[1]["H", 1]["SymbolicMeasure", 1]},
        <|
            "Head"          -> Head[psSym],
            "Phase"         -> psSym["Phase"],
            "FreshSymbols"  -> DeleteDuplicates @ Cases[psSym["Phase"], _\[FormalS], Infinity]
        |>
    ]
];

block["SubstituteOutcomes-roundtrip",
    Module[{psSym, sym, ps0, ps1Reg},
        psSym = PauliStabilizer[1]["H", 1]["SymbolicMeasure", 1];
        sym = First @ DeleteDuplicates @ Cases[psSym["Phase"], _\[FormalS], Infinity];
        ps0 = psSym["SubstituteOutcomes", sym -> 0];
        ps1Reg = psSym["SubstituteOutcomes", sym -> 1];
        <|
            "outcome 0 stabilizers" -> ps0["Stabilizers"],
            "outcome 1 stabilizers" -> ps1Reg["Stabilizers"]
        |>
    ]
];

block["SampleOutcomes-multiple",
    Block[{},
        SeedRandom[42];
        With[{psSym = PauliStabilizer[1]["H", 1]["SymbolicMeasure", 1]},
            #["Stabilizers"] & /@ psSym["SampleOutcomes", 10]
        ]
    ]
];


(* ============================================================================ *)
(* StabilizerFrame -- block order matches API.md doc order:                     *)
(*   Constructors / FromList   -> Clifford-distributes (line 668 in MD)         *)
(*   Constructors / FromPS     -> from-PS              (line 681)               *)
(*   Constructors / FromT      -> from-T               (line 694)               *)
(*   Methods / Distribute      -> Clifford-distributes (line 725, dup of 668)   *)
(*   Methods / TT              -> TT                   (line 741)               *)
(*   Methods / Materialization -> T-on-zero            (line 757)               *)
(* ============================================================================ *)

block["StabilizerFrame-Constructor-FromList",
    (* MD section: Constructors / From a list of {coeff, PS} pairs.              *)
    (* Same code as StabilizerFrame-Method-CliffordDistributes below; in MD this  *)
    (* example serves both as a constructor demo and as a Clifford-distribution  *)
    (* demo, so the verifier exercises it twice with different IDs.              *)
    With[{f = StabilizerFrame[{{1, PauliStabilizer[1]}, {1, PauliStabilizer[1]["H", 1]}}]},
        With[{fH = f["H", 1]},
            <|"Length" -> fH["Length"], "Stabilizers-each" -> #["Stabilizers"] & /@ fH["Stabilizers"]|>
        ]
    ]
];

block["StabilizerFrame-from-PS",
    (* Build a frame from a single PauliStabilizer (coefficient = 1) *)
    With[{f = StabilizerFrame[PauliStabilizer[1]]},
        <|"Length" -> f["Length"], "Components" -> f["Components"]|>
    ]
];

block["StabilizerFrame-from-T",
    (* T gate produces 2-component frame *)
    With[{psT = PauliStabilizer[1]["T", 1]},
        <|
            "Head"         -> Head[psT],
            "Length"       -> psT["Length"],
            "Coefficients" -> psT["Coefficients"]
        |>
    ]
];

block["StabilizerFrame-Method-CliffordDistributes",
    (* MD section: Methods / Clifford gates distribute over components.          *)
    (* Identical code to StabilizerFrame-Constructor-FromList above.             *)
    With[{f = StabilizerFrame[{{1, PauliStabilizer[1]}, {1, PauliStabilizer[1]["H", 1]}}]},
        With[{fH = f["H", 1]},
            <|"Length" -> fH["Length"], "Stabilizers-each" -> #["Stabilizers"] & /@ fH["Stabilizers"]|>
        ]
    ]
];

block["StabilizerFrame-TT",
    (* T^2 produces a 4-component frame *)
    With[{psT2 = PauliStabilizer[1]["H", 1]["T", 1]["T", 1]},
        {Head[psT2], psT2["Length"]}
    ]
];

block["StabilizerFrame-T-on-zero-equals-zero",
    (* T|0> materialization: should equal |0> after FullSimplify *)
    Module[{psT = PauliStabilizer[1]["T", 1], vec},
        vec = Normal @ psT["StateVector"];
        Chop @ N @ FullSimplify[vec - {1, 0}]
    ]
];


(* ============================================================================ *)
(* InnerProduct method (Phase 6.5: was StabilizerInnerProduct)                  *)
(* ============================================================================ *)

block["InnerProduct-self",
    With[{psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]},
        psBell["InnerProduct", psBell]
    ]
];

block["InnerProduct-orthogonal",
    Module[{psPhiPlus, psPhiMinus},
        psPhiPlus = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
        psPhiMinus = psPhiPlus["Z", 1];
        Chop @ N @ psPhiPlus["InnerProduct", psPhiMinus]
    ]
];

block["InnerProduct-frame",
    (* Inner product between StabilizerFrame and PauliStabilizer *)
    With[{psT = PauliStabilizer[1]["T", 1], ps0 = PauliStabilizer[1]},
        FullSimplify @ psT["InnerProduct", ps0]
    ]
];


(* ============================================================================ *)
(* Expectation method (Phase 6.5: was StabilizerExpectation)                    *)
(* ============================================================================ *)

block["Expectation-stabilizer-element",
    With[{psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]},
        <|
            "<XX>" -> psBell["Expectation", "XX"],
            "<ZZ>" -> psBell["Expectation", "ZZ"],
            "<YY>" -> psBell["Expectation", "YY"]
        |>
    ]
];

block["Expectation-anticommuting",
    With[{psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]},
        psBell["Expectation", "XI"]
    ]
];

block["Expectation-5Q-syndromes",
    With[{ps5 = PauliStabilizer["5QubitCode"]},
        ps5["Expectation", #] & /@ {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"}
    ]
];


(* ============================================================================ *)
(* GraphState / LocalComplement                                                 *)
(* ============================================================================ *)

block["GraphState-from-Graph",
    With[{gs = GraphState[Graph[Range[3], {1 \[UndirectedEdge] 2, 2 \[UndirectedEdge] 3}]]},
        <|
            "Vertices"    -> gs["VertexCount"],
            "Edges"       -> gs["EdgeCount"],
            "Stabilizers" -> gs["Stabilizers"]
        |>
    ]
];

block["GraphState-cluster5",
    GraphState[Graph[Range[5], Table[i \[UndirectedEdge] (i + 1), {i, 4}]]]["Stabilizers"]
];

block["GraphState-to-PauliStabilizer",
    With[{gs = GraphState[Graph[Range[3], {1 \[UndirectedEdge] 2, 2 \[UndirectedEdge] 3}]]},
        gs["PauliStabilizer"]["Stabilizers"]
    ]
];

block["LocalComplement-star",
    With[{g = Graph[{1, 2, 3, 4}, {1 \[UndirectedEdge] 2, 1 \[UndirectedEdge] 3, 1 \[UndirectedEdge] 4}]},
        Sort @ EdgeList @ LocalComplement[g, 1]
    ]
];

block["LocalComplement-involutive",
    With[{g = Graph[Range[5], Table[i \[UndirectedEdge] (i + 1), {i, 4}]]},
        Sort @ EdgeList @ LocalComplement[LocalComplement[g, 3], 3] === Sort @ EdgeList[g]
    ]
];



Print["=== ALL BLOCKS COMPLETE ==="];
