#!/usr/bin/env wolframscript

(* Verifier for OngoingProjects/Stabilizer/Documentation/synthesis-implementation.md.
   Each block-id heading marks a code example referenced in the MD; the verifier
   re-runs every block and prints labeled output. The MD's output blocks should
   match this script's output verbatim (modulo formatting).

   Run from repo root:
     wolframscript -f OngoingProjects/Stabilizer/Documentation/verify-synthesis-implementation.wls
*)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

block[id_, expr_] := Module[{result},
    result = expr;
    Print["=== ", id, " ==="];
    Print[ToString[result, InputForm]];
    Print[]
];


(* ============================================================================ *)
(* \[Section]1.1 \[Dash] Pauli string                                           *)
(* ============================================================================ *)

block["1.1-symplectic-encoding",
    (* Bell state (n=2 stabilizers on 2 qubits) -- no auto-pad ellipsis *)
    With[{ps = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]},
        Association[
            "Stabilizers" -> ps["Stabilizers"],
            "X-bits-shape" -> Dimensions[ps["X"]],
            "X-bits" -> ps["X"],
            "Z-bits" -> ps["Z"],
            "Phase" -> ps["Phase"]
        ]
    ]
];

block["1.1-multiplication-via-symplectic",
    (* X*Z = -iY: in symplectic encoding, (1,0) XOR (0,1) = (1,1) = Y bits.
       Demonstrate by composing two single-qubit stabilizers: X and Z stabilize disjoint
       eigenspaces; their composition's tableau bits XOR. *)
    With[{psX = PauliStabilizer[{"X"}], psZ = PauliStabilizer[{"Z"}]},
        {psX["Stabilizers"], psZ["Stabilizers"], psX["Tableau"][[All, 1, 1]] + psZ["Tableau"][[All, 1, 1]]}
    ]
];


(* ============================================================================ *)
(* \[Section]1.2 \[Dash] Tableau with destabilizers                             *)
(* ============================================================================ *)

block["1.2-tableau-shape",
    With[{ps = PauliStabilizer["5QubitCode"]},
        <|
            "Qubits" -> ps["Qubits"],
            "GeneratorCount" -> ps["GeneratorCount"],
            "TableauDimensions" -> Dimensions[ps["Tableau"]],
            "MatrixDimensions" -> Dimensions[ps["Matrix"]],
            "DestabilizerExample" -> ps["Destabilizers"][[1]],
            "StabilizerExample" -> ps["Stabilizers"][[1]]
        |>
    ]
];

block["1.2-AG-symplectic-invariant",
    (* AarGot04 Prop 2 invariant: M . Omega . M^T == Omega (mod 2) over F_2 *)
    With[{
        ps = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]],
        n = 2,
        omega = ArrayFlatten[{
            {ConstantArray[0, {2, 2}], IdentityMatrix[2]},
            {IdentityMatrix[2], ConstantArray[0, {2, 2}]}
        }]
    },
        With[{m = ps["Matrix"]},
            Mod[m . omega . Transpose[m] - omega, 2]
        ]
    ]
];


(* ============================================================================ *)
(* \[Section]1.3 \[Dash] Graph state representation (Anders & Briegel)          *)
(* ============================================================================ *)

block["1.3-graph-state-stabilizers",
    With[{gs = GraphState[Graph[Range[5], Table[i \[UndirectedEdge] (i + 1), {i, 4}]]]},
        <|
            "Vertices" -> gs["VertexCount"],
            "Edges" -> gs["EdgeCount"],
            "Stabilizers" -> gs["Stabilizers"]
        |>
    ]
];

block["1.3-graph-state-equals-cluster-circuit",
    (* GraphState stabilizers should match the cluster-state circuit
       (apply H to all qubits, then CZ between adjacent qubits). *)
    Module[{gs, ps, psFromCirc, n = 4},
        gs = GraphState[Graph[Range[n], Table[i \[UndirectedEdge] (i + 1), {i, n - 1}]]];
        ps = gs["PauliStabilizer"];
        psFromCirc = PauliStabilizer[QuantumCircuitOperator[
            Join[Table["H" -> i, {i, n}], Table["CZ" -> {i, i + 1}, {i, n - 1}]]
        ]];
        ps["Stabilizers"] === psFromCirc["Stabilizers"]
    ]
];


(* ============================================================================ *)
(* \[Section]1.5 \[Dash] Stabilizer frame                                       *)
(* ============================================================================ *)

block["1.5-stabilizer-frame-from-T-gate",
    With[{psT = PauliStabilizer[1]["T", 1]},
        <|
            "Head" -> Head[psT],
            "Length" -> psT["Length"],
            "Coefficients" -> psT["Coefficients"]
        |>
    ]
];

block["1.5-frame-closes-under-Clifford",
    (* T followed by H stays in StabilizerFrame; the frame still has 2 components. *)
    Module[{psPlus, psTH},
        psPlus = PauliStabilizer[1]["H", 1];
        psTH = psPlus["T", 1]["H", 1];
        {Head[psTH], psTH["Length"]}
    ]
];

block["1.5-frame-materialization-T-on-zero",
    (* T|0> = |0> exactly (eigenstate). Materialized vector after FullSimplify
       should be {1, 0}. *)
    Module[{psT, vec},
        psT = PauliStabilizer[1]["T", 1];
        vec = Normal @ psT["StateVector"];
        Chop @ N @ FullSimplify[vec - {1, 0}]
    ]
];


(* ============================================================================ *)
(* \[Section]2.1 \[Dash] Clifford gate updates                                  *)
(* ============================================================================ *)

block["2.1-Heisenberg-H-on-Z",
    (* H takes Z to X via Heisenberg conjugation. Starting from |0> (stab Z),
       applying H gives stabilizer X. *)
    PauliStabilizer[1]["H", 1]["Stabilizers"]
];

block["2.1-Heisenberg-CNOT-XI-to-XX",
    (* CNOT takes X*I to X*X. Starting from H|00> (stab {XI, IZ}), applying CNOT_12
       gives stab {XX, ZZ} (= Bell state). *)
    Sort @ PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["Stabilizers"]
];


(* ============================================================================ *)
(* \[Section]2.3 \[Dash] Local complementation                                  *)
(* ============================================================================ *)

block["2.3-local-complement-star-to-wheel",
    (* LC at the center of a star K_{1,3} produces the complete graph on the leaves
       (joined to the center). *)
    Module[{g = Graph[{1, 2, 3, 4}, {1 \[UndirectedEdge] 2, 1 \[UndirectedEdge] 3, 1 \[UndirectedEdge] 4}]},
        Sort @ EdgeList @ LocalComplement[g, 1]
    ]
];

block["2.3-local-complement-involutive",
    (* AndBri05 Def 1: LC is involutive. *)
    Module[{g = Graph[Range[5], Table[i \[UndirectedEdge] (i + 1), {i, 4}]]},
        Sort @ EdgeList @ LocalComplement[LocalComplement[g, 3], 3] === Sort @ EdgeList[g]
    ]
];


(* ============================================================================ *)
(* \[Section]2.4 \[Dash] Measurement (Z-basis and Pauli-string)                 *)
(* ============================================================================ *)

block["2.4-Z-measurement-Bell",
    (* Bell['M', 1] returns Association of conditional outcomes. *)
    With[{psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]},
        Sort @ Keys @ psBell["M", 1]
    ]
];

block["2.4-Pauli-string-measurement-5Q",
    (* 5-qubit code state has all 5 stabilizer measurements deterministic = +1
       (outcome bit = 0). *)
    With[{ps5Q = PauliStabilizer["5QubitCode"]},
        {
            Sort @ Keys @ ps5Q["M", "XZZXI"],
            Sort @ Keys @ ps5Q["M", "IXZZX"],
            Sort @ Keys @ ps5Q["M", "XIXZZ"],
            Sort @ Keys @ ps5Q["M", "ZXIXZ"]
        }
    ]
];


(* ============================================================================ *)
(* \[Section]2.5 \[Dash] Inner products                                         *)
(* ============================================================================ *)

block["2.5-inner-product-orthogonal",
    Module[{psPhiPlus, psPhiMinus},
        psPhiPlus = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
        psPhiMinus = psPhiPlus["Z", 1];
        Chop @ N @ StabilizerInnerProduct[psPhiPlus, psPhiMinus]
    ]
];

block["2.5-inner-product-self",
    With[{psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]},
        StabilizerInnerProduct[psBell, psBell]
    ]
];

block["2.5-stabilizer-expectation-on-Bell",
    With[{psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]},
        <|
            "<XX>" -> StabilizerExpectation[psBell, "XX"],
            "<ZZ>" -> StabilizerExpectation[psBell, "ZZ"],
            "<YY>" -> StabilizerExpectation[psBell, "YY"],
            "<XI>" -> StabilizerExpectation[psBell, "XI"]
        |>
    ]
];


(* ============================================================================ *)
(* \[Section]2.7 \[Dash] Counting and enumeration                               *)
(* ============================================================================ *)

block["2.7-stabilizer-state-count",
    (* AarGot04 Prop 1: N(n) = 2^n Product[2^(n-k) + 1, {k, 0, n-1}] *)
    Function[n, 2^n Product[2^(n - k) + 1, {k, 0, n - 1}]] /@ Range[1, 4]
];

block["2.7-clifford-group-order",
    (* KoeSmo14 Eq 2: |C_n| = 2^(n^2 + 2n) Product[4^j - 1, {j, 1, n}] *)
    Function[n, 2^(n^2 + 2 n) Product[4^j - 1, {j, 1, n}]] /@ Range[1, 3]
];


(* ============================================================================ *)
(* \[Section]2.8 \[Dash] Random Clifford                                        *)
(* ============================================================================ *)

block["2.8-random-clifford",
    Block[{},
        SeedRandom[20260430];
        With[{ps = RandomClifford[3]},
            <|"Qubits" -> ps["Qubits"], "Stabilizers" -> ps["Stabilizers"]|>
        ]
    ]
];

block["2.8-random-clifford-uniformity",
    (* 200 samples on n=1 should hit most of the 24 elements. *)
    Block[{},
        SeedRandom[12345];
        Length @ DeleteDuplicates @ Table[
            With[{r = RandomClifford[1]}, {r["Stabilizers"], r["Signs"]}],
            {200}
        ]
    ]
];


(* ============================================================================ *)
(* \[Section]2.12 \[Dash] QEC code extraction                                   *)
(* ============================================================================ *)

block["2.12-5Q-syndromes",
    (* For each of 15 single-qubit Pauli errors {X_i, Y_i, Z_i : i=1..5}, the
       4-bit syndrome (against the first 4 stabilizers) must be unique. *)
    Module[{n = 5, gens, errors, syndromes, sympIP},
        sympIP[v1_, v2_] := Mod[v1[[;; n]] . v2[[n + 1 ;; 2 n]] + v2[[;; n]] . v1[[n + 1 ;; 2 n]], 2];
        gens = (Module[{xs, zs},
                {xs, zs} = Transpose @ Replace[Characters[#],
                    {"I" -> {0, 0}, "X" -> {1, 0}, "Y" -> {1, 1}, "Z" -> {0, 1}}, {1}];
                Join[xs, zs]
            ]) & /@ Take[PauliStabilizer["5QubitCode"]["Stabilizers"], 4];
        errors = Flatten[Table[
            With[{e = ConstantArray[0, 2 n]},
                Switch[op,
                    "X", ReplacePart[e, i -> 1],
                    "Y", ReplacePart[e, {i -> 1, n + i -> 1}],
                    "Z", ReplacePart[e, n + i -> 1]
                ]
            ],
            {i, n}, {op, {"X", "Y", "Z"}}], 1];
        syndromes = Table[Table[sympIP[err, g], {g, gens}], {err, errors}];
        <|"NumErrors" -> Length[errors], "DistinctSyndromes" -> Length @ Union[syndromes]|>
    ]
];


(* ============================================================================ *)
(* \[Section]3.1-3.2 \[Dash] Symbolic Pauli arithmetic + symbolic measurement   *)
(* ============================================================================ *)

block["3.2-symbolic-measurement-allocates-fresh-symbol",
    (* StabilizerMeasure on H|0> allocates a fresh symbol \[FormalS][k]. *)
    Module[{psPlus = PauliStabilizer[1]["H", 1], psSym, syms},
        psSym = StabilizerMeasure[psPlus, 1];
        syms = DeleteDuplicates @ Cases[psSym["Phase"], _\[FormalS], Infinity];
        <|
            "Head" -> Head[psSym],
            "FreshSymbols" -> syms,
            "Phase" -> psSym["Phase"]
        |>
    ]
];

block["3.2-substitute-outcomes-roundtrip",
    (* Substituting symbol -> 0 or 1 reproduces the regular ["M"] outcome branches. *)
    Module[{psPlus, psSym, sym, psSub0, psSub1, ps0, ps1},
        psPlus = PauliStabilizer[1]["H", 1];
        psSym = StabilizerMeasure[psPlus, 1];
        sym = First @ DeleteDuplicates @ Cases[psSym["Phase"], _\[FormalS], Infinity];
        psSub0 = SubstituteOutcomes[psSym, sym -> 0];
        psSub1 = SubstituteOutcomes[psSym, sym -> 1];
        ps0 = psPlus["M", 1][0];
        ps1 = psPlus["M", 1][1];
        <|
            "outcome 0 stabilizers match" -> (psSub0["Stabilizers"] === ps0["Stabilizers"]),
            "outcome 1 stabilizers match" -> (psSub1["Stabilizers"] === ps1["Stabilizers"])
        |>
    ]
];


(* ============================================================================ *)
(* Integration: QuantumCircuitOperator routing                                  *)
(* ============================================================================ *)

block["integration-method-stabilizer",
    QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}][Method -> "Stabilizer"]["Stabilizers"]
];

block["integration-named-circuit-GHZ",
    QuantumCircuitOperator["GHZ"[3]][Method -> "Stabilizer"]["Stabilizers"]
];


Print["=== ALL BLOCKS COMPLETE ==="];
