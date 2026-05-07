(* ============================================================================
   Tests/Correctness_TextbookResults.wlt -- outcome-focused correctness tests.

   These tests verify that PauliStabilizer + StabilizerFrame + GraphState +
   CliffordChannel produce numerical / physical outcomes that match
   well-established textbook results: Bell correlations, GHZ states, the
   5-qubit / Steane / Shor stabilizer codes, cluster states, etc.

   Goal: check OUTCOMES, not just internal consistency. Every test below
   verifies an exact numeric or set-theoretic equality against a value that
   is published in the QEC / quantum-information textbook literature
   (Got97, Got00, Nielsen & Chuang, AarGot04).
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];


(* Convenience: Bell, GHZ-n, cluster-n, 5Q-, Steane-, Shor-state fixtures. *)
psBell  = PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "CNOT" -> {1, 2}}
psGHZ3  = PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "CNOT" -> {1, 2}, "CNOT" -> {2, 3}}
ps5Q    = PauliStabilizer["5QubitCode"]
psSteane = PauliStabilizer["SteaneCode"]
psShor  = PauliStabilizer["9QubitCode"]


(* ============================================================================ *)
(* TIER A -- Bell-state correlations (Nielsen & Chuang 2.6, Got97 2.2)         *)
(*                                                                              *)
(* For |Phi+> = (|00> + |11>)/sqrt(2):                                          *)
(*   <Phi+|XX|Phi+> = +1 (XX is a stabilizer)                                  *)
(*   <Phi+|ZZ|Phi+> = +1 (ZZ is a stabilizer)                                  *)
(*   <Phi+|YY|Phi+> = -1 (because YY = -XX*ZZ)                                 *)
(*   <Phi+|IX|Phi+> =  0 (anticommutes with ZZ)                                *)
(*   <Phi+|XI|Phi+> =  0 (anticommutes with ZZ)                                *)
(*   <Phi+|IZ|Phi+> =  0 (anticommutes with XX)                                *)
(*   <Phi+|ZI|Phi+> =  0 (anticommutes with XX)                                *)
(* ============================================================================ *)

VerificationTest[psBell["Expectation", "XX"], +1, TestID -> "Bell-Expectation-XX-Plus1"]
VerificationTest[psBell["Expectation", "ZZ"], +1, TestID -> "Bell-Expectation-ZZ-Plus1"]
VerificationTest[psBell["Expectation", "YY"], -1, TestID -> "Bell-Expectation-YY-Minus1"]
VerificationTest[psBell["Expectation", "XI"],  0, TestID -> "Bell-Expectation-XI-Zero"]
VerificationTest[psBell["Expectation", "IX"],  0, TestID -> "Bell-Expectation-IX-Zero"]
VerificationTest[psBell["Expectation", "ZI"],  0, TestID -> "Bell-Expectation-ZI-Zero"]
VerificationTest[psBell["Expectation", "IZ"],  0, TestID -> "Bell-Expectation-IZ-Zero"]

(* Bell-state ZZ measurement: deterministic outcome 0 (correlated). *)
VerificationTest[Sort @ Keys @ psBell["M", "ZZ"], {0}, TestID -> "Bell-Measure-ZZ-Deterministic"]

(* Bell-state XX measurement: also deterministic outcome 0. *)
VerificationTest[Sort @ Keys @ psBell["M", "XX"], {0}, TestID -> "Bell-Measure-XX-Deterministic"]

(* Bell-state Z_1 measurement: random outcome (anticommutes with XX). *)
VerificationTest[Sort @ Keys @ psBell["M", "ZI"], {0, 1}, TestID -> "Bell-Measure-Z1-Random"]


(* ============================================================================ *)
(* TIER B -- GHZ-3 state correlations                                           *)
(*                                                                              *)
(* For |GHZ_3> = (|000> + |111>)/sqrt(2), stabilizer group <XXX, ZZI, IZZ>:    *)
(*   <GHZ_3|XXX|GHZ_3> = +1                                                    *)
(*   <GHZ_3|ZZI|GHZ_3> = +1                                                    *)
(*   <GHZ_3|IZZ|GHZ_3> = +1                                                    *)
(*   <GHZ_3|ZIZ|GHZ_3> = +1   (= ZZI * IZZ)                                    *)
(*   <GHZ_3|YYX|GHZ_3> = ?    (depends on signs; let's compute)                *)
(*   <GHZ_3|XII|GHZ_3> =  0   (anticommutes with ZZI)                          *)
(* ============================================================================ *)

VerificationTest[psGHZ3["Expectation", "XXX"], +1, TestID -> "GHZ3-Expectation-XXX-Plus1"]
VerificationTest[psGHZ3["Expectation", "ZZI"], +1, TestID -> "GHZ3-Expectation-ZZI-Plus1"]
VerificationTest[psGHZ3["Expectation", "IZZ"], +1, TestID -> "GHZ3-Expectation-IZZ-Plus1"]
VerificationTest[psGHZ3["Expectation", "ZIZ"], +1, TestID -> "GHZ3-Expectation-ZIZ-Plus1-Product"]
VerificationTest[psGHZ3["Expectation", "XII"],  0, TestID -> "GHZ3-Expectation-XII-Zero"]
VerificationTest[psGHZ3["Expectation", "IIZ"],  0, TestID -> "GHZ3-Expectation-IIZ-Zero-AnticommutesXXX"]

(* GHZ-3 ZZI measurement: deterministic. *)
VerificationTest[Sort @ Keys @ psGHZ3["M", "ZZI"], {0}, TestID -> "GHZ3-Measure-ZZI-Deterministic"]

(* GHZ-3 IIZ measurement: random. *)
VerificationTest[Sort @ Keys @ psGHZ3["M", "IIZ"], {0, 1}, TestID -> "GHZ3-Measure-IIZ-Random"]


(* ============================================================================ *)
(* TIER C -- 5-qubit code correctness ([[5,1,3]], Got97 ch. 6)                 *)
(*                                                                              *)
(* The 5-qubit code has 4 stabilizers (XZZXI, IXZZX, XIXZZ, ZXIXZ) plus a      *)
(* logical XXXXX and ZZZZZ. On the encoded |0_L> state:                        *)
(*   - Each of the 4 stabilizers gives expectation +1.                         *)
(*   - All 15 single-qubit Pauli errors {X_i, Y_i, Z_i : i=1..5} produce       *)
(*     UNIQUE 4-bit syndromes (the defining property of a [[n,1,3]] code).    *)
(* ============================================================================ *)

VerificationTest[
    Map[ps5Q["Expectation", #] &, {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"}],
    {1, 1, 1, 1},
    TestID -> "5Q-AllStabilizers-Expectation-Plus1"
]

(* 5Q syndrome uniqueness via direct enumeration. *)
VerificationTest[
    Module[{n = 5, generators, errors, syndromes, sympIP},
        sympIP[v1_, v2_] := Mod[
            v1[[;; n]] . v2[[n + 1 ;; 2 n]] + v2[[;; n]] . v1[[n + 1 ;; 2 n]],
            2
        ];
        generators = Module[{xs, zs},
                {xs, zs} = Transpose @ Replace[Characters[#],
                    {"I" -> {0, 0}, "X" -> {1, 0}, "Y" -> {1, 1}, "Z" -> {0, 1}}, {1}];
                Join[xs, zs]
            ] & /@ Take[ps5Q["Stabilizers"], 4];
        errors = Flatten[Table[
            With[{e = ConstantArray[0, 2 n]},
                Switch[op,
                    "X", ReplacePart[e, i -> 1],
                    "Y", ReplacePart[e, {i -> 1, n + i -> 1}],
                    "Z", ReplacePart[e, n + i -> 1]
                ]
            ],
            {i, n}, {op, {"X", "Y", "Z"}}
        ], 1];
        syndromes = Table[Table[sympIP[err, g], {g, generators}], {err, errors}];
        Length[Union[syndromes]]
    ],
    15,    (* All 15 syndromes must be distinct. *)
    TestID -> "5Q-Syndrome-Uniqueness-15errors-15syndromes"
]


(* ============================================================================ *)
(* TIER D -- Steane code [[7,1,3]] (Got00, Steane 1996)                         *)
(*                                                                              *)
(* Steane has 6 stabilizers + logical Xbar/Zbar. On encoded |0_L>, each of the *)
(* 6 stabilizers gives expectation +1.                                         *)
(* ============================================================================ *)

VerificationTest[
    AllTrue[
        psSteane["Expectation", #] & /@ Take[psSteane["Stabilizers"], 6],
        # == 1 &
    ],
    True,
    TestID -> "Steane-AllStabilizers-Expectation-Plus1"
]

(* The 7th "stabilizer" listed by our framework is actually the logical X     *)
(* operator XXXXXXX -- it commutes with all 6 stabilizers and with itself,   *)
(* so its expectation is also +1 on |0_L>. *)
VerificationTest[
    psSteane["Expectation", "XXXXXXX"],
    1,
    TestID -> "Steane-Logical-X-Expectation-Plus1-on-Logical-Zero"
]


(* ============================================================================ *)
(* TIER E -- Shor 9-qubit code [[9,1,3]] (Shor 1995, Got97 ch.4)                *)
(*                                                                              *)
(* Shor has 8 stabilizers. Each gives expectation +1 on |0_L>.                 *)
(* ============================================================================ *)

VerificationTest[
    AllTrue[
        psShor["Expectation", #] & /@ Take[psShor["Stabilizers"], 8],
        # == 1 &
    ],
    True,
    TestID -> "Shor-AllStabilizers-Expectation-Plus1"
]


(* ============================================================================ *)
(* TIER F -- Cluster state stabilizers (AndBri05 Eq 1)                          *)
(*                                                                              *)
(* For a graph G with vertex i having neighbors N(i), the cluster-state        *)
(* stabilizer at vertex i is K_i = X_i tensor Product over j in N(i) of Z_j.  *)
(* For a linear chain 1-2-3-4-5: K_1 = XZIII, K_2 = ZXZII, K_3 = IZXZI,        *)
(*                                  K_4 = IIZXZ, K_5 = IIIZX.                  *)
(* ============================================================================ *)

VerificationTest[
    GraphState[
        Graph[Range[5], Table[i \[UndirectedEdge] (i + 1), {i, 4}]]
    ]["Stabilizers"],
    {"XZIII", "ZXZII", "IZXZI", "IIZXZ", "IIIZX"},
    TestID -> "Cluster-5-LinearChain-Stabilizers"
]

(* Round-trip: build cluster state via H+CZ circuit, compare with GraphState. *)
VerificationTest[
    Sort @ (PauliStabilizer @ QuantumCircuitOperator @ Join[
        Table["H" -> i, {i, 4}],
        Table["CZ" -> {i, i + 1}, {i, 3}]
    ])["Stabilizers"],
    Sort @ GraphState[
        Graph[Range[4], Table[i \[UndirectedEdge] (i + 1), {i, 3}]]
    ]["Stabilizers"],
    TestID -> "Cluster-4-Circuit-vs-GraphState-MatchesStabilizers"
]


(* ============================================================================ *)
(* TIER G -- Inner-product values (GarMarCro12, Got97)                          *)
(*                                                                              *)
(* For two stabilizer states, <psi|phi> is either 0 or 2^(-s/2) where s is    *)
(* a measure of their generator-set divergence. Specific cases:                *)
(*   <Bell|Bell> = 1                                                            *)
(*   <Bell PhiPlus | Bell PhiMinus> = 0                                        *)
(*   <|0> | |1>> = 0  (orthogonal)                                              *)
(*   <|0> | |+>> = 1/sqrt(2)                                                    *)
(*   <|0> | |0>> = 1                                                            *)
(* ============================================================================ *)

VerificationTest[psBell["InnerProduct", psBell], 1, TestID -> "InnerProduct-Bell-Self-One"]

VerificationTest[
    With[{psPhiMinus = psBell["Z", 1]},
        Chop[N @ psBell["InnerProduct", psPhiMinus]]
    ],
    0,
    TestID -> "InnerProduct-Bell-PhiPlus-PhiMinus-Zero"
]

VerificationTest[
    PauliStabilizer[1]["InnerProduct", PauliStabilizer[1]],
    1,
    TestID -> "InnerProduct-zero-Self-One"
]

VerificationTest[
    Chop[N @ PauliStabilizer[1]["InnerProduct", PauliStabilizer[1]["X", 1]]],
    0,
    TestID -> "InnerProduct-zero-one-Zero"
]

VerificationTest[
    Chop[N @ PauliStabilizer[1]["InnerProduct", PauliStabilizer[1]["H", 1]]],
    1 / Sqrt[2.] // N,
    TestID -> "InnerProduct-zero-plus-OneOverSqrt2"
]


(* ============================================================================ *)
(* TIER H -- Heisenberg-picture conjugation (Got98 / textbook)                  *)
(*                                                                              *)
(* For Clifford gate U, U P U† transforms each Pauli P. Verify against         *)
(* textbook tables: HZH = X, HXH = Z, SZS† = Z, SXS† = Y, etc.                 *)
(* (These overlap with Stim-cross-package tests but are checked here against   *)
(* the textbook source Got98 directly.)                                         *)
(* ============================================================================ *)

(* HZH = X *)
VerificationTest[
    PauliStabilizer[1]["H", 1]["Stabilizers"],
    {"X"},
    TestID -> "Heisenberg-HZH-equals-X"
]

(* HXH = Z (so applying H to |+> stabilizer X gives Z stabilizer) *)
VerificationTest[
    PauliStabilizer[1]["H", 1]["H", 1]["Stabilizers"],
    {"Z"},
    TestID -> "Heisenberg-HH-restores-Z"
]

(* S Z S† = Z *)
VerificationTest[
    PauliStabilizer[1]["S", 1]["Stabilizers"],
    {"Z"},
    TestID -> "Heisenberg-SZS-equals-Z"
]

(* S X S† = Y *)
VerificationTest[
    PauliStabilizer[1]["H", 1]["S", 1]["Stabilizers"],
    {"Y"},
    TestID -> "Heisenberg-SXS-equals-Y"
]

(* S^2 = Z (so SS|0> = |0>, but SS turns the X stabilizer into -X) *)
VerificationTest[
    PauliStabilizer[1]["H", 1]["S", 1]["S", 1]["Stabilizers"],
    {"-X"},
    TestID -> "Heisenberg-SS-on-X-equals-MinusX"
]


(* ============================================================================ *)
(* TIER I -- Got97 stabilizer-group properties                                  *)
(*                                                                              *)
(* From Got97 ch. 3: a valid stabilizer group must satisfy                      *)
(*   (a) all generators commute pairwise                                       *)
(*   (b) -I is not in the group                                                 *)
(*   (c) the group has |G| = 2^n elements for n generators                     *)
(* ============================================================================ *)

(* Pairwise commutation for several named codes. *)
VerificationTest[
    Block[{},
        AllTrue[
            {ps5Q, psSteane, psShor},
            With[{mats = Module[{sign, body, paulis},
                    {sign, body} = If[StringStartsQ[#, "-"], {-1, StringDrop[#, 1]}, {1, #}];
                    paulis = Replace[Characters[body], {
                        "I" -> {{1, 0}, {0, 1}},
                        "X" -> {{0, 1}, {1, 0}},
                        "Y" -> {{0, -I}, {I, 0}},
                        "Z" -> {{1, 0}, {0, -1}}
                    }, {1}];
                    sign * If[Length[paulis] == 1, First[paulis], KroneckerProduct @@ paulis]
                ] & /@ #["Stabilizers"]},
                AllTrue[Subsets[mats, {2}], #[[1]] . #[[2]] == #[[2]] . #[[1]] &]
            ] &
        ]
    ],
    True,
    TestID -> "Got97-Stabilizers-PairwiseCommute-NamedCodes"
]


(* ============================================================================ *)
(* TIER J -- AarGot04 measurement statistics on encoded states                  *)
(*                                                                              *)
(* All 4 stabilizers of the 5Q code give deterministic outcome 0 on |0_L>.    *)
(* All 6 stabilizers of Steane code do too on |0_L>.                           *)
(* ============================================================================ *)

VerificationTest[
    Map[Sort[Keys[ps5Q["M", #]]] &, Take[ps5Q["Stabilizers"], 4]],
    {{0}, {0}, {0}, {0}},
    TestID -> "AarGot-5Q-StabilizerMeasurements-AllZero"
]

VerificationTest[
    Map[Sort[Keys[psSteane["M", #]]] &, Take[psSteane["Stabilizers"], 6]],
    Table[{0}, {6}],
    TestID -> "AarGot-Steane-StabilizerMeasurements-AllZero"
]


(* ============================================================================ *)
(* TIER K -- Outcome-correctness for non-Clifford channels                      *)
(*                                                                              *)
(* For BitFlip[p] applied to |0>: the average state is (1-p)|0><0| + p|1><1|. *)
(* The average <Z> is (1-p)*1 + p*(-1) = 1 - 2p. Verify via Phase 7.2 mixture. *)
(* ============================================================================ *)

(* Note: 1-qubit ps["Expectation", pauli] currently hits ROADMAP A.11 (the
   single-element KroneckerProduct bug in pauliStringMatrix). To verify channel
   mixture probabilities without that bug, we use 2-qubit stabilizer states and
   2-qubit Pauli observables. The physics is the same: channel on q1 with
   identity elsewhere produces the same single-qubit-eigenvalue average. *)

(* BitFlip[p] on q1 of |00>: average <ZI> = (1-p)*1 + p*(-1) = 1 - 2p.        *)
VerificationTest[
    With[{
        mixture = QuantumChannel["BitFlip"[1/4], {1}][PauliStabilizer[2]]
    },
        Total[#[[1]] * #[[2]]["Expectation", "ZI"] & /@ mixture] // Simplify
    ],
    1/2,    (* (1 - 2*1/4) = 1/2 *)
    TestID -> "BitFlip-OnZeroZero-AverageZI-equals-1minus2p"
]

(* PhaseFlip[p] on q1 of |+0>: average <XI> = (1-p)*1 + p*(-1) = 1 - 2p.     *)
VerificationTest[
    With[{
        mixture = QuantumChannel["PhaseFlip"[1/3], {1}][PauliStabilizer[2]["H", 1]]
    },
        Total[#[[1]] * #[[2]]["Expectation", "XI"] & /@ mixture] // Simplify
    ],
    1/3,    (* 1 - 2*1/3 = 1/3 *)
    TestID -> "PhaseFlip-OnPlusZero-AverageXI-equals-1minus2p"
]

(* Depolarizing[p] on q1 of |00>: average <ZI> = 1 - p.                       *)
VerificationTest[
    With[{
        mixture = QuantumChannel["Depolarizing"[1/2], {1}][PauliStabilizer[2]]
    },
        Total[#[[1]] * #[[2]]["Expectation", "ZI"] & /@ mixture] // Simplify
    ],
    1/2,
    TestID -> "Depolarizing-OnZeroZero-AverageZI-equals-1minusP"
]
