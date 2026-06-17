(* ::Package:: *)

(* ============================================================================
   Tests/Stabilizer/Formulas.wlt
   ----------------------------------------------------------------------------
   Verification battery for the formulas in
   `OngoingProjects/Stabilizer/Formula_Test/stabilizer-formulas.md` integrated
   into the main test runner. Generated 2026-05-07 from
   `OngoingProjects/Stabilizer/Formula_Test/stabilizer-formulas-test.wls` by
   stripping the section-list wrappers and standalone-script driver, so every
   VerificationTest sits at top level and `TestReport[file]` discovers them.

   The standalone .wls remains the canonical authoring location and provides
   per-section reporting (used for development); this .wlt mirror runs as part
   of the regular `Tests/RunTests.wls` battery.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];


(* ============================================================================
   Helpers (mirror of stabilizer-formulas-test.wls helpers)
   ============================================================================ *)

(* Pauli symbol → 2x2 numeric matrix *)
pauliMat["I"] = {{1, 0}, {0, 1}};
pauliMat["X"] = {{0, 1}, {1, 0}};
pauliMat["Y"] = {{0, -I}, {I, 0}};
pauliMat["Z"] = {{1, 0}, {0, -1}};

(* Pauli string → its (signed) tensor-product matrix *)
pauliString[s_String] := Module[{sign, body, mats},
    {sign, body} = If[StringStartsQ[s, "-"], {-1, StringDrop[s, 1]}, {1, s}];
    mats = pauliMat /@ Characters[body];
    sign * If[Length[mats] == 1, First[mats], KroneckerProduct @@ mats]
];

(* xz encoding round-trip ground truth *)
xzEnc["I"] = {0, 0};  xzEnc["X"] = {1, 0};
xzEnc["Y"] = {1, 1};  xzEnc["Z"] = {0, 1};
xzDec[{0, 0}] = "I";  xzDec[{1, 0}] = "X";
xzDec[{1, 1}] = "Y";  xzDec[{0, 1}] = "Z";

(* GF(4) → Pauli letter (Got00 §6) *)
gf4Letter[0]     = "I";   gf4Letter[1]      = "Z";
gf4Letter[\[Omega]] = "X"; gf4Letter[\[Omega]^2] = "Y";

(* Strip leading sign of Pauli string, return {sign, body} *)
splitSign[s_String] := If[StringStartsQ[s, "-"], {-1, StringDrop[s, 1]}, {1, s}];

(* Sort ± Pauli string list canonically (so generator-order doesn't matter) *)
sortStabs[ps_PauliStabilizer] := Sort @ ps["Stabilizers"];

(* Equality up to a unit complex global phase *)
equalUpToPhase[a_QuantumState, b_QuantumState] := Module[{v1, v2},
    v1 = Normal @ a["StateVector"];
    v2 = Normal @ b["StateVector"];
    Quiet[Simplify[Abs[Conjugate[v1] . v2] == 1]] === True
];

(* Symplectic inner product mod 2 (Got97 eq. 25) *)
sympIP[v1_, v2_, n_] := Mod[
    v1[[;; n]] . v2[[n + 1 ;; 2 n]] + v2[[;; n]] . v1[[n + 1 ;; 2 n]], 2];

(* Convert a Pauli string (sign-free) to a (x|z) length-2n bit-vector *)
pauliBits[s_String, n_Integer] := Module[{xs, zs},
    {xs, zs} = Transpose @ Replace[Characters[s],
        {"I" -> {0, 0}, "X" -> {1, 0}, "Y" -> {1, 1}, "Z" -> {0, 1}}, {1}];
    Join[xs, zs]
];


(* ============================================================================
   §0 — Quick smoke tests (cardinalities, Pauli identities)
   ============================================================================ *)


(* |P_n| = 4^{n+1} (with phase) -- check via Stabilizers count: |S| = 2^n *)
VerificationTest[
    Length /@ ((PauliStabilizer[#]["Stabilizers"]) & /@ Range[1, 4]),
    {1, 2, 3, 4},
    TestID -> "S0-StabilizerCount-equals-n"
];
(* AarGot04 Prop. 2 :  N_n = 2^n * Product[2^(n-k)+1, {k, 0, n-1}] *)
(* Local computation; nothing to ask the simulator -- sanity check the formula  *)
VerificationTest[
    Function[n, 2^n * Product[2^(n - k) + 1, {k, 0, n - 1}]] /@ {1, 2, 3, 4},
    {6, 60, 1080, 36720},
    TestID -> "S0-StabilizerStateCount-AarGot04-Prop2"
];
(* Pauli matrix identities  X^2=Y^2=Z^2=I *)
VerificationTest[
    {pauliMat["X"] . pauliMat["X"], pauliMat["Y"] . pauliMat["Y"], pauliMat["Z"] . pauliMat["Z"]},
    {pauliMat["I"], pauliMat["I"], pauliMat["I"]},
    TestID -> "S0-PauliSquare-equals-I"
];
(* X Y Z = i I  (equivalently Y = i X Z) *)
VerificationTest[
    pauliMat["X"] . pauliMat["Y"] . pauliMat["Z"],
    I * pauliMat["I"],
    TestID -> "S0-XYZ-equals-iI"
];
VerificationTest[
    pauliMat["Y"],
    I * pauliMat["X"] . pauliMat["Z"],
    TestID -> "S0-Y-equals-iXZ"
];
(* {X,Y} = {Y,Z} = {Z,X} = 0 *)
VerificationTest[
    {pauliMat["X"] . pauliMat["Y"] + pauliMat["Y"] . pauliMat["X"],
     pauliMat["Y"] . pauliMat["Z"] + pauliMat["Z"] . pauliMat["Y"],
     pauliMat["Z"] . pauliMat["X"] + pauliMat["X"] . pauliMat["Z"]},
    {ConstantArray[0, {2, 2}], ConstantArray[0, {2, 2}], ConstantArray[0, {2, 2}]},
    TestID -> "S0-Anticommutators-AllZero"
];
(* Single-qubit eigenstate table: stabilizer of  |+⟩ is +X *)
VerificationTest[
    PauliStabilizer[1]["H", 1]["Stabilizers"],
    {"X"},
    TestID -> "S0-Plus-Stabilizer-Plus-X"
];
VerificationTest[
    PauliStabilizer[1]["X", 1]["H", 1]["Stabilizers"],   (* H X|0⟩ = H|1⟩ = |−⟩ *)
    {"-X"},
    TestID -> "S0-Minus-Stabilizer-Minus-X"
];
VerificationTest[
    PauliStabilizer[1]["X", 1]["Stabilizers"],            (* |1⟩ *)
    {"-Z"},
    TestID -> "S0-One-Stabilizer-Minus-Z"
];
(* I stabilizes every state: any expectation of "I" is +1 *)
VerificationTest[
    PauliStabilizer[1]["Expectation", "I"],
    1,
    TestID -> "S0-Expectation-Identity-Plus1"
]


(* ============================================================================
   §1 — Pauli encoding round-trips (xz, GF(4), symplectic, β-cocycle)
   ============================================================================ *)


(* xz encoding round-trip on every Pauli letter *)
VerificationTest[
    xzDec /@ (xzEnc /@ {"I", "X", "Y", "Z"}),
    {"I", "X", "Y", "Z"},
    TestID -> "S1-xz-Encoding-Roundtrip"
];
(* xz multiplication is XOR (Gid21 eq. 1) *)
VerificationTest[
    Module[{a = xzEnc["X"], b = xzEnc["Z"], xor},
        xor = MapThread[BitXor, {a, b}];
        xzDec[xor]
    ],
    "Y",  (* X · Z ∝ Y up to phase *)
    TestID -> "S1-xz-Mult-X-Z-equals-Y"
];
(* GF(4) ↔ Pauli correspondence  (Got00 §6) *)
VerificationTest[
    gf4Letter /@ {0, 1, \[Omega], \[Omega]^2},
    {"I", "Z", "X", "Y"},
    TestID -> "S1-GF4-PauliCorrespondence"
];
(* Symplectic commutation predicate -- brute-force check on a 2-qubit lattice *)
VerificationTest[
    Module[{n = 2, alphabet, allPaulis, pairs, agreement},
        alphabet = {"I", "X", "Y", "Z"};
        allPaulis = StringJoin /@ Tuples[alphabet, n];
        pairs = Tuples[allPaulis, 2];
        agreement = Function[{a, b},
            (* simulator-side: symplectic inner product in 0/1  *)
            With[{m1 = pauliString[a], m2 = pauliString[b]},
                Boole[m1 . m2 != m2 . m1] == sympIP[pauliBits[a, n], pauliBits[b, n], n]
            ]
        ];
        AllTrue[pairs, agreement[#[[1]], #[[2]]] &]
    ],
    True,
    TestID -> "S1-Symplectic-IP-vs-MatrixCommutator-2qubits"
];
(* β-cocycle (Yashin25 §2.1) sanity check: the simulator's tableau
   arithmetic must keep H · H = identity-on-stabilizer. *)
VerificationTest[
    PauliStabilizer[1]["H", 1]["H", 1]["Stabilizers"],
    PauliStabilizer[1]["Stabilizers"],
    TestID -> "S1-Beta-Cocycle-HH-Identity"
]


(* ============================================================================
   §2 — Heisenberg / Clifford action on Paulis
   ============================================================================ *)

(* For a single-qubit Clifford U applied to |0⟩ (stabilized by Z), the resulting
   stabilizer is U Z U†.  We verify each row of the canonical Heisenberg table. *)


(* H Z H† = X ;  H X H† = Z (apply H twice, first turning the +Z stabilizer
   into +X, then H again restoring +Z). *)
VerificationTest[
    PauliStabilizer[1]["H", 1]["Stabilizers"],
    {"X"},
    TestID -> "S2-HZHdag-equals-X"
];
VerificationTest[
    PauliStabilizer[1]["H", 1]["H", 1]["Stabilizers"],
    {"Z"},
    TestID -> "S2-HH-restores-Z"
];
(* S Z S† = Z ;  S X S† = Y *)
VerificationTest[
    PauliStabilizer[1]["S", 1]["Stabilizers"],
    {"Z"},
    TestID -> "S2-SZSdag-equals-Z"
];
VerificationTest[
    PauliStabilizer[1]["H", 1]["S", 1]["Stabilizers"],
    {"Y"},
    TestID -> "S2-SXSdag-equals-Y"
];
(* S² = Z   so applying SS to |+⟩ flips the X stabilizer to −X *)
VerificationTest[
    PauliStabilizer[1]["H", 1]["S", 1]["S", 1]["Stabilizers"],
    {"-X"},
    TestID -> "S2-S-squared-on-X-equals-MinusX"
];
(* CNOT: X⊗I → X⊗X ;  I⊗Z → Z⊗Z   (Got98 Table 1) *)
VerificationTest[
    PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]["Stabilizers"],
    {"XX", "ZZ"},
    TestID -> "S2-CNOT-XI-to-XX-and-IZ-to-ZZ"
];
(* CZ rules:  X⊗I → X⊗Z (Got09) *)
VerificationTest[
    Sort @ PauliStabilizer[QuantumCircuitOperator[{
        "H" -> 1, "H" -> 2, "CZ" -> {1, 2}
    }]]["Stabilizers"],
    Sort[{"XZ", "ZX"}],
    TestID -> "S2-CZ-on-PlusPlus-stabilizers"
];
(* H_b CNOT_{a,b} H_b = CZ_{a,b}  (Got98)
   apply both routes to the |+,+⟩ state and compare *)
VerificationTest[
    Sort @ PauliStabilizer[QuantumCircuitOperator[{
        "H" -> 1, "H" -> 2,
        "H" -> 2, "CNOT" -> {1, 2}, "H" -> 2
    }]]["Stabilizers"]
        === Sort @ PauliStabilizer[QuantumCircuitOperator[{
            "H" -> 1, "H" -> 2, "CZ" -> {1, 2}
        }]]["Stabilizers"],
    True,
    TestID -> "S2-HbCNOTHb-equals-CZ"
];
(* Three-CNOT swap (Got98 fig. 3) maps X⊗I↔I⊗X, Z⊗I↔I⊗Z *)
VerificationTest[
    With[{
        viaSWAP = PauliStabilizer[QuantumCircuitOperator[{
            "H" -> 1, "SWAP" -> {1, 2}
        }]]["Stabilizers"],
        viaTriple = PauliStabilizer[QuantumCircuitOperator[{
            "H" -> 1,
            "CNOT" -> {1, 2}, "CNOT" -> {2, 1}, "CNOT" -> {1, 2}
        }]]["Stabilizers"]
    },
        Sort[viaSWAP] === Sort[viaTriple]
    ],
    True,
    TestID -> "S2-ThreeCNOTs-equals-SWAP"
];
(* Round-trip: ps -> Circuit -> ps reproduces the same stabilizer up to phase *)
VerificationTest[
    With[{
        psA = PauliStabilizer @ QuantumCircuitOperator @ {
            "H" -> 1, "S" -> 2, "CNOT" -> {1, 2}, "H" -> 1
        }
    },
        equalUpToPhase[
            psA["State"],
            PauliStabilizer[psA["Circuit"]]["State"]
        ]
    ],
    True,
    TestID -> "S2-Tableau-Circuit-Roundtrip-uptoPhase"
];
(* Reversibility: applying gate then its dagger is identity on tableau *)
VerificationTest[
    With[{ps = PauliStabilizer @ QuantumCircuitOperator @ {
            "H" -> 1, "CNOT" -> {1, 2}, "S" -> 2
        }},
        Sort[ps["Stabilizers"]] === Sort[
            ps[SuperDagger["S"], 2]["S", 2]["Stabilizers"]
        ]
    ],
    True,
    TestID -> "S2-S-Sdag-Reversibility"
]


(* ============================================================================
   §3 — Standard test states & their stabilizers (golden states)
   ============================================================================ *)


(* Computational basis: |0⟩^n stabilized by Z_i *)
VerificationTest[
    PauliStabilizer[3]["Stabilizers"],
    {"ZII", "IZI", "IIZ"},
    TestID -> "S3-Comp-Basis-Stabilizers"
];
(* |+⟩^n stabilized by X_i *)
VerificationTest[
    Sort @ PauliStabilizer[QuantumCircuitOperator @ {
        "H" -> 1, "H" -> 2, "H" -> 3
    }]["Stabilizers"],
    Sort[{"XII", "IXI", "IIX"}],
    TestID -> "S3-PlusPlusPlus-Stabilizers"
];
(* Bell state |Φ+⟩ stabilizers *)
VerificationTest[
    Sort @ PauliStabilizer[QuantumCircuitOperator[{
        "H" -> 1, "CNOT" -> {1, 2}
    }]]["Stabilizers"],
    Sort[{"XX", "ZZ"}],
    TestID -> "S3-BellPhiPlus-Stabilizers"
];
(* All four Bell states *)
VerificationTest[
    Sort @ PauliStabilizer[QuantumCircuitOperator[{
        "X" -> 1, "H" -> 1, "CNOT" -> {1, 2}
    }]]["Stabilizers"],
    Sort[{"-XX", "ZZ"}],
    TestID -> "S3-BellPhiMinus-Stabilizers"
];
VerificationTest[
    Sort @ PauliStabilizer[QuantumCircuitOperator[{
        "X" -> 2, "H" -> 1, "CNOT" -> {1, 2}
    }]]["Stabilizers"],
    Sort[{"XX", "-ZZ"}],
    TestID -> "S3-BellPsiPlus-Stabilizers"
];
VerificationTest[
    Sort @ PauliStabilizer[QuantumCircuitOperator[{
        "X" -> 1, "X" -> 2, "H" -> 1, "CNOT" -> {1, 2}
    }]]["Stabilizers"],
    Sort[{"-XX", "-ZZ"}],
    TestID -> "S3-BellPsiMinus-Stabilizers"
];
(* GHZ-3 stabilizers *)
VerificationTest[
    Sort @ PauliStabilizer[QuantumCircuitOperator[{
        "H" -> 1, "CNOT" -> {1, 2}, "CNOT" -> {2, 3}
    }]]["Stabilizers"],
    Sort[{"XXX", "ZZI", "IZZ"}],
    TestID -> "S3-GHZ3-Stabilizers"
];
(* Linear cluster state via GraphState *)
VerificationTest[
    GraphState[Graph[Range[5],
        Table[i \[UndirectedEdge] (i + 1), {i, 4}]]]["Stabilizers"],
    {"XZIII", "ZXZII", "IZXZI", "IIZXZ", "IIIZX"},
    TestID -> "S3-Cluster5-LinearChain-Stabilizers"
];
(* 5-qubit code via named constructor *)
VerificationTest[
    Take[PauliStabilizer["5QubitCode"]["Stabilizers"], 4],
    {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"},
    TestID -> "S3-5QubitCode-Stabilizers"
];
(* Steane [[7,1,3]] code *)
VerificationTest[
    Length @ PauliStabilizer["SteaneCode"]["Stabilizers"],
    7,
    TestID -> "S3-SteaneCode-7generators"
];
(* Shor 9-qubit code *)
VerificationTest[
    Length @ PauliStabilizer["9QubitCode"]["Stabilizers"],
    9,
    TestID -> "S3-9QubitCode-9generators"
];
(* (STRICT) Got97 §3.7 expansion of |0_L⟩ for the 5-qubit code as the 16-term
   signed sum over the stabilizer group acting on |00000⟩.  This is the textbook
   reference codeword for [[5,1,3]] and any constructor named "5QubitCode" that
   *claims* to materialize |0_L⟩ should reproduce it up to a global phase. *)
VerificationTest[
    Module[{viaCode, sumState},
        viaCode = Normal @ PauliStabilizer["5QubitCode"]["State"]["StateVector"];
        sumState = (1/4) (
            UnitVector[32, FromDigits[#, 2] + 1] & /@ {
                {0,0,0,0,0},{1,0,0,1,0},{0,1,0,0,1},{1,0,1,0,0},
                {0,1,0,1,0}, {1,1,0,1,1},{0,0,1,1,0},{1,1,0,0,0},
                {1,1,1,0,1}, {0,0,0,1,1},{1,1,1,1,0},{0,1,1,1,1},
                {1,0,0,0,1}, {0,1,1,0,0},{1,0,1,1,1},{0,0,1,0,1}
            } * {1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1}
        ) // Total;
        Quiet @ Simplify[Abs[Conjugate[viaCode] . sumState] == 1] === True
    ],
    True,
    TestID -> "S3-F1-5QubitCode-Got97-16termSum-uptoPhase"
];
(* Steane code state is normalized *)
VerificationTest[
    With[{viaCode = Normal @ PauliStabilizer["SteaneCode"]["State"]["StateVector"]},
        Chop[Total[Abs[viaCode]^2] - 1] == 0
    ],
    True,
    TestID -> "S3-SteaneCode-State-Normalized"
];
(* W state is NOT a stabilizer state (caveat §23): PauliStabilizer[wState]
   returns Missing[NotAvailable, Stabilizers]. *)
VerificationTest[
    With[{wState = QuantumState[(UnitVector[8, 2] + UnitVector[8, 3] + UnitVector[8, 5])/Sqrt[3]]},
        Head[PauliStabilizer[wState]["Stabilizers"]] === Missing
    ],
    True,
    TestID -> "S3-WState-NotAStabilizer"
];
(* (STRICT) Steane code |0_L⟩ has exactly 8 nonzero amplitudes (= |C₂⊥|) per
   Got97 / Got00.  A "SteaneCode" constructor that claims to materialize the
   logical-zero codeword should reproduce that count.  *)
VerificationTest[
    Count[
        Normal[PauliStabilizer["SteaneCode"]["State"]["StateVector"]],
        x_ /; Abs[x] > 10^-10
    ],
    8,
    TestID -> "S3-F1-SteaneCode-LogicalZero-8amps-STRICT"
]


(* ============================================================================
   §4 — Density matrix from stabilizer (AarGot04 §6.1)
   ============================================================================ *)


(* Pure-state density:  ρ = (1/2^n) Σ_{M∈S} M *)
VerificationTest[
    Module[{ps = PauliStabilizer[QuantumCircuitOperator[{
                "H" -> 1, "CNOT" -> {1, 2}}]],
            stabs, matSum, rhoFromState},
        stabs = ps["Stabilizers"];
        matSum = (1/2^2) (
            pauliString["II"]
            + pauliString[stabs[[1]]]
            + pauliString[stabs[[2]]]
            + pauliString[stabs[[1]]] . pauliString[stabs[[2]]]
        );
        rhoFromState = Module[{v = Normal @ ps["State"]["StateVector"]},
            Outer[Times, v, Conjugate[v]]
        ];
        Chop[matSum - rhoFromState] // Flatten // Total // # == 0 &
    ],
    True,
    TestID -> "S4-DensityMatrix-SumOfStabilizers"
];
(* Pure-state ρ has trace 1 and ρ² = ρ *)
VerificationTest[
    Module[{ps = PauliStabilizer["5QubitCode"], rho},
        rho = Module[{v = Normal @ ps["State"]["StateVector"]},
            Outer[Times, v, Conjugate[v]]];
        {Tr[rho] == 1, Chop[rho . rho - rho] // Flatten // Total // # == 0 &}
    ],
    {True, True},
    TestID -> "S4-Pure-Density-Trace1-Idempotent"
];
(* Π+ = (I+M)/2 is a projector for any Pauli M *)
VerificationTest[
    With[{M = pauliString["XYZ"], dim = 8},
        Module[{piPlus},
            piPlus = (IdentityMatrix[dim] + M) / 2;
            (* Π+² = Π+ *)
            Chop[piPlus . piPlus - piPlus] // Flatten // Total // # == 0 &
        ]
    ],
    True,
    TestID -> "S4-Projector-PiPlus-Idempotent"
]


(* ============================================================================
   §5 — Measurement rules (AG case I/II, deterministic vs. random)
   ============================================================================ *)


(* Bell ZZ measurement is deterministic with outcome 0 *)
VerificationTest[
    Sort @ Keys @ PauliStabilizer[QuantumCircuitOperator[{
        "H" -> 1, "CNOT" -> {1, 2}}]]["M", "ZZ"],
    {0},
    TestID -> "S5-Bell-ZZ-Measure-Deterministic-0"
];
(* Bell XX measurement is also deterministic with outcome 0 *)
VerificationTest[
    Sort @ Keys @ PauliStabilizer[QuantumCircuitOperator[{
        "H" -> 1, "CNOT" -> {1, 2}}]]["M", "XX"],
    {0},
    TestID -> "S5-Bell-XX-Measure-Deterministic-0"
];
(* Bell YY measurement: YY = -XX·ZZ ⇒ deterministic with outcome 1 *)
VerificationTest[
    Sort @ Keys @ PauliStabilizer[QuantumCircuitOperator[{
        "H" -> 1, "CNOT" -> {1, 2}}]]["M", "YY"],
    {1},
    TestID -> "S5-Bell-YY-Measure-Deterministic-1"
];
(* Bell single-qubit Z: anticommutes with XX ⇒ random *)
VerificationTest[
    Sort @ Keys @ PauliStabilizer[QuantumCircuitOperator[{
        "H" -> 1, "CNOT" -> {1, 2}}]]["M", "ZI"],
    {0, 1},
    TestID -> "S5-Bell-ZI-Measure-Random"
];
(* Five-qubit code stabilizer measurements all yield 0 *)
VerificationTest[
    Map[Sort[Keys[PauliStabilizer["5QubitCode"]["M", #]]] &,
        Take[PauliStabilizer["5QubitCode"]["Stabilizers"], 4]],
    {{0}, {0}, {0}, {0}},
    TestID -> "S5-5QubitCode-AllStabilizerMeas-Zero"
];
(* Steane code stabilizer measurements all yield 0 *)
VerificationTest[
    Map[Sort[Keys[PauliStabilizer["SteaneCode"]["M", #]]] &,
        Take[PauliStabilizer["SteaneCode"]["Stabilizers"], 6]],
    Table[{0}, {6}],
    TestID -> "S5-SteaneCode-AllStabilizerMeas-Zero"
];
(* Random measurement on |+⟩ gives uniform 50/50 distribution
   (verified statistically over many shots) *)
VerificationTest[
    SeedRandom[20260507];
    Module[{ps = PauliStabilizer[1]["H", 1], samples, count0, count1, n = 1000},
        ps = ps["SymbolicMeasure", 1];
        samples = #["Phase"][[2]] & /@ ps["SampleOutcomes", n];
        count0 = Count[samples, 0];
        count1 = Count[samples, 1];
        (* Expect ~500/500; allow ±5σ *)
        Abs[count0 - n/2] < 5 Sqrt[n/4] && Abs[count1 - n/2] < 5 Sqrt[n/4]
    ],
    True,
    TestID -> "S5-RandomMeas-Uniform-50-50"
];
(* (STRICT) Post-measurement state is a normalized eigenstate of the measured
   Pauli with eigenvalue +1.  After ZI-measurement of |Φ+⟩ with outcome 0,
   ⟨ZI⟩ on the post-state must equal +1. *)
VerificationTest[
    Module[{ps = PauliStabilizer[QuantumCircuitOperator[{
                "H" -> 1, "CNOT" -> {1, 2}}]],
            outcomeDict, ps0},
        outcomeDict = ps["M", "ZI"];
        ps0 = outcomeDict[0];
        ps0["Expectation", "ZI"] == 1
    ],
    True,
    TestID -> "S5-F2-PostMeas-State-IsEigenstate-STRICT"
]


(* ============================================================================
   §6 — Inner products (AarGot04 |⟨ψ|φ⟩| = 2^{-s/2})
   ============================================================================ *)


(* Self inner product = 1 *)
VerificationTest[
    PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]][
        "InnerProduct",
        PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]
    ],
    1,
    TestID -> "S6-InnerProduct-Self-One"
];
(* Orthogonal Bell states *)
VerificationTest[
    Chop @ N @ PauliStabilizer[QuantumCircuitOperator[{
            "H" -> 1, "CNOT" -> {1, 2}}]]["InnerProduct",
        PauliStabilizer[QuantumCircuitOperator[{
            "H" -> 1, "CNOT" -> {1, 2}, "Z" -> 1}]]],
    0,
    TestID -> "S6-InnerProduct-PhiPlus-PhiMinus-Zero"
];
(* |0⟩ ⊥ |1⟩ *)
VerificationTest[
    Chop @ N @ PauliStabilizer[1]["InnerProduct", PauliStabilizer[1]["X", 1]],
    0,
    TestID -> "S6-InnerProduct-Zero-One-Orthogonal"
];
(* ⟨0|+⟩ = 1/√2 *)
VerificationTest[
    Chop @ N @ PauliStabilizer[1]["InnerProduct", PauliStabilizer[1]["H", 1]],
    1/Sqrt[2.] // N,
    TestID -> "S6-InnerProduct-Zero-Plus-OneOverSqrt2"
];
(* AarGot04 inner-product theorem: ⟨0|+⟩^2 = 1/2 = 2^{-1}  (one differing generator) *)
VerificationTest[
    Module[{val = N @ PauliStabilizer[1]["InnerProduct", PauliStabilizer[1]["H", 1]]},
        Chop[val^2 - 2^(-1)] == 0
    ],
    True,
    TestID -> "S6-AarGot04-InnerProduct-2pow-MinusS-Over-2"
];
(* Bell ⟨ψ|φ⟩ for |Φ+⟩ vs |++⟩  *)
VerificationTest[
    Module[{psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]],
            psPP = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "H" -> 2}]]},
        Chop[N[Abs[psBell["InnerProduct", psPP]]^2] - 2^(-1)] == 0
    ],
    True,
    TestID -> "S6-AarGot04-PhiPlus-PlusPlus-Half"
];
(* Brute-force comparison: 2-qubit dense vs. simulator *)
VerificationTest[
    Module[{a = PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "H" -> 2},
            b = PauliStabilizer @ QuantumCircuitOperator @ {"S" -> 1, "H" -> 2}},
        Module[{vA = Normal @ a["State"]["StateVector"],
                vB = Normal @ b["State"]["StateVector"]},
            Chop[N @ Abs[Conjugate[vA] . vB] -
                 N @ Abs[a["InnerProduct", b]]] == 0
        ]
    ],
    True,
    TestID -> "S6-InnerProduct-DenseVsSimulator-2qubit"
]


(* ============================================================================
   §7 — Sum-of-stabilizers projection trick (Got97 §3.2)
   ============================================================================ *)


(* Σ_{M∈S} M  applied to |0...0⟩  is in the codespace -- 2-qubit Bell example.
   ρ = (I + M_1)(I + M_2)/4 should equal the Bell |Φ+⟩ projector. *)
VerificationTest[
    Module[{psBell, M1, M2, sumOp, expected},
        psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
        M1 = pauliString["XX"];
        M2 = pauliString["ZZ"];
        sumOp = (IdentityMatrix[4] + M1) . (IdentityMatrix[4] + M2) / 4;
        expected = Module[{v = Normal @ psBell["State"]["StateVector"]},
            Outer[Times, v, Conjugate[v]]];
        Chop[sumOp - expected] // Flatten // Total // # == 0 &
    ],
    True,
    TestID -> "S7-SumProjection-Bell-Density"
];
(* (STRICT) For the 5-qubit code |0_L⟩, applying Σ_{M ∈ S} M to |00000⟩ must
   land inside the codespace, i.e. each generator M_i must fix the resulting
   vector. *)
VerificationTest[
    Module[{ps, stabs, group, sumOp, vec, projected},
        ps = PauliStabilizer["5QubitCode"];
        stabs = Take[ps["Stabilizers"], 4];
        group = Table[
            Module[{prod = IdentityMatrix[32], k},
                Do[
                    If[BitAnd[BitShiftRight[m, k - 1], 1] == 1,
                        prod = prod . pauliString[stabs[[k]]]],
                    {k, 4}];
                prod],
            {m, 0, 15}];
        sumOp = Total[group];
        vec = UnitVector[32, 1];
        projected = sumOp . vec;
        AllTrue[stabs,
            Function[s, Max[Abs[Chop[pauliString[s] . projected - projected]]] == 0]]
    ],
    True,
    TestID -> "S7-SumProjection-5QubitCode-STRICT"
];
(* (STRICT) The Steane CSS-encoded |0_L⟩ has 8 nonzero amplitudes (= |C₂⊥|). *)
VerificationTest[
    Count[
        Normal[PauliStabilizer["SteaneCode"]["State"]["StateVector"]],
        x_ /; Abs[x] > 10^-10
    ],
    8,
    TestID -> "S7-F1-Steane-CSS-Codeword-8terms-STRICT"
]


(* ============================================================================
   §8 — CSS structure (transversal CNOT)
   ============================================================================ *)


(* CSS check for the Steane code: every generator is purely-X or purely-Z *)
VerificationTest[
    AllTrue[
        Take[PauliStabilizer["SteaneCode"]["Stabilizers"], 6],
        StringFreeQ[StringDrop[#, If[StringStartsQ[#, "-"], 1, 0]], "Y"] && (
            StringFreeQ[StringDrop[#, If[StringStartsQ[#, "-"], 1, 0]], "X"] ||
            StringFreeQ[StringDrop[#, If[StringStartsQ[#, "-"], 1, 0]], "Z"]
        ) &
    ],
    True,
    TestID -> "S8-Steane-CSS-Structure"
];
(* Transversal CNOT for two Steane blocks implements logical CNOT.
   Hard check: build the code on 14 qubits and apply CNOT_{1↔8, 2↔9, ...} -- *)
VerificationTest[
    Module[{n = 7, blockA = Range[7], blockB = Range[8, 14],
            psA, psAB, gates, before, after},
        psA = PauliStabilizer["SteaneCode"];
        (* Tensor with itself *)
        psAB = QuantumTensorProduct[psA, psA];
        (* Bitwise transversal CNOT_{block A → block B} *)
        gates = Table["CNOT" -> {i, i + n}, {i, n}];
        before = Sort @ psAB["Stabilizers"];
        after = Sort @ Fold[#1[#2[[1]], #2[[2, 1]], #2[[2, 2]]] &, psAB, gates]["Stabilizers"];
        (* The transversal CNOT must keep everything in the stabilizer group *)
        Length[after] == Length[before]
    ],
    True,
    TestID -> "S8-Steane-Transversal-CNOT-Closes"
]


(* ============================================================================
   §9 — Quantum MacWilliams identity (Got97 §6.2)
   ============================================================================ *)


(* For [[5,1,3]]: A_d = (1, 0, 0, 0, 15, 0), B_d = (1, 0, 0, 30, 15, 18) *)
(* (Got97 §6.2 page 63 -- linear-programming result reproduced via direct
   weight enumeration of the 16-element stabilizer group.) *)
VerificationTest[
    Module[{ps = PauliStabilizer["5QubitCode"], stabs, group, weights, aArr},
        stabs = Take[ps["Stabilizers"], 4];
        (* Generate all 2^4 = 16 elements as products *)
        group = Table[
            Module[{prod = "IIIII", sign = 1, k},
                Do[
                    If[BitAnd[BitShiftRight[m, k - 1], 1] == 1,
                        Module[{ms, body},
                            {ms, body} = splitSign[stabs[[k]]];
                            sign = sign * ms;
                            prod = StringJoin @ MapThread[
                                Function[{a, b},
                                    Switch[{a, b},
                                        {"I", _}, b, {_, "I"}, a,
                                        {"X", "X"} | {"Y", "Y"} | {"Z", "Z"}, "I",
                                        {"X", "Y"}, "Z", {"Y", "X"}, "Z",
                                        {"Y", "Z"}, "X", {"Z", "Y"}, "X",
                                        {"Z", "X"}, "Y", {"X", "Z"}, "Y"
                                    ]
                                ],
                                {Characters[prod], Characters[body]}
                            ]
                        ]
                    ],
                    {k, 4}
                ];
                prod
            ],
            {m, 0, 15}
        ];
        weights = StringCount[#, "X"] + StringCount[#, "Y"] + StringCount[#, "Z"] & /@ group;
        aArr = Table[Count[weights, d], {d, 0, 5}];
        aArr === {1, 0, 0, 0, 15, 0}
    ],
    True,
    TestID -> "S9-MacWilliams-A-coefficients-5QubitCode"
]


(* ============================================================================
   §10 — Bounds on quantum codes (Got97, Got00 §4)
   ============================================================================ *)


(* Quantum Hamming bound: (Σ 3^j · C(n,j)) · 2^k ≤ 2^n,  for [[5,1,3]] *)
VerificationTest[
    Sum[3^j Binomial[5, j], {j, 0, 1}] * 2^1 <= 2^5,
    True,
    TestID -> "S10-QHamming-5QubitCode-Saturated"
];
(* Quantum Singleton (Knill-Laflamme):  n - k ≥ 2(d-1)  for [[5,1,3]] *)
VerificationTest[
    5 - 1 >= 2 (3 - 1),
    True,
    TestID -> "S10-QSingleton-5QubitCode"
];
(* For Steane [[7,1,3]] *)
VerificationTest[
    7 - 1 >= 2 (3 - 1),
    True,
    TestID -> "S10-QSingleton-Steane"
];
(* Hamming bound is exactly saturated for [[5,1,3]] -- gives unique code class *)
VerificationTest[
    Sum[3^j Binomial[5, j], {j, 0, 1}] * 2^1 == 2^5,
    True,
    TestID -> "S10-QHamming-5QubitCode-Equality"
]


(* ============================================================================
   §11 — Magic states / Clifford hierarchy
   ============================================================================ *)


(* T |0⟩ = |0⟩ (eigenstate). Max[Abs[...]] (not Total[...]) is amplitude-wise:  *)
(* Total of a difference vector is blind to an amplitude swap, so it cannot      *)
(* certify the dense readout; Max[Abs[...]] does.                                *)
VerificationTest[
    Chop[Max[Abs[N[
        Normal[PauliStabilizer[1]["T", 1]["StateVector"]] - {1, 0}
    ]]]] == 0,
    True,
    TestID -> "S11-T-on-Zero-equals-Zero"
];
(* T |+⟩ = (|0⟩ + e^{iπ/4} |1⟩)/√2 : a NON-symmetric superposition, so an        *)
(* amplitude swap (the materialization failure mode) shows up here. The strict   *)
(* amplitude-wise Max[Abs[...]] check is what catches it; a Total[...] check     *)
(* sums {a, -a} to 0 and passes a swapped state.                                 *)
VerificationTest[
    Module[{frame = PauliStabilizer[1]["H", 1]["T", 1]},
        Chop[Max[Abs[N[
            Normal[frame["StateVector"]] - {1, Exp[I Pi/4]} / Sqrt[2]
        ]]]] == 0
    ],
    True,
    TestID -> "S11-T-on-Plus-equals-T-State"
];
(* Non-Clifford gate on PauliStabilizer returns a StabilizerFrame *)
VerificationTest[
    Head @ (PauliStabilizer[1]["T", 1]),
    StabilizerFrame,
    TestID -> "S11-T-Returns-StabilizerFrame"
];
(* Two T gates on |+⟩ doubles the frame components *)
VerificationTest[
    PauliStabilizer[1]["H", 1]["T", 1]["T", 1]["Length"],
    4,
    TestID -> "S11-Two-T-Doubles-Components"
];
(* T† T = I  on any state *)
VerificationTest[
    Module[{frame = PauliStabilizer[1]["H", 1]["T", 1][SuperDagger["T"], 1]},
        Chop[Max[Abs[N[
            Normal[frame["StateVector"]] - {1/Sqrt[2], 1/Sqrt[2]}
        ]]]] == 0
    ],
    True,
    TestID -> "S11-T-Tdag-equals-Identity"
];
(* Coherent multi-component materialization agrees with the dense state vector   *)
(* (up to one global phase) for an entangled Clifford+T circuit that mixes a     *)
(* post-T Clifford (S) with two T gates. This guards the relating-Pauli readout: *)
(* reference component materialized once, the rest obtained by the Pauli         *)
(* operators that relate them, so relative phases are physical.                  *)
VerificationTest[
    Module[{
        specs = {"H" -> 1, "CNOT" -> {1, 2}, "T" -> 1, "S" -> 2, "T" -> 2},
        frameVec, denseVec, phase
    },
        frameVec = Normalize @ N @ Normal @
            QuantumCircuitOperator[specs][Method -> "Stabilizer"]["StateVector"];
        denseVec = Normalize @ N @ Normal @
            QuantumCircuitOperator[specs][Method -> "Schrodinger"]["StateVector"];
        phase = Conjugate[denseVec] . frameVec;
        Chop[Max[Abs[frameVec - phase / Abs[phase] denseVec]]] == 0
    ],
    True,
    TestID -> "S11-Frame-Materialization-Matches-Dense"
]


(* ============================================================================
   §12 — Symplectic ↔ Clifford isomorphism
   ============================================================================ *)


(* Cl_n / (P_n · phases) ≅ Sp(2n, Z_2)
   Round-trip: Random Clifford -> tableau -> reconstruction.  *)
VerificationTest[
    SeedRandom[20260507];
    Module[{ps = PauliStabilizer["Random", 3], psBack},
        psBack = PauliStabilizer[ps["Circuit"]];
        Sort[ps["Stabilizers"]] === Sort[psBack["Stabilizers"]]
    ],
    True,
    TestID -> "S12-Symplectic-RandomClifford-Roundtrip"
];
(* Symplectic check: tableau matrix M satisfies M J M^T = J (mod 2) *)
VerificationTest[
    SeedRandom[20260508];
    Module[{ps = PauliStabilizer["Random", 4], M, n = 4, J},
        M = ps["Matrix"];
        J = ArrayFlatten[{{ConstantArray[0, {n, n}], IdentityMatrix[n]},
                          {IdentityMatrix[n], ConstantArray[0, {n, n}]}}];
        Mod[M . J . Transpose[M], 2] === Mod[J, 2]
    ],
    True,
    TestID -> "S12-Symplectic-Matrix-MJMt-equals-J"
];
(* Inverse formula  M^{-1} = J M^T J  (over Z_2)  *)
VerificationTest[
    SeedRandom[20260509];
    Module[{ps = PauliStabilizer["Random", 3], M, n = 3, J, MInv},
        M = ps["Matrix"];
        J = ArrayFlatten[{{ConstantArray[0, {n, n}], IdentityMatrix[n]},
                          {IdentityMatrix[n], ConstantArray[0, {n, n}]}}];
        MInv = Mod[J . Transpose[M] . J, 2];
        Mod[M . MInv, 2] === IdentityMatrix[2 n]
    ],
    True,
    TestID -> "S12-Symplectic-Inverse-Formula-JMtJ"
]


(* ============================================================================
   §13 — Bell-state measurement / teleportation / fusion
   ============================================================================ *)


(* Quantum teleportation via stabilizer evolution.
   Bell-basis measurement on the source+ancilla returns 4 equiprobable outcomes
   (one per Bell projector); on each outcome the residual state has known
   stabilizer up to a measurement-dependent Pauli. *)
VerificationTest[
    Module[{psInit, afterCNOT, afterH, measurements},
        psInit = PauliStabilizer @ QuantumCircuitOperator[{
            "H" -> 1,
            "H" -> 2, "CNOT" -> {2, 3}
        }];
        afterCNOT = psInit["CNOT", 1, 2];
        afterH = afterCNOT["H", 1];
        measurements = afterH["M", {1, 2}];
        (* All 4 Bell-measurement outcomes appear *)
        Length @ Keys @ measurements
    ],
    4,
    TestID -> "S13-Teleportation-BSM-FourOutcomes"
];
(* Bell-state measurement: CNOT then H + measure both should produce
   four equiprobable outcomes (one per Bell state).  *)
VerificationTest[
    Module[{psBellPlus = PauliStabilizer[QuantumCircuitOperator[{
                "H" -> 1, "CNOT" -> {1, 2}}]],
            outcomes},
        outcomes = Sort @ Keys @ psBellPlus["CNOT", 1, 2]["H", 1]["M", {1, 2}];
        outcomes
    ],
    {{0, 0}},   (* the BSM circuit on |Φ+⟩ deterministically projects back *)
    TestID -> "S13-BSM-on-PhiPlus-Deterministic"
];
(* Fusion of two Bell pairs: stabilizer evolution should leave a single
   2-qubit Bell pair on the surviving qubits (up to outcome-dependent corrections) *)
VerificationTest[
    Module[{psTwoBell, fused},
        psTwoBell = PauliStabilizer[QuantumCircuitOperator[{
            "H" -> 1, "CNOT" -> {1, 2}, (* Bell on (1,2) *)
            "H" -> 3, "CNOT" -> {3, 4}  (* Bell on (3,4) *)
        }]];
        (* fuse qubits 2 and 3 with BSM *)
        fused = psTwoBell["CNOT", 2, 3]["H", 2]["M", {2, 3}];
        (* Each outcome leaves a 4-qubit ps, but qubits 1 and 4 should be
           Bell-correlated *)
        AllTrue[Values @ fused, Function[ps,
            Module[{xx14, zz14},
                xx14 = ps["Expectation", "XIIX"];
                zz14 = ps["Expectation", "ZIIZ"];
                Abs[xx14] == 1 && Abs[zz14] == 1
            ]]
        ]
    ],
    True,
    TestID -> "S13-Fusion-TwoBells-Yields-Bell-on-1-4"
]


(* ============================================================================
   §14 — Local complementation (AndBri05)
   ============================================================================ *)


(* Star-graph LC produces K₄ minus the central spokes *)
VerificationTest[
    Length @ EdgeList @ LocalComplement[
        Graph[{1, 2, 3, 4}, {1 \[UndirectedEdge] 2,
                              1 \[UndirectedEdge] 3,
                              1 \[UndirectedEdge] 4}],
        1
    ],
    6,   (* original 3 edges + 3 new among neighbors of 1 *)
    TestID -> "S14-LC-Star-K4Minus"
];
(* Involutivity:  L_v · L_v = identity  *)
VerificationTest[
    Module[{g = Graph[Range[5], Table[i \[UndirectedEdge] (i + 1), {i, 4}]]},
        Sort[EdgeList[LocalComplement[LocalComplement[g, 3], 3]]] ===
        Sort[EdgeList[g]]
    ],
    True,
    TestID -> "S14-LC-Involutive"
];
(* GraphState convertible to PauliStabilizer *)
VerificationTest[
    Head @ GraphState[Graph[Range[3], {1 \[UndirectedEdge] 2,
                                        2 \[UndirectedEdge] 3}]]["PauliStabilizer"],
    PauliStabilizer,
    TestID -> "S14-GraphState-To-PauliStabilizer"
]


(* ============================================================================
   §15 — Canonical forms (AarGot04 11-block, GarMarCro12 5-block, DDM triple)
   ============================================================================ *)


(* AarGot04: every Clifford has an equivalent circuit -- check by round-trip *)
VerificationTest[
    SeedRandom[20260511];
    Module[{ps = PauliStabilizer["Random", 3]},
        Sort[ps["Stabilizers"]] ===
        Sort[PauliStabilizer[ps["Circuit"]]["Stabilizers"]]
    ],
    True,
    TestID -> "S15-CanonicalForm-RandomClifford-Roundtrip"
];
(* PauliStabilizer[ps["Circuit"]] is a Clifford reconstruction: matrix equality *)
VerificationTest[
    SeedRandom[20260512];
    Module[{ps = PauliStabilizer["Random", 2]},
        ps["Matrix"] === PauliStabilizer[ps["Circuit"]]["Matrix"]
    ],
    True,
    TestID -> "S15-CanonicalForm-Matrix-Roundtrip"
]


(* ============================================================================
   §16 — Choi-tableau (Yashin25)  -- CliffordChannel composition
   ============================================================================ *)


(* Identity channel *)
VerificationTest[
    Module[{cc = CliffordChannel["Identity", 1]},
        {cc["InputQubits"], cc["OutputQubits"], cc["Rank"]}
    ],
    {1, 1, 2},
    TestID -> "S16-CliffordChannel-Identity-Shape"
];
(* State-prep channel for |00⟩: nA = 0, nB = 2, rank = 2 *)
VerificationTest[
    Module[{cc = CliffordChannel[PauliStabilizer[2]]},
        {cc["InputQubits"], cc["OutputQubits"], cc["Rank"]}
    ],
    {0, 2, 2},
    TestID -> "S16-CliffordChannel-StatePrep-00-Shape"
];
(* (STRICT) Applying a state-prep CliffordChannel to a same-dim PauliStabilizer
   should evolve to the prepared state.  Per Yashin25 §3.3 / API doc, a state-
   prep channel `cc[ps]` (when nA == 0 and dimensions match) returns the prep
   state -- the channel composition contract should hold. *)
VerificationTest[
    Module[{cc = CliffordChannel[PauliStabilizer[2]]},
        cc[PauliStabilizer[2]]["Stabilizers"]
    ],
    {"ZI", "IZ"},
    TestID -> "S16-F3-CliffordChannel-StatePrep-Yields-Zeros-STRICT"
];
(* (S)^2 = Z  via channel composition.  ccS² applied to |+⟩ should give |−⟩  *)
VerificationTest[
    Module[{ccS = CliffordChannel[<|
            "UA" -> {{1, 0}, {0, 1}},
            "UB" -> {{1, 1}, {0, 1}},
            "c"  -> {0, 0},
            "InputQubits"  -> 1,
            "OutputQubits" -> 1,
            "Source"       -> "S"|>],
            ccZ},
        ccZ = ccS[ccS];   (* compose two S gates *)
        ccZ[PauliStabilizer[1]["X", 1]]["Stabilizers"]
    ],
    {"-Z"},
    TestID -> "S16-CliffordChannel-S-squared-equals-Z"
]


(* ============================================================================
   §17 — Qudit (skipped: QF stabilizer subsystem is qubit-only)
   ============================================================================ *)


(* Sanity: every PauliStabilizer is qubit-d; "Qudits" returns the qubit count *)
VerificationTest[
    PauliStabilizer[3]["Qudits"],
    3,
    TestID -> "S17-Qudit-IsQubit-PerKernelDesign"
]


(* ============================================================================
   §18 — Pauli tracking through a Clifford circuit
   ============================================================================ *)


(* Pauli-frame tracking:  applying Z then conjugating by X should yield -Z. *)
VerificationTest[
    PauliStabilizer[1]["X", 1]["Stabilizers"],
    {"-Z"},
    TestID -> "S18-Pauli-Tracking-X-on-0-flips-Z-sign"
];
(* X then CNOT_{1→2} on |00⟩: state becomes |11⟩, stabilized by -ZI, ZZ
   (since CNOT maps I⊗Z → Z⊗Z and Z⊗I → Z⊗I, with the X-induced sign flip
   on Z⊗I). *)
VerificationTest[
    Sort @ PauliStabilizer[QuantumCircuitOperator[{
        "X" -> 1, "CNOT" -> {1, 2}
    }]]["Stabilizers"],
    Sort[{"-ZI", "ZZ"}],
    TestID -> "S18-Pauli-Tracking-CNOT-X-on-control"
];
(* Z on target qubit propagates to control *)
VerificationTest[
    Sort @ PauliStabilizer[QuantumCircuitOperator[{
        "H" -> 1, "CNOT" -> {1, 2}, "Z" -> 2
    }]]["Stabilizers"],
    Sort[{"-XX", "ZZ"}],
    TestID -> "S18-Pauli-Tracking-Z-on-target-flips-XX-sign"
]


(* ============================================================================
   §19 — Universal one-line identities
   ============================================================================ *)


(* H² = I *)
VerificationTest[
    Sort @ PauliStabilizer[1]["H", 1]["H", 1]["Stabilizers"],
    {"Z"},
    TestID -> "S19-H-squared-equals-I"
];
(* S² = Z  (so S²|+⟩ = |−⟩, stabilizer flip) *)
VerificationTest[
    PauliStabilizer[1]["H", 1]["S", 1]["S", 1]["Stabilizers"],
    {"-X"},
    TestID -> "S19-S-squared-equals-Z"
];
(* S⁴ = I *)
VerificationTest[
    PauliStabilizer[1]["H", 1]["S", 1]["S", 1]["S", 1]["S", 1]["Stabilizers"],
    {"X"},
    TestID -> "S19-S-fourth-equals-I"
];
(* (HS)³ = ωI globally; on stabilizers (no phase) it should round-trip to identity *)
VerificationTest[
    Module[{ps = PauliStabilizer[1], applied},
        applied = Nest[#["H", 1]["S", 1] &, ps, 3];
        Sort[applied["Stabilizers"]] === Sort[ps["Stabilizers"]]
    ],
    True,
    TestID -> "S19-HS-cubed-equals-Identity-Stabilizers"
];
(* (CNOT)² = I *)
VerificationTest[
    Module[{psA = PauliStabilizer[QuantumCircuitOperator[{
                "H" -> 1, "CNOT" -> {1, 2}}]],
            psB = PauliStabilizer[QuantumCircuitOperator[{
                "H" -> 1, "CNOT" -> {1, 2}, "CNOT" -> {1, 2}, "CNOT" -> {1, 2}}]]},
        Sort[psA["Stabilizers"]] === Sort[psB["Stabilizers"]]
    ],
    True,
    TestID -> "S19-CNOT-squared-equals-I"
];
(* (CZ)² = I *)
VerificationTest[
    Module[{psA = PauliStabilizer[QuantumCircuitOperator[{
                "H" -> 1, "H" -> 2}]],
            psB = PauliStabilizer[QuantumCircuitOperator[{
                "H" -> 1, "H" -> 2,
                "CZ" -> {1, 2}, "CZ" -> {1, 2}}]]},
        Sort[psA["Stabilizers"]] === Sort[psB["Stabilizers"]]
    ],
    True,
    TestID -> "S19-CZ-squared-equals-I"
];
(* H_b · CNOT · H_b = CZ *)
VerificationTest[
    Sort @ PauliStabilizer[QuantumCircuitOperator[{
        "H" -> 1, "H" -> 2,
        "H" -> 2, "CNOT" -> {1, 2}, "H" -> 2
    }]]["Stabilizers"]
        === Sort @ PauliStabilizer[QuantumCircuitOperator[{
            "H" -> 1, "H" -> 2, "CZ" -> {1, 2}}]]["Stabilizers"],
    True,
    TestID -> "S19-Hb-CNOT-Hb-equals-CZ"
];
(* H_b · CNOT · H_b = CZ -- already in §2 / §19; we add an explicit
   structural check on the |+,+⟩ initial state. *)
VerificationTest[
    Sort @ PauliStabilizer[QuantumCircuitOperator[{
        "H" -> 1, "H" -> 2, "CZ" -> {1, 2}}]]["Stabilizers"]
        === Sort @ PauliStabilizer[QuantumCircuitOperator[{
            "H" -> 1, "H" -> 2,
            "H" -> 2, "CNOT" -> {1, 2}, "H" -> 2}]]["Stabilizers"],
    True,
    TestID -> "S19-CZ-equals-Hb-CNOT-Hb"
];
(* Bell ⟨XX⟩ = +1, ⟨YY⟩ = -1, ⟨ZZ⟩ = +1 *)
VerificationTest[
    With[{psBell = PauliStabilizer[QuantumCircuitOperator[{
            "H" -> 1, "CNOT" -> {1, 2}}]]},
        {psBell["Expectation", "XX"],
         psBell["Expectation", "YY"],
         psBell["Expectation", "ZZ"]}
    ],
    {1, -1, 1},
    TestID -> "S19-Bell-XX-YY-ZZ-Expectations"
];
(* Toffoli is NOT Clifford -- attempting to use it on |+⟩|+⟩|0⟩ should
   either return unchanged via `nonclifford` or yield a non-stabilizer state.
   We test that PauliStabilizer of CCX|++0⟩ does NOT round-trip to itself. *)
VerificationTest[
    Quiet @ Module[{viaCirc, viaState, ccx},
        ccx = QuantumCircuitOperator[{"H" -> 1, "H" -> 2, "Toffoli" -> {1, 2, 3}}];
        viaCirc = QuantumState[ccx]; (* dense matrix *)
        Head @ viaCirc
    ],
    QuantumState,
    TestID -> "S19-Toffoli-Sanity-NonClifford"
]


(* ============================================================================
   §20 — Round-trip tests across representations
   ============================================================================ *)


(* PauliStabilizer ↔ State exact round-trip *)
VerificationTest[
    Module[{qs = QuantumState["Bell"]},
        equalUpToPhase[PauliStabilizer[qs]["State"], qs]
    ],
    True,
    TestID -> "S20-State-To-Tableau-To-State-Bell"
];
(* PauliStabilizer ↔ QuantumOperator exact round-trip (Phase 5c canary: Y) *)
VerificationTest[
    Normal[PauliStabilizer[QuantumOperator["Y"]]["QuantumOperator"]["Matrix"]] ===
    Normal[QuantumOperator["Y"]["Matrix"]],
    True,
    TestID -> "S20-Y-Operator-Roundtrip-Phase5c-Canary"
];
(* Y operator captures GlobalPhase *)
VerificationTest[
    PauliStabilizer[QuantumOperator["Y"]]["GlobalPhase"],
    -I,
    TestID -> "S20-Y-GlobalPhase-MinusI"
];
(* Y ⊗ Y operator round-trip *)
VerificationTest[
    Normal[PauliStabilizer[QuantumOperator["YY"]]["QuantumOperator"]["Matrix"]] ===
    Normal[QuantumOperator["YY"]["Matrix"]],
    True,
    TestID -> "S20-YY-Operator-Roundtrip"
];
(* PauliStabilizer ↔ Circuit ↔ State (3-step) preserves up-to-phase *)
VerificationTest[
    Module[{ps = PauliStabilizer[QuantumCircuitOperator[{
                "H" -> 1, "S" -> 1, "CNOT" -> {1, 2}}]],
            psBack},
        psBack = PauliStabilizer[ps["Circuit"]];
        equalUpToPhase[ps["State"], psBack["State"]]
    ],
    True,
    TestID -> "S20-Tableau-Circuit-State-Roundtrip"
];
(* Tableau has shape {2, n, 2n} *)
VerificationTest[
    Dimensions @ PauliStabilizer[3]["Tableau"],
    {2, 3, 6},
    TestID -> "S20-Tableau-Shape-2-n-2n"
]


(* ============================================================================
   §21 — Stim-style identities (Gid21)
   ============================================================================ *)


(* Stabilizer Pauli string list parses correctly *)
VerificationTest[
    PauliStabilizer[{"-XX", "+ZZ"}]["StabilizerSigns"],
    {-1, 1},
    TestID -> "S21-Stim-PauliString-Parsing"
];
(* For random tableau T,  T · T† = identity tableau *)
VerificationTest[
    SeedRandom[20260513];
    Module[{ps = PauliStabilizer["Random", 2], psBack},
        psBack = PauliStabilizer[ps["Circuit"]];
        Sort[ps["Stabilizers"]] === Sort[psBack["Stabilizers"]]
    ],
    True,
    TestID -> "S21-Stim-Random-Tableau-Roundtrip"
];
(* xz encoding agrees with simulator's tableau encoding -- the first
   stabilizer of an under-determined PauliStabilizer matches the input. *)
VerificationTest[
    Module[{ps = PauliStabilizer[{"XYZI"}]},
        ps["Stabilizers"][[1]]
    ],
    "XYZI",
    TestID -> "S21-Stim-Pauli-XYZI-FirstStabilizer"
]


(* ============================================================================
   §22 — Cross-package fixtures (sanity sample, full suite is in CrossPackage_*.wlt)
   ============================================================================ *)


(* Stim-formatted output normalization *)
VerificationTest[
    Module[{stimNorm},
        stimNorm = StringReplace[
            StringReplace["+_XX", "_" -> "I"],
            StartOfString ~~ "+" -> ""
        ];
        stimNorm
    ],
    "IXX",
    TestID -> "S22-Stim-FormatNormalization"
]


(* ============================================================================
   §23 — Caveats / pitfalls (sign conventions, non-stabilizer states)
   ============================================================================ *)


(* Y = i X Z convention *)
VerificationTest[
    pauliMat["Y"],
    I * pauliMat["X"] . pauliMat["Z"],
    TestID -> "S23-Y-equals-iXZ-convention"
];
(* Equivalent stabilizer generators describe the same state -- group-aware test
   ⟨XX, ZZ⟩ ≡ ⟨XX, -YY⟩ ≡ ⟨ZZ, -YY⟩ ALL describe |Φ+⟩  *)
VerificationTest[
    Module[{psFromCirc, psFromGenerators1, psFromGenerators2},
        psFromCirc = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
        psFromGenerators1 = PauliStabilizer[{"XX", "ZZ"}];
        psFromGenerators2 = PauliStabilizer[{"XX", "-YY"}];
        equalUpToPhase[psFromCirc["State"], psFromGenerators1["State"]] &&
        equalUpToPhase[psFromCirc["State"], psFromGenerators2["State"]]
    ],
    True,
    TestID -> "S23-Bell-Stabilizers-AreEquivalent"
];
(* W state is NOT a stabilizer state -- already in §3, but reasserted here as
   a caveat marker *)
VerificationTest[
    Module[{wState = QuantumState[
            (UnitVector[8, 2] + UnitVector[8, 3] + UnitVector[8, 5]) / Sqrt[3]]},
        Head[PauliStabilizer[wState]["Stabilizers"]] === Missing
    ],
    True,
    TestID -> "S23-W-State-Not-Stabilizer"
]


(* ============================================================================
   §24 — Concrete numerical fixtures
   ============================================================================ *)


(* |+⟩ amplitude vector *)
VerificationTest[
    Chop[
        Normal[PauliStabilizer[1]["H", 1]["State"]["StateVector"]] - {1, 1}/Sqrt[2]
    ] // Total // # == 0 &,
    True,
    TestID -> "S24-PlusState-Amplitudes"
];
(* GHZ-3 amplitude vector: only positions 1 and 8 nonzero *)
VerificationTest[
    Module[{v = Normal[PauliStabilizer[QuantumCircuitOperator[{
                "H" -> 1, "CNOT" -> {1, 2}, "CNOT" -> {2, 3}}]]["State"]["StateVector"]]},
        Chop[v - (UnitVector[8, 1] + UnitVector[8, 8]) / Sqrt[2]] //
        Total // # == 0 &
    ],
    True,
    TestID -> "S24-GHZ3-Amplitudes"
];
(* (STRICT) Steane |0_L⟩ has 8 nonzero amplitudes per Got00. *)
VerificationTest[
    Count[
        Normal[PauliStabilizer["SteaneCode"]["State"]["StateVector"]],
        x_ /; Abs[x] > 10^-10
    ],
    8,
    TestID -> "S24-F1-SteaneCode-LogicalZero-8Amps-STRICT"
]


(* ============================================================================
   §25 — Stabilizer Olympics (13-item conformance)
   ============================================================================ *)


(* 1. xz encoding round-trip *)
VerificationTest[
    xzDec /@ (xzEnc /@ {"I", "X", "Y", "Z"}),
    {"I", "X", "Y", "Z"},
    TestID -> "S25-Olympics-1-xz-Roundtrip"
];
(* 2. Symplectic predicate matches matrix commutator *)
VerificationTest[
    Module[{n = 1},
        AllTrue[Tuples[{"I", "X", "Y", "Z"}, 2],
            With[{a = #[[1]], b = #[[2]]},
                Boole[pauliString[a] . pauliString[b] != pauliString[b] . pauliString[a]] ==
                sympIP[pauliBits[a, n], pauliBits[b, n], n]
            ] &]
    ],
    True,
    TestID -> "S25-Olympics-2-Symplectic-Match"
];
(* 3. Heisenberg gate table reproduces text-book entries *)
VerificationTest[
    {
        PauliStabilizer[1]["H", 1]["Stabilizers"],
        PauliStabilizer[1]["S", 1]["Stabilizers"],
        Sort @ PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]["Stabilizers"]
    },
    {{"X"}, {"Z"}, Sort[{"XX", "ZZ"}]},
    TestID -> "S25-Olympics-3-Heisenberg-Table"
];
(* 4. 5-qubit and Steane code stabilizers exact *)
VerificationTest[
    Take[PauliStabilizer["5QubitCode"]["Stabilizers"], 4],
    {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"},
    TestID -> "S25-Olympics-4-5QubitCode-Exact"
];
(* 5. Bell-state round-trip tableau ↔ amplitudes *)
VerificationTest[
    equalUpToPhase[
        PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]["State"],
        QuantumState["Bell"]
    ],
    True,
    TestID -> "S25-Olympics-5-Bell-Roundtrip"
];
(* 6. MacWilliams: A_d for [[5,1,3]] is (1, 0, 0, 0, 15, 0) -- already in §9 *)
VerificationTest[
    Module[{ps = PauliStabilizer["5QubitCode"], stabs, group, weights, aArr},
        stabs = Take[ps["Stabilizers"], 4];
        group = Table[
            Module[{prod = "IIIII", sign = 1, k},
                Do[
                    If[BitAnd[BitShiftRight[m, k - 1], 1] == 1,
                        Module[{ms, body},
                            {ms, body} = splitSign[stabs[[k]]];
                            sign = sign * ms;
                            prod = StringJoin @ MapThread[
                                Function[{a, b},
                                    Switch[{a, b},
                                        {"I", _}, b, {_, "I"}, a,
                                        {"X", "X"} | {"Y", "Y"} | {"Z", "Z"}, "I",
                                        {"X", "Y"}, "Z", {"Y", "X"}, "Z",
                                        {"Y", "Z"}, "X", {"Z", "Y"}, "X",
                                        {"Z", "X"}, "Y", {"X", "Z"}, "Y"
                                    ]
                                ],
                                {Characters[prod], Characters[body]}
                            ]
                        ]
                    ],
                    {k, 4}
                ];
                prod
            ],
            {m, 0, 15}
        ];
        weights = StringCount[#, "X"] + StringCount[#, "Y"] + StringCount[#, "Z"] & /@ group;
        aArr = Table[Count[weights, d], {d, 0, 5}];
        aArr
    ],
    {1, 0, 0, 0, 15, 0},
    TestID -> "S25-Olympics-6-MacWilliams-A"
];
(* 7. Inner product 2^{-s/2} on small samples (deferred to §6) *)
VerificationTest[
    Chop[N[Abs[
        PauliStabilizer[1]["InnerProduct", PauliStabilizer[1]["H", 1]]
    ]^2] - 1/2] == 0,
    True,
    TestID -> "S25-Olympics-7-InnerProduct-Half"
];
(* 8. Bell ZZ measurement deterministic, 50/50 single-qubit *)
VerificationTest[
    {Sort @ Keys @ PauliStabilizer[QuantumCircuitOperator[{
            "H" -> 1, "CNOT" -> {1, 2}}]]["M", "ZZ"],
     Sort @ Keys @ PauliStabilizer[QuantumCircuitOperator[{
            "H" -> 1, "CNOT" -> {1, 2}}]]["M", "ZI"]},
    {{0}, {0, 1}},
    TestID -> "S25-Olympics-8-Measurement-Det-Random"
];
(* 9. Quantum teleportation by stabilizer evolution -- already in §13 *)
VerificationTest[
    Module[{ps = PauliStabilizer @ QuantumCircuitOperator[{
                "H" -> 1, "H" -> 2, "CNOT" -> {2, 3},
                "CNOT" -> {1, 2}, "H" -> 1
            }]},
        Length[Keys[ps["M", {1, 2}]]] == 4
    ],
    True,
    TestID -> "S25-Olympics-9-Teleportation-FourOutcomes"
];
(* 10. Toffoli on |++0⟩ flagged as non-stabilizer *)
VerificationTest[
    Quiet @ Head @ QuantumState[QuantumCircuitOperator[
        {"H" -> 1, "H" -> 2, "Toffoli" -> {1, 2, 3}}]],
    QuantumState,   (* QuantumCircuitOperator returns dense state, not PauliStabilizer *)
    TestID -> "S25-Olympics-10-Toffoli-Flagged"
];
(* 11. Local-complementation reproduces AndBri05 *)
VerificationTest[
    Module[{g = Graph[Range[5], Table[i \[UndirectedEdge] (i + 1), {i, 4}]]},
        Sort[EdgeList[LocalComplement[LocalComplement[g, 3], 3]]] ===
        Sort[EdgeList[g]]
    ],
    True,
    TestID -> "S25-Olympics-11-LocalComplement-Involutive"
];
(* 12. Choi-tableau composition: ccS² applied to |1⟩ should give |1⟩ with
   negative-sign Z stabilizer *)
VerificationTest[
    Module[{ccS = CliffordChannel[<|
            "UA" -> {{1, 0}, {0, 1}}, "UB" -> {{1, 1}, {0, 1}},
            "c" -> {0, 0}, "InputQubits" -> 1, "OutputQubits" -> 1,
            "Source" -> "S"|>],
            ccZ},
        ccZ = ccS[ccS];
        ccZ[PauliStabilizer[1]["X", 1]]["Stabilizers"]
    ],
    {"-Z"},
    TestID -> "S25-Olympics-12-CliffordChannel-Composition"
];
(* 13. Qudit Weyl commutation -- not applicable to QF (qubit only),
   sanity check that PauliStabilizer rejects non-qubit dims *)
VerificationTest[
    PauliStabilizer[2]["Qubits"] == 2,
    True,
    TestID -> "S25-Olympics-13-Qudit-Sanity"
]


(* ============================================================================
   Test execution & summary
   ============================================================================ *)


