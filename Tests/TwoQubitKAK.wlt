(* Tests for the native two-qubit Cartan (KAK) decomposition: qo["KAK"] and the
   qo["Decompose", gateset] umbrella. The decomposition factors any 2-qubit unitary as
   U = e^{i phi} (a1 (x) a2) . exp(i (c1 XX + c2 YY + c3 ZZ)) . (b1 (x) b2) and emits a
   QuantumCircuitOperator of single-qubit "U" gates and "CX" gates.

   Ordering follows the round-trip-first convention: exact matrix round-trips (including
   global phase) first, then unitarity, gate count, structured families, the umbrella, and
   the symbolic guard. Numeric tolerance is 10^-6: most gates round-trip exactly, but
   X-heavy factors lose ~10^-8 through ArcSin[Abs[b]] near theta = Pi in the shared
   UnitaryEulerAngles helper (a pre-existing limitation, also present in ZYZ). *)

(* ---- helpers (top-level, evaluated by TestReport before the tests use them) ---- *)
kakX = PauliMatrix[1]; kakY = PauliMatrix[2]; kakZ = PauliMatrix[3]; kakId = IdentityMatrix[2];
kakKP = KroneckerProduct;
kakCircuit[m_] := QuantumOperator[m, {1, 2}]["KAK"];
kakReconErr[m_] := Chop[Norm[Flatten[N @ (QuantumOperator @ kakCircuit[m])["MatrixRepresentation"] - N[m]], Infinity], 10^-6];
kakUnitErr[m_] := With[{mat = N @ (QuantumOperator @ kakCircuit[m])["MatrixRepresentation"]},
    Chop[Norm[Flatten[ConjugateTranspose[mat] . mat - IdentityMatrix[4]], Infinity], 10^-6]];
kakCX[m_] := Count[kakCircuit[m]["Elements"], _?(MatchQ[#["Label"], Subscript["C", ___][___]] &)];
kakGates = <|
    "Identity"  -> IdentityMatrix[4],
    "localX"    -> kakKP[kakX, kakId],
    "localY"    -> kakKP[kakId, kakY],
    "localHH"   -> kakKP[(kakX + kakZ) / Sqrt[2], (kakX + kakZ) / Sqrt[2]],
    "CNOT"      -> {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}},
    "CNOTrev"   -> {{1, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}, {0, 1, 0, 0}},
    "CZ"        -> DiagonalMatrix[{1, 1, 1, -1}],
    "CPhase"    -> DiagonalMatrix[{1, 1, 1, Exp[I 7/10]}],
    "SWAP"      -> {{1, 0, 0, 0}, {0, 0, 1, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}},
    "iSWAP"     -> {{1, 0, 0, 0}, {0, 0, I, 0}, {0, I, 0, 0}, {0, 0, 0, 1}},
    "DCX"       -> {{1, 0, 0, 0}, {0, 0, 0, 1}, {0, 1, 0, 0}, {0, 0, 1, 0}},
    "RXX"       -> MatrixExp[-I 1/2 kakKP[kakX, kakX]],
    "RYY"       -> MatrixExp[-I 9/10 kakKP[kakY, kakY]],
    "RZZ"       -> MatrixExp[-I 13/10 kakKP[kakZ, kakZ]],
    "Can"       -> MatrixExp[I (3/10 kakKP[kakX, kakX] + 6/10 kakKP[kakY, kakY] + 3/20 kakKP[kakZ, kakZ])],
    "Magic"     -> 1 / Sqrt[2] {{1, 0, 0, I}, {0, I, 1, 0}, {0, I, -1, 0}, {1, 0, 0, -I}}
|>;


BeginTestSection["TwoQubitKAK - round-trip (named gates, global phase included)"]

(* exact matrix round-trip for every named gate; e^{i phi} is recovered by projection,
   so global phase is reproduced, not just the local-equivalence class *)
VerificationTest[
    Max[N @ kakReconErr /@ Values[kakGates]],
    0,
    SameTest -> (Abs[#1 - #2] < 10^-6 &),
    TestID -> "KAK-roundtrip-named-gates"
]

(* the gates that are exactly reproduced (no X-heavy ArcSin loss) round-trip to machine zero *)
VerificationTest[
    kakReconErr /@ {kakGates["Identity"], kakGates["CZ"], kakGates["iSWAP"], kakGates["RXX"], kakGates["RZZ"], kakGates["Magic"]},
    {0, 0, 0, 0, 0, 0},
    TestID -> "KAK-roundtrip-exact-zero"
]

(* a state-preparation-style global-phase guard: e^{i 6/5} times a gate must round-trip *)
VerificationTest[
    kakReconErr[Exp[I 6/5] kakGates["CNOT"]],
    0,
    SameTest -> (Abs[#1 - #2] < 10^-6 &),
    TestID -> "KAK-roundtrip-global-phase"
]

(* separator-resonance guard: a gate in the c1 = ArcTan[Sqrt[2]]/2 class with nontrivial local
   factors merges two M2 eigenvalues under the default Sqrt[2] separator (reconErr ~ 0.085 without
   the guard); kakEigenbasis re-picks a separator so the decomposition stays exact *)
VerificationTest[
    BlockRandom[
        SeedRandom[5];
        With[{l = RandomVariate[CircularUnitaryMatrixDistribution[2], 4]},
            kakReconErr[
                kakKP[l[[1]], l[[2]]] .
                N @ MatrixExp[I (ArcTan[Sqrt[2]] / 2 kakKP[kakX, kakX] + 3/10 kakKP[kakY, kakY] + 1/10 kakKP[kakZ, kakZ])] .
                kakKP[l[[3]], l[[4]]]
            ]
        ]
    ],
    0,
    SameTest -> (Abs[#1 - #2] < 10^-6 &),
    TestID -> "KAK-separator-resonance-guard"
]

EndTestSection[]


BeginTestSection["TwoQubitKAK - random SU(4)"]

(* 50 Haar-random two-qubit unitaries must each round-trip *)
VerificationTest[
    SeedRandom[20260606];
    Max[Table[kakReconErr[RandomVariate[CircularUnitaryMatrixDistribution[4]]], {50}]],
    0,
    SameTest -> (Abs[#1 - #2] < 10^-6 &),
    TestID -> "KAK-roundtrip-random-50"
]

(* the compiled circuit operator is unitary for random input *)
VerificationTest[
    SeedRandom[1];
    Max[Table[kakUnitErr[RandomVariate[CircularUnitaryMatrixDistribution[4]]], {10}]],
    0,
    SameTest -> (Abs[#1 - #2] < 10^-6 &),
    TestID -> "KAK-unitary-random-10"
]

EndTestSection[]


BeginTestSection["TwoQubitKAK - structured families"]

(* 20 random diagonal unitaries (commuting, ZZ-class) round-trip exactly *)
VerificationTest[
    SeedRandom[3];
    Max[Table[kakReconErr[DiagonalMatrix[Exp[I RandomReal[{0, 2 Pi}, 4]]]], {20}]],
    0,
    SameTest -> (Abs[#1 - #2] < 10^-6 &),
    TestID -> "KAK-roundtrip-diagonal-20"
]

(* arbitrary-precision input is accepted and round-trips *)
VerificationTest[
    kakReconErr[N[MatrixExp[I (3/10 kakKP[kakX, kakX] + 6/10 kakKP[kakY, kakY] + 3/20 kakKP[kakZ, kakZ])], 30]],
    0,
    SameTest -> (Abs[#1 - #2] < 10^-6 &),
    TestID -> "KAK-roundtrip-precision30"
]

(* permutation gates *)
VerificationTest[
    Max[kakReconErr /@ {
        {{0, 1, 0, 0}, {1, 0, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}},
        {{0, 0, 1, 0}, {0, 1, 0, 0}, {1, 0, 0, 0}, {0, 0, 0, 1}}
    }],
    0,
    SameTest -> (Abs[#1 - #2] < 10^-6 &),
    TestID -> "KAK-roundtrip-permutations"
]

EndTestSection[]


BeginTestSection["TwoQubitKAK - gate count"]

(* a local A (x) B needs no entangling gate *)
VerificationTest[
    kakCX[kakGates["localHH"]],
    0,
    TestID -> "KAK-cx-local-zero"
]

(* CNOT-class gates use at most the optimal 3 CX (currently 2) *)
VerificationTest[
    kakCX[kakGates["CNOT"]] <= 3,
    True,
    TestID -> "KAK-cx-cnot-bound"
]

(* every 2-qubit gate is at most 6 CX with the current generic template (optimal target is 3) *)
VerificationTest[
    SeedRandom[5];
    Max[Table[kakCX[RandomVariate[CircularUnitaryMatrixDistribution[4]]], {10}]] <= 6,
    True,
    TestID -> "KAK-cx-generic-bound"
]

EndTestSection[]


BeginTestSection["TwoQubitKAK - Decompose umbrella"]

(* 1-qubit Decompose routes to ZYZ (which returns a bare QuantumOperator when the global phase
   is zero, or a circuit otherwise); check it reconstructs the gate either way *)
VerificationTest[
    With[{m = N @ {{0, 1}, {1, 0}}},
        Chop[Norm[Flatten[N @ (QuantumOperator @ QuantumOperator[m, {1}]["Decompose"])["MatrixRepresentation"] - m], Infinity], 10^-6]
    ],
    0,
    SameTest -> (Abs[#1 - #2] < 10^-6 &),
    TestID -> "KAK-decompose-1q-routes-zyz"
]

(* 2-qubit Decompose equals KAK *)
VerificationTest[
    SeedRandom[7];
    With[{m = RandomVariate[CircularUnitaryMatrixDistribution[4]]},
        Chop[Norm[Flatten[
            N @ (QuantumOperator @ QuantumOperator[m, {1, 2}]["Decompose"])["MatrixRepresentation"] -
            N @ (QuantumOperator @ QuantumOperator[m, {1, 2}]["KAK"])["MatrixRepresentation"]
        ], Infinity], 10^-6]
    ],
    0,
    SameTest -> (Abs[#1 - #2] < 10^-6 &),
    TestID -> "KAK-decompose-2q-equals-kak"
]

(* "KAK" and "Decompose" are advertised in the property list *)
VerificationTest[
    ContainsAll[QuantumOperator["X"]["Properties"], {"KAK", "Decompose"}],
    True,
    TestID -> "KAK-properties-listed"
]

EndTestSection[]


BeginTestSection["TwoQubitKAK - symbolic guard"]

(* a symbolic / parametric operator is not decomposed; it returns the operator undecomposed
   (a one-gate circuit) with a message, and does not hang or recurse *)
VerificationTest[
    Head @ QuantumOperator[MatrixExp[-I \[Theta] / 2 kakKP[kakX, kakX]], {1, 2}]["KAK"],
    QuantumCircuitOperator,
    {QuantumCircuitOperator::kaksymbolic},
    TestID -> "KAK-symbolic-guard"
]

EndTestSection[]
