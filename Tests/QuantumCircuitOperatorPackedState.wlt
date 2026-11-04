(* A circuit applied to a state contracts a tensor network whose first tensor is the state.
   A full state of machine numbers enters that network as a packed array, 16 bytes per
   amplitude, instead of a rank-n SparseArray, which stores n positions with every
   amplitude. Exact, symbolic and sparse states enter as before.

   References are built by hand: the circuit's unitary as a product of Kronecker products of
   2x2 and 4x4 gate matrices, wire 1 the most significant, applied to the amplitudes.

   Regimes covered:
     numerical reference  random machine states on 6 and 10 qubits under layers of R_X and
                          CNOT, against the unitary built by hand
     exact                H, H, CNOT on |00> stays exact: (1, 1, 1, 1)/2
     symbolic             R_X(theta) on a machine state keeps theta, and agrees with the
                          reference at a value of theta
     sparse               a machine GHZ state on 8 qubits, sparse enough to stay sparse
     mixed                a machine density matrix, against U rho U^dagger
     memory               a full 18-qubit state is applied in a few times its own size *)

rxMatrix[t_] := {{Cos[t / 2], -I Sin[t / 2]}, {-I Sin[t / 2], Cos[t / 2]}}

cnotMatrix = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}};

onWire[g_, q_, n_] := KroneckerProduct[IdentityMatrix[2 ^ (q - 1)], g, IdentityMatrix[2 ^ (n - q - Log2[Length[g]] + 1)]]

(* the gates of one layer: R_X(angle q) on each wire q, then CNOT on (1, 2), (3, 4), ... *)
layerSpec[n_, angle_] := Join[Table["RX"[angle[q]] -> q, {q, n}], Table["CNOT" -> {q, q + 1}, {q, 1, n - 1, 2}]]

layerMatrix[n_, angle_] := Dot @@ Reverse @ Join[
    Table[onWire[rxMatrix[angle[q]], q, n], {q, n}],
    Table[onWire[cnotMatrix, q, n], {q, 1, n - 1, 2}]
]

randomState[n_] := (SeedRandom[n]; Normalize[RandomComplex[{-1 - I, 1 + I}, 2 ^ n]])

Scan[
    Function[n,
        VerificationTest[
            With[{
                v = randomState[n],
                angles = {0.3 # + 0.1 &, 0.7 - 0.05 # &}
            },
                Max @ Abs[
                    Normal[QuantumCircuitOperator[Catenate[layerSpec[n, #] & /@ angles]][QuantumState[v]]["StateVector"]] -
                        Fold[layerMatrix[n, #2] . #1 &, v, angles]
                ]
            ],
            _ ? (# < 10^-12 &),
            SameTest -> MatchQ,
            TestID -> "PackedState-random-machine-" <> ToString[n]
        ]
    ],
    {6, 10}
]

VerificationTest[
    Normal[QuantumCircuitOperator[{"H" -> 1, "H" -> 2, "CNOT" -> {1, 2}}][QuantumState["00"]]["StateVector"]],
    {1, 1, 1, 1} / 2,
    TestID -> "PackedState-exact-state-stays-exact"
]

VerificationTest[
    With[{v = randomState[2]},
        With[{out = Normal[QuantumCircuitOperator[{"RX"[\[Theta]] -> 1}][QuantumState[v]]["StateVector"]]},
            {! FreeQ[out, \[Theta]], Max @ Abs[(out - onWire[rxMatrix[\[Theta]], 1, 2] . v) /. \[Theta] -> 0.7]}
        ]
    ],
    {True, _ ? (# < 10^-12 &)},
    SameTest -> MatchQ,
    TestID -> "PackedState-symbolic-gate-on-machine-state"
]

VerificationTest[
    With[{v = N @ Normal @ SparseArray[{1 -> 1 / Sqrt[2], 256 -> 1 / Sqrt[2]}, 256]},
        Max @ Abs[
            Normal[QuantumCircuitOperator[layerSpec[8, 0.2 # &]][QuantumState[v]]["StateVector"]] - layerMatrix[8, 0.2 # &] . v
        ]
    ],
    _ ? (# < 10^-12 &),
    SameTest -> MatchQ,
    TestID -> "PackedState-sparse-GHZ"
]

VerificationTest[
    With[{rho = (SeedRandom[3]; With[{a = RandomComplex[{-1 - I, 1 + I}, {4, 4}]}, # / Tr[#] & [a . ConjugateTranspose[a]]]), u = layerMatrix[2, 0.4 &]},
        Max @ Abs[Normal[QuantumCircuitOperator[layerSpec[2, 0.4 &]][QuantumState[rho]]["DensityMatrix"]] - u . rho . ConjugateTranspose[u]]
    ],
    _ ? (# < 10^-12 &),
    SameTest -> MatchQ,
    TestID -> "PackedState-density-matrix"
]

(* A rank-18 SparseArray of the state alone is about 8 times the state, and applying the
   circuit through it peaks near 10 times; through a packed tensor it peaks near 4 times. *)
VerificationTest[
    With[{qs = QuantumState[randomState[18]], qc = QuantumCircuitOperator[layerSpec[18, 0.1 # &]]},
        MaxMemoryUsed[qc[qs]] / ByteCount[qs] < 7
    ],
    True,
    TestID -> "PackedState-memory"
]
