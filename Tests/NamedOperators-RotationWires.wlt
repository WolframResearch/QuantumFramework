(* The rotation gates R_k(theta) = exp(-i theta/2 sigma_k) built in one step on the wire
   they are asked for, and a single-qudit gate on an integer wire n built there, as on
   {n}.

   References are independent of the route under test: MatrixExp of the Pauli matrix
   written out by hand for a qubit, and MatrixExp of the generalized Pauli matrix for a
   qutrit and a ququart.

   Regimes covered:
     general symbolic     R_X, R_Y, R_Z of a symbol theta, on wires 1 and 3
     exactly solvable     theta = Pi/7, Pi/2, 0; R_Z diagonal e^(-+i theta/2);
                          R(theta1) R(theta2) = R(theta1 + theta2)
     numerical reference  theta = 0.3 and 1.5 + 0.5 i against MatrixExp
     qudit                R_X, R_Y, R_Z of theta in dimensions 3 and 4
     wires                n and {n} agree for every single-qudit gate name; gates on
                          more than one wire keep their meaning for an integer n (the
                          first wire moved to n)
     scaling              a distinct symbolic rotation costs less than one general
                          2x2 operator *)

cm[qo_] := Normal[qo["Computational"]["Matrix"]]

zeroQ[m_] := AllTrue[Flatten[FullSimplify[m]], # === 0 &]

numZeroQ[m_] := Max[Abs[N[Flatten[m]]]] < 10^-12

$pauli = <|"RX" -> PauliMatrix[1], "RY" -> PauliMatrix[2], "RZ" -> PauliMatrix[3]|>;

$generator = <|"RX" -> "PauliX", "RY" -> "PauliY", "RZ" -> "PauliZ"|>;

KeyValueMap[
    Function[{gate, sigma},
        Scan[
            Function[angle,
                VerificationTest[
                    With[{d = cm[QuantumOperator[gate[angle], {3}]] - MatrixExp[- I angle / 2 sigma]},
                        If[MachineNumberQ[angle] || MachineNumberQ[Re[angle]], numZeroQ[d], zeroQ[d]]
                    ],
                    True,
                    TestID -> "RotationWires-" <> gate <> "-matrix-" <> ToString[angle, InputForm]
                ]
            ],
            {\[Theta], Pi / 7, Pi / 2, 0, 0.3, 1.5 + 0.5 I}
        ]
    ],
    $pauli
];

KeyValueMap[
    Function[{gate, sigma},
        VerificationTest[
            zeroQ[cm[QuantumOperator[gate[a]] @ QuantumOperator[gate[b]]] - MatrixExp[- I (a + b) / 2 sigma]],
            True,
            TestID -> "RotationWires-" <> gate <> "-group-law"
        ]
    ],
    $pauli
];

Scan[
    Function[{gd},
        (* the symbolic entries are compared at several angles: FullSimplify does not
           reduce every difference of the two closed forms to 0 in dimension 4 *)
        VerificationTest[
            With[{d = cm[QuantumOperator[gd[[1]][\[Theta], gd[[2]]], {2}]] -
                    MatrixExp[- I \[Theta] / 2 cm[QuantumOperator[$generator[gd[[1]]][gd[[2]]]]]]},
                AllTrue[{1 / 3, 7 / 5, -2, 11}, numZeroQ[d /. \[Theta] -> #] &]
            ],
            True,
            TestID -> "RotationWires-" <> gd[[1]] <> "-qudit-" <> ToString[gd[[2]]]
        ]
    ],
    Tuples[{Keys[$pauli], {3, 4}}]
];

(* a qudit rotation is stored as Simplify returns the operator exponential
   exp(-i theta/2 sigma_k), and is the same object whether it is built on {2}, placed on {2}
   after being built without a wire, or given as name -> 2 in a circuit. *)
Scan[
    Function[{gd},
        VerificationTest[
            With[{
                op = gd[[1]][\[Theta], gd[[2]]],
                exp = Simplify[QuantumOperator[Exp[- I \[Theta] / 2 QuantumOperator[$generator[gd[[1]]][gd[[2]]]]]]]
            },
                With[{qo = QuantumOperator[op, {2}]},
                    {
                        Normal[qo["State"]["StateVector"]] === Normal[exp["State"]["StateVector"]],
                        qo === QuantumOperator[QuantumOperator[op], {2}],
                        qo === QuantumCircuitOperator[{op -> 2}]["Elements"][[1]]
                    }
                ]
            ],
            {True, True, True},
            TestID -> "RotationWires-" <> gd[[1]] <> "-qudit-form-" <> ToString[gd[[2]]]
        ]
    ],
    Tuples[{Keys[$pauli], {3, 4}}]
];

VerificationTest[
    With[{qo = QuantumOperator["RZ"[\[Theta]], {3}]},
        {qo["Order"], qo["Label"], zeroQ[cm[qo] - DiagonalMatrix[{Exp[- I \[Theta] / 2], Exp[I \[Theta] / 2]}]]}
    ],
    {{{3}, {3}}, Subscript["R", "Z"][\[Theta]], True},
    TestID -> "RotationWires-RZ-order-label-diagonal"
]

(* every gate built on wire n by name is the gate built on {n} *)
Scan[
    Function[g,
        VerificationTest[
            With[{q = QuantumCircuitOperator[{g -> 3}]["Elements"][[1]], l = QuantumCircuitOperator[{g -> {3}}]["Elements"][[1]]},
                {q === l, q["Order"]}
            ],
            {True, {{3}, {3}}},
            TestID -> "RotationWires-integer-wire-" <> ToString[g, InputForm]
        ]
    ],
    {"X", "Y", "Z", "H", "S", "T", "P"[\[Phi]], "PhaseShift"[2], "RX"[\[Theta]], "RY"[0.3], "RZ"[Pi / 5],
        "RX"[\[Theta], 3], "X"[3], "H"[3], "P"[\[Phi], 3], "XRotation"[\[Theta]], "ZRotation"[\[Theta]]}
]

VerificationTest[
    {
        QuantumCircuitOperator[{"SWAP" -> 3}]["Elements"][[1]]["Order"],
        QuantumCircuitOperator[{"CNOT" -> 3}]["Elements"][[1]]["Order"],
        Head[FromOperatorShorthand["RandomHermitian" -> 3]]
    },
    {{{3, 2}, {3, 2}}, {{3, 2}, {3, 2}}, QuantumMeasurementOperator},
    TestID -> "RotationWires-integer-wire-multi-qudit-meaning"
]

(* A distinct symbolic rotation is one 2x2 matrix placed on one wire, so it costs less than
   a general 2x2 operator built through the generic constructor. Built through the operator
   exponential and moved to its wire it cost about twice as much. *)
VerificationTest[
    With[{k = 300},
        With[
            {
                general = Min @ Table[ClearSystemCache[]; First @ AbsoluteTiming[Table[QuantumOperator[{{Cos[Unique[]], 0}, {0, 1}}], k]], 3],
                rotation = Min @ Table[ClearSystemCache[]; First @ AbsoluteTiming[QuantumCircuitOperator[Table["RX"[Unique[]] -> Mod[j, 10] + 1, {j, k}]]], 3]
            },
            rotation / general < 1
        ]
    ],
    True,
    TestID -> "RotationWires-distinct-rotation-cost"
]
