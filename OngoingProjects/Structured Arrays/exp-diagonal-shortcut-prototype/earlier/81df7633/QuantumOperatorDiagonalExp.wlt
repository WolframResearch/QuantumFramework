(* The exponential of a diagonal operator is the exponential of each diagonal entry.
   The references are built without the route under test: Exp of the list of
   diagonal entries, a closed form, or MatrixExp on the matrix itself.

   Regimes covered:
     general symbolic     exp of diag(a, b, c, d); a declared parameter gamma in
                          exp(-i gamma H) for a 12-qubit Ising H
     exactly solvable     exp of an exact diagonal (i pi, log 2, 0); dephasing of a
                          3-qubit register, whose coherences decay in closed form
     limiting             16-qubit Ising H; gamma = 0 gives the identity
     numerical reference  MatrixExp of the same matrix, for a coupling at roundoff
                          next to a large entry, where the exponential is not
                          diagonal
     failure / edge       underflow of e^(-2000); an off-diagonal entry that is zero
                          only after simplification; collision of two parameters *)

BeginTestSection["QuantumOperator - exponential of a diagonal operator"]

deIsing[n_] := With[{
    couplings = RandomReal[{-1, 1}, {n, n}],
    spins = 1 - 2 Tuples[{0, 1}, n]
},
    Total[(spins . UpperTriangularize[couplings, 1]) spins, {2}]
]

deH16 = BlockRandom[deIsing[16], RandomSeeding -> 7];

(* exp(-i 0.3 H) for a 16-qubit Ising H is one exponential per basis state; its
   diagonal is Exp of the list, to the last bit. *)
VerificationTest[
    With[{u = Exp[-I 0.3 QuantumOperator[SparseArray[Band[{1, 1}] -> deH16], Range[16]]]},
        Max[Abs[Normal[Diagonal[u["Matrix"]]] - Exp[-I 0.3 deH16]]] == 0
    ],
    True,
    TimeConstraint -> 30,
    TestID -> "DiagonalExp-16-qubit-Ising"
]

(* The same exponential with the time as a declared parameter, read at 25 values. *)
VerificationTest[
    With[{
        h = BlockRandom[deIsing[12], RandomSeeding -> 5],
        pts = Subdivide[0.05, 1.25, 24]
    },
        With[{u = Exp[QuantumOperator[SparseArray[Band[{1, 1}] -> -I deGamma h], Range[12], "Parameters" -> {deGamma}]]},
            {
                Map[Max[Abs[Normal[Diagonal[u[#]["Matrix"]]] - Exp[-I # h]]] == 0 &, pts],
                Normal[u[0]["Matrix"]] == IdentityMatrix[4096]
            }
        ]
    ],
    {ConstantArray[True, 25], True},
    TimeConstraint -> 30,
    TestID -> "DiagonalExp-parameter-scan-12-qubit"
]

(* Symbolic and exact diagonals keep their closed form. *)
VerificationTest[
    {
        Normal[Exp[QuantumOperator[DiagonalMatrix[{deA, deB, deC, deD}], {1, 2}]]["Matrix"]],
        Normal[Exp[QuantumOperator[DiagonalMatrix[{I Pi, Log[2], 0, -I Pi / 2}], {1, 2}]]["Matrix"]]
    },
    {DiagonalMatrix[Exp[{deA, deB, deC, deD}]], DiagonalMatrix[{-1, 2, 1, -I}]},
    TestID -> "DiagonalExp-symbolic-and-exact"
]

(* An off-diagonal entry that simplifies to zero leaves a diagonal exponential, with
   no pole where the two diagonal entries coincide. *)
VerificationTest[
    Normal[Exp[QuantumOperator[{{deA, Sin[deY]^2 + Cos[deY]^2 - 1}, {0, deB}}]]["Matrix"]] /. deB -> deA,
    {{E^deA, 0}, {0, E^deA}},
    TestID -> "DiagonalExp-hidden-zero-no-pole"
]

(* e^(-2000) underflows to zero, and e^0 stays exactly the machine number 1. *)
VerificationTest[
    With[{m = Normal[Exp[-2000. QuantumOperator[DiagonalMatrix[{1., 0.}]]]["Matrix"]]},
        {m[[1, 1]] == 0, m[[2, 2]] - 1 == 0, Precision[m[[2, 2]]]}
    ],
    {True, True, MachinePrecision},
    TestID -> "DiagonalExp-underflow"
]

(* A machine diagonal with exact entries among them, integers or constants such as
   Pi, is exponentiated in machine numbers, and an underflowing entry gives zero
   without a message. *)
VerificationTest[
    With[{d = Normal[Diagonal[Exp[QuantumOperator[#]]["Matrix"]]] & /@ {{{1., 0}, {0, 2}}, {{-800., 0}, {0, 0}}, {{1., 0}, {0, Pi}}, {{-800., 0}, {0, Pi}}}},
        {Precision /@ d, Max[MapThread[Abs[#1 - #2] / Max[1, Abs[#2]] &, {Flatten[d], N[{E, E^2, 0, 1, E, E^Pi, 0, E^Pi}]}]] < 10^-15}
    ],
    {ConstantArray[MachinePrecision, 4], True},
    TestID -> "DiagonalExp-mixed-exact-and-machine"
]

(* A matrix-type object with no qudits keeps its 1 x 1 matrix. *)
VerificationTest[
    Normal[QuantumPartialTrace[QuantumState[{{3/10, 0}, {0, 7/10}}]]["Matrix"]],
    {{1}},
    TestID -> "DiagonalExp-rank-zero-matrix"
]

(* Pure dephasing of a 3-qubit register: every coherence rho_xy picks up the phase
   e^(-i (h_x - h_y) t) and decays at 2 gamma per differing qubit, and the
   populations do not move. *)
VerificationTest[
    With[{
        h = BlockRandom[deIsing[3], RandomSeeding -> 5],
        rho = BlockRandom[QuantumState["RandomMixed"[3]], RandomSeeding -> 3],
        t = 0.7, gamma = 0.2
    },
        With[{
            out = Exp[t QuantumOperator["Liouvillian"[
                QuantumOperator[SparseArray[Band[{1, 1}] -> h], Range[3]],
                Table[QuantumOperator[SparseArray[KroneckerProduct @@ ReplacePart[ConstantArray[IdentityMatrix[2], 3], k -> {{1, 0}, {0, -1}}]], Range[3]], {k, 3}],
                ConstantArray[gamma, 3]
            ]]][rho]["DensityMatrix"],
            r = Normal[rho["DensityMatrix"]]
        },
            Max[Abs[Normal[out] - Table[
                r[[x, y]] Exp[-I (h[[x]] - h[[y]]) t - 2 gamma HammingDistance[IntegerDigits[x - 1, 2, 3], IntegerDigits[y - 1, 2, 3]] t],
                {x, 8}, {y, 8}
            ]]] < 10^-14
        ]
    ],
    True,
    TestID -> "DiagonalExp-dephasing-closed-form"
]

(* A coupling of 10^-9 between two degenerate levels, beside an entry 10^6, is at
   roundoff relative to the largest entry, yet over t = 10^6 it rotates the pair by
   10^-3. The exponential must not treat the matrix as diagonal. *)
VerificationTest[
    With[{m = {{1.*^6, 0., 0.}, {0., 0., 1.*^-9}, {0., 1.*^-9, 0.}}, t = 1.*^6},
        Abs[Abs[Normal[Exp[-I t QuantumOperator[m, 3]]["Matrix"]][[2, 3]]] - Sin[10^-3]] < 10^-4
    ],
    True,
    TestID -> "DiagonalExp-resonant-roundoff-coupling-kept"
]

(* Two parameters that coincide: the exponential of diag(a, b) at a = b = 1 is e I. *)
VerificationTest[
    Normal[Exp[QuantumOperator[DiagonalMatrix[{deA, deB}], "Parameters" -> {deA, deB}]][1, 1]["Matrix"]],
    E IdentityMatrix[2],
    TestID -> "DiagonalExp-parameter-collision"
]

EndTestSection[]
