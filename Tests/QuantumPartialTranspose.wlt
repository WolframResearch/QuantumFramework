BeginTestSection["QuantumPartialTranspose - computational frame"]

(* Anchors the partial-transpose convention: swap the ket and bra legs of the
   chosen qudit. *)
VerificationTest[
    Normal[QuantumPartialTranspose[QuantumState["Bell"], {2}]["DensityMatrix"]],
    {{1, 0, 0, 0}, {0, 0, 1, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}}/2,
    TestID -> "Bell-computational-T2"
]

EndTestSection[]


BeginTestSection["QuantumPartialTranspose - complex-element bases"]

(* comp(rho^T2) must equal the computational partial transpose of comp(rho) in
   any frame. Complex frames (PauliY, Fourier) are the discriminating cases: a
   real frame cannot see a dropped conjugation, and PT spectra are invariant
   under local frame changes, so only matrix-element comparisons discriminate. *)

VerificationTest[
    With[{qs = QuantumState[Normalize[{1, 2, 3 I, 4}], QuantumBasis[{"PauliY", "PauliY"}]]},
        With[{
            rho = Normal[QuantumState[qs, QuantumBasis[4]]["DensityMatrix"]],
            pt = Normal[QuantumState[QuantumPartialTranspose[qs, {2}], QuantumBasis[4]]["DensityMatrix"]]
        },
            Simplify[pt - ArrayReshape[Transpose[ArrayReshape[rho, {2, 2, 2, 2}], {1, 4, 3, 2}], {4, 4}]]
        ]
    ],
    ConstantArray[0, {4, 4}],
    TestID -> "PauliY-product-T2-frame-faithful"
]

(* General symbolic amplitudes (real parameters, complex structure from the
   frame and the explicit I): the identity holds as a closed form. *)
VerificationTest[
    With[{qs = QuantumState[{a, b I, c, d}, QuantumBasis[{"PauliY", "PauliY"}]]},
        With[{
            rho = Normal[QuantumState[qs, QuantumBasis[4]]["DensityMatrix"]],
            pt = Normal[QuantumState[QuantumPartialTranspose[qs, {2}], QuantumBasis[4]]["DensityMatrix"]]
        },
            Simplify[
                pt - ArrayReshape[Transpose[ArrayReshape[rho, {2, 2, 2, 2}], {1, 4, 3, 2}], {4, 4}],
                Assumptions -> (a | b | c | d) \[Element] Reals
            ]
        ]
    ],
    ConstantArray[0, {4, 4}],
    TestID -> "PauliY-symbolic-T2-frame-faithful"
]

(* Mixed dimensions: qutrit Fourier frame times qubit PauliY frame, transposing
   the qutrit leg. *)
VerificationTest[
    With[{qs = QuantumState[Normalize[{1, I, 2, 3 I, 4, 5}], QuantumBasis[{"Fourier"[3], "PauliY"}]]},
        With[{
            rho = Normal[QuantumState[qs, QuantumBasis[6]]["DensityMatrix"]],
            pt = Normal[QuantumState[QuantumPartialTranspose[qs, {1}], QuantumBasis[6]]["DensityMatrix"]]
        },
            Simplify[pt - ArrayReshape[Transpose[ArrayReshape[rho, {3, 2, 3, 2}], {3, 2, 1, 4}], {6, 6}]]
        ]
    ],
    ConstantArray[0, {6, 6}],
    TestID -> "Fourier3-PauliY-T1-frame-faithful"
]

(* Tripartite, transposing the two non-adjacent complex legs. *)
VerificationTest[
    With[{qs = QuantumState[Normalize[{1, I, 2, 3, 5 I, 8, 13, 21 I}], QuantumBasis[{"PauliY", "Fourier"[2], "PauliY"}]]},
        With[{
            rho = Normal[QuantumState[qs, QuantumBasis[8]]["DensityMatrix"]],
            pt = Normal[QuantumState[QuantumPartialTranspose[qs, {1, 3}], QuantumBasis[8]]["DensityMatrix"]]
        },
            Simplify[pt - ArrayReshape[Transpose[ArrayReshape[rho, ConstantArray[2, 6]], Cycles[{{1, 4}, {3, 6}}]], {8, 8}]]
        ]
    ],
    ConstantArray[0, {8, 8}],
    TestID -> "Tripartite-T13-frame-faithful"
]

(* A state carrying input legs: transposing an input qudit conjugates the input
   side of the basis and is frame-faithful through the matrix conversion route
   (both sides of the comparison go through it, matrix-form state in, matrix-form
   state out). *)
VerificationTest[
    With[{qs = QuantumOperator[QuantumOperator["CNOT"], QuantumBasis[{"PauliY", "PauliY"}]]["State"]["MatrixState"]},
        With[{rho = Normal[qs["Computational"]["DensityMatrix"]]},
            Simplify[
                Normal[qs["Transpose", {3}]["Computational"]["DensityMatrix"]] -
                ArrayReshape[Transpose[ArrayReshape[rho, ConstantArray[2, 8]], Cycles[{{3, 7}}]], {16, 16}]
            ]
        ]
    ],
    ConstantArray[0, {16, 16}],
    TestID -> "InputLeg-T3-matrix-route-frame-faithful"
]

(* T2 is an involution in any frame. *)
VerificationTest[
    With[{qs = QuantumState[Normalize[{1, 2, 3 I, 4}], QuantumBasis[{"PauliY", "PauliY"}]]},
        Simplify[
            Normal[QuantumState[QuantumPartialTranspose[QuantumPartialTranspose[qs, {2}], {2}], QuantumBasis[4]]["DensityMatrix"]] -
            Normal[QuantumState[qs, QuantumBasis[4]]["DensityMatrix"]]
        ]
    ],
    ConstantArray[0, {4, 4}],
    TestID -> "PauliY-T2-involution"
]

(* The partial transpose of a Hermitian state stays Hermitian and unit-trace in
   the computational frame. *)
VerificationTest[
    With[{pt = Normal[QuantumState[
            QuantumPartialTranspose[QuantumState[Normalize[{1, 2, 3 I, 4}], QuantumBasis[{"PauliY", "PauliY"}]], {2}],
            QuantumBasis[4]]["DensityMatrix"]]},
        {Simplify[pt - ConjugateTranspose[pt]], Simplify[Tr[pt]]}
    ],
    {ConstantArray[0, {4, 4}], 1},
    TestID -> "PauliY-T2-hermitian-unit-trace"
]

EndTestSection[]
