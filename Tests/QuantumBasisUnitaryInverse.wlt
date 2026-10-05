(* A basis matrix that is exactly unitary inverts by its conjugate transpose, B^-1 = B^dagger.
   QF's Fourier basis matrix is the exact FourierMatrix[d], whose entries are roots of unity,
   and the general exact Inverse of it does not finish in two minutes at d = 16. Each change of
   basis below is checked against F^dagger written out from FourierMatrix: a pure state goes to
   F^dagger v, a mixed state to F^dagger rho F, and an operator M whose input basis is the
   Fourier basis acts as M F^dagger. A matrix that is not exactly unitary, a machine matrix and
   a symbolic matrix keep the general Inverse. *)

BeginTestSection["QuantumBasis - change of basis with a unitary basis matrix"]

uinvF = FourierMatrix[16];
uinvV = Array[uinvA, 16];
uinvRho = DiagonalMatrix[Range[16]] / 136 + SparseArray[{{1, 2} -> 1 / 100, {2, 1} -> 1 / 100}, {16, 16}];
uinvClose[a_, b_] := Max[Abs[N[Normal[a]] - N[Normal[b]]]] < 10 ^ -12;

VerificationTest[
    Simplify[Normal[QuantumState[QuantumState[uinvV, 16], QuantumBasis["Fourier"[16]]]["StateVector"]] - ConjugateTranspose[uinvF] . uinvV],
    ConstantArray[0, 16],
    TestID -> "UnitaryInverse-Fourier16-pure-state",
    TimeConstraint -> 10
]

VerificationTest[
    Simplify[Normal[QuantumState[QuantumState[QuantumState[uinvV, 16], QuantumBasis["Fourier"[16]]], QuantumBasis[16]]["StateVector"]] - uinvV],
    ConstantArray[0, 16],
    TestID -> "UnitaryInverse-Fourier16-round-trip",
    TimeConstraint -> 10
]

VerificationTest[
    uinvClose[QuantumState[QuantumState[uinvRho, 16], QuantumBasis["Fourier"[16]]]["DensityMatrix"], ConjugateTranspose[uinvF] . uinvRho . uinvF],
    True,
    TestID -> "UnitaryInverse-Fourier16-mixed-state",
    TimeConstraint -> 10
]

VerificationTest[
    With[{qb = QuditBasis["Fourier"[16]]}, uinvClose[qb["Inverse"]["ReducedMatrix"] . qb["ReducedMatrix"], IdentityMatrix[16]]],
    True,
    TestID -> "UnitaryInverse-Fourier16-QuditBasis-Inverse",
    TimeConstraint -> 10
]

VerificationTest[
    With[{m = Table[Mod[3 i + 5 j, 7] - 3, {i, 16}, {j, 16}], w = Range[16] - 8},
        uinvClose[
            QuantumCircuitOperator[{QuantumOperator[m, {1}, QuantumBasis[QuditBasis[16], QuditBasis["Fourier"[16]]]]}][QuantumState[w, 16]]["StateVector"],
            m . ConjugateTranspose[uinvF] . w
        ]
    ],
    True,
    TestID -> "UnitaryInverse-Fourier16-input-basis-in-a-circuit",
    TimeConstraint -> 10
]

VerificationTest[
    With[{m = SparseArray[QuditBasis["Fourier"[16]]["ReducedMatrix"]]}, MatchQ[MatrixInverse[m], _SparseArray]],
    True,
    TestID -> "UnitaryInverse-sparse-matrix-stays-sparse",
    TimeConstraint -> 10
]

VerificationTest[
    With[{m = QuditBasis["Tetrahedron"]["ReducedMatrix"]}, MatrixInverse[m] === Inverse[m]],
    True,
    TestID -> "UnitaryInverse-frame-basis-keeps-Inverse"
]

VerificationTest[
    With[{m = N[FourierMatrix[16]]}, MatrixInverse[m] === Inverse[m]],
    True,
    TestID -> "UnitaryInverse-machine-matrix-keeps-Inverse"
]

VerificationTest[
    With[{m = {{Cos[uinvT], -Sin[uinvT]}, {Sin[uinvT], Cos[uinvT]}}}, MatrixInverse[m] === Inverse[m]],
    True,
    TestID -> "UnitaryInverse-symbolic-matrix-keeps-Inverse"
]

EndTestSection[]
