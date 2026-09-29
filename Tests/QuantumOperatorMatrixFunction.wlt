(* A numeric function applied to an operator, f[qo], is the matrix function of the
   stored matrix. Every reference below is built independently of the route under
   test: either f of a spectrum whose eigenbasis is known by construction, or a
   closed form. The numeric operators are large enough, and their spectra clustered
   enough, that a route interpolating f by a polynomial in the matrix loses all
   accuracy: the Ising spectra here have 32 distinct eigenvalues, each twice
   degenerate under the global spin flip.

   Physics. f(H) = sum over eigenspaces of f(lambda) P_lambda for normal H, and the
   Jordan form for a defective one; for an operator with parameters f is applied to
   the matrix with the values in place, which is regular where two eigenvalues
   collide. Invariants checked: f(H) Hermitian for real f and Hermitian H,
   exp(-i t H) unitary, cos^2 + sin^2 = 1, [f(H), H] = 0, sqrt(rho) positive,
   exp(A (x) 1 + 1 (x) B) = exp(A) (x) exp(B), exp(-i t H) = 1 - i t H + O(t^2).
   Refused base case: a single qubit at generic parameter values.

   Regimes covered:
     general symbolic     cos(t X), sin(t ZZ), t(XX + YY), [[a, b], [b, -a]],
                          cos(t sum_i Z_i) on n = 2 .. 6 qubits, cos(t H) for the
                          4-qubit Heisenberg chain, exp(-i H) of the 3-qubit
                          transverse-field Ising chain in (J, g), "Parameters",
                          exact input
     exactly solvable     Jordan block, exact Log, Hadamard-rotated Ising, |A| of a
                          diagonalizable non-normal A, sqrt of a projector
     limiting             [[1, 1], [0, 1 + d]] for d = 10^-1 .. 10^-15 (separated to
                          defective), eigenvalue 10^-12 next to 1, coupling 10^-9,
                          decoupled subsystems, short time t -> 0
     numerical reference  64-dim Ising spectrum in a random eigenbasis, 600 random
                          Hermitian draws against a 30-digit diagonalization,
                          cos^2 + sin^2 = 1 and [cos H, H] = 0
     failure / edge       Log of a singular operator (numeric, diagonal, exact,
                          defective, or reached by substitution) is a Failure;
                          degenerate spectra of multiplicity 8 and 4; parameter
                          values at and next to a collision of two eigenvalues

   The zero base, 0^M, the limit of b^M as b -> 0 (the projector onto the null space
   of M along its range), has its own rows. Refused base case: one diagonal mode with
   a rank-1 kernel.
     general symbolic     0^(w n), 0^(w M) for every w > 0, 0^(c M) = 0^M for every
                          c > 0, a parametric Jordan block, the XXZ pair on both
                          sides of its level crossing in d, Lindblad generators for
                          every rate: damping, the driven qubit (resonance
                          fluorescence, through its exceptional point) and the
                          thermal generator, the ring of 4 sites in a field h for
                          every h > 0
     exactly solvable     vacuum projector of one and three modes, ground spaces of
                          the XX chain, the Heisenberg triangle and the
                          ferromagnetic ring for n = 4, 6, 8, the entangled kernel
                          of B^dagger B,
                          steady states of damping, dephasing, collective decay,
                          the driven qubit and the thermal generator (ThermalState),
                          Louisell's identity, the total-loss channel,
                          ThermalState at nbar = 0
     limiting             thermal state to first order in q, exp(-beta (H - E0))
                          -> 0^(H - E0) with leading term exp(-beta gap) times the
                          first-excited projector, steady state as the
                          t -> Infinity limit of exp(t L), {{e, 1}, {0, 0}} as
                          e -> 0, where the limits in b and e do not commute, and
                          the ferromagnetic ring in a field h, of rank 1 for every
                          h > 0 and n + 1 at h = 0, where those in b and h do not
     numerical reference  exact against machine and 40 digits, a machine number
                          operator, diagonal and in a rotated basis, a rotated number
                          operator of dimension 128, projectors built by
                          construction, P.P = P, M.P = P.M = 0, Tr P = nullity,
                          Hermitian P for Hermitian M, Choi positivity and trace
                          preservation for steady-state projectors
     failure / edge       eigenvalue -1, a unitary generator (imaginary spectrum),
                          a dark coherence that rotates, Jordan blocks at zero
                          (exact, machine, conjugated, symbolic), machine matrices
                          within roundoff of one, exact spectra over 150 orders of
                          magnitude and an eigenvalue -10^-130, spectra {I, 1} and
                          {-10^-40, 1} in a basis of condition 10^40 and a Jordan
                          block at 10^-16 beside the null space, entries whose
                          value cancels over hundreds of digits, Hermitian or not,
                          complex values written without I, as ArcCos[2], and
                          precisions past 307 digits, all decided exactly, and a
                          zero and a real part that are 0 but not written as 0,
                          alone, inside an entry or on a branch cut, and huge
                          entries. The exact decisions end where a number that
                          Element does not prove algebraic lies within 10^-3200 of
                          zero, or where telling it from zero takes more working
                          precision than N is given (a cancellation between terms
                          far beyond 10^10000 in size): it is read as zero, and in
                          the second case N says so. Past $MinNumber or $MaxNumber
                          N reports an underflow or an overflow (with N::meprec): a
                          value that underflows is read as zero, and one that
                          overflows as nonzero with the sign Sign reads from the
                          expression, Indeterminate where Sign cannot. The exact
                          decisions hold between those two *)

Needs["Wolfram`QuantumFramework`SecondQuantization`"]

BeginTestSection["QuantumOperator - matrix functions"]

mfIsing[n_] := With[{
    couplings = RandomReal[{-1, 1}, {n, n}],
    spins = 1 - 2 Tuples[{0, 1}, n]
},
    Total[(spins . UpperTriangularize[couplings, 1]) spins, {2}]
]

mfDistance[a_, b_] := Max[Abs[Normal[a] - Normal[b]]]

mfSpectrum = BlockRandom[mfIsing[6], RandomSeeding -> 7];

mfUnitary = BlockRandom[Orthogonalize[RandomComplex[{-1 - I, 1 + I}, {64, 64}]], RandomSeeding -> 11];

(* The Hermitian matrix with spectrum spec whose eigenvectors are the columns of u^†. *)
mfRotated[u_, spec_] := ConjugateTranspose[u] . DiagonalMatrix[spec] . u

(* Six qubits, diagonal in the computational basis: f acts entry by entry. *)
mfIsingDiagonal = QuantumOperator[SparseArray[Band[{1, 1}] -> mfSpectrum], Range[6]];

VerificationTest[
    Map[
        mfDistance[#[mfIsingDiagonal]["Matrix"], DiagonalMatrix[# /@ mfSpectrum]] &,
        {Cos, Sin, Sqrt[# + 10] &, Log[# + 10] &}
    ],
    {_ ? (# < 10^-12 &) ..},
    SameTest -> MatchQ,
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-Ising64-diagonal"
]

(* The same spectrum in the Hadamard-rotated basis, H^n diag(h) H^n: a dense real
   symmetric matrix whose eigenbasis is known exactly. *)
mfHadamard = Normal[KroneckerProduct @@ ConstantArray[{{1., 1.}, {1., -1.}} / Sqrt[2], 6]];

VerificationTest[
    With[{qo = QuantumOperator[mfHadamard . DiagonalMatrix[mfSpectrum] . mfHadamard, Range[6]]},
        Map[
            mfDistance[#[qo]["Matrix"], mfHadamard . DiagonalMatrix[# /@ mfSpectrum] . mfHadamard] &,
            {Cos, Sin}
        ]
    ],
    {_ ? (# < 10^-12 &) ..},
    SameTest -> MatchQ,
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-Ising64-hadamard-basis"
]

(* A complex Hermitian operator with the Ising spectrum in a random eigenbasis. *)
VerificationTest[
    With[{qo = QuantumOperator[mfRotated[mfUnitary, mfSpectrum], Range[6]]},
        mfDistance[Cos[qo]["Matrix"], mfRotated[mfUnitary, Cos[mfSpectrum]]]
    ],
    _ ? (# < 10^-12 &),
    SameTest -> MatchQ,
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-Ising64-random-eigenbasis"
]

(* Highly degenerate spectrum: eigenvalue 1 eight times, -1/2 four times. *)
VerificationTest[
    With[{
        u = Orthogonalize[mfUnitary[[;; 16, ;; 16]]],
        spectrum = Join[ConstantArray[1., 8], ConstantArray[-0.5, 4], {2., 2., 3., 0.1}]
    },
        mfDistance[
            Sin[QuantumOperator[mfRotated[u, spectrum], Range[4]]]["Matrix"],
            mfRotated[u, Sin[spectrum]]
        ]
    ],
    _ ? (# < 10^-12 &),
    SameTest -> MatchQ,
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-degenerate-spectrum"
]

(* A non-analytic f on a Hermitian operator is its spectral function: |H| = U |h| U^†. *)
VerificationTest[
    With[{qo = QuantumOperator[mfRotated[mfUnitary, mfSpectrum], Range[6]]},
        mfDistance[Abs[qo]["Matrix"], mfRotated[mfUnitary, Abs[mfSpectrum]]]
    ],
    _ ? (# < 10^-12 &),
    SameTest -> MatchQ,
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-Abs-Hermitian"
]

(* A non-normal operator with a known eigenbasis V: f(A) = V f(D) V^-1. *)
VerificationTest[
    With[{
        v = {{1, 2, 0, 1}, {0, 1, 1, 0}, {1, 0, 1, 1}, {0, 1, 0, 2}},
        d = {0.3, -1.2, 2.5, 0.7}
    },
        mfDistance[
            Cos[QuantumOperator[v . DiagonalMatrix[d] . Inverse[v], {1, 2}]]["Matrix"],
            v . DiagonalMatrix[Cos[d]] . Inverse[v]
        ]
    ],
    _ ? (# < 10^-12 &),
    SameTest -> MatchQ,
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-non-normal"
]

(* Symbolic angles: cos(t X) = cos t I, sin(t ZZ) = sin t ZZ. *)
VerificationTest[
    Simplify[Normal[Cos[mfT QuantumOperator["X"]]["Matrix"]] - Cos[mfT] IdentityMatrix[2]],
    ConstantArray[0, {2, 2}],
    TestID -> "MatrixFunction-symbolic-Cos-X"
]

VerificationTest[
    Simplify[Normal[Sin[mfT QuantumOperator["ZZ"]]["Matrix"]] - Sin[mfT] DiagonalMatrix[{1, -1, -1, 1}]],
    ConstantArray[0, {4, 4}],
    TestID -> "MatrixFunction-symbolic-Sin-ZZ"
]

(* XX + YY has spectrum {0, 0, 2, -2}: it vanishes on |00>, |11> and acts as 2 X
   on span{|01>, |10>}. *)
VerificationTest[
    With[{qo = mfT (QuantumOperator["XX"] + QuantumOperator["YY"])},
        Simplify[{
            Normal[Cos[qo]["Matrix"]] - {{1, 0, 0, 0}, {0, Cos[2 mfT], 0, 0}, {0, 0, Cos[2 mfT], 0}, {0, 0, 0, 1}},
            Normal[Sin[qo]["Matrix"]] - {{0, 0, 0, 0}, {0, 0, Sin[2 mfT], 0}, {0, Sin[2 mfT], 0, 0}, {0, 0, 0, 0}}
        }]
    ],
    ConstantArray[0, {2, 4, 4}],
    TestID -> "MatrixFunction-symbolic-degenerate-XX+YY"
]

(* On the non-orthonormal PauliX basis R stores diag(1, -1) and represents X, so
   cos(t R) represents cos t I and sin(t R) represents sin t X. *)
VerificationTest[
    With[{r = QuantumOperator[{{1, 0}, {0, -1}}, "PauliX"]},
        Simplify[{
            Normal[Cos[mfT r]["MatrixRepresentation"]] - Cos[mfT] IdentityMatrix[2],
            Normal[Sin[mfT r]["MatrixRepresentation"]] - Sin[mfT] PauliMatrix[1]
        }]
    ],
    ConstantArray[0, {2, 2, 2}],
    TestID -> "MatrixFunction-nonorthonormal-basis"
]

VerificationTest[
    With[{qo = QuantumOperator[{{0, mfTheta}, {mfTheta, 0}}, "Parameters" -> {mfTheta}]},
        {Cos[qo]["Parameters"], Simplify[Normal[Cos[qo]["Matrix"]] - Cos[mfTheta] IdentityMatrix[2]], Normal[Cos[qo][Pi / 3]["Matrix"]]}
    ],
    {{mfTheta}, ConstantArray[0, {2, 2}], IdentityMatrix[2] / 2},
    TestID -> "MatrixFunction-parameters"
]

(* Exact input stays exact. *)
VerificationTest[
    With[{m = Log[QuantumOperator[{{2, 1}, {1, 2}}]]["Matrix"]},
        {Precision[m], Simplify[Normal[m] - Log[3] / 2 {{1, 1}, {1, 1}}]}
    ],
    {Infinity, ConstantArray[0, {2, 2}]},
    TestID -> "MatrixFunction-exact-Log"
]

VerificationTest[
    Simplify[Normal[Log[2, QuantumOperator[{{2, 1}, {1, 2}}]]["Matrix"]] - Log[2, 3] / 2 {{1, 1}, {1, 1}}],
    ConstantArray[0, {2, 2}],
    TestID -> "MatrixFunction-exact-two-argument-Log"
]

VerificationTest[
    Normal[Cos[QuantumOperator[DiagonalMatrix[{1/2, 1/2, 1/3, 2}], {1, 2}]]["Matrix"]],
    DiagonalMatrix[Cos[{1/2, 1/2, 1/3, 2}]],
    TestID -> "MatrixFunction-exact-degenerate-diagonal"
]

(* A Jordan block needs the derivative: cos [[1, 1], [0, 1]] = [[cos 1, -sin 1], [0, cos 1]]. *)
VerificationTest[
    Simplify[Normal[Cos[QuantumOperator[{{1, 1}, {0, 1}}]]["Matrix"]] - {{Cos[1], -Sin[1]}, {0, Cos[1]}}],
    ConstantArray[0, {2, 2}],
    TestID -> "MatrixFunction-exact-Jordan-block"
]

VerificationTest[
    Cos[QuantumOperator[mfHadamard[[;; 4, ;; 4]], {2, 5}]]["Order"],
    {{2, 5}, {2, 5}},
    TestID -> "MatrixFunction-keeps-order"
]

(* On a state f acts on the density matrix: sqrt(rho) is the positive Hermitian
   root and squares back to rho, here for a Gibbs-like state on the degenerate Ising
   spectrum in a random eigenbasis. *)
VerificationTest[
    With[{rho = QuantumState[mfRotated[mfUnitary, Normalize[Exp[-mfSpectrum], Total]], 2^6]},
        With[{root = Normal[Sqrt[rho]["DensityMatrix"]]},
            {mfDistance[root . root, rho["DensityMatrix"]], Min[Re[Eigenvalues[root]]], root === ConjugateTranspose[root]}
        ]
    ],
    {_ ? (# < 10^-12 &), _ ? (# > -10^-12 &), True},
    SameTest -> MatchQ,
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-state-Sqrt"
]

(* Nearly normal: a Jordan-like block with a 10^-9 coupling. Its Schur form is not
   diagonal, so the coupling must survive: cos [[1, e], [0, 1]] = [[cos 1, -e sin 1], [0, cos 1]]. *)
VerificationTest[
    mfDistance[
        Cos[QuantumOperator[{{1., 10.^-9}, {0., 1.}}]]["Matrix"],
        {{Cos[1.], -10.^-9 Sin[1.]}, {0., Cos[1.]}}
    ],
    _ ? (# < 10^-15 &),
    SameTest -> MatchQ,
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-nearly-normal"
]

(* A pure state is a projector, so its square root is itself: the zero eigenvalues
   must map to zero, not to the square root of their roundoff. *)
VerificationTest[
    With[{rho = QuantumState[BlockRandom[Normalize[RandomComplex[{-1 - I, 1 + I}, 16]], RandomSeeding -> 5]]},
        mfDistance[Sqrt[rho]["DensityMatrix"], rho["DensityMatrix"]]
    ],
    _ ? (# < 10^-12 &),
    SameTest -> MatchQ,
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-pure-state-Sqrt"
]

(* Small Hermitian matrices, many draws: Abs must never fall off the spectral route.
   The reference diagonalizes at 30 digits, and |H|^2 = H^2 checks it independently. *)
mfAbsReference[h_] := With[
    {es = Eigensystem[SetPrecision[h, 30]]},
    {v = Orthogonalize[es[[2]]]},
    N[Transpose[v] . DiagonalMatrix[Abs[es[[1]]]] . Conjugate[v]]
]

VerificationTest[
    BlockRandom[
        Max @ Table[
            With[{m = RandomComplex[{-1 - I, 1 + I}, {d, d}]},
                With[{h = (m + ConjugateTranspose[m]) / 2},
                    With[{abs = Normal[Abs[QuantumOperator[h]]["Matrix"]]},
                        Max[mfDistance[abs, mfAbsReference[h]], mfDistance[abs . abs, h . h]]
                    ]
                ]
            ],
            {d, {2, 4}}, {300}
        ],
        RandomSeeding -> 17
    ],
    _ ? (# < 10^-12 &),
    SameTest -> MatchQ,
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-Abs-Hermitian-sweep"
]

VerificationTest[
    BlockRandom[
        Max @ Table[
            With[{rho = QuantumState[Normalize[RandomComplex[{-1 - I, 1 + I}, 2]]]},
                mfDistance[Sqrt[rho]["DensityMatrix"], rho["DensityMatrix"]]
            ],
            {300}
        ],
        RandomSeeding -> 19
    ],
    _ ? (# < 10^-12 &),
    SameTest -> MatchQ,
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-qubit-pure-state-Sqrt-sweep"
]

(* Log of a singular operator does not exist: the result is a Failure, not a matrix
   of infinities, on the dense and on the diagonal route alike. *)
VerificationTest[
    Head /@ {Log[QuantumOperator[{{0.5, 0.5}, {0.5, 0.5}}]], Log[QuantumOperator[{{0., 0.}, {0., 1.}}]]},
    {Failure, Failure},
    TestID -> "MatrixFunction-singular-Log-fails"
]

(* Symbolic entries with a machine coefficient keep the compact closed form. *)
VerificationTest[
    Simplify[Normal[Cos[0.5 mfT QuantumOperator["X"]]["Matrix"]] - Cos[0.5 mfT] IdentityMatrix[2]],
    ConstantArray[0, {2, 2}],
    SameTest -> (Chop[#1] === #2 &),
    TestID -> "MatrixFunction-inexact-symbolic"
]

(* A full-rank state with a small but resolved eigenvalue is not singular: Log stays
   finite and matches R.log(p).R^T. *)
VerificationTest[
    With[{r = RotationMatrix[0.3], p = {1 - 10.^-12, 10.^-12}},
        With[{result = Log[QuantumOperator[r . DiagonalMatrix[p] . Transpose[r]]]},
            {Head[result], mfDistance[result["Matrix"], r . DiagonalMatrix[Log[p]] . Transpose[r]] < 10^-3}
        ]
    ],
    {QuantumOperator, True},
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-small-resolved-eigenvalue"
]

(* A diagonalizable non-normal matrix takes any f through its eigenbasis:
   A = [[1, 1], [0, -2]] = V.diag(1, -2).V^-1 with V = [[1, 1], [0, -3]], so
   |A| = V.diag(1, 2).V^-1 = [[1, -1/3], [0, 2]]. *)
VerificationTest[
    mfDistance[Abs[QuantumOperator[{{1., 1.}, {0., -2.}}]]["Matrix"], {{1, -1/3}, {0, 2}}],
    _ ? (# < 10^-14 &),
    SameTest -> MatchQ,
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-Abs-non-normal-diagonalizable"
]

(* Real symmetric input gives a real symmetric result. *)
VerificationTest[
    With[{m = Normal[Log[QuantumOperator[mfHadamard . DiagonalMatrix[mfSpectrum + 10] . mfHadamard, Range[6]]]["Matrix"]]},
        {FreeQ[m, _Complex], m === Transpose[m]}
    ],
    {True, True},
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-real-symmetric-stays-real"
]

(* Functional-calculus identities on a complex Hermitian H with the degenerate Ising
   spectrum: cos^2 H + sin^2 H = 1, [cos H, H] = 0, and cos H is exactly Hermitian. *)
VerificationTest[
    With[{h = mfRotated[mfUnitary, mfSpectrum]},
        With[{c = Normal[Cos[QuantumOperator[h, Range[6]]]["Matrix"]], s = Normal[Sin[QuantumOperator[h, Range[6]]]["Matrix"]]},
            {mfDistance[c . c + s . s, IdentityMatrix[64]], mfDistance[c . h, h . c], c === ConjugateTranspose[c]}
        ]
    ],
    {_ ? (# < 10^-12 &), _ ? (# < 10^-12 &), True},
    SameTest -> MatchQ,
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-functional-calculus-identities"
]

(* Eigenvalues that depend on the parameters: [[a, b], [b, -a]] squares to r^2 with
   r = sqrt(a^2 + b^2), so cos of it is cos r times the identity. *)
VerificationTest[
    Simplify[
        Normal[Cos[QuantumOperator[{{mfA, mfB}, {mfB, -mfA}}]]["Matrix"]] - Cos[Sqrt[mfA^2 + mfB^2]] IdentityMatrix[2],
        Element[{mfA, mfB}, Reals]
    ],
    ConstantArray[0, {2, 2}],
    TestID -> "MatrixFunction-symbolic-parameter-dependent-spectrum"
]

(* Across the route switches: [[1, 1], [0, 1 + d]] runs from well separated to
   defective as d -> 0, and cos of it is [[cos 1, (cos(1 + d) - cos 1)/d], [0, cos(1 + d)]]. *)
VerificationTest[
    Max @ Table[
        With[{d = 10^-k},
            mfDistance[
                Cos[QuantumOperator[N[{{1, 1}, {0, 1 + d}}]]]["Matrix"],
                N[{{Cos[1], (Cos[1 + d] - Cos[1]) / d}, {0, Cos[1 + d]}}, 30]
            ]
        ],
        {k, 15}
    ],
    _ ? (# < 10^-12 &),
    SameTest -> MatchQ,
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-near-defective-sweep"
]

(* Log of an exact singular operator and of a defective numeric one is a Failure. *)
VerificationTest[
    Head /@ {Log[QuantumOperator[{{1, 1}, {1, 1}} / 2]], Log[QuantumOperator[{{0., 1.}, {0., 0.}}]]},
    {Failure, Failure},
    {Infinity::indet},
    TestID -> "MatrixFunction-singular-Log-exact-and-defective-fails"
]

(* A normal, non-Hermitian operator: a random three-qubit unitary U. Its principal
   square root is unitary and squares back to U, and exp(log U) = U. *)
VerificationTest[
    With[{u = Orthogonalize[mfUnitary[[;; 8, ;; 8]]]},
        With[{root = Normal[Sqrt[QuantumOperator[u, Range[3]]]["Matrix"]], log = Normal[Log[QuantumOperator[u, Range[3]]]["Matrix"]]},
            {mfDistance[root . root, u], mfDistance[root . ConjugateTranspose[root], IdentityMatrix[8]], mfDistance[MatrixExp[log], u]}
        ]
    ],
    {_ ? (# < 10^-12 &) ..},
    SameTest -> MatchQ,
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-unitary"
]

(* -Tr rho log2 rho through the operator logarithm is the von Neumann entropy. *)
VerificationTest[
    With[{rho = QuantumState[mfRotated[mfUnitary, Normalize[Exp[-mfSpectrum], Total]], 2^6]},
        Abs[-Tr[rho["DensityMatrix"] . Log[2, rho]["DensityMatrix"]] - QuantityMagnitude[rho["VonNeumannEntropy"]]]
    ],
    _ ? (# < 10^-12 &),
    SameTest -> MatchQ,
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-von-Neumann-entropy"
]

(* f commutes with a change of basis near a zero eigenvalue: a roundoff-level
   eigenvalue is zero on the diagonal and in a rotated basis alike. *)
VerificationTest[
    With[{r = RotationMatrix[0.3], d = DiagonalMatrix[{1., -10.^-17}]},
        {
            Head /@ {Log[QuantumOperator[d]], Log[QuantumOperator[r . d . Transpose[r]]]},
            mfDistance[r . Normal[Sqrt[QuantumOperator[d]]["Matrix"]] . Transpose[r], Sqrt[QuantumOperator[r . d . Transpose[r]]]["Matrix"]]
        }
    ],
    {{Failure, Failure}, _ ? (# < 10^-15 &)},
    SameTest -> MatchQ,
    TestID -> "MatrixFunction-diagonal-and-rotated-agree"
]

(* A non-analytic f of a defective matrix is undefined: a Failure, with the
   built-in message saying why. *)
VerificationTest[
    Head[Abs[QuantumOperator[{{1., 1.}, {0., 1.}}]]],
    Failure,
    {MatrixFunction::drvnnum},
    TestID -> "MatrixFunction-Abs-defective-fails"
]

(* An interacting exact Hamiltonian: the open 4-qubit Heisenberg chain
   H = sum_i (X_i X_i+1 + Y_i Y_i+1 + Z_i Z_i+1), spectrum 3 (5 times), -1 (3 times),
   -1 +- 2 sqrt 2 (3 times each), -3 +- 2 sqrt 3. The reference is the sum over
   eigenspaces of cos(t lambda) times the orthogonal projector, built from an exact
   Eigensystem; the two closed forms must agree identically in t. *)
mfHeisenberg[n_] := Total[Flatten[Table[QuantumOperator[p, {i, i + 1}], {i, n - 1}, {p, {"XX", "YY", "ZZ"}}]]]

VerificationTest[
    With[
        {h = Normal[mfHeisenberg[4]["Matrix"]]},
        {eigenspaces = GatherBy[Transpose[Eigensystem[h]], First]},
        {reference = Total[
            With[{v = Orthogonalize[#[[All, 2]]]}, Cos[mfT #[[1, 1]]] Transpose[v] . Conjugate[v]] & /@ eigenspaces
        ]},
        Simplify[Normal[Cos[mfT mfHeisenberg[4]]["Matrix"]] - reference]
    ],
    ConstantArray[0, {16, 16}],
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-symbolic-Heisenberg-chain"
]

(* M = [[a, b], [b, -a]] squares to r^2 with r = sqrt(a^2 + b^2), so its sine series
   sums to sin M = sinc(r) M, regular at r = 0 where the two eigenvalues +-r
   collide. The closed form f[qo] shows is sin(r)/r M, generic in (a, b). *)
mfCollision = QuantumOperator[{{mfA, mfB}, {mfB, -mfA}}, "Parameters" -> {mfA, mfB}];

(* Close to the collision the closed form is still exact: it equals sinc(r) M
   identically at a = 10^-8, b = 2 10^-8. *)
VerificationTest[
    With[{a = 10^-8, b = 2 10^-8},
        FullSimplify[Normal[Sin[mfCollision][a, b]["Matrix"]] - Sinc[Sqrt[a^2 + b^2]] {{a, b}, {b, -a}}]
    ],
    ConstantArray[0, {2, 2}],
    TestID -> "MatrixFunction-parameter-near-collision"
]

(* At the collision itself, substituting the parameters before applying f gives the
   right answer, sin 0 = 0. *)
VerificationTest[
    Normal[Sin[mfCollision[0, 0]]["Matrix"]],
    ConstantArray[0, {2, 2}],
    TestID -> "MatrixFunction-parameter-collision-substitute-first"
]

(* Substituting values into f[qo] applies f to the substituted matrix rather than
   evaluating the closed form sin(r)/r M at r = 0, so at the collision the result is
   sin 0 = 0 with no messages, whether the values arrive at once, one parameter at a
   time, or positionally. *)
VerificationTest[
    Normal /@ {
        Sin[mfCollision][<|mfA -> 0, mfB -> 0|>]["Matrix"],
        Sin[mfCollision][<|mfA -> 0|>][<|mfB -> 0|>]["Matrix"],
        Sin[mfCollision][0, 0]["Matrix"]
    },
    ConstantArray[0, {3, 2, 2}],
    TestID -> "MatrixFunction-parameter-collision-substitute-after"
]

(* A collision that makes the matrix defective: [[a, 1], [0, b]] at a = b is a
   Jordan block, and sin of it needs the derivative, [[sin b, cos b], [0, sin b]],
   also when a -> b is substituted symbolically. *)
VerificationTest[
    With[{jordan = QuantumOperator[{{mfA, 1}, {0, mfB}}, "Parameters" -> {mfA, mfB}]},
        Simplify[{
            Normal[Sin[jordan][1, 1]["Matrix"]] - {{Sin[1], Cos[1]}, {0, Sin[1]}},
            Normal[Sin[jordan][<|mfA -> mfB|>]["Matrix"]] - {{Sin[mfB], Cos[mfB]}, {0, Sin[mfB]}}
        }]
    ],
    ConstantArray[0, {2, 2, 2}],
    TestID -> "MatrixFunction-parameter-collision-defective"
]

(* Near, not at, a collision with machine values: the closed form's divided
   difference (sin a - sin b)/(a - b) cancels catastrophically at a - b = 10^-12;
   applying f to the substituted matrix is accurate. The reference is the exact
   matrix at 30 digits. *)
VerificationTest[
    With[{jordan = QuantumOperator[{{mfA, 1}, {0, mfB}}, "Parameters" -> {mfA, mfB}]},
        Max[Abs[
            Normal[Sin[jordan][1., 1. + 10.^-12]["Matrix"]] -
                N[MatrixFunction[Sin, {{1, 1}, {0, 1 + 10^-12}}], 30]
        ]]
    ],
    _ ? (# < 10^-12 &),
    SameTest -> MatchQ,
    TestID -> "MatrixFunction-parameter-near-collision-machine"
]

(* Where f itself is not defined at the substituted matrix, the substitution fails:
   Log of [[a, b], [b, -a]] at a = b = 0 is the log of the zero matrix. *)
VerificationTest[
    Head[Log[mfCollision][0, 0]],
    Failure,
    TestID -> "MatrixFunction-parameter-substitution-Log-fails"
]

(* A substitution fails when it leaves without a value an amplitude that had one, read
   amplitude by amplitude: in a state whose first amplitude is infinite for s < -5, the
   value s = -1 leaves the second one Indeterminate and fails, while s = 2 gives both a
   value; and s = -10 leaves the first one infinite, which had a value before the
   substitution chose its case, and fails. A substitution that leaves a Piecewise
   amplitude with an infinite case still unchosen, as t = 2 in
   Piecewise[{{Infinity, s < -5}}, t], leaves it a value and gives the state. *)
VerificationTest[
    With[{qs = QuantumState[{Piecewise[{{Infinity, mfS < -5}}, 1], Piecewise[{{1, mfS > 0}}, Indeterminate]}, "Parameters" -> {mfS}]},
        {
            FailureQ[qs[-1]], Normal[qs[2]["StateVector"]], FailureQ[QuantumState[{Piecewise[{{Infinity, mfS < -5}}, 1], 1}, "Parameters" -> {mfS}][-10]],
            Head[QuantumState[{Piecewise[{{Infinity, mfS < -5}}, mfT], 1}, "Parameters" -> {mfS, mfT}][<|mfT -> 2|>]]
        }
    ],
    {True, {1, 1}, True, QuantumState},
    TestID -> "MatrixFunction-substitution-guard-per-amplitude"
]

(* Parameters that cannot be Function variables (indexed th[1], th[2]) keep the
   closed form: still sin(r)/r M, and correct away from the collision. *)
VerificationTest[
    With[{q = QuantumOperator[{{th[1], th[2]}, {th[2], -th[1]}}, "Parameters" -> {th[1], th[2]}]},
        Simplify[{
            Normal[Sin[q]["Matrix"]] - Sin[Sqrt[th[1]^2 + th[2]^2]] / Sqrt[th[1]^2 + th[2]^2] {{th[1], th[2]}, {th[2], -th[1]}},
            Normal[Sin[q][3/10, 4/10]["Matrix"]] - 2 Sin[1/2] {{3/10, 4/10}, {4/10, -3/10}}
        }]
    ],
    ConstantArray[0, {2, 2, 2}],
    TestID -> "MatrixFunction-parameter-indexed-parameters"
]

(* A non-square operator has no matrix function, with or without parameters. *)
VerificationTest[
    Head[Sin[QuantumOperator[{{mfA, 1}, {0, mfA}, {1, 0}, {0, 1}}, "Parameters" -> {mfA}]]],
    Failure,
    TestID -> "MatrixFunction-parameter-nonsquare-fails"
]

(* Abs is not differentiable, and a Jordan block needs the derivative: substituting
   a = b = 1 into Abs of [[a, 1], [0, b]] fails rather than returning an unevaluated
   Derivative[1][Abs][1]. *)
VerificationTest[
    Head[Abs[QuantumOperator[{{mfA, 1}, {0, mfB}}, "Parameters" -> {mfA, mfB}]][1, 1]],
    Failure,
    TestID -> "MatrixFunction-parameter-Abs-defective-fails"
]

(* A symbolic second argument keeps a symbolic result: log base b of [[2, 1], [1, 2]]
   (eigenvalues 3 and 1) is log 3/(2 log b) [[1, 1], [1, 1]], with or without
   parameters in the operator. *)
VerificationTest[
    Simplify[{
        Normal[Log[mfB, QuantumOperator[{{2, 1}, {1, 2}}]]["Matrix"]] - Log[3] / (2 Log[mfB]) {{1, 1}, {1, 1}},
        Normal[Log[mfB, QuantumOperator[{{mfA, 1}, {1, mfA}}, "Parameters" -> {mfA}]][2]["Matrix"]] - Log[3] / (2 Log[mfB]) {{1, 1}, {1, 1}}
    }],
    ConstantArray[0, {2, 2, 2}],
    TestID -> "MatrixFunction-symbolic-second-argument"
]

(* f of f: exp(sin M) at the collision a = b = 0 is exp 0 = 1, where exp of the closed
   form sin(r)/r M is 0/0; one parameter at a time too, and away from the collision
   it matches the built-in MatrixExp of MatrixFunction[Sin, m]. The same holds when
   the inner result is placed on the reversed order {2, 1} of a two-qubit
   operator a Z(x)X + b X(x)I, whose order the outer function sorts. *)
VerificationTest[
    With[{
        m = {{1/3, 1/5}, {1/5, -1/3}},
        reordered = QuantumOperator[
            Sin[QuantumOperator[mfA KroneckerProduct[PauliMatrix[3], PauliMatrix[1]] + mfB KroneckerProduct[PauliMatrix[1], IdentityMatrix[2]], {1, 2}, "Parameters" -> {mfA, mfB}]],
            {2, 1}
        ]
    },
        Simplify[{
            Normal[Exp[Sin[mfCollision]][0, 0]["Matrix"]] - IdentityMatrix[2],
            Normal[Exp[Sin[mfCollision]][<|mfA -> 0|>][<|mfB -> 0|>]["Matrix"]] - IdentityMatrix[2],
            Normal[Exp[Sin[mfCollision]][1/3, 1/5]["Matrix"]] - MatrixExp[MatrixFunction[Sin, m]],
            Normal[Exp[reordered][0, 0]["Matrix"]][[;; 2, ;; 2]] - IdentityMatrix[2],
            Normal[Exp[reordered][0, 0]["Matrix"]][[3 ;;, 3 ;;]] - IdentityMatrix[2]
        }]
    ],
    ConstantArray[0, {5, 2, 2}],
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-parameter-nested-collision"
]

(* The exponential spellings take the same route: exp M, e^M, 2^M and MatrixExp at the
   collision a = b = 0 are the identity, and exp of [[a, 1], [0, b]] at a = b = 1 is
   e [[1, 1], [0, 1]]. *)
VerificationTest[
    {
        Normal /@ {
            Exp[mfCollision][0, 0]["Matrix"],
            (E ^ mfCollision)[0, 0]["Matrix"],
            (2 ^ mfCollision)[0, 0]["Matrix"],
            MatrixExp[mfCollision][0, 0]["Matrix"]
        },
        Normal[Exp[QuantumOperator[{{mfA, 1}, {0, mfB}}, "Parameters" -> {mfA, mfB}]][1, 1]["Matrix"]]
    },
    {ConstantArray[IdentityMatrix[2], 4], E {{1, 1}, {0, 1}}},
    TestID -> "MatrixFunction-parameter-exponential-collision"
]

(* A zero eigenvalue that does not depend on the parameter: log of [[a, 1], [0, 0]]
   exists for no a, so it fails at construction. *)
VerificationTest[
    Head[Log[QuantumOperator[{{mfA, 1}, {0, 0}}, "Parameters" -> {mfA}]]],
    Failure,
    TestID -> "MatrixFunction-parameter-undefined-everywhere-fails"
]

(* Unitarity and Hermiticity on the parametric route: for H = [[w, g], [g, -w]],
   exp(-i t H) is unitary and cos H Hermitian identically in real (t, w, g), and at
   the collision w = g = 0 the propagator is the identity. *)
VerificationTest[
    With[{h = QuantumOperator[{{mfA, mfB}, {mfB, -mfA}}, "Parameters" -> {mfT, mfA, mfB}]},
        With[
            {u = Normal[Exp[-I mfT h]["Matrix"]], c = Normal[Cos[h]["Matrix"]]},
            {
                FullSimplify[ComplexExpand[u . ConjugateTranspose[u]], Element[{mfT, mfA, mfB}, Reals]],
                FullSimplify[ComplexExpand[c - ConjugateTranspose[c]], Element[{mfA, mfB}, Reals]],
                Normal[Exp[-I mfT h][<|mfT -> 1, mfA -> 0, mfB -> 0|>]["Matrix"]]
            }
        ]
    ],
    {IdentityMatrix[2], ConstantArray[0, {2, 2}], IdentityMatrix[2]},
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-parameter-unitarity-hermiticity"
]

(* A many-body level crossing: the 3-qubit transverse-field Ising chain
   H = J (Z1 Z2 + Z2 Z3) + g (X1 + X2 + X3). Its closed-form propagator has about
   9 million leaves and is not finite at g = 0; the lazy propagator equals the
   built-in MatrixExp at the crossing (J, g) = (1, 0) and at a generic point. *)
mfTFIM = With[{z = PauliMatrix[3], x = PauliMatrix[1], i2 = IdentityMatrix[2]},
    mfA (KroneckerProduct[z, z, i2] + KroneckerProduct[i2, z, z]) +
        mfB (KroneckerProduct[x, i2, i2] + KroneckerProduct[i2, x, i2] + KroneckerProduct[i2, i2, x])
];

VerificationTest[
    With[{u = Exp[QuantumOperator[-I mfTFIM, {1, 2, 3}, "Parameters" -> {mfA, mfB}]]},
        {
            Simplify[Normal[u[1, 0]["Matrix"]] - MatrixExp[-I mfTFIM /. {mfA -> 1, mfB -> 0}]],
            Max[Abs[Normal[u[0.7, 0.3]["Matrix"]] - MatrixExp[-I mfTFIM /. {mfA -> 0.7, mfB -> 0.3}]]] < 10^-12
        }
    ],
    {ConstantArray[0, {8, 8}], True},
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-parameter-transverse-Ising-crossing"
]

(* Decoupled subsystems: exp(a X (x) 1 + b 1 (x) Z) = exp(a X) (x) exp(b Z), in closed
   form and after substitution, including b = 0 where the Z factor is the identity. *)
VerificationTest[
    With[{q = QuantumOperator[mfA KroneckerProduct[PauliMatrix[1], IdentityMatrix[2]] + mfB KroneckerProduct[IdentityMatrix[2], PauliMatrix[3]], {1, 2}, "Parameters" -> {mfA, mfB}]},
        Simplify[{
            Normal[Exp[q]["Matrix"]] - KroneckerProduct[MatrixExp[mfA PauliMatrix[1]], MatrixExp[mfB PauliMatrix[3]]],
            Normal[Exp[q][1/3, 0]["Matrix"]] - KroneckerProduct[MatrixExp[PauliMatrix[1] / 3], IdentityMatrix[2]]
        }]
    ],
    ConstantArray[0, {2, 4, 4}],
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-parameter-decoupled-subsystems"
]

(* Short time: exp(-i t H) = 1 - i t H + O(t^2) for H = [[a, b], [b, -a]]. *)
VerificationTest[
    Simplify[
        Normal[Series[Normal[Exp[-I mfT QuantumOperator[{{mfA, mfB}, {mfB, -mfA}}, "Parameters" -> {mfT, mfA, mfB}]]["Matrix"]], {mfT, 0, 1}]] -
            (IdentityMatrix[2] - I mfT {{mfA, mfB}, {mfB, -mfA}})
    ],
    ConstantArray[0, {2, 2}],
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-parameter-short-time"
]

(* cos(t sum_i Z_i) over n = 2 .. 6 qubits: diagonal with cos(t m), m the
   magnetization of each basis state. *)
VerificationTest[
    Table[
        Simplify[
            Normal[Cos[mfT Total[QuantumOperator["Z", {#}] & /@ Range[n]]]["Matrix"]] -
                DiagonalMatrix[Cos[mfT Total[1 - 2 Tuples[{0, 1}, n], {2}]]]
        ] === ConstantArray[0, {2^n, 2^n}],
        {n, 2, 6}
    ],
    ConstantArray[True, 5],
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-symbolic-magnetization-sweep"
]

(* Known limitation: the lazy route stops at operations applied afterwards (a
   product with another operator, a sum, an action on a state); those are built from
   the closed form, so values at a collision give a Failure rather than Indeterminate
   amplitudes, and substituting first gives the answer, sin 0 X = 0. When these
   operations become lazy, the Failure should become the zero matrix. *)
VerificationTest[
    {
        Head[(Sin[mfCollision] @ QuantumOperator["X"])[0, 0]],
        Normal[(Sin[mfCollision[0, 0]] @ QuantumOperator["X"])["Matrix"]]
    },
    {Failure, ConstantArray[0, {2, 2}]},
    {Power::infy, Infinity::indet, Power::infy, Infinity::indet, Power::infy, General::stop, Infinity::indet, General::stop},
    TestID -> "MatrixFunction-parameter-after-composition-fails"
]

(* Known limitation, near a collision: a product built from the closed form keeps
   its divided difference, so sin([[a, 1], [0, b]]) X at a = 1., b = 1. + 10^-12 is
   off by about 10^-4 with no message, while sin of the substituted matrix times X is
   accurate. When products become lazy, the first error should drop to roundoff. *)
VerificationTest[
    With[{
        jordan = QuantumOperator[{{mfA, 1}, {0, mfB}}, "Parameters" -> {mfA, mfB}],
        reference = N[MatrixFunction[Sin, {{1, 1}, {0, 1 + 10^-12}}] . PauliMatrix[1], 30]
    },
        {
            Max[Abs[Normal[(Sin[jordan] @ QuantumOperator["X"])[1., 1. + 10.^-12]["Matrix"]] - reference]] > 10^-6,
            Max[Abs[Normal[(Sin[jordan[1., 1. + 10.^-12]] @ QuantumOperator["X"])["Matrix"]] - reference]] < 10^-12
        }
    ],
    {True, True},
    TestID -> "MatrixFunction-parameter-product-near-collision-known-limitation"
]

(* A real matrix with a real-analytic f gives a real result, also off the Hermitian
   branch: cos of a rotation matrix (normal, complex eigenvalues) and of the
   non-normal [[0, 1], [-2, 1/2]]. *)
VerificationTest[
    {
        FreeQ[Normal[Cos[QuantumOperator[N[RotationMatrix[3/10]]]]["Matrix"]], _Complex],
        FreeQ[Normal[Cos[QuantumOperator[{{0., 1.}, {-2., 0.5}}]]["Matrix"]], _Complex],
        Max[Abs[Normal[Cos[QuantumOperator[N[RotationMatrix[3/10]]]]["Matrix"]] - N[MatrixFunction[Cos, RotationMatrix[3/10]], 30]]] < 10^-14
    },
    {True, True, True},
    TestID -> "MatrixFunction-real-input-real-result"
]

(* The operator as base, M^p, is a matrix power: X^p = (1 + (-1)^p)/2 1 + (1 - (-1)^p)/2 X,
   X^(1/2) is the square root, and on the lazy route q^(1/2) at the collision of
   q = [[a, b], [b, -a]] is 0. *)
VerificationTest[
    {
        Simplify[Normal[(QuantumOperator["X"] ^ mfT)["Matrix"]] - ((1 + (-1)^mfT) / 2 IdentityMatrix[2] + (1 - (-1)^mfT) / 2 PauliMatrix[1])],
        (QuantumOperator["X"] ^ (1/2))["Matrix"] === Sqrt[QuantumOperator["X"]]["Matrix"],
        Normal[(mfCollision ^ (1/2))[0, 0]["Matrix"]]
    },
    {ConstantArray[0, {2, 2}], True, ConstantArray[0, {2, 2}]},
    TestID -> "MatrixFunction-operator-base-power"
]

(* 0^M is the limit of b^M = MatrixExp[Log[b] M] as b -> 0: 1 on a semisimple zero
   eigenvalue and 0 on eigenvalues with positive real part, so the spectral projector
   onto the null space of M along its range. Each test returns its physical claims by
   name, and a failure is read by the tag of the Failure that 0^M itself returns. *)
mfZeroBaseTag[Failure[tag_, _]] := tag

mfZeroBaseTag[_] := None

(* The matrix of 0^M for an operator or a matrix M. *)
mfZero[m_] := Normal[(0 ^ If[MatrixQ[m], QuantumOperator[m], m])["Matrix"]]

(* MatrixPower[a, 0] of the singular annihilator stays unevaluated with
   MatrixPower::sing, so its powers are built by Dot. *)
mfPower[x_, k_] := Nest[x . # &, IdentityMatrix[Length[x]], k]

(* The Choi matrix of a superoperator acting on row-major vectorized d x d matrices. *)
mfChoi[p_, d_] := Flatten[Transpose[ArrayReshape[p, {d, d, d, d}], {1, 3, 2, 4}], {{1, 2}, {3, 4}}]

(* For a number operator n, 0^n is the vacuum projector, the value the literature
   uses: the thermal state (1 - q) q^n as q -> 0, ThermalState at nbar = 0 included,
   which leaves the vacuum at first order in q; the total-loss channel, whose Kraus
   operators carry eta^(n/2), at eta = 0, where every state goes to the vacuum; and
   Louisell's normal-ordered sum_k (-1)^k a^dagger^k a^k / k!, the lambda = 1 case
   of (1 - lambda)^n. Also for three modes through their total number, for a machine
   number operator in a rotated basis, and for a parameter base sent to zero. *)
VerificationTest[
    With[
        {a = AnnihilationOperator[5], a3 = AnnihilationOperator[3], u = {{1, 1}, {1, -1}} / Sqrt[2], vacuum = DiagonalMatrix[{1, 0, 0, 0, 0}]},
        {number = a["Dagger"] @ a, annihilator = Normal[a["Matrix"]], state = Outer[Times, {0, 1, 0, 1, 0}, {0, 1, 0, 1, 0}] / 2},
        {
            thermal = (1 - mfQ) mfQ ^ QuantumOperator[number, "Parameters" -> {mfQ}],
            base = mfT ^ QuantumOperator[number, "Parameters" -> {mfT}],
            loss = Normal[(mfEta ^ QuantumOperator[number / 2, "Parameters" -> {mfEta}])[0]["Matrix"]],
            rotated = mfZero[N[u . DiagonalMatrix[{0, 1}] . u]],
            normalOrdered = Function[lambda, Total[Table[(-lambda)^k mfPower[ConjugateTranspose[annihilator], k] . mfPower[annihilator, k] / k!, {k, 0, 4}]]]
        },
        {kraus = Table[loss . mfPower[annihilator, k] / Sqrt[k!], {k, 0, 4}], closedThermal = Normal[thermal["Matrix"]]},
        <|
            "one mode" -> mfZero[number] == vacuum,
            "three modes" -> mfZero[Total[QuantumOperator[a3["Dagger"] @ a3, {#}] & /@ {1, 2, 3}]] == Outer[Times, UnitVector[27, 1], UnitVector[27, 1]],
            "rotated machine" -> Max[Abs[rotated - N[u . DiagonalMatrix[{1, 0}] . u]]] < 10^-14 && rotated === Transpose[rotated],
            "machine" -> mfZero[N[number]] === N[vacuum],
            "parameter base" -> Normal[base[0]["Matrix"]] == vacuum && Normal[base[1/2]["Matrix"]] == DiagonalMatrix[2^-Range[0, 4]],
            "ThermalState at nbar = 0" -> Normal[ThermalState[0, 5]["DensityMatrix"]] == vacuum,
            "thermal limit" -> Normal[thermal[0]["Matrix"]] == Limit[closedThermal, mfQ -> 0] == vacuum,
            "thermal first order" -> Limit[(closedThermal - vacuum) / mfQ, mfQ -> 0] == DiagonalMatrix[{-1, 1, 0, 0, 0}],
            "total loss" -> Total[ConjugateTranspose[#] . # & /@ kraus] == IdentityMatrix[5] && Total[# . state . ConjugateTranspose[#] & /@ kraus] == vacuum,
            "Louisell" -> normalOrdered[1] == vacuum,
            "Louisell (1 - lambda)^n" -> Simplify[normalOrdered[mfLambda] - Normal[((1 - mfLambda) ^ QuantumOperator[number, "Parameters" -> {mfLambda}])["Matrix"]]] ==
                ConstantArray[0, {5, 5}]
        |>
    ],
    <|
        "one mode" -> True, "three modes" -> True, "rotated machine" -> True, "machine" -> True, "parameter base" -> True, "ThermalState at nbar = 0" -> True,
        "thermal limit" -> True, "thermal first order" -> True, "total loss" -> True, "Louisell" -> True, "Louisell (1 - lambda)^n" -> True
    |>,
    TestID -> "MatrixFunction-zero-base-vacuum-projector"
]

(* 0^(H - E0) is the projector onto the ground space, degenerate and entangled for
   an interacting H, and Hermitian with H.P = P.H = 0, which tells it from an oblique
   projector onto the same null space: the XX chain -(X1 X2 + X2 X3), whose ground
   space is spanned by |+++> and |--->; the Heisenberg triangle, whose ground space
   is total spin 1/2, P = (15/4 - S^2)/3 = (3 - H)/6; the ferromagnetic ring of n
   sites, the sum of bond singlet projectors (1 - SWAP)/2, whose ground space is the
   symmetric subspace of dimension n + 1, for n = 4, 6, 8; and B^dagger B for
   B = (a1 + a2)/Sqrt[2], whose kernel holds the antisymmetric (|01> - |10>)/Sqrt[2].
   The XXZ pair XX + YY + d ZZ has a level crossing at d = -1: its ground space is the
   singlet for every d > -1 and the doublet {|00>, |11>} for every d < -1, both in
   closed form, and with the ground energy Min[d, -2 - d] as the shift, the values
   of d give the ranks 1, 3, 2, 2 at d = 0, -1, -2, -3. The ground-space projector
   is the zero-temperature limit of the Boltzmann operator exp(-beta (H - E0)),
   approached as exp(-beta gap) times the projector onto the first excited level,
   itself 0^((H - E0 - gap)^2). In a field h on every site the ring keeps only the
   fully polarized |0...0> for every h > 0, so the limits in b and h do not
   commute, in closed form for every h > 0 at n = 4; at n = 8 and h = 10^-12 that gap
   lies inside the roundoff bound of the machine eigenvalues, and the precision is
   raised until it clears it. *)
VerificationTest[
    With[
        {
            xx = -(QuantumOperator["XX", {1, 2}] + QuantumOperator["XX", {2, 3}]),
            heisenberg = Total[QuantumOperator[#, {1, 2}] + QuantumOperator[#, {2, 3}] + QuantumOperator[#, {1, 3}] & /@ {"XX", "YY", "ZZ"}],
            ring = Function[n, Total[Table[(QuantumOperator[IdentityMatrix[4], {i, Mod[i, n] + 1}] - QuantumOperator["SWAP", {i, Mod[i, n] + 1}]) / 2, {i, n}]]],
            field = Total[QuantumOperator[{{0, 0}, {0, 1}}, {#}] & /@ Range[8]],
            b = (AnnihilationOperator[3, {1}] + AnnihilationOperator[3, {2}]) / Sqrt[2],
            xxz = QuantumOperator["XX"] + QuantumOperator["YY"] + mfD QuantumOperator["ZZ"],
            singlet = Outer[Times, {0, 1, -1, 0}, {0, 1, -1, 0}] / 2,
            doublet = DiagonalMatrix[{1, 0, 0, 1}]
        },
        {
            ground = mfZero[xx + 2],
            shifted = Normal[(xx + 2)["Matrix"]],
            orthogonalQ = Function[{h, p}, p == ConjugateTranspose[p] && p . p == p && Normal[h["Matrix"]] . p == 0 p && p . Normal[h["Matrix"]] == 0 p],
            crossing = 0 ^ QuantumOperator[xxz - Min[mfD, -2 - mfD], "Parameters" -> {mfD}]
        },
        {gap = Min[DeleteCases[Eigenvalues[shifted], 0]]},
        {excited = mfZero[QuantumOperator[(shifted - gap IdentityMatrix[8]) . (shifted - gap IdentityMatrix[8]), {1, 2, 3}]]},
        <|
            "XX chain" -> ground == Total[Normal[QuantumState[#]["DensityMatrix"]] & /@ {"+++", "---"}] && orthogonalQ[xx + 2, ground],
            "Heisenberg triangle" -> mfZero[heisenberg + 3] == (3 IdentityMatrix[8] - Normal[heisenberg["Matrix"]]) / 6,
            "ferromagnetic ring" -> (With[{h = ring[#]}, {p = mfZero[ring[#]]}, orthogonalQ[h, p] && Tr[p] == # + 1] & /@ {4, 6, 8}),
            "ring in a field" -> (With[{p = mfZero[ring[8] + # field]}, {Tr[p], p[[1, 1]]}] & /@ {0, 10^-6, 10^-12}) == {{9, 1}, {1, 1}, {1, 1}},
            "ring in a symbolic field" -> Simplify[mfZero[ring[4] + mfH Total[QuantumOperator[{{0, 0}, {0, 1}}, {#}] & /@ Range[4]]], mfH > 0] ==
                Outer[Times, UnitVector[16, 1], UnitVector[16, 1]],
            "entangled kernel" -> With[{bb = b["Dagger"] @ b}, {p = mfZero[b["Dagger"] @ b], antisymmetric = UnitVector[9, 2] - UnitVector[9, 4]},
                orthogonalQ[bb, p] && Tr[p] == 3 && p . antisymmetric == antisymmetric],
            "XXZ crossing" -> Simplify[mfZero[xxz + 2 + mfD], mfD > -1] == singlet && Simplify[mfZero[xxz - mfD], mfD < -1] == doublet &&
                (MatrixRank[Normal[crossing[#]["Matrix"]]] & /@ {0, -1, -2, -3}) == {1, 3, 2, 2} &&
                Normal[crossing[-1]["Matrix"]] == singlet + doublet && Normal[crossing[-3]["Matrix"]] == doublet,
            "zero temperature" -> Simplify[Limit[MatrixExp[-mfBeta shifted], mfBeta -> Infinity] - ground] == ConstantArray[0, {8, 8}],
            "approach as exp(-beta gap)" -> Tr[excited] == 4 && Simplify[Limit[Exp[gap mfBeta] (MatrixExp[-mfBeta shifted] - ground), mfBeta -> Infinity]] == excited
        |>
    ],
    <|
        "XX chain" -> True, "Heisenberg triangle" -> True, "ferromagnetic ring" -> {True, True, True}, "ring in a field" -> True,
        "ring in a symbolic field" -> True, "entangled kernel" -> True,
        "XXZ crossing" -> True, "zero temperature" -> True, "approach as exp(-beta gap)" -> True
    |>,
    TestID -> "MatrixFunction-zero-base-ground-space"
]

(* For a Lindblad generator L, b^(-L) = exp(t L) with b = exp(-t), so 0^(-L) is the
   limit of the evolution as t -> Infinity: the projector P onto the steady states,
   with P.P = P and L.P = P.L = 0, trace preserving, and completely positive (its
   Choi matrix is positive), in closed form for every rate for the driven qubit and
   the thermal generator. Amplitude damping sends every state to |0><0| (the rows
   read vec(|0><0|) Tr) for every rate g > 0; dephasing keeps the populations and
   removes the coherences. The resonantly driven qubit, H = Omega X / 2 with decay at
   rate gamma, has in closed form the resonance-fluorescence steady state,
   rho_11 = Omega^2 / (gamma^2 + 2 Omega^2) and rho_01 = I gamma Omega / (gamma^2 +
   2 Omega^2), where the condition that the other eigenvalues (gamma/2 and
   (3 gamma +- Sqrt[gamma^2 - 16 Omega^2])/4) have positive real part holds on both
   sides of the exceptional point gamma = 4 Omega and at it, where two of them
   collide and the closed form still gives rho_11 = 1/18. The thermal generator,
   decay at rate nbar + 1 and excitation at rate nbar, has ThermalState as its steady
   state for every nbar > 0, and the vacuum at nbar = 0. Collective decay of two
   qubits, J = s1^- + s2^-, leaves a noiseless qubit on {|00>, |S>} with the singlet
   |S> dark, so P is oblique (P != P^dagger) and the steady state depends on the
   initial one: |01><01| goes to (|00><00| + |S><S|)/2 and |11><11| to |00><00|.
   For three qubits the states that J annihilates span three dimensions, so the
   steady states span nine. With H = Z1 + Z2 the coherence |00><S| rotates forever
   and there is no limit, nor for a unitary generator, whose spectrum is
   imaginary. *)
VerificationTest[
    With[
        {
            lowering = QuantumOperator[{{0, 1}, {0, 0}}],
            collectiveJump = QuantumOperator[{{0, 1}, {0, 0}}, {1}] + QuantumOperator[{{0, 1}, {0, 0}}, {2}],
            a3 = AnnihilationOperator[3],
            density = Normal[QuantumState[#]["DensityMatrix"]] &
        },
        {
            damping = QuantumOperator["Liouvillian"[None, {lowering}, {1}]],
            dampingRate = QuantumOperator["Liouvillian"[None, {lowering}, {mfG}]],
            dephasing = QuantumOperator["Liouvillian"[None, {QuantumOperator["Z"]}, {1}]],
            driven = QuantumOperator["Liouvillian"[mfOmega QuantumOperator["X"] / 2, {lowering}, {mfGamma}]],
            thermalGenerator = QuantumOperator["Liouvillian"[None, {a3, a3["Dagger"]}, {mfN + 1, mfN}]],
            collective = QuantumOperator["Liouvillian"[None, {collectiveJump}, {1}]],
            collective3 = QuantumOperator["Liouvillian"[None, {Total[QuantumOperator[{{0, 1}, {0, 0}}, {#}] & /@ {1, 2, 3}]}, {1}]],
            singletState = With[{s = Normal[QuantumState["01"]["StateVector"] - QuantumState["10"]["StateVector"]] / Sqrt[2]}, Outer[Times, s, s]],
            steady = Function[l, mfZero[-l]],
            everyRateQ = Function[{l, d, assumptions},
                With[{p = mfZero[-l] /. Piecewise[{{value_, _}}, _] :> value, lm = Normal[l["Matrix"]]},
                    Simplify[{p . p - p, lm . p, p . lm}, assumptions] === ConstantArray[0, {3, d^2, d^2}] &&
                        Simplify[Flatten[IdentityMatrix[d]] . p - Flatten[IdentityMatrix[d]], assumptions] === ConstantArray[0, d^2] &&
                        TrueQ[Simplify[And @@ Thread[Eigenvalues[mfChoi[p, d]] >= 0], assumptions]]
                ]
            ],
            apply = Function[{p, rho}, ArrayReshape[p . Flatten[rho], Dimensions[rho]]]
        },
        {
            projectorQ = Function[{l, d}, With[{p = steady[l], lm = Normal[l["Matrix"]]},
                p . p == p && lm . p == 0 lm && p . lm == 0 lm && Flatten[IdentityMatrix[d]] . p == Flatten[IdentityMatrix[d]] &&
                    PositiveSemidefiniteMatrixQ[mfChoi[p, d]] && Simplify[Limit[MatrixExp[mfTime lm], mfTime -> Infinity] - p] == 0 lm
            ]],
            fluorescence = apply[steady[driven], density["0"]],
            thermalSteady = apply[steady[thermalGenerator], DiagonalMatrix[{1, 0, 0}]]
        },
        <|
            "projector invariants" -> {projectorQ[damping, 2], projectorQ[dephasing, 2], projectorQ[collective, 4]},
            "projector invariants for every rate" -> {everyRateQ[driven, 2, mfOmega > 0 && mfGamma > 0], everyRateQ[thermalGenerator, 3, mfN > 0]},
            "amplitude damping" -> steady[damping] == {{1, 0, 0, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}} &&
                Simplify[steady[dampingRate], mfG > 0] == steady[damping],
            "dephasing" -> steady[dephasing] == DiagonalMatrix[{1, 0, 0, 1}],
            "resonance fluorescence" -> (fluorescence /. Piecewise[{{value_, _}}, _] :> value) ==
                {{mfGamma^2 + mfOmega^2, I mfGamma mfOmega}, {-I mfGamma mfOmega, mfOmega^2}} / (mfGamma^2 + 2 mfOmega^2),
            "condition holds for every rate" -> Reduce[ForAll[{mfOmega, mfGamma}, mfOmega > 0 && mfGamma >= 4 mfOmega, 3 mfGamma > Sqrt[mfGamma^2 - 16 mfOmega^2]], Reals] &&
                Refine[Re[Sqrt[mfGamma^2 - 16 mfOmega^2]], 0 < mfGamma < 4 mfOmega] == 0,
            "exceptional point" -> (fluorescence /. {mfOmega -> 1, mfGamma -> 4}) == {{17/18, 2 I / 9}, {-2 I / 9, 1/18}},
            "thermal steady state" -> Simplify[thermalSteady - Normal[ThermalState[mfN, 3]["DensityMatrix"]], mfN > 0] == ConstantArray[0, {3, 3}] &&
                (thermalSteady /. mfN -> 0) == DiagonalMatrix[{1, 0, 0}],
            "collective decay" -> With[{p = steady[collective]},
                Tr[p] == 4 && p =!= ConjugateTranspose[p] && apply[p, density["01"]] == (density["00"] + singletState) / 2 && apply[p, density["11"]] == density["00"]],
            "collective decay of three qubits" -> With[{p = steady[collective3], lm = Normal[collective3["Matrix"]]},
                Tr[p] == 9 && p . p == p && lm . p == 0 lm && p . lm == 0 lm && Flatten[IdentityMatrix[8]] . p == Flatten[IdentityMatrix[8]] &&
                    PositiveSemidefiniteMatrixQ[mfChoi[p, 8]]],
            "no limit" -> mfZeroBaseTag /@ {
                0 ^ -QuantumOperator["Liouvillian"[QuantumOperator["Z", {1}] + QuantumOperator["Z", {2}], {collectiveJump}, {1}]],
                0 ^ -QuantumOperator["Liouvillian"[QuantumOperator["X"], {}, {}]]
            }
        |>
    ],
    <|
        "projector invariants" -> {True, True, True}, "projector invariants for every rate" -> {True, True}, "amplitude damping" -> True, "dephasing" -> True,
        "resonance fluorescence" -> True, "condition holds for every rate" -> True, "exceptional point" -> True,
        "thermal steady state" -> True, "collective decay" -> True, "collective decay of three qubits" -> True,
        "no limit" -> {"ZeroBasePowerNoLimit", "ZeroBasePowerNoLimit"}
    |>,
    TestID -> "MatrixFunction-zero-base-steady-state"
]

(* The projector P = 0^M of M = S.(0_k (+) D).S^-1, with Re D > 0 and D holding a
   Jordan block or a complex eigenvalue, is S.(1_k (+) 0).S^-1: P.P = P,
   M.P = P.M = 0, Tr P = k. It is the same for c M at every c > 0, in closed form, it
   factors over a sum A (x) 1 + 1 (x) B of commuting parts, and the machine and
   40-digit routes agree. *)
VerificationTest[
    With[
        {
            s4 = {{1, 0, 0, 0}, {1, 1, 0, 0}, {0, 1, 1, 0}, {1, 0, 1, 1}} . {{1, 2, 0, 1}, {0, 1, 1, 0}, {0, 0, 1, 2}, {0, 0, 0, 1}},
            s5 = {{1, 0, 0, 0, 0}, {0, 1, 0, 0, 0}, {1, 1, 1, 0, 0}, {0, 1, 0, 1, 0}, {1, 0, 1, 1, 1}} .
                {{1, 1, 0, 2, 0}, {0, 1, 1, 0, 1}, {0, 0, 1, 1, 0}, {0, 0, 0, 1, 1}, {0, 0, 0, 0, 1}}
        },
        {
            cases = {
                {s4, 2, ArrayFlatten[{{ConstantArray[0, {2, 2}], 0}, {0, {{2, 1}, {0, 2}}}}]},
                {s5, 3, ArrayFlatten[{{ConstantArray[0, {3, 3}], 0}, {0, {{1 + I, 0}, {0, 1/2}}}}]}
            },
            invariants = Function[{sim, k, block},
                With[{m = sim . block . Inverse[sim]}, {p = mfZero[sim . block . Inverse[sim]]},
                    <|
                        "idempotent" -> p . p == p,
                        "annihilates M" -> m . p == 0 m && p . m == 0 m,
                        "trace is nullity" -> Tr[p] == k,
                        "projector by construction" -> p == sim . DiagonalMatrix[Join[ConstantArray[1, k], ConstantArray[0, 2]]] . Inverse[sim],
                        "scale invariance" -> Simplify[mfZero[mfC m] - p, mfC > 0] == 0 m,
                        "machine agrees" -> Max[Abs[mfZero[N[m]] - p]] < 10^-12,
                        "40 digits agree" -> Max[Abs[mfZero[N[m, 40]] - p]] < 10^-35
                    |>
                ]
            ]
        },
        {a = s4 . cases[[1, 3]] . Inverse[s4]},
        Append[
            Merge[invariants @@@ cases, Identity],
            "factorizes" -> mfZero[KroneckerProduct[a, IdentityMatrix[2]] + KroneckerProduct[IdentityMatrix[4], DiagonalMatrix[{0, 3}]]] ==
                KroneckerProduct[mfZero[a], DiagonalMatrix[{1, 0}]]
        ]
    ],
    <|
        "idempotent" -> {True, True}, "annihilates M" -> {True, True}, "trace is nullity" -> {True, True},
        "projector by construction" -> {True, True}, "scale invariance" -> {True, True}, "machine agrees" -> {True, True}, "40 digits agree" -> {True, True},
        "factorizes" -> True
    |>,
    TestID -> "MatrixFunction-zero-base-invariants"
]

(* Approaching a Jordan block, M = {{e, 1}, {0, 0}} has 0^M = {{0, -1/e}, {0, 1}} for
   every e > 0, whose norm grows as 1/e, while at fixed b the e -> 0 limit of b^M is
   {{1, Log[b]}, {0, 1}}: the limits b -> 0 and e -> 0 do not commute. The machine
   route follows the closed form while e is resolvable, and fails as defective once
   the matrix is within roundoff of the block. *)
VerificationTest[
    With[
        {closed = mfZero[{{mfE, 1}, {0, 0}}], machine = Function[e, 0 ^ QuantumOperator[N[{{e, 1}, {0, 0}}]]]},
        <|
            "closed form" -> Simplify[closed, mfE > 0] == {{0, -1 / mfE}, {0, 1}},
            "limits do not commute" -> Limit[MatrixExp[Log[mfB] {{mfE, 1}, {0, 0}}], mfE -> 0] == {{1, Log[mfB]}, {0, 1}},
            "machine follows" -> AllTrue[{10^-1, 10^-3, 10^-5}, Max[Abs[Normal[machine[#]["Matrix"]] - {{0, -1 / #}, {0, 1}}]] # < 10^-10 &],
            "machine at roundoff" -> mfZeroBaseTag[machine[10^-9]]
        |>
    ],
    <|"closed form" -> True, "limits do not commute" -> True, "machine follows" -> True, "machine at roundoff" -> "ZeroBasePowerDefective"|>,
    TestID -> "MatrixFunction-zero-base-jordan-approach"
]

(* With every eigenvalue of positive real part the limit is the zero matrix, which
   is also what MatrixFunction[0^# &, m] gives there. A nonzero base is the matrix
   exponential. *)
VerificationTest[
    With[{m = {{3/2, 1/2}, {1/2, 3/2}}},
        <|
            "diagonal" -> mfZero[DiagonalMatrix[{1, 2}]] == ConstantArray[0, {2, 2}],
            "as MatrixFunction" -> mfZero[m] == MatrixFunction[0^# &, m],
            "complex" -> mfZero[DiagonalMatrix[{1 + I, 2}]] == ConstantArray[0, {2, 2}],
            "machine" -> Max[Abs[Normal[(0. ^ QuantumOperator[N[m]])["Matrix"]]]] == 0,
            "parametric" -> Normal[(0 ^ QuantumOperator[mfT DiagonalMatrix[{1, 2}], "Parameters" -> {mfT}])[1/2]["Matrix"]] == ConstantArray[0, {2, 2}],
            "nonzero base" -> Normal[(2 ^ QuantumOperator[DiagonalMatrix[{1, 2}]])["Matrix"]] == {{2, 0}, {0, 4}}
        |>
    ],
    <|"diagonal" -> True, "as MatrixFunction" -> True, "complex" -> True, "machine" -> True, "parametric" -> True, "nonzero base" -> True|>,
    TestID -> "MatrixFunction-zero-base-positive-spectrum"
]

(* The limit does not exist for an eigenvalue that is neither zero nor of positive
   real part (b^-1 diverges, b^I oscillates) or for a Jordan block at zero, where
   b^M = 1 + Log[b] N diverges. The result is the Failure that names the reason,
   without messages: a Jordan block at zero exact or machine, triangular or
   conjugated, and a machine matrix within roundoff of one, are all defective, and a
   real part that is 0 but not written as 0 has no limit. For exact input the
   message names the eigenvalue exactly when it is a root of a factor of the
   characteristic polynomial with real coefficients or of a linear or quadratic
   factor, as for a root of z^3 - z - 1, or, for a Hermitian matrix of Gaussian
   rationals, a rational or quadratic irrational that its computed value at 64 digits
   identifies, as (1 - Sqrt[5])/2, -10^-12, -10^-130, -10^-400 and
   -123456789/1000000007 are. Otherwise a Hermitian
   matrix names it to the digits its roundoff bound leaves, as for a root of
   x^3 - x^2 - 2 x + 1, and any other matrix says only that one exists, as for the
   roots of z^3 + I z + 1 and of z^3 + ArcCos[2] z + 1, whose coefficients are complex
   by value. Past the range of machine numbers it is named the same way,
   as the smallest eigenvalues of {{1, 1}, {1, 1 - 10^-400}} and
   {{10^400, 1}, {1, 10^-400 - 10^-1000}} are. A
   machine diagonal is held to the roundoff tolerance, so -10^-12 beside 1 has no
   limit. Other matrix functions return the Failure that names their reason too,
   as Log at a zero eigenvalue and Abs at a Jordan block, which needs its derivative,
   do. *)
VerificationTest[
    With[
        {
            r3 = {{1, 2, 2}, {2, 1, -2}, {2, -2, 1}} / 3,
            s20 = {{1, 10^20}, {1, 10^20 + 1}},
            named = (0 ^ QuantumOperator[#])["MessageParameters"] &
        },
        Join[
            mfZeroBaseTag /@ <|
                "eigenvalue -1" -> 0 ^ QuantumOperator[DiagonalMatrix[{0, -1}]],
                "eigenvalues +-I" -> 0 ^ QuantumOperator[{{0, 1}, {-1, 0}}],
                "eigenvalue -10^-130" -> 0 ^ QuantumOperator[r3 . DiagonalMatrix[{0, -10^-130, 1}] . Transpose[r3]],
                "Jordan block" -> 0 ^ QuantumOperator[{{0, 1}, {0, 0}}],
                "machine Jordan block" -> 0 ^ QuantumOperator[N[{{0, 1}, {0, 0}}]],
                "conjugated Jordan block" -> 0 ^ QuantumOperator[{{3, -1}, {9, -3}}],
                "machine conjugated Jordan block" -> 0 ^ QuantumOperator[N[{{3, -1}, {9, -3}}]],
                "machine nilpotent" -> 0 ^ QuantumOperator[N[{{-1, 1}, {-1, 1}}]],
                "nearly defective" -> 0 ^ QuantumOperator[{{10.^-17, 1.}, {0., 0.}}],
                "real part 0 not written as 0" -> 0 ^ QuantumOperator[DiagonalMatrix[{I + (Sqrt[2] + Sqrt[3])^2 - 5 - 2 Sqrt[6], 1}]],
                "eigenvalue -(Sqrt[2] - 1)^100 expanded" -> 0 ^ QuantumOperator[DiagonalMatrix[{-Expand[(Sqrt[2] - 1)^100], 1}]],
                "machine diagonal -10^-12" -> 0 ^ QuantumOperator[DiagonalMatrix[N[{-10^-12, 1}]]],
                "eigenvalue I, basis of condition 10^40" -> 0 ^ QuantumOperator[s20 . DiagonalMatrix[{I, 1}] . Inverse[s20]],
                "eigenvalue -10^-40, basis of condition 10^40" -> 0 ^ QuantumOperator[s20 . DiagonalMatrix[{-10^-40, 1}] . Inverse[s20]],
                "Log at a zero eigenvalue" -> Log[QuantumOperator[DiagonalMatrix[{0, 1}]]],
                "Abs at a Jordan block" -> Abs[QuantumOperator[{{1, 1}, {0, 1}}]]
            |>,
            <|
                "eigenvalues named" -> named /@ {
                    {{0, 1}, {-1, 0}}, {{0, 1}, {2, 0}}, {{1, 2}, {3, 4}}, {{0, 1}, {1, 0}},
                    s20 . DiagonalMatrix[{I, 1}] . Inverse[s20], s20 . DiagonalMatrix[{-10^-40, 1}] . Inverse[s20],
                    {{0, 0, 1}, {1, 0, 1}, {0, 1, 0}}, {{0, 0, -1}, {1, 0, -I}, {0, 1, 0}}, {{0, 0, -1}, {1, 0, -ArcCos[2]}, {0, 1, 0}},
                    {{0, 1}, {1, 1}}, r3 . DiagonalMatrix[{0, -10^-12, 1}] . Transpose[r3], r3 . DiagonalMatrix[{0, -10^-130, 1}] . Transpose[r3],
                    r3 . DiagonalMatrix[{0, -10^-400, 1}] . Transpose[r3], r3 . DiagonalMatrix[{0, -123456789/1000000007, 1}] . Transpose[r3],
                    {{Exp[1000], 1}, {2, 0}}
                },
                "Hermitian, to its digits" -> With[{v = First[named[{{0, 1, 0}, {1, 0, 1}, {0, 1, 1}}]]},
                    InexactNumberQ[v] && Abs[v - Root[#^3 - #^2 - 2 # + 1 &, 1]] < 10^-14],
                "Hermitian, past the machine range" -> Function[mat,
                    With[{v = First[named[mat]], smallest = Det[mat] / ((Tr[mat] + Sqrt[Tr[mat]^2 - 4 Det[mat]]) / 2)}, InexactNumberQ[v] && Abs[v / smallest - 1] < 10^-14]
                ] /@ {{{1, 1}, {1, 1 - 10^-400}}, {{10^400, 1}, {1, 10^-400 - 10^-1000}}}
            |>
        ]
    ],
    <|
        "eigenvalue -1" -> "ZeroBasePowerNoLimit", "eigenvalues +-I" -> "ZeroBasePowerNoLimit", "eigenvalue -10^-130" -> "ZeroBasePowerNoLimit",
        "Jordan block" -> "ZeroBasePowerDefective", "machine Jordan block" -> "ZeroBasePowerDefective",
        "conjugated Jordan block" -> "ZeroBasePowerDefective", "machine conjugated Jordan block" -> "ZeroBasePowerDefective",
        "machine nilpotent" -> "ZeroBasePowerDefective", "nearly defective" -> "ZeroBasePowerDefective",
        "real part 0 not written as 0" -> "ZeroBasePowerNoLimit", "eigenvalue -(Sqrt[2] - 1)^100 expanded" -> "ZeroBasePowerNoLimit",
        "machine diagonal -10^-12" -> "ZeroBasePowerNoLimit",
        "eigenvalue I, basis of condition 10^40" -> "ZeroBasePowerNoLimit",
        "eigenvalue -10^-40, basis of condition 10^40" -> "ZeroBasePowerNoLimit", "Log at a zero eigenvalue" -> "NonFiniteMatrixFunction",
        "Abs at a Jordan block" -> "NonFiniteMatrixFunction",
        "eigenvalues named" -> {
            {-I}, {-Sqrt[2]}, {(5 - Sqrt[33]) / 2}, {-1}, {I}, {-10^-40}, {Root[-1 - #1 + #1^3 &, 2, 0]}, Missing["NotAvailable", "MessageParameters"], Missing["NotAvailable", "MessageParameters"],
            {(1 - Sqrt[5]) / 2}, {-10^-12}, {-10^-130}, {-10^-400}, {-123456789/1000000007}, {-4 / (E^1000 + Sqrt[8 + E^2000])}
        },
        "Hermitian, to its digits" -> True,
        "Hermitian, past the machine range" -> {True, True}
    |>,
    TestID -> "MatrixFunction-zero-base-undefined-fails"
]

(* The null space of an inexact matrix is read off its singular values, within 100 n
   eps of the norm, and its zero eigenvalues may move by that tolerance times their
   condition number: a non-normal machine matrix with eigenvalues 3, 2, 0 gives the
   projector of its exact form, so does a rotated {{10^-4, 1}, {0, 0}}, whose zero
   eigenvalue has condition number near 10^4, a machine number operator of dimension
   128 in a random orthogonal basis gives its vacuum projector (trace 1, exactly
   symmetric), and a 40-digit matrix keeps its precision. A positive spectrum gives
   machine zeros for machine input. An exact matrix that is not Hermitian is decided
   exactly from its characteristic polynomial, and a Hermitian one from eigenvalues
   whose precision is raised until the roundoff bound decides, so the projector does
   not depend on the basis however wide the spectrum or ill-conditioned the basis:
   eigenvalues 0, 10^-12, 10^20 and 0, 1, 10^150, non-normal or Hermitian, 0, 10^-130,
   1, and 10^-40, 1 or 10^-20 + I, 1, 0 in bases of condition 10^40, and a Jordan
   block of order 3 at 10^-16 beside the null space, whose eigenvalue a numerical
   computation moves by 10^-10. Entries whose values cancel over tens to hundreds of
   digits are evaluated with the working precision the cancellation needs: a
   Hermitian matrix with entries 14142135623730951 - 10^16 Sqrt[2] and the like is
   positive definite or indefinite as the exact entries say, (Sqrt[2] - 1)^100 on the
   diagonal is positive, and conjugated eigenvalues Sqrt[10^40 + 1] - 10^20 and
   (Sqrt[2] - 1)^30 beside zero give the projector, while -(Sqrt[10^400 + 1] - 10^200)
   has no limit and is named exactly. A complex value written without I, ArcCos[2] =
   I Log[2 + Sqrt[3]], is told from a real one by value: {{2, a}, {a, a^2/2}} for
   a = ArcCos[2] is not Hermitian and its projector is the oblique x x^T / (x^T x);
   an eigenvalue -1/10 + 2 a has no limit; a Jordan block at 1/10 + a has one; and a
   nilpotent matrix of such entries is defective. The real part b = Im[ArcCos[2]] of
   such a value is real however it is written: {{1 + I b, 1}, {1, (1 - I b)/(1 + b^2)}}
   is not Hermitian and gets its oblique projector, eigenvalues I b and I (b - 1/b)
   have no limit, and {{1, I b}, {-I b, b^2}} is Hermitian. A number is read as
   algebraic when it is after FunctionExpand, so 10^-4000 Sin[ArcCos[1/3] / 2] =
   10^-4000 / Sqrt[3] is positive, not zero; a zero on a branch cut, as
   Sqrt[-1 - g I] - I or Log[-1 + g I] - I Pi for g = GoldenRatio - (1 + Sqrt[5]) / 2,
   is read as zero, alone or inside an entry, and so is Log[6] - Log[2] - Log[3]
   inside entries; a part read as zero by the convention, as -Exp[-10^8], is zero
   throughout, so r3.diag(0, -Exp[-10^8], 1).r3^T has a null space of dimension 2;
   Exp[10^7] is read without all of its digits; and spectra 0, 1 + Sqrt[2],
   2 + Sqrt[3], 5, 4 + Sqrt[5] and 0, 0, 2 + Sqrt[3], 5, 4 + Sqrt[5], 5 + Sqrt[6] in
   integer bases give their exact projectors. The roundoff bound holds past 307
   digits, where 10^-p no longer fits a machine number, and an entry that holds a zero
   not written as 0 is evaluated to an accuracy, which N reaches. *)
VerificationTest[
    With[
        {
            m = {{-2, 3, 4}, {6, -3, -6}, {-8, 6, 10}},
            ill = {{3/5, -4/5}, {4/5, 3/5}} . {{1/10000, 1}, {0, 0}} . {{3/5, 4/5}, {-4/5, 3/5}},
            u = BlockRandom[SeedRandom[3]; Orthogonalize[RandomReal[NormalDistribution[], {128, 128}]]],
            s3 = {{1, 0, 0}, {1, 1, 0}, {0, 1, 1}} . {{1, 2, 0}, {0, 1, 1}, {0, 0, 1}},
            r3 = {{1, 2, 2}, {2, 1, -2}, {2, -2, 1}} / 3,
            s20 = {{1, 10^20}, {1, 10^20 + 1}},
            s20Null = {{1, 10^20, 0}, {1, 10^20 + 1, 0}, {0, 1, 1}},
            sr = {{3, 1, 0, 2}, {1, 7, 1, 0}, {0, 2, 11, 1}, {1, 0, 1, 13}} / 3,
            zeroBasePower = Wolfram`QuantumFramework`PackageScope`zeroBasePower
        },
        {
            illExact = mfZero[ill],
            vacuum = mfZero[Transpose[u] . DiagonalMatrix[N[Range[0, 127]]] . u],
            digits40 = zeroBasePower[N[{{1, 1}, {1, 1}}, 40]],
            conjugated = Function[{sim, spectrum}, mfZero[sim . DiagonalMatrix[spectrum] . Inverse[sim]] == sim . DiagonalMatrix[Boole[PossibleZeroQ[spectrum]]] . Inverse[sim]]
        },
        <|
            "exact non-normal" -> mfZero[m] == {{1, -1, -1}, {-2, 2, 2}, {2, -2, -2}},
            "machine non-normal" -> Max[Abs[mfZero[N[m]] - {{1, -1, -1}, {-2, 2, 2}, {2, -2, -2}}]] < 10^-13,
            "ill-conditioned" -> Max[Abs[mfZero[N[ill]] - illExact]] < 10^-8 Max[Abs[illExact]],
            "rotated number operator" -> Max[Abs[vacuum - Outer[Times, u[[1]], u[[1]]]]] < 10^-12 && Abs[Tr[vacuum] - 1] < 10^-12 && vacuum === Transpose[vacuum],
            "40 digits" -> Precision[digits40] > 35 && Max[Abs[Normal[digits40] - {{1, -1}, {-1, 1}} / 2]] < 10^-35,
            "machine zeros" -> Precision[zeroBasePower[N[{{3/2, 1/2}, {1/2, 3/2}}]]] === MachinePrecision,
            "wide spectra, any basis" -> {
                conjugated[s3, {0, 10^-12, 10^20}], conjugated[s3, {0, 1, 10^150}],
                conjugated[r3, {0, 1, 10^150}], conjugated[r3, {0, 10^-130, 1}]
            },
            "ill-conditioned basis" -> {conjugated[s20, {10^-40, 1}], conjugated[s20, {10^-20 + I, 1}], conjugated[s20Null, {10^-20 + I, 1, 0}]},
            "defective beside the null space" -> With[{block = ArrayFlatten[{{{{0}}, 0}, {0, 10^-16 IdentityMatrix[3] + DiagonalMatrix[{1, 1}, 1]}}]},
                mfZero[sr . block . Inverse[sr]] == sr . DiagonalMatrix[{1, 0, 0, 0}] . Inverse[sr]],
            "cancellations" -> With[
                {
                    epsilon = Sqrt[10^400 + 1] - 10^200,
                    sameQ = Function[{a, b}, AllTrue[Flatten[RootReduce[a - b]], # === 0 &]],
                    positiveDefinite = {{14142135623730951 - 10^16 Sqrt[2], 1}, {1, (11/10) (14142135623730951 + 10^16 Sqrt[2]) / (14142135623730951^2 - 2 10^32)}},
                    indefinite = {{10^16 Sqrt[2] - 14142135623730950, 1}, {1, (9/10) (10^16 Sqrt[2] + 14142135623730950) / (2 10^32 - 14142135623730950^2)}}
                },
                {
                    mfZero[positiveDefinite] == ConstantArray[0, {2, 2}],
                    mfZeroBaseTag[0 ^ QuantumOperator[indefinite]],
                    mfZero[DiagonalMatrix[{(Sqrt[2] - 1)^100, 1}]] == ConstantArray[0, {2, 2}],
                    sameQ[mfZero[s3 . DiagonalMatrix[{0, Sqrt[10^40 + 1] - 10^20, 1}] . Inverse[s3]], s3 . DiagonalMatrix[{1, 0, 0}] . Inverse[s3]],
                    sameQ[mfZero[s3 . DiagonalMatrix[{0, (Sqrt[2] - 1)^30, 1}] . Inverse[s3]], s3 . DiagonalMatrix[{1, 0, 0}] . Inverse[s3]],
                    RootReduce[First[(0 ^ QuantumOperator[s3 . DiagonalMatrix[{0, -epsilon, 1}] . Inverse[s3]])["MessageParameters"]] + epsilon] === 0
                }
            ],
            "complex values written without I" -> With[{a = ArcCos[2], lambda = 1/10 + ArcCos[2]},
                {
                    Simplify[mfZero[{{2, a}, {a, a^2 / 2}}] - Outer[Times, {a, -2}, {a, -2}] / (a^2 + 4)] == ConstantArray[0, {2, 2}],
                    mfZeroBaseTag[0 ^ QuantumOperator[s3 . DiagonalMatrix[{0, -1/10 + 2 a, 1 - 2 a}] . Inverse[s3]]],
                    Simplify[mfZero[{{0, 1, 0}, {0, lambda, 1}, {0, 0, lambda}}] - {{1, -1 / lambda, 1 / lambda^2}, {0, 0, 0}, {0, 0, 0}}] == ConstantArray[0, {3, 3}],
                    mfZeroBaseTag[0 ^ QuantumOperator[{{ArcCosh[2], a}, {a, -ArcCosh[2]}}]]
                }
            ],
            "real parts of complex values" -> With[{b = Im[ArcCos[2]], x = {1, -(1 + I Im[ArcCos[2]])}, y = {-I Im[ArcCos[2]], 1}},
                {
                    Max[Abs[N[mfZero[{{1 + I b, 1}, {1, (1 - I b) / (1 + b^2)}}] - Outer[Times, x, x] / (x . x), {Infinity, 30}]]] < 10^-25,
                    mfZeroBaseTag[0 ^ QuantumOperator[s3 . DiagonalMatrix[{0, I b, 1 - I b}] . Inverse[s3]]],
                    mfZeroBaseTag[0 ^ QuantumOperator[{{I b, 1}, {1, - I / b}}]],
                    Max[Abs[N[mfZero[{{1, I b}, {- I b, b^2}}] - Outer[Times, y, Conjugate[y]] / (1 + b^2), {Infinity, 30}]]] < 10^-25
                }
            ],
            "algebraic after FunctionExpand" -> With[{alg = 10^-4000 Sin[ArcCos[1/3] / 2]},
                {mfZero[DiagonalMatrix[{alg, 1}]] == ConstantArray[0, {2, 2}], mfZeroBaseTag[0 ^ QuantumOperator[DiagonalMatrix[{- alg, 1}]]]}
            ],
            "zero on a branch cut" -> With[{g0 = GoldenRatio - (1 + Sqrt[5]) / 2},
                {
                    mfZero[DiagonalMatrix[{Sqrt[-1 - g0 I] - I, 1}]] == {{1, 0}, {0, 0}},
                    Max[Abs[N[mfZero[r3 . DiagonalMatrix[{Log[-1 + g0 I] - I Pi, 1, 2}] . Transpose[r3]] - Outer[Times, r3[[All, 1]], r3[[All, 1]]], {Infinity, 30}]]] < 10^-25,
                    mfZero[{{1, Sqrt[-1 + g0 I]}, {-I, 1}}] == {{1, -I}, {I, 1}} / 2
                }
            ],
            "zero written with Log" -> With[{z = Log[6] - Log[2] - Log[3]}, mfZero[{{3 z - 2, 2 - 2 z}, {3 z - 3, 3 - 2 z}}] == {{3, -2}, {3, -2}}],
            "eight square roots" -> With[{r = (Sqrt[2] + Sqrt[3] + Sqrt[5] + Sqrt[7] + Sqrt[11] + Sqrt[13] + Sqrt[17] + Sqrt[19])^2},
                {
                    Max[Abs[N[mfZero[{{r, r}, {Expand[r], Expand[r]}}] - {{1, -1}, {-1, 1}} / 2, {Infinity, 30}]]] < 10^-25,
                    mfZeroBaseTag[0 ^ QuantumOperator[DiagonalMatrix[{I + r - Expand[r], 1}]]],
                    mfZeroBaseTag[0 ^ QuantumOperator[{{r, 1}, {-1, - Expand[r]}}]]
                }
            ],
            "symbolic beside a part read as zero" -> Simplify[mfZero[{{mfW, 1}, {0, Exp[-10^8]}}] /. mfW -> 2] == {{0, -1/2}, {0, 1}},
            "inexact base that is not zero" -> Max[Abs[(Normal[((1.`20*^-400) ^ QuantumOperator[DiagonalMatrix[{0, -1}]])["Matrix"]] - {{1, 0}, {0, 10^400}}) / {{1, 1}, {1, 10^400}}]] < 10^-15,
            "part read as zero" -> mfZero[r3 . DiagonalMatrix[{0, - Exp[-10^8], 1}] . Transpose[r3]] == r3 . DiagonalMatrix[{1, 1, 0}] . Transpose[r3],
            "huge entry" -> mfZero[DiagonalMatrix[{Exp[10^7], 0}]] == {{0, 0}, {0, 1}},
            "algebraic spectrum" -> With[
                {
                    s = {{-2, -2, -2, 1, -1}, {-1, 1, -2, 2, 0}, {2, -1, 2, 2, 0}, {-2, -1, 0, -1, 2}, {1, -2, 0, -1, 0}},
                    s6 = {{-2, 0, 0, 1, -1, -2}, {1, 2, 0, 1, 0, -1}, {1, -1, -1, 1, 0, 2}, {-1, -1, 1, -1, 0, 0}, {1, 1, 0, 1, -2, 2}, {-1, 2, 0, -2, 2, -1}}
                },
                {
                    AllTrue[Flatten[RootReduce[mfZero[s . DiagonalMatrix[{0, 1 + Sqrt[2], 2 + Sqrt[3], 5, 4 + Sqrt[5]}] . Inverse[s]] - s . DiagonalMatrix[{1, 0, 0, 0, 0}] . Inverse[s]]], # === 0 &],
                    Max[Abs[N[mfZero[s6 . DiagonalMatrix[{0, 0, 2 + Sqrt[3], 5, 4 + Sqrt[5], 5 + Sqrt[6]}] . Inverse[s6]] - s6 . DiagonalMatrix[{1, 1, 0, 0, 0, 0}] . Inverse[s6], {Infinity, 30}]]] < 10^-25
                }
            ],
            "precision past 307 digits" -> {
                mfZero[r3 . DiagonalMatrix[{0, 10^-400, 3}] . Transpose[r3]] == Outer[Times, r3[[All, 1]], r3[[All, 1]]],
                Max[Abs[Normal[zeroBasePower[N[{{3/5, -4/5}, {4/5, 3/5}} . DiagonalMatrix[{0, 1}] . {{3/5, 4/5}, {-4/5, 3/5}}, 400]]] - N[{{9, 12}, {12, 16}} / 25, 400]]] < 10^-390
            },
            "entry holding a zero" -> mfZero[r3 . DiagonalMatrix[{GoldenRatio - (1 + Sqrt[5]) / 2, 1, 2}] . Transpose[r3]] == Outer[Times, r3[[All, 1]], r3[[All, 1]]]
        |>
    ],
    <|
        "exact non-normal" -> True, "machine non-normal" -> True, "ill-conditioned" -> True,
        "rotated number operator" -> True, "40 digits" -> True, "machine zeros" -> True,
        "wide spectra, any basis" -> {True, True, True, True},
        "ill-conditioned basis" -> {True, True, True}, "defective beside the null space" -> True,
        "cancellations" -> {True, "ZeroBasePowerNoLimit", True, True, True, True},
        "complex values written without I" -> {True, "ZeroBasePowerNoLimit", True, "ZeroBasePowerDefective"},
        "real parts of complex values" -> {True, "ZeroBasePowerNoLimit", "ZeroBasePowerNoLimit", True},
        "algebraic after FunctionExpand" -> {True, "ZeroBasePowerNoLimit"}, "zero on a branch cut" -> {True, True, True},
        "zero written with Log" -> True, "eight square roots" -> {True, "ZeroBasePowerNoLimit", "ZeroBasePowerNoLimit"},
        "symbolic beside a part read as zero" -> True, "inexact base that is not zero" -> True, "part read as zero" -> True,
        "huge entry" -> True, "algebraic spectrum" -> {True, True},
        "precision past 307 digits" -> {True, True}, "entry holding a zero" -> True
    |>,
    TestID -> "MatrixFunction-zero-base-inexact-zero-eigenvalue"
]

(* A symbolic diagonal entry t gives the limit in closed form, 1 at t == 0, 0 for
   Re[t] > 0 and Indeterminate on the half-plane where there is no limit: w n is the
   identity at w = 0, the vacuum projector at w = 1 and Indeterminate at w = -1. A
   symbolic matrix that is not diagonal gives the projector from its null spaces at
   generic values, carrying the condition that the other eigenvalues have positive
   real part: 0^(w M) for M with eigenvalues 0 and 2 is the null-space projector for
   every w > 0; the parametric Jordan block {{a, 1}, {0, a}} has no null space and
   gives 0 wherever Re a > 0; a block at 2 beside a symbolic eigenvalue w and a zero
   one gives diag(0, 0, 0, 1) at w = 3, and at w = 0, where the null space grows, the
   closed form has no value while the declared parameter reaches diag(0, 0, 1, 1). A
   symbolic block at zero, or a numeric eigenvalue -1 beside a symbolic one, fails
   with its reason. With w declared, values reach the limit directly, also through a
   further matrix function, and an indexed parameter substitutes into the closed
   form. A base that is identically zero but not written as 0 reads the same way, and
   so does a diagonal entry, GoldenRatio - (1 + Sqrt[5]) / 2. A numeric block beside
   w whose eigenvalue is 10^-100 is decided exactly, so the condition is on w alone. *)
VerificationTest[
    With[
        {
            m = {{1, 1}, {1, 1}},
            null = {{1, -1}, {-1, 1}} / 2,
            block = {{2, 1, 0, 0}, {0, 2, 0, 0}, {0, 0, mfW, 0}, {0, 0, 0, 0}},
            limit = Function[t, Piecewise[{{1, t == 0}, {0, Re[t] > 0}}, Indeterminate]]
        },
        {
            wn = mfZero[mfW DiagonalMatrix[{0, 1, 2}]],
            wm = 0 ^ QuantumOperator[mfW m, "Parameters" -> {mfW}],
            jordan = 0 ^ QuantumOperator[{{mfA, 1}, {0, mfA}}, "Parameters" -> {mfA}],
            blockClosed = mfZero[block],
            indexed = 0 ^ QuantumOperator[mfX[1] DiagonalMatrix[{0, 1, 2}], "Parameters" -> {mfX[1]}]
        },
        <|
            "closed form" -> wn === DiagonalMatrix[{1, limit[mfW], limit[2 mfW]}],
            "values of w" -> {wn /. mfW -> 0, wn /. mfW -> 1, wn /. mfW -> -1} ===
                {IdentityMatrix[3], DiagonalMatrix[{1, 0, 0}], DiagonalMatrix[{1, Indeterminate, Indeterminate}]},
            "non-diagonal" -> Simplify[mfZero[mfW m], mfW > 0] == null,
            "parametric" -> Normal[wm[1]["Matrix"]] == null && Normal[wm[0]["Matrix"]] == IdentityMatrix[2] && mfZeroBaseTag[wm[-1]] === "ZeroBasePowerNoLimit",
            "Jordan block" -> Normal[jordan["Matrix"]] === ConstantArray[Piecewise[{{0, Re[mfA] > 0}}, Indeterminate], {2, 2}] &&
                Normal[jordan[1]["Matrix"]] == ConstantArray[0, {2, 2}] && mfZeroBaseTag[jordan[0]] === "ZeroBasePowerDefective",
            "block beside w" -> (blockClosed /. mfW -> 3) == DiagonalMatrix[{0, 0, 0, 1}] && ! FreeQ[blockClosed /. mfW -> 0, Indeterminate] &&
                Normal[(0 ^ QuantumOperator[block, "Parameters" -> {mfW}])[0]["Matrix"]] == DiagonalMatrix[{0, 0, 1, 1}],
            "symbolic failures" -> {mfZeroBaseTag[0 ^ QuantumOperator[{{0, mfW}, {0, 0}}]], mfZeroBaseTag[0 ^ QuantumOperator[{{mfW, 1}, {0, -1}}]]} ==
                {"ZeroBasePowerDefective", "ZeroBasePowerNoLimit"},
            "nested" -> Normal[Sqrt[wm][1]["Matrix"]] == null,
            "indexed parameter" -> Normal[indexed[0]["Matrix"]] == IdentityMatrix[3] && FailureQ[indexed[-1]],
            "zero base" -> Normal[((Sin[mfTheta]^2 + Cos[mfTheta]^2 - 1) ^ QuantumOperator[DiagonalMatrix[{0, 1}]])["Matrix"]] == {{1, 0}, {0, 0}},
            "zero eigenvalue not written as 0" -> mfZero[DiagonalMatrix[{GoldenRatio, 3, 5/2}] - (1 + Sqrt[5]) / 2 IdentityMatrix[3]] == DiagonalMatrix[{1, 0, 0}],
            "numeric eigenvalue 10^-100 beside w" -> (mfZero[{{mfW, 1, 0, 0}, {0, 0, 0, 0}, {0, 0, 1, 1}, {0, 0, 1, 1 + 10^-100}}] /. mfW -> 1) ==
                {{0, -1, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}
        |>
    ],
    <|
        "closed form" -> True, "values of w" -> True, "non-diagonal" -> True, "parametric" -> True, "Jordan block" -> True,
        "block beside w" -> True, "symbolic failures" -> True, "nested" -> True, "indexed parameter" -> True, "zero base" -> True,
        "zero eigenvalue not written as 0" -> True, "numeric eigenvalue 10^-100 beside w" -> True
    |>,
    TestID -> "MatrixFunction-zero-base-symbolic"
]

EndTestSection[]
