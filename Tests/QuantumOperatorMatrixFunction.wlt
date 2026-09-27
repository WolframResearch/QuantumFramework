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
                          values at and next to a collision of two eigenvalues *)

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

EndTestSection[]
