(* A numeric function applied to an operator, f[qo], is the matrix function of the
   stored matrix. Every reference below is built independently of the route under
   test: either f of a spectrum whose eigenbasis is known by construction, or a
   closed form. The numeric operators are large enough, and their spectra clustered
   enough, that a route interpolating f by a polynomial in the matrix loses all
   accuracy: the Ising spectra here have 32 distinct eigenvalues, each twice
   degenerate under the global spin flip.

   Regimes covered:
     general symbolic     cos(t X), sin(t ZZ), t(XX + YY), [[a, b], [b, -a]],
                          cos(t sum_i Z_i) on 6 qubits, cos(t H) for the 4-qubit
                          Heisenberg chain, "Parameters", exact input
     exactly solvable     Jordan block, exact Log, Hadamard-rotated Ising, |A| of a
                          diagonalizable non-normal A, sqrt of a projector
     limiting             [[1, 1], [0, 1 + d]] for d = 10^-1 .. 10^-15 (separated to
                          defective), eigenvalue 10^-12 next to 1, coupling 10^-9
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

(* n qubits: cos(t sum_i Z_i) is diagonal with cos(t m), m the magnetization of each
   basis state. *)
VerificationTest[
    With[{n = 6},
        Simplify[
            Normal[Cos[mfT Total[QuantumOperator["Z", {#}] & /@ Range[n]]]["Matrix"]] -
                DiagonalMatrix[Cos[mfT Total[1 - 2 Tuples[{0, 1}, n], {2}]]]
        ]
    ],
    ConstantArray[0, {64, 64}],
    TimeConstraint -> 120,
    TestID -> "MatrixFunction-symbolic-n-qubit-magnetization"
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
   difference (sin a - sin b)/(a - b) cancels catastrophically at a - b = 10^-12
   (it was off by 1e-4); applying f to the substituted matrix is accurate. The
   reference is the exact matrix at 30 digits. *)
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

EndTestSection[]
