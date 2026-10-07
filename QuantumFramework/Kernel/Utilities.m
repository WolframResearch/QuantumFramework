Package["Wolfram`QuantumFramework`"]

PackageImport["Wolfram`Arrays`"]

PackageScope["ToList"]
PackageScope["SymbolicQ"]
PackageScope["basisMultiplicity"]
PackageScope["nameQ"]
PackageScope["nameString"]
PackageScope["propQ"]
PackageScope["propName"]
PackageScope["stateQ"]
PackageScope["orderQ"]
PackageScope["autoOrderQ"]
PackageScope["targetQ"]
PackageScope["targetsQ"]
PackageScope["measurementReprQ"]
PackageScope["emptyTensorQ"]

PackageScope["primeFactors"]
PackageScope["powerPrimeFactors"]

PackageScope["normalizeMatrix"]
PackageScope["tensorToVector"]
PackageScope["toTensor"]
PackageScope["tensorDimensions"]
PackageScope["tensorRank"]
PackageScope["identityMatrix"]
PackageScope["kroneckerProduct"]
PackageScope["projector"]
PackageScope["MatrixPartialTrace"]
PackageScope["blockDiagonalMatrix"]
PackageScope["eigenvalues"]
PackageScope["eigenvectors"]
PackageScope["eigensystem"]
PackageScope["roundoffChop"]
PackageScope["pauliMatrix"]
PackageScope["spinMatrix"]
PackageScope["fanoMatrix"]
PackageScope["GellMannMatrices"]
PackageScope["RegularSimplex"]
PackageScope["GramMatrix"]
PackageScope["GramDual"]

PackageScope["toggleSwap"]
PackageScope["toggleShift"]
PackageScope["alignDimensions"]

PackageScope["MatrixInverse"]
PackageScope["matrixFunction"]
PackageScope["zeroBasePower"]
PackageScope["diagonalMatrixQ"]
PackageScope["matrixExponential"]
PackageScope["valuelessEntriesQ"]
PackageScope["exactZeroQ"]
PackageScope["zeroBaseQ"]
PackageScope["SetPrecisionNumeric"]
PackageScope["TranscendentalRecognize"]

PackageScope["QuditAdjacencyMatrix"]

PackageScope["$QuantumFrameworkProfile"]
PackageScope["profile"]
PackageScope["cacheProperty"]



ToList = Developer`ToList


(* *)

SymbolicQ[expr_] := ! FreeQ[expr, sym_Symbol ? Developer`SymbolQ /; ! MemberQ[{Rational, Complex}, sym] && ! MemberQ[Attributes[sym], NumericFunction | Constant]]

basisMultiplicity[dim_, size_] := If[
    size === 1, 1,
    Replace[Ceiling @ Log[size, dim], Except[_Integer ? Positive] -> 1]
]


(* test functions *)

nameQ[name_] := MatchQ[name, _String | {_String, ___} | _String[___]]

nameString[name_] := Replace[name, s_String | {s_String, ___} | s_String[___] :> s]

propQ[prop_] := MatchQ[prop, _String | {_String, ___}]

propName[prop_] := Replace[prop, name_String | {name_String, ___} :> name]

(* The rank contract for a state's amplitude container: a vector of amplitudes
   or a square density matrix.  The shape comes from the container rather than
   from VectorQ and SquareMatrixQ, which answer for an explicit array only and
   would reject a lazy or symbolic one that knows its own shape perfectly
   well. *)

stateQ[state_] := MatchQ[ArrayDimensions[state], {_Integer} | {d_Integer, d_}]

orderQ[order_] := VectorQ[order, IntegerQ] && DuplicateFreeQ[order]

autoOrderQ[order_] := MatchQ[order, _ ? orderQ | Automatic | {_ ? orderQ | Automatic, _ ? orderQ | Automatic}]

targetQ[target_] := VectorQ[target, IntegerQ] && AllTrue[target, Positive]

targetsQ[targets_] := VectorQ[targets, targetQ]

measurementReprQ[state_] := TensorQ[state] && MemberQ[{2, 3}, tensorRank[state]]


(* numbers *)

primeFactors[n_] := Catenate[Table @@@ FactorInteger[n]]

powerPrimeFactors[n_] := Power @@@ FactorInteger[n]


(* Matrix tools *)

tensorToVector[t_ ? TensorQ] := Flatten[t]

(* scalar *)
tensorToVector[t_] := {t}


toTensor[t_ ? TensorQ] := t

toTensor[t_] := {t}


tensorDimensions[t_] := Replace[TensorDimensions[t], Except[_List] :> {}]

tensorRank[t_] := Length @ tensorDimensions[t]


identityMatrix[0 | {_, 0} | {0, _}] := {{}}

identityMatrix[n_] := IdentityMatrix[n, SparseArray]


normalizeMatrix[matrix_] := With[{tr = Tr[matrix]}, If[NumericQ[tr] && tr == 0, matrix, matrix / tr]]


kroneckerProduct[ts___] := Fold[If[ArrayQ[#1] && ArrayQ[#2], KroneckerProduct[##], Times[##]] &, {ts}]


projector[v_] := KroneckerProduct[v, Conjugate[v]]


MatrixPartialTrace[matrix_, trace_, dimensions_] := ArrayReshape[
    toTensor @ TensorContract[
        ArrayReshape[matrix, Join[dimensions, dimensions]], Thread[{trace, trace + Length[dimensions]}]
    ],
    Table[Times @@ Delete[dimensions, List /@ trace], 2]
]


emptyTensorQ[t_] := MatchQ[tensorDimensions[t], {___, 0}]

matrixQ[m_] := MatrixQ[m] && ! emptyTensorQ[m]


blockDiagonalMatrix[ms : {__ ? MatrixQ}] := MyBlockDiagonalMatrix[DeleteCases[ms, _ ? emptyTensorQ]]


(* The eigenvector v scaled to unit norm and, where the rules below apply, its phase
   fixed by first dividing v by one of its entries.

   An exact v with an entry that is not a Gaussian rational is divided by its first
   entry that is not zero, and each entry that is zero is written 0; a v whose entries
   are all zero comes back with its entries written 0. Such an entry can be zero without
   being written 0, as the first entry of an eigenvector of the Fourier transform can,
   and Unequal then leaves x != 0 undecided. An entry Unequal cannot tell from zero is
   zero when RootReduce writes it 0, for an algebraic entry, and when exactZeroQ reads it
   as zero otherwise; on the eigenvectors of the 7 x 7 Fourier transform the exact test
   exactZeroQ applies gives up and assumes zero, with a message. Normalize divides every
   entry by the norm of v, and for entries that are sums of roots of unity that norm is a
   nested radical about as large as the whole vector, so the normalized vector is many
   times the size of v. RootReduce writes the norm of a vector of algebraic numbers as a
   rational, a quadratic radical or one Root object, whose value N finds by isolating a
   root of an integer polynomial; Simplify can instead return a smaller expression whose
   terms cancel, which loses digits under machine N, and on some norms it runs for
   minutes. The vector is divided by that form of its norm when its entries are
   algebraic and this gives the smaller vector.

   Any other v, inexact, symbolic, or of Gaussian rationals, is divided by its first
   entry when that entry is a number Unequal finds nonzero. *)
normalizedEigenvector[v_] := If[
    VectorQ[v, NumericQ] && Precision[v] === Infinity && ! VectorQ[v, gaussianRationalQ],
    exactNormalizedEigenvector[Normal[v]],
    Normalize[If[NumericQ[First[v]] && First[v] != 0, v / First[v], v]]
]

exactNormalizedEigenvector[v_] := With[
    {z = Replace[v, _ ? zeroEntryQ -> 0, {1}]},
    {k = FirstPosition[z, Except[0], {0}, {1}, Heads -> False][[1]]},
    If[ k == 0,
        z,
        With[{w = z / z[[k]]},
            {plain = Normalize[w]},
            {reduced = If[AllTrue[w, TrueQ[Element[#, Algebraics]] &], w / RootReduce[Norm[w]], plain]},
            If[LeafCount[reduced] < LeafCount[plain], reduced, plain]
        ]
    ]
]

zeroEntryQ[x_] := ! TrueQ[x != 0] && If[TrueQ[Element[x, Algebraics]], RootReduce[x] === 0, exactZeroQ[x]]


Options[eigensystem] = {"Sort" -> False, "Normalize" -> False, "Orthogonalize" -> False, Chop -> False}

eigensystem[matrix_, OptionsPattern[]] := Module[{values, vectors},
    {values, vectors} = chopEigensystem[matrix, Simplify @ Enclose[
        ConfirmBy[
            If[ TrueQ[OptionValue[Chop]],
                Which[
                    MatrixQ[matrix, NumericQ] && Precision[matrix] === MachinePrecision,
                    machineEigensystem[N[Normal[matrix]]],
                    Precision[matrix] === MachinePrecision,
                    Quiet[
                        Check[
                            Eigensystem[matrix, ZeroTest -> (Chop[N[#1]] == 0 &)],
                            Eigensystem[matrix],
                            Eigensystem::eivec0
                        ],
                        Eigensystem::eivec0
                    ],
                    True,
                    Eigensystem[roundoffChop[matrix]]
                ],
                Eigensystem[matrix]
            ],
            MatchQ[{_ ? ListQ, _ ? MatrixQ}]
        ],
        Eigensystem[matrix] &
    ]];
    If[ ! MatchQ[OptionValue["Sort"], False | None] && AllTrue[values, NumericQ],
        With[{ordering = OrderingBy[values, Replace[OptionValue["Sort"],
                True | Automatic :> If[Length[values] > 2 && ContainsOnly[Arg[values], {0, Pi}], Identity, {Mod[Arg[#], 2 Pi], Abs[#]} &]
            ]]
        },
            values = values[[ordering]];
            vectors = vectors[[ordering]]
        ]
    ];
    If[ TrueQ[OptionValue["Normalize"]], vectors = normalizedEigenvector /@ vectors];
    (* Eigensystem returns an arbitrary (generally non-orthonormal) basis within a degenerate
       eigenspace, so the eigenvectors need not resolve the identity. Gram-Schmidt in the
       (sorted) eigenvalue order repairs each degenerate block; distinct-eigenvalue vectors of
       a Hermitian matrix are already orthogonal, so it leaves a well-conditioned basis
       essentially unchanged. Numeric bases only. *)
    If[ TrueQ[OptionValue["Orthogonalize"]] && ArrayQ[vectors, 2, NumericQ], vectors = Orthogonalize[vectors]];

    {values, vectors}
]

(* A machine matrix is diagonalized with no ZeroTest: an explicit ZeroTest takes
   it off the numerical eigensolvers, although the Eigensystem documentation
   limits that option to exact and symbolic matrices. Of the two solvers the
   default method chooses between, only the Hermitian one guarantees orthonormal
   eigenvectors inside a degenerate eigenspace, and it runs only when
   HermitianMatrixQ holds at its default tolerance, which roundoff alone can
   fail; so a matrix Hermitian to roundoff is made exactly Hermitian first. Any
   other matrix goes to the general solver as it is. *)
machineEigensystem[m_] := Eigensystem[If[nearlyHermitianQ[m, roundoffTolerance[m]], (m + ConjugateTranspose[m]) / 2, m]]

(* The eigenvalues and eigenvectors of the matrix m with each number below the roundoff of
   the eigensystem, relative to its scale, written 0. For an inexact numeric n x n matrix at
   precision p, with tol the smaller of 100 n 10^-p and 10^-10, that is an eigenvalue, or its
   real or imaginary part, smaller in magnitude than tol times the largest eigenvalue in
   magnitude, and the real or imaginary part of an eigenvector entry smaller than tol, against
   the unit norm of the eigenvectors Eigensystem returns for an inexact matrix (it pairs an
   extra eigenvalue of a defective one with a zero vector). The bound 10^-10 acts only below
   about 12 digits, where 100 n 10^-p reaches the digits the numbers carry: at 3 digits it
   would write 0 for the eigenvalue -0.372 of {{1, 2}, {3, 4}} and for the eigenvalue 1 of
   diag(1, 2, 3, 4). An eigenvalue below tol is not told from roundoff in a
   general basis, though a diagonal matrix gives it exactly: at machine precision
   diag(10^10, 10^-5) has the eigenvalues 10^10 and 0. A matrix that is roundoff throughout,
   as the commutator of two commuting matrices computed at machine precision, has no scale of
   its own to measure that roundoff against, and its eigenvalues are that roundoff. Any other
   matrix gets Chop's threshold, which leaves exact numbers as they are. *)
chopEigensystem[m_ ? inexactNumericArrayQ, {values_, vectors_}] := With[{tol = Min[roundoffTolerance[m], 10^-10]},
    {chopAtScale[values, tol, Max[Abs[values]]], chopAtScale[vectors, tol, 1.]}
]

chopEigensystem[_, es_] := Chop[es]

(* roundoffChop[a, n, ref] writes 0 for each approximate number in the numeric array a, or for
   its real or imaginary part, that roundoff alone could have left there, a being computed from
   ref, the entries of an n x n matrix at precision p: a number smaller in magnitude than tol
   times the largest entry of ref. tol is the smaller of 100 n 10^-p, what a computation at
   precision p leaves relative to the scale of such a matrix, with room to spare
   (roundoffTolerance), and 10^-10, a bound that acts only below about 12 digits, where
   100 n 10^-p reaches the digits the numbers carry: at 3 digits it would write 0 for the
   entry 1 of diag(1, 2, 3, 4). ref is a itself when it is not given, and n the number of rows
   of a matrix a when neither is. Chop's own threshold, 10^-10, does not scale with the
   numbers: it writes 0 for every entry of an array whose entries are all smaller, as those of
   hbar/2 sigma_z in SI units, and for the entries below it that a matrix known to 30 digits
   resolves. An array of exact numbers keeps them, and one that is not numeric, with no scale
   to measure roundoff against, gets Chop's threshold. *)
roundoffChop[a_ ? (ArrayQ[#, _, NumericQ] &), n_Integer, ref_ ? inexactNumericArrayQ] :=
    chopAtScale[a, Min[100 n roundoff[ref], 10^-10], Max[Abs[ref]]]

roundoffChop[a_, _Integer, _] := Chop[a]

roundoffChop[a_, n_Integer] := roundoffChop[a, n, a]

roundoffChop[m_] := roundoffChop[m, Length[m]]

(* A numeric array with an approximate number in it other than a zero known only to an
   accuracy, so that its entries have a precision and roundoff has a size (entryPrecision). *)
inexactNumericArrayQ[a_] := ArrayQ[a, _, NumericQ] && entryPrecision[a] < Infinity

(* Chop at tol times scale, as a machine number where it is one, since Chop compares machine
   numbers with a machine threshold faster. For a machine scale of at least 10^-290 and a tol
   of at least 10^-16, as machine precision gives, the product is formed in machine arithmetic
   and stays above the smallest machine number; otherwise it is formed at 20 digits, so that
   it does not underflow, and N is taken of it only inside the machine range. It is 0 only
   when the scale is, every number then being zero, and Chop at its own threshold writes those
   0. *)
chopAtScale[x_, tol_, scale_] := Which[
    ! TrueQ[scale > 0], Chop[x],
    MachineNumberQ[scale] && scale >= 1.*^-290 && tol >= 1.*^-16, Chop[x, N[tol] scale],
    True, With[{delta = SetPrecision[tol, 20] SetPrecision[scale, 20]},
        Chop[x, If[$MinMachineNumber <= delta <= $MaxMachineNumber, N[delta], delta]]
    ]
]

Options[eigenvalues] = Options[eigensystem]

eigenvalues[matrix_, opts : OptionsPattern[]] := First @ eigensystem[matrix, opts] 

Options[eigenvectors] = Options[eigensystem]

eigenvectors[matrix_, opts : OptionsPattern[]] := Last @ eigensystem[matrix, opts] 



pauliMatrix[n_] := pauliMatrix[n, 2]

spinMatrix[n_] := spinMatrix[n, 2]

pauliMatrix[0, dimension_] := identityMatrix[dimension]

spinMatrix[0, dimension_] := identityMatrix[dimension]

pauliMatrix[1, dimension_] := SparseArray[Table[{n, Mod[n - 1, dimension]} + 1 -> 1, {n, 0, dimension - 1}], {dimension, dimension}]

spinMatrix[1, dimension_] := With[{
    s = (dimension - 1) / 2
},
    SparseArray[
        {a_, b_} :> (KroneckerDelta[a, b + 1] + KroneckerDelta[a + 1, b]) Sqrt[(s + 1) (a + b - 1) - a b],
        {dimension, dimension}
    ]
]

pauliMatrix[2, dimension_] := - I pauliMatrix[3, dimension] . pauliMatrix[1, dimension]

spinMatrix[2, dimension_] := With[{
    s = (dimension - 1) / 2
},
    SparseArray[
        {a_, b_} :> I (KroneckerDelta[a, b + 1] -  KroneckerDelta[a + 1, b]) Sqrt[(s + 1) (a + b - 1) - a b],
        {dimension, dimension}
    ]
]

pauliMatrix[3, dimension_] := SparseArray[Table[{n, n} + 1 -> Exp[- 2 Pi I n / dimension], {n, 0, dimension - 1}], {dimension, dimension}]

spinMatrix[3, dimension_] := With[{
    s = (dimension - 1) / 2
},
    SparseArray[
        {a_, b_} :> 2 (s + 1 - a) KroneckerDelta[a, b],
        {dimension, dimension}
    ]
]

fanoMatrix[d_, q_, p_, x_ : Automatic, z_ : Automatic] :=
    Simplify @
    Exp[I Pi q p / d] *
        MatrixPower[FourierMatrix[d], 2] .
            MatrixPower[Replace[z, Automatic :> ConjugateTranspose[pauliMatrix[3, d]]], q] .
                    MatrixPower[Replace[x, Automatic :> ConjugateTranspose[pauliMatrix[1, d]]], p]


h[d_, k_] := Piecewise[{
    {IdentityMatrix[d, SparseArray], k == 1},
    {BlockDiagonalMatrix[{h[d - 1, k], {{0}}}], 1 < k < d},
    {Sqrt[2 / d / (d - 1)] BlockDiagonalMatrix[{IdentityMatrix[d - 1, SparseArray], {{1 - d}}}], k == d}
}]

f[d_, k_, j_] := Piecewise[{
    {SparseArray[{{k, j} -> 1, {j, k} -> 1}, {d, d}], k < j},
    {h[d, k], k == j},
    {SparseArray[{{k, j} -> I, {j, k} -> -I}, {d, d}], k > j}}
]

GellMannMatrix[n_Integer, i_Integer] /; n > 1 && 1 <= i <= n ^ 2 - 1 := Block[{k, l, j},
    k = Ceiling[Sqrt[i + 1]];
    l = i - (k - 1) ^ 2 + 1;
    j = Ceiling[l / 2];
    If[OddQ[l], f[n, j, k], f[n, k, j]]
]


GellMannMatrices[d_Integer ? Positive] := Table[GellMannMatrix[d, i], {i, 1, d ^ 2 - 1}]


RegularSimplex[d_Integer ? Positive] := 1 / Sqrt[d + 1] (Append[identityMatrix[d], ConstantArray[1 / d (1 + Sqrt[d + 1]), d]] - 1 / d (1 + 1 / Sqrt[d + 1]))


GramMatrix[x_] := Inverse[Outer[Tr @* Dot, x, x, 1]]

GramDual[x_] := GramMatrix[x] . x


(* optimization *)

(* An exactly unitary matrix, such as the Fourier basis matrix of roots of unity, inverts by
   its conjugate transpose; the general exact inverse of the same matrix does not finish at
   dimension 16. Machine and symbolic matrices keep the general inverse. *)
exactUnitaryMatrixQ[matrix_] := SquareMatrixQ[matrix] && Precision[matrix] === Infinity && MatrixQ[matrix, NumericQ] && UnitaryMatrixQ[matrix]

MatrixInverse[matrix_ ? exactUnitaryMatrixQ] := ConjugateTranspose[matrix]

MatrixInverse[matrix_] := If[
    SquareMatrixQ[matrix],
    Quiet[Check[Inverse[matrix], PseudoInverse[matrix], Inverse::sing], Inverse::sing],
    PseudoInverse[matrix]
]


(* A scalar function of a matrix. A diagonal matrix takes f entry by entry on its
   diagonal. An exact or symbolic matrix goes through its minimal polynomial
   (ComputeMatrixFunction), which is sound in exact arithmetic; its closed form is
   generic in any parameters, valid where the eigenvalues it separates stay
   distinct. In floating point that polynomial interpolation is ill-conditioned in
   the number of distinct eigenvalues, so an inexact matrix is diagonalized instead,
   and f of a diagonalized matrix needs no derivatives, so non-analytic f such as
   Abs works too: a normal matrix by its Schur decomposition m = q.t.q^†, whose t
   is diagonal to roundoff exactly when m is normal, and a non-normal one with a
   well-conditioned eigenbasis v by v.f(d).v^-1. Any other goes to MatrixFunction (Schur-Parlett), which needs f
   differentiable at the eigenvalues. An eigenvalue within the error bound of the
   decomposition from the negative real axis is put on it, so that Sqrt and Log take
   their principal value on the branch cut there, as for the exact matrix. The bound is
   the eigenvalue's own on the normal and eigenvector routes; on the Schur route of a
   nearly defective matrix it is one for the whole matrix, and it can take in an
   eigenvalue the decomposition has resolved. *)
matrixFunction[f : Plus | Minus | Times | Conjugate, mat_, {left___}, {right___}, ___] := f[left, mat, right]

(* The power M^p of an inexact matrix for a real exponent p that is not an Integer
   (inexactMatrixPower). MatrixPower keeps an Integer exponent, an exact or symbolic
   matrix, and the power applied to a vector. *)
matrixFunction[Power, mat_, {}, {p_}] /; realNonIntegerQ[p] && SquareMatrixQ[mat] && inexactMatrixQ[mat] :=
    inexactMatrixPower[mat, Re[p]]

matrixFunction[Power, mat_, {left___}, {right___}, opts : OptionsPattern[]] := MatrixPower[mat, left, right, opts]

matrixFunction[f_, mat_, {left___}, {right___}, opts : OptionsPattern[]] := scalarMatrixFunction[f[left, #, right] &, mat, opts]

(* A number with no imaginary part that is not an Integer, as 1/2, 0.5, 2. and Pi are. *)
realNonIntegerQ[p_] := NumericQ[p] && ! IntegerQ[p] && TrueQ[Im[p] == 0]

(* An exponent of integer value, as 2. or -1., gives the power MatrixPower computes for
   that Integer, from products and the inverse. *)
inexactMatrixPower[mat_, p_] /; TrueQ[p == Round[p]] := MatrixPower[mat, Round[p]]

(* For any other p, x^p has the branch cut of Sqrt and Log on the negative real axis, and
   M^p takes their route and their principal value there, unless the matrix has a Jordan
   block at zero. Such a block of size s needs the derivatives of x^p at 0 up to order
   s - 1: they exist, all 0, when p > s - 1, and one of them does not when p < s - 1. That
   route would take derivatives of x^p near zero, which grow without bound, on the
   eigenvalues roundoff splits apart there, so the block is read first, as 0^m reads it,
   with bound = zeroJordanBlockBound: for p < bound - 1 M^p is the Failure, as Sqrt is
   there, and for larger p it is MatrixPower's, off on the block by up to the order of
   ||m||^p just above bound - 1 and by far less well above it. The bound can
   exceed s, when a small nonzero eigenvalue lies as close to zero as the eigenvalues of a
   longer block would be spread, and M^p is then the Failure although it exists; and
   beside a Jordan block on the negative real axis, MatrixPower's M^p takes x^p there from
   both sides of the cut. *)
inexactMatrixPower[mat_, p_] := With[{bound = powerZeroBlockBound[mat]},
    Which[
        bound < 2, scalarMatrixFunction[# ^ p &, mat],
        p < bound - 1, derivativeFailure,
        True, MatrixPower[mat, p]
    ]
]

(* zeroJordanBlockBound of the inexact mat, 0 for a diagonal or Hermitian mat, which has no
   Jordan block. *)
powerZeroBlockBound[mat_] /; DiagonalMatrixQ[mat] := 0

powerZeroBlockBound[mat_] := With[{m = Normal[mat]},
    If[nearlyHermitianQ[m, roundoffTolerance[m]], 0, zeroJordanBlockBound[m, roundoffZeroEigenvalue[m, Eigenvalues[m]]]]
]

(* f must have a finite value at each numeric eigenvalue: Log at a zero eigenvalue is a
   Failure, not a matrix of infinities. At an inexact eigenvalue the value must pass
   numberValueQ, not only NumericQ: a function defined on real arguments only, as CubeRoot
   or Surd, stays unevaluated at a complex eigenvalue, CubeRoot[-1. + 0.5 I], and NumericQ
   of that is True. A caller that does not use the values, and asks only that they be
   finite, gives valueQ = NumericQ. *)
spectralValues[f_, eigenvalues_, valueQ_ : numberValueQ] := With[{values = f /@ eigenvalues},
    If[ ! VectorQ[eigenvalues, NumericQ] || And @@ MapThread[If[InexactNumberQ[#1], valueQ[#2], NumericQ[#2]] &, {eigenvalues, values}],
        values,
        Failure["NonFiniteMatrixFunction", <|"MessageTemplate" -> "The function has no finite value at an eigenvalue."|>]
    ]
]

(* A value of f at an inexact argument: a number, or an exact numeric value, as
   Arg[-1.] = Pi; not an expression left unevaluated at the argument, as
   CubeRoot[-1. + 0.5 I], whose precision is that of the argument. *)
numberValueQ[v_] := NumberQ[v] || NumericQ[v] && Precision[v] === Infinity

(* An inexact diagonal is held to the same roundoff rules as a dense matrix, so f
   commutes with a change of basis near a zero eigenvalue and near the negative real
   axis too. *)
scalarMatrixFunction[f_, mat_, ___] /; SquareMatrixQ[mat] && DiagonalMatrixQ[mat] := Enclose @ With[
    {eigenvalues = Normal[Diagonal[mat]]},
    SparseArray[
        Band[{1, 1}] -> Confirm[spectralValues[f,
            If[ MatrixQ[mat, NumericQ] && Precision[mat] < Infinity,
                onNegativeAxis[roundoffEigenvalues[eigenvalues, roundoff[mat]], 16 roundoff[mat] Max[Abs[eigenvalues]]],
                eigenvalues
            ]
        ]],
        Dimensions[mat]
    ]
]

scalarMatrixFunction[f_, mat_, opts___] /; SquareMatrixQ[mat] && MatrixQ[mat, NumericQ] && Precision[mat] < Infinity :=
    With[{m = Normal[mat]},
        inexactMatrixFunction[f, m, roundoff[m], roundoffTolerance[m], opts]
    ]

(* A derivative at numeric arguments with no numeric value, Derivative[1][Abs][1] at
   a Jordan block, means f is not differentiable where the matrix needs it; one that
   is merely unevaluated, Derivative[1][Zeta][2], has a value. *)
scalarMatrixFunction[f_, mat_, opts___] := With[{values = ResourceFunction["ComputeMatrixFunction"][f, mat, opts]},
    If[! valuelessEntriesQ[values] && FreeQ[values, (d : Derivative[__][_][__ ? NumericQ] /; ! NumericQ[N[d]])], values, derivativeFailure]
]

(* An entry holds Indeterminate, Undefined or an infinity somewhere other than the
   default of a Piecewise. A Piecewise whose symbols have not chosen a case, as in
   the closed form of 0^m, has a value; once values send it to its default, the
   default is the entry, and the check sees it. A SparseArray is atomic, so its
   values are read out. *)
valuelessEntriesQ[array_SparseArray] := valuelessEntriesQ[Append[array["NonzeroValues"], array["Background"]]]

valuelessEntriesQ[array_] := ! FreeQ[array /. HoldPattern[Piecewise[cases_, _]] :> Piecewise[cases], Indeterminate | Undefined | _DirectedInfinity]

(* eps is the relative precision of the entries; within tol = 100 n eps the matrix
   counts as Hermitian and the strict upper triangle of t as zero. An eigenvalue
   within 10 eps of zero, relative to the largest, is set to zero: roundoff alone
   puts it there, and f may be singular at zero (Sqrt, Log). *)
roundoffEigenvalues[eigenvalues_, eps_] := Chop[eigenvalues, 10 eps Max[Abs[eigenvalues]]]

(* The roundoff a decomposition of an n x n matrix at precision p leaves, relative to
   its norm, with room to spare: 100 n 10^-p. *)
roundoffTolerance[mat_] := 100 Length[mat] roundoff[mat]

(* The error a Schur decomposition of an n x n matrix at precision p leaves in q.t.q^† - m,
   with room to spare: 10 n 10^-p times the scale of the matrix, a tenth of
   roundoffTolerance. *)
backwardError[mat_, scale_] := 10 Length[mat] roundoff[mat] scale

(* The unit roundoff 10^-p of a matrix of precision p, computed at 20 digits so that it
   does not underflow for p above 307; the machine value is computed once. *)
roundoff[mat_] := roundoffAt[entryPrecision[mat]]

roundoffAt[MachinePrecision] = 10 ^ -SetPrecision[MachinePrecision, 20]

roundoffAt[p_] := 10 ^ -SetPrecision[p, 20]

(* The precision of the entries of mat other than zeros known to an accuracy, Infinity
   when there are none. A zero known to an accuracy, as 0``40, has precision 0, which
   says nothing about the roundoff in the other entries; a machine 0. does. The precision
   of mat is either 0 or already that of its other entries, so they are gone through one
   at a time only when it is 0, which for a large matrix takes longer than its
   eigensystem. *)
entryPrecision[mat_] := With[{p = Precision[mat]}, If[p > 0, p, nonzeroEntryPrecision[mat]]]

nonzeroEntryPrecision[mat_SparseArray] := Precision[Select[Append[mat["NonzeroValues"], mat["Background"]], ! accuracyZeroQ[#] &]]

nonzeroEntryPrecision[mat_] := Precision[Select[Flatten[mat], ! accuracyZeroQ[#] &]]

accuracyZeroQ[x_] := TrueQ[x == 0] && Precision[x] == 0

(* The precision of a result computed from the inexact mat: that of its nonzero entries,
   and machine precision when it has none. *)
resultPrecision[mat_] := Replace[entryPrecision[mat], Infinity -> MachinePrecision]

(* HermitianMatrixQ's Tolerance t sets entries below t to zero and compares the others
   to relative precision t one entry at a time, which does not bound m - m^† against
   the largest entry, so the test is written out: the largest entry of m - m^† against
   tol times the largest absolute entry s. For s below 1 it divides by s instead, so
   that the bound does not underflow for a matrix near the smallest machine numbers. *)
nearlyHermitianQ[mat_, tol_] := With[{s = Max[Abs[mat]], d = Max[Abs[mat - ConjugateTranspose[mat]]]},
    s == 0 || If[s < 1, d / s <= tol, d <= tol s]
]

(* The Schur factor q is unitary even inside a degenerate eigenspace, where the
   eigenvectors Eigensystem returns need not be orthonormal, so f(m) = q.f(t).q^†
   when t is diagonal to roundoff. q.(x q^†) is q.DiagonalMatrix[x].q^† without the
   dense diagonal product. A Hermitian matrix keeps real eigenvalues and, for real
   f values, gives an exactly Hermitian result, real for real input. The eigenvalue
   t_ii of a normal matrix is known to within the residual of its Schur vector. *)
inexactMatrixFunction[f_, mat_, eps_, tol_, opts___] := Enclose @ With[
    {qt = SchurDecomposition[mat, RealBlockDiagonalForm -> False]},
    {q = First[qt], t = Last[qt], hermitianQ = nearlyHermitianQ[mat, tol]},
    If[ hermitianQ || Max[Abs[UpperTriangularize[t, 1]]] <= tol Max[Abs[t]],
        With[
            {eigenvalues = With[{e = roundoffEigenvalues[If[hermitianQ, Re, Identity][Diagonal[t]], eps]},
                If[hermitianQ, e, onNegativeAxis[e, eigenvalueError[mat, q, Diagonal[t], Max[Abs[t]]]]]
            ]},
            {values = Confirm[spectralValues[f, eigenvalues]]},
            {result = q . (values ConjugateTranspose[q])},
            Which[
                ! hermitianQ, realOnRealInput[mat, eigenvalues, values, tol, result],
                ! FreeQ[values, _Complex], result,
                FreeQ[mat, _Complex], Re[(result + Transpose[result]) / 2],
                True, (result + ConjugateTranspose[result]) / 2
            ]
        ],
        nonNormalMatrixFunction[f, mat, qt, Eigensystem[mat], eps, opts]
    ]
]

(* An eigenvalue x left of the imaginary axis whose distance |Im x| from the negative
   real axis, where Sqrt and Log have their branch cut, is within bound (one number, or
   one per eigenvalue), the error the decomposition may have made in it, is put on that
   axis: which side of the cut it landed on may be roundoff. f of the real number is the
   principal value, the value the exact matrix gives when x is real. *)
onNegativeAxis[eigenvalues_, bound_] := MapThread[
    If[Re[#1] < 0 && Abs[Im[#1]] <= #2, Re[#1], #1] &,
    {eigenvalues, If[ListQ[bound], bound, ConstantArray[bound, Length[eigenvalues]]]}
]

(* How far each computed eigenvalue x_i with the unit vector v_i (the columns of vectors)
   can lie from an eigenvalue of m, with room to spare: four times its residual
   ||m.v_i - x_i v_i|| plus 16 10^-p scale for the roundoff in computing it. That bounds
   the error for a normal m; for any other it is multiplied by the condition number of
   the eigenvalue. *)
eigenvalueError[mat_, vectors_, eigenvalues_, scale_] :=
    4 (Norm /@ Transpose[mat . vectors - (# eigenvalues &) /@ vectors]) + 16 roundoff[mat] scale

(* A real matrix has eigenvalues in conjugate pairs, and f of it is real when f takes
   each pair to conjugate values, as Cos and Exp do, and Log off the negative axis.
   Each eigenvalue is paired with the one nearest its conjugate, which is itself when
   it is real, so that f of a real eigenvalue must be real; the imaginary part of the
   result is then roundoff and is dropped. *)
realOnRealInput[mat_, eigenvalues_, values_, tol_, result_] := If[
    FreeQ[mat, _Complex] && VectorQ[values, NumericQ] &&
        Max[Abs[values[[Nearest[eigenvalues -> "Index", Conjugate[eigenvalues]][[All, 1]]]] - Conjugate[values]]] <=
            tol Max[1, Max[Abs[values]]],
    Re[result],
    result
]

(* An eigenbasis v whose condition number stays below eps^(-1/4) keeps the error of
   v.f(d).v^-1 near eps^(3/4). The rows of v^-1 are the left eigenvectors scaled to the
   unit right ones, and the length of a row is the condition number of its eigenvalue:
   to first order, a perturbation of size e moves the eigenvalue by up to e times it. A
   defective or nearly defective matrix fails the condition and goes to the Schur form
   (cutMatrixFunction). *)
nonNormalMatrixFunction[f_, mat_, _, {eigenvalues_, vectors_}, eps_, ___] /;
    With[{sv = SingularValueList[vectors]}, Length[sv] == Length[vectors] && Max[sv] <= eps ^ (-1/4) Min[sv]] :=
    Enclose @ With[
        {inverse = Inverse[Transpose[vectors]], tol = roundoffTolerance[mat]},
        {placed = onNegativeAxis[eigenvalues, eigenvalueError[mat, Transpose[vectors], eigenvalues, Norm[Flatten[mat]]] (Norm /@ inverse)]},
        {values = Confirm[spectralValues[f, placed]]},
        realOnRealInput[mat, placed, values, tol, Transpose[vectors] . (values inverse)]
    ]

(* MatrixFunction takes divided differences of f between eigenvalues that roundoff has
   split apart around zero, and those stay finite where f has no derivative at zero, so
   Sqrt and Log of a matrix within roundoff of a nilpotent one can come out as large
   numbers that mean nothing. So the eigenvalue 0 is read to within roundoff as 0^m reads
   it, and f fails when it has no finite value there or, as its series shows, no finite
   derivative of an order below a bound on the size of the largest Jordan block there.
   Otherwise the Schur form goes to cutMatrixFunction. Of f's values at these eigenvalues
   only finiteness is asked: MatrixFunction takes f at eigenvalues it computes itself, and
   builtinMatrixFunction guards those calls; a Failure it returns is the result. *)
nonNormalMatrixFunction[f_, mat_, qt_, {eigenvalues_, _}, _, opts___] := Enclose[
    Confirm[spectralValues[f, eigenvalues, NumericQ]];
    If[ lacksDerivativeAtZeroQ[f, zeroJordanBlockBound[mat, roundoffZeroEigenvalue[mat, eigenvalues]]],
        derivativeFailure,
        Replace[cutMatrixFunction[f, mat, qt, opts], Except[_Failure | _ ? (MatrixQ[#, numberValueQ] &)] -> derivativeFailure]
    ]
]

(* f of the Schur form q.t.q^† of a defective or nearly defective m. Roundoff splits an
   eigenvalue with a Jordan block into several about it, and on the negative real axis
   they land on both sides of the cut of Sqrt and Log. Taken there, f comes from both
   sides, and the divided differences MatrixFunction takes between them, of order one
   over the split, mean nothing. Each eigenvalue that roundoff may have moved off the
   axis, at a point where f jumps across it (cutMove), is put back on it, at its real
   part. An imaginary part within the backward error is dropped; a larger one is restored
   by expanding f(t) about that t' in the imaginary parts d, since f(t) can be sensitive
   to it far beyond its size, as when a genuine eigenvalue lies close by: the blocks
   along the top row of f of the block matrix with t' on the diagonal and
   DiagonalMatrix[d] above it are the terms (1/j!) d^j/ds^j f(t' + s d) at s = 0. f is
   differentiated only on the axis, so every eigenvalue split from one on the cut takes
   the value from above the cut, as the exact matrix does. The terms shrink at least as
   fast as (|d| / |Re x|)^j; a matrix that would need more than four terms at that rate,
   or has nothing to move, goes to builtinMatrixFunction unchanged, and so does one that
   leaves beside the restored eigenvalues another whose side of the cut the decomposition
   has not determined (undecidedQ). *)
cutMatrixFunction[f_, mat_, {q_, t_}, opts___] := With[
    {eigenvalues = Diagonal[t], scale = Norm[Flatten[t]]},
    {bound = backwardError[mat, scale], mirrors = If[FreeQ[mat, _Complex], Nearest[eigenvalues -> "Index", Conjugate[eigenvalues]], None]},
    {moves = MapIndexed[cutMove[f, t, #1, First[#2], bound, scale, mirrors, roundoff[mat]] &, eigenvalues]},
    {placed = MapThread[If[#2 > 0, Re[#1], #1] &, {eigenvalues, moves}]},
    {displacement = MapThread[If[#3 == 2, #1 - #2, 0] &, {eigenvalues, placed, moves}]},
    {order = If[Max[Abs[displacement]] == 0, 0,
        Ceiling[Log[roundoff[mat]] / Log[Max[MapThread[If[#3 == 2, Abs[#1] / Abs[#2], 0] &, {displacement, placed, moves}]]]]
    ]},
    If[ Max[moves] == 0 || order > 4 ||
            MemberQ[moves, 2] && AnyTrue[Pick[eigenvalues, moves, 0], undecidedQ[f, t, #, bound, scale, roundoff[mat]] &],
        realResult[mat, builtinMatrixFunction[f, mat, opts]],
        seriesMatrixFunction[f, mat, q, t + DiagonalMatrix[placed - eigenvalues], displacement, order, opts]
    ]
]

(* Whether f may need its value from above the cut at the eigenvalue x: Re x < 0, Im x not
   zero, since a real eigenvalue is on the axis already, and at most half the distance
   -Re x from the branch point 0, so that an expansion about Re x reaches it, and f jumps
   across the axis at Re x (jumpQ), which Cos and Exp do not. *)
nearCutQ[f_, x_, eps_] := Re[x] < 0 && Im[x] != 0 && Abs[Im[x]] <= - Re[x] / 2 && jumpQ[f, Re[x], eps]

(* Whether the eigenvalue x near the cut, left where it is by cutMove, may nevertheless lie
   on the axis: it is farther from it than any split of a Jordan block of order up to 4,
   yet a perturbation of t within the backward error puts an eigenvalue at the point
   halfway from x to the axis. Its side of the cut is then not determined, and no series
   restores it. *)
undecidedQ[f_, t_, x_, bound_, scale_, eps_] := nearCutQ[f, x, eps] && Abs[Im[x]] > 2 (bound / scale) ^ (1 / 4) scale &&
    halfwaySingularValue[t, x] <= bound

(* The smallest singular value of t at the point halfway from x to the real axis. *)
halfwaySingularValue[t_, x_] := First[SingularValueList[t - (Re[x] + I Im[x] / 2) IdentityMatrix[Length[t]], -1, Tolerance -> 0]]

(* How the eigenvalue x of t, the i-th, is put on the negative real axis: 0 when it stays,
   1 when it moves to Re x, 2 when it moves and its imaginary part is restored by the
   series. It stays unless it is near the cut (nearCutQ). It just moves when Im x is
   within the backward error. It moves and is restored when, for a real m (mirrors
   lists, for each eigenvalue, those nearest its conjugate), x is its own mirror image,
   which in a real matrix only an eigenvalue on the real axis, or split from one there,
   can be; or when a perturbation of t within the backward error puts an eigenvalue at
   the point halfway from x to the axis: the smallest singular value of t there is at
   most bound. The pair a +- b i of a real m is each other's mirror image and passes only
   that last test, which a pair farther apart than roundoff can split fails. The singular
   value is computed only within the widest split roundoff gives a Jordan block of order
   up to 4, 2 (bound / scale)^(1/4) scale. *)
cutMove[f_, t_, x_, i_, bound_, scale_, mirrors_, eps_] := Which[
    ! nearCutQ[f, x, eps],
        0,
    Abs[Im[x]] <= bound,
        1,
    ListQ[mirrors] && MemberQ[mirrors[[i]], i] ||
        Abs[Im[x]] <= 2 (bound / scale) ^ (1 / 4) scale && halfwaySingularValue[t, x] <= bound,
        2,
    True,
        0
]

(* Whether f is discontinuous across the real axis at a, where it has a derivative for the
   expansion to use (which Arg does not): its values a distance h above and below a differ
   by far more than the 2 h |f'(a)| a smooth f moves over the gap, and by more than
   roundoff in their size. Values that are not numbers, as CubeRoot of a complex number,
   count as no jump. *)
jumpQ[f_, a_, eps_] := With[{h = Sqrt[eps] Abs[a]}, {above = f[a + I h], below = f[a - I h], slope = Derivative[1][f][a]},
    NumberQ[above] && NumberQ[below] && NumberQ[slope] &&
        Abs[above - below] > 10 ^ 4 2 h Abs[slope] + 10 ^ -4 (Abs[above] + Abs[below])
]

(* f(t' + d) for the Schur factor t' and the displacements d, through the terms of
   order up to order in d, read off the top block row of f of the block matrix with t'
   on the diagonal and DiagonalMatrix[d] above it, and returned in the frame q; at
   order 0, when nothing is restored, it is f(t'). The series has converged when its
   last term is within roundoff of the sum. Otherwise the order doubles, up to 8: an
   eigenvalue near the moved ones slows the series. One that has not converged by then,
   or that MatrixFunction cannot compute, gives way to builtinMatrixFunction of m. *)
seriesMatrixFunction[f_, mat_, q_, t_, displacement_, order_, opts___] := With[
    {n = Length[t]},
    {blocks = MatrixFunction[f,
        KroneckerProduct[IdentityMatrix[order + 1], t] +
            If[order > 0, KroneckerProduct[DiagonalMatrix[ConstantArray[1, order], 1], DiagonalMatrix[displacement]], 0],
        opts
    ]},
    {terms = If[MatrixQ[blocks, NumericQ], ArrayReshape[blocks[[;; n]], {n, order + 1, n}], None]},
    Which[
        terms === None,
            realResult[mat, builtinMatrixFunction[f, mat, opts]],
        order == 0 || Max[Abs[terms[[All, -1]]]] <= roundoffTolerance[mat] Max[Abs[Total[terms, {2}]]],
            realResult[mat, q . Total[terms, {2}] . ConjugateTranspose[q]],
        order < 8,
            seriesMatrixFunction[f, mat, q, t, displacement, Min[2 order, 8], opts],
        True,
            realResult[mat, builtinMatrixFunction[f, mat, opts]]
    ]
]

(* f of a real matrix whose imaginary part is roundoff is real. *)
realResult[mat_, result_] := If[
    FreeQ[mat, _Complex] && MatrixQ[result, NumericQ] && Max[Abs[Im[result]]] <= roundoffTolerance[mat] Max[Abs[result]],
    Re[result],
    result
]

(* The built-in MatrixFunction. Its default method for an inexact matrix, Schur-Parlett,
   takes f at the mean of each cluster of close eigenvalues of its own complex Schur form:
   at complex numbers at machine precision, x + 0. I for a real eigenvalue x, and above it
   at the eigenvalues that roundoff moves off the real axis. A function defined on real
   arguments only, as CubeRoot, Surd or RealAbs, stays unevaluated at a complex number, and
   MatrixFunction then crashes the kernel at machine precision and, above it, returns a
   matrix typically off by order 1. So f is first taken once at each distinct entry of the
   diagonal of SchurDecomposition's complex form, read as complex at machine precision
   when the matrix is to be balanced first, whose Schur form is then another one, and the
   result is a Failure that names the entry farthest from the real axis at which f has no
   value. That diagonal can be complex where the cluster means are real, as for a complex
   matrix with a real eigenvalue, so the check can fail where MatrixFunction would have
   been right. Above machine precision MatrixFunction also returns huge numbers for exact
   values of f, as RealSign[-1.`30] = -1, which inexactValue gives at the precision of the
   argument. Another method, or an option MatrixFunction does not know, is left to
   MatrixFunction, and a result with an entry that is not a value is the Failure. *)
builtinMatrixFunction[f_, mat_, opts___] := With[
    {method = OptionValue[MatrixFunction, FilterRules[{opts}, Options[MatrixFunction]], Method]},
    {points = If[
        MatchQ[method, Automatic | "Schur" | {Automatic | "Schur", ___}] && FilterRules[{opts}, Except[Options[MatrixFunction]]] === {},
        DeleteDuplicates @ SortBy[
            Diagonal[Last[SchurDecomposition[mat, RealBlockDiagonalForm -> False]]] +
                If[MemberQ[method, (("Balanced" -> v_) | ("Balanced" :> v_)) /; TrueQ[v]], 0. I, 0],
            - Abs[Im[#]] &
        ],
        {}
    ]},
    {checked = Reap[SelectFirst[points, With[{v = f[#]}, Sow[v, builtinMatrixFunction]; ! numberValueQ[v]] &], builtinMatrixFunction]},
    Replace[First[checked], {
        _Missing :> Replace[
            MatrixFunction[
                If[Precision[mat] =!= MachinePrecision && MemberQ[Catenate[Last[checked]], v_ /; Precision[v] === Infinity], inexactValue[f, #] &, f],
                mat,
                opts
            ],
            Except[_ ? (MatrixQ[#, numberValueQ] &)] -> derivativeFailure
        ],
        z_ :> complexValueFailure[z]
    }]
]

(* f at an inexact argument, with an exact value given by N at the precision of the
   argument, which leaves 0 exact; at any other argument, a symbolic one included, f
   itself, so that MatrixFunction takes f's derivatives. *)
inexactValue[f_, x_ ? InexactNumberQ] := With[{v = f[x]}, If[Precision[v] === Infinity, N[v, Precision[x]], v]]

inexactValue[f_, x_] := f[x]

(* A bound on the size of the largest Jordan block of the eigenvalue 0 of the inexact m to
   within roundoff, given zero = roundoffZeroEigenvalue[m, eigenvalues]: 0 when 0 is not
   an eigenvalue, 1 when it is a semisimple one, and when it is defective a - k + 1, and at
   least 2, for its k null vectors and the largest a for which the a eigenvalues of
   smallest modulus lie within tol^(1/a) ||m|| of zero, as far as a perturbation of size
   tol ||m|| spreads the eigenvalues of a Jordan block of size a. Small eigenvalues that
   roundoff cannot bring to zero can count in a too, so the bound can exceed the size. *)
zeroJordanBlockBound[m_, zero_] := Which[
    zero["Defective"], Max[2, Max[0, Select[Range[Length[m]], Abs[zero["Eigenvalues"][[#]]] <= zero["Tolerance"] ^ (1 / #) zero["Norm"] &]] - zero["Nullity"] + 1],
    zero["Nullity"] > 0, 1,
    True, 0
]

(* True when f has no finite value at zero, for s at least 1, or, for s above 1, its series
   at zero shows that it has no finite derivative there of some order below s: a series
   that comes back as a SeriesData with a nonzero term of power at most s - 1 whose
   coefficient is not numeric (a logarithm) or whose power is not a nonnegative integer (a
   pole, or a branch point as Sqrt has), or as a Piecewise, an expansion that depends on
   the side of zero (Surd, RealAbs). f is read with exact parameters, so BesselJ[1., x]
   expands as BesselJ[1, x] does, and its value is taken at an exact 0, so where it has
   none the kernel may say so, as for 1/x; the series may warn as well, as for
   AlternatingFactorial. A derivative whose formula is 0/0 at zero, as that of
   SinhIntegral, is finite. A series of any other shape (Abs and Round do not expand,
   LogIntegral comes back as a multiple of a SeriesData), or one not read within a second
   (DirichletEta at a block of size 64), leaves the decision to MatrixFunction. *)
lacksDerivativeAtZeroQ[_, s_] /; s < 1 := False

lacksDerivativeAtZeroQ[f_, s_] := With[{exact = SetPrecision[f, Infinity]},
    ! NumericQ[N[exact[0]]] || s > 1 && TimeConstrained[nonTaylorQ[Series[exact[\[FormalX]], {\[FormalX], 0, s - 1}], s - 1], 1, False]
]

(* SeriesData[x, 0, c, nmin, nmax, den] is the sum of c[[i]] x^((nmin + i - 1)/den) plus
   terms of order x^(nmax/den). A coefficient is tested for zero last, and only on a term
   that is not a Taylor term, since that test can take seconds. *)
nonTaylorQ[series_SeriesData, order_] := With[{c = series[[3]], nmin = series[[4]], den = series[[6]]},
    AnyTrue[
        Transpose[{c, (nmin + Range[0, Length[c] - 1]) / den}],
        #[[2]] <= order && ! (NumericQ[#[[1]]] && IntegerQ[#[[2]]] && #[[2]] >= 0) && ! TrueQ[#[[1]] == 0] &
    ]
]

nonTaylorQ[_Piecewise, _] := True

nonTaylorQ[_, _] := False

derivativeFailure = Failure["NonFiniteMatrixFunction", <|
    "MessageTemplate" -> "The function has no finite value, or no derivative where a Jordan block of the matrix needs one, at an eigenvalue."
|>]

complexValueFailure[z_] := Failure["NonFiniteMatrixFunction", <|
    "MessageTemplate" -> "The function has no finite value at the eigenvalue `1` of the complex Schur form, from which the matrix function of a matrix without a well-conditioned eigenbasis is computed.",
    "MessageParameters" -> {z}
|>]

(* 0^m is the limit of b^m = MatrixExp[Log[b] m] as b -> 0. On an eigenvalue t the
   limit of b^t is 1 at t = 0 and 0 for Re t > 0, and there is none for any other t.
   A Jordan block adds the terms Log[b]^k b^t, which vanish for Re t > 0 and diverge
   at t = 0. So the limit exists exactly when every eigenvalue is zero or has
   positive real part and the zero eigenvalue is semisimple, and it is then the
   spectral projector onto the null space of m along its range. For a number
   operator it is the vacuum projector, for H - E0 the projector onto the ground
   space, and for minus a Lindblad generator L the limit of exp(t L) as t -> Infinity,
   the projector onto its steady states. The kind of matrix picks the computation. *)
zeroBasePower[mat_ ? SquareMatrixQ] /; DiagonalMatrixQ[mat] && ! inexactMatrixQ[mat] := diagonalZeroBasePower[mat]

zeroBasePower[mat_ ? SquareMatrixQ] /; DiagonalMatrixQ[mat] && inexactMatrixQ[mat] := inexactDiagonalZeroBasePower[mat]

zeroBasePower[mat_ ? SquareMatrixQ] /; ! DiagonalMatrixQ[mat] && ! MatrixQ[mat, NumericQ] := symbolicZeroBasePower[zeroesWritten[Normal[mat]]]

zeroBasePower[mat_ ? SquareMatrixQ] /; ! DiagonalMatrixQ[mat] && inexactMatrixQ[mat] := inexactZeroBasePower[Normal[mat]]

zeroBasePower[mat_ ? SquareMatrixQ] /; ! DiagonalMatrixQ[mat] && MatrixQ[mat, NumericQ] && ! inexactMatrixQ[mat] :=
    With[{m = zeroesWritten[Normal[mat]]}, Which[
        DiagonalMatrixQ[m], diagonalZeroBasePower[m],
        exactHermitianQ[m] && FreeQ[settledValue[#, Abs] & /@ Select[DeleteDuplicates[Flatten[m]], ! gaussianRationalQ[#] &], Missing["Overflow" | "Underflow" | "NotANumber"]], exactHermitianZeroBasePower[m],
        True, exactZeroBasePower[m]
    ]]

inexactMatrixQ[mat_] := MatrixQ[mat, NumericQ] && Precision[mat] < Infinity

(* m with every numeric part that is zero but not written as 0 written as 0, decided as
   exactZeroQ decides, so that a number the convention there reads as zero, as
   Exp[-10^8], is zero throughout what follows; each distinct part is decided once, and
   a part that is zero is written as 0 whole. An entry takes the rewrite only when
   what it changes by is a number that is zero too: in -10^4344 Exp[-10^4], about
   -11.35, the factor Exp[-10^4] stays, so the entry keeps its value in any basis, and
   so it does in 10^4344 Exp[-10^4] w for a symbol w. *)
zeroesWritten[m_] := With[
    {zero = AssociationMap[exactZeroQ, DeleteDuplicates[Cases[m, x_ /; ! AtomQ[x] && NumericQ[x], {2, Infinity}]]]},
    Map[
        With[{w = # /. x_ /; TrueQ[zero[x]] :> 0}, If[w === # || NumericQ[# - w] && exactZeroQ[# - w], w, #]] &,
        m,
        {2}
    ]
]

(* An exact m is Hermitian when m minus its conjugate transpose is zero entry by entry,
   decided exactly, as the real part of m_ij - m_ji and the imaginary part of
   m_ij + m_ji for i <= j. (HermitianMatrixQ compares the entries of a SparseArray
   numerically, and calls a matrix with long algebraic entries Hermitian that is
   not; Conjugate decides realness numerically, and hits the precision limit on an
   entry that holds a zero on a branch cut.) *)
exactHermitianQ[m_ ? (MatrixQ[#, gaussianRationalQ] &)] := m === ConjugateTranspose[m]

exactHermitianQ[m_] := AllTrue[
    Flatten[{UpperTriangularize[m - Transpose[m]], UpperTriangularize[I (m + Transpose[m])]}],
    realPartSign[#] === 0 &
]

(* A diagonal m that is exact or symbolic: each entry t is an eigenvalue with its own
   null vector. The limit of b^t is 1 where t is zero and 0 where t is a number of
   positive real part, both decided exactly however t is written, and for a symbolic t
   it is the closed form Piecewise[{{1, t == 0}, {0, Re[t] > 0}}, Indeterminate],
   Indeterminate on the half-plane where there is no limit. (Not Undefined: Piecewise
   turns an Undefined default into a ConditionalExpression, and Normal of a SparseArray
   drops that.) A numeric entry with no limit is the Failure. *)
diagonalZeroBasePower[mat_] := With[
    {limits = zeroBaseLimit /@ Normal[Diagonal[mat]]},
    Replace[FirstCase[limits, Missing["NoLimit", t_] :> t], {
        _Missing :> SparseArray[Band[{1, 1}] -> limits, Dimensions[mat]],
        t_ :> noLimitFailure[t]
    }]
]

zeroBaseLimit[t_ ? NumericQ] := Which[exactZeroQ[t], 1, realPartSign[t] === 1, 0, True, Missing["NoLimit", t]]

zeroBaseLimit[t_ ? PossibleZeroQ] := 1

zeroBaseLimit[t_] := Piecewise[{{1, t == 0}, {0, Re[t] > 0}}, Indeterminate]

(* An inexact diagonal is held to the tolerance of a dense inexact matrix, so that 0^m
   commutes with a change of basis: an entry within it of zero counts as zero, and an
   entry whose real part does not exceed it has no limit. A machine m gives machine
   zeros and ones, and an arbitrary-precision m ones at its precision and exact zeros,
   as N does. *)
inexactDiagonalZeroBasePower[mat_] := With[
    {eigenvalues = Normal[Diagonal[mat]]},
    {scale = roundoffTolerance[mat] Max[Abs[eigenvalues]]},
    {zero = UnitStep[scale - Abs[eigenvalues]]},
    {noLimit = (1 - zero) UnitStep[scale - Re[eigenvalues]]},
    If[ Total[noLimit] > 0,
        noLimitFailure[First[Pick[eigenvalues, noLimit, 1]]],
        SparseArray[Band[{1, 1}] -> N[zero, resultPrecision[mat]], Dimensions[mat], N[0, resultPrecision[mat]]]
    ]
]

(* An exact Hermitian m: its zero eigenvalue is semisimple and the projector is the
   orthogonal one onto its null space, from the null vectors of m and their conjugates,
   which are the null vectors of its transpose. (Those are taken when an entry is not
   a Gaussian rational: Conjugate decides realness numerically, and on an entry that
   holds a zero on a branch cut hits the precision limit.) The exact nullity k and
   numerical eigenvalues decide the signs: each eigenvalue computed at precision p from
   the entries of m rounded to p digits lies within the roundoff bound
   100 n 10^-p ||m|| of an exact one (Weyl), so m has a negative eigenvalue when the
   smallest lies below minus the bound, and none when the (k+1)-th smallest lies above
   it. That exact eigenvalue is not zero, so doubling p until one of the two holds
   ends; p starts at machine precision, or at twice it when an entry of m is too small
   or too large for a machine number, and after 8 doublings the route for any exact
   matrix decides instead. *)
exactHermitianZeroBasePower[m_] := With[
    {x = NullSpace[m]},
    {spectrum = NestWhile[
        hermitianSpectrum[m, 2 #["Precision"]] &,
        hermitianSpectrum[m, If[machineRangeQ[m], MachinePrecision, 2 MachinePrecision]],
        undecidedSpectrumQ[#, Length[x]] &,
        1,
        8
    ]},
    Which[
        undecidedSpectrumQ[spectrum, Length[x]], exactZeroBasePower[m],
        First[spectrum["Eigenvalues"]] < - spectrum["Bound"], noLimitFailure[smallestEigenvalue[m, spectrum]],
        x === {}, SparseArray[{}, Dimensions[m]],
        True, nullSpaceProjector[Transpose[x], If[MatrixQ[m, gaussianRationalQ], Conjugate[x], NullSpace[Transpose[m]]]]
    ]
]

(* True when every entry of m that is not 0 lies between 10^-290 and 10^290 in size, so
   that m, its eigenvalues and their roundoff bound round to machine numbers without
   underflow. *)
machineRangeQ[m_] := AllTrue[
    DeleteCases[DeleteDuplicates[Flatten[m]], 0],
    With[{v = If[gaussianRationalQ[#], Abs[#], settledValue[#, Abs]]}, ! MissingQ[v] && 10 ^ -290 < v < 10 ^ 290] &
]

(* The eigenvalues of the exact Hermitian m at precision p in increasing order, and the
   bound within which each lies of an exact one. *)
hermitianSpectrum[m_, p_] := With[
    {mp = roundedEntries[m, p]},
    {eigenvalues = Sort[Re[Eigenvalues[mp]]]},
    <|"Eigenvalues" -> eigenvalues, "Bound" -> roundoffTolerance[mp] Max[Abs[eigenvalues]], "Precision" -> p|>
]

(* The entries of m, whose zeros are written as 0, rounded to p digits of its largest
   entry. An entry that is not a Gaussian rational is first valued to 3 more digits
   than p, and to 23 at least, of its own precision or of accuracy relative to the
   largest entry, whichever N reaches first, so that each entry is off by rounding
   alone. (The
   accuracy goal is reached on an entry that holds a zero not written as 0, such as
   4/3 + (GoldenRatio - (1 + Sqrt[5]) / 2), and the precision goal on a huge entry, as
   Exp[10^7], without all of its digits; neither goal is negative.) *)
roundedEntries[m_, p_] := With[
    {entries = DeleteDuplicates[Cases[m, Except[_ ? gaussianRationalQ], {2}]]},
    {digits = Max[N[p], 20] + 3},
    {accuracy = digits - Floor[Log10[Max[Abs[Cases[m, _ ? gaussianRationalQ, {2}]], Replace[settledValue[#, Abs] & /@ entries, _Missing -> 0, {1}], 10 ^ -3200]]]},
    N[Replace[m, Dispatch[(# -> valueWithin[#, Identity, {digits, Max[accuracy, digits]}]) & /@ entries], {2}], p]
]

undecidedSpectrumQ[spectrum_, k_] := With[
    {eigenvalues = spectrum["Eigenvalues"], bound = spectrum["Bound"]},
    First[eigenvalues] >= - bound && k < Length[eigenvalues] && eigenvalues[[k + 1]] <= bound
]

(* The smallest eigenvalue t of the exact Hermitian m, from its spectrum computed at 64
   digits at least. For m of Gaussian rationals it is exact when the rational or
   quadratic irrational r that Rationalize or RootApproximant finds within the bound of
   t, both applied to t scaled by a power of 10 to order one, is an eigenvalue, decided
   on f(m) for the minimal polynomial f of r, a matrix of Gaussian rationals. Otherwise
   it is t to machine precision, or to the fewer digits the bound leaves it, and t as
   computed when it lies within the bound of zero. *)
smallestEigenvalue[m_, spectrum_] := With[
    {refined = If[spectrum["Precision"] < 64, hermitianSpectrum[m, 64], spectrum]},
    {t = First[refined["Eigenvalues"]], bound = refined["Bound"]},
    If[ Abs[t] > bound,
        With[{digits = SetPrecision[t, Min[Log10[Abs[t] / bound], $MachinePrecision]]},
            If[MatrixQ[m, gaussianRationalQ], exactEigenvalueNear[m, t, bound, refined["Eigenvalues"], digits], digits]
        ],
        t
    ]
]

(* The rational that Rationalize finds within the bound of t, or else the quadratic
   irrational that RootApproximant finds from the digits of t the bound leaves, at most
   50, both applied to t scaled by a power of 10 to order one, when it is an eigenvalue
   of m; default otherwise. (Given the digits past the bound, RootApproximant fits a
   quadratic to the roundoff in them and evaluates its roots past the precision
   limit.) *)
exactEigenvalueNear[m_, t_, bound_, eigenvalues_, default_] := With[
    {scale = 10 ^ -Round[Log10[Abs[t]]], acceptQ = Abs[valueAt[#, Identity, 10 - Log10[bound]] - t] <= bound && eigenvalueNearQ[m, #, eigenvalues, bound] &},
    With[{r = Rationalize[scale t, scale bound] / scale},
        If[acceptQ[r], r, With[{q = RootApproximant[SetPrecision[scale t, Min[Log10[Abs[t] / bound], 50]], 2] / scale}, If[acceptQ[q], q, default]]]
    ]
]

(* True when the rational or quadratic irrational r is an eigenvalue of the Hermitian m
   of Gaussian rationals, given its eigenvalues each computed within the bound: the
   nullity of f(m), for the minimal polynomial f of r, counts the eigenvalues equal to
   a root of f, and it equals the number computed within the bound of one only when
   each of those is a root. *)
eigenvalueNearQ[m_, r_, eigenvalues_, bound_] := With[
    {f = MinimalPolynomial[r, \[FormalX]]},
    {c = CoefficientList[f, \[FormalX]], roots = valueAt[#, Identity, 10 - Log10[bound]] & /@ SolveValues[f == 0, \[FormalX]]},
    Length[m] - MatrixRank[Total[MapThread[#1 #2 &, {c, NestList[m . # &, IdentityMatrix[Length[m], SparseArray], Length[c] - 1]}]]] ==
        Count[eigenvalues, e_ /; AnyTrue[roots, Abs[e - #] <= bound &]]
]

(* An exact m. Its core-nilpotent decomposition m = t.(c (+) n).t^-1 separates the
   nonsingular core c from the nilpotent part n on the generalized null space, so the
   zero eigenvalue is semisimple exactly when n = 0, that is when the rank of m equals
   the size of c (the entries of n, zero, can be written as expressions too long to
   evaluate), the other eigenvalues are those of c, and the projector is
   t.(0 (+) 1).t^-1, the last columns of t against the last rows of its inverse. *)
exactZeroBasePower[m_] := coreNilpotentZeroBasePower[m, CoreNilpotentDecomposition[m]]

(* True when the zero eigenvalue of the exact m, of algebraic multiplicity k, is
   semisimple, that is when m has k independent null vectors. At k <= 1 it is: a
   nilpotent part of size 1 is 0. Otherwise the null vectors NullSpace finds are each
   checked exactly, m.x = 0 entry by entry, and when that fails the rank of m with
   every pivot decided as exactZeroQ decides settles it. (That rank alone is slow: its
   pivots grow unsimplified.) *)
semisimpleZeroQ[_, k_ /; k <= 1] := True

semisimpleZeroQ[m_ ? (MatrixQ[#, gaussianRationalQ] &), k_] := MatrixRank[m] == Length[m] - k

semisimpleZeroQ[m_, k_] := With[{x = NullSpace[m]},
    If[ Length[x] == k && AllTrue[Flatten[m . Transpose[x]], exactZeroQ],
        True,
        MatrixRank[m, ZeroTest -> exactZeroQ] == Length[m] - k
    ]
]

coreNilpotentZeroBasePower[m_, {t_, core_, nilpotent_}] := If[
    semisimpleZeroQ[m, Length[m] - Length[core]],
    Replace[If[core === {}, Missing[], coreNoLimit[corePolynomial[m, Length[m] - Length[core]]]], _Missing :> Which[
        nilpotent === {}, SparseArray[{}, Dimensions[t]],
        True, t[[All, Length[core] + 1 ;;]] . Inverse[t][[Length[core] + 1 ;;]]
    ]],
    defectiveZeroFailure
]

(* The characteristic polynomial of the core, that of m divided by (-z)^k for the size
   k of the nilpotent part, which is zero here: det(m - z) = det(c - z) (-z)^k, so its
   k lowest coefficients vanish, however they are written, and are dropped. From m its
   coefficients stay small, where the entries of the core the decomposition returns
   can grow by orders of magnitude in size. *)
corePolynomial[m_, k_] := keepingPowers[
    With[{c = CoefficientList[CharacteristicPolynomial[#, \[FormalZ]], \[FormalZ]]}, Together[(-1) ^ k Drop[c, k]] . \[FormalZ] ^ Range[0, Length[c] - k - 1]] &,
    m
]

(* Missing[] when every root of the characteristic polynomial p of the exact core has
   positive real part, and otherwise the Failure naming one that has not. *)
coreNoLimit[p_] := If[rightHalfPlaneQ[p], Missing[], noLimitFailure[namedNoLimitEigenvalue[p]]]

(* True when every root of the polynomial p has positive real part, decided exactly and
   without the roots. A real p is kept; otherwise p times the polynomial of conjugate
   coefficients is real, and its roots are those of p and their conjugates. Either way
   the roots all have positive real part exactly when the roots of that polynomial at
   -s all have negative real part, which the Routh-Hurwitz test decides from the
   coefficients. *)
rightHalfPlaneQ[p_] := With[
    {q = realCoefficients[CoefficientList[p, \[FormalZ]]]},
    hurwitzStableQ[Reverse[togetherKeepingPowers[q] (-1) ^ Range[0, Length[q] - 1]]]
]

(* The coefficients a of a polynomial when they are real, as for the characteristic
   polynomial of every Lindblad generator, and otherwise those of its product with the
   polynomial of conjugate coefficients. *)
realCoefficients[a_] := With[{conjugate = conjugateEntry /@ a}, If[conjugate === a, a, ListConvolve[a, conjugate, {1, -1}, 0]]]

(* The complex conjugate of the exact matrix c. An entry whose imaginary part is zero,
   decided by value, is kept. In any other entry I goes to -I, a Root object to its
   conjugate, and to their conjugates as a whole a power with an exponent that is not
   an integer and a base whose real part is not positive, and any other function value
   whose imaginary part is not zero, as ArcCos[2] or Log[-2]; a function value whose
   imaginary part is zero, as Im[ArcCos[2]], is kept whole, since conjugating inside it
   would change it. So c is real when this leaves it unchanged. (Conjugate itself
   decides whether a long real expression is real numerically, and hits the precision
   limit on one with cancellations.) *)
exactConjugate[c_] := Map[conjugateEntry, c, {2}]

conjugateEntry[x_ ? gaussianRationalQ] := Conjugate[x]

conjugateEntry[x_] := If[realPartSign[I x] === 0, x, x /. {
    Complex[a_, b_] :> Complex[a, -b],
    y : Power[b_, e_] /; ! IntegerQ[e] && realPartSign[b] =!= 1 :> Conjugate[y],
    r_Root :> Conjugate[r],
    y : f_Symbol[__] /; NumericQ[y] && ! MemberQ[{Plus, Times, Power}, f] :> If[realPartSign[I y] === 0, y, Conjugate[y]]
}]

(* True when every root of the real polynomial with coefficients a, from the highest
   power down, has negative real part: every leading entry of the Routh array is
   positive. Each row of the array is the row two above it less the multiple of the
   row above that clears its leading entry, and the rows stop at the first leading
   entry that is not positive, so the roots are all in the left half-plane when a[[1]]
   is positive and the rows reach the last one with a positive leading entry. *)
hurwitzStableQ[a_List] := With[
    {width = Ceiling[Length[a] / 2] + 1},
    {pairs = NestWhileList[routhStep, {PadRight[a[[1 ;; ;; 2]], width], PadRight[a[[2 ;; ;; 2]], width]}, lowerLeadingPositiveQ, 1, Length[a] - 2]},
    realPartSign[First[a]] === 1 && Length[pairs] == Length[a] - 1 && lowerLeadingPositiveQ[Last[pairs]]
]

routhStep[{upper_, lower_}] := {lower, togetherKeepingPowers[Append[Rest[upper] - First[upper] / First[lower] Rest[lower], 0]]}

(* f of x with every integer power beyond the 64th of a compound number held as a
   symbol, and put back after, for f that holds whatever the values of those powers,
   as Together, the characteristic polynomial, a factorization and the roots of a
   linear or quadratic factor do: Together, CharacteristicPolynomial, Exponent and
   FactorList expand (Sqrt[2] - 1)^1000000 into a polynomial in Sqrt[2], for
   minutes. A power of an atom, as E^1000, cannot be expanded and is left as it
   is. *)
keepingPowers[f_, x_] := With[
    {powers = DeleteDuplicates[Cases[x, _ ? heldPowerQ, {0, Infinity}]]},
    {held = Array[\[FormalCapitalP], Length[powers]]},
    f[x /. Thread[powers -> held]] /. Thread[held -> powers]
]

heldPowerQ[x_] := MatchQ[x, Power[b_, n_Integer] /; Abs[n] > 64 && ! AtomQ[b] && NumericQ[b]]

togetherKeepingPowers[x_] := keepingPowers[Together, x]

(* The degree in z of the polynomial p, whose leading coefficient is not zero, with
   the powers keepingPowers holds held: Exponent expands them. *)
polynomialDegree[p_] := keepingPowers[Exponent[#, \[FormalZ]] &, p]

lowerLeadingPositiveQ[{_, lower_}] := realPartSign[First[lower]] === 1

(* An eigenvalue without positive real part among the roots of the characteristic
   polynomial p, named exactly when it is a root of a linear or quadratic factor, or
   of a factor with real coefficients, decided by value, that holds no power that
   keepingPowers holds, and Missing[] when such eigenvalues are roots only of other
   factors. (The kernel crashes building the exact roots of a factor with complex
   coefficients of degree 32, as Eigenvalues of a dense complex 32 x 32 matrix does,
   and a Root object expands such a power.) *)
namedNoLimitEigenvalue[p_] := SelectFirst[
    Catenate[keepingPowers[SolveValues[# == 0, \[FormalZ]] &, #] & /@ Select[
        polynomialFactors[p],
        polynomialDegree[#] >= 1 && (
            polynomialDegree[#] <= 2 ||
                FreeQ[#, _ ? heldPowerQ] && AllTrue[CoefficientList[#, \[FormalZ]], realPartSign[I #] === 0 &]
        ) &
    ]],
    realPartSign[#] =!= 1 &
]

(* The factors of the polynomial p. A p of degree 3 or more is factored, with the
   powers keepingPowers holds held, when Element proves its coefficients algebraic
   and they hold at most 3 other distinct radicals or Root objects; FactorList fails on
   others, as Exp[10^7], or takes minutes, as with 8 square roots, and p is then its
   own factor. (Element is not asked about an integer exponent beyond 10^9, on some
   of which it crashes the kernel.) *)
polynomialFactors[p_] := If[
    polynomialDegree[p] >= 3 && FreeQ[p, Power[_, n_Integer /; Abs[n] > 10 ^ 9]] &&
        AllTrue[CoefficientList[p, \[FormalZ]], TrueQ[Element[#, Algebraics]] &],
    keepingPowers[
        If[ Length[DeleteDuplicates[Cases[CoefficientList[#, \[FormalZ]], Power[_Integer | _Rational, _Rational] | _Root | _AlgebraicNumber, {0, Infinity}]]] <= 3,
            FactorList[#, Extension -> Automatic][[All, 1]],
            {#}
        ] &,
        p
    ],
    {p}
]

(* Exact numbers are decided by value. The value of x at 50, then 400 and 3200 digits of
   precision or of accuracy, whichever N reaches first, tells a nonzero x from zero, and
   a sign read off it is exact. (The accuracy goal is reached on a number that holds a
   zero not written as 0, and the precision goal on a huge number without all of its
   digits.) What stays within 10^-3200 of zero is zero when Element proves it
   algebraic, after FunctionExpand, and algebraicZeroQ, which is exact, finds it zero,
   as for GoldenRatio - (1 + Sqrt[5]) / 2; it is taken as zero when Element does not, as
   for Log[6] - Log[2] - Log[3], and when algebraicZeroQ does not settle it in its time,
   so a number that is not zero, lies within 10^-3200 of it and is not proved
   algebraic, as Tanh[10^4] - 1, is read as zero. A real part proved algebraic and
   found nonzero is read to digits doubling from 6400 to 102400, and past that its sign
   is read from RootReduce within 5 s; it is Indeterminate, which counts as not
   positive, when RootReduce does not end in time. With these readings the decisions
   are exact for numbers between $MinNumber and $MaxNumber, the range of
   arbitrary-precision numbers. Beyond it N reports an underflow or an overflow with
   its own messages: a value that underflows is read as zero, proved algebraic or not
   (PossibleZeroQ with the method ExactAlgebraics crashes the kernel on some, as
   (Sqrt[2] - 1)^(10^20) - (Sqrt[2] - 1)^(10^20 + 1)), and one that overflows as
   nonzero, with the sign of its real part as Sign reads it from the expression, as for
   -Exp[Exp[100]], and Indeterminate when Sign cannot. *)
exactZeroQ[x_ ? gaussianRationalQ] := x == 0

exactZeroQ[x_] := Replace[settledValue[x, Abs], {
    Missing["Overflow"] -> False,
    Missing["Underflow"] -> True,
    Missing["NotANumber"] :> With[{y = simplifiedForm[x]}, If[y =!= x, exactZeroQ[y], False]],
    _Missing :> Replace[provenAlgebraic[x], {_Missing -> True, y_ :> Replace[algebraicZeroQ[y], _Missing -> True]}],
    _ -> False
}]

(* x as FullSimplify writes it, within 5 s, for a number N cannot value as written,
   because a part of it overflows or underflows, as
   Exp[Exp[100]] / (1 + Exp[Exp[100]]); x itself when that fails. *)
simplifiedForm[x_] := TimeConstrained[FullSimplify[x], 5, x]

(* The sign of the real part of x when Sign can read it from the expression, as for
   -Exp[Exp[100]], and Indeterminate when it cannot. *)
expressionSign[x_] := Replace[Sign[Re[x]], Except[-1 | 0 | 1] -> Indeterminate]

(* True when the algebraic number y is zero and False when it is not, decided exactly,
   and Missing[] when neither decision ends in its time. PossibleZeroQ with the method
   ExactAlgebraics is fast on sums of many radicals and slow on high powers, as
   (1 + Sqrt[2])^5225 - (Sqrt[2] - 1)^-5225 (minutes); RootReduce the other way round.
   So the first is given 1 s and the second 5 s after it; both take minutes on the sum
   of that difference and S^2 - Expand[S^2] for S the sum of the square roots of the
   primes up to 19. *)
algebraicZeroQ[y_] := TimeConstrained[
    PossibleZeroQ[y, Method -> "ExactAlgebraics"],
    1,
    TimeConstrained[RootReduce[y] === 0, 5, Missing[]]
]

(* True when the number base of a power base^op is zero, without the conventions of
   exactZeroQ, which read a tiny nonzero number such as Exp[-10^4] as zero: a
   Gaussian rational that is 0, and a number that stays within 10^-3200 of zero and
   is zero exactly, by algebraicZeroQ when Element proves it algebraic and otherwise
   when FullSimplify reduces it to 0, as Log[6] - Log[2] - Log[3]. One it does not
   reduce, as Erf[100] - 1, or that algebraicZeroQ does not settle in its time, is taken
   as nonzero, and base^op is computed from it. (PossibleZeroQ assumes such a number
   zero, with a message.) *)
zeroBaseQ[x_ ? gaussianRationalQ] := x == 0

zeroBaseQ[x_] := Replace[settledValue[x, Abs], {
    Missing["Overflow"] -> False,
    Missing["NotANumber"] :> With[{y = simplifiedForm[x]}, If[y =!= x, zeroBaseQ[y], False]],
    _Missing :> Replace[provenAlgebraic[x], {_Missing :> TimeConstrained[FullSimplify[x] === 0, 5, False], y_ :> Replace[algebraicZeroQ[y], _Missing -> False]}],
    _ -> False
}]

(* The sign of the real part of the exact number x. *)
realPartSign[x_ ? gaussianRationalQ] := Sign[Re[x]]

realPartSign[x_] := Replace[settledValue[x, Re], {
    Missing["Overflow"] :> expressionSign[x],
    Missing["Underflow"] -> 0,
    Missing["NotANumber"] :> With[{y = simplifiedForm[x]}, If[y =!= x, realPartSign[y], expressionSign[x]]],
    _Missing :> With[{re = Replace[provenAlgebraic[x], {_Missing :> Re[x], y_ :> Re[y]}]},
        If[exactZeroQ[re], 0, tinySign[x, re]]
    ],
    v_ :> Sign[v]
}]

(* The sign of re, the real part of x, which is not zero and lies within 10^-3200 of it:
   read at digits doubling from 6400, 4 times, and past that from RootReduce, which is
   exact, within 5 s; Indeterminate when neither tells. *)
tinySign[x_, re_] := Replace[digitsValue[x, Re, NestWhile[2 # &, 6400, TrueQ[digitsValue[x, Re, #] == 0] &, 1, 4]], {
    v_ /; NumberQ[v] && TrueQ[v != 0] :> Sign[v],
    _ :> Replace[TimeConstrained[reducedSign[RootReduce[re]], 5, Indeterminate], Except[-1 | 0 | 1] -> Indeterminate]
}]

(* The sign of the algebraic number r as RootReduce writes it, read by Sign with 4 times
   the digits of its integers as extra working precision: a quadratic irrational
   a + b Sqrt[2] that is not zero, with a and b of d digits, lies no closer than about
   10^-d to zero, as 1 / (1 + (1 + Sqrt[2])^(10^6)) does, of 382776 digits, where
   Sign under the default limit gives up. *)
reducedSign[r_] := Block[{$MaxExtraPrecision = 10000 + 4 integerDigits[r]}, Sign[r]]

(* x, or x as FunctionExpand writes it, when Element proves it algebraic, as
   Sin[ArcCos[1/3] / 2] = 1 / Sqrt[3]; Missing[] when Element does not. Element is
   asked first, since FunctionExpand takes minutes on 1 / (1 + (1 + Sqrt[2])^(10^6)),
   which Element proves at once; and neither is asked about a number with an integer
   exponent beyond 10^9, on some of which Element crashes the kernel, as
   (Sqrt[2] - 1)^(10^20) / ((Sqrt[2] - 1)^(10^20) + (Sqrt[2] - 1)^(10^20 + 1)). *)
provenAlgebraic[x_] /; ! FreeQ[x, Power[_, n_Integer /; Abs[n] > 10 ^ 9]] := Missing[]

provenAlgebraic[x_] := If[TrueQ[Element[x, Algebraics]],
    x,
    With[{y = TimeConstrained[FunctionExpand[x], 1, x]}, If[y =!= x && TrueQ[Element[y, Algebraics]], y, Missing[]]]
]

(* The part of the value of x, as read by part, at the first of 50, 400 and 3200 digits
   that tells it from zero; Missing[] when none does, Missing["Underflow"] or
   Missing["Overflow"] when the value underflows or overflows, and
   Missing["NotANumber"] when N gives no number, as when a part of x overflows while
   x does not. *)
settledValue[x_, part_] := firstNonzeroValue[x, part, {50, 400, 3200}]

firstNonzeroValue[_, _, {}] := Missing[]

firstNonzeroValue[x_, part_, {digits_, rest___}] := With[{v = digitsValue[x, part, digits]},
    Which[
        v === Overflow[], Missing["Overflow"],
        v === Underflow[], Missing["Underflow"],
        ! NumberQ[v], Missing["NotANumber"],
        TrueQ[v != 0], v,
        True, firstNonzeroValue[x, part, {rest}]
    ]
]

(* The part of x, as part picks it, valued to the given precision or accuracy, whichever
   N reaches first (Infinity for a goal not set). The part is taken before N, so that
   the goal applies to it: the real part of 5 + I Exp[10^7] is 5, which the digits of
   the whole number would not show, and N reads Abs[I Exp[10^7]] where it hits the
   precision limit on I Exp[10^7]. *)
valueWithin[x_, part_, {precision_, accuracy_}] :=
    Block[{$MaxExtraPrecision = extraPrecision[x, Max[Select[{precision, accuracy}, NumericQ]]]}, N[part[x], {precision, accuracy}]]

digitsValue[x_, part_, digits_] := valueWithin[x, part, {digits, digits}]

(* The part of x valued to the given accuracy. *)
valueAt[x_, part_, accuracy_] := valueWithin[x, part, {Infinity, accuracy}]

(* The working precision N may add to reach its goal: twice the digits asked for, plus
   200, plus twice the size in digits of the largest integer or part of x, as its
   structure shows it, up to 10^5 digits, enough for the cancellation between such
   parts (as (1 + Sqrt[2])^52250 - (Sqrt[2] - 1)^-52250 or
   E^(10^5 Sqrt[2]) - Cosh[10^5 Sqrt[2]] - Sinh[10^5 Sqrt[2]]). A part of more than
   10^6 digits does not count: no working precision given here reaches a cancellation
   between such parts, and N would spend it in vain, for minutes on
   (Sqrt[2] - 1)^(10^20) / ((Sqrt[2] - 1)^(10^20) + (Sqrt[2] - 1)^(10^20 + 1)), whose
   parts underflow. The limit is finite: under an unlimited one
   N never finishes on a number that holds an exact zero where a function jumps, as
   Sqrt[-1 + (GoldenRatio - (1 + Sqrt[5]) / 2) I] - I does on the branch cut of Sqrt.
   And it is no larger than that: on a number that is zero, N raises its precision to
   the limit before it gives the zero, so the limit sets the time, minutes for
   Gamma[1/3] Gamma[2/3] - 2 Pi / Sqrt[3] under a limit of tens of thousands of digits.
   A number whose cancellation needs more, as between terms far beyond 10^100000 in size
   or between parts whose size the structure does not show, as BesselI[0, 10^5], is not
   told from zero: it is read as zero, as one within 10^-3200 of it is, and N says so
   with its own message. *)
extraPrecision[x_, digits_] := 2 Max[digits, 0] + 200 + 2 Min[10 ^ 5, Max[0, Select[
    Join[integerLengths[x], Abs[sizeDigits /@ Level[x, {0, Infinity}]]],
    # <= 10 ^ 6 &
]]]

(* The number of digits of each integer in x, of the numerator and denominator of each
   rational, and of the largest integer of all. *)
integerLengths[x_] := Cases[x, r : _Integer | _Rational | _Complex :> IntegerLength[Max[Abs[Numerator[{Re[r], Im[r]}]], Denominator[{Re[r], Im[r]}]]], {0, Infinity}]

integerDigits[x_] := Max[0, integerLengths[x]]

(* About log10 of the size of the number x, read from its structure alone, so that no
   part of x is valued (a value can underflow, or hit the precision limit on a zero
   hidden on a branch cut): log10 itself for integers, rationals and constants such as
   Pi; for a power, the exponent times the size of the base, and for Cosh or Sinh of y,
   |y| log10(e), where the exponent or y is a rational of moderate size, and otherwise
   the same with 10^s, s the size of the exponent or of y up to 12, in place of its
   magnitude (its sign unread, so a power with a negative exponent counts as large);
   the sum for a product, the largest term plus log10 of their number for a sum, and 0
   for anything else. *)
sizeDigits[x : _Integer | _Rational] := If[x == 0, 0, N[Log10[Abs[x]]]]

sizeDigits[Complex[a_, b_]] := Max[sizeDigits[a], sizeDigits[b]]

sizeDigits[x_Symbol ? NumericQ] := Log10[Abs[N[x]]]

sizeDigits[Power[b_, e : _Integer | _Rational]] /; 10 ^ -6 < Abs[e] < 10 ^ 12 := e sizeDigits[b]

sizeDigits[Power[b_, e_]] /; unsizedNumberQ[e] := 10 ^ Min[sizeDigits[e], 12] Abs[sizeDigits[b]]

sizeDigits[(Cosh | Sinh)[y : _Integer | _Rational]] /; Abs[y] < 10 ^ 12 := Abs[y] Log10[N[E]]

sizeDigits[(Cosh | Sinh)[y_]] /; unsizedNumberQ[y] := 10 ^ Min[sizeDigits[y], 12] Log10[N[E]]

(* True for a number whose magnitude the rules above do not take as it is: one that is
   not a rational, or a rational of 10^12 or more. *)
unsizedNumberQ[y_] := NumericQ[y] && ! (MatchQ[y, _Integer | _Rational] && Abs[y] < 10 ^ 12)

sizeDigits[x_Times] := Total[sizeDigits /@ List @@ x]

sizeDigits[x_Plus] := Max[sizeDigits /@ List @@ x] + Log10[N[Length[x]]]

sizeDigits[_] := 0

gaussianRationalQ[x_] := MatchQ[x, _Integer | _Rational | Complex[_Integer | _Rational, _Integer | _Rational]]

(* An inexact m. The limit exists when its zero eigenvalue is semisimple to within
   roundoff (roundoffZeroEigenvalue) and every other eigenvalue has real part above
   scale = roundoffTolerance[m] ||m||, and it is then the projector from the null spaces.
   A Hermitian m gives an exactly Hermitian projector. *)
inexactZeroBasePower[m_] := With[
    {zero = roundoffZeroEigenvalue[m, Eigenvalues[m]]},
    {noLimit = SelectFirst[Drop[zero["Eigenvalues"], zero["Nullity"]], Re[#] <= zero["Tolerance"] zero["Norm"] &]},
    Which[
        zero["Defective"], defectiveZeroFailure,
        ! MissingQ[noLimit], noLimitFailure[noLimit],
        zero["Nullity"] == 0, SparseArray[{}, Dimensions[m], N[0, resultPrecision[m]]],
        True, With[{p = nullSpaceProjector @@ zero["NullSpaces"]}, If[nearlyHermitianQ[m, roundoffTolerance[m]], (p + ConjugateTranspose[p]) / 2, p]]
    ]
]

(* The eigenvalue 0 of an inexact m to within roundoff, read from the singular value
   decomposition of m and its eigenvalues, which are kept sorted by modulus. Roundoff
   moves each singular value by at most of order eps ||m||, so those within scale =
   roundoffTolerance[m] ||m|| count as zero and give the nullity k, and the singular
   vectors give x and y. The zero eigenvalue is semisimple to within roundoff when the
   smallest singular value of y.x, the cosine of the largest angle between the null
   spaces of m and m^†, stays above the tolerance; its computed eigenvalues then lie
   within scale over that cosine of zero, so they are the k eigenvalues of smallest
   modulus, and every other eigenvalue must lie outside that bound, or the zero
   eigenvalue is defective to within roundoff. *)
roundoffZeroEigenvalue[m_, eigenvalues_] := roundoffZeroEigenvalue[m, SingularValueDecomposition[m], SortBy[eigenvalues, Abs]]

roundoffZeroEigenvalue[m_, {u_, s_, v_}, eigenvalues_] := With[
    {sigma = Diagonal[s], tol = roundoffTolerance[m]},
    {scale = tol Max[sigma]},
    {k = Total[UnitStep[scale - sigma]]},
    {x = Take[v, All, -k], y = ConjugateTranspose[Take[u, All, -k]]},
    {cosine = If[k == 0, 1, Min[SingularValueList[y . x, Tolerance -> 0]]]},
    <|
        "Nullity" -> k, "NullSpaces" -> {x, y}, "Eigenvalues" -> eigenvalues, "Tolerance" -> tol, "Norm" -> Max[sigma],
        "Defective" -> (cosine <= tol || Max[Abs[Take[eigenvalues, k]], 0] > scale / cosine || Min[Abs[Drop[eigenvalues, k]], Infinity] <= scale / cosine)
    |>
]

(* The spectral projector onto the null space along the range, from columns x
   spanning the null space and rows y spanning it from the left. *)
nullSpaceProjector[x_, y_] := x . Inverse[y . x] . y

(* A symbolic m that is not diagonal: the projector from the null spaces of m and of
   its transpose for generic values of the symbols, rational in them, where every
   other eigenvalue has positive real part. Each nonzero entry of the projector, or
   every entry when it is zero, carries that condition as
   Piecewise[{{p, condition}}, Indeterminate]. The other eigenvalues enter only
   through the condition, so values where two of them collide, as at an exceptional
   point of a Lindblad generator, leave the projector finite. Where a generically
   nonzero eigenvalue vanishes, the null space grows and the closed form is
   Indeterminate; with the symbols declared as parameters the values reach the limit
   there directly. A zero eigenvalue with fewer null vectors than its multiplicity at
   generic values, or a numeric eigenvalue with no limit, is the Failure. *)
symbolicZeroBasePower[m_] := With[
    {x = NullSpace[m], others = DeleteCases[Eigenvalues[m], t_ /; If[NumericQ[t], exactZeroQ[t], PossibleZeroQ[t]]]},
    {
        noLimit = SelectFirst[others, NumericQ[#] && realPartSign[#] =!= 1 &],
        condition = And @@ DeleteDuplicates[Re[#] > 0 & /@ Select[others, ! NumericQ[#] &]]
    },
    Which[
        Length[x] + Length[others] < Length[m], defectiveZeroFailure,
        ! MissingQ[noLimit], noLimitFailure[noLimit],
        x === {}, ConstantArray[Piecewise[{{0, condition}}, Indeterminate], Dimensions[m]],
        True, Map[
            If[PossibleZeroQ[#], 0, Piecewise[{{#, condition}}, Indeterminate]] &,
            Together[nullSpaceProjector[Transpose[x], NullSpace[Transpose[m]]]],
            {2}
        ]
    ]
]

noLimitFailure[_Missing] := Failure["ZeroBasePowerNoLimit", <|
    "MessageTemplate" -> "0^m is the limit of b^m as b -> 0, which does not exist: m has an eigenvalue neither zero nor with positive real part."
|>]

noLimitFailure[t_] := Failure["ZeroBasePowerNoLimit", <|
    "MessageTemplate" -> "0^m is the limit of b^m as b -> 0, which does not exist: m has the eigenvalue `1`, neither zero nor with positive real part.",
    "MessageParameters" -> {t}
|>]

defectiveZeroFailure = Failure["ZeroBasePowerDefective", <|
    "MessageTemplate" -> "0^m is the limit of b^m as b -> 0, which does not exist: the zero eigenvalue of m is defective, where b^m grows as Log[b]."
|>]


(* For the exponential, a matrix counts as diagonal when it is square and every
   off-diagonal entry is zero, with no tolerance: a coupling at roundoff relative to
   the largest entry still mixes two degenerate levels over a long time, and
   MatrixExp resolves that mixing. A symbolic off-diagonal entry counts as zero when
   DiagonalMatrixQ's zero test proves it zero; that test is best effort, and an
   identically zero entry it cannot prove sends the matrix to MatrixExp. *)
diagonalMatrixQ[mat_] := SquareMatrixQ[mat] && DiagonalMatrixQ[mat, Tolerance -> 0]

(* The exponential of a diagonal matrix is the diagonal of the exponentials of its
   entries, and its action on a vector multiplies the vector by them. Any other
   matrix, and a diagonal with an entry that has no value (infinite, indeterminate),
   goes to MatrixExp, which fails on the latter. *)
matrixExponential[mat_ ? diagonalMatrixQ, v___] := With[{d = Normal[Diagonal[mat]]},
    diagonalAction[diagonalExp[d], v] /; ! valuelessEntriesQ[d]
]
matrixExponential[mat_, v___] := MatrixExp[mat, v]

diagonalAction[values_] := DiagonalMatrix[values, TargetStructure -> "Sparse"]
diagonalAction[values_, v_] := values v

(* A numeric diagonal at machine precision, exact entries included, is exponentiated
   in machine numbers, as MatrixExp does, and as one packed array. On a packed array
   Exp returns an exponential below the smallest normalized machine number as a
   subnormal number or zero, and one above the largest as an arbitrary-precision
   number, without a message; on an unpacked list it raises General::munfl. A numeric
   diagonal at a higher finite precision is exponentiated at that precision. *)
diagonalExp[d_ ? machineVectorQ] := Exp[packedVector[N[d]]]
diagonalExp[d_ ? (VectorQ[#, NumericQ] && Precision[#] < Infinity &)] := Exp[N[d, Precision[d]]]
diagonalExp[d_] := Exp[d]

machineVectorQ[d_] := VectorQ[d, NumericQ] && Precision[d] === MachinePrecision

packedVector[y_ ? (FreeQ[#, _Complex] &)] := Developer`ToPackedArray[y, Real]
packedVector[y_] := Developer`ToPackedArray[y, Complex]


SetPrecisionNumeric[x_ /; NumericQ[x] || ArrayQ[x, _, NumericQ]] := SetPrecision[x, $MachinePrecision - 3]

SetPrecisionNumeric[x_] := x


(* helpers *)

toggleSwap[xs : {_Integer...}, n_Integer] := MapIndexed[(#1 > n) != (First[#2] > n) &, xs]

toggleShift[xs : {_Integer...}, n_Integer] := n - Subtract @@ Total /@ TakeDrop[Boole @ toggleSwap[xs, n], n]


alignDimensions[xs_, {}] := {{xs}, {xs}}

alignDimensions[{}, ys_] := {{ys}, {ys}}

alignDimensions[xs : {_Integer..}, ys : {_Integer..}] := Module[{
    as = FoldList[Times, xs], bs = FoldList[Times, ys], p, first, second
},
    p = Min[Intersection[as, bs]];
    If[ IntegerQ[p],
        first = TakeDrop[xs, First @ FirstPosition[as, p]];
        second = TakeDrop[ys, First @ FirstPosition[bs, p]];
        DeleteCases[{}] /@ MapThread[
            ReverseApplied[Prepend],
            {
                {first, second}[[All, 1]],
                alignDimensions @@ {first, second}[[All, 2]]
            }
        ],

        Missing[]
    ]
]


TranscendentalRecognize[num_ ? NumericQ, basis : _ ? VectorQ : {Pi}] := Enclose[
    Block[{lr, ans},
        (* TODO: identify exact FindIntegerNullVector::* messages emitted on no-relation cases; until then, ConfirmBy -> Enclose handles failure path *)
        lr = ConfirmBy[FindIntegerNullVector[Prepend[N[basis, Precision[num]], num]], ListQ];
        ans = Rest[lr] . basis / First[lr];
        If[Numerator[ans] > 1*^3 || Denominator[ans] > 1*^3,
            num,
            Sign[N[ans]] Sign[num] ans
        ]
    ],
    num &
]


$QuantumFrameworkProfile = False

profile[label_] := Function[{expr}, If[TrueQ[$QuantumFrameworkProfile], EchoTiming[expr, label], expr], HoldFirst]

ReverseHalf[list_] /; Divisible[Length[list], 2] := Catenate @ Reverse[TakeDrop[list, Length[list] / 2]]

HadamardGrayRowPermutation[n_Integer ? Positive] := FindPermutation @ Nest[Insert[#, Splice @ ReverseHalf[# + Length[#]], Length[#] / 2 + 1] &, {1, 2}, n - 1]

(* Safe store for the *Prop property caches. Replaces the fragile
   `Quiet[HeadProp[obj, prop, args] = result, Rule::rhs]` (and bare Set) used in the object
   wrappers: it commits the cache DownValue only when the held key is free of pattern
   constructs. A pattern-bearing key - only ever a symbolic / meta call - would otherwise create
   an over-matching DownValue (the latent bug the old code masked with Quiet) and emit Rule::rhs;
   such keys are simply left uncached. The key is small (obj + prop + args), so the FreeQ check
   is cheap, and no Quiet is needed. Returns result either way. *)
SetAttributes[cacheProperty, HoldFirst];
cacheProperty[prop_, result_] := (
    If[ FreeQ[Unevaluated[prop], _Pattern | _Blank | _BlankSequence | _BlankNullSequence | _Optional],
        prop = result
    ];
    result
)

MyBlockDiagonalMatrix[{m1_ ? matrixQ, m2_ ? matrixQ}] := Block[{
    r1, r2, c1, c2
},
  	{r1, c1} = tensorDimensions @ m1;
  	{r2, c2} = tensorDimensions @ m2;
    Join[
        Join[m1, ConstantArray[0, {r1, c2}, SparseArray], 2],
        Join[ConstantArray[0, {r2, c1}, SparseArray], m2, 2]
    ]
]

MyBlockDiagonalMatrix[mm : {__ ? matrixQ}] := Fold[MyBlockDiagonalMatrix[{##}] &, mm]

MyBlockDiagonalMatrix[{}] := {{}}


ActivateTensor[expr_] := Activate[
    Activate[expr /. GeneralizedPower[TensorProduct, t_, n_Integer] :> Inactive[TensorProduct] @@ ConstantArray[Activate[t, TensorContract], n], TensorContract]
]


QuditAdjacencyMatrix[x_, d_ : Automatic] := {{Abs[Normalize[x]] Mod[Arg[x], 2 Pi] / (Pi / (2 Replace[d, Automatic -> 2]))}}
QuditAdjacencyMatrix[xs_ ? VectorQ, dim_ : Automatic] := With[{d = Replace[dim, Automatic :> Length[xs]]},
	Abs[xs] DiagonalMatrix[Mod[Arg[xs], 2 Pi] / (Pi / (2 d)) + 1] + If[d > 1, SparseArray[Band[{2, 1}] -> 1, {d, d}], 0]
]

