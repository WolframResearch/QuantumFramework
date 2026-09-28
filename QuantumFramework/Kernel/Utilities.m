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
PackageScope["SetPrecisionNumeric"]
PackageScope["TranscendentalRecognize"]

PackageScope["QuditAdjacencyMatrix"]

PackageScope["$QuantumFrameworkProfile"]
PackageScope["profile"]
PackageScope["Memoize"]
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



Options[eigensystem] = {"Sort" -> False, "Normalize" -> False, "Orthogonalize" -> False, Chop -> False}

eigensystem[matrix_, OptionsPattern[]] := Module[{values, vectors},
    {values, vectors} = Chop @ Simplify @ Enclose[
        ConfirmBy[
            If[ TrueQ[OptionValue[Chop]],
                If[ Precision[matrix] === MachinePrecision,
                    Quiet[
                        Check[
                            Eigensystem[matrix, ZeroTest -> (Chop[N[#1]] == 0 &)],
                            Eigensystem[matrix],
                            Eigensystem::eivec0
                        ],
                        Eigensystem::eivec0
                    ],
                    Eigensystem[Chop @ matrix]
                ],
                Eigensystem[matrix]
            ],
            MatchQ[{_ ? ListQ, _ ? MatrixQ}]
        ],
        Eigensystem[matrix] &
    ];
    If[ ! MatchQ[OptionValue["Sort"], False | None] && AllTrue[values, NumericQ],
        With[{ordering = OrderingBy[values, Replace[OptionValue["Sort"],
                True | Automatic :> If[Length[values] > 2 && ContainsOnly[Arg[values], {0, Pi}], Identity, {Mod[Arg[#], 2 Pi], Abs[#]} &]
            ]]
        },
            values = values[[ordering]];
            vectors = vectors[[ordering]]
        ]
    ];
    If[ TrueQ[OptionValue["Normalize"]], vectors = Normalize[If[NumericQ[First[#]] && First[#] != 0, # / First[#], #]] & /@ vectors];
    (* Eigensystem returns an arbitrary (generally non-orthonormal) basis within a degenerate
       eigenspace, so the eigenvectors need not resolve the identity. Gram-Schmidt in the
       (sorted) eigenvalue order repairs each degenerate block; distinct-eigenvalue vectors of
       a Hermitian matrix are already orthogonal, so it leaves a well-conditioned basis
       essentially unchanged. Numeric bases only. *)
    If[ TrueQ[OptionValue["Orthogonalize"]] && ArrayQ[vectors, 2, NumericQ], vectors = Orthogonalize[vectors]];

    {values, vectors}
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
   differentiable at the eigenvalues. *)
matrixFunction[f : Plus | Minus | Times | Conjugate, mat_, {left___}, {right___}, ___] := f[left, mat, right]

matrixFunction[Power, mat_, {left___}, {right___}, opts : OptionsPattern[]] := MatrixPower[mat, left, right, opts]

matrixFunction[f_, mat_, {left___}, {right___}, opts : OptionsPattern[]] := scalarMatrixFunction[f[left, #, right] &, mat, opts]

(* f must be finite on numeric eigenvalues: Log at a zero eigenvalue is a Failure,
   not a matrix of infinities. *)
spectralValues[f_, eigenvalues_] := With[{values = f /@ eigenvalues},
    If[ ! VectorQ[eigenvalues, NumericQ] || VectorQ[values, NumericQ],
        values,
        Failure["NonFiniteMatrixFunction", <|"MessageTemplate" -> "The function is not finite at an eigenvalue."|>]
    ]
]

(* An inexact diagonal is held to the same roundoff rule as a dense matrix, so f
   commutes with a change of basis near a zero eigenvalue too. *)
scalarMatrixFunction[f_, mat_, ___] /; SquareMatrixQ[mat] && DiagonalMatrixQ[mat] := Enclose @ With[
    {eigenvalues = Normal[Diagonal[mat]]},
    SparseArray[
        Band[{1, 1}] -> Confirm[spectralValues[f,
            If[MatrixQ[mat, NumericQ] && Precision[mat] < Infinity, roundoffEigenvalues[eigenvalues, 10 ^ -Precision[mat]], eigenvalues]
        ]],
        Dimensions[mat]
    ]
]

scalarMatrixFunction[f_, mat_, opts___] /; SquareMatrixQ[mat] && MatrixQ[mat, NumericQ] && Precision[mat] < Infinity :=
    With[{m = Normal[mat], eps = 10 ^ -Precision[mat]},
        inexactMatrixFunction[f, m, eps, 100 Length[m] eps, opts]
    ]

(* A derivative at numeric arguments with no numeric value, Derivative[1][Abs][1] at
   a Jordan block, means f is not differentiable where the matrix needs it; one that
   is merely unevaluated, Derivative[1][Zeta][2], has a value. *)
scalarMatrixFunction[f_, mat_, opts___] := Enclose @ ConfirmBy[
    ResourceFunction["ComputeMatrixFunction"][f, mat, opts],
    FreeQ[#, Indeterminate | _DirectedInfinity | (d : Derivative[__][_][__ ? NumericQ] /; ! NumericQ[N[d]])] &,
    "The function is not finite, or not differentiable where the matrix needs it, at an eigenvalue."
]

(* eps is the relative precision of the entries; within tol = 100 n eps the matrix
   counts as Hermitian and the strict upper triangle of t as zero. An eigenvalue
   within 10 eps of zero, relative to the largest, is set to zero: roundoff alone
   puts it there, and f may be singular at zero (Sqrt, Log). *)
roundoffEigenvalues[eigenvalues_, eps_] := Chop[eigenvalues, 10 eps Max[Abs[eigenvalues]]]

(* The Schur factor q is unitary even inside a degenerate eigenspace, where the
   eigenvectors Eigensystem returns need not be orthonormal, so f(m) = q.f(t).q^†
   when t is diagonal to roundoff. q.(x q^†) is q.DiagonalMatrix[x].q^† without the
   dense diagonal product. A Hermitian matrix keeps real eigenvalues and, for real
   f values, gives an exactly Hermitian result, real for real input. HermitianMatrixQ's
   Tolerance zeroes small entries rather than bounding m - m^†, so the test is
   written out. *)
inexactMatrixFunction[f_, mat_, eps_, tol_, opts___] := Enclose @ With[
    {qt = SchurDecomposition[mat, RealBlockDiagonalForm -> False]},
    {q = First[qt], t = Last[qt], hermitianQ = Max[Abs[mat - ConjugateTranspose[mat]]] <= tol Max[Abs[mat]]},
    If[ hermitianQ || Max[Abs[UpperTriangularize[t, 1]]] <= tol Max[Abs[t]],
        With[
            {values = Confirm[spectralValues[f, roundoffEigenvalues[If[hermitianQ, Re, Identity][Diagonal[t]], eps]]]},
            {result = q . (values ConjugateTranspose[q])},
            Which[
                ! hermitianQ, realOnRealInput[f, mat, Diagonal[t], values, tol, result],
                ! FreeQ[values, _Complex], result,
                FreeQ[mat, _Complex], Re[(result + Transpose[result]) / 2],
                True, (result + ConjugateTranspose[result]) / 2
            ]
        ],
        nonNormalMatrixFunction[f, mat, Eigensystem[mat], eps, opts]
    ]
]

(* A real matrix has eigenvalues in conjugate pairs, and when f maps each conjugate
   eigenvalue to the conjugate value (Cos, Exp, Log off the negative axis), f of it
   is real: its imaginary part is roundoff and is dropped. *)
realOnRealInput[f_, mat_, eigenvalues_, values_, tol_, result_] := If[
    FreeQ[mat, _Complex] && With[{conjugateValues = f /@ Conjugate[eigenvalues]},
        VectorQ[conjugateValues, NumericQ] &&
            Max[Abs[conjugateValues - Conjugate[values]]] <= tol Max[1, Max[Abs[values]]]
    ],
    Re[result],
    result
]

(* An eigenbasis v whose condition number stays below eps^(-1/4) keeps the error of
   v.f(d).v^-1 near eps^(3/4). A defective or nearly defective matrix fails this and
   goes to MatrixFunction (Schur-Parlett), which needs f differentiable there. *)
nonNormalMatrixFunction[f_, mat_, {eigenvalues_, vectors_}, eps_, ___] /;
    With[{sv = SingularValueList[vectors]}, Length[sv] == Length[vectors] && Max[sv] <= eps ^ (-1/4) Min[sv]] :=
    Enclose @ With[{values = Confirm[spectralValues[f, eigenvalues]]},
        realOnRealInput[f, mat, eigenvalues, values, 100 Length[mat] eps,
            Transpose[vectors] . (values Inverse[Transpose[vectors]])
        ]
    ]

nonNormalMatrixFunction[f_, mat_, {eigenvalues_, _}, _, opts___] :=
    Enclose[Confirm[spectralValues[f, eigenvalues]]; MatrixFunction[f, mat, opts]]

(* 0^m is the limit of b^m = MatrixExp[Log[b] m] as b -> 0. On an eigenvalue t the
   limit of b^t is 1 at t = 0 and 0 for Re t > 0, and there is none for any other t.
   A Jordan block adds the terms Log[b]^k b^t, which vanish for Re t > 0 and diverge
   at t = 0. So the limit exists exactly when every eigenvalue is zero or has
   positive real part and the zero eigenvalue is semisimple, with as many zero
   eigenvalues as null vectors, and it is then the spectral projector onto the null
   space of m along its range: x.(y.x)^-1.y for columns x spanning the null space of
   m and rows y spanning the null space of m acting from the left. For a number
   operator it is the vacuum projector. Each entry of a diagonal m is an eigenvalue
   with its own null vector, so there the limit is taken entry by entry. *)
zeroBasePower[mat_] /; SquareMatrixQ[mat] && DiagonalMatrixQ[mat] && MatrixQ[mat, NumericQ] := With[
    {diagonal = Normal[Diagonal[mat]]},
    {masks = zeroBaseMasks[diagonal, zeroBaseScale[mat, Max[Abs[diagonal]]]]},
    If[ Total[Last[masks]] > 0,
        noLimitFailure[First[Pick[diagonal, Last[masks], 1]]],
        zeroBaseMatrix[mat, Band[{1, 1}] -> atPrecisionOf[mat, First[masks]]]
    ]
]

zeroBasePower[mat_] /; SquareMatrixQ[mat] && MatrixQ[mat, NumericQ] && Precision[mat] === Infinity := With[
    {m = Normal[mat]},
    {eigenvalues = Eigenvalues[m], x = Transpose[NullSpace[m]], y = Conjugate[NullSpace[ConjugateTranspose[m]]]},
    {masks = zeroBaseMasks[eigenvalues, 0]},
    Which[
        Total[First[masks]] != Length[y], defectiveZeroFailure,
        Total[Last[masks]] > 0, noLimitFailure[First[Pick[eigenvalues, Last[masks], 1]]],
        y === {}, SparseArray[{}, Dimensions[m]],
        True, x . Inverse[y . x] . y
    ]
]

(* An inexact m. Roundoff moves each singular value by at most of order eps ||m||, so
   those within scale = 100 n eps ||m|| count as zero and give the nullity k, and the
   singular vectors give x and y. The zero eigenvalue is semisimple to within
   roundoff when the smallest singular value of y.x (the cosine of the largest angle
   between the two null spaces) stays above 100 n eps; the computed zero eigenvalues
   then lie within scale over that cosine of zero, so they are the k eigenvalues of
   smallest modulus. *)
zeroBasePower[mat_] /; SquareMatrixQ[mat] && MatrixQ[mat, NumericQ] := With[
    {m = Normal[mat], tol = 100 Length[mat] 10 ^ -Precision[mat]},
    {usv = SingularValueDecomposition[m]},
    {scale = tol Max[Diagonal[usv[[2]]]]},
    {k = Count[Diagonal[usv[[2]]], _ ? (# <= scale &)]},
    {x = Take[usv[[3]], All, -k], y = ConjugateTranspose[Take[usv[[1]], All, -k]]},
    {cosine = If[k == 0, 1, Min[SingularValueList[y . x, Tolerance -> 0]]]},
    If[ cosine <= tol,
        defectiveZeroFailure,
        inexactZeroBaseProjector[m, x, y, SortBy[Eigenvalues[m], Abs], k, scale, scale / cosine]
    ]
]

(* Every eigenvalue besides the k zero ones must lie farther from zero than the bound
   the zero ones keep to, or the zero eigenvalue is defective to within roundoff, and
   must have real part above scale, or the limit does not exist. A Hermitian m gives
   an exactly Hermitian projector. *)
inexactZeroBaseProjector[m_, x_, y_, eigenvalues_, k_, scale_, bound_] := With[
    {zeros = Take[eigenvalues, k], others = Drop[eigenvalues, k]},
    {noLimit = SelectFirst[others, Re[#] <= scale &]},
    Which[
        Max[Abs[zeros], 0] > bound || Min[Abs[others], Infinity] <= bound, defectiveZeroFailure,
        ! MissingQ[noLimit], noLimitFailure[noLimit],
        k == 0, zeroBaseMatrix[m, {}],
        True, With[{p = x . Inverse[y . x] . y},
            If[Max[Abs[m - ConjugateTranspose[m]]] <= scale, (p + ConjugateTranspose[p]) / 2, p]
        ]
    ]
]

(* A symbolic diagonal entry t keeps 0^t, the limit wherever Re t > 0, and one that is
   identically zero gives 1. *)
zeroBasePower[mat_] /; SquareMatrixQ[mat] && DiagonalMatrixQ[mat] := With[
    {diagonal = Normal[Diagonal[mat]]},
    {limits = zeroBaseLimitFunction /@ diagonal},
    If[ FreeQ[limits, Indeterminate],
        SparseArray[Band[{1, 1}] -> (limits /. zeroBaseLimitFunction[t_] :> 0 ^ t), Dimensions[mat]],
        noLimitFailure[First[Pick[diagonal, limits, Indeterminate]]]
    ]
]

zeroBasePower[mat_] := matrixFunction[zeroBaseLimitFunction, mat, {}, {}] /.
    (zeroBaseLimitFunction | zeroBaseLimitDerivative)[t_] :> 0 ^ t

(* 1 where an eigenvalue counts as zero, and 1 where it is neither zero nor of positive
   real part: exactly when scale is 0, and to within scale otherwise. An inexact
   diagonal is held to the same scale as a dense matrix, so that 0^m commutes with a
   change of basis. *)
zeroBaseMasks[eigenvalues_, 0] := With[{zero = Boole[PossibleZeroQ /@ eigenvalues]},
    {zero, (1 - zero) (1 - Boole[Positive[Re[eigenvalues]]])}
]

zeroBaseMasks[eigenvalues_, scale_] := With[{zero = UnitStep[scale - Abs[eigenvalues]]},
    {zero, (1 - zero) UnitStep[scale - Re[eigenvalues]]}
]

zeroBaseScale[mat_, norm_] := If[Precision[mat] === Infinity, 0, 100 Length[mat] 10 ^ -Precision[mat] norm]

(* A symbolic m goes through ComputeMatrixFunction, which needs the scalar function
   and, at a repeated root of the minimal polynomial, its derivatives. The
   derivatives of the limit are the limits of Log[b]^k b^t, 0 for Re t > 0 and none
   at t = 0 or elsewhere, which differentiating 0^t does not give: its derivative
   -Infinity 0^t has no value at any numeric t. So the function and its derivatives
   stay unevaluated at a symbolic t while ComputeMatrixFunction runs, and become 0^t
   afterwards, the limit wherever Re t > 0. At a numeric eigenvalue with no limit the
   closed form holds Indeterminate and the route fails as the other symbolic matrix
   functions do; a symbolic eigenvalue keeps 0^t, which has no value where the limit
   does not exist. *)
zeroBaseLimitFunction[t_ ? NumericQ] := Which[PossibleZeroQ[t], 1, TrueQ[Re[t] > 0], 0, True, Indeterminate]

zeroBaseLimitFunction[t_] /; TrueQ[PossibleZeroQ[t]] := 1

zeroBaseLimitFunction /: Derivative[_Integer ? Positive][zeroBaseLimitFunction] := zeroBaseLimitDerivative

zeroBaseLimitDerivative[t_ ? NumericQ] := If[TrueQ[Re[t] > 0] && ! PossibleZeroQ[t], 0, Indeterminate]

zeroBaseLimitDerivative[t_] /; TrueQ[PossibleZeroQ[t]] := Indeterminate

zeroBaseLimitDerivative /: Derivative[_Integer ? Positive][zeroBaseLimitDerivative] := zeroBaseLimitDerivative

(* The limit as a SparseArray at the precision of m: a machine m gives machine zeros,
   and an arbitrary-precision m keeps exact zeros, as N does. *)
zeroBaseMatrix[m_, entries_] := SparseArray[entries, Dimensions[m], atPrecisionOf[m, 0]]

atPrecisionOf[m_, x_] := If[MatrixQ[m, NumericQ] && Precision[m] < Infinity, N[x, Precision[m]], x]

noLimitFailure[t_] := Failure["ZeroBasePowerNoLimit", <|
    "MessageTemplate" -> "0^m is the limit of b^m as b -> 0, which does not exist: m has the eigenvalue `1`, neither zero nor with positive real part.",
    "MessageParameters" -> {t}
|>]

defectiveZeroFailure = Failure["ZeroBasePowerDefective", <|
    "MessageTemplate" -> "0^m is the limit of b^m as b -> 0, which does not exist: the zero eigenvalue of m is defective, where b^m grows as Log[b]."
|>]


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

(* Memoize[f] rewrites f so its results are cached. Results are stored in a private
   Association keyed on Hash[Hold[args]] (with the held key kept alongside for collision
   safety), instead of appending literal f[args] -> res DownValues. Because the cache key
   is an integer hash and never a pattern, this never trips Rule::rhs even when arguments
   contain Blank/Pattern (the failure the old QuantumOperatorProp cache had to Quiet).

   Memoize[f, "CacheQ" -> pred] only caches calls for which pred @@ args is True; pred is
   evaluated once on a cache miss (never on a hit), so an expensive cacheability predicate
   does not tax repeated reads. The original definitions are preserved on a shadow symbol
   whose bodies still call f, so recursion stays memoized.
   TODO: WFR. *)
Options[Memoize] = {"CacheQ" -> Automatic};

Memoize[f_ ? Developer`SymbolQ, opts : OptionsPattern[]] := With[{g = Unique[f], cache = Unique[f]},
    cache = <||>;
    SetAttributes[g, HoldAll];
    DownValues[g] = MapAt[ReplaceAll[f -> g], DownValues[f], {All, 1}];
    With[{cacheQ = Replace[OptionValue[Memoize, {opts}, "CacheQ"], Automatic -> (True &)]},
        ResourceFunction["BlockProtected"][{f},
            DownValues[f] = {
                HoldPattern[f[args___]] :> With[{held = Hold[args]},
                    Module[{key = Hash[held], hit},
                        hit = Lookup[cache, key, Missing[]];
                        If[ ! MissingQ[hit] && First[hit] === held,
                            Last[hit],
                            With[{res = g[args]},
                                (* a failed result (Failure, $Failed, or an abort) is
                                   never cached: its Message side effects would be
                                   swallowed on a hit, and an abort is not an answer *)
                                If[ TrueQ[cacheQ[args]] && ! FailureQ[res], cache[key] = {held, res}];
                                res
                            ]
                        ]
                    ]
                ]
            }
        ]
    ]
]

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

