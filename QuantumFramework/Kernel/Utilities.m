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
PackageScope["valuelessEntriesQ"]
PackageScope["exactZeroQ"]
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

(* A machine matrix is diagonalized with no ZeroTest: an explicit ZeroTest takes
   it off the numerical eigensolvers, although the Eigensystem documentation
   limits that option to exact and symbolic matrices. Of the two solvers the
   default method chooses between, only the Hermitian one guarantees orthonormal
   eigenvectors inside a degenerate eigenspace, and it runs only when
   HermitianMatrixQ holds at its default tolerance, which roundoff alone can
   fail; so a matrix Hermitian to roundoff is made exactly Hermitian first. Any
   other matrix goes to the general solver as it is. *)
machineEigensystem[m_] := Eigensystem[If[nearlyHermitianQ[m, roundoffTolerance[m]], (m + ConjugateTranspose[m]) / 2, m]]

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
            If[MatrixQ[mat, NumericQ] && Precision[mat] < Infinity, roundoffEigenvalues[eigenvalues, roundoff[mat]], eigenvalues]
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

(* The unit roundoff 10^-p of a matrix of precision p, computed at 20 digits so that it
   does not underflow for p above 307. *)
roundoff[mat_] := 10 ^ -SetPrecision[Precision[mat], 20]

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
   f values, gives an exactly Hermitian result, real for real input. *)
inexactMatrixFunction[f_, mat_, eps_, tol_, opts___] := Enclose @ With[
    {qt = SchurDecomposition[mat, RealBlockDiagonalForm -> False]},
    {q = First[qt], t = Last[qt], hermitianQ = nearlyHermitianQ[mat, tol]},
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
        realOnRealInput[f, mat, eigenvalues, values, roundoffTolerance[mat],
            Transpose[vectors] . (values Inverse[Transpose[vectors]])
        ]
    ]

(* MatrixFunction takes divided differences of f between eigenvalues that roundoff has
   split apart around zero, and those stay finite where f has no derivative at zero, so
   Sqrt and Log of a matrix within roundoff of a nilpotent one can come out as large
   numbers that mean nothing. So the eigenvalue 0 is read to within roundoff as 0^m reads
   it, and f fails when it has no finite value there or, as its series shows, no finite
   derivative of an order below a bound on the size of the largest Jordan block there. *)
nonNormalMatrixFunction[f_, mat_, {eigenvalues_, _}, _, opts___] := Enclose[
    Confirm[spectralValues[f, eigenvalues]];
    If[ lacksDerivativeAtZeroQ[f, zeroJordanBlockBound[mat, roundoffZeroEigenvalue[mat, eigenvalues]]],
        derivativeFailure,
        Replace[MatrixFunction[f, mat, opts], Except[_ ? (MatrixQ[#, NumericQ] &)] -> derivativeFailure]
    ]
]

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
        exactHermitianQ[m] && FreeQ[settledValue[#, Identity] & /@ Select[DeleteDuplicates[Flatten[m]], ! gaussianRationalQ[#] &], Missing["Overflow"]], exactHermitianZeroBasePower[m],
        True, exactZeroBasePower[m]
    ]]

inexactMatrixQ[mat_] := MatrixQ[mat, NumericQ] && Precision[mat] < Infinity

(* m with every numeric part that is zero but not written as 0 written as 0, decided as
   exactZeroQ decides, so that a number the convention there reads as zero, as
   Exp[-10^8], is zero throughout what follows; each distinct part is decided once, and
   a part that is zero is written as 0 whole. *)
zeroesWritten[m_] := With[
    {zero = AssociationMap[exactZeroQ, DeleteDuplicates[Cases[m, x_ /; ! AtomQ[x] && NumericQ[x], {2, Infinity}]]]},
    m /. x_ /; TrueQ[zero[x]] :> 0
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
        SparseArray[Band[{1, 1}] -> N[zero, Precision[mat]], Dimensions[mat], N[0, Precision[mat]]]
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
    With[{v = If[gaussianRationalQ[#], #, settledValue[#, Identity]]}, ! MissingQ[v] && 10 ^ -290 < Abs[v] < 10 ^ 290] &
]

(* The eigenvalues of the exact Hermitian m at precision p in increasing order, and the
   bound within which each lies of an exact one. *)
hermitianSpectrum[m_, p_] := With[
    {mp = roundedEntries[m, p]},
    {eigenvalues = Sort[Re[Eigenvalues[mp]]]},
    <|"Eigenvalues" -> eigenvalues, "Bound" -> roundoffTolerance[mp] Max[Abs[eigenvalues]], "Precision" -> p|>
]

(* The entries of m, whose zeros are written as 0, rounded to p digits of its largest
   entry. An entry that is not a Gaussian rational is first evaluated to an accuracy of
   3 more digits than that, with as much working precision as its cancellations need,
   so that each entry is off by rounding alone. (An accuracy goal is reached on an
   entry that holds a zero not written as 0, such as 4/3 + (GoldenRatio -
   (1 + Sqrt[5]) / 2), where a precision goal is not.) *)
roundedEntries[m_, p_] := With[
    {entries = DeleteDuplicates[Cases[m, Except[_ ? gaussianRationalQ], {2}]]},
    {accuracy = Max[N[p], 20] + 3 - Floor[Log10[Max[Abs[Cases[m, _ ? gaussianRationalQ, {2}]], Abs[Replace[settledValue[#, Identity] & /@ entries, _Missing -> 0, {1}]], 10 ^ -3200]]]},
    N[Replace[m, Dispatch[(# -> valueAt[#, Identity, accuracy]) & /@ entries], {2}], p]
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
corePolynomial[m_, k_] := With[{c = CoefficientList[CharacteristicPolynomial[m, \[FormalZ]], \[FormalZ]]},
    Together[(-1) ^ k Drop[c, k]] . \[FormalZ] ^ Range[0, Length[c] - k - 1]
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
    hurwitzStableQ[Reverse[Together[q] (-1) ^ Range[0, Length[q] - 1]]]
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

routhStep[{upper_, lower_}] := {lower, Together[Append[Rest[upper] - First[upper] / First[lower] Rest[lower], 0]]}

lowerLeadingPositiveQ[{_, lower_}] := realPartSign[First[lower]] === 1

(* An eigenvalue without positive real part among the roots of the characteristic
   polynomial p, named exactly when it is a root of a factor with real coefficients,
   decided by value, or of a linear or quadratic one, and Missing[] when such
   eigenvalues are roots only of factors with complex coefficients of degree 3 or more.
   (The kernel crashes building the exact roots of such a factor of degree 32, as
   Eigenvalues of a dense complex 32 x 32 matrix does.) A p of degree 3 or more is
   factored when Element proves its coefficients algebraic and they hold at most 3
   distinct radicals or Root objects; FactorList fails on others, as Exp[10^7], or
   takes minutes, as with 8 square roots, and p is then its own factor. *)
namedNoLimitEigenvalue[p_] := SelectFirst[
    Catenate[SolveValues[# == 0, \[FormalZ]] & /@ Select[
        If[ Exponent[p, \[FormalZ]] >= 3 && AllTrue[CoefficientList[p, \[FormalZ]], TrueQ[Element[#, Algebraics]] &] &&
                Length[DeleteDuplicates[Cases[CoefficientList[p, \[FormalZ]], Power[_Integer | _Rational, _Rational] | _Root | _AlgebraicNumber, {0, Infinity}]]] <= 3,
            FactorList[p, Extension -> Automatic][[All, 1]],
            {p}
        ],
        Exponent[#, \[FormalZ]] >= 1 && (Exponent[#, \[FormalZ]] <= 2 || AllTrue[CoefficientList[#, \[FormalZ]], realPartSign[I #] === 0 &]) &
    ]],
    realPartSign[#] =!= 1 &
]

(* Exact numbers are decided by value. The value of x at 50, then 400 and 3200 digits of
   precision or of accuracy, whichever N reaches first, tells a nonzero x from zero, and
   a sign read off it is exact. (The accuracy goal is reached on a number that holds a
   zero not written as 0, and the precision goal on a huge number without all of its
   digits.) What stays within 10^-3200 of zero is zero when Element proves it
   algebraic, after FunctionExpand, and PossibleZeroQ with the method ExactAlgebraics,
   which is exact for algebraic numbers, finds it zero, as for
   GoldenRatio - (1 + Sqrt[5]) / 2, and is taken as zero when Element does not, as
   for Log[6] - Log[2] - Log[3]; a number that is not zero, lies within 10^-3200 of
   it and is not proved algebraic, as Tanh[10^4] - 1, is read as zero. A real part
   proved algebraic and found nonzero is read to digits that keep doubling until they
   tell, as they must, up to 12 doublings. These decisions are exact for numbers
   between $MinNumber and $MaxNumber, the range of arbitrary-precision numbers. Beyond
   it N reports an underflow or an overflow with its own messages: a value that
   underflows is read as zero, proved algebraic or not (PossibleZeroQ with the method
   ExactAlgebraics crashes the kernel on some, as (Sqrt[2] - 1)^(10^20) -
   (Sqrt[2] - 1)^(10^20 + 1)), and one that overflows as nonzero, with the sign of
   its real part as Sign reads it from the expression, as for -Exp[Exp[100]], and
   Indeterminate when Sign cannot. *)
exactZeroQ[x_ ? gaussianRationalQ] := x == 0

exactZeroQ[x_] := Replace[settledValue[x, Identity], {
    Missing["Overflow"] -> False,
    Missing["Underflow"] -> True,
    _Missing :> Replace[provenAlgebraic[x], {_Missing -> True, y_ :> PossibleZeroQ[y, Method -> "ExactAlgebraics"]}],
    _ -> False
}]

(* The sign of the real part of the exact number x. *)
realPartSign[x_ ? gaussianRationalQ] := Sign[Re[x]]

realPartSign[x_] := Replace[settledValue[x, Re], {
    Missing["Overflow"] :> Replace[Sign[Re[x]], Except[-1 | 0 | 1] -> Indeterminate],
    Missing["Underflow"] -> 0,
    _Missing :> With[{re = Replace[provenAlgebraic[x], {_Missing :> Re[x], y_ :> Re[y]}]},
        If[exactZeroQ[re], 0, Sign[digitsValue[x, Re, NestWhile[2 # &, 6400, digitsValue[x, Re, #] == 0 &, 1, 12]]]]
    ],
    v_ :> Sign[v]
}]

(* x as FunctionExpand writes it when Element then proves it algebraic, as
   Sin[ArcCos[1/3] / 2] = 1 / Sqrt[3]; Missing[] when Element does not. *)
provenAlgebraic[x_] := With[{y = FunctionExpand[x]}, If[TrueQ[Element[y, Algebraics]], y, Missing[]]]

(* The part of the value of x, as read by part, at the first of 50, 400 and 3200 digits
   that tells it from zero; Missing[] when none does, and Missing["Underflow"] or
   Missing["Overflow"] when the value underflows or overflows. *)
settledValue[x_, part_] := firstNonzeroValue[x, part, {50, 400, 3200}]

firstNonzeroValue[_, _, {}] := Missing[]

firstNonzeroValue[x_, part_, {digits_, rest___}] := With[{v = digitsValue[x, part, digits]},
    Which[
        v === Overflow[], Missing["Overflow"],
        v === Underflow[], Missing["Underflow"],
        v != 0, v,
        True, firstNonzeroValue[x, part, {rest}]
    ]
]

(* The part of the value of x at the given digits of precision or of accuracy, whichever
   N reaches first. *)
digitsValue[x_, part_, digits_] := Block[{$MaxExtraPrecision = extraPrecision[x, digits]}, part[N[x, {digits, digits}]]]

(* The part of the value of x to the given accuracy. *)
valueAt[x_, part_, accuracy_] := Block[{$MaxExtraPrecision = extraPrecision[x, accuracy]}, part[N[x, {Infinity, accuracy}]]]

(* The working precision N may add to reach its goal: enough for the cancellation that
   integers of d digits in x can cause, and finite. Under an unlimited one N never
   finishes on a number that holds an exact zero where a function jumps, as
   Sqrt[-1 + (GoldenRatio - (1 + Sqrt[5]) / 2) I] - I does on the branch cut of
   Sqrt. A number whose cancellation needs more, as between terms far beyond
   10^10000 in size, is not told from zero: it is read as zero, as one within
   10^-3200 of it is, and N says so with its own message. *)
extraPrecision[x_, digits_] := 10 Max[digits, 0] + 10000 + 2 Max[0, Cases[x,
    r : _Integer | _Rational | _Complex :> IntegerLength[Max[Abs[Numerator[{Re[r], Im[r]}]], Denominator[{Re[r], Im[r]}]]],
    {0, Infinity}
]]

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
        zero["Nullity"] == 0, SparseArray[{}, Dimensions[m], N[0, Precision[m]]],
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

