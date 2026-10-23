Package["Wolfram`QuantumFramework`"]

PackageExport["QuantumEntangledQ"]
PackageExport["QuantumEntanglementMonotone"]
PackageScope["$QuantumEntanglementMonotones"]
PackageScope["numericStateNotPSDQ"]



$QuantumEntanglementMonotones = {
    "Concurrence", "Negativity", "LogNegativity", "EntanglementEntropy", "RenyiEntropy", "Realignment",
    "MutualInformationI", "MutualInformationJ", "Discord"
}

QuantumEntanglementMonotone::mixedentropy =
    "The entanglement entropy of a subsystem measures entanglement only for a pure state; the reduced von Neumann or Renyi entropy of a mixed state is not an entanglement measure. Use Negativity, LogNegativity, or Concurrence for a mixed state."

QuantumEntanglementMonotone::notphysical =
    "The input is not a positive-semidefinite density matrix; an entanglement monotone of a non-physical state may be meaningless."

(* The physicality test reads numbers only, so a symbolic state is not tested. An exact state held in another frame is
   tested in the computational frame, simplified: the change of basis can leave an exact zero as a sum of radicals or
   roots of unity, on which the test would raise N::meprec. *)
warnUnphysicalMonotone[qs_] := If[
    qs["NumericQ"] && numericStateNotPSDQ[If[exactNumericStateQ[qs] && ! qs["ComputationalQ"], simplifiedExactState[qs["Computational"]], qs]],
    Message[QuantumEntanglementMonotone::notphysical]
]

QuantumEntangledQ::undecided =
    "The sign of `1` could not be decided exactly within `2` seconds, so the answer is Indeterminate."

QuantumEntangledQ::realparameters =
    "The parameters `1` are taken to be real, and only their values for which the state is a density matrix count: the answer holds for all of them, or it is Indeterminate."


(* QuantumEntangledQ answers True when the method's test proves the two-part state entangled, False when a
   separability proof holds for it, and Indeterminate when neither does. Automatic tests the partial transpose in
   every size, for a negative eigenvalue, and beyond 2 x 3, where no efficient test is complete, the realignment
   value ||R(rho)||_1 - 1 as well, which detects some states with a positive partial transpose. Its separability
   proofs: a positive partial transpose where that test is exact (pptSeparableQ), and a purity Tr[rho^2] of at most
   1 / (dA dB - 1), the largest ball of states around the maximally mixed one, all of them separable (Gurvits and
   Barnum, PRA 66, 062311 (2002)). A named method answers with its own test alone: Negativity and LogNegativity
   with the partial transpose and the proofs it gives, the mutual informations and the discord, positive on
   separable states too, never, and any other method with a positive value read from QuantumEntanglementMonotone
   (the realignment value within the time limit that exact and symbolic input get), or a zero value where that
   method is exact. A symbolic input raises a message naming its parameters, taken to be real; a symbolic state whose
   physical parameter values cannot be found is Indeterminate. *)
QuantumEntangledQ[qs_ ? QuantumStateQ, biPartition_ : Automatic, method : _String | Automatic : Automatic] /;
        method === Automatic || MemberQ[$QuantumEntanglementMonotones, method] :=
    Enclose @ With[{q = exactSymbolicState[qs]},
        If[! q["NumericQ"], Message[QuantumEntangledQ::realparameters, DeleteDuplicates @ symbolicParameters[Normal @ q["State"]]]];
        With[{bp = ConfirmBy[q["Bipartition", biPartition]["Normalized"], QuantumStateQ[#] && #["Qudits"] == 2 &]},
            With[{c = numberContext[bp]},
                If[c["Kind"] === "Symbolic" && ! c["Exists"], Indeterminate, entangledAnswer[method, q, biPartition, bp, c]]
            ]
        ]
    ]

(* A symbolic state with approximate numbers is read with them at their exact binary values, so that a product of
   small coefficients neither underflows nor drops below round-off before the symbolic decision sees it. *)
exactSymbolicState[qs_] := If[! qs["NumericQ"] && Precision[qs["State"]] < Infinity, QuantumState[Rationalize[Normal @ qs["State"], 0], qs["Basis"]], qs]

(* A separability proof on an approximate density matrix rests on cuts at round-off, so there the realignment, a proof
   of entanglement, is read first. *)
entangledAnswer[Automatic, qs_, _, bp_, c_] := (
    warnUnphysicalMonotone[qs];
    With[{npt = partialTransposeNegativeQ[bp, c]},
        Which[
            TrueQ[npt], True,
            c["Kind"] === "Approximate" && bp["MatrixQ"],
            Which[realignmentDetectsQ[bp, c], True, separabilityProofQ[npt, bp, c], False, True, Indeterminate],
            separabilityProofQ[npt, bp, c], False,
            realignmentDetectsQ[bp, c], True,
            True, Indeterminate
        ]
    ]
)

separabilityProofQ[npt_, bp_, c_] := npt === False && pptSeparableQ[bp, c] || purityBoundQ[bp, c]

realignmentDetectsQ[bp_, c_] := Times @@ bp["Dimensions"] > 6 && criterionSign[realignmentCriterion[bp, c], c] === 1

entangledAnswer["Negativity" | "LogNegativity", qs_, _, bp_, c_] := (
    warnUnphysicalMonotone[qs];
    With[{npt = partialTransposeNegativeQ[bp, c]},
        Which[TrueQ[npt], True, npt === False && pptSeparableQ[bp, c], False, True, Indeterminate]
    ]
)

(* The correlation measures are positive on separable states too, so their values prove nothing. *)
entangledAnswer["MutualInformationI" | "MutualInformationJ" | "Discord", ___] := Indeterminate

(* The realignment value is read with the time limit realignmentCriterion gives exact and symbolic input; zero proves
   separability on a pure state only. *)
entangledAnswer["Realignment", qs_, _, bp_, c_] := (
    warnUnphysicalMonotone[qs];
    With[{sign = criterionSign[realignmentCriterion[bp, c], c]},
        Which[sign === 1, True, sign === 0 && pureTwoPartQ[bp, c], False, True, Indeterminate]
    ]
)

(* A zero value proves separability where the method is exact: the concurrence on two qubits and on pure states,
   the realignment value and the entropies on pure states. An entropy is read in bits. *)
entangledAnswer[method_, qs_, biPartition_, bp_, c_] := With[{
    sign = criterionSign[Replace[QuantumEntanglementMonotone[qs, biPartition, method], Quantity[x_, "Bits"] :> x], c]
},
    Which[
        sign === 1, True,
        sign === 0 && (method === "Concurrence" && bp["Dimensions"] === {2, 2} || pureTwoPartQ[bp, c]), False,
        True, Indeterminate
    ]
]

(* Whether the partial transpose of the two-part state bp has a negative eigenvalue: True, False, or Indeterminate
   when that cannot be decided. A state vector has one exactly when its Schmidt rank exceeds 1: approximate input
   reads its negativity Sum[s_i s_j, i < j] off the Schmidt coefficients s, the singular values of its amplitudes as
   a dA x dB matrix, and exact input the sign of Tr[rhoA]^2 - Tr[rhoA^2], 0 on the zero vector. A density matrix is
   read through its partial transpose: the negativity of approximate input, the sign of every eigenvalue of exact
   input, and for symbolic input the coefficients of det(x + rho^T_B), all of them nonnegative exactly when no
   eigenvalue is negative, for every physical value of the parameters. *)
partialTransposeNegativeQ[bp_, c_] /; bp["VectorQ"] && c["Kind"] =!= "Symbolic" :=
    With[{m = ArrayReshape[Normal @ bp["StateVector"], bp["Dimensions"]]},
        If[ c["Kind"] === "Approximate",
            With[{s = SingularValueList[m, Tolerance -> 0]}, Chop[(Total[s] ^ 2 - s . s) / 2] > 0],
            Replace[quantitySign[With[{a = m . ConjugateTranspose[m]}, Tr[a] ^ 2 - Total[a Transpose[a], 2]]], {1 -> True, 0 -> False, _ -> Indeterminate}]
        ]
    ]

partialTransposeNegativeQ[bp_, c_] := With[{pt = partialTransposeMatrix[Normal @ bp["DensityMatrix"], bp["Dimensions"]]},
    Switch[c["Kind"],
        "Approximate",
        Chop[(Total[SingularValueList[pt]] - 1) / 2] > 0,
        "Exact",
        With[{signs = quantitySign /@ Eigenvalues[pt]}, Which[MemberQ[signs, -1], True, SubsetQ[{0, 1}, signs], False, True, Indeterminate]],
        "Symbolic",
        With[{e = positiveSemidefiniteCoefficients[pt, c["Parameters"]]},
            Which[! ListQ[e], Indeterminate, forAllQ[c, Or @@ Thread[e < 0]], True, forAllQ[c, And @@ Thread[e >= 0]], False, True, Indeterminate]
        ]
    ]
]

partialTransposeMatrix[m_, {dA_, dB_}] := ArrayReshape[Transpose[ArrayReshape[m, {dA, dB, dA, dB}], {1, 4, 3, 2}], {dA dB, dA dB}]

(* A positive partial transpose proves separability where the test is exact: on 2 x 2 and 2 x 3 (Horodecki,
   Horodecki and Horodecki, PLA 223, 1 (1996)), on a pure state, and on a density matrix whose rank is at most the
   larger of the ranks of its two marginals (Horodecki, Lewenstein, Vidal and Cirac, PRA 62, 032310 (2000), Thm 1;
   the rank against the local dimensions alone proves nothing, since a state can sit in a larger space). Approximate
   ranks are counted so that an eigenvalue of uncertain size never makes the proof: the rank of rho counts every
   eigenvalue above round-off, tol = n 10^(2 - p) at precision p, and the ranks of the marginals only those above
   10^-10. The eigenvalues rho drops lift an eigenvalue of a marginal by at most n tol, so the proof is given only
   where n tol <= 10^-10, which keeps them from adding to a marginal's rank: up to n = 100 at machine precision. An
   approximate density matrix is read as pure when its impurity 1 - Tr[rho^2] is at most tol, and only where
   tol <= 10^-10. The precision p is rounded up to a
   whole digit, so that the fraction of a digit a normalization costs never moves the cut, and a higher p only makes
   either proof harder to give. A symbolic density matrix is not ranked. *)
pptSeparableQ[bp_, c_] := Times @@ bp["Dimensions"] <= 6 || pureTwoPartQ[bp, c] || lowRankQ[bp, c]

pureTwoPartQ[bp_, c_] := bp["VectorQ"] || With[{m = Normal @ bp["DensityMatrix"]},
    Switch[c["Kind"],
        "Symbolic", forAllQ[c, realForm[Total[m Transpose[m], 2], c] == 1],
        "Exact", pureDensityMatrixQ[m],
        "Approximate", With[{tol = roundoffCut[m]}, tol <= 10 ^ -10 && Abs[1 - Total[m Transpose[m], 2]] <= tol]
    ]
]

lowRankQ[bp_, c_] := bp["MatrixQ"] && c["Kind"] =!= "Symbolic" && With[{m = Normal @ bp["DensityMatrix"], dims = bp["Dimensions"]},
    With[{t = ArrayReshape[m, Join[dims, dims]]},
        If[ c["Kind"] === "Exact",
            MatrixRank[m] <= Max[MatrixRank[TensorContract[t, {{2, 4}}]], MatrixRank[TensorContract[t, {{1, 3}}]]],
            With[{tol = roundoffCut[m]},
                Length[m] tol <= 10 ^ -10 &&
                    eigenvalueCount[m, tol] <= Max[eigenvalueCount[TensorContract[t, {{2, 4}}], 10 ^ -10], eigenvalueCount[TensorContract[t, {{1, 3}}], 10 ^ -10]]
            ]
        ]
    ]
]

eigenvalueCount[m_, tol_] := Count[Re @ Eigenvalues[(m + ConjugateTranspose[m]) / 2], x_ /; x > tol]

roundoffCut[m_] := Length[m] 10 ^ (2 - Ceiling[Precision[m]])

purityBoundQ[bp_, c_] := bp["MatrixQ"] && With[{m = Normal @ bp["DensityMatrix"]},
    nonnegativeQ[1 / (Times @@ bp["Dimensions"] - 1) - Total[m Transpose[m], 2], c]
]

(* The realignment value of the two-part state, given a time limit on exact and symbolic input, whose singular values
   may not come out in closed form. *)
realignmentCriterion[bp_, c_] := If[c["Kind"] === "Approximate", realignmentValue[bp], TimeConstrained[realignmentValue[bp], $entanglementSignTime, Indeterminate]]


(* The number type of the two-part state bp, with what deciding a sign on it needs: "Approximate", "Exact", or
   "Symbolic" with the parameters, taken as real, the condition that bp is a density matrix for their values:
   Hermitian, and every coefficient of det(x + rho) nonnegative, so that no eigenvalue is negative, and whether any
   value meets it. The condition is $Failed when it is not found within $entanglementSignTime seconds. *)
numberContext[bp_] := With[{a = If[bp["VectorQ"], bp["StateVector"], bp["DensityMatrix"]]},
    Which[
        ! ArrayQ[a, _, NumericQ], symbolicContext[Normal @ bp["DensityMatrix"]],
        Precision[a] === Infinity, <|"Kind" -> "Exact"|>,
        True, <|"Kind" -> "Approximate"|>
    ]
]

symbolicContext[rho_] := With[{ps = DeleteDuplicates @ symbolicParameters[rho]},
    With[{physical = With[{e = positiveSemidefiniteCoefficients[rho, ps]},
        If[ ListQ[e],
            TimeConstrained[
                Simplify[
                    And @@ Join[Thread[Flatten[ComplexExpand[{Re[#], Im[#]}] & [rho - ConjugateTranspose[rho]]] == 0], Thread[e >= 0]],
                    Element[ps, Reals]
                ],
                $entanglementSignTime,
                $Failed
            ],
            $Failed
        ]
    ]},
        <|"Kind" -> "Symbolic", "Parameters" -> ps, "Physical" -> physical, "Exists" -> resolvedQ[Exists[ps, physical]]|>
    ]
]

(* The parameters of a symbolic expression: its outermost parts that are not numeric and on which no numeric function
   acts, plain symbols and indexed ones such as x[1] or Subscript[p, 1]. *)
symbolicParameters[_ ? NumericQ] := {}
symbolicParameters[h_Symbol[args___]] /; h === List || MemberQ[Attributes[h], NumericFunction] := Catenate[symbolicParameters /@ {args}]
symbolicParameters[x_] := {x}

(* The coefficients e_k of det(x + h) = Sum[e_k x^(n - k)]: the elementary symmetric functions of the eigenvalues of
   h, all of them nonnegative exactly when no eigenvalue of the Hermitian h is negative; $Failed when they are not
   found within $entanglementSignTime seconds. *)
positiveSemidefiniteCoefficients[h_, ps_] := Module[{x},
    TimeConstrained[
        Simplify[ComplexExpand[(-1) ^ Length[h] CoefficientList[CharacteristicPolynomial[-h, x], x]], Element[ps, Reals]],
        $entanglementSignTime,
        $Failed
    ]
]

(* Whether the condition cond holds for every real value of the parameters at which the symbolic state is a density
   matrix. *)
forAllQ[c_, cond_] := c["Exists"] && With[{ps = c["Parameters"], physical = c["Physical"]}, resolvedQ[ForAll[ps, physical, cond]]]

(* Whether a statement about real parameters is true, decided by Resolve within $entanglementSignTime seconds: False
   when it is not decided, when a step that built it ran out of time ($Failed), or when it holds a nonreal constant.
   Its approximate numbers are taken at their exact binary values, none of them dropped however small, so that an
   imaginary part 0. becomes 0 and is not read as a nonreal constant. *)
resolvedQ[q_] := TrueQ @ TimeConstrained[
    With[{r = Rationalize[ComplexExpand[q], 0]}, FreeQ[r, $Failed | _Complex] && Resolve[r, Reals]],
    $entanglementSignTime,
    False
]

realForm[v_, c_] := TimeConstrained[Simplify[ComplexExpand[v], Element[c["Parameters"], Reals]], $entanglementSignTime, $Failed]


(* The sign, 1, 0 or -1, of a real criterion quantity, or Indeterminate when it cannot be decided. An approximate
   value within 10^-10 of zero counts as zero. An exact value is read in 60-digit significance arithmetic, by
   SetPrecision, which raises no N::meprec: N does, with a precision goal on an exact zero, and with an accuracy goal
   on a sum of large terms that cancel, which it can then read as zero. The reading gives the sign when any of its
   digits survives the cancellation; a value that keeps none is decided within $entanglementSignTime seconds: by
   readings at 200 and 1000 digits; then as zero when FullSimplify reduces x to 0 within a fifth of that time; then by
   a 5000-digit reading, which for a constant such as Khinchin's takes longer than the whole limit; then, when x is
   algebraic, by the sign of RootReduce[x], 0 for an exact zero. A number built from rationals and Root objects by
   sums, products, rational powers and Abs is algebraic, as is any that Element proves algebraic. Any other is
   Indeterminate, with a message, which shows x unevaluated, each integer or rational of more than 20 digits in it
   replaced by its 6-digit value. The quantity of a symbolic state has a sign when it has that sign for every
   physical value of the parameters. *)
$entanglementSignTime = 5;

criterionSign[v_, c_] := If[c["Kind"] === "Symbolic", symbolicSign[realForm[v, c], c], quantitySign[v]]

quantitySign[x_ ? NumericQ] /; Precision[x] < Infinity := Sign[Chop[Re[x]]]
quantitySign[x_ ? NumericQ] := With[{v = Re[SetPrecision[x, 60]]}, If[Precision[v] > 0, Sign[v], exactSign[x]]]
quantitySign[_] := Indeterminate

exactSign[x_] := Replace[
    TimeConstrained[
        With[{digits = SelectFirst[{200, 1000}, Precision[Re[SetPrecision[x, #]]] > 0 &]},
            Which[
                IntegerQ[digits], Sign[Re[SetPrecision[x, digits]]],
                TimeConstrained[FullSimplify[x], $entanglementSignTime / 5, x] === 0, 0,
                True, With[{v = Re[SetPrecision[x, 5000]]},
                    Which[
                        Precision[v] > 0, Sign[v],
                        algebraicNumberQ[x] || TrueQ[Element[x, Algebraics]], Block[{$MaxExtraPrecision = 5000}, Sign[RootReduce[x]]],
                        True, Indeterminate
                    ]
                ]
            ]
        ],
        $entanglementSignTime,
        Indeterminate
    ],
    Except[-1 | 0 | 1] :> (Message[QuantumEntangledQ::undecided, Short[HoldForm[x] /. q : _Integer | _Rational /; IntegerLength[Numerator[q]] + IntegerLength[Denominator[q]] > 21 :> With[{a = N[q, 6]}, a /; True], 2], $entanglementSignTime]; Indeterminate)
]

algebraicNumberQ[_Integer | _Rational | _Complex | _Root | _AlgebraicNumber] := True
algebraicNumberQ[(Plus | Times | Abs)[args__]] := AllTrue[{args}, algebraicNumberQ]
algebraicNumberQ[Power[a_, _Integer | _Rational]] := algebraicNumberQ[a]
algebraicNumberQ[_] := False

symbolicSign[v_, c_] := If[FreeQ[v, Alternatives @@ c["Parameters"]],
    quantitySign[v],
    Which[forAllQ[c, v > 0], 1, forAllQ[c, v == 0], 0, forAllQ[c, v < 0], -1, True, Indeterminate]
]

nonnegativeQ[v_, c_] := If[c["Kind"] === "Symbolic", forAllQ[c, realForm[v, c] >= 0], MatchQ[quantitySign[v], 0 | 1]]



(* The default method is decided on the two-part state, after the split and the partial trace: the concurrence for
   a pure state and for two qubits, where it is exact, and the negativity for a mixed state beyond two qubits, where
   the concurrence is only a lower bound that a local unitary can change. A method reaches this head either as a
   bare string in the second slot (the optional biPartition then defaults) or as the {"RenyiEntropy", alpha} list;
   both are excluded here so a method spec is never mistaken for a bipartition, which is always Automatic, an
   integer, or a list of qudit indices, never a method name. *)
QuantumEntanglementMonotone[qs_ ? QuantumStateQ,
    biPartition : Except[_String | {"RenyiEntanglementEntropy" | "RenyiEntropy", _}] : Automatic] :=
    Enclose @ QuantumEntanglementMonotone[qs, biPartition, defaultMonotone[ConfirmBy[qs["Bipartition", biPartition], QuantumStateQ]]]

defaultMonotone[bp_] := If[
    bp["Qudits"] == 2 && (bp["Dimensions"] === {2, 2} || bp["VectorQ"] || pureDensityMatrixQ[Normal @ bp["Normalized"]["DensityMatrix"]]),
    "Concurrence",
    "Negativity"
]


y[{j_Integer, k_Integer}, n_Integer] /; 1 <= j < k <= n := SparseArray[{{j, k} -> -I, {k, j} -> I}, {n, n}]

Y[n_] := Y[n] = Catenate @ Table[y[{j, k}, n], {k, 2, n}, {j, k - 1}]

(* Rungta-Buzek-Caves-Hillery-Milburn I-concurrence (PRA 64, 042315 (2001)): the concurrence-vector
   components are the per-generator-pair Wootters combinations lambda1 - Sum[rest] of the singular
   values of Sqrt[rho].Sqrt[rho-tilde], with the spin-flipped rho-tilde = o.Conjugate[rho].o built
   from the su(d) generator pair o = yn (x) ym. Those singular values are the square roots of the
   eigenvalues of rho.rho-tilde. Norm of the vector is exact for a pure state of any dimension
   (= Sqrt[2 (1 - Tr[rhoA^2])]) and exact for two qubits (the textbook Wootters concurrence); for a
   d1 d2 > 4 mixed state it is a lower bound on the convex-roof concurrence, which has no closed form,
   and (unlike the two-qubit Wootters concurrence) it is basis-dependent: not invariant under local
   unitaries, though it stays below the LU-invariant convex-roof value in every basis.

   An exact or symbolic reduction, in any dimension, reads the Wootters values off the spectrum of
   rho.rho-tilde, which is similar to the positive-semidefinite Sqrt[rho].rho-tilde.Sqrt[rho] and so real
   and non-negative. Eigenvalues needs only the characteristic polynomial, where the operator square root
   Sqrt[rho] (MatrixPower) needs the eigenvectors, which blow up into huge unsimplified Root objects on
   exact input and stall on symbolic input; exact stays exact. The real part of each eigenvalue is kept;
   a numeric one is clipped at zero and a symbolic one is left as it is, so that the closed form simplifies
   under real assumptions and evaluates to a real number at numeric parameter values.

   On approximate numbers the square roots of computed eigenvalues are not taken: round-off of order eps
   in a zero eigenvalue becomes order Sqrt[eps] in the value, about 1e-8 on machine numbers, which reads a
   product state as entangled. With rho = W.W^dagger, the eigenvalues of rho.rho-tilde are the squared
   singular values of the symmetric W^T.o.W, so the Wootters values are read as those singular values,
   whose round-off stays of order eps, two qubits included. *)

(* lambda1 - Sum[rest]: largest minus the rest, order-independent (the singular values are not guaranteed sorted for
   symbolic input). An exact value numerically near zero is reduced to one algebraic number first: square roots of Root
   objects that cancel to exactly zero, as for a separable mixture of three product states, cannot be compared with
   zero numerically even at $MaxExtraPrecision. *)
wootterCombination[lambda_] := With[{x = 2 Max[lambda] - Total[lambda]},
    Max[0, If[NumericQ[x] && Precision[x] === Infinity && Abs[N[x]] < 10^-6, RootReduce[x], x]]
]

approximateNumericArrayQ[a_] := ArrayQ[a, _, NumericQ] && Precision[a] < Infinity

(* A factor W with rho = W.W^dagger: the orthonormal eigenvectors of rho scaled by the square roots of its
   eigenvalues. An eigenvalue below n 10^-p of the largest, the eigensolver's own error at precision p, is set to
   zero: a zero eigenvalue of a rank-deficient rho comes out as round-off of order eps, and its square root, of
   order Sqrt[eps], would enter the Wootters values linearly. *)
densityMatrixFactor[mat_] := With[{es = Eigensystem[(mat + ConjugateTranspose[mat]) / 2]},
	With[{ev = Re @ es[[1]]},
		Transpose[Sqrt[ev UnitStep[ev - Length[mat] Max[ev] 10 ^ -Max[Precision /@ ev]]] es[[2]]]
	]
]

ConcurrenceVector[qs_ ? QuantumStateQ, biPartition_ : Automatic] :=
    Enclose @ concurrenceVector[ConfirmBy[qs["Bipartition", biPartition]["Normalized"], QuantumStateQ[#] && #["Qudits"] == 2 &]["Operator"]]

(* The concurrence vector of a normalized two-qudit operator rho. *)
concurrenceVector[rho_] := Block[{component, d1, d2, y1, y2},
	{d1, d2} = rho["OutputDimensions"];
	y1 = Y[d1];
	y2 = Y[d2];
	component = If[ approximateNumericArrayQ[rho["Matrix"]],
		With[{w = densityMatrixFactor[Normal @ rho["Matrix"]]},
			Function[{ya, yb}, wootterCombination @ SingularValueList[Transpose[w] . KroneckerProduct[ya, yb] . w, Tolerance -> 0]]
		],
		With[{mat = Normal @ rho["Matrix"]}, With[{cmat = Conjugate[mat]},
			Function[{ya, yb}, With[{omat = Normal @ KroneckerProduct[ya, yb]},
				wootterCombination @ Sqrt[If[NumericQ[#], Max[0, #], #] & /@ Re[Eigenvalues[mat . omat . cmat . omat]]]
			]]
		]]
	];
	Catenate @ Table[component[yn, ym], {yn, y1}, {ym, y2}]
]

QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, "ConcurrenceVector"] :=
    (warnUnphysicalMonotone[qs]; ConcurrenceVector[qs, biPartition])

(* The concurrence depends only on the state's direction, not its overall scale: the bipartition is
   normalized before the reduced purity is read, so an input with trace or vector-norm != 1 gives the
   same value as its normalized form. Every monotone in this file shares that contract (the mixed route
   normalizes inside ConcurrenceVector). Without it the reduced purity picks up the scale, and the
   Max[0, ...] clamp on the now-shifted 2 (1 - Purity) silently returns 0 for a scaled pure state.

   Whether the pure-state formula applies is decided on the two-part state, after the split and the partial
   trace over the qudits outside it, not on the input: tracing qudits out of a pure state leaves the pair
   mixed in general, as for any two qubits of a GHZ state, and a mixed pair takes the mixed-state route. The
   pair is pure when its purity Tr[rho^2], read as a sum over the entries with no size limit, is 1: exactly
   for exact entries, and within n 10^(2 - p) for approximate entries on an n-dimensional pair, p being the
   precision of the purity itself: an entry the partial trace cancels to zero has precision 0 and would make
   the precision of the matrix 0. That tolerance cannot tell a pure pair from one mixed by less than it, so an
   approximate density matrix of up to 64 dimensions always takes the mixed-state route, which is accurate on
   pure states as well; the purity test serves the larger pairs, where the (d1 d2)^2 / 4 generator pairs of
   that route are out of reach.

   A pure two-part state given by approximate numbers is read through its Schmidt coefficients s, the
   singular values of its amplitudes arranged as a d1 x d2 matrix: C = 2 Sqrt[Sum[s_i^2 s_j^2, i < j]], a
   sum of non-negative terms. 2 (1 - Purity) cancels to round-off of order eps on a product state, and its
   square root lifts that to order Sqrt[eps]; the Schmidt form keeps it of order eps. The zero vector, whose
   normalized form is still zero, has concurrence 0. A pure pair held as a
   density matrix rho = psi.psi^dagger gives its amplitudes, up to a global phase, as the column of rho with
   the largest diagonal entry, divided by the square root of that entry. *)
schmidtConcurrence[v_, dims_] := With[{p = SingularValueList[ArrayReshape[v, dims], Tolerance -> 0] ^ 2},
    If[Total[p] == 0, 0, 2 Sqrt[Rest[p] . Most[Accumulate[p]]] / Total[p]]
]

pureDensityMatrixQ[m_] := With[{p = Total[m Transpose[m], 2]},
    If[ approximateNumericArrayQ[m],
        Abs[1 - p] <= Length[m] 10 ^ (2 - Precision[p]),
        TrueQ[Simplify[p] == 1]
    ]
]

densityMatrixAmplitudes[m_] := With[{k = First @ Ordering[Re @ Normal @ Diagonal[m], -1]},
    Normal[m[[All, k]]] / Sqrt[Re @ m[[k, k]]]
]

QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, "Concurrence"] := (
    warnUnphysicalMonotone[qs];
    Enclose @ With[{bp = ConfirmBy[qs["Bipartition", biPartition]["Normalized"], QuantumStateQ[#] && #["Qudits"] == 2 &]},
        Which[
            bp["VectorQ"] && approximateNumericArrayQ[bp["StateVector"]],
            schmidtConcurrence[Normal @ bp["Computational"]["StateVector"], bp["Dimensions"]],
            bp["MatrixQ"] && approximateNumericArrayQ[bp["DensityMatrix"]],
            If[ bp["Dimension"] > 64 && pureDensityMatrixQ[bp["DensityMatrix"]],
                schmidtConcurrence[densityMatrixAmplitudes[bp["Computational"]["DensityMatrix"]], bp["Dimensions"]],
                Norm @ concurrenceVector[bp["Operator"]]
            ],
            bp["VectorQ"] || pureDensityMatrixQ[bp["DensityMatrix"]],
            If[ bp["VectorQ"] && Normal[bp["StateVector"]] === ConstantArray[0, bp["Dimension"]],
                0,
                With[{val = 2 (1 - (QuantumPartialTrace[bp, {1}] ^ 2)["Norm"])},
                    Sqrt[If[NumericQ[val], Max[0, Re[val]], val]]
                ]
            ],
            True,
            Norm @ concurrenceVector[bp["Operator"]]
        ]
    ]
)


QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, "Negativity"] := (
    warnUnphysicalMonotone[qs];
    Enclose[(ConfirmBy[qs["Bipartition", biPartition]["Normalized"], QuantumStateQ[#] && #["Qudits"] == 2 &]["Transpose", {2}]["TraceNorm"] - 1) / 2]
)


QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, "LogNegativity"] := (
    warnUnphysicalMonotone[qs];
    Enclose @ Log2 @ ConfirmBy[qs["Bipartition", biPartition]["Normalized"], QuantumStateQ[#] && #["Qudits"] == 2 &]["Transpose", {2}]["TraceNorm"]
)


(* Entanglement entropy is the von Neumann entropy of a reduced state, an entanglement measure only for
   a pure global state. A pure state reads it through its Schmidt weights (vector) or its reduced state
   (density-matrix form); a genuinely mixed state's reduced von Neumann entropy mixes classical ignorance
   into the count, so it is not returned as entanglement. *)
QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, "EntanglementEntropy"] := (
    warnUnphysicalMonotone[qs];
    Enclose @ With[{
        bp = ConfirmBy[qs["Bipartition", biPartition]["Normalized"], QuantumStateQ[#] && #["Qudits"] == 2 &]
    },
        Which[
            bp["VectorQ"],
                Quantity[Total[-# Log2[#] & @ Select[Confirm @ bp["SchmidtBasis"]["Probability"], If[NumericQ[#], # > 0, True] &]], "Bits"],
            TrueQ[bp["PureStateQ"]],
                QuantumPartialTrace[bp, {1}]["VonNeumannEntropy"],
            True,
                Message[QuantumEntanglementMonotone::mixedentropy]; Indeterminate
        ]
    ]
)

QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, "RenyiEntanglementEntropy" | "RenyiEntropy"] :=
    QuantumEntanglementMonotone[qs, biPartition, {"RenyiEntanglementEntropy", 1 / 2}]

(* The Renyi entanglement entropy is the Renyi entropy of the reduced state, an entanglement measure only
   for a pure global state, so a genuinely mixed input is guarded exactly as EntanglementEntropy is. *)
QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, {"RenyiEntanglementEntropy" | "RenyiEntropy", alpha_}] := (
    warnUnphysicalMonotone[qs];
    Enclose @ With[{
        bp = ConfirmBy[qs["Bipartition", biPartition]["Normalized"], QuantumStateQ[#] && #["Qudits"] == 2 &]
    },
        If[ bp["VectorQ"] || TrueQ[bp["PureStateQ"]],
            With[{val = (1 / (1 - alpha)) Log[2, Tr @ MatrixPower[QuantumPartialTrace[bp, {1}]["DensityMatrix"], alpha]]},
                If[NumericQ[val], Re[val], val]
            ],
            Message[QuantumEntanglementMonotone::mixedentropy]; Indeterminate
        ]
    ]
)

QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, "Realignment"] := (
    warnUnphysicalMonotone[qs];
    Enclose @ realignmentValue @ ConfirmBy[qs["Bipartition", biPartition]["Normalized"], QuantumStateQ[#] && #["Qudits"] == 2 &]
)

realignmentValue[bqs_] := Total @ SingularValueList @ ArrayReshape[Transpose[bqs["Bend"]["Tensor"], 2 <-> 3], bqs["Dimensions"] ^ 2] - 1


(* Quantum Discord *)

MutualInformationI[rho_QuantumState, biPartition_ : Automatic] := Enclose @ With[
    {s = ConfirmBy[rho["Bipartition", biPartition], QuantumStateQ[#] && #["Qudits"] == 2 &]},

	QuantumPartialTrace[s, {1}]["Entropy"] + QuantumPartialTrace[s, {2}]["Entropy"] - s["Entropy"]
]

MutualInformationJ[rho_QuantumState, qm : _QuantumMeasurementOperator | Automatic : Automatic,biPartition_ : Automatic] := Enclose @ With[
    {s = ConfirmBy[rho["Bipartition", biPartition], QuantumStateQ[#] && #["Qudits"] == 2 &]},
    {m = QuantumMeasurementOperator[Replace[qm, Automatic :> QuantumMeasurementOperator[]], {2}][s]}, 

    QuantumPartialTrace[s, {2}]["Entropy"] - m["ProbabilitiesList"] . (Simplify[QuantumPartialTrace[#, {2}]]["Entropy"] & /@ m["States"])
]

QuantumDiscord[rho_QuantumState, qm : _QuantumMeasurementOperator | Automatic : Automatic, biPartition_ : Automatic] :=
	Enclose[Confirm[MutualInformationI[rho, biPartition]] - Confirm[MutualInformationJ[rho, qm, biPartition]]]

QuantumEntanglementMonotone[qs_ ? QuantumStateQ, biPartition_ : Automatic, "MutualInformationI"] :=
    (warnUnphysicalMonotone[qs]; MutualInformationI[qs, biPartition])

QuantumEntanglementMonotone[qs_ ? QuantumStateQ, qm : _QuantumMeasurementOperator | Automatic : Automatic, biPartition_ : Automatic, "MutualInformationJ"] :=
    (warnUnphysicalMonotone[qs]; MutualInformationJ[qs, qm, biPartition])

QuantumEntanglementMonotone[qs_ ? QuantumStateQ, qm : _QuantumMeasurementOperator | Automatic : Automatic, biPartition_ : Automatic, "Discord"] :=
    (warnUnphysicalMonotone[qs]; QuantumDiscord[qs, qm, biPartition])

