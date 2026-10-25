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

QuantumEntanglementMonotone::realignment =
    "The realignment value ||R(rho)||_1 - 1 is the statistic of a separability test, not an entanglement monotone: a positive value proves the state entangled, but the value is negative on mixed product states, and local operations can raise it."

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

QuantumEntangledQ::args = "The arguments `1` after the state are not a split followed by a method."

QuantumEntanglementMonotone::args =
    "The arguments `1` after the state are not a split followed by a method; MutualInformationJ and Discord also take a QuantumMeasurementOperator before the split."

QuantumEntangledQ::method = QuantumEntanglementMonotone::method = "`1` is not a method; the methods are Automatic and `2`."

QuantumEntangledQ::correlation =
    "`1` measures all correlations, classical ones included, and is positive on some separable states, so it is no test of entanglement."

QuantumEntangledQ::split = QuantumEntanglementMonotone::split =
    "`1` is not a split of a state of `2` qudits: give Automatic, a list of distinct qudits that leaves at least one out, or two disjoint nonempty lists of qudits."

QuantumEntangledQ::quditsplit = QuantumEntanglementMonotone::quditsplit =
    "`1` is not a split of a single qudit of dimension `2`: give Automatic or the dimension of the first part, one of `3`."

QuantumEntangledQ::nosplit = QuantumEntanglementMonotone::nosplit =
    "A state of dimension `1` on a single qudit has no split into two parts of dimension greater than 1."

QuantumEntangledQ::dependent = QuantumEntanglementMonotone::dependent =
    "The elements of the state's basis do not form a basis of each qudit's space: they are linearly dependent, too few, or too close to dependent to invert at the precision of their entries."

QuantumEntangledQ::order = QuantumEntanglementMonotone::order =
    "`1` is not a Renyi order: give a finite nonnegative number, or a real expression in parameters that is not negative for all of them."

QuantumEntangledQ::nodensitymatrix = QuantumEntanglementMonotone::nodensitymatrix =
    "The state is held in the phase-space picture, which has no density matrix in the computational basis."

QuantumEntangledQ::state = QuantumEntanglementMonotone::state = "The first argument must be a QuantumState; it is `1`."

QuantumEntangledQ::zero = QuantumEntanglementMonotone::zero =
    "The state has no normalized form: it is a state vector of norm 0, or a density matrix whose trace is 0, for approximate entries within the round-off of its diagonal."


(* The methods each function takes besides Automatic: QuantumEntanglementMonotone every monotone and the concurrence
   vector, QuantumEntangledQ the tests of entanglement, which the correlation measures are not. "RenyiEntanglementEntropy"
   is another name for "RenyiEntropy", and {"RenyiEntropy", alpha} gives its order. *)
$correlationMeasures = {"MutualInformationI", "MutualInformationJ", "Discord"};

entanglementMethods[QuantumEntanglementMonotone] := Append[$QuantumEntanglementMonotones, "ConcurrenceVector"]
entanglementMethods[QuantumEntangledQ] := DeleteCases[$QuantumEntanglementMonotones, Alternatives @@ $correlationMeasures]

methodName["RenyiEntanglementEntropy" | {"RenyiEntanglementEntropy" | "RenyiEntropy", _}] := "RenyiEntropy"
methodName[name_String] := name
methodName[_] := None

(* A Renyi order is a nonnegative number, or an expression in parameters of the user's (symbolicParameters), none of
   them a system symbol such as $Failed, None or Missing, that is not negative, or nonreal, for all real values of
   them: x + I is refused, 2 + I x, real at x = 0, is read. *)
renyiOrderQ[alpha_ ? NumericQ] := TrueQ[NonNegative[alpha]]
renyiOrderQ[alpha_] := FreeQ[alpha, _String | _List] && With[{ps = symbolicParameters[alpha]},
    ps =!= {} && AllTrue[ps, FreeQ[#, s_Symbol /; Context[s] === "System`" && s =!= Subscript] &] &&
        ! TrueQ[TimeConstrained[Simplify[alpha < 0 || Im[alpha] != 0, Element[ps, Reals]], 1]]
]

(* The arguments after the state, read as {split, method}, or as {measurement, split, method} for MutualInformationJ
   and Discord, which measure the second part. The method is always last: a string, or a list headed by one, is a
   method, and a split never is. The method, the split and the state are then checked, in that order, and the first
   problem gives one message and a Failure. *)
entanglementArguments[head_, qs_, args_List] := With[{
    read = Replace[args, {
        {qm : _QuantumMeasurementOperator | Automatic : Automatic, split_ : Automatic, method : "MutualInformationJ" | "Discord"} :>
            {qm, split, method},
        {_QuantumMeasurementOperator, ___} -> $Failed,
        {split : Except[_String | {_String, ___}] : Automatic, method_ : Automatic} :> {split, method},
        {split_, method_} :> {split, method},
        _ -> $Failed
    }]
},
    With[{method = If[read === $Failed, None, Last[read]], split = If[read === $Failed, None, read[[-2]]]},
        Which[
            read === $Failed,
            entanglementFailure[head, "InvalidArguments", "args", {args}],
            head === QuantumEntangledQ && MemberQ[$correlationMeasures, methodName[method]],
            entanglementFailure[head, "InvalidMethod", "correlation", {method}],
            MatchQ[method, {"RenyiEntanglementEntropy" | "RenyiEntropy", alpha_} /; ! renyiOrderQ[alpha]],
            entanglementFailure[head, "InvalidMethod", "order", {Last[method]}],
            method =!= Automatic && ! MemberQ[entanglementMethods[head], methodName[method]],
            entanglementFailure[head, "InvalidMethod", "method", {method, entanglementMethods[head]}],
            ! TrueQ[splitQ[qs, split]],
            splitFailure[head, qs, split],
            qs["Picture"] === "PhaseSpace",
            entanglementFailure[head, "NoDensityMatrix", "nodensitymatrix", {}],
            ! TrueQ[independentFrameQ[qs]],
            entanglementFailure[head, "DependentBasis", "dependent", {}],
            TrueQ[zeroStateQ[qs]],
            entanglementFailure[head, "ZeroState", "zero", {}],
            True,
            read
        ]
    ]
]

entanglementFailure[head_, tag_, name_, params_] := (
    Message[MessageName[head, name], Sequence @@ params];
    Failure[tag, <|"MessageTemplate" :> MessageName[head, name], "MessageParameters" -> params|>]
)

(* A split: Automatic, the first qudit against the rest; a list of distinct qudits against all the others, of which
   there is at least one; two disjoint nonempty lists of qudits, every other qudit traced out; or, for a single qudit
   of dimension D, the dimension d of the first part, a divisor of D with 1 < d < D, which Automatic takes too. *)
splitQ[qs_, Automatic] := qs["Qudits"] > 1 || quditSplits[qs] =!= {}
splitQ[qs_, d_Integer] := MemberQ[quditSplits[qs], d]
splitQ[qs_, qudits : {__Integer}] := Length[qudits] < qs["Qudits"] && DuplicateFreeQ[qudits] && SubsetQ[Range[qs["Qudits"]], qudits]
splitQ[qs_, {a : {__Integer}, b : {__Integer}}] := DuplicateFreeQ[Join[a, b]] && SubsetQ[Range[qs["Qudits"]], Join[a, b]]
splitQ[_, _] := False

quditSplits[qs_] := If[qs["Qudits"] == 1, Select[Divisors[qs["Dimension"]], 1 < # < qs["Dimension"] &], {}]

splitFailure[head_, qs_, split_] := Which[
    qs["Qudits"] > 1, entanglementFailure[head, "InvalidSplit", "split", {split, qs["Qudits"]}],
    quditSplits[qs] === {}, entanglementFailure[head, "InvalidSplit", "nosplit", {qs["Dimension"]}],
    True, entanglementFailure[head, "InvalidSplit", "quditsplit", {split, qs["Dimension"], quditSplits[qs]}]
]

(* Whether the elements of the state's basis form a basis of each qudit's space, a square matrix of full rank
   (invertibleFrameQ), as in the computational frame, and, for an approximate state, one well enough conditioned that
   the change to the computational frame, which multiplies its round-off by kappa for a state vector and kappa^2 for a
   density matrix (frameContext), leaves it a digit. A zero state is left to zeroStateQ. *)
independentFrameQ[qs_] := qs["ComputationalQ"] || invertibleFrameQ[qs["Basis"]] && ! roundoffSwampedQ[qs]

roundoffSwampedQ[qs_] := qs["NumericQ"] && Precision[qs["State"]] < Infinity && With[{m = Normal @ qs["State"]},
    Max[Abs[m]] > 0 && TrueQ[frameCondition[qs["Basis"]] ^ If[qs["VectorQ"], 1, 2] qs["Dimension"] 10 ^ -wholeDigits[m] >= 1]
]

(* A state with no normalized form: a state vector whose coefficients are all 0, or a density matrix of trace 0. Both
   are read in the frame the state is stored in, which is invertible by now, and for approximate entries conditioned
   well enough to keep a digit, so they are 0 exactly when the state is.
   An exact value is 0 when its sign is; an approximate vector must be 0, and an approximate trace lie within the
   round-off of the diagonal it sums, read relative to its largest entry s: the trace divided by s lies between 0 and
   n, and falls below the range of machine numbers only when it is 0 relative to s anyway. A symbolic density matrix
   is 0 when its trace simplifies to 0; a symbolic state vector, whose coefficients would evaluate to 0, is not
   tested. A machine vector all of whose coefficients lie below the normal range of machine numbers has no norm that
   can be divided out either. The state is read as stored, a sparse array kept sparse. *)
zeroStateQ[qs_] := With[{a = qs["State"]},
    Which[
        ! qs["NumericQ"], qs["MatrixQ"] && TimeConstrained[Simplify[Tr[a]], $entanglementSignTime / 5, Tr[a]] === 0,
        Precision[a] === Infinity, quantitySign[If[qs["VectorQ"], Conjugate[a] . a, Abs[Tr[a]]]] === 0,
        qs["VectorQ"], Max[Abs[a]] == 0 || Precision[a] === MachinePrecision && Max[Abs[a]] < $MinMachineNumber,
        True, With[{s = Max[Abs[Diagonal[a]]]}, s == 0 || Abs[Tr[a]] / s <= criterionCut[a]]
    ]
]


(* QuantumEntangledQ answers True when the method's test proves the two-part state entangled, False when a
   separability proof holds for it, and Indeterminate when neither does. Automatic tests the partial transpose in
   every size, for a negative eigenvalue, and beyond 2 x 3, where no efficient test is complete, the realignment
   value ||R(rho)||_1 - 1 as well, which detects some states with a positive partial transpose. Its separability
   proofs: a positive partial transpose where that test is exact (pptSeparableQ), and a purity Tr[rho^2] of at most
   1 / (dA dB - 1), the largest ball of states around the maximally mixed one, all of them separable (Gurvits and
   Barnum, PRA 66, 062311 (2002)). A named method answers with its own test alone: Negativity and LogNegativity
   with the partial transpose and the proofs it gives; the entropies and the concurrence of a pure two-part state with
   its Schmidt rank; the concurrence of an arbitrary-precision mixed state beyond two qubits with its machine-number
   value; and any other method with a positive value of its monotone (the realignment value within the time limit
   that exact and symbolic input get), or a zero value where that method is exact. A symbolic input raises a message naming its parameters, taken to be real; a symbolic state whose
   physical parameter values cannot be found is Indeterminate. Invalid arguments and states, and the mutual
   informations and the discord, which are positive on separable states too, are refused before any of this, by
   entanglementArguments. *)
QuantumEntangledQ[qs_ ? QuantumStateQ, args___] := With[{a = entanglementArguments[QuantumEntangledQ, qs, {args}]},
    If[FailureQ[a], a, entangledQ[qs, Sequence @@ a]]
]

QuantumEntangledQ[args___] /; ! MatchQ[{args}, {_ ? QuantumStateQ, ___}] :=
    entanglementFailure[QuantumEntangledQ, "InvalidState", "state", {If[{args} === {}, "missing", First[{args}]]}]

entangledQ[qs_, biPartition_, method_] :=
    Enclose @ With[{q = exactSymbolicState[qs]},
        If[! q["NumericQ"], Message[QuantumEntangledQ::realparameters, DeleteDuplicates @ symbolicParameters[Normal @ q["State"]]]];
        With[{bp = ConfirmBy[q["Bipartition", biPartition]["Normalized"], QuantumStateQ[#] && #["Qudits"] == 2 &]},
            With[{c = frameContext[numberContext[bp], q]},
                If[c["Kind"] === "Symbolic" && ! c["Exists"], Indeterminate, entangledAnswer[method, q, biPartition, bp, c]]
            ]
        ]
    ]

(* An approximate state stored in a frame that is not unitary reaches the computational frame through the frame and its
   inverse, which multiply its round-off by the frame's condition number kappa, the product over the qudits of the
   ratio of the largest to the smallest singular value of their elements: kappa for a state vector, and kappa^2 for a
   density matrix, changed by similarity. Its cut is scaled by that factor, so that amplified round-off proves no sign,
   and no value is read as zero (a floor of -1): amplified round-off carries a criterion into the floor as readily as
   out of it. A frame whose kappa is within 10^-6 of 1, unitary up to a scale, changes nothing. *)
frameContext[c_, qs_] := If[
    c["Kind"] === "Approximate" && ! qs["ComputationalQ"],
    With[{k = frameCondition[qs["Basis"]] ^ If[qs["VectorQ"], 1, 2]},
        If[k < 1 + 10 ^ -6, c, <|c, "Cut" -> k c["Cut"], "Floor" -> -1|>]
    ],
    c
]

(* An exact frame with a nonzero entry outside 10^-290 to 10^290 in size (machineRangeQ) is read with 20 digits; a
   frame that is not numeric gives 1. *)
frameCondition[qb_] := Times @@ Map[
    With[{e = Normal @ #["ReducedMatrix"]},
        If[ ArrayQ[e, 2, NumericQ],
            With[{m = If[machineRangeQ[e], N[e], N[e, 20]]},
                With[{s = SingularValueList[m, Tolerance -> 0]}, If[Length[s] < Length[m] || Min[s] == 0, Infinity, Max[s] / Min[s]]]
            ],
            1
        ]
    ] &,
    Join[qb["Output"]["Decompose"], qb["Input"]["Decompose"]]
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

(* The realignment value of a pure state, (Sum[s_i])^2 - 1, is twice its Schmidt negativity, already read; it is
   computed for density matrices only. *)
realignmentDetectsQ[bp_, c_] := bp["MatrixQ"] && Times @@ bp["Dimensions"] > 6 && criterionSign[realignmentCriterion[bp, c], c] === 1

entangledAnswer["Negativity" | "LogNegativity", qs_, _, bp_, c_] := (
    warnUnphysicalMonotone[qs];
    With[{npt = partialTransposeNegativeQ[bp, c]},
        Which[TrueQ[npt], True, npt === False && pptSeparableQ[bp, c], False, True, Indeterminate]
    ]
)

(* The realignment value is read with the time limit realignmentCriterion gives exact and symbolic input; zero proves
   separability on a pure state only. *)
entangledAnswer["Realignment", qs_, _, bp_, c_] := (
    warnUnphysicalMonotone[qs];
    With[{sign = criterionSign[realignmentCriterion[bp, c], c]},
        Which[sign === 1, True, sign === 0 && pureTwoPartQ[bp, c], False, True, Indeterminate]
    ]
)

(* An entropy or the concurrence of a pure two-part state is positive exactly when its Schmidt rank exceeds 1, which
   the partial transpose reads: exactly on exact input, for every physical parameter value on symbolic input, and to
   first order in the smaller Schmidt coefficients on approximate input. The entropies of order 1 and above are of
   second order in them, below the cut of the state's numbers long before they are zero: the von Neumann entropy of
   the ket |00> + 10^-11 |11> is about 10^-20. A Renyi entropy of an approximate order comes out approximate on exact
   input too. *)
entangledAnswer[method_, qs_, _, bp_, c_] /;
        MatchQ[methodName[method], "EntanglementEntropy" | "RenyiEntropy" | "Concurrence"] && TrueQ[pureTwoPartQ[bp, c]] := (
    warnUnphysicalMonotone[qs];
    Replace[schmidtRankAboveOneQ[bp, c], Except[True | False] -> Indeterminate]
)

(* Whether a pure two-part state has Schmidt rank above 1. An exact density matrix is read through the single scalar
   Tr[rho]^2 - Tr[rhoA^2], as an exact state vector is, rather than the eigenvalues of its partial transpose. *)
schmidtRankAboveOneQ[bp_, c_] /; bp["MatrixQ"] && c["Kind"] === "Exact" := With[{m = Normal @ bp["DensityMatrix"], d = bp["Dimensions"]},
    With[{a = TensorContract[ArrayReshape[m, Join[d, d]], {{2, 4}}]},
        Replace[quantitySign[Tr[m] ^ 2 - Total[a Transpose[a], 2]], {1 -> True, 0 -> False, _ -> Indeterminate}]
    ]
]
schmidtRankAboveOneQ[bp_, c_] := partialTransposeNegativeQ[bp, c]

(* Beyond two qubits, the mixed-state concurrence of an arbitrary-precision density matrix keeps an error of order
   10^-18 near the separable states, whatever its precision, while claiming many more digits: isotropic states at their
   separability boundary read as entangled. So it is read from the matrix in machine numbers, against the larger of
   their cut and the input's own, and a value within it proves nothing. The matrix is read raw: a QuantumState of it would equal the
   input, whose cached properties keep the input's precision. The Wootters value of two qubits keeps its digits. *)
entangledAnswer["Concurrence", qs_, _, bp_, c_] /;
        c["Kind"] === "Approximate" && bp["MatrixQ"] && Times @@ bp["Dimensions"] > 4 && Precision[bp["DensityMatrix"]] =!= MachinePrecision &&
            ! pureTwoPartQ[bp, c] := (
    warnUnphysicalMonotone[qs];
    With[{m = N[Normal @ bp["DensityMatrix"]]},
        If[ roundoffSign[Norm @ concurrenceMatrixVector[m, bp["Dimensions"]], With[{mc = approximateContext[m]}, <|mc, "Cut" -> Max[mc["Cut"], c["Cut"]]|>]] === 1,
            True,
            Indeterminate
        ]
    ]
)

(* A zero value proves separability where the method is exact: the concurrence on two qubits and on pure states,
   the realignment value and the entropies on pure states. An entropy is read in bits. *)
entangledAnswer[method_, qs_, biPartition_, bp_, c_] := With[{
    sign = criterionSign[Replace[entanglementMonotone[qs, biPartition, method], Quantity[x_, "Bits"] :> x], c]
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
   eigenvalue is negative, for every physical value of the parameters. The approximate negativity sums every singular
   value of the partial transpose of the Hermitian part of rho, small ones included, which the default tolerance would
   drop: round-off that breaks Hermiticity would add singular values of its own size. It proves a negative
   eigenvalue above the cut of numberContext, and none at most its floor, and between the two decides nothing. A state
   vector's is read off its min(dA, dB) Schmidt coefficients as Sum[s_i s_j, i < j], with no cancellation, and their
   round-off grows far more slowly than n, so its cut and floor take min(dA, dB) Sqrt[max(dA, dB)] in place of n:
   product kets in 2 x 4096 reach 4 machine epsilons. *)
partialTransposeNegativeQ[bp_, c_] /; bp["VectorQ"] && c["Kind"] =!= "Symbolic" :=
    With[{m = ArrayReshape[Normal @ bp["StateVector"], bp["Dimensions"]]},
        If[ c["Kind"] === "Approximate",
            With[{s = SingularValueList[m, Tolerance -> 0]}, roundoffNegativeQ[Rest[s] . Most[Accumulate[s]], c, 1 / Sqrt[Max[Dimensions[m]]]]],
            Replace[quantitySign[With[{a = m . ConjugateTranspose[m]}, Tr[a] ^ 2 - Total[a Transpose[a], 2]]], {1 -> True, 0 -> False, _ -> Indeterminate}]
        ]
    ]

partialTransposeNegativeQ[bp_, c_] := With[{pt = partialTransposeMatrix[Normal @ bp["DensityMatrix"], bp["Dimensions"]]},
    Switch[c["Kind"],
        "Approximate",
        roundoffNegativeQ[(Total[SingularValueList[hermitianPart[pt], Tolerance -> 0]] - 1) / 2, c, 1],
        "Exact",
        With[{signs = quantitySign /@ Eigenvalues[pt]}, Which[MemberQ[signs, -1], True, SubsetQ[{0, 1}, signs], False, True, Indeterminate]],
        "Symbolic",
        With[{e = positiveSemidefiniteCoefficients[pt, c["Parameters"]]},
            Which[! ListQ[e], Indeterminate, forAllQ[c, Or @@ Thread[e < 0]], True, forAllQ[c, And @@ Thread[e >= 0]], False, True, Indeterminate]
        ]
    ]
]

roundoffNegativeQ[negativity_, c_, scale_] := Replace[roundoffSign[negativity / scale, c], {1 -> True, 0 | -1 -> False, _ -> Indeterminate}]

hermitianPart[m_] := (m + ConjugateTranspose[m]) / 2

partialTransposeMatrix[m_, {dA_, dB_}] := ArrayReshape[Transpose[ArrayReshape[m, {dA, dB, dA, dB}], {1, 4, 3, 2}], {dA dB, dA dB}]

(* A positive partial transpose proves separability where the test is exact: on 2 x 2 and 2 x 3 (Horodecki,
   Horodecki and Horodecki, PLA 223, 1 (1996)), on a pure state, and on a density matrix whose rank is at most the
   larger of the ranks of its two marginals (Horodecki, Lewenstein, Vidal and Cirac, PRA 62, 032310 (2000), Thm 1;
   the rank against the local dimensions alone proves nothing, since a state can sit in a larger space). Approximate
   ranks are counted so that an eigenvalue of uncertain size never makes the proof: the rank of rho counts every
   eigenvalue above round-off, tol = n 10^(2 - p) (roundoffCut) with p the digits its entries are known to, and the
   ranks of the marginals only those above 10^-10. The eigenvalues rho drops lift an eigenvalue of a marginal by at
   most n tol, so the proof is given only where n tol <= 10^-10, which keeps them from adding to a marginal's rank: up
   to n = 100 at machine precision. An approximate density matrix is read as pure when its impurity 1 - Tr[rho^2] is
   at most tol, and only where tol <= 10^-10. The digits p are rounded up to a whole one, so that the fraction of a
   digit a normalization costs never moves the cut, and a higher p only makes either proof harder to give. A symbolic
   density matrix is not ranked. *)
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

roundoffCut[m_] := Length[m] 10 ^ (2 - wholeDigits[m])

(* The cut below which an approximate criterion of an n-dimensional state counts as zero: roundoffCut's n 10^(2 - p),
   about 45 n machine epsilons, for machine numbers, whose arithmetic adds round-off of a few n machine epsilons of
   its own, and n 10^-p above machine precision, where the arithmetic tracks its own error and the partial transpose
   and the realignment of separable states stay far below it. The mixed-state concurrence does not, and is read in
   machine numbers. *)
criterionCut[m_] := If[Precision[m] === MachinePrecision, roundoffCut[m], Length[m] 10 ^ -wholeDigits[m]]

(* The context of an approximate array m: its cut, and the floor within which a value is zero, n machine epsilons for
   machine numbers and the cut above machine precision. *)
approximateContext[m_] := With[{cut = criterionCut[m]},
    <|"Kind" -> "Approximate", "Cut" -> cut, "Floor" -> If[Precision[m] === MachinePrecision, Length[m] $MachineEpsilon, cut]|>
]

(* The digits to which the entries of an approximate array m are known, relative to its largest entry. They are read
   from its accuracy: an entry that cancelled to zero carries accuracy but no precision, and would make Precision[m]
   0. The cuts round them up to a whole digit, past the float noise of 10^-6 that this sum can leave on a whole
   number. *)
arrayDigits[m_] := Accuracy[m] + Log10[Max[Abs[m]]]

wholeDigits[m_] := Ceiling[arrayDigits[m] - 10 ^ -6]

purityBoundQ[bp_, c_] := bp["MatrixQ"] && With[{m = Normal @ bp["DensityMatrix"]},
    nonnegativeQ[1 / (Times @@ bp["Dimensions"] - 1) - Total[m Transpose[m], 2], c]
]

(* The realignment value of the two-part state, given a time limit on exact and symbolic input, whose singular values
   may not come out in closed form, and read from the Hermitian part of an approximate density matrix. *)
realignmentCriterion[bp_, c_] := Which[
    c["Kind"] =!= "Approximate", TimeConstrained[realignmentValue[bp], $entanglementSignTime, Indeterminate],
    bp["MatrixQ"], With[{d = bp["Dimensions"]},
        Total[SingularValueList[ArrayReshape[Transpose[ArrayReshape[hermitianPart[Normal @ bp["DensityMatrix"]], Join[d, d]], {1, 3, 2, 4}], d ^ 2], Tolerance -> 0]] - 1
    ],
    True, realignmentValue[bp]
]


(* The number type of the two-part state bp, with what deciding a sign on it needs: "Approximate" with the round-off
   cut of its numbers, "Exact", or "Symbolic" with the parameters, taken as real, the condition that bp is a density
   matrix for their values: Hermitian, and every coefficient of det(x + rho) nonnegative, so that no eigenvalue is
   negative, and whether any value meets it. The condition is $Failed when it is not found within
   $entanglementSignTime seconds. The cut is criterionCut's: each eigenvalue or singular value a criterion sums
   carries round-off of order 10^-p at precision p, and an n-dimensional criterion sums at most n of them. A value
   proves a sign beyond the cut, and is zero within the floor below it: n machine epsilons for machine numbers, whose
   separable states keep their criteria at about a third of that at most, and the cut itself above machine
   precision. *)
numberContext[bp_] := With[{a = If[bp["VectorQ"], bp["StateVector"], bp["DensityMatrix"]]},
    Which[
        ! ArrayQ[a, _, NumericQ], symbolicContext[Normal @ bp["DensityMatrix"]],
        Precision[a] === Infinity, <|"Kind" -> "Exact"|>,
        True, approximateContext[Normal @ a]
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
   value has a sign beyond the round-off cut of the state's numbers and is zero within its floor, both from
   numberContext, and is Indeterminate between the two; one read off an exact state is Indeterminate. An exact value is read in 60-digit
   significance arithmetic, by SetPrecision, which raises no N::meprec: N does, with a precision goal on an exact
   zero, and with an accuracy goal on a sum of large terms that cancel, which it can then read as zero. The reading
   gives the sign when any of its digits survives the cancellation; a value that keeps none is decided within
   $entanglementSignTime seconds: by readings at 200 and 1000 digits; then as zero when FullSimplify reduces x to 0
   within a fifth of that time; then by a 5000-digit reading, which for a constant such as Khinchin's takes longer
   than the whole limit; then, when x is algebraic, by the sign of RootReduce[x], 0 for an exact zero. A number built
   from rationals and Root objects by sums, products, rational powers and Abs is algebraic, as is any that Element
   proves algebraic. Any other is Indeterminate, with a message, which shows x unevaluated, each integer or rational
   of more than 20 digits in it replaced by its 6-digit value. The quantity of a symbolic state has a sign when it has
   that sign for every physical value of the parameters. *)
$entanglementSignTime = 5;

criterionSign[v_, c_] := Switch[c["Kind"],
    "Symbolic", symbolicSign[realForm[v, c], c],
    "Approximate", roundoffSign[v, c],
    _, quantitySign[v]
]

(* The value is compared by its nominal value, the binary number it holds: the comparison of a low-precision number
   with the cut would otherwise be tolerant, and read a value above it as within it. *)
roundoffSign[x_ ? NumericQ, c_] := With[{r = SetPrecision[Re[x], Infinity]}, Which[Abs[r] <= c["Floor"], 0, Abs[r] > c["Cut"], Sign[r], True, Indeterminate]]
roundoffSign[_, _] := Indeterminate

quantitySign[x_ ? NumericQ] /; Precision[x] === Infinity := With[{v = Re[SetPrecision[x, 60]]}, If[Precision[v] > 0, Sign[v], exactSign[x]]]
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

(* A margin that lies between the floor and the cut, either side of 0, decides nothing. *)
nonnegativeQ[v_, c_] := If[c["Kind"] === "Symbolic", forAllQ[c, realForm[v, c] >= 0], MatchQ[criterionSign[v, c], 0 | 1]]



(* The arguments are read and checked once, by entanglementArguments, and each method below receives a valid split.
   The default method is decided on the two-part state, after the split and the partial trace: the concurrence for
   a pure state and for two qubits, where it is exact, and the negativity for a mixed state beyond two qubits, where
   the concurrence is only a lower bound that a local unitary can change. *)
QuantumEntanglementMonotone[qs_ ? QuantumStateQ, args___] := With[{a = entanglementArguments[QuantumEntanglementMonotone, qs, {args}]},
    If[FailureQ[a], a, entanglementMonotone[qs, Sequence @@ a]]
]

QuantumEntanglementMonotone[args___] /; ! MatchQ[{args}, {_ ? QuantumStateQ, ___}] :=
    entanglementFailure[QuantumEntanglementMonotone, "InvalidState", "state", {If[{args} === {}, "missing", First[{args}]]}]

entanglementMonotone[qs_, biPartition_, Automatic] :=
    Enclose @ entanglementMonotone[qs, biPartition, defaultMonotone[ConfirmBy[qs["Bipartition", biPartition], QuantumStateQ]]]

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
concurrenceVector[rho_] := concurrenceMatrixVector[Normal @ rho["Matrix"], rho["OutputDimensions"]]

(* The concurrence vector of a normalized two-qudit density matrix mat of dimensions {d1, d2}. *)
concurrenceMatrixVector[mat_, {d1_, d2_}] := With[{
	component = If[ approximateNumericArrayQ[mat],
		With[{w = densityMatrixFactor[mat]},
			Function[{ya, yb}, wootterCombination @ SingularValueList[Transpose[w] . KroneckerProduct[ya, yb] . w, Tolerance -> 0]]
		],
		With[{cmat = Conjugate[mat]},
			Function[{ya, yb}, With[{omat = Normal @ KroneckerProduct[ya, yb]},
				wootterCombination @ Sqrt[If[NumericQ[#], Max[0, #], #] & /@ Re[Eigenvalues[mat . omat . cmat . omat]]]
			]]
		]
	]
},
	Flatten[Outer[component, Y[d1], Y[d2], 1], 1]
]

entanglementMonotone[qs_, biPartition_, "ConcurrenceVector"] :=
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

entanglementMonotone[qs_, biPartition_, "Concurrence"] := (
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


entanglementMonotone[qs_, biPartition_, "Negativity"] := (
    warnUnphysicalMonotone[qs];
    Enclose[(ConfirmBy[qs["Bipartition", biPartition]["Normalized"], QuantumStateQ[#] && #["Qudits"] == 2 &]["Transpose", {2}]["TraceNorm"] - 1) / 2]
)


entanglementMonotone[qs_, biPartition_, "LogNegativity"] := (
    warnUnphysicalMonotone[qs];
    Enclose @ Log2 @ ConfirmBy[qs["Bipartition", biPartition]["Normalized"], QuantumStateQ[#] && #["Qudits"] == 2 &]["Transpose", {2}]["TraceNorm"]
)


(* The entanglement entropy and the Renyi entanglement entropies are entropies of the reduced state, entanglement
   measures only for a pure two-part state: the reduced entropy of a mixed state mixes classical ignorance into the
   count, so it is not returned as entanglement. They are read off the Schmidt weights, in bits, and the entanglement
   entropy, the von Neumann entropy, is the Renyi entropy of order 1. *)
entanglementMonotone[qs_, biPartition_, "EntanglementEntropy"] := entanglementMonotone[qs, biPartition, {"RenyiEntropy", 1}]

entanglementMonotone[qs_, biPartition_, "RenyiEntanglementEntropy" | "RenyiEntropy"] :=
    entanglementMonotone[qs, biPartition, {"RenyiEntanglementEntropy", 1 / 2}]

entanglementMonotone[qs_, biPartition_, {"RenyiEntanglementEntropy" | "RenyiEntropy", alpha_}] := (
    warnUnphysicalMonotone[qs];
    Enclose @ With[{
        bp = ConfirmBy[qs["Bipartition", biPartition]["Normalized"], QuantumStateQ[#] && #["Qudits"] == 2 &]
    },
        With[{impurity = pureImpurity[qs, bp]},
            If[ NumericQ[impurity],
                Quantity[pureStateRenyiEntropy[qs, bp, alpha, impurity], "Bits"],
                Message[QuantumEntanglementMonotone::mixedentropy]; Indeterminate
            ]
        ]
    ]
)

(* The impurity 1 - Tr[rho^2] up to which the two-part state bp is taken as pure, or $Failed when it is not taken
   as pure: 0 for a state vector and an exact density matrix, and for an approximate density matrix tol, the round-off
   of its purity, criterionCut times frameRoundoff, within which its computational purity must be 1. "PureStateQ"
   compares the purity with 1 by Equal, which accepts impurities well above tol, so a density matrix it reads as pure
   must meet tol too. It reads the purity only below 2^10 dimensions, and of the matrix as stored: a numeric density
   matrix of 2^10 dimensions or more, or split from a state stored in another frame, that it does not read as pure is
   pure when its computational purity is 1, exactly, or within tol where tol <= 10^-10, as for pureTwoPartQ: a larger
   tol would take mixed states as pure. *)
pureImpurity[qs_, bp_] := Which[
    bp["VectorQ"],
    0,
    TrueQ[bp["PureStateQ"]],
    With[{m = Normal @ computationalForm[bp]["DensityMatrix"]},
        If[ approximateNumericArrayQ[m],
            With[{tol = frameRoundoff[qs] criterionCut[m]}, If[Abs[1 - Total[m Transpose[m], 2]] <= tol, tol, $Failed]],
            0
        ]
    ],
    bp["NumericQ"] && (bp["Dimension"] >= 2 ^ 10 || ! qs["ComputationalQ"]),
    With[{m = Normal @ computationalForm[bp]["DensityMatrix"]},
        If[ approximateNumericArrayQ[m],
            With[{e = Abs[1 - Total[m Transpose[m], 2]]},
                If[e > 10 ^ -10, $Failed, With[{tol = frameRoundoff[qs] criterionCut[m]}, If[tol <= 10 ^ -10 && e <= tol, tol, $Failed]]]
            ],
            If[pureDensityMatrixQ[m], 0, $Failed]
        ]
    ],
    True,
    $Failed
]

(* The state in the computational frame, which a split state usually is already. *)
computationalForm[qs_] := If[qs["ComputationalQ"], qs, qs["Computational"]]

(* The factor by which the change from the frame of the approximate state qs to the computational one multiplies its
   round-off, as in frameContext: kappa for a state vector and kappa^2 for a density matrix, 1 in the computational
   frame. *)
frameRoundoff[qs_] := If[qs["ComputationalQ"], 1, frameCondition[qs["Basis"]] ^ If[qs["VectorQ"], 1, 2]]

(* The Renyi entropy of order alpha of the weights p, in bits: Log2[Total[p^alpha]] / (1 - alpha), and at order 1,
   where that is 0 / 0, its limit, the von Neumann entropy -p . Log2[p], to the precision of the order. An approximate
   order is read as the number it holds to its own precision, SetPrecision[alpha, Infinity], which is how the callers
   raise it (1.00001 given with 3 digits is 1), and the formula is taken whenever that number is not 1, with the
   digits its 1 - alpha cancels added by the callers (orderDigits). A symbolic order gives both, so that the value has
   that limit when 1 is substituted for the order. *)
renyiEntropy[p_, alpha_ ? NumericQ] :=
    If[SetPrecision[alpha, Infinity] == 1, N[- p . Log2[p], Precision[alpha]], Log2[Total[p ^ alpha]] / (1 - alpha)]
renyiEntropy[p_, alpha_] := Piecewise[{{- p . Log2[p], alpha == 1}}, Log2[Total[p ^ alpha]] / (1 - alpha)]

(* The Renyi entropy of the pure two-part state bp split from qs, taken as pure up to the impurity given, read in the
   computational frame off its Schmidt weights, the nonzero eigenvalues of its reduced state. *)
pureStateRenyiEntropy[qs_, bp_, alpha_, impurity_] := With[{c = computationalForm[bp], d = bp["Dimensions"]},
    With[{a = Normal @ If[c["VectorQ"], c["StateVector"], c["DensityMatrix"]]},
        Which[
            approximateNumericArrayQ[a],
            approximateRenyiEntropy[a, d, alpha, frameRoundoff[qs], impurity],
            c["VectorQ"], exactRenyiEntropy[With[{m = ArrayReshape[a, d]}, m . ConjugateTranspose[m]], alpha],
            True, exactRenyiEntropy[TensorContract[ArrayReshape[a, Join[d, d]], {{2, 4}}], alpha]
        ]
    ]
]

(* The weights of an approximate state are its squared Schmidt coefficients s, the singular values of its amplitudes
   as a dA x dB matrix: those of a state vector, or those a pure density matrix holds in its column with the largest
   diagonal entry (densityMatrixAmplitudes), which are known to the same digits p relative to the largest. Each
   amplitude is then off by at most 10^-p times the largest, and a singular value by at most Sqrt[n] 10^-p times the
   largest of them, n = dA dB, times the factor k = frameRoundoff of the frame the state was given in. A density
   matrix taken as pure up to an impurity e (pureImpurity) may hold a mixture whose smaller part, of weight at most
   e / 2, puts at most Sqrt[n] e / 2 into its column. A coefficient within 10 Sqrt[n] Max[k 10^-p, e] times the
   largest is dropped: an order below 1 would raise its weight w to w^alpha, far above round-off, and read a product
   state given with 10 digits, or a separable mixture taken as pure, as entangled. Order 0 thus counts the
   coefficients above that cut. The entropy is evaluated at twice the digits p, and as many more as 1 - alpha cancels
   (orderDigits), so that a weight of 10^-24 is not lost in the 1 + 10^-24 it enters. The value is a machine number
   when the state or the order is one. Otherwise it is labeled by the errors it can carry (labeledEntropy): the
   relative error of the smallest weight kept, twice that of its coefficient, 2 Sqrt[n] k Max[s] / Min[kept] 10^-p,
   which no entropy amplifies at any order; the entropy that the min(dA, dB) - Length[kept] coefficients below the
   cut could hold, at most cut / Max[s] each (hiddenEntropy); and the change of the value across the order's own
   uncertainty (orderSpread). *)
approximateRenyiEntropy[a_, d_, alpha_, k_, impurity_] := With[{
    s = SingularValueList[ArrayReshape[If[VectorQ[a], a, densityMatrixAmplitudes[a]], d], Tolerance -> 0],
    p = If[Precision[a] === MachinePrecision, MachinePrecision, arrayDigits[a]]
},
    With[{top = Max[s], cut = 10 Sqrt[Times @@ d] Max[k 10 ^ -p, impurity] Max[s], wp = 2 p + orderDigits[alpha]},
        With[{kept = Select[s, # == top || # > cut &], order = If[Precision[alpha] < Infinity, SetPrecision[SetPrecision[alpha, Infinity], wp], alpha]},
            With[{w = # / Total[#] & [SetPrecision[kept, wp] ^ 2]},
                With[{v = SetPrecision[renyiEntropy[w, order], wp]},
                    Which[
                        MemberQ[{p, Precision[alpha]}, MachinePrecision],
                        N[v, MachinePrecision],
                        ! NumericQ[alpha],
                        N[v, p - Log10[2 Sqrt[Times @@ d] k top / Min[kept]]],
                        True,
                        labeledEntropy[
                            v,
                            hiddenEntropy[Min[d] - Length[kept], cut / top, order, wp] + orderSpread[w, alpha, v, wp],
                            p - Log10[2 Sqrt[Times @@ d] k top / Min[kept]]
                        ]
                    ]
                ]
            ]
        ]
    ]
]

(* The value v, off by at most delta and known to at most prec digits, labeled with its error err, the larger of
   twice delta (lowering a number's precision also rounds it, by up to half its new error) and |v| 10^-prec: as a
   precision when v exceeds err; otherwise as an accuracy that covers v itself, 2 (|v| + err), since the number can
   round to 0 at that accuracy; a zero with no error to prec digits. *)
labeledEntropy[v_, delta_, prec_] := With[{err = Max[2 delta, Abs[v] 10 ^ -prec]},
    Which[
        err == 0, SetAccuracy[v, prec],
        Abs[v] <= err, SetAccuracy[v, - Log10[2 (Abs[v] + err)]],
        True, SetPrecision[v, Log10[Abs[v] / err]]
    ]
]

(* The Renyi entropy that h Schmidt coefficients of at most c times the largest could hold, the bound on what the
   coefficients dropped as round-off leave out, with their weights at most equal to the largest one's; 0 at order 0,
   which counts the coefficients kept, and when none was dropped. *)
hiddenEntropy[h_, c_, order_, wp_] := If[ h < 1 || SetPrecision[order, Infinity] == 0,
    0,
    With[{w = Min[c ^ 2, 1 / (h + 1)]}, Abs[renyiEntropy[SetPrecision[Join[{1 - h w}, ConstantArray[w, h]], wp], order]]]
]

(* How far the entropy v of the weights w moves when an approximate order alpha moves by its own relative
   uncertainty 10^-Precision[alpha] either way, evaluated with as many more digits as that cancels; 0 for an exact or
   symbolic order. The entropy changes with the order by about alpha |Log[w]| times the order's relative change, for
   the smallest weight w, which can far exceed the order's own uncertainty. *)
orderSpread[w_, alpha_, v_, wp_] := If[ ! NumericQ[alpha] || Precision[alpha] === Infinity,
    0,
    With[{x = SetPrecision[alpha, Infinity], e = 10 ^ -Precision[alpha], wq = wp + Ceiling[Precision[alpha]]},
        Max[Abs[renyiEntropy[SetPrecision[w, wq], SetPrecision[x (1 + #), wq]] - v] & /@ {e, -e}]
    ]
]

(* The weights of an exact or symbolic state are the eigenvalues of its reduced state r that are not zero; an integer
   order above 1 takes Tr[r^alpha], a polynomial in the entries of r, without them. An approximate order is applied,
   as the number it holds, to the weights at twice its digits, as many more as 1 - alpha cancels, and as many as the
   smallest weight needs, so that it is not lost in the 1 + w^alpha it enters; the value is a machine number for a
   machine order, and otherwise is labeled by its change across the order's own uncertainty (orderSpread). *)
exactRenyiEntropy[r_, alpha_] := Which[
    IntegerQ[alpha] && alpha > 1,
    Log2[Tr[MatrixPower[r, alpha]]] / (1 - alpha),
    Precision[alpha] < Infinity,
    With[{w = Select[Eigenvalues[r], ! PossibleZeroQ[#] &]},
        With[{wp = 2 Precision[alpha] + orderDigits[alpha] + If[VectorQ[w, NumericQ], Max[0, Ceiling[- Log10[Min[N[w]]]]], 0]},
            With[{wn = N[w, wp]},
                With[{v = SetPrecision[renyiEntropy[wn, SetPrecision[SetPrecision[alpha, Infinity], wp]], wp]},
                    If[ Precision[alpha] === MachinePrecision || ! VectorQ[w, NumericQ],
                        N[v, Precision[alpha]],
                        labeledEntropy[v, orderSpread[wn, alpha, v, wp], wp]
                    ]
                ]
            ]
        ]
    ],
    True,
    renyiEntropy[Select[Eigenvalues[r], ! PossibleZeroQ[#] &], alpha]
]

(* The digits that 1 - alpha cancels, read off the number a numeric order holds to its own precision, for an order
   other than 1. *)
orderDigits[alpha_] := With[{x = SetPrecision[alpha, Infinity]}, If[NumericQ[x] && x != 1, Max[0, Ceiling[- Log10[Abs[1 - x]]]], 0]]

(* The realignment value is the statistic of a separability test, not an entanglement monotone, which a direct call
   says in a message; QuantumEntangledQ reads it through realignmentCriterion, without one. *)
entanglementMonotone[qs_, biPartition_, "Realignment"] := (
    warnUnphysicalMonotone[qs];
    Enclose @ With[{bp = ConfirmBy[qs["Bipartition", biPartition]["Normalized"], QuantumStateQ[#] && #["Qudits"] == 2 &]},
        Message[QuantumEntanglementMonotone::realignment];
        realignmentValue[bp]
    ]
)

(* Every singular value counts: for approximate input the default tolerance drops those below about 100 machine
   epsilons of the largest, and with them a small realignment value. *)
realignmentValue[bqs_] := Total[SingularValueList[ArrayReshape[Transpose[bqs["Bend"]["Tensor"], 2 <-> 3], bqs["Dimensions"] ^ 2], Tolerance -> 0]] - 1


(* Quantum Discord *)

MutualInformationI[rho_QuantumState, biPartition_ : Automatic] := Enclose @ With[
    {s = ConfirmBy[rho["Bipartition", biPartition], QuantumStateQ[#] && #["Qudits"] == 2 &]},

	QuantumPartialTrace[s, {1}]["Entropy"] + QuantumPartialTrace[s, {2}]["Entropy"] - s["Entropy"]
]

(* J measures the second part of the split, by default in its computational basis, in its own dimension. It is summed
   in bits and returned in them, so that symbolic outcome probabilities do not leave a sum of quantities. *)
MutualInformationJ[rho_QuantumState, qm : _QuantumMeasurementOperator | Automatic : Automatic,biPartition_ : Automatic] := Enclose @ With[
    {s = ConfirmBy[rho["Bipartition", biPartition], QuantumStateQ[#] && #["Qudits"] == 2 &]},
    {m = QuantumMeasurementOperator[Replace[qm, Automatic :> QuantumBasis[Last @ s["Dimensions"]]], {2}][s]},

    Quantity[
        QuantityMagnitude[QuantumPartialTrace[s, {2}]["Entropy"]] -
            m["ProbabilitiesList"] . (QuantityMagnitude[Simplify[QuantumPartialTrace[#, {2}]]["Entropy"]] & /@ m["States"]),
        "Bits"
    ]
]

QuantumDiscord[rho_QuantumState, qm : _QuantumMeasurementOperator | Automatic : Automatic, biPartition_ : Automatic] :=
	Enclose[Confirm[MutualInformationI[rho, biPartition]] - Confirm[MutualInformationJ[rho, qm, biPartition]]]

entanglementMonotone[qs_, biPartition_, "MutualInformationI"] :=
    (warnUnphysicalMonotone[qs]; MutualInformationI[qs, biPartition])

entanglementMonotone[qs_, qm_, biPartition_, "MutualInformationJ"] :=
    (warnUnphysicalMonotone[qs]; MutualInformationJ[qs, qm, biPartition])

entanglementMonotone[qs_, qm_, biPartition_, "Discord"] :=
    (warnUnphysicalMonotone[qs]; QuantumDiscord[qs, qm, biPartition])

