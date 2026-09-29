Package["Wolfram`QuantumFramework`"]

PackageExport["QuantumDistance"]
PackageExport["QuantumSimilarity"]

PackageScope["$QuantumDistances"]
PackageScope["numericStateNotPSDQ"]



$QuantumDistances = {"Fidelity", "RelativeEntropy", "RelativePurity", "Trace", "Bures", "BuresAngle", "HilbertSchmidt", "Bloch"}


(* Whether a state is not physical, for the entanglement monotones, by the rule distanceMatrix applies to an
   input. A monotone computes from a state divided by its trace (a state vector divided by its norm): a
   trace that traceClass finds to be 0, negative, or not real is never a state's, a positive trace is
   divided out before notPhysicalMatrixQ tests the matrix, and a trace whose sign the kernel cannot
   decide, or a symbolic one, leaves the matrix to be tested as given. Multiplying a state by a positive
   number therefore leaves the verdict as it is when the kernel can decide the sign of both the original
   and the new trace, and neither the machine entries nor a state vector's squared norm leave the normal
   range of machine numbers. A state vector is physical unless it is zero, which
   zeroVectorQ reads from its dense computational vector without forming its density matrix. The density
   matrix of a mixed state is made dense as well: traceClass reads a machine trace against Norm, which can
   return a wrong value, or crash the kernel, on a SparseArray whose default element is not 0, and
   QuantumState can store an input holding unsimplified exact zeros in such an array. *)
numericStateNotPSDQ[qs_ ? QuantumStateQ] := If[
    qs["VectorQ"],
    zeroVectorQ[Normal @ qs["Computational"]["StateVector"]],
    With[{rho = Normal @ qs["Computational"]["DensityMatrix"]},
        With[{t = Tr[rho]},
            Switch[traceClass[t, rho],
                "Zero" | "Other", True,
                "Positive", notPhysicalMatrixQ[rho / Re[t]],
                _, notPhysicalMatrixQ[rho]
            ]
        ]
    ]
]
numericStateNotPSDQ[_] := False

(* Whether a state vector is 0, from its squared norm: through reducedTrace for exact entries, as traceClass
   compares an exact trace with 0, and for machine entries when the squared norm is 0., as it is too when
   the squared norm underflows, which distanceMatrix then reads as a trace of 0. *)
zeroVectorQ[v_] := With[{t = Conjugate[v] . v},
    Which[
        ! NumericQ[t], False,
        Precision[t] === Infinity, TrueQ[reducedTrace[t] == 0],
        True, t == 0
    ]
]

QuantumDistance::notphysical =
    "An input is not a positive-semidefinite density matrix; a distance or similarity involving a non-physical state may be meaningless."

QuantumDistance::notnormalized =
    "Input `1` has trace `2`, not 1 (for a state vector, the trace is its squared norm); it is divided by its trace."

QuantumDistance::bloch =
    "The \"Bloch\" distance is defined for single-qubit states; these states have dimension `1`. \"Trace\" and \"HilbertSchmidt\" are defined for any dimension."

(* The tolerance below which a machine number is read as round-off: a trace this close to 1, or to 0 in
   units of the Frobenius norm of its matrix, a negative eigenvalue or a Hermiticity defect this small in a
   state of unit trace, and a fidelity distance this far below 0. *)
$distanceTolerance = 1.*^-8


(* The two density matrices a measure is computed from, each input's in the computational basis. An input
   whose trace (for a state vector, its squared norm) is a positive number other than 1 stands for the
   state rho / Tr[rho], and is divided by its trace: an exact trace other than 1 with a notnormalized
   message naming the input, and a machine trace in any case, since round-off leaves it a few units in the
   last place from 1, with the message only when it is more than $distanceTolerance from 1. The trace is
   not simplified: a symbolic input is rescaled only when its trace evaluates to a number, as that of
   {{2 p, 0}, {0, 2 - 2 p}} does, and is otherwise used as given, as is an exact input whose trace has a
   sign the kernel cannot decide. An input whose trace traceClass finds to be 0, negative, or not real, or
   a numeric one that is not Hermitian or has a negative eigenvalue once divided by its trace, is not
   physical, and one notphysical message covers the call. *)
distanceMatrices[qs1_, qs2_] := With[{inputs = MapThread[distanceMatrix, {{qs1, qs2}, {1, 2}}]},
    If[Or @@ inputs[[All, 2]], Message[QuantumDistance::notphysical]];
    inputs[[All, 1]]
]

distanceMatrix[qs_, i_] := With[{rho = qs["Computational"]["DensityMatrix"]},
    With[{t = Tr[rho]}, inputMatrix[traceClass[t, rho], qs, rho, t, i]]
]

(* {the matrix the measures read, whether the input is not physical}; a trace that is 0, negative, or not
   real is never a state's, whatever the scale of the input, and such an input is read as given *)
inputMatrix["Zero" | "Other", _, rho_, _, _] := {rho, True}
inputMatrix["Positive", qs_, rho_, t_, i_] := withPhysicality[qs, rescaled[rho, Re[t], i]]
inputMatrix[_, qs_, rho_, _, _] := withPhysicality[qs, rho]

withPhysicality[qs_, m_] := {m, ! qs["VectorQ"] && notPhysicalMatrixQ[m]}

rescaled[rho_, t_, i_] := (
    If[Precision[t] === Infinity || Abs[t - 1] > $distanceTolerance, Message[QuantumDistance::notnormalized, i, t]];
    rho / t
)

(* The trace of an input as "Zero", "Unit", "Positive" (real and positive, other than an exact 1), "Other"
   (negative or not real), "Undecided" (an exact trace whose sign the kernel cannot decide), or "Symbolic"
   (one that does not evaluate to a number). An exact trace is compared with 0 and 1 through reducedTrace,
   and classed by its sign through signClass. A machine trace is judged against the Frobenius norm of its
   matrix, so a state vector of norm 10^-5 is rescaled like any other. *)
traceClass[t_ ? NumericQ, _] /; Precision[t] === Infinity := With[{u = reducedTrace[t]},
    Which[
        TrueQ[u == 0], "Zero",
        TrueQ[u == 1], "Unit",
        True, signClass[u, Positive[u]]
    ]
]
traceClass[t_ ? NumericQ, rho_] := With[{tolerance = $distanceTolerance frobeniusNorm[rho]},
    Which[
        Abs[t] <= tolerance, "Zero",
        Abs[Im[t]] <= tolerance && Re[t] > tolerance, "Positive",
        True, "Other"
    ]
]
traceClass[_, _] := "Symbolic"

(* An exact trace to compare with 0 and 1: one that is algebraic and whose machine value is near 0 or 1 is
   reduced with RootReduce, since Equal cannot always decide a sum of nested radicals such as
   Expand[(Sqrt[3 + 2 Sqrt[2]] - Sqrt[2])^2], which is 1. An exact rational or complex rational is
   compared exactly as it is, and is not converted to a machine number, which for one below about 10^-308
   raises General::munfl. *)
reducedTrace[t_ ? ExactNumberQ] := t
reducedTrace[t_] := If[Min[Abs[N[t]], Abs[N[t] - 1]] < 1.*^-6 && TrueQ[Element[t, Algebraics]], RootReduce[t], t]

(* The class of a nonzero exact trace u from s = Positive[u]. When Positive cannot decide an algebraic
   trace, RootReduce decides it when the sign of the reduced form can be decided: (1 + I Sqrt[2])
   (1 - I Sqrt[2]) reduces to 3, and a Root object is known exactly to be real or not, even when Positive
   cannot decide the sign of one with a tiny imaginary part; the square of Sqrt[2 + 10^-200] - Sqrt[2]
   stays undecided. A trace whose sign the kernel cannot decide within $MaxExtraPrecision, as that of an
   input multiplied by E^(10^-200) - 1, is "Undecided", and the input is used as given. *)
signClass[_, True] := "Positive"
signClass[_, False] := "Other"
signClass[u_, _] /; TrueQ[Element[u, Algebraics]] := With[{r = RootReduce[u]},
    Which[
        Element[r, Reals] === False, "Other",
        r === u, "Undecided",
        True, Replace[Positive[r], {True -> "Positive", False -> "Other", _ -> "Undecided"}]
    ]
]
signClass[_, _] := "Undecided"

(* The Frobenius norm as the norm of the flattened entries, which keeps a SparseArray sparse and is many
   times faster than Norm[rho, "Frobenius"] on one; 1. for a matrix with a symbolic entry, whose trace is
   then judged on its own. *)
frobeniusNorm[rho_] := With[{f = N @ Norm[Flatten[rho]]}, If[NumericQ[f], f, 1.]]

(* A numeric matrix that is not Hermitian, or that has an eigenvalue below -$distanceTolerance. The check
   runs on the matrix after it is divided by a positive trace, so the tolerance is relative to the state; a
   matrix whose trace is undecided or symbolic is checked as given, against the same tolerance. The
   eigenvalues are machine numbers: an exact negative eigenvalue closer to 0 than the tolerance is not
   seen, which is the price of not deciding the spectrum of an exact matrix exactly, and an exact matrix
   whose entries machine arithmetic cannot evaluate is not checked at all, as DiagonalMatrix[{t, -t/2}]
   divided by its trace, with t = Sqrt[2 + 10^-200] - Sqrt[2]: its entries evaluate to 0./0., with the
   messages that division raises. A matrix of exact rational entries is checked for Hermiticity exactly;
   any other is checked in its dense form, since HermitianMatrixQ with a Tolerance reads a SparseArray more
   strictly than the same matrix dense, and fails a round-off entry opposite a structural zero. *)
notPhysicalMatrixQ[m_] := With[{dm = Normal @ N[m]},
    MatrixQ[dm, NumericQ] && Or[
        ! If[ArrayQ[m, 2, ExactNumberQ], HermitianMatrixQ[m], HermitianMatrixQ[dm, Tolerance -> $distanceTolerance]],
        Min[Re @ Eigenvalues[dm]] < - $distanceTolerance
    ]
]


QuantumDistance[qs1_ ? QuantumStateQ, qs2_ ? QuantumStateQ] := QuantumDistance[qs1, qs2, "Fidelity"]

QuantumDistance[qs1_ ? QuantumStateQ, qs2_ ? QuantumStateQ, "Bloch"] /; qs1["Dimension"] == qs2["Dimension"] != 2 :=
    nonQubitBlochFailure[qs1["Dimension"]]

QuantumDistance[qs1_ ? QuantumStateQ, qs2_ ? QuantumStateQ, measure : Alternatives @@ $QuantumDistances] /; qs1["Dimension"] == qs2["Dimension"] :=
    Apply[
        If[equalMatricesQ[##], equalStateDistance[measure, #1], matrixDistance[measure][##]] &,
        distanceMatrices[qs1, qs2]
    ]

(* Each measure as a function of the two density matrices. *)
matrixDistance["Fidelity"] = fidelityDistance;
matrixDistance["RelativeEntropy"] = relativeEntropy;
matrixDistance["RelativePurity"] = 1 - Chop[Tr[#1 . #2]] &;
(* For symbolic states that do not commute, the matrix square root in "Trace" leaves an expression that can
   evaluate to Indeterminate at numeric parameter values; FullSimplify under assumptions on the parameters
   recovers the value. *)
matrixDistance["Trace"] = With[{m = #1 - #2}, Re @ Tr[MatrixPower[ConjugateTranspose[m] . m, 1 / 2]] / 2] &;
matrixDistance["Bures"] = Sqrt[2 fidelityDistance[##]] &;
matrixDistance["BuresAngle"] = Re @ ArcCos[1 - fidelityDistance[##]] &;
matrixDistance["HilbertSchmidt"] = Norm[#1 - #2, "Frobenius"] &;
matrixDistance["Bloch"] = EuclideanDistance[blochVector[#1], blochVector[#2]] / 2 &;

(* Two matrices are equal when they are identical, or when they are exact and every entry of their
   difference reduces to 0 under RootReduce, as Sqrt[3 + 2 Sqrt[2]] - 1 - Sqrt[2] does. Equal states are at
   distance 0 in every measure, exact or machine as the matrix is, and at 0 bits of relative entropy;
   "RelativePurity" gives 1 - Tr[rho^2], which is 0 only for a pure state. Computed from the formulas
   instead, the distance of an exact state with irrational eigenvalues to itself is a sum of nested
   radicals that is 0 but is never simplified, on which Re and ArcCos cannot decide. *)
equalMatricesQ[r_, s_] := r === s || Precision[{r, s}] === Infinity && With[{nr = N[r], ns = N[s]},
    MatrixQ[nr, NumericQ] && MatrixQ[ns, NumericQ] && Max[Abs[Flatten[nr - ns]]] < $distanceTolerance
] && AllTrue[Flatten[Normal[r - s]], RootReduce[#] === 0 &]

equalStateDistance["RelativePurity", r_] := matrixDistance["RelativePurity"][r, r]
equalStateDistance["RelativeEntropy", r_] := Quantity[zeroLike[r], "Bits"]
equalStateDistance[_, r_] := zeroLike[r]

zeroLike[r_] := If[Precision[r] === Infinity, 0, 0.]

(* One minus the fidelity Tr[Sqrt[Sqrt[r] . s . Sqrt[r]]], computed as Tr[Sqrt[r . s]], which has the same
   eigenvalues under the square root. For machine input the value can land a few units in the last place
   below 0, which the square root in "Bures" would turn into an imaginary part, and such a value is read as
   0.; a value further below 0, from a non-physical input, and every symbolic value, are returned as they
   are. *)
fidelityDistance[r_, s_] := roundoffZero[1 - Re[Tr[MatrixPower[r . s, 1 / 2]]]]

roundoffZero[d_Real] /; - $distanceTolerance < d < 0 := 0.
roundoffZero[d_] := d

(* Umegaki relative entropy, S(s||t) = Tr[s log s] - Tr[s log t], where the
   cross term is taken as Sum_j <t_j|s|t_j> log q_j over the eigenbasis of t
   alone - a density matrix is routinely singular and 0 log 0 = 0 is what
   decides the answer there. A zero eigenvalue of s drops out of the entropy,
   which keeps a pure s against a full-rank t finite; a zero q_j carrying weight
   sends its term to -Infinity, so a violated support condition
   supp(s) subset of supp(t) surfaces as Quantity[Infinity, "Bits"], which
   QuantumSimilarity carries through as 2^-Infinity == 0. Chop supplies the
   tolerance an inexact spectrum needs and leaves exact input untouched. The
   sum only telescopes to Tr[s log t] on an orthonormal eigenbasis, which
   Eigensystem does not supply inside a degenerate eigenspace, hence
   "Orthogonalize". *)
relativeEntropy[s_, t_] := With[{positive = ! TrueQ[Re[#] <= 0] &},
    Apply[
        Function[{vals, vecs},
            Quantity[
                Chop[(
                    With[{p = Select[Chop @ eigenvalues[Normal @ s], positive]}, Total[p Log[p]]] -
                    Total @ MapThread[If[positive[#2], #2 Log[#1], 0] &, Chop @ {
                        If[positive[#], #, 0] & /@ vals,
                        Conjugate[#] . s . # & /@ vecs
                    }]
                ) / Log[2]],
                "Bits"
            ]
        ],
        eigensystem[Normal @ t, "Normalize" -> True, "Orthogonalize" -> True]
    ]
]

(* The Bloch vector {Tr[m . X], Tr[m . Y], Tr[m . Z]} of a qubit density matrix m, in the generators
   GellMannMatrices[2] that "BlochVector" also uses; half the Euclidean distance between two Bloch vectors
   is the trace distance of the two states. In a higher dimension the generalized Bloch vector gives only a
   multiple of the Hilbert-Schmidt distance, so the measure is defined for a single qubit, and any other
   dimension gives a Failure, as QuantumQASM does for a circuit that is not on qubits. *)
blochVector[m_] := Tr[m . #] & /@ GellMannMatrices[2]

nonQubitBlochFailure[d_] := Failure["NonQubitBloch", <|
    "MessageTemplate" :> QuantumDistance::bloch,
    "MessageParameters" -> {d},
    "Dimension" -> d
|>]


QuantumSimilarity[qs1_ ? QuantumStateQ, qs2_ ? QuantumStateQ, measure : Alternatives @@ $QuantumDistances : "Fidelity"] :=
    With[{d = QuantumDistance[qs1, qs2, measure]}, If[FailureQ[d], d, unitInterval @ similarity[measure, d]]]

(* Each similarity from its distance d: 1 - d for the measures bounded by 1, 1 - d / Sqrt[2] and
   1 - d / (Pi / 2) for those bounded by Sqrt[2] and Pi / 2, and an exponential decay for the unbounded
   relative entropy. A state has similarity 1 with itself in every measure but "RelativePurity", whose
   similarity Tr[rho . sigma] is then the purity Tr[rho^2], 1 only for a pure state. *)
similarity["Bures" | "HilbertSchmidt", d_] := 1 - d / Sqrt[2]
similarity["BuresAngle", d_] := 1 - d / (Pi / 2)
similarity["RelativeEntropy", d_] := 2 ^ (-Replace[d, q_Quantity :> QuantityMagnitude[q]])
similarity[_, d_] := 1 - d

(* A machine similarity that round-off leaves just outside [0, 1], as the self RelativePurity similarity of
   a pure state a unit in the last place above 1, is returned at the end of the interval. The differences
   are compared with 0, since a comparison of two machine numbers such as 1. + 2.^-52 > 1 ignores the last
   bits; a value further outside, from a non-physical input, is returned as it is. *)
unitInterval[x_Real] /; - $distanceTolerance < x < 0 := 0.
unitInterval[x_Real] /; 0 < x - 1 < $distanceTolerance := 1.
unitInterval[x_] := x

