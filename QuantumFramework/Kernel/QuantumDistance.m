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
matrixDistance["Trace"] = traceDistance;
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

(* One minus the fidelity of the two matrices. Round-off can leave the value just below 0, which the square
   root in "Bures" would turn into an imaginary part: a Real less than $distanceTolerance below 0 is read as
   the machine number 0., whatever its precision; a value further below 0, from a non-physical input, and
   every symbolic value, are returned as they are. *)
fidelityDistance[r_, s_] := roundoffZero[1 - fidelity[r, s]]

(* Below machine precision the arbitrary-precision linear algebra has too few digits to work with, and
   SingularValueList can fail to converge or run for minutes, as it does on some pairs of 4 x 4 states at
   precision 2. The fidelity and the trace distance are then computed from the machine numbers of the input
   and given the accuracy the input and the machine computation support. With every entry known to within
   10^-a, each state is known to within d^(3/2) 10^-a in trace norm. The trace distance moves by no more
   than that. The fidelity, by the Powers-Stormer inequality, which bounds the squared Hilbert-Schmidt
   distance of the square roots of two states by their trace-norm distance, moves by no more than
   2 (d^(3/2) 10^-a)^(1/2), so it keeps about half the digits of its input; the eigenvalues the machine
   computation drops below its cut add at most 2 Sqrt[(d - 1) cut], and its round-off about
   d $MachineEpsilon. *)
fidelity[r_, s_] /; lowPrecisionMatricesQ[r, s] := With[{a = workingPrecision[{r, s}], d = Length[r]},
    SetAccuracy[
        fidelity[N[r], N[s]],
        -Log10[2 Sqrt[d^(3 / 2) 10^-a] + 2 Sqrt[(d - 1) N[100 10^-MachinePrecision]] + d $MachineEpsilon]
    ]
]

(* The fidelity Tr[Sqrt[Sqrt[r] . s . Sqrt[r]]] of two matrices of inexact numbers, from the eigensystem of
   each. With r = Sum_i a_i |u_i><u_i| and s = Sum_j b_j |v_j><v_j|, it is the sum of the singular values of
   Sqrt[r] . Sqrt[s], and so of M_ij = Sqrt[a_i] <u_i|v_j> Sqrt[b_j], which is the same matrix between
   factors with orthonormal columns. Round-off moves an eigenvalue of a Hermitian matrix, and a singular
   value, no further than the norm of the round-off, however small the overlap of the two states. The
   eigenvalues of r . s, which is not Hermitian, have no such bound: for two nearly orthogonal states the
   round-off in them grows, relative to the small eigenvalues that carry the fidelity, as the overlap
   shrinks, and the square roots of those that are round-off of 0 put an error far above the round-off into
   their sum. An input that is not physical enters as the positive part of its Hermitian part. The entries
   that are exactly 0 in an arbitrary-precision input stay exact under N, and a fidelity made of them alone,
   as that of two orthogonal states in the computational basis, is an exact 0, which is given the accuracy
   of the input. *)
fidelity[r_, s_] /; inexactMatricesQ[r, s] := With[{p = workingPrecision[{r, s}]},
    With[{f = inexactFidelity[N[r, p], N[s, p], p]}, If[Precision[f] === Infinity && p =!= MachinePrecision, SetAccuracy[f, p], f]]
]

(* The fidelity of exact or symbolic matrices, from the eigenvalues of r . s, which are those of
   Sqrt[r] . s . Sqrt[r]: the trace of a square root is the sum of the square roots of the eigenvalues, so
   no square root of a matrix is formed. An exact matrix has radicals or Root objects for eigenvalues, and
   the sum of their square roots stays small, where the square root of the matrix carries its eigenvectors
   in every entry and grows far faster with the dimension; a symbolic sum has no difference of two
   eigenvalues in a denominator, so it stays finite at parameter values where they coincide; and a product
   that is not diagonalizable, which only a non-physical input gives, can have no square root, as the
   nilpotent |+><+| . Z has none, while its eigenvalues always exist. *)
fidelity[r_, s_] := Re[Total[Sqrt[Eigenvalues[r . s]]]]

inexactMatricesQ[r_, s_] := Precision[{r, s}] < Infinity && MatrixQ[r, NumericQ] && MatrixQ[s, NumericQ]

lowPrecisionMatricesQ[r_, s_] := inexactMatricesQ[r, s] && workingPrecision[{r, s}] < $MachinePrecision

(* The working precision of an inexact matrix, or of a list of matrices: machine precision for machine
   numbers, and otherwise the least accuracy of the entries, the number of digits known after the decimal
   point. The entries of a density matrix are no larger than 1, so their accuracy measures what is known of
   the state; the precision of an entry falls as the entry shrinks, to 0 for a zero known only to an
   accuracy, as the 0``29.8 that Cos[N[Pi, 30] / 2] gives, and says little about the other entries. *)
workingPrecision[m_] := With[{p = Precision[m]}, If[p === MachinePrecision, p, Accuracy[m]]]

(* Two diagonal matrices commute, and their fidelity is the sum over the diagonal of the square roots of the
   products of their populations, taken as Sqrt[p] Sqrt[q] so that a product of two small populations does
   not leave the range of machine numbers; a negative population, from a non-physical input, counts as 0. *)
inexactFidelity[r_, s_, _] /; diagonalQ[r] && diagonalQ[s] :=
    Total[Sqrt[Ramp[Re[Normal[Diagonal[r]]]]] Sqrt[Ramp[Re[Normal[Diagonal[s]]]]]]
inexactFidelity[r_, s_, p_] := supportFidelity[hermitianSupport[r], hermitianSupport[s], p]

(* A matrix whose every entry off the diagonal is 0, exact or a zero known only to an accuracy, as 0``29.8;
   DiagonalMatrixQ with its default tolerance also takes a machine matrix with small nonzero entries off the
   diagonal for diagonal. *)
diagonalQ[m_] := DiagonalMatrixQ[m, Tolerance -> 0]

(* The sum of the singular values of M from the two supports, {eigenvalues, eigenvectors as rows}; with no
   eigenvalue left in a support, which only a non-physical input gives, the fidelity is 0. *)
supportFidelity[{a_, u_}, {b_, v_}, p_] := If[a === {} || b === {},
    N[0, p],
    Total[SingularValueList[KroneckerProduct[Sqrt[a], Sqrt[b]] Normal[Conjugate[u] . Transpose[v]], Tolerance -> 0]]
]

(* The positive eigenvalues of a diagonal matrix, its diagonal entries as stored, with the unit vectors as
   eigenvectors. No eigensolver adds round-off to them, so each is kept however small: the fidelity of
   diag(1 - 10^-20, 10^-20) with the pure state {10^-12, Sqrt[1 - 10^-24]} lies almost all in the
   population 10^-20. *)
hermitianSupport[m_] /; diagonalQ[m] := With[{l = Re[Normal[Diagonal[m]]]},
    With[{keep = Thread[l > 0]},
        {Pick[l, keep], IdentityMatrix[Length[l], SparseArray][[Pick[Range[Length[l]], keep]]]}
    ]
]

(* The eigenvalues of the Hermitian part of an inexact matrix that are above 100 10^-p times the largest in
   magnitude, at precision p, with their eigenvectors as rows; the others, negative ones included, are
   dropped. Round-off moves the eigenvalues by about 10^-p times the largest, so the cut keeps the square
   roots of round-off of 0 out of the fidelity. A genuine eigenvalue below the cut is dropped too, which
   changes the fidelity by no more than the square root of the cut for each one dropped, and only by that
   much when the other state has its weight on the eigenvector; the precision an arbitrary-precision result
   claims does not account for that change, which can exceed its stated uncertainty at any precision. A
   diagonal matrix keeps every positive entry instead, so two states whose overlap lies in eigenvalues below
   the cut get a fidelity that depends, by up to that much, on whether the basis they are given in makes
   them diagonal. The eigenvectors must be orthonormal. Eigensystem gives the eigenvectors of an inexact
   matrix normalized and independent, but not orthogonal, and at higher precision two eigenvectors of a
   repeated eigenvalue can be orthogonal to far fewer digits than the working precision; there they are the
   columns of q in the Schur decomposition q . t . ConjugateTranspose[q], where q is unitary and t, for a
   Hermitian matrix, is diagonal with the eigenvalues. The kept eigenvectors are then made orthonormal with
   Orthogonalize by Householder reflections, which keeps the span of each leading set of them, and so each
   eigenspace; vectors that are already orthonormal come back changed only by round-off and by a factor of
   modulus 1 on each, which leaves the singular values of M as they are. *)
hermitianSupport[m_] := With[{e = hermitianEigensystem[m]},
    With[{keep = Thread[Re[e[[1]]] > N[100 10^-workingPrecision[m]] Max[Abs[e[[1]]]]]},
        {Re[Pick[e[[1]], keep]], orthonormalRows[Pick[e[[2]], keep]]}
    ]
]

orthonormalRows[v_] := Orthogonalize[v, Method -> "Householder"]

(* The eigensystem of the Hermitian part of m. A SparseArray of machine numbers with both real and
   complex entries gives an unpacked list, on which forming the Hermitian part is many times slower than
   on a packed one, so the list is packed as complex first; only machine numbers are packed, since packing
   would round arbitrary-precision entries to machine numbers. *)
hermitianEigensystem[m_] /; Precision[m] === MachinePrecision :=
    Eigensystem[hermitianPart[With[{n = Normal[m]}, If[Developer`PackedArrayQ[n], n, Developer`ToPackedArray[n, Complex]]]]]
hermitianEigensystem[m_] := With[{qt = SchurDecomposition[hermitianPart[Normal[m]]]}, {Diagonal[qt[[2]]], Transpose[qt[[1]]]}]

hermitianPart[m_] := (m + ConjugateTranspose[m]) / 2

(* Below machine precision, as for the fidelity, the trace distance is computed from the machine numbers of
   the input and given the accuracy its bound d^(3/2) 10^-a, plus the round-off of the machine computation,
   supports. *)
traceDistance[r_, s_] /; lowPrecisionMatricesQ[r, s] := With[{a = workingPrecision[{r, s}], d = Length[r]},
    SetAccuracy[traceDistance[N[r], N[s]], -Log10[d^(3 / 2) 10^-a + d $MachineEpsilon]]
]

(* Two diagonal matrices: the singular values of their difference are the absolute differences of their
   populations. *)
traceDistance[r_, s_] /; inexactMatricesQ[r, s] && diagonalQ[r] && diagonalQ[s] := Total[Abs[Normal[Diagonal[r] - Diagonal[s]]]] / 2

(* Half the trace norm of r - s, the sum of its singular values. SingularValueList finds them without
   forming ConjugateTranspose[r - s] . (r - s), so a singular value of the order of the round-off stays that
   small, where the square root of an eigenvalue of that product would be of the order of the square root
   of the round-off. For inexact input it is given the dense difference, since a SparseArray with few nonzero
   entries, which is how QuantumState keeps a sparse density matrix, is converted to a dense one above
   dimension 100 with a SingularValueList::arh message; and Tolerance -> 0, since its default drops singular
   values below 100 10^-p times the largest, which are well above round-off. The trace norm of a difference
   that is not Hermitian, which only a non-physical input gives, is still the sum of its singular values, not
   of the absolute values of its eigenvalues. *)
traceDistance[r_, s_] /; inexactMatricesQ[r, s] := Total[SingularValueList[Normal[r - s], Tolerance -> 0]] / 2

(* For exact and symbolic input the default tolerance is 0 and only exact zeros are dropped. The singular
   values are radicals or Root objects with no difference of two eigenvalues in a denominator, so they stay
   finite where eigenvalues of ConjugateTranspose[r - s] . (r - s) coincide, as they do at every real point
   for two qubit states; Re removes the zero imaginary part such a value can pick up when it is evaluated
   numerically. *)
traceDistance[r_, s_] := Re @ Total[SingularValueList[r - s]] / 2

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

