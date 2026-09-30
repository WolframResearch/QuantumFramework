Package["Wolfram`QuantumFramework`"]

PackageExport["QuantumDistance"]
PackageExport["QuantumSimilarity"]

PackageScope["$QuantumDistances"]
PackageScope["numericStateNotPSDQ"]



$QuantumDistances = {"Fidelity", "RelativeEntropy", "RelativePurity", "Trace", "Bures", "BuresAngle", "HilbertSchmidt", "Bloch"}


(* Whether a state is not physical, for the entanglement monotones. QuantumDistance takes its notphysical
   verdict on an input from the same functions on the same matrix, so the two warn on the same states. A
   monotone computes from a state divided by its trace (a state vector divided by its norm), and a state
   vector is physical unless it is 0, which zeroVectorQ reads without forming its density matrix; a mixed
   state is judged by densityReading. Multiplying a machine state by a positive machine number, or a state
   of exact rational or complex rational entries by a positive rational number, leaves the verdict as it
   is when the entries, before and after the division by the trace, and a state vector's squared norm lie
   in the normal range of machine numbers, and no quantity compared with a tolerance lies within round-off
   of it. Another multiplier can change the Hermiticity test from exact to one with the tolerance: a
   machine number makes an exact state machine, and an irrational number can leave a trace that does not
   cancel from the divided entries, as the trace 1/(3 Sqrt[2]) + Sqrt[2]/3 of diag(1/3, 2/3) / Sqrt[2]
   does. The state is read in its dense form: QuantumState can store an input holding several unsimplified
   exact zeros in a SparseArray whose default element is such a zero, and Norm, which traceScale takes of a
   matrix, can crash the kernel on such an array. *)
numericStateNotPSDQ[qs_ ? QuantumStateQ] := If[
    qs["VectorQ"],
    zeroVectorQ[Normal @ qs["Computational"]["StateVector"]],
    Last @ densityReading[Normal @ qs["Computational"]["DensityMatrix"]]
]
numericStateNotPSDQ[_] := False

(* {the class of the trace of a density matrix, the trace, whether the matrix is not physical}: a trace that
   traceClass finds to be 0, negative, or not real is never a state's, a positive trace is divided out
   before notPhysicalMatrixQ tests the matrix, and a trace whose sign the kernel cannot decide, or a
   symbolic one, leaves the matrix to be tested as given. Tr raises General::munfl when machine entries in
   the normal range cancel to a nonzero trace below it, as those of diag(5 10^-301, -4.99999999999 10^-301)
   do. *)
densityReading[rho_] := With[{t = Tr[rho]},
    With[{class = traceClass[t, rho]},
        {class, t, Switch[class,
            "Zero" | "Other", True,
            "Positive", notPhysicalMatrixQ[rho / Re[t]],
            _, notPhysicalMatrixQ[rho]
        ]}
    ]
]

(* Whether a state vector is 0, from its squared norm: through reducedTrace for exact entries, as traceClass
   compares an exact trace with 0, and for machine entries when the squared norm is 0., as it is too when
   the squared norm underflows. *)
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
   units of the scale traceScale gives its matrix, a negative eigenvalue this small in a state of unit
   trace, and a fidelity distance this far below 0. It is also the Tolerance of the Hermiticity test, under
   which HermitianMatrixQ takes an entry of absolute value at most 10^-8 for 0 and compares larger entries
   except for their last Log2[10^-8 / $MachineEpsilon] bits, about 25.4, so that two of them may differ
   by a relative 5 10^-9. *)
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

(* {the matrix the measures read, whether the input is not physical}. The verdict is numericStateNotPSDQ's:
   zeroVectorQ's for a state vector, and for a mixed state densityReading's of the dense matrix, which also
   gives the trace and its class that decide the division. The measures read the matrix as the state
   stores it, divided by its trace when that is positive; any other is read as given. *)
distanceMatrix[qs_, i_] := With[{rho = qs["Computational"]["DensityMatrix"]},
    inputMatrix[rho, i, If[
        qs["VectorQ"],
        With[{t = Tr[rho]}, {traceClass[t, rho], t, numericStateNotPSDQ[qs]}],
        densityReading[Normal @ rho]
    ]]
]

inputMatrix[rho_, i_, {"Positive", t_, notPhysical_}] := {rescaled[rho, Re[t], i], notPhysical}
inputMatrix[rho_, _, {_, _, notPhysical_}] := {rho, notPhysical}

rescaled[rho_, t_, i_] := (
    If[Precision[t] === Infinity || Abs[t - 1] > $distanceTolerance, Message[QuantumDistance::notnormalized, i, t]];
    rho / t
)

(* The trace of an input as "Zero", "Unit", "Positive" (real and positive, other than an exact 1), "Other"
   (negative or not real), "Undecided" (an exact trace whose sign the kernel cannot decide), or "Symbolic"
   (one that does not evaluate to a number). An exact trace is compared with 0 and 1 through reducedTrace,
   and classed by its sign through signClass. A machine trace is judged against the scale traceScale gives
   its matrix, so a state vector of norm 10^-5 is rescaled like any other. *)
traceClass[t_ ? NumericQ, _] /; Precision[t] === Infinity := With[{u = reducedTrace[t]},
    Which[
        TrueQ[u == 0], "Zero",
        TrueQ[u == 1], "Unit",
        True, signClass[u, Positive[u]]
    ]
]
traceClass[t_ ? NumericQ, rho_] := If[t == 0, "Zero", With[{s = traceScale[t, rho]},
    If[s < 1, scaledTraceClass[t / s, 1], scaledTraceClass[t, s]]
]]
traceClass[_, _] := "Symbolic"

(* The class of a machine trace t against the scale s: 0 within $distanceTolerance s, and positive when its
   real part exceeds that and its imaginary part does not. traceClass passes t / s against 1 when s is
   below 1, so that, for a trace in the normal range of machine numbers, neither that quotient nor the
   product $distanceTolerance s falls below it. A comparison that does not evaluate to True, as one with
   Indeterminate, leaves the trace "Other". *)
scaledTraceClass[t_, s_] := Which[
    TrueQ[Abs[t] <= $distanceTolerance s], "Zero",
    TrueQ[Abs[Im[t]] <= $distanceTolerance s && Re[t] > $distanceTolerance s], "Positive",
    True, "Other"
]

(* An exact trace to compare with 0 and 1: one that is algebraic and lies within 10^-6 of 0 or 1 is reduced
   with RootReduce, since Equal cannot always decide a sum of nested radicals such as
   Expand[(Sqrt[3 + 2 Sqrt[2]] - Sqrt[2])^2], which is 1. Less compares the trace with the exact bound
   10^-6 without the General::munfl that converting a trace below about 10^-308 to a machine number
   raises; like Positive, it raises Less::meprec for a trace closer to the bound than $MaxExtraPrecision
   resolves, as 10^-6 + Sqrt[2 + 10^-200] - Sqrt[2], which is then not reduced. An exact rational or
   complex rational is compared exactly as it is. *)
reducedTrace[t_ ? ExactNumberQ] := t
reducedTrace[t_] := If[
    (TrueQ[Abs[t] < 10^-6] || TrueQ[Abs[t - 1] < 10^-6]) && TrueQ[Element[t, Algebraics]],
    RootReduce[t],
    t
]

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

(* The scale a nonzero machine trace t is judged against: the larger of the absolute value of t and the
   Frobenius norm of the numeric entries of its matrix, all of them for a numeric matrix, taken as the norm
   of the flattened entries. Both scale with the matrix, so the class of t does not depend on its scale: t
   is 0 when its absolute value is at most $distanceTolerance times that norm, positive when its real part
   exceeds, and its imaginary part is at most, $distanceTolerance times the scale, and "Other" otherwise.
   The norm is not converted to a machine number, which for arbitrary-precision entries below the range of
   machine numbers would underflow to 0. *)
traceScale[t_, rho_] := Max[Abs[t], numericEntriesNorm[Flatten[rho]]]
numericEntriesNorm[entries_] :=
    If[VectorQ[entries, NumericQ], Norm[entries], numericNorm[Select[Normal[entries], NumericQ]]]
numericNorm[{}] := 0
numericNorm[entries_] := Norm[entries]

(* A numeric matrix that is not Hermitian, or that has an eigenvalue below -$distanceTolerance. The check
   runs on the matrix after it is divided by a positive trace, so the tolerance is relative to the state; a
   matrix whose trace is undecided or symbolic is checked as given, against the same tolerance. The
   eigenvalues are those of N of the matrix, so an exact negative eigenvalue closer to 0 than the
   tolerance is not seen, which is the price of not deciding the spectrum of an exact matrix exactly. N
   turns an exact or arbitrary-precision entry below the range of machine numbers into 0., with
   General::munfl, and an exact one above that range into an arbitrary-precision number, whose eigenvalues
   raise General::munfl, as the entries 10^400 of {{10^-400, 1}, {1, 0}} divided by its trace do. An
   exact matrix whose entries machine arithmetic cannot evaluate is not checked at all, as
   DiagonalMatrix[{t, -t/2}] divided by its trace, with t = Sqrt[2 + 10^-200] - Sqrt[2]: its entries
   evaluate to 0./0., with the Power::infy and Infinity::indet messages of that division. A matrix whose
   entries, as divided, are exact rational or complex rational numbers is checked for Hermiticity exactly;
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

(* The fidelity Tr[Sqrt[Sqrt[r] . s . Sqrt[r]]] of two matrices of inexact numbers, from the eigensystems
   of the two states, or, above machine precision, for a state that is pure to within the error of its
   entries, from a column of its matrix. With r = Sum_i a_i |u_i><u_i| and s = Sum_j b_j |v_j><v_j|, it is
   the sum of the singular values of Sqrt[r] . Sqrt[s], and so of M_ij = Sqrt[a_i] <u_i|v_j> Sqrt[b_j],
   which is the same matrix between factors with orthonormal columns. Round-off moves an eigenvalue of a
   Hermitian matrix, and a singular value, no further than the norm of the round-off, however small the
   overlap of the two states. The eigenvalues of r . s, which is not Hermitian, have no such bound: for two
   nearly orthogonal states the round-off in them grows, relative to the small eigenvalues that carry the
   fidelity, as the overlap shrinks, and the square roots of those that are round-off of 0 put an error far
   above the round-off into their sum. An input that is not physical enters as the positive part of its
   Hermitian part. *)
fidelity[r_, s_] /; inexactMatricesQ[r, s] := With[{p = workingPrecision[{r, s}]},
    withClaimedAccuracy[inexactFidelity[N[r, p], N[s, p], p], p]
]

(* The fidelity from {fidelity, bound}, the bound being on the change that the eigenvalues dropped below the
   cut and resolved above their uncertainty make in it. A machine fidelity claims no accuracy. Above machine
   precision a fidelity claims an uncertainty of twice the larger of the one significance arithmetic gives
   it and the bound: a number read exactly, as SetPrecision[x, Infinity] reads it, is rounded at the
   accuracy it carries, by up to half of its uncertainty, and the second half leaves room for that
   rounding. The entries that are exactly 0 in an arbitrary-precision input stay exact under N, and a
   fidelity made of them alone, as that of two orthogonal states in the computational basis, is an exact 0,
   which is given the accuracy of the input first. Significance arithmetic follows the round-off of the
   computation, not the sensitivity of the square root to the input, and an eigenvalue that the arithmetic
   cannot tell from its uncertainty d 10^-p, up to a few times that, is taken as 0. Where one state has such
   an eigenvalue, or a kept one near 0, that the other state weighs, the fidelity can move by more than the
   accuracy it claims when the entries move within theirs: by up to about the square root of a few times
   d 10^-p for an eigenvalue taken as 0, and, for a 30-digit state vector against a mixed state whose
   weight on it is 10^-16, by about 10^-23 when the entries of the mixed state move by 10^-30. Bounding that
   for every input, as the Powers-Stormer inequality does below machine precision, would leave every
   arbitrary-precision fidelity about half the digits of its input. *)
withClaimedAccuracy[{f_, _}, MachinePrecision] := f
withClaimedAccuracy[{f_, b_}, p_] := With[{g = If[Precision[f] === Infinity, SetAccuracy[f, p], f]},
    SetAccuracy[g, Min[Accuracy[g], If[TrueQ[b > 0], - Log10[b], Infinity]] - Log10[2]]
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
    {Total[Sqrt[Ramp[Re[Normal[Diagonal[r]]]]] Sqrt[Ramp[Re[Normal[Diagonal[s]]]]]], 0}
inexactFidelity[r_, s_, p : MachinePrecision] := {supportFidelity[hermitianSupport[r], hermitianSupport[s], p], 0}
inexactFidelity[r_, s_, p_] := formFidelity[stateForm[r], stateForm[s], r, s, p]

(* Above machine precision the eigensystem is the costly step, and a state that is pure to within the error
   of its entries needs none. For m = l |u><u|, the column c of m through its largest diagonal entry m_kk is
   Sqrt[l m_kk] |u> up to a phase, and c c^dagger / m_kk is m. When that rank-one matrix lies within d 10^-p
   of m in Frobenius norm, every other eigenvalue of m lies within d 10^-p of 0, whatever its sign, since no
   matrix of rank one is closer to m than the root of the sum of their squares: none of them is resolved,
   and they are taken as 0. The cut drops them too while d 10^-p lies below it, that is for d below about
   100 times the largest eigenvalue; above that they include eigenvalues the cut would keep. The distance is
   computed at precision p, and is known only to a few times d 10^-p, so the comparison, which holds when
   the two agree to within that, lets through eigenvalues up to a few times d 10^-p. A state with more
   weight off u, a resolved small eigenvalue or a negative one of a non-physical input included, takes its
   eigensystem. As for the eigensystem, m is read through its Hermitian part. A diagonal matrix keeps its
   populations, as below. *)
stateForm[m_] /; diagonalQ[m] := hermitianSupport[m]
stateForm[m_] := With[{n = hermitianPart[Normal[m]]},
    With[{k = First[Ordering[Re[Diagonal[n]], -1]]},
        With[{c = n[[All, k]], mkk = Re[n[[k, k]]]},
            If[TrueQ[mkk > 0 && Norm[Flatten[n - Outer[Times, c, Conjugate[c]] / mkk]] <= N[Length[n] 10^-workingPrecision[m]]],
                rankOne[c, mkk],
                hermitianSupport[m]
            ]
        ]
    ]
]

(* {fidelity, bound} from the forms of the two states. For a pure state l |u><u| and a state s the fidelity
   is Sqrt[l <u|s|u>], taken as the norm of the overlaps of u with the eigenvectors of s weighted by the
   square roots of their eigenvalues: each overlap is computed directly, and the sum of their squares has
   no cancellation, where <u|s|u> summed over the entries of s can cancel to far below the entries. For two
   pure states it is Sqrt[l l'] |<u|u'>|, the overlap of their columns. A pure state drops no resolved
   eigenvalue, so only the eigenvalues a support drops add to the bound. *)
formFidelity[rankOne[c_, k_], rankOne[c2_, k2_], _, _, p_] := {Abs[uniformAccuracy[Conjugate[c] . c2 / Sqrt[k k2], p]], 0}
formFidelity[rankOne[c_, k_], {b_, v_, eb_, vb_}, r_, _, p_] :=
    {If[b === {}, N[0, p], Norm[uniformAccuracy[Sqrt[b] (Conjugate[v] . c) / Sqrt[k], p]]], droppedBound[eb, vb, r]}
formFidelity[support_List, pure_rankOne, r_, s_, p_] := formFidelity[pure, support, s, r, p]
formFidelity[{a_, u_, ea_, ua_}, {b_, v_, eb_, vb_}, r_, s_, p_] :=
    {supportFidelity[{a, u}, {b, v}, p], droppedBound[ea, ua, s] + droppedBound[eb, vb, r]}

(* Dropping eigenvalues e_i with eigenvectors u_i from r changes its fidelity with s by at most
   Sum_i Sqrt[e_i <u_i|s|u_i>], the trace norm of (Sqrt[r] - Sqrt[r']) . Sqrt[s] bounded term by term. *)
droppedBound[{}, _, _] := 0
droppedBound[e_, u_, s_] := Total[Sqrt[e Ramp[Re[Total[Conjugate[u] (u . Transpose[Normal[s]]), {2}]]]]]

(* A matrix whose every entry off the diagonal is 0, exact or a zero known only to an accuracy, as 0``29.8;
   DiagonalMatrixQ with its default tolerance also takes a machine matrix with small nonzero entries off the
   diagonal for diagonal. *)
diagonalQ[m_] := DiagonalMatrixQ[m, Tolerance -> 0]

(* The sum of the singular values of M from the two supports, {eigenvalues, eigenvectors as rows, ...}; with
   no eigenvalue left in a support, which only a non-physical input gives, the fidelity is 0. *)
supportFidelity[{a_, u_, ___}, {b_, v_, ___}, p_] := If[a === {} || b === {},
    N[0, p],
    Total[SingularValueList[zeroUnresolvedParts[uniformAccuracy[KroneckerProduct[Sqrt[a], Sqrt[b]] Normal[Conjugate[u] . Transpose[v]], p]], Tolerance -> 0]]
]

(* Above machine precision, a sum that cancels can leave a number whose imaginary part lies far below its
   accuracy and has few digits, as Complex[0``30, -1.8`13.5*^-58], and SingularValueList and Norm then return
   a result with about as few: SingularValueList gives 0``0.7 for the singular value 1/4 of a 2 x 2 matrix
   holding one, and Norm[{0.5`30, Complex[0``30, -1.8`13.5*^-58]}] is 0.5`1.95. A zero known to an accuracy
   does no such harm. Every entry is given one accuracy, the least of theirs and that of the input, which
   turns such an imaginary part into a zero known to that accuracy and keeps every value larger than its
   uncertainty; the entries are the terms of the fidelity itself, divided by Sqrt[m_kk] for a pure state,
   so that this accuracy is theirs. That is enough for Norm and Abs. SingularValueList still collapses, or
   does not return, when a real or imaginary part lies within a few times its uncertainty without being 0,
   and so has less than one digit: every part within ten times the uncertainty 10^-a of the entries is made
   a zero known to the accuracy a - 1, whose uncertainty covers the value it replaces. *)
uniformAccuracy[x_, MachinePrecision] := x
uniformAccuracy[x_, p_] := SetAccuracy[x, Min[Accuracy[x], p]]

zeroUnresolvedParts[m_] /; Precision[m] === MachinePrecision := m
zeroUnresolvedParts[m_] := With[{a = Accuracy[m]},
    With[{zeroed = Function[x, If[Abs[x] <= 10 10^-a, SetAccuracy[0, a - 1], x]]},
        Map[zeroed[Re[#]] + I zeroed[Im[#]] &, m, {2}]
    ]
]

(* The positive eigenvalues of a diagonal matrix, its diagonal entries as stored, with the unit vectors as
   eigenvectors. No eigensolver adds round-off to them, so each is kept however small: the fidelity of
   diag(1 - 10^-20, 10^-20) with the pure state {10^-12, Sqrt[1 - 10^-24]} lies almost all in the
   population 10^-20. *)
hermitianSupport[m_] /; diagonalQ[m] := With[{l = Re[Normal[Diagonal[m]]]},
    With[{keep = Thread[l > 0]},
        {Pick[l, keep], IdentityMatrix[Length[l], SparseArray][[Pick[Range[Length[l]], keep]]], {}, {}}
    ]
]

(* {kept eigenvalues, their eigenvectors as rows, bounds on the resolved dropped eigenvalues, their
   eigenvectors as rows} of the Hermitian part of an inexact matrix at precision p. An eigenvalue is kept
   when it is above 100 10^-p times the largest in magnitude; the others, negative ones included, are
   dropped. Round-off moves the eigenvalues by about 10^-p times the largest, so the cut keeps the square
   roots of round-off of 0 out of the fidelity. A genuine eigenvalue below the cut is dropped too, which
   changes the fidelity by no more than the square root of the cut for each one dropped, and only by that
   much when the other state has its weight on the eigenvector. A dropped eigenvalue that is above d 10^-p,
   the bound on how far an error of 10^-p in each entry moves an eigenvalue, by more than its own
   uncertainty is resolved: it is kept aside with d 10^-p added, to bound the change it makes; any other is
   taken as 0. A diagonal matrix keeps
   every positive entry instead, so two states whose overlap lies in eigenvalues below the cut get a
   fidelity that depends, by up to that much, on whether the basis they are given in makes them diagonal.
   The eigenvectors are orthonormal. At machine precision they come from Eigensystem, whose documentation
   promises normalized, independent eigenvectors only; on a matrix that is Hermitian exactly, as the
   Hermitian part of a machine matrix is, it runs its Hermitian solver, which returns them orthonormal. At
   higher precision Eigensystem leaves two eigenvectors of a repeated eigenvalue orthogonal to far fewer
   digits than the working precision, and they are the columns of q in the Schur decomposition
   q . t . ConjugateTranspose[q] instead, where q is unitary and t, for a Hermitian matrix, is diagonal with
   the eigenvalues. *)
hermitianSupport[m_] := With[{e = hermitianEigensystem[m], p = workingPrecision[m]},
    With[{l = Re[e[[1]]], delta = N[Length[m] 10^-p]},
        With[{keep = Thread[l > N[100 10^-p] Max[Abs[e[[1]]]]]},
            With[{resolved = Boole[Thread[l > delta]] - Boole[keep]},
                {Pick[l, keep], Pick[e[[2]], keep], Pick[l, resolved, 1] + delta, Pick[e[[2]], resolved, 1]}
            ]
        ]
    ]
]

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

