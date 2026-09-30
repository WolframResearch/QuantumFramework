(* ::Package:: *)

BeginTestSection["QuantumDistance"]


(* ========== Setup ========== *)

qs0 = QuantumState[{1, 0}];
qs1 = QuantumState[{0, 1}];
qsPlus = QuantumState[{1, 1} / Sqrt[2]];
qsMinus = QuantumState[{1, -1} / Sqrt[2]];
qsMixed = QuantumState[{{1/2, 0}, {0, 1/2}}];


(* ========== QuantumDistance: basic behavior ========== *)

(* Default metric is Fidelity *)
VerificationTest[
    QuantumDistance[qs0, qs1],
    QuantumDistance[qs0, qs1, "Fidelity"],
    TestID -> "Distance-DefaultIsFidelity"
]

(* Distance to self is 0 for all metrics *)
VerificationTest[
    Chop /@ (QuantumDistance[qs0, qs0, #] & /@ {"Fidelity", "Trace", "Bures", "BuresAngle", "HilbertSchmidt", "Bloch"}),
    {0, 0, 0, 0, 0, 0},
    TestID -> "Distance-SelfIsZero"
]

(* Distance is non-negative for all metrics *)
VerificationTest[
    AllTrue[
        QuantumDistance[qs0, qs1, #] & /@ {"Fidelity", "Trace", "Bures", "BuresAngle", "HilbertSchmidt", "Bloch"},
        # >= 0 &
    ],
    True,
    TestID -> "Distance-NonNegative"
]


(* ========== Fidelity ========== *)

(* Orthogonal pure states -> max fidelity distance = 1 *)
VerificationTest[
    QuantumDistance[qs0, qs1, "Fidelity"],
    1,
    TestID -> "Fidelity-OrthogonalStates"
]

(* Same state -> fidelity distance = 0 *)
VerificationTest[
    QuantumDistance[qs0, qs0, "Fidelity"] // Chop,
    0,
    TestID -> "Fidelity-SameState"
]

(* Non-orthogonal states -> between 0 and 1 *)
VerificationTest[
    0 < QuantumDistance[qs0, qsPlus, "Fidelity"] < 1,
    True,
    TestID -> "Fidelity-NonOrthogonal"
]


(* ========== Trace distance ========== *)

(* Orthogonal pure states -> trace distance = 1 *)
VerificationTest[
    QuantumDistance[qs0, qs1, "Trace"],
    1,
    TestID -> "Trace-OrthogonalStates"
]

VerificationTest[
    QuantumDistance[qs0, qs0, "Trace"] // Chop,
    0,
    TestID -> "Trace-SameState"
]

(* Trace distance is bounded [0, 1] *)
VerificationTest[
    0 <= QuantumDistance[qs0, qsPlus, "Trace"] <= 1,
    True,
    TestID -> "Trace-Bounded"
]


(* ========== Bures distance ========== *)

(* Bures distance range is [0, Sqrt[2]] *)
VerificationTest[
    QuantumDistance[qs0, qs1, "Bures"],
    Sqrt[2],
    TestID -> "Bures-MaxDistance"
]

VerificationTest[
    QuantumDistance[qs0, qs0, "Bures"] // Chop,
    0,
    TestID -> "Bures-SameState"
]

(* Intermediate value *)
VerificationTest[
    0 < QuantumDistance[qs0, qsPlus, "Bures"] < Sqrt[2],
    True,
    TestID -> "Bures-Intermediate"
]


(* ========== BuresAngle ========== *)

(* BuresAngle range is [0, Pi/2] *)
VerificationTest[
    QuantumDistance[qs0, qs1, "BuresAngle"],
    Pi / 2,
    TestID -> "BuresAngle-MaxDistance"
]

VerificationTest[
    QuantumDistance[qs0, qs0, "BuresAngle"] // Chop,
    0,
    TestID -> "BuresAngle-SameState"
]


(* ========== HilbertSchmidt ========== *)

VerificationTest[
    QuantumDistance[qs0, qs1, "HilbertSchmidt"],
    Sqrt[2],
    TestID -> "HilbertSchmidt-OrthogonalStates"
]

VerificationTest[
    QuantumDistance[qs0, qs0, "HilbertSchmidt"] // Chop,
    0,
    TestID -> "HilbertSchmidt-SameState"
]


(* ========== Bloch distance ========== *)

(* Bloch distance between antipodal states on Bloch sphere *)
VerificationTest[
    QuantumDistance[qs0, qs1, "Bloch"],
    1,
    TestID -> "Bloch-AntipodalStates"
]

VerificationTest[
    QuantumDistance[qs0, qs0, "Bloch"] // Chop,
    0,
    TestID -> "Bloch-SameState"
]


(* ========== RelativePurity ========== *)

(* Tr[rho . rho] = 1 for pure state with itself *)
VerificationTest[
    QuantumDistance[qs0, qs0, "RelativePurity"],
    0,
    TestID -> "RelativePurity-SameState"
]

(* Orthogonal pure states -> Tr[rho . sigma] = 0 *)
VerificationTest[
    QuantumDistance[qs0, qs1, "RelativePurity"],
    1,
    TestID -> "RelativePurity-OrthogonalStates"
]

(* Non-orthogonal -> between 0 and 1 *)
VerificationTest[
    0 < QuantumDistance[qs0, qsPlus, "RelativePurity"] < 1,
    True,
    TestID -> "RelativePurity-Intermediate"
]


(* ========== RelativeEntropy ========== *)

(* Relative entropy DIVERGES against a pure state unless the states are equal:
   supp(s) must lie inside supp(t), and a pure t has a one-dimensional support. *)
VerificationTest[
    {
        QuantumDistance[qsMixed, qs0, "RelativeEntropy"],
        QuantumDistance[qs0, qs1, "RelativeEntropy"],
        QuantumDistance[qs0, qs0, "RelativeEntropy"]
    },
    {Quantity[Infinity, "Bits"], Quantity[Infinity, "Bits"], Quantity[0, "Bits"]},
    TestID -> "RelativeEntropy-ToPureState"
]

(* A PURE first argument against a full-rank second is finite: the zero
   eigenvalue of s contributes nothing under 0 log 0 = 0, which is why MatrixLog
   cannot be used here at all - Log is not analytic at 0.
   S(|0><0| || I/2) = log2(2) = 1 bit. *)
VerificationTest[
    Chop[QuantityMagnitude[QuantumDistance[qs0, qsMixed, "RelativeEntropy"]] - 1, 1*^-8],
    0,
    TestID -> "RelativeEntropy-PureAgainstFullRankIsFinite"
]

(* Exact input gives an exact answer, as it does for the other distances: the
   tolerance and the Re that an inexact spectrum needs would both spend
   exactness the input still has. *)
(* S(diag(3/4,1/4) || I/2) = 1 - H(3/4) against the maximally mixed state. *)
VerificationTest[
    With[{
        d = QuantumDistance[QuantumState[{{3/4, 0}, {0, 1/4}}], qsMixed, "RelativeEntropy"],
        ref = 1 + (3/4) Log[2, 3/4] + (1/4) Log[2, 1/4]
    },
        {
            Precision[d],
            Simplify[QuantityMagnitude[d] - ref],
            Precision @ QuantumDistance[qs0, qsMixed, "RelativeEntropy"]
        }
    ],
    {Infinity, 0, Infinity},
    TestID -> "RelativeEntropy-ExactInputStaysExact"
]

(* Machine input stays machine, and agrees with the exact answer. *)
VerificationTest[
    With[{
        d = QuantumDistance[QuantumState[N @ {{3/4, 0}, {0, 1/4}}], QuantumState[N @ {{1/2, 0}, {0, 1/2}}], "RelativeEntropy"],
        ref = N[1 + (3/4) Log[2, 3/4] + (1/4) Log[2, 1/4]]
    },
        {Precision[d], Chop[QuantityMagnitude[d] - ref, 1*^-10]}
    ],
    {MachinePrecision, 0},
    TestID -> "RelativeEntropy-MachineInputStaysMachine"
]

(* Eigensystem hands back an arbitrary basis inside a degenerate eigenspace, and
   for exact input that basis is generally not orthonormal, so Sum_j |t_j><t_j|
   fails to resolve the identity and the cross term stops telescoping to
   Tr[s log t]. Here t has eigenvalues {1/2, 1/4, 1/4} in a rotated basis: the
   exact path read ~0.07 bits high while the machine path was right to 1e-16,
   which is why only an EXACT degenerate second argument catches it. *)
VerificationTest[
    With[{
        t = With[{p = Outer[Times, {1, 1, 1} / Sqrt[3], {1, 1, 1} / Sqrt[3]]},
            p / 2 + (IdentityMatrix[3] - p) / 4],
        s = With[{a = {{2 + I, 1, 3}, {0, 1 - 2 I, 1}, {1, 2, 1 + I}}},
            # / Tr[#] & [a . ConjugateTranspose[a]]]
    },
        Chop[
            Re @ N @ QuantityMagnitude @ QuantumDistance[QuantumState[s], QuantumState[t], "RelativeEntropy"] -
                Re[Tr[N[s] . MatrixLog[N[s]]] - Tr[N[s] . MatrixLog[N[t]]]] / Log[2],
            1*^-8
        ]
    ],
    0,
    TestID -> "RelativeEntropy-DegenerateExactSecondArgument"
]

(* Umegaki value against a reference computed here from the definition, both
   directions, since the quantity is not symmetric. Both operands are full rank,
   which is the case where MatrixLog is legitimate, so the reference is
   independent of how QuantumDistance computes it. *)
VerificationTest[
    With[{am = {{0.7, 0.1}, {0.1, 0.3}}, bm = {{0.6, 0.}, {0., 0.4}}},
        With[{umegaki = Function[{x, y},
            Re[Tr[x . MatrixLog[x, Method -> "Jordan"]] - Tr[x . MatrixLog[y, Method -> "Jordan"]]] / Log[2]]},
            Chop[{
                QuantityMagnitude[QuantumDistance[QuantumState[am], QuantumState[bm], "RelativeEntropy"]] - umegaki[am, bm],
                QuantityMagnitude[QuantumDistance[QuantumState[bm], QuantumState[am], "RelativeEntropy"]] - umegaki[bm, am]
            }, 1*^-6]
        ]
    ],
    {0, 0},
    TestID -> "RelativeEntropy-matches-Umegaki-both-directions"
]

(* Infinite distance means zero similarity, which the exponential mapping in
   QuantumSimilarity gives without a special case. *)
VerificationTest[
    QuantumSimilarity[qs0, qs1, "RelativeEntropy"],
    0,
    TestID -> "Similarity-RelativeEntropy-OrthogonalPuresAreZero"
]


(* ========== QuantumSimilarity: metric-aware normalization ========== *)

(* Similarity to self should be 1 for bounded metrics *)
VerificationTest[
    DeleteDuplicates[Chop /@ (QuantumSimilarity[qs0, qs0, #] & /@ {"Fidelity", "Trace", "Bures", "BuresAngle", "HilbertSchmidt", "Bloch"})],
    {1},
    TestID -> "Similarity-SelfIsOne"
]

(* Similarity of maximally distant states should be 0 *)
VerificationTest[
    DeleteDuplicates[Chop /@ (QuantumSimilarity[qs0, qs1, #] & /@ {"Fidelity", "Trace", "Bures", "BuresAngle", "HilbertSchmidt", "Bloch"})],
    {0},
    TestID -> "Similarity-OrthogonalIsZero"
]

(* Similarity is in [0, 1] for intermediate states, all bounded metrics *)
VerificationTest[
    AllTrue[
        QuantumSimilarity[qs0, qsPlus, #] & /@ {"Fidelity", "Trace", "Bures", "BuresAngle", "HilbertSchmidt", "Bloch"},
        0 < # < 1 &
    ],
    True,
    TestID -> "Similarity-IntermediateInUnitInterval"
]

(* RelativePurity similarity: same as Tr[rho . sigma] directly *)
VerificationTest[
    QuantumSimilarity[qs0, qs0, "RelativePurity"],
    1,
    TestID -> "Similarity-RelativePurity-SelfIsOne"
]

VerificationTest[
    QuantumSimilarity[qs0, qs1, "RelativePurity"],
    0,
    TestID -> "Similarity-RelativePurity-OrthogonalIsZero"
]

(* RelativeEntropy similarity: pure state self -> 2^0 = 1 *)
VerificationTest[
    QuantumSimilarity[qs0, qs0, "RelativeEntropy"],
    1,
    TestID -> "Similarity-RelativeEntropy-PureSelf"
]


(* ========== Symmetry: d(a,b) = d(b,a) for symmetric metrics ========== *)

VerificationTest[
    AllTrue[
        {"Fidelity", "Trace", "Bures", "BuresAngle", "HilbertSchmidt", "Bloch"},
        Chop[QuantumDistance[qs0, qsPlus, #] - QuantumDistance[qsPlus, qs0, #]] == 0 &
    ],
    True,
    TestID -> "Distance-Symmetry"
]


(* ========== Triangle inequality for true metrics ========== *)

VerificationTest[
    AllTrue[
        {"Trace", "Bures", "HilbertSchmidt", "Bloch"},
        QuantumDistance[qs0, qs1, #] <=
            QuantumDistance[qs0, qsPlus, #] + QuantumDistance[qsPlus, qs1, #] + 10^-10 &
    ],
    True,
    TestID -> "Distance-TriangleInequality"
]


(* ========== Mixed state tests ========== *)

(* Mixed state distances should be numeric *)
VerificationTest[
    AllTrue[
        QuantumDistance[qsMixed, qs0, #] & /@ {"Fidelity", "Trace", "Bures", "BuresAngle", "HilbertSchmidt"},
        NumericQ
    ],
    True,
    TestID -> "Distance-MixedStateNumeric"
]

(* Maximally mixed state is equidistant from |0> and |1> *)
VerificationTest[
    AllTrue[
        {"Fidelity", "Trace", "Bures", "BuresAngle", "HilbertSchmidt"},
        Chop[QuantumDistance[qsMixed, qs0, #] - QuantumDistance[qsMixed, qs1, #]] == 0 &
    ],
    True,
    TestID -> "Distance-MixedEquidistant"
]


(* ========== Non-physical input warning ========== *)

(* a Hermitian matrix with a negative eigenvalue is not a physical state; the distance still computes but
   warns *)
VerificationTest[
    QuantumDistance[QuantumState[{{0.5, 0.6}, {0.6, 0.5}}], QuantumState[{{0.5, 0.6}, {0.6, 0.5}}], "Fidelity"],
    _ ? NumericQ,
    {QuantumDistance::notphysical},
    SameTest -> MatchQ,
    TestID -> "Distance-NonPSD-Warns"
]

(* a genuine physical state does not warn *)
VerificationTest[
    QuantumDistance[qs0, qsMixed, "Fidelity"],
    _ ? NumericQ,
    {},
    SameTest -> MatchQ,
    TestID -> "Distance-PSD-NoWarn"
]

(* A matrix that is not Hermitian is not a state even when its eigenvalues are those of a state:
   {{1/2, 1}, {0, 1/2}} has the eigenvalues of I/2, and warns. *)
VerificationTest[
    QuantumDistance[QuantumState[{{1/2, 1}, {0, 1/2}}], qsMixed, "Trace"],
    1/2,
    {QuantumDistance::notphysical},
    TestID -> "NonHermitian-NotPhysical"
]

(* A matrix of exact rational entries is checked for Hermiticity exactly, so an asymmetry of 10^-12, below
   any machine tolerance, is seen. *)
VerificationTest[
    QuantumDistance[QuantumState[{{1/2, 10^-12}, {0, 1/2}}], qsMixed, "Trace"],
    1/2000000000000,
    {QuantumDistance::notphysical},
    TestID -> "NonHermitian-Exact-NotPhysical"
]

(* A machine matrix is Hermitian up to a tolerance of 10^-8 on its entries: an asymmetry of 10^-6 warns, one
   of 10^-10 is round-off. *)
VerificationTest[
    {
        QuantumDistance[QuantumState[{{0.5, 0.2 + 1.*^-6}, {0.2, 0.5}}], qsMixed, "Trace"],
        QuantumDistance[QuantumState[{{0.5, 0.2 + 1.*^-10}, {0.2, 0.5}}], qsMixed, "Trace"]
    },
    {_Real, _Real},
    {QuantumDistance::notphysical},
    SameTest -> MatchQ,
    TestID -> "NonHermitian-Machine-Tolerance"
]

(* So is an asymmetry of 10^-17 opposite a structural zero, which QuantumState stores in a SparseArray. *)
VerificationTest[
    QuantumDistance[QuantumState[{{0.5, 1.*^-17}, {0., 0.5}}], qsMixed, "Trace"],
    _Real,
    {},
    SameTest -> MatchQ,
    TestID -> "NonHermitian-Machine-SparseRoundoff"
]

(* A zero state has no trace to divide by; it is not physical, and the measure is computed as given. *)
VerificationTest[
    QuantumDistance[QuantumState[{{0, 0}, {0, 0}}], qs0, "Trace"],
    1/2,
    {QuantumDistance::notphysical},
    TestID -> "ZeroTrace-NotPhysical"
]

(* Nor is a trace that is negative or not real a state's, at any scale: -10^-9 diag(1/4, 3/4) has no
   eigenvalue below -10^-8, and every entry of 10^-9 diag(1/4 + I, 3/4) lies below the 10^-8 Hermiticity
   tolerance, yet neither is physical. Neither is divided by its trace, so the measure is computed as given. *)
VerificationTest[
    {
        QuantumDistance[QuantumState[-10^-9 {{1/4, 0}, {0, 3/4}}], qsMixed, "Trace"],
        QuantumDistance[QuantumState[1.*^-9 {{0.25 + I, 0.}, {0., 0.75}}], qsMixed, "Trace"]
    },
    {1000000001/2000000000, _Real ? (Abs[# - 0.4999999995] < 1.*^-12 &)},
    {QuantumDistance::notphysical, QuantumDistance::notphysical},
    SameTest -> MatchQ,
    TestID -> "NegativeOrComplexTrace-NotPhysical"
]


(* ========== Input whose trace is not 1 ========== *)

(* An input whose trace (for a state vector, its squared norm) is a positive number other than 1 is
   divided by its trace before any measure is computed, with a notnormalized message naming it: 2|0> is
   read as |0>, and the identity matrix, with trace 2, as the maximally mixed state. *)
VerificationTest[
    QuantumSimilarity[qs0, QuantumState[{2, 0}]],
    1,
    {QuantumDistance::notnormalized},
    TestID -> "Unnormalized-Ket-Rescaled"
]

VerificationTest[
    QuantumDistance[QuantumState[{{1/4, 0}, {0, 3/4}}], QuantumState[{{1, 0}, {0, 1}}]] ===
        QuantumDistance[QuantumState[{{1/4, 0}, {0, 3/4}}], qsMixed],
    True,
    {QuantumDistance::notnormalized},
    TestID -> "Unnormalized-Matrix-Rescaled"
]

(* Every measure gives 2|+> and 3 diag(1/4, 3/4) the value it gives |+> and diag(1/4, 3/4), a pair at
   which each distance lies inside its range, and in this order the relative entropy is finite. Each
   input is rescaled, so each call issues two messages. *)
Scan[
    Function[measure,
        VerificationTest[
            QuantumDistance[QuantumState[{Sqrt[2], Sqrt[2]}], QuantumState[{{3/4, 0}, {0, 9/4}}], measure],
            QuantumDistance[qsPlus, QuantumState[{{1/4, 0}, {0, 3/4}}], measure],
            {QuantumDistance::notnormalized, QuantumDistance::notnormalized},
            TestID -> "Unnormalized-ScaleInvariant-" <> measure
        ]
    ],
    {"Fidelity", "RelativeEntropy", "RelativePurity", "Trace", "Bures", "BuresAngle", "HilbertSchmidt", "Bloch"}
]

(* The same for a qutrit pair, 5 diag(1/2, 1/3, 1/6) and 2 (1, 1, 1) against diag(1/2, 1/3, 1/6) and
   (1, 1, 1) / Sqrt[3], for every measure defined beyond a qubit. *)
Scan[
    Function[measure,
        VerificationTest[
            QuantumDistance[QuantumState[{{5/2, 0, 0}, {0, 5/3, 0}, {0, 0, 5/6}}], QuantumState[{2, 2, 2}, 3], measure],
            QuantumDistance[QuantumState[{{1/2, 0, 0}, {0, 1/3, 0}, {0, 0, 1/6}}], QuantumState[{1, 1, 1} / Sqrt[3], 3], measure],
            {QuantumDistance::notnormalized, QuantumDistance::notnormalized},
            TestID -> "Unnormalized-Qutrit-ScaleInvariant-" <> measure
        ]
    ],
    {"Fidelity", "RelativeEntropy", "RelativePurity", "Trace", "Bures", "BuresAngle", "HilbertSchmidt"}
]

(* A two-qubit state with a symbolic parameter: QuantumState["Werner"[w]] puts weight 1 - w on the singlet
   and w/3 on each triplet state, |Phi+> among them, so it is diagonal in the Bell basis with spectrum
   {1 - w, w/3, w/3, w/3}. Scaled by 3 and rescaled, its fidelity distance to |Phi+> is 1 - Sqrt[w/3], and
   its trace and Hilbert-Schmidt distances follow from the spectrum {w/3 - 1, w/3, w/3, 1 - w} of the
   difference with |Phi+><Phi+|. *)
Scan[
    Function[rule,
        VerificationTest[
            FullSimplify[
                QuantumDistance[
                    QuantumState[3 Normal[QuantumState["Werner"[w]]["DensityMatrix"]], QuantumBasis[{2, 2}]],
                    QuantumState["PhiPlus"],
                    First[rule]
                ] - Last[rule],
                0 < w < 1
            ],
            0,
            {QuantumDistance::notnormalized},
            TestID -> "Unnormalized-Werner-" <> First[rule]
        ]
    ],
    {"Fidelity" -> 1 - Sqrt[w / 3], "Trace" -> 1 - w / 3, "HilbertSchmidt" -> Sqrt[2 + 4 w (w - 2) / 3]}
]

(* However small, a positive trace is a state's: 10^-5 |0> and a matrix with trace 10^-9 are rescaled. These
   exact traces are classified exactly. *)
VerificationTest[
    {
        QuantumDistance[QuantumState[{10^-5, 0}], qs0],
        QuantumDistance[QuantumState[{{10^-9, 0}, {0, 0}}], qs0, "Trace"]
    },
    {0, 0},
    {QuantumDistance::notnormalized, QuantumDistance::notnormalized},
    TestID -> "Unnormalized-TinyTrace-Rescaled"
]

(* An algebraic trace whose sign Positive cannot decide is reduced with RootReduce, and once found positive
   is divided out like any other: the ket {1 + I Sqrt[2], 1} has the squared norm
   (1 + I Sqrt[2]) (1 - I Sqrt[2]) + 1, which is 4, and is read as its normalized form. *)
VerificationTest[
    QuantumDistance[QuantumState[{1 + I Sqrt[2], 1}], qsPlus, "Trace"],
    1/2,
    {QuantumDistance::notnormalized},
    TestID -> "Unnormalized-UndecidedTrace-Rescaled"
]

(* A machine trace is judged against the Frobenius norm of its matrix: the machine ket 10^-5 |0> is rescaled
   too, while the zero machine ket, and a machine matrix whose trace is 10^-9 against a norm near 1, have a
   zero trace and are not physical. *)
VerificationTest[
    Chop[
        {
            QuantumDistance[QuantumState[{1.*^-5, 0.}], qs0, "Trace"],
            QuantumDistance[QuantumState[{0., 0.}], qs0, "Trace"],
            QuantumDistance[QuantumState[{{0.5, 0.}, {0., -0.5 + 1.*^-9}}], qs0, "Trace"]
        } - {0, 1/2, 1/2},
        1.*^-8
    ],
    {0, 0, 0},
    {QuantumDistance::notnormalized, QuantumDistance::notphysical, QuantumDistance::notphysical},
    TestID -> "MachineTrace-JudgedAgainstNorm"
]

(* An exact trace is classified exactly, however small: {{p, 0}, {0, 1 - p}} / 10^9 has trace 10^-9 and is
   read as {{p, 0}, {0, 1 - p}}. *)
VerificationTest[
    QuantumDistance[QuantumState[{{p / 10^9, 0}, {0, 1 / 10^9 - p / 10^9}}], qs0],
    1 - Re[Sqrt[p]],
    {QuantumDistance::notnormalized},
    TestID -> "Symbolic-TinyExactTrace-Rescaled"
]

(* An exact trace written as a sum of nested radicals is reduced before it is compared: the ket
   (Sqrt[3 + 2 Sqrt[2]] - Sqrt[2]) |0> has norm exactly 1 and issues no message, and a matrix whose trace is
   Sqrt[2 + 10^-200] - Sqrt[2], positive but beyond machine resolution, is rescaled to |0><0|. *)
VerificationTest[
    {
        QuantumDistance[QuantumState[{Sqrt[3 + 2 Sqrt[2]] - Sqrt[2], 0}], qs0, "Trace"],
        QuantumDistance[QuantumState[{{Sqrt[2 + 10^-200] - Sqrt[2], 0}, {0, 0}}], qs0, "Trace"]
    },
    {0, 0},
    {QuantumDistance::notnormalized},
    TestID -> "ExactTrace-NestedRadicals"
]

(* An exact trace is compared with 1 exactly: the squared norm of {1, 7/10^5} is 1 + 49/10^10, so the ket
   is rescaled, and its similarity to itself is 1, not a number above 1. *)
VerificationTest[
    With[{v = QuantumState[{1, 7 / 10^5}]}, QuantumSimilarity[v, v, "RelativePurity"]],
    1,
    {QuantumDistance::notnormalized, QuantumDistance::notnormalized},
    TestID -> "Unnormalized-ExactNearUnitTrace-Rescaled"
]

(* A machine trace is always divided, since round-off leaves it a few units in the last place from 1, and
   within 10^-8 of 1 without a message. *)
VerificationTest[
    Chop[QuantumDistance[QuantumState[{{0.25 + 1.*^-12, 0.}, {0., 0.75}}], qsMixed, "Trace"] - 0.25, 1.*^-10],
    0,
    {},
    TestID -> "Normalized-MachineRounding-NoMessage"
]

(* Divided by its machine trace, the same ket is |0> up to a unit in the last place, so its fidelity with
   |0> can land a unit in the last place below 1, which the Bures distance turns into Sqrt[$MachineEpsilon];
   the distance is bounded by that of a fidelity two units in the last place below 1. *)
VerificationTest[
    With[{v = QuantumState[{Sqrt[1 - 9.*^-9], 0.}]}, {QuantumDistance[v, v, "Bures"], QuantumDistance[v, qs0, "Bures"] <= Sqrt[2 $MachineEpsilon]}],
    {0., True},
    {},
    TestID -> "Normalized-MachineNearUnitTrace-NoMessage"
]

(* A symbolic input is rescaled only when its trace evaluates to a number. {{p, 0}, {0, 1 - p}} has trace 1
   and is used as given; {{2 p, 0}, {0, 2 - 2 p}} has trace 2 and is read as the same state. *)
VerificationTest[
    QuantumDistance[QuantumState[{{p, 0}, {0, 1 - p}}], qs0],
    1 - Re[Sqrt[p]],
    {},
    TestID -> "Symbolic-UnitTrace-NoMessage"
]

VerificationTest[
    QuantumDistance[QuantumState[{{2 p, 0}, {0, 2 - 2 p}}], qs0],
    1 - Re[Sqrt[p]],
    {QuantumDistance::notnormalized},
    TestID -> "Symbolic-NumericTrace-Rescaled"
]

(* The trace of a ket with symbolic amplitudes is not a number, and the ket is used as given: {1, x} is read
   through its unnormalized density matrix, at distance 0 from |0> for every x, while its "Normalized" state
   gives 1 - 1/Sqrt[1 + Abs[x]^2]. A ket of unit norm for real parameters needs no rescaling. *)
VerificationTest[
    {
        QuantumDistance[QuantumState[{1, x}], qs0],
        FullSimplify[QuantumDistance[QuantumState[{1, x}]["Normalized"], qs0] - (1 - 1 / Sqrt[1 + Abs[x]^2])]
    },
    {0, 0},
    {},
    TestID -> "Symbolic-Ket-UsedAsGiven"
]

VerificationTest[
    FullSimplify[
        QuantumDistance[QuantumState[{Cos[\[Theta] / 2], E^(I \[Phi]) Sin[\[Theta] / 2]}], qsPlus],
        {\[Theta], \[Phi]} \[Element] Reals
    ],
    1 - Sqrt[2 + 2 Cos[\[Phi]] Sin[\[Theta]]] / 2,
    {},
    TestID -> "Symbolic-Ket-NoMessage"
]


(* ========== Equal states are at distance 0 ========== *)

randomDensityMatrix[d_] := With[{a = RandomComplex[{-1 - I, 1 + I}, {d, d}]}, # / Re[Tr[#]] &[a . ConjugateTranspose[a]]]

(* Identical machine states are at fidelity, Bures and BuresAngle distance 0., whichever side of 0 the
   round-off of the formula would have landed on: twenty random qutrit states against themselves. *)
VerificationTest[
    BlockRandom[
        Union @ Flatten @ Table[
            With[{r = QuantumState[randomDensityMatrix[3]]}, QuantumDistance[r, r, #] & /@ {"Fidelity", "Bures", "BuresAngle"}],
            20
        ],
        RandomSeeding -> 3
    ],
    {0.},
    TestID -> "SameMachineState-ExactlyZero"
]

(* Two machine states equal up to round-off, one rotated by a random unitary and back: the fidelity distance
   of the pair lands a few units in the last place from 0 on either side, and the Bures distance is a real
   number of the order of the square root of the machine epsilon, never an imaginary one. *)
VerificationTest[
    BlockRandom[
        With[{bures = Table[
            With[{m = randomDensityMatrix[3], u = RandomVariate[CircularUnitaryMatrixDistribution[3]]},
                QuantumDistance[QuantumState[m], QuantumState[ConjugateTranspose[u] . (u . m . ConjugateTranspose[u]) . u], "Bures"]
            ],
            20
        ]},
            {Union[Head /@ bures], Max[bures] < 1.*^-6}
        ],
        RandomSeeding -> 4
    ],
    {{Real}, True},
    TestID -> "NearEqualMachineStates-BuresReal"
]

(* An exact state against itself, and two equal exact matrices written differently, Sqrt[3 + 2 Sqrt[2]]
   against 1 + Sqrt[2], which the kernel does not denest: RootReduce finds every entry of their difference
   to be 0, so every measure gives the distance of a state to itself, 0 or 0 bits, and "RelativePurity"
   gives 1 - Tr[rho^2], with no message. From the formulas the distances would be sums of nested radicals
   equal to 0, on which Re and ArcCos cannot decide. *)
VerificationTest[
    With[{
        e = QuantumState[{{1/4, 1/4}, {1/4, 3/4}}],
        a = QuantumState[{{1/4, (1 + Sqrt[2]) / 8}, {(1 + Sqrt[2]) / 8, 3/4}}],
        b = QuantumState[{{1/4, Sqrt[3 + 2 Sqrt[2]] / 8}, {Sqrt[3 + 2 Sqrt[2]] / 8, 3/4}}]
    },
        {
            QuantumDistance[e, e, #] & /@ {"Fidelity", "Bures", "BuresAngle"},
            QuantumDistance[a, b, #] & /@ {"Fidelity", "RelativeEntropy", "Trace", "Bures", "BuresAngle", "HilbertSchmidt", "Bloch"},
            RootReduce[QuantumDistance[a, b, "RelativePurity"] - (1 - Tr[MatrixPower[Normal[a["DensityMatrix"]], 2]])],
            QuantumSimilarity[a, b, #] & /@ {"Fidelity", "Trace", "Bures", "BuresAngle", "HilbertSchmidt", "Bloch"}
        }
    ],
    {{0, 0, 0}, {0, Quantity[0, "Bits"], 0, 0, 0, 0, 0}, 0, {1, 1, 1, 1, 1, 1}},
    {},
    TestID -> "EqualExactStates-ExactlyZero"
]


(* ========== Near-equal states ========== *)

(* Two pure states at angle d on the Bloch sphere: the BuresAngle distance is d/2, the Bures distance
   2 Sin[d/4] and the fidelity distance 2 Sin[d/4]^2, all going to 0 with d, where equal states give 0.
   In machine precision, for the same pair in a random frame, the fidelity distance is of order d^2 and
   carries an absolute error of order the machine epsilon, so the BuresAngle distance has a relative error
   of order $MachineEpsilon / d^2: a few 10^-9 at d = 10^-3, and no correct digit once d is near 10^-8. *)
VerificationTest[
    {
        FullSimplify[QuantumDistance[qs0, QuantumState[{Cos[d / 2], Sin[d / 2]}], #] & /@ {"BuresAngle", "Bures", "Fidelity"}, 0 < d < Pi],
        BlockRandom[
            With[{u = RandomVariate[CircularUnitaryMatrixDistribution[2]]},
                Function[dd,
                    Abs[QuantumDistance[QuantumState[u . {1., 0.}], QuantumState[u . {Cos[dd / 2], Sin[dd / 2]}], "BuresAngle"] / (dd / 2) - 1] <
                        100 $MachineEpsilon / dd^2
                ] /@ {1.*^-3, 1.*^-4}
            ],
            RandomSeeding -> 7
        ]
    },
    {{d / 2, 2 Sin[d / 4], 2 Sin[d / 4]^2}, {True, True}},
    TestID -> "SmallAngle-BuresAngleIsHalfTheAngle"
]

(* Commuting states are classical: for diag(p, 1 - p) and diag(q, 1 - q) every measure reduces to its
   classical counterpart on the distributions {p, 1 - p} and {q, 1 - q}, the Bhattacharyya coefficient for
   the fidelity family, |p - q| for "Trace" and "Bloch", and the Kullback-Leibler divergence in bits for
   the relative entropy. *)
VerificationTest[
    With[{bhattacharyya = Sqrt[p q] + Sqrt[(1 - p) (1 - q)]},
        FullSimplify[
            Replace[
                QuantumDistance[QuantumState[{{p, 0}, {0, 1 - p}}], QuantumState[{{q, 0}, {0, 1 - q}}], First[#]],
                d_Quantity :> QuantityMagnitude[d]
            ] - Last[#],
            0 < p < 1 && 0 < q < 1
        ] & /@ {
            "Fidelity" -> 1 - bhattacharyya,
            "RelativeEntropy" -> (p Log[p / q] + (1 - p) Log[(1 - p) / (1 - q)]) / Log[2],
            "RelativePurity" -> 1 - p q - (1 - p) (1 - q),
            "Trace" -> Abs[p - q],
            "Bures" -> Sqrt[2 (1 - bhattacharyya)],
            "BuresAngle" -> ArcCos[bhattacharyya],
            "HilbertSchmidt" -> Sqrt[2] Abs[p - q],
            "Bloch" -> Abs[p - q]
        }
    ],
    ConstantArray[0, 8],
    TestID -> "ClassicalLimit-CommutingStates"
]


(* ========== Independent references ========== *)

(* The fidelity of random qutrit states against the nuclear norm of Sqrt[rho] . Sqrt[sigma], computed
   without QuantumDistance, and the Fuchs-van de Graaf inequalities 1 - F <= T <= Sqrt[1 - F^2] between the
   fidelity F and the trace distance T. *)
VerificationTest[
    BlockRandom[
        Union @ Flatten @ Table[
            With[{r = randomDensityMatrix[3], s = randomDensityMatrix[3]},
                With[{
                    f = 1 - QuantumDistance[QuantumState[r], QuantumState[s]],
                    t = QuantumDistance[QuantumState[r], QuantumState[s], "Trace"]
                },
                    {
                        Abs[f - Total[SingularValueList[MatrixPower[r, 1/2] . MatrixPower[s, 1/2]]]] < 1.*^-10,
                        1 - f <= t + 1.*^-12,
                        t <= Sqrt[1 - f^2] + 1.*^-12
                    }
                ]
            ],
            10
        ],
        RandomSeeding -> 5
    ],
    {True},
    TestID -> "Fidelity-NuclearNorm-FuchsVanDeGraaf"
]

(* Every measure is invariant under a unitary change of frame applied to both states, and all but the
   relative entropy are symmetric in the two states: random mixed qutrit pairs that do not commute. *)
VerificationTest[
    BlockRandom[
        Max @ Flatten @ Table[
            With[{r = randomDensityMatrix[3], s = randomDensityMatrix[3], u = RandomVariate[CircularUnitaryMatrixDistribution[3]]},
                {
                    Abs[QuantumDistance[QuantumState[r], QuantumState[s], #] - QuantumDistance[QuantumState[s], QuantumState[r], #]] & /@
                        {"Fidelity", "RelativePurity", "Trace", "Bures", "BuresAngle", "HilbertSchmidt"},
                    Abs[QuantityMagnitude[QuantumDistance[QuantumState[u . r . ConjugateTranspose[u]], QuantumState[u . s . ConjugateTranspose[u]], #]] -
                        QuantityMagnitude[QuantumDistance[QuantumState[r], QuantumState[s], #]]] & /@
                        {"Fidelity", "RelativeEntropy", "RelativePurity", "Trace", "Bures", "BuresAngle", "HilbertSchmidt"}
                }
            ],
            10
        ],
        RandomSeeding -> 8
    ] < 1.*^-12,
    True,
    TestID -> "Symmetry-UnitaryInvariance"
]

(* A channel cannot make two states more distinguishable: tracing out one qubit of random two-qubit mixed
   states never increases the fidelity, relative entropy, trace, Bures and BuresAngle distances. (The
   Hilbert-Schmidt distance is not contractive in general, and is left out.) *)
VerificationTest[
    BlockRandom[
        Max @ Flatten @ Table[
            With[{x = QuantumState[randomDensityMatrix[4], QuantumBasis[{2, 2}]], y = QuantumState[randomDensityMatrix[4], QuantumBasis[{2, 2}]]},
                Replace[
                    QuantumDistance[QuantumPartialTrace[x, {2}], QuantumPartialTrace[y, {2}], #] - QuantumDistance[x, y, #],
                    d_Quantity :> QuantityMagnitude[d]
                ] & /@ {"Fidelity", "RelativeEntropy", "Trace", "Bures", "BuresAngle"}
            ],
            20
        ],
        RandomSeeding -> 10
    ] <= 1.*^-12,
    True,
    TestID -> "Contraction-UnderPartialTrace"
]

(* The trace, Bures, BuresAngle and Hilbert-Schmidt distances are metrics: the triangle inequality holds on
   random mixed qutrit triples. The relative entropy bounds the trace distance T by Pinsker's inequality,
   S >= 2 T^2 / Log[2] in bits. *)
VerificationTest[
    BlockRandom[
        {
            Max @ Flatten @ Table[
                With[{x = QuantumState[randomDensityMatrix[3]], y = QuantumState[randomDensityMatrix[3]], z = QuantumState[randomDensityMatrix[3]]},
                    QuantumDistance[x, z, #] - QuantumDistance[x, y, #] - QuantumDistance[y, z, #] & /@ {"Trace", "Bures", "BuresAngle", "HilbertSchmidt"}
                ],
                30
            ] <= 1.*^-12,
            Min @ Table[
                With[{x = QuantumState[randomDensityMatrix[3]], y = QuantumState[randomDensityMatrix[3]]},
                    QuantityMagnitude[QuantumDistance[x, y, "RelativeEntropy"]] - 2 QuantumDistance[x, y, "Trace"]^2 / Log[2]
                ],
                30
            ] >= 0
        },
        RandomSeeding -> 11
    ],
    {True, True},
    TestID -> "Triangle-Pinsker"
]

(* A symbolic pair that does not commute, a qubit state on the z axis against one tilted by t toward x:
   the root fidelity F of two qubit states satisfies F^2 = Tr[rho . sigma] + 2 Sqrt[Det[rho] Det[sigma]]. *)
VerificationTest[
    With[{
        r = (IdentityMatrix[2] + a PauliMatrix[3]) / 2,
        s = (IdentityMatrix[2] + b (Cos[t] PauliMatrix[3] + Sin[t] PauliMatrix[1])) / 2
    },
        FullSimplify[
            (1 - QuantumDistance[QuantumState[r], QuantumState[s]])^2 - (Tr[r . s] + 2 Sqrt[Det[r] Det[s]]),
            0 < a < 1 && 0 < b < 1 && 0 < t < Pi
        ]
    ],
    0,
    TestID -> "Symbolic-NonCommutingQubits-FidelityFormula"
]

(* The same pair in the other measures, with Bloch vectors {0, 0, a} and b {Sin[t], 0, Cos[t]}: the relative
   entropy from the two spectra and the overlap Tr[rho . Log[sigma]], the relative purity from
   Tr[rho . sigma] = (1 + a b Cos[t])/2, and the Hilbert-Schmidt and Bloch distances from the distance
   between the Bloch vectors. *)
VerificationTest[
    With[{
        r = (IdentityMatrix[2] + a PauliMatrix[3]) / 2,
        s = (IdentityMatrix[2] + b (Cos[t] PauliMatrix[3] + Sin[t] PauliMatrix[1])) / 2,
        blochDistance = Sqrt[a^2 + b^2 - 2 a b Cos[t]]
    },
        FullSimplify[
            Replace[QuantumDistance[QuantumState[r], QuantumState[s], First[#]], d_Quantity :> QuantityMagnitude[d]] - Last[#],
            0 < a < 1 && 0 < b < 1 && 0 < t < Pi
        ] & /@ {
            "RelativeEntropy" -> (
                (1 + a) / 2 Log[(1 + a) / 2] + (1 - a) / 2 Log[(1 - a) / 2] - Log[(1 - b^2) / 4] / 2 - a Cos[t] Log[(1 + b) / (1 - b)] / 2
            ) / Log[2],
            "RelativePurity" -> 1 - (1 + a b Cos[t]) / 2,
            "HilbertSchmidt" -> blochDistance / Sqrt[2],
            "Bloch" -> blochDistance / 2
        }
    ],
    {0, 0, 0, 0},
    TestID -> "Symbolic-NonCommutingQubits-ClosedForms"
]

(* The Werner state QuantumState["Werner"[p, d]] of two qudits, with weight p on the symmetric subspace,
   against the maximally entangled state |Phi+> = Sum_i |ii> / Sqrt[d], which lies in that subspace: its
   overlap is 2 p / (d (d + 1)), so the fidelity distance is 1 - Sqrt[2 p / (d (d + 1))] and the trace
   distance 1 - 2 p / (d (d + 1)), for symbolic p and d = 2, 3, 4. *)
VerificationTest[
    Table[
        With[{w = QuantumState["Werner"[p, dim]], phi = QuantumState[Flatten[IdentityMatrix[dim]] / Sqrt[dim], QuantumBasis[{dim, dim}]]},
            FullSimplify[
                {QuantumDistance[w, phi] - (1 - Sqrt[2 p / (dim (dim + 1))]), QuantumDistance[w, phi, "Trace"] - (1 - 2 p / (dim (dim + 1)))},
                0 < p < 1
            ]
        ],
        {dim, 2, 4}
    ],
    ConstantArray[{0, 0}, 3],
    TestID -> "Symbolic-WernerQudits-ClosedForms"
]

(* For a qubit, half the Euclidean distance between Bloch vectors is the trace distance, and both are read
   from the same matrices. *)
VerificationTest[
    BlockRandom[
        Max @ Table[
            With[{a = QuantumState[randomDensityMatrix[2]], b = QuantumState[randomDensityMatrix[2]]},
                Abs[QuantumDistance[a, b, "Trace"] - QuantumDistance[a, b, "Bloch"]]
            ],
            10
        ] < 1.*^-12,
        RandomSeeding -> 6
    ],
    True,
    TestID -> "Qubit-TraceEqualsBloch"
]

(* A machine similarity stays in [0, 1] at round-off: the self RelativePurity similarity of a random pure
   state, which the formula can put a unit in the last place above 1, and the Trace and HilbertSchmidt
   similarities of orthogonal random states, which it can put just below 0. The differences are compared
   with 0, since a comparison of two machine numbers ignores the last bits. *)
VerificationTest[
    BlockRandom[
        With[{
            self = Table[
                With[{k = Normalize[RandomComplex[{-1 - I, 1 + I}, 4]]}, QuantumSimilarity[QuantumState[k], QuantumState[k], "RelativePurity"]],
                50
            ],
            orthogonal = Flatten @ Table[
                With[{basis = Orthogonalize[RandomComplex[{-1 - I, 1 + I}, {2, 4}]]},
                    QuantumSimilarity[QuantumState[First[basis]], QuantumState[Last[basis]], #] & /@ {"Trace", "HilbertSchmidt"}
                ],
                50
            ]
        },
            {Max[self] - 1 <= 0, Min[orthogonal] >= 0}
        ],
        RandomSeeding -> 9
    ],
    {True, True},
    TestID -> "Similarity-MachineRoundoff-InUnitInterval"
]


(* ========== Measure names ========== *)

(* A name that is not a measure leaves QuantumDistance and QuantumSimilarity unevaluated. *)
VerificationTest[
    {Head[QuantumDistance[qs0, qsPlus, "NoSuchMeasure"]], Head[QuantumSimilarity[qs0, qsPlus, "NoSuchMeasure"]]},
    {QuantumDistance, QuantumSimilarity},
    TestID -> "UnknownMeasure-Unevaluated"
]


(* ========== "Bloch" is defined for a single qubit ========== *)

(* Any other dimension, a qutrit or two qubits, gives a Failure naming the dimension, with no message, and
   QuantumSimilarity passes the Failure through. *)
VerificationTest[
    {
        QuantumDistance[QuantumState[{1, 0, 0}, 3], QuantumState[{0, 1, 0}, 3], "Bloch"],
        QuantumDistance[QuantumState["00"], QuantumState["11"], "Bloch"],
        QuantumSimilarity[QuantumState[{1, 0, 0}, 3], QuantumState[{0, 1, 0}, 3], "Bloch"]
    },
    {
        Failure["NonQubitBloch", KeyValuePattern["Dimension" -> 3]],
        Failure["NonQubitBloch", KeyValuePattern["Dimension" -> 4]],
        Failure["NonQubitBloch", KeyValuePattern["Dimension" -> 3]]
    },
    {},
    SameTest -> MatchQ,
    TestID -> "Bloch-NonQubit-Failure"
]


(* ========== Fidelity and trace distance from spectra ========== *)

(* A pair of exact mixed states of dimension 4 in general position: the fidelity is a sum of the square
   roots of the eigenvalues of r . s, which are Root objects, and the trace distance half a sum of
   singular values, so each result stays under a thousand leaves, the bound checked here, and both agree to
   25 digits with the nuclear norm of Sqrt[r] . Sqrt[s] and the singular values of r - s at 40 digits. *)
VerificationTest[
    BlockRandom[
        With[{
            r = With[{a = RandomInteger[{-4, 4}, {4, 4}] + I RandomInteger[{-4, 4}, {4, 4}]}, # / Tr[#] & [a . ConjugateTranspose[a]]],
            s = With[{a = RandomInteger[{-4, 4}, {4, 4}] + I RandomInteger[{-4, 4}, {4, 4}]}, # / Tr[#] & [a . ConjugateTranspose[a]]]
        },
            With[{f = QuantumDistance[QuantumState[r], QuantumState[s], "Fidelity"], t = QuantumDistance[QuantumState[r], QuantumState[s], "Trace"]},
                {
                    LeafCount[f] < 1000,
                    LeafCount[t] < 1000,
                    Abs[N[1 - f, 30] - Total[SingularValueList[MatrixPower[N[r, 40], 1/2] . MatrixPower[N[s, 40], 1/2], Tolerance -> 0]]] < 10^-25,
                    Abs[N[t, 30] - Total[SingularValueList[N[r - s, 40], Tolerance -> 0]] / 2] < 10^-25
                }
            ]
        ],
        RandomSeeding -> 9
    ],
    {True, True, True, True},
    {},
    TestID -> "Exact-Spectral-CompactAndCorrect"
]

(* Machine pure states in dimension 64: the fidelity distance is 1 - |<a|b>| and the trace distance
   Sqrt[1 - |<a|b>|^2] to round-off, for a random pair and for a pair with |<a|b>| = 10^-6. For that pair
   r . s has one eigenvalue |<a|b>|^2 = 10^-12 and 63 that are round-off of 0, and the round-off in its
   eigenvalues grows, relative to |<a|b>|^2, as the overlap shrinks, so the square roots of those 63 would
   put an error far above the round-off into the fidelity; r and s each have one eigenvalue 1 and round-off
   of 0, and their eigenvectors give the overlap to round-off. *)
VerificationTest[
    BlockRandom[
        With[{a = Normalize[RandomComplex[{-1 - I, 1 + I}, 64]], c = RandomComplex[{-1 - I, 1 + I}, 64]},
            With[{
                b = Normalize[RandomComplex[{-1 - I, 1 + I}, 64]],
                bNearlyOrthogonal = 1.*^-6 a + Sqrt[1 - 1.*^-12] Normalize[c - (Conjugate[a] . c) a]
            },
                Flatten @ Map[
                    Function[v, With[{overlap = Abs[Conjugate[a] . v]},
                        {
                            Abs[QuantumDistance[QuantumState[a], QuantumState[v], "Fidelity"] - (1 - overlap)] < 1.*^-13,
                            Abs[QuantumDistance[QuantumState[a], QuantumState[v], "Trace"] - Sqrt[1 - overlap^2]] < 1.*^-13
                        }
                    ]],
                    {b, bNearlyOrthogonal}
                ]
            ]
        ],
        RandomSeeding -> 10
    ],
    {True, True, True, True},
    {},
    TestID -> "Machine-PureStates-NoSqrtOfRoundoff"
]

(* Pure states at 30 digits in dimension 8, a random pair and a pair with |<a|b>| = 10^-10: the fidelity
   distance keeps its digits to 10^-25 for both. For the second pair the eigenvalue |<a|b>|^2 = 10^-20 of r . s is
   not far above the round-off in the eigenvalues of that product, and a fidelity taken from them keeps far
   fewer of its digits. *)
VerificationTest[
    BlockRandom[
        With[{
            a = Normalize[RandomComplex[{-1 - I, 1 + I}, 8, WorkingPrecision -> 30]],
            b = Normalize[RandomComplex[{-1 - I, 1 + I}, 8, WorkingPrecision -> 30]],
            c = RandomComplex[{-1 - I, 1 + I}, 8, WorkingPrecision -> 30]
        },
            With[{bNearlyOrthogonal = N[10^-10, 30] a + N[Sqrt[1 - 10^-20], 30] Normalize[c - (Conjugate[a] . c) a]},
                Map[
                    Abs[QuantumDistance[QuantumState[a], QuantumState[#]] - (1 - Abs[Conjugate[a] . #])] < 10^-25 &,
                    {b, bNearlyOrthogonal}
                ]
            ]
        ],
        RandomSeeding -> 11
    ],
    {True, True},
    {},
    TestID -> "ArbitraryPrecision-PureStates-KeepDigits"
]

(* A rank-2 state with a repeated eigenvalue, 1/2 twice, against a rank-4 state, at 60 digits. With the
   eigenvectors the columns of u and v, the fidelity is the sum of the singular values of
   M_ij = Sqrt[1/2] <u_i|v_j> Sqrt[q_j], which needs no eigensolver, and the distance agrees with it to
   10^-50; two eigenvectors of the repeated eigenvalue that are orthogonal to far fewer digits than the
   working precision put an error far above 10^-50 into the distance. *)
VerificationTest[
    BlockRandom[
        With[{
            u = Orthogonalize[RandomComplex[{-1 - I, 1 + I}, {6, 6}, WorkingPrecision -> 90]],
            v = Orthogonalize[RandomComplex[{-1 - I, 1 + I}, {6, 6}, WorkingPrecision -> 90]],
            p = {1/2, 1/2, 0, 0, 0, 0},
            q = {1/10, 2/10, 3/10, 4/10, 0, 0}
        },
            Abs[
                QuantumDistance[
                    QuantumState[N[u . DiagonalMatrix[p] . ConjugateTranspose[u], 60]],
                    QuantumState[N[v . DiagonalMatrix[q] . ConjugateTranspose[v], 60]]
                ] -
                (1 - Total[SingularValueList[KroneckerProduct[Sqrt[p[[1 ;; 2]]], Sqrt[q[[1 ;; 4]]]] (ConjugateTranspose[u][[1 ;; 2]] . v[[All, 1 ;; 4]])]])
            ] < 10^-50
        ],
        RandomSeeding -> 7
    ],
    True,
    {},
    TestID -> "ArbitraryPrecision-RepeatedEigenvalue-KeepDigits"
]

(* Two commuting machine states with small populations, 10^-8 and 2 10^-8: the fidelity takes
   Sqrt[10^-8 2 10^-8] from the product of the populations, an eigenvalue of r . s at the level of the
   round-off of its largest, which a cut for round-off in r . s drops, while the eigenvalues of r and s
   themselves are far above round-off and keep it. *)
VerificationTest[
    Abs[
        QuantumDistance[QuantumState[N[DiagonalMatrix[{1 - 10^-8, 10^-8}]]], QuantumState[N[DiagonalMatrix[{1 - 2 10^-8, 2 10^-8}]]]] -
        N[1 - Sqrt[(1 - 10^-8) (1 - 2 10^-8)] - Sqrt[10^-8 2 10^-8]]
    ] < 1.*^-15,
    True,
    {},
    TestID -> "Machine-CommutingStates-KeepSmallPopulations"
]

(* Two diagonal machine states with swapped populations, diag(1 - e, e) and diag(e, 1 - e), for e = 10^-14
   and 10^-20: the fidelity 2 Sqrt[e (1 - e)] lies entirely in the small population e. The populations are
   the eigenvalues, stored exactly, and each is kept however small, so the distance is
   1 - 2 Sqrt[e (1 - e)] to round-off; a cut relative to the largest eigenvalue of each state would drop e
   and give 1. *)
VerificationTest[
    Map[
        Abs[
            QuantumDistance[QuantumState[N[DiagonalMatrix[{1 - #, #}]]], QuantumState[N[DiagonalMatrix[{#, 1 - #}]]]] -
            N[1 - 2 Sqrt[# (1 - #)]]
        ] < 1.*^-15 &,
        {10^-14, 10^-20}
    ],
    {True, True},
    {},
    TestID -> "Machine-SwappedPopulations-KeepSmallPopulations"
]

(* Two orthogonal machine states: every measure built on the fidelity, and the trace distance, is a machine
   number, as the input is, and none is an exact number. *)
VerificationTest[
    QuantumDistance[QuantumState[{1., 0.}], QuantumState[{0., 1.}], #] & /@ {"Fidelity", "Bures", "BuresAngle", "Trace"},
    {1., Sqrt[2.], N[Pi / 2], 1.},
    {},
    TestID -> "Machine-OrthogonalStates-MachineResults"
]

(* A diagonal machine state against a pure state that is not diagonal, diag(1 - 10^-20, 10^-20) against
   {10^-12, Sqrt[1 - 10^-24]}: the fidelity Sqrt[(1 - 10^-20) 10^-24 + 10^-20 (1 - 10^-24)] lies almost all
   in the population 10^-20, which the diagonal state gives as stored, far below a cut relative to its
   largest eigenvalue. *)
VerificationTest[
    Abs[
        QuantumDistance[QuantumState[N[DiagonalMatrix[{1 - 10^-20, 10^-20}]]], QuantumState[N[{10^-12, Sqrt[1 - 10^-24]}]]] -
        N[1 - Sqrt[(1 - 10^-20) 10^-24 + 10^-20 (1 - 10^-24)]]
    ] < 1.*^-15,
    True,
    {},
    TestID -> "Machine-DiagonalAgainstPure-KeepSmallPopulation"
]

(* 30-digit states whose zero entries are known only to an accuracy, as the 0``29.8 of Cos[N[Pi, 30] / 2]:
   such a zero has precision 0, and so does a matrix that holds one, while the nonzero entries carry the 30
   digits. |1> built that way is at distance 0 from |1> and at 1 - 1/Sqrt[2] from |+>, and so is the
   maximally mixed state built with such zeros from |+>. *)
VerificationTest[
    With[{
        one = QuantumState[{Cos[N[Pi, 30] / 2], Sin[N[Pi, 30] / 2]}],
        mixed = QuantumState[{{N[1/2, 30], Cos[N[Pi, 30] / 2]}, {Cos[N[Pi, 30] / 2], N[1/2, 30]}}],
        plus = QuantumState[N[{1, 1} / Sqrt[2], 30]]
    },
        {
            Abs[QuantumDistance[one, QuantumState[N[{0, 1}, 30]]]] < 10^-25,
            Abs[QuantumDistance[one, plus] - (1 - 1 / Sqrt[2])] < 10^-25,
            Abs[QuantumDistance[mixed, plus] - (1 - 1 / Sqrt[2])] < 10^-25
        }
    ],
    {True, True, True},
    {},
    TestID -> "ArbitraryPrecision-AccuracyOnlyZeros-KeepPrecision"
]

(* Two pairs of 4 x 4 states at precision 2. In arbitrary-precision arithmetic, the fidelity of the first
   pair computed from the eigensystems of the two states comes back with an accuracy it does not have, and
   the singular values of the difference of the second pair run for minutes without converging. Computed from
   the machine numbers of the input, the fidelity distance of the first pair and the trace distance of the
   second return at once, and the exact distance of the states the input rounds lies within the uncertainty
   each result states. *)
VerificationTest[
    With[{
        dm = Function[{k, seed},
            BlockRandom[With[{a = RandomInteger[{-4, 4}, {4, k}] + I RandomInteger[{-4, 4}, {4, k}]}, # / Tr[#] &[a . ConjugateTranspose[a]]], RandomSeeding -> seed]]
    },
        Map[
            Function[c, With[{r = dm @@ c[[1]], s = dm @@ c[[2]]},
                With[{d = TimeConstrained[QuantumDistance[QuantumState[N[r, 2]], QuantumState[N[s, 2]], c[[3]]], 30, $Aborted]},
                    NumberQ[d] && Abs[N[d] - N[QuantumDistance[QuantumState[r], QuantumState[s], c[[3]]]]] <= 10^-Accuracy[d]
                ]
            ]],
            {{{4, 911}, {2, 961}, "Fidelity"}, {{2, 909}, {3, 959}, "Trace"}}
        ]
    ],
    {True, True},
    {},
    TestID -> "LowPrecision-NoHang-HonestAccuracy"
]

(* Near machine precision the machine computation drops an eigenvalue below its cut, and the stated
   accuracy accounts for it. A qubit pair whose fidelity 10^-7 lies entirely in an eigenvalue 10^-14 of the
   first state, and a qutrit pair whose trace distance 3/10 includes a singular value 5 10^-15 of the
   difference, both rotated out of the computational basis, at precision 15.9: the exact distance lies
   within the uncertainty each result states. At machine precision the trace distance keeps that singular
   value, to 10^-15. *)
VerificationTest[
    With[{
        rot2 = {{3/5, -4/5}, {4/5, 3/5}},
        rot3 = {{2/3, -2/3, 1/3}, {2/3, 1/3, -2/3}, {1/3, 2/3, 2/3}}
    },
        With[{
            qubits = {rot2 . DiagonalMatrix[{1 - 10^-14, 10^-14}] . Transpose[rot2], rot2 . DiagonalMatrix[{0, 1}] . Transpose[rot2]},
            qutrits = {rot3 . DiagonalMatrix[{1/2, 1/2, 0}] . Transpose[rot3], rot3 . DiagonalMatrix[{1/5, 4/5 - 5 10^-15, 5 10^-15}] . Transpose[rot3]}
        },
            {
                With[{d = QuantumDistance @@ Append[QuantumState[N[#, 159/10]] & /@ qubits, "Fidelity"]},
                    Abs[N[d] - N[1 - 10^-7]] <= 10^-Accuracy[d]],
                With[{d = QuantumDistance @@ Append[QuantumState[N[#, 159/10]] & /@ qutrits, "Trace"]},
                    Abs[N[d] - 3/10] <= 10^-Accuracy[d]],
                Abs[QuantumDistance @@ Append[QuantumState[N[#]] & /@ qutrits, "Trace"] - 3/10] < 1.*^-15
            }
        ]
    ],
    {True, True, True},
    {},
    TestID -> "NearMachinePrecision-StatedAccuracyHolds"
]

(* Two 256 x 256 machine states that QuantumState keeps as SparseArray, one banded and one diagonal: the
   trace distance is computed on the dense difference, with no SingularValueList::arh message, and matches
   half the sum of the absolute eigenvalues of the difference. *)
VerificationTest[
    With[{
        banded = With[{b = IdentityMatrix[256] + 0.3 (DiagonalMatrix[ConstantArray[1., 255], 1] + DiagonalMatrix[ConstantArray[1., 255], -1])}, b / Tr[b]],
        thermal = With[{w = Exp[-5. Range[0, 255] / 256]}, DiagonalMatrix[w / Total[w]]]
    },
        Abs[QuantumDistance[QuantumState[banded], QuantumState[thermal], "Trace"] - Total[Abs[Eigenvalues[banded - thermal]]] / 2] < 1.*^-13
    ],
    True,
    {},
    TestID -> "Machine-SparseStates-TraceNoMessage"
]

(* Two orthogonal 30-digit states, |0> and |1>: the zero entries of N[{1, 0}, 30] are exact, so their
   fidelity is made of exact zeros alone, and given the accuracy of the input it keeps every measure an
   inexact number. *)
VerificationTest[
    With[{d = QuantumDistance[QuantumState[N[{1, 0}, 30]], QuantumState[N[{0, 1}, 30]], #] & /@ {"Fidelity", "Bures", "BuresAngle", "Trace"}},
        {AllTrue[d, InexactNumberQ], Max[Abs[d - {1, Sqrt[2], Pi / 2, 1}]] < 10^-25}
    ],
    {True, True},
    {},
    TestID -> "ArbitraryPrecision-OrthogonalStates-InexactResults"
]

(* A symbolic qubit pair, r on the z axis with Bloch length a and s at angle t with length b, whose
   distances are evaluated at parameter values after they are computed. A closed form that divides by a
   difference of two eigenvalues is 0/0 wherever they coincide: the eigenvalues of
   ConjugateTranspose[r - s] . (r - s) coincide at every real point, and those of r . s at points such as
   a = b = 0, two maximally mixed states, and a = b, t = Pi, two states with swapped populations. The
   spectral sums have no such denominator: the trace distance is half the Bloch distance
   Sqrt[a^2 + b^2 - 2 a b Cos[t]] at a = 3/10, b = 6/10, t = 11/10, as a closed form, and at a pure pair,
   and the fidelity distance is 0 and 1 - Sqrt[1 - a^2] at the two points where the two eigenvalues of
   r . s are equal. *)
VerificationTest[
    With[{
        r = (IdentityMatrix[2] + a PauliMatrix[3]) / 2,
        s = (IdentityMatrix[2] + b (Cos[t] PauliMatrix[3] + Sin[t] PauliMatrix[1])) / 2,
        bloch = Sqrt[a^2 + b^2 - 2 a b Cos[t]] / 2
    },
        With[{trace = QuantumDistance[QuantumState[r], QuantumState[s], "Trace"], fidelity = QuantumDistance[QuantumState[r], QuantumState[s]]},
            {
                Chop[N[trace /. {a -> 3/10, b -> 6/10, t -> 11/10}] - N[bloch /. {a -> 3/10, b -> 6/10, t -> 11/10}]],
                FullSimplify[trace - bloch, 0 < a < 1 && 0 < b < 1 && 0 < t < Pi],
                Chop[N[QuantumDistance[QuantumState[{1, 0}], QuantumState[{Cos[d / 2], Sin[d / 2]}], "Trace"] /. d -> 11/10] - N[Sin[11/20]]],
                N[fidelity /. {a -> 0, b -> 0, t -> 1}],
                Chop[N[fidelity /. {a -> 3/10, b -> 3/10, t -> Pi}] - N[1 - Sqrt[1 - (3/10)^2]]]
            }
        ]
    ],
    {0, 0, 0, 0., 0},
    {},
    TestID -> "Symbolic-Distances-EvaluateAtCoincidentEigenvalues"
]

(* |+> against the traceless Pauli Z, whose expectation in |+> is 0: the product of the two matrices is a
   nonzero 2 x 2 nilpotent matrix, which has no square root, while its eigenvalues are 0 and the fidelity is
   0. The input is not physical, and the distance comes with the notphysical message. *)
VerificationTest[
    QuantumDistance[QuantumState["+"], QuantumState[PauliMatrix[3]]],
    1,
    {QuantumDistance::notphysical},
    TestID -> "NonPhysical-NilpotentProduct-HasFidelity"
]

(* 30-digit states with an eigenvalue below the cut 100 10^-p but several times its uncertainty d 10^-p, on
   whose eigenvector the other state has its weight: dropping it changes the fidelity by its square root,
   about 10^-15, and the exact distance lies within the uncertainty the result states. The states: rank 2
   in dimension 2, 3 and 8, pure but for that eigenvalue, against a pure state; rank 3 against a pure state
   and against a mixed state. *)
VerificationTest[
    With[{
        rot2 = {{3/5, -4/5}, {4/5, 3/5}},
        rot3 = {{2/3, -2/3, 1/3}, {2/3, 1/3, -2/3}, {1/3, 2/3, 2/3}},
        rot8 = IdentityMatrix[8] - ConstantArray[1/4, {8, 8}]
    },
        With[{
            rankThree = rot3 . DiagonalMatrix[{1/2, 1/2 - 3 10^-29, 3 10^-29}] . Transpose[rot3],
            rankTwo = Function[{u, x}, u . DiagonalMatrix[PadRight[{1 - x, x}, Length[u]]] . Transpose[u]],
            pure = Function[{u, j}, Outer[Times, u[[All, j]], u[[All, j]]]]
        },
            Map[
                Function[c, With[{d = QuantumDistance[QuantumState[N[c[[1]], 30]], QuantumState[N[c[[2]], 30]]]},
                    Abs[N[d] - N[c[[3]]]] <= 10^-Accuracy[d]
                ]],
                {
                    {rankTwo[rot2, 10^-29], pure[rot2, 2], 1 - Sqrt[10^-29]},
                    {rankTwo[rot3, 2 10^-29], pure[rot3, 2], 1 - Sqrt[2 10^-29]},
                    {rankTwo[rot8, 3 10^-29], pure[rot8, 2], 1 - Sqrt[3 10^-29]},
                    {rankThree, pure[rot3, 3], 1 - Sqrt[3 10^-29]},
                    {rankThree, rot3 . DiagonalMatrix[{1/2, 0, 1/2}] . Transpose[rot3], 1/2 - Sqrt[3 10^-29 / 2]}
                }
            ]
        ]
    ],
    {True, True, True, True, True},
    {},
    TestID -> "ArbitraryPrecision-DroppedEigenvalue-HonestAccuracy"
]

(* A 30-digit input of trace 1 that is pure but for the eigenvalues 10^-9 and -10^-9, with eigenvectors u_2
   and u_3, which are too small for the notphysical message: it enters as the positive part of its Hermitian
   part, in which -10^-9 is 0. Against the pure state u_2 its fidelity is 10^-(9/2), and against the pure
   state (u_2 + u_3) / Sqrt[2], which also weighs the negative eigenvalue, it is (10^-9 / 2)^(1/2). *)
VerificationTest[
    With[{rot3 = {{2/3, -2/3, 1/3}, {2/3, 1/3, -2/3}, {1/3, 2/3, 2/3}}},
        With[{m = QuantumState[N[rot3 . DiagonalMatrix[{1, 10^-9, -10^-9}] . Transpose[rot3], 30]]},
            Map[
                Function[c, Abs[QuantumDistance[m, QuantumState[N[Outer[Times, c[[1]], c[[1]]], 30]]] - (1 - Sqrt[c[[2]]])] < 10^-20],
                {{rot3[[All, 2]], 10^-9}, {(rot3[[All, 2]] + rot3[[All, 3]]) / Sqrt[2], 10^-9 / 2}}
            ]
        ]
    ],
    {True, True},
    {},
    TestID -> "ArbitraryPrecision-NearlyPureNonPhysical-PositivePart"
]

(* 30-digit states in dimension 8 whose entries are dyadic, so that their 30-digit numbers are exact, built
   on the columns e_j of a unitary matrix of dyadic entries. A pure state a = e_1 against a state of rank 2
   holding it, and against a state whose eigenvector of eigenvalue 1/2 has overlap 2^-33 with a while its
   other eigenvectors are orthogonal to it, so that <a|s|a> = 2^-67: a pure state needs no eigensystem, and
   its fidelity with s comes from the overlaps of a with the eigenvectors of s, so the distance keeps its
   digits to 10^-25, and agrees to 10^-25 with the states in either order. Taken at 30 digits from entries
   of order 1/10, <a|s|a> cancels to 2^-67, and its square root keeps far fewer digits. Two mixed states of
   rank 2 with fidelity 1/4: overlaps of their eigenvectors that are 0 come out with an imaginary part far
   below their accuracy and with few digits, and the distance is still 3/4 to 10^-25. *)
VerificationTest[
    With[{w = DiagonalMatrix[{1, I, 1, I, 1, I, 1, I}] . (IdentityMatrix[8] - ConstantArray[1/4, {8, 8}])},
        With[{a = w[[All, 1]], e = Transpose[w], proj = Outer[Times, #, Conjugate[#]] &},
            Map[
                Function[c, With[{r30 = N[c[[1]], 30], s30 = N[c[[2]], 30]},
                    With[{d = QuantumDistance[QuantumState[r30], QuantumState[s30]]},
                        {Abs[d - c[[3]]] < 10^-25, Abs[d - QuantumDistance[QuantumState[s30], QuantumState[r30]]] < 10^-25}
                    ]
                ]],
                {
                    {a, proj[a + e[[2]] + e[[3]]] / 4 + proj[e[[4]]] / 4, 1/2},
                    {a, proj[2^-33 a + e[[2]]] / 2 + proj[e[[3]]] / 4 + (1/4 - 2^-67) proj[e[[4]]], 1 - 2^-33 / Sqrt[2]},
                    {proj[a + e[[2]] + e[[3]]] / 4 + proj[e[[4]]] / 4, proj[e[[2]] + e[[5]]] / 4 + proj[a - e[[3]]] / 4, 3/4}
                }
            ]
        ]
    ],
    {{True, True}, {True, True}, {True, True}},
    {},
    TestID -> "ArbitraryPrecision-DyadicStates-KeepDigits"
]

(* Three pairs of 30-digit mixed states of rank 3 in dimension 8, with entries in multiples of 1/64 or
   1/128, so that their 30-digit numbers are exact. The matrix M of each pair holds entries with a real or
   imaginary part within a few times its uncertainty but not 0, and SingularValueList on it returned no
   digit, with Divide::indet, or did not return within minutes. Each distance comes back within a minute,
   with no message and more than 20 digits, and holds the exact distance within its uncertainty. *)
VerificationTest[
    Map[
        Function[pair, With[{r = pair[[1]] / 64, s = pair[[2]] / 64},
            With[{d = TimeConstrained[QuantumDistance[QuantumState[N[r, 30]], QuantumState[N[s, 30]]], 60, $Aborted]},
                NumberQ[d] && Accuracy[d] > 20 &&
                    Abs[SetPrecision[d, Infinity] - (1 - Re[Total[Sqrt[Eigenvalues[N[r . s, 120]]]]])] <= 10^-Accuracy[d]
            ]
        ]],
        {
            {
                {{9, 2 + I, 5, 2 - 3 I, 3 + 2 I, 2 + I, 1 - 4 I, 4 + 7 I}, {2 - I, 5, 2 - 5 I, 1, I, -3, -2 - I, 7 + 2 I},
                 {5, 2 + 5 I, 9, 2 + I, -1 + 2 I, 2 - 3 I, -3 - 4 I, 4 + 11 I}, {2 + 3 I, 1, 2 - I, 5, 5 I, 1, -2 + 3 I, 3 + 2 I},
                 {3 - 2 I, -I, -1 - 2 I, -5 I, 5, -I, 3 + 2 I, 2 - 3 I}, {2 - I, -3, 2 + 3 I, 1, I, 5, -2 - I, -1 + 2 I},
                 {1 + 4 I, -2 + I, -3 + 4 I, -2 - 3 I, 3 - 2 I, -2 + I, 9, -8 - I}, {4 - 7 I, 7 - 2 I, 4 - 11 I, 3 - 2 I, 2 + 3 I, -1 - 2 I, -8 + I, 17}},
                {{13, 3 I, -3, 11 I, -3, -I, -3, 3 I}, {-3 I, 5, 5 I, -3, 5 I, 1, 5 I, 5}, {-3, -5 I, 5, 3 I, 5, -I, 5, -5 I},
                 {-11 I, -3, -3 I, 21, -3 I, -7, -3 I, -3}, {-3, -5 I, 5, 3 I, 5, -I, 5, -5 I}, {I, 1, I, -7, I, 5, I, 1},
                 {-3, -5 I, 5, 3 I, 5, -I, 5, -5 I}, {-3 I, 5, 5 I, -3, 5 I, 1, 5 I, 5}}
            },
            {
                {{7, I, I, 2 + I, -1 - 6 I, -2 - 3 I, 4 - 3 I, -I}, {-I, 7, 3 + 4 I, -1 + 2 I, 6 - I, 3 - 2 I, -1, 1},
                 {-I, 3 - 4 I, 5, 5, -3 I, 1 - 4 I, 1, -1 - 2 I}, {2 - I, -1 - 2 I, 5, 15, -8 - I, 3 - 4 I, 5, -3 - 6 I},
                 {-1 + 6 I, 6 + I, 3 I, -8 + I, 15, 4 - 3 I, 3 I, 2 + 3 I}, {-2 + 3 I, 3 + 2 I, 1 + 4 I, 3 + 4 I, 4 + 3 I, 7, 1 + 4 I, 1 - 2 I},
                 {4 + 3 I, -1, 1, 5, -3 I, 1 - 4 I, 5, -1 - 2 I}, {I, 1, -1 + 2 I, -3 + 6 I, 2 - 3 I, 1 + 2 I, -1 + 2 I, 3}},
                {{11, -2 - I, -1 + 4 I, 4 - 7 I, -5 - 8 I, -2 - 5 I, -15 - 2 I, -2 - I}, {-2 + I, 7, 2 + 5 I, 5 - 2 I, -2 + 9 I, 3, -9 I, 7},
                 {-1 - 4 I, 2 - 5 I, 11, -8 - 3 I, 7 + 12 I, 2 - I, -11 + 2 I, 2 - 5 I}, {4 + 7 I, 5 + 2 I, -8 + 3 I, 19, -4 - 9 I, 1 - 6 I, -2 - 11 I, 5 + 2 I},
                 {-5 + 8 I, -2 - 9 I, 7 - 12 I, -4 + 9 I, 27, 6 - 5 I, -7 - 2 I, -2 - 9 I}, {-2 + 5 I, 3, 2 + I, 1 + 6 I, 6 + 5 I, 7, -13 I, 3},
                 {-15 + 2 I, 9 I, -11 - 2 I, -2 + 11 I, -7 + 2 I, 13 I, 39, 9 I}, {-2 + I, 7, 2 + 5 I, 5 - 2 I, -2 + 9 I, 3, -9 I, 7}} / 2
            },
            {
                {{5, -4 - I, -3 + 4 I, -I, 1, -I, -3, -I}, {-4 + I, 13, 8 - 7 I, -3 + 4 I, -4 - 3 I, -3 + 4 I, -4 - 7 I, -3 + 4 I},
                 {-3 - 4 I, 8 + 7 I, 13, -4 - I, 1 - 4 I, -4 - I, -3 - 4 I, -4 - I}, {I, -3 - 4 I, -4 + I, 5, 5 I, 5, I, 5},
                 {1, -4 + 3 I, 1 + 4 I, -5 I, 5, -5 I, 1, -5 I}, {I, -3 - 4 I, -4 + I, 5, 5 I, 5, I, 5},
                 {-3, -4 + 7 I, -3 + 4 I, -I, 1, -I, 13, -I}, {I, -3 - 4 I, -4 + I, 5, 5 I, 5, I, 5}},
                {{2, -2, 2, -2 I, -2, -2 I, 2 I, -2 I}, {-2, 6, -2, 2 I, -2 - 4 I, 2 I, 4 - 2 I, 2 I}, {2, -2, 18, -2 I, -18, -2 I, 2 I, -2 I},
                 {2 I, -2 I, 2 I, 2, -2 I, 2, -2, 2}, {-2, -2 + 4 I, -18, 2 I, 26, 2 I, -4 + 2 I, 2 I}, {2 I, -2 I, 2 I, 2, -2 I, 2, -2, 2},
                 {-2 I, 4 + 2 I, -2 I, -2, -4 - 2 I, -2, 6, -2}, {2 I, -2 I, 2 I, 2, -2 I, 2, -2, 2}}
            }
        }
    ],
    {True, True, True},
    {},
    TestID -> "ArbitraryPrecision-DyadicMixedStates-NoCollapse"
]


EndTestSection[]
