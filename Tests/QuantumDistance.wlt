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


EndTestSection[]
