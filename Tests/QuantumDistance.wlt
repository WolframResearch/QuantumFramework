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

(* Relative entropy to pure state returns entropy of first state *)
VerificationTest[
    QuantumDistance[qsMixed, qs0, "RelativeEntropy"],
    qsMixed["Entropy"],
    TestID -> "RelativeEntropy-ToPureState"
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
