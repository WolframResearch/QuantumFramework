BeginTestSection["QuantumPositiveTransform - sign-sector transform"]

(* QuantumPositiveTransform on a framework object must dispatch to the typed rule
   for that object.  An unrestricted array rule was tried ahead of the object
   rules, so it captured the object, ran ArrayReshape on it, and returned an
   unevaluated SparseArray (ArrayReshape::listrp) or a Failure.  Each call below
   must return the matching head with no message. *)

VerificationTest[
    Head @ QuantumPositiveTransform[QuantumState["0", 3]],
    QuantumState,
    TestID -> "Positive-on-Schrodinger-state"
]

VerificationTest[
    Head @ QuantumPositiveTransform[QuantumWignerTransform[QuantumState["0", 3]]],
    QuantumState,
    TestID -> "Positive-on-phase-space-state"
]

VerificationTest[
    Head @ QuantumPositiveTransform[QuditBasis["Wigner"[3]]],
    QuditBasis,
    TestID -> "Positive-on-qudit-basis"
]

VerificationTest[
    Head @ QuantumPositiveTransform[QuantumOperator["X"[3]]],
    QuantumOperator,
    TestID -> "Positive-on-operator"
]

(* The array form is the mathematical content: it splits a quasiprobability vector
   into a nonnegative pair (w+, w-) whose difference is the original vector. *)
VerificationTest[
    With[{w = QuantumPhaseSpaceTransform[QuantumState["0", 3]]["StateVector"]},
        With[{p = Normal @ QuantumPositiveTransform[w]},
            {Dimensions[p], p[[1]] - p[[2]] == Normal[w], AllTrue[Flatten[p], # >= 0 &]}
        ]
    ],
    {{2, 9}, True, True},
    TestID -> "Positive-array-splits-into-nonnegative-pair"
]

EndTestSection[]
