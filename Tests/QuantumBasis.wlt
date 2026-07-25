BeginTestSection["QuantumBasis - constructors"]

VerificationTest[QuantumBasis["Bell"]["Dimensions"], {4}, TestID -> "Bell-bare"]

VerificationTest[QuantumBasis["Fourier"[2]]["Dimensions"], {2}, TestID -> "Fourier-2"]

VerificationTest[QuantumBasis[QuditBasis[2], QuditBasis[3]]["Dimensions"], {2, 3}, TestID -> "Output-Input-explicit"]

VerificationTest[QuantumBasis[QuantumBasis["Computational"], 3]["Dimension"], 8, TestID -> "Multiplicity-3"]

VerificationTest[
    QuantumBasis[{QuditBasis[2], QuditBasis[3]}]["Dimensions"],
    {2, 3},
    TestID -> "Tensor-list"
]

EndTestSection[]


BeginTestSection["QuantumBasis - properties"]

VerificationTest[QuantumBasis["Bell"]["Qudits"], 1, TestID -> "Bell-Qudits"]

VerificationTest[QuantumBasis["Fourier"[3]]["OutputDimension"], 3, TestID -> "Fourier-OutputDim"]

VerificationTest[QuantumBasis[QuditBasis[2]]["Label"], None, TestID -> "QuditBasis-Label"]

EndTestSection[]


BeginTestSection["QuantumBasis - equality"]

VerificationTest[
    QuantumBasis[QuditBasis[2], QuditBasis[2]] == QuantumBasis[QuditBasis[2], QuditBasis[2]],
    True,
    TestID -> "Equal-same"
]

VerificationTest[
    QuantumBasis[QuditBasis[2], QuditBasis[3]] == QuantumBasis[QuditBasis[3], QuditBasis[2]],
    False,
    TestID -> "Unequal-different"
]

EndTestSection[]


BeginTestSection["QuantumBasis - pictures"]

(* every registered picture is accepted as a positional argument, not routed into QuditBasis.
   Swept over the registry rather than sampled, so a name added to $QuantumBasisPictures is
   covered here without editing this file *)

VerificationTest[
    AssociationMap[QuantumBasis[#]["Picture"] &, Wolfram`QuantumFramework`PackageScope`$QuantumBasisPictures],
    AssociationThread[
        Wolfram`QuantumFramework`PackageScope`$QuantumBasisPictures ->
            Wolfram`QuantumFramework`PackageScope`$QuantumBasisPictures
    ],
    TestID -> "Picture-bare-every-registered-picture"
]

VerificationTest[
    AssociationMap[QuantumBasis["Computational", #]["Picture"] &, Wolfram`QuantumFramework`PackageScope`$QuantumBasisPictures],
    AssociationThread[
        Wolfram`QuantumFramework`PackageScope`$QuantumBasisPictures ->
            Wolfram`QuantumFramework`PackageScope`$QuantumBasisPictures
    ],
    TestID -> "Picture-after-Computational-every-registered-picture"
]

(* a bare picture says nothing about the basis, so the basis stays at the default *)

VerificationTest[QuantumBasis["Heisenberg"]["Dimensions"], {2}, TestID -> "Picture-bare-default-dimensions"]

VerificationTest[QuantumBasis["Schrodinger"] === QuantumBasis[], True, TestID -> "Picture-Schrodinger-is-default"]

(* the option form, and its agreement with each positional form *)

VerificationTest[
    QuantumBasis["Computational", "Picture" -> "Heisenberg"]["Picture"],
    "Heisenberg",
    TestID -> "Picture-option-form"
]

VerificationTest[
    QuantumBasis["Heisenberg"] === QuantumBasis["Picture" -> "Heisenberg"],
    True,
    TestID -> "Picture-bare-matches-option"
]

VerificationTest[
    QuantumBasis["Computational", "Heisenberg"] === QuantumBasis["Computational", "Picture" -> "Heisenberg"],
    True,
    TestID -> "Picture-positional-matches-option"
]

(* a picture may sit on either side of the basis specification *)

VerificationTest[QuantumBasis["Computational", "Heisenberg"]["Picture"], "Heisenberg", TestID -> "Picture-after-Computational"]

VerificationTest[QuantumBasis["PauliX", "Interaction"]["Picture"], "Interaction", TestID -> "Picture-after-PauliX"]

VerificationTest[QuantumBasis[3, "Heisenberg"]["Dimensions"], {3}, TestID -> "Picture-after-dimension"]

VerificationTest[
    QuantumBasis["Heisenberg", "PauliX"] === QuantumBasis["PauliX", "Heisenberg"],
    True,
    TestID -> "Picture-position-symmetric"
]

VerificationTest[QuantumBasis["Heisenberg", "Label" -> "x"]["Label"], "x", TestID -> "Picture-with-option-tail"]

(* a picture also survives at the tail of the forms that reach the option rule as
   QuantumBasis["Output" -> ..., "Input" -> ..., picture] rather than as a bare name *)

VerificationTest[QuantumBasis["X", "Y", "Heisenberg"]["Picture"], "Heisenberg", TestID -> "Picture-trailing-after-two-names"]

VerificationTest[QuantumBasis["X", "Y", "Heisenberg"]["Dimensions"], {2, 2}, TestID -> "Picture-trailing-keeps-both-qudits"]

VerificationTest[
    QuantumBasis[QuditBasis[2], QuditBasis[3], "Heisenberg"]["Dimensions"],
    {2, 3},
    TestID -> "Picture-trailing-after-QuditBasis-pair"
]

VerificationTest[QuantumBasis[{"X"}, {2}, "Heisenberg"]["Picture"], "Heisenberg", TestID -> "Picture-trailing-after-list-forms"]

VerificationTest[QuantumBasis[QuditBasis[2], "Heisenberg"]["Picture"], "Heisenberg", TestID -> "Picture-trailing-after-QuditBasis"]

(* whichever picture is written FIRST wins, positional or option alike: the arguments are
   folded over Reverse @ {args}, so the earliest one is applied last *)

VerificationTest[
    QuantumBasis["Picture" -> "Interaction", "Heisenberg"]["Picture"],
    "Interaction",
    TestID -> "Picture-first-wins-option-then-positional"
]

VerificationTest[
    QuantumBasis["Heisenberg", "Picture" -> "Interaction"]["Picture"],
    "Heisenberg",
    TestID -> "Picture-first-wins-positional-then-option"
]

VerificationTest[
    QuantumBasis["Heisenberg", "Interaction"]["Picture"],
    "Heisenberg",
    TestID -> "Picture-first-wins-two-positionals"
]

VerificationTest[
    QuantumBasis["Picture" -> "PhaseSpace", "Picture" -> "Heisenberg"]["Picture"],
    "PhaseSpace",
    TestID -> "Picture-first-wins-two-options"
]

(* an Integer tail is still a multiplicity, not a dimension *)

VerificationTest[QuantumBasis["Heisenberg", 3]["Dimensions"], {2, 2, 2}, TestID -> "Picture-integer-tail-is-multiplicity"]

(* a negative tail is not a multiplicity, and fails exactly as it does without a picture *)

VerificationTest[
    FailureQ @ QuantumBasis["Heisenberg", -1],
    True,
    {QuditBasis::invalidArgs},
    TestID -> "Picture-negative-integer-tail-fails"
]

(* QuantumBasis[{}] answers with a QuditBasis rather than a QuantumBasis, which is its own
   pre-existing defect. QuantumBasis[{}, picture] used to be a Failure instead, because the
   picture reached QuditBasis as a name; excluding it here sends the pair to the list rule,
   so the pictured form now agrees with the unpictured one. That agreement is what is pinned,
   not that either head is right. Restricting the list rule to a non-empty list would fix the
   head and does break other suites, so the head is left where it was *)

VerificationTest[
    Head @ QuantumBasis[{}, "Heisenberg"],
    Head @ QuantumBasis[{}],
    TestID -> "Picture-empty-list-matches-unpictured"
]

(* A phase-space basis derives "PhaseSpace" from its own elements, so renaming its picture
   would leave the elements and the tag disagreeing. QuantumState reads the phase-space route
   off that tag, so a demoted Wigner basis sends an already-transformed basis back through
   QuantumWignerTransform and the quasi-probability comes out doubled. The closed form for a
   symbolic qubit amplitude vector is the reference: it must not depend on a picture name. *)

VerificationTest[
    Normal @ QuantumState[{a, b, c, d}, QuantumBasis["Wigner"]]["PhaseSpace"],
    {{a, b, a, b}, {c, d, -c, -d}, {a, -b, a, -b}, {c, -d, -c, d}},
    TestID -> "PhaseSpace-Wigner-closed-form-is-the-reference"
]

(* the invariant the closed form exists to protect: naming a picture must not change the
   quasi-probability. This is what catches a demotion that slips past the guard by any
   route, rather than only the routes enumerated below *)

VerificationTest[
    Normal @ QuantumState[{a, b, c, d}, QuantumBasis["Wigner"]]["PhaseSpace"] ===
        Normal @ QuantumState[{a, b, c, d}, QuantumBasis[QuantumBasis["Wigner"], "PhaseSpace"]]["PhaseSpace"],
    True,
    TestID -> "PhaseSpace-independent-of-picture-name"
]

VerificationTest[
    FailureQ @ QuantumBasis[QuantumBasis["Heisenberg"], "Wigner"],
    True,
    {QuantumBasis::phaseSpacePicture},
    TestID -> "PhaseSpace-basis-refuses-demotion-through-a-merged-basis"
]

(* the qubit grid is 2d x 2d because d is even; an odd d is the d x d case, and it is the
   one that shows the 4x4 above is a doubled grid rather than a d^2 reshape *)

VerificationTest[
    Dimensions @ QuantumState[Array[x, 9], QuantumBasis["Wigner"[3]]]["PhaseSpace"],
    {3, 3},
    TestID -> "PhaseSpace-odd-dimension-grid-is-d-by-d"
]

VerificationTest[
    FailureQ @ QuantumBasis["Wigner", "Heisenberg"],
    True,
    {QuantumBasis::phaseSpacePicture},
    TestID -> "PhaseSpace-basis-refuses-demotion-trailing"
]

VerificationTest[
    FailureQ @ QuantumBasis["Heisenberg", "Wigner"],
    True,
    {QuantumBasis::phaseSpacePicture},
    TestID -> "PhaseSpace-basis-refuses-demotion-leading"
]

(* swept over the whole phase-space registry, not sampled: PhaseSpace on a phase-space basis
   stays a no-op. The refusing direction is pinned one message at a time above, because a
   sweep of it would emit one message per name and trip General::stop *)

VerificationTest[
    Union[QuantumBasis[#, "PhaseSpace"]["Picture"] & /@ Wolfram`QuantumFramework`PackageScope`$QuditPhaseSpaceBasisNames],
    {"PhaseSpace"},
    TestID -> "PhaseSpace-registry-accepts-its-own-picture"
]

(* a Failure handed in as the specification reports once, not twice *)

VerificationTest[
    FailureQ @ QuantumBasis["Heisenberg", QuditBasis["NoSuchThing"]],
    True,
    {QuditBasis::invalidName},
    TestID -> "Picture-failed-spec-reports-once"
]

(* a tail that specifies no basis at all fails by name rather than silently *)

VerificationTest[
    FailureQ @ QuantumBasis["Heisenberg", Sqrt[2]],
    True,
    {QuantumBasis::invalidSpec},
    TestID -> "Picture-uninterpretable-tail-names-itself"
]

(* the two-name rule the picture guard sits on must keep reading its second argument as an input basis *)

VerificationTest[QuantumBasis["Computational", "Pauli"]["Dimensions"], {2, 4}, TestID -> "Picture-guard-keeps-Computational-Pauli"]

VerificationTest[QuantumBasis["X", "Y"]["Dimensions"], {2, 2}, TestID -> "Picture-guard-keeps-X-Y"]

VerificationTest[QuantumBasis[{"X"}, {2}]["Dimensions"], {2, 2}, TestID -> "Picture-guard-keeps-list-forms"]

VerificationTest[QuantumBasis[3, 2]["Dimensions"], {3, 3}, TestID -> "Picture-guard-keeps-integer-pair"]

(* an unrecognized basis name still fails, and the message names the basis, never the picture;
   one message per test, since three copies of one message symbol would trip General::stop *)

VerificationTest[
    FailureQ @ QuantumBasis["NoSuchBasis"],
    True,
    {QuditBasis::invalidName},
    TestID -> "Picture-unknown-name-still-fails"
]

VerificationTest[
    FailureQ @ QuantumBasis["Heisenberg", "NoSuchBasis"],
    True,
    {QuditBasis::invalidName},
    TestID -> "Picture-unknown-name-after-picture-fails"
]

EndTestSection[]
