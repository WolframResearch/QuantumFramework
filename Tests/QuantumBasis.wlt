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
   Swept over the registry rather than sampled. Note the sweep quantifies over the registry as
   it stands at load: the rules below match against Alternatives @@ $QuantumBasisPictures,
   which is expanded into the pattern when the package loads, so a name appended to the
   registry at runtime is not picked up by them and recurses *)

VerificationTest[
    AssociationMap[QuantumBasis[#]["Picture"] &, $QuantumBasisPictures],
    AssociationThread[
        $QuantumBasisPictures ->
            $QuantumBasisPictures
    ],
    TestID -> "Picture-bare-every-registered-picture"
]

VerificationTest[
    AssociationMap[QuantumBasis["Computational", #]["Picture"] &, $QuantumBasisPictures],
    AssociationThread[
        $QuantumBasisPictures ->
            $QuantumBasisPictures
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

(* "PhaseSpace" is set for the names in $QuditPhaseSpaceBasisNames, at one place in the
   constructor. It is a tag, not something read back off the elements, and nothing ties the
   two together, so elements and tag can be put in disagreement by several routes that these
   tests do not close. What is closed here is the positional route on an already-tagged
   basis: renaming it sends an already-transformed basis back through QuantumWignerTransform
   and the quasi-probability comes out doubled. The array below is the value at d = 2, kept
   as a tripwire on that doubling; it is not a normalized quasi-probability of any state,
   and being read off the reshaped state vector it does not discriminate the elements. *)

VerificationTest[
    Normal @ QuantumState[{a, b, c, d}, QuantumBasis["Wigner"]]["PhaseSpace"],
    {{a, b, a, b}, {c, d, -c, -d}, {a, -b, a, -b}, {c, -d, -c, d}},
    TestID -> "PhaseSpace-Wigner-closed-form-is-the-reference"
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
    Union[QuantumBasis[#, "PhaseSpace"]["Picture"] & /@ $QuditPhaseSpaceBasisNames],
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


BeginTestSection["QuantumBasis - parameter substitution"]

(* Substituting a parameter used to rebuild the basis from its representations
   unconditionally.  A basis declares its parameters in ParameterSpec whether or
   not its elements depend on them, so the common case rebuilt a basis identical
   to the one it started from.  The rebuild is now skipped when the substitution
   turns out to have changed nothing; these tests pin both branches. *)

(* Elements that genuinely depend on the parameter must still substitute. *)

$rotationBasis := QuantumBasis[
    QuditBasis[<|
        QuditName["a"] -> {Cos[\[FormalX]], Sin[\[FormalX]]},
        QuditName["b"] -> {-Sin[\[FormalX]], Cos[\[FormalX]]}
    |>],
    "ParameterSpec" -> {{\[FormalX], 0, 1}}
]

VerificationTest[
    Normal /@ Values @ $rotationBasis[<|\[FormalX] -> 0|>]["Output"]["Representations"],
    {{1, 0}, {0, 1}},
    TestID -> "ParameterSubstitution-rotates-elements-at-zero"
]

VerificationTest[
    Normal /@ Values @ $rotationBasis[<|\[FormalX] -> Pi / 2|>]["Output"]["Representations"],
    {{0, 1}, {-1, 0}},
    TestID -> "ParameterSubstitution-rotates-elements-at-quarter-turn"
]

(* Substituting spends the parameter either way. *)
VerificationTest[
    $rotationBasis[<|\[FormalX] -> 0|>]["ParameterArity"],
    0,
    TestID -> "ParameterSubstitution-spends-the-parameter"
]

(* A basis that declares a parameter its elements do not use takes the fast
   path, and must give exactly what the rebuild gave: the same elements, the
   same picture and label, one fewer parameter. *)

$declaredOnly := QuantumBasis[QuditBasis[2], "ParameterSpec" -> {{\[FormalX], 0, 1}}]

VerificationTest[
    Block[{qb = $declaredOnly, sub},
        sub = qb[<|\[FormalX] -> 0.5|>];
        {
            sub["Output"]["Representations"] === qb["Output"]["Representations"],
            sub["Input"]["Representations"] === qb["Input"]["Representations"],
            sub["Picture"] === qb["Picture"],
            sub["Label"] === qb["Label"],
            sub["ParameterArity"]
        }
    ],
    {True, True, True, True, 0},
    TestID -> "ParameterSubstitution-unused-parameter-keeps-the-basis"
]

(* An evolved state carries its parameter on the basis while the amplitudes
   alone depend on it, so binding it must leave the basis untouched. *)
VerificationTest[
    Block[{psi, bound},
        psi = QuantumEvolve[
            QuantumOperator["PauliX"], {} -> {},
            QuantumState["0"],
            {\[FormalT], 0, 1}
        ];
        bound = psi[0.5];
        {
            bound["Basis"]["Output"]["Representations"] === psi["Basis"]["Output"]["Representations"],
            bound["Basis"]["ParameterArity"],
            psi["Basis"]["ParameterArity"]
        }
    ],
    {True, 0, 1},
    TestID -> "ParameterSubstitution-evolved-state-keeps-its-basis"
]

EndTestSection[]
