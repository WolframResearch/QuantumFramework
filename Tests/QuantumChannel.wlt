BeginTestSection["QuantumChannel - constructors"]

VerificationTest[QuantumChannel["BitFlip"]["Dimensions"], {2, 2, 2}, TestID -> "BitFlip-bare"]

VerificationTest[QuantumChannel["BitFlip"[0.1]]["Dimensions"], {2, 2, 2}, TestID -> "BitFlip-call"]

VerificationTest[QuantumChannel["PhaseFlip"[0.2]]["Dimensions"], {2, 2, 2}, TestID -> "PhaseFlip"]

VerificationTest[QuantumChannel["BitPhaseFlip"[0.2]]["Dimensions"], {2, 2, 2}, TestID -> "BitPhaseFlip"]

VerificationTest[QuantumChannel["Depolarizing"[0.2]]["Dimensions"], {4, 2, 2}, TestID -> "Depolarizing"]

VerificationTest[QuantumChannel["AmplitudeDamping"[0.5]]["Dimensions"], {2, 2, 2}, TestID -> "AmplitudeDamping"]

VerificationTest[
    QuantumChannel["GeneralizedAmplitudeDamping"[0.3, 0.4]]["Dimensions"],
    {4, 2, 2},
    TestID -> "GeneralizedAmplitudeDamping"
]

VerificationTest[QuantumChannel["PhaseDamping"[0.3]]["Dimensions"], {2, 2, 2}, TestID -> "PhaseDamping"]

VerificationTest[QuantumChannel["ResetError"[0.1, 0.1]]["Dimensions"], {5, 2, 2}, TestID -> "ResetError"]

EndTestSection[]


BeginTestSection["QuantumChannel - action on states"]

(* BitFlip[1] always flips: |0> -> |1>, fidelity 0 with |0> *)
VerificationTest[
    Total @ Flatten @ Chop[QuantumChannel["BitFlip"[1]][QuantumState["0"]]["MatrixRepresentation"] - {{0, 0}, {0, 1}}],
    0,
    TestID -> "BitFlip-1-on-0"
]

(* BitFlip[0] is identity *)
VerificationTest[
    Total @ Flatten[QuantumChannel["BitFlip"[0]][QuantumState["0"]]["MatrixRepresentation"] - {{1, 0}, {0, 0}}],
    0,
    TestID -> "BitFlip-0-identity"
]

(* Depolarizing[1] => maximally mixed *)
VerificationTest[
    Total @ Flatten @ Chop[QuantumChannel["Depolarizing"[1]][QuantumState["0"]]["MatrixRepresentation"] - IdentityMatrix[2]/2],
    0,
    TestID -> "Depolarizing-1-mixed"
]

EndTestSection[]


BeginTestSection["QuantumChannel - failure"]

VerificationTest[
    Quiet @ QuantumChannel["NotAChannel"[0.1]],
    Failure["InvalidName", _],
    SameTest -> MatchQ,
    TestID -> "InvalidName-call-form"
]

VerificationTest[
    Quiet @ QuantumChannel["NotAChannel"],
    Failure["InvalidName", _],
    SameTest -> MatchQ,
    TestID -> "InvalidName-bare"
]

VerificationTest[
    Quiet @ QuantumChannel["BitFlip"["bad-arg-1", "bad-arg-2", "bad-arg-3"]],
    Failure["InvalidArguments", _],
    SameTest -> MatchQ,
    TestID -> "InvalidArgs-BitFlip"
]

EndTestSection[]


BeginTestSection["QuantumChannel - trace preservation and unitality"]

(* "Trace" is the trace-preservation witness Sum_i K_i^dagger K_i and must be the
   identity for every named channel. It previously composed the operator the other
   way round and returned Sum_i K_i K_i^dagger, so the two channels with a
   non-Hermitian jump operator were reported as not trace preserving. *)
VerificationTest[
    QuantumChannel[#]["TracePreservingQ"] & /@ {
        "BitFlip"[1/5], "PhaseFlip"[1/5], "BitPhaseFlip"[1/5], "Depolarizing"[1/5],
        "PhaseDamping"[1/5], "AmplitudeDamping"[1/5],
        "GeneralizedAmplitudeDamping"[1/5, 1/3], "ResetError"[1/5, 1/5]
    },
    ConstantArray[True, 8],
    TestID -> "Channel-TracePreserving-AllNamed"
]

(* Amplitude damping is the smallest case that separates the two witnesses:
   Sum K^dag K = I while Sum K K^dag = diag(1 + gamma, 1 - gamma). *)
VerificationTest[
    With[{qc = QuantumChannel["AmplitudeDamping"[1/5]]},
        {Normal @ qc["Adjoint"]["Unitality"]["MatrixRepresentation"], Normal @ qc["Unitality"]["MatrixRepresentation"]}
    ],
    {{{1, 0}, {0, 1}}, {{6/5, 0}, {0, 4/5}}},
    TestID -> "Channel-AdjointUnitality-vs-Unitality-AmplitudeDamping"
]

(* Unitality is a genuinely different question, so it must still say No where the
   channel is not unital. Amplitude damping drives every state toward |0>, so it
   cannot fix the identity; generalized amplitude damping is unital only at p = 1/2. *)
VerificationTest[
    QuantumChannel[#]["UnitalQ"] & /@ {
        "BitFlip"[1/5], "PhaseFlip"[1/5], "BitPhaseFlip"[1/5], "Depolarizing"[1/5],
        "PhaseDamping"[1/5], "AmplitudeDamping"[1/5],
        "GeneralizedAmplitudeDamping"[1/5, 1/2], "GeneralizedAmplitudeDamping"[1/5, 1/4]
    },
    {True, True, True, True, True, False, True, False},
    TestID -> "Channel-UnitalQ-AllNamed"
]

(* "Unitality" is N(I), so it must agree with applying the channel to the
   maximally mixed state and rescaling, which never touches either witness. *)
VerificationTest[
    With[{names = {"PhaseDamping"[1/5], "AmplitudeDamping"[1/5], "GeneralizedAmplitudeDamping"[1/5, 1/4]}},
        And @@ (
            With[{qc = QuantumChannel[#]},
                Chop[Normal[qc["Unitality"]["MatrixRepresentation"]] -
                     2 Normal[qc[QuantumState[IdentityMatrix[2] / 2]]["DensityMatrix"]]] == {{0, 0}, {0, 0}}
            ] & /@ names
        )
    ],
    True,
    TestID -> "Channel-Unitality-equals-channel-of-identity"
]

(* A Kraus set scaled off unit norm is not trace preserving, and a symbolic
   parameter leaves both witnesses undecidable rather than silently True. *)
VerificationTest[
    {
        QuantumChannel[{2 {{1, 0}, {0, Sqrt[4/5]}}, 2 {{0, Sqrt[1/5]}, {0, 0}}}]["TracePreservingQ"],
        QuantumChannel["AmplitudeDamping"[g]]["TracePreservingQ"],
        QuantumChannel["AmplitudeDamping"[g]]["UnitalQ"]
    },
    {False, Undefined, Undefined},
    TestID -> "Channel-NonTracePreserving-and-symbolic"
]

(* More than one Kraus qudit, and a channel wider than one wire. *)
VerificationTest[
    {
        QuantumChannel["AmplitudeDamping"[1/5], {1, 2}]["TracePreservingQ"],
        QuantumTensorProduct[QuantumChannel["AmplitudeDamping"[1/5]], QuantumChannel["PhaseDamping"[3/10]]]["TraceQudits"],
        QuantumTensorProduct[QuantumChannel["AmplitudeDamping"[1/5]], QuantumChannel["PhaseDamping"[3/10]]]["TracePreservingQ"],
        QuantumChannel[{
            {{1, 0, 0}, {0, Sqrt[3/4], 0}, {0, 0, Sqrt[1/2]}},
            {{0, Sqrt[1/4], 0}, {0, 0, 0}, {0, 0, 0}},
            {{0, 0, 0}, {0, 0, Sqrt[1/2]}, {0, 0, 0}}}]["TracePreservingQ"]
    },
    {True, 2, True, True},
    TestID -> "Channel-TracePreserving-MultiQuditShapes"
]

(* Both witnesses are advertised, so they are discoverable rather than only
   reachable by knowing the name. *)
VerificationTest[
    MemberQ[QuantumChannel["Properties"], #] & /@ {"Unitality", "UnitalQ", "TracePreservingQ"},
    {True, True, True},
    TestID -> "Channel-witness-properties-are-advertised"
]

EndTestSection[]


BeginTestSection["QuantumChannel - action on operators and circuits"]

(* Applying a channel to a framework operator forms the two-element circuit
   {op, channel} and compiles it.  It used to forward an undeclared "Trace" option
   into the circuit compiler, which raised OptionValue::nodef.  A default
   VerificationTest fails on any emitted message, so these pin down that the calls
   are now message-free and return an operator. *)

VerificationTest[
    Head @ QuantumChannel["BitFlip"[1/10]][QuantumOperator["Z"]],
    QuantumOperator,
    TestID -> "Channel-on-operator-message-free"
]

VerificationTest[
    Head @ QuantumChannel["BitFlip"[1/10]][QuantumCircuitOperator[{"H"}]],
    QuantumOperator,
    TestID -> "Channel-on-circuit-message-free"
]

(* The compiled value is unchanged: it equals compiling the explicit nested
   circuit {circuit, channel}. *)
VerificationTest[
    Normal @ QuantumChannel["BitFlip"[1/10]][QuantumCircuitOperator[{"H"}]]["Matrix"] ===
        Normal @ QuantumCircuitOperator[{QuantumCircuitOperator[{"H"}], QuantumChannel["BitFlip"[1/10]]}]["QuantumOperator"]["Matrix"],
    True,
    TestID -> "Channel-on-circuit-value-unchanged"
]

EndTestSection[]


BeginTestSection["QuantumChannel - single Kraus operator"]

(* A one-element Kraus list is the deterministic channel rho -> M rho M^dagger.
   The general list path stacked the operators against an environment qudit of
   dimension Length[ops], which for one operator is 1: it collapsed and left the
   lone output wire on the non-positive environment label, so quantumChannelQ
   rejected the object.  "ValidQ" guards application and every property downvalue,
   so the channel could neither be applied to a state nor answer "TracePreservingQ".
   The one-operator case now routes to the single-operator constructor. *)

VerificationTest[
    {
        QuantumChannel[{IdentityMatrix[2]}]["ValidQ"],
        QuantumChannel[{PauliMatrix[1]}]["ValidQ"],
        QuantumChannel[{QuantumOperator["X"]}]["ValidQ"]
    },
    {True, True, True},
    TestID -> "Channel-single-Kraus-ValidQ"
]

(* The identity Kraus is the identity channel: |0><0| stays |0><0|. *)
VerificationTest[
    Normal @ QuantumChannel[{IdentityMatrix[2]}][QuantumState["0"]]["DensityMatrix"],
    {{1, 0}, {0, 0}},
    TestID -> "Channel-single-Kraus-identity-applies"
]

(* A single X Kraus, matrix or named operator, flips |0> to |1>. *)
VerificationTest[
    {
        Normal @ QuantumChannel[{PauliMatrix[1]}][QuantumState["0"]]["DensityMatrix"],
        Normal @ QuantumChannel[{QuantumOperator["X"]}][QuantumState["0"]]["DensityMatrix"]
    },
    {{{0, 0}, {0, 1}}, {{0, 0}, {0, 1}}},
    TestID -> "Channel-single-Kraus-X-applies"
]

(* M^dagger M = I makes the one-Kraus channel trace preserving. *)
VerificationTest[
    QuantumChannel[{IdentityMatrix[2]}]["TracePreservingQ"],
    True,
    TestID -> "Channel-single-Kraus-TracePreserving"
]

(* Validity is structural, not a CPTP check: a non-isometric single Kraus still
   builds a valid, applicable object, but it is not trace preserving. *)
VerificationTest[
    {
        QuantumChannel[{2 IdentityMatrix[2]}]["ValidQ"],
        QuantumChannel[{2 IdentityMatrix[2]}]["TracePreservingQ"]
    },
    {True, False},
    TestID -> "Channel-single-Kraus-nonTracePreserving"
]

(* The single-element list matches the single-operator constructor bit for bit. *)
VerificationTest[
    Normal @ QuantumChannel[{QuantumOperator["X"]}]["QuantumOperator"]["Matrix"] ===
        Normal @ QuantumChannel[QuantumOperator["X"]]["QuantumOperator"]["Matrix"],
    True,
    TestID -> "Channel-single-Kraus-matches-single-operator"
]

EndTestSection[]
