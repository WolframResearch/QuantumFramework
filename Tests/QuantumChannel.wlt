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

EndTestSection[]
