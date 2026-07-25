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

(* each picture name is accepted as a positional argument, not routed into QuditBasis *)

VerificationTest[QuantumBasis["Schrodinger"]["Picture"], "Schrodinger", TestID -> "Picture-bare-Schrodinger"]

VerificationTest[QuantumBasis["Heisenberg"]["Picture"], "Heisenberg", TestID -> "Picture-bare-Heisenberg"]

VerificationTest[QuantumBasis["Interaction"]["Picture"], "Interaction", TestID -> "Picture-bare-Interaction"]

VerificationTest[QuantumBasis["PhaseSpace"]["Picture"], "PhaseSpace", TestID -> "Picture-bare-PhaseSpace"]

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

(* an Integer tail is still a multiplicity, not a dimension *)

VerificationTest[QuantumBasis["Heisenberg", 3]["Dimensions"], {2, 2, 2}, TestID -> "Picture-integer-tail-is-multiplicity"]

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
