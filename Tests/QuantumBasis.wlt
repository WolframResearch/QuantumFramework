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
