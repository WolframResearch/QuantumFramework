BeginTestSection["QuantumMeasurementOperator - constructors"]

VerificationTest[QuantumMeasurementOperator["RandomHermitian"[2]]["Dimensions"], {2, 2}, TestID -> "RandomHermitian-2"]

VerificationTest[QuantumMeasurementOperator["GellMannMICPOVM"[2]]["Dimensions"], {4, 2, 2}, TestID -> "GellMannMICPOVM-2"]

VerificationTest[QuantumMeasurementOperator["TetrahedronSICPOVM"]["Dimensions"], {4, 2, 2}, TestID -> "TetrahedronSICPOVM-bare"]

VerificationTest[QuantumMeasurementOperator["QBismSICPOVM"[2]]["Dimensions"], {4, 2, 2}, TestID -> "QBismSICPOVM-2"]

VerificationTest[QuantumMeasurementOperator["HesseSICPOVM"]["Dimensions"], {9, 3, 3}, TestID -> "HesseSICPOVM-bare"]

VerificationTest[QuantumMeasurementOperator["HoggarSICPOVM"]["Dimensions"], {64, 8, 8}, TestID -> "HoggarSICPOVM-bare"]

VerificationTest[QuantumMeasurementOperator[1]["Targets"], {{1}}, TestID -> "Integer-target"]

EndTestSection[]


BeginTestSection["QuantumMeasurementOperator - basic action"]

(* projector measurement: probabilities of |+> in computational basis *)
VerificationTest[
    Chop @ Values @ QuantumMeasurementOperator[QuantumOperator["I"]][QuantumState["+"]]["Probabilities"],
    {1/2, 1/2},
    TestID -> "Plus-projector-probabilities"
]

EndTestSection[]


BeginTestSection["QuantumMeasurementOperator - failure"]

VerificationTest[
    Quiet @ QuantumMeasurementOperator["NotAMeasurement"[2]],
    Failure["InvalidName", _],
    SameTest -> MatchQ,
    TestID -> "InvalidName-call-form"
]

EndTestSection[]
