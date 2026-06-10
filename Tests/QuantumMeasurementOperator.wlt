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


BeginTestSection["QuantumMeasurementOperator - Computational property"]

(* gh#3: ["Computational"] used to silently corrupt non-Z measurement statistics
   by rebuilding the QMO with the operator re-seated into the computational frame,
   stripping the eigenvalue-tagged pointer labels. It should now leave a non-
   computational measurement untouched. *)
VerificationTest[
    Values @ QuantumMeasurementOperator["X", {1}]["Computational"][QuantumState["+"]]["Probabilities"],
    {1, 0},
    TestID -> "Computational-X-plus"
]

VerificationTest[
    Values @ QuantumMeasurementOperator["X", {1}]["Computational"][QuantumState["-"]]["Probabilities"],
    {0, 1},
    TestID -> "Computational-X-minus"
]

VerificationTest[
    Values @ QuantumMeasurementOperator["Y", {1}]["Computational"][QuantumState["+"]]["Probabilities"],
    {1/2, 1/2},
    TestID -> "Computational-Y-plus"
]

(* a Z measurement is already in the computational frame and continues to be
   rebuilt through the standard path *)
VerificationTest[
    Values @ QuantumMeasurementOperator["Z", {1}]["Computational"][QuantumState["0"]]["Probabilities"],
    {1, 0},
    TestID -> "Computational-Z-zero"
]

(* circuit-level: QuantumCircuitOperator[...]["Computational"] maps the property
   over each element; the non-Z measurement inside must stay faithful *)
VerificationTest[
    Values @ QuantumCircuitOperator[{"H" -> 1, QuantumMeasurementOperator["X", {1}]}]["Computational"][QuantumState["0"]]["Probabilities"],
    {1, 0},
    TestID -> "Computational-circuit-H-then-X"
]

EndTestSection[]


BeginTestSection["QuantumMeasurementOperator - failure"]

VerificationTest[
    Quiet @ QuantumMeasurementOperator["NotAMeasurement"[2]],
    Failure["InvalidName", _],
    SameTest -> MatchQ,
    TestID -> "InvalidName-call-form"
]

VerificationTest[
    Quiet @ QuantumMeasurementOperator["HesseSICPOVM"["bad-arg"]],
    Failure["InvalidArguments", _],
    SameTest -> MatchQ,
    TestID -> "InvalidArgs-HesseSICPOVM"
]

EndTestSection[]
