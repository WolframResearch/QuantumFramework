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


BeginTestSection["QuantumMeasurement - Entropy"]

(* Shannon entropy of the outcome distribution, H = -Sum p_i Log_b p_i, in closed form
   (mirrors QuantumState["VonNeumannEntropy", logBase]). Computing it via
   Information[CategoricalDistribution[..]] is unsafe: on symbolic / non-numeric
   probabilities (any parametric circuit) that path falls back to Monte-Carlo sampling
   and emits RandomVariate::unsdst + Extract::psl1, including on plain display. *)

(* --- known numeric values --- *)

VerificationTest[
    QuantumCircuitOperator[{"H", "CNOT", {1, 2}}][]["Entropy"],
    Quantity[1, "Bits"],
    TestID -> "Entropy-Bell-1bit"
]

VerificationTest[
    QuantumCircuitOperator[{"H", {1}}][]["Entropy"],
    Quantity[1, "Bits"],
    TestID -> "Entropy-PlusState-1bit"
]

(* deterministic outcome carries no information *)
VerificationTest[
    QuantumMeasurementOperator["Z"][QuantumState["0"]]["Entropy"],
    Quantity[0, "Bits"],
    TestID -> "Entropy-Deterministic-0bit"
]

(* biased {1/4, 3/4}: H = 2 - (3/4) Log2[3], exact *)
VerificationTest[
    Simplify[QuantumMeasurementOperator["Z"][QuantumState[{1, Sqrt[3]} / 2]]["Entropy"] - Quantity[2 - 3/4 Log2[3], "Bits"]],
    Quantity[0, "Bits"],
    TestID -> "Entropy-Biased-ExactShannon"
]

(* --- logBase argument: parity with QuantumState["Entropy", logBase] (QM lacked it) --- *)

(* base e returns the entropy in nats *)
VerificationTest[
    QuantumCircuitOperator[{"H", "CNOT", {1, 2}}][]["Entropy", E],
    Log[2],
    TestID -> "Entropy-LogBase-Nats"
]

(* {"Entropy", logBase} list-form dispatches to the bare base-b value *)
VerificationTest[
    QuantumCircuitOperator[{"H", "CNOT", {1, 2}}][][{"Entropy", 2}],
    1,
    TestID -> "Entropy-ListForm-Base2"
]

(* base conversion is consistent on a non-trivial distribution: nats = bits * Ln[2] *)
VerificationTest[
    With[{qm = QuantumMeasurementOperator["Z"][QuantumState[{1, Sqrt[3]} / 2]]},
        FullSimplify[qm["Entropy", E] - qm["Entropy", 2] Log[2]]
    ],
    0,
    TestID -> "Entropy-LogBase-Consistency"
]

(* --- symbolic robustness (the reported regression) --- *)

(* a bare parametric rotation measured in Z keeps symbolic probabilities (theta is not
   assumed real); qm["Entropy"] must stay a Quantity and emit NO messages. The empty
   expected-message list makes any leaked message fail the test. *)
VerificationTest[
    Module[{t}, Head @ QuantumCircuitOperator[{QuantumState[{1, 0}], "RY"[t], {1}}][]["Entropy"]],
    Quantity,
    {},
    TestID -> "Entropy-Symbolic-NoMessages"
]

(* the summary box evaluates N @ qm["Entropy"] on display, so the display path must be
   quiet too *)
VerificationTest[
    Module[{t}, ToBoxes @ QuantumCircuitOperator[{QuantumState[{1, 0}], "RY"[t], {1}}][]; True],
    True,
    {},
    TestID -> "Entropy-Symbolic-Display-NoMessages"
]

(* symbolic entropy is correct, not merely quiet: at theta = Pi/2 the Z outcomes are
   equiprobable, so H = 1 bit *)
VerificationTest[
    Module[{t}, Simplify[QuantumCircuitOperator[{QuantumState[{1, 0}], "RY"[t], {1}}][]["Entropy"] /. t -> Pi / 2]],
    Quantity[1, "Bits"],
    TestID -> "Entropy-Symbolic-CorrectAtPiOver2"
]

EndTestSection[]
