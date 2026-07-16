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


BeginTestSection["QuantumMeasurement - Symbolic Probabilities"]

(* Outcome ordering ("TopProbabilities") and Monte-Carlo sampling ("SimulatedMeasurement",
   "SimulatedCounts", "SimulatedStateMeasurement"), together with the generic
   "DistributionInformation" passthrough, are defined only when every outcome probability
   is numeric. On a symbolic / parametric measurement the old code fed a non-numeric
   distribution into Information / RandomVariate and leaked CategoricalDistribution::elmntavsl,
   MultinomialDistribution::vprobprm, Extract::psl1, RandomVariate::unsdst, KeyMap::invak, etc.,
   including on plain display. These now return Indeterminate (mirroring "Entropy"). The empty
   expected-message list makes any leaked message fail the test. *)

(* --- the reported regression: no leaked messages, Indeterminate value --- *)

(* TopProbabilities needs to order the weights (was CategoricalDistribution::elmntavsl) *)
VerificationTest[
    Module[{t}, QuantumCircuitOperator[{QuantumState[{1, 0}], "RY"[t], {1}}][]["TopProbabilities"]],
    Indeterminate,
    {},
    TestID -> "SymProb-TopProbabilities"
]

(* SimulatedCounts samples a MultinomialDistribution (was MultinomialDistribution::vprobprm) *)
VerificationTest[
    Module[{t}, QuantumCircuitOperator[{QuantumState[{1, 0}], "RY"[t], {1}}][]["SimulatedCounts", 10]],
    Indeterminate,
    {},
    TestID -> "SymProb-SimulatedCounts"
]

(* SimulatedMeasurement samples the CategoricalDistribution; previously returned an
   unevaluated RandomVariate[...] instead of a clean value *)
VerificationTest[
    QuantumMeasurement[<|0 -> p, 1 -> 1 - p|>]["SimulatedMeasurement"],
    Indeterminate,
    {},
    TestID -> "SymProb-SimulatedMeasurement"
]

VerificationTest[
    QuantumMeasurement[<|0 -> p, 1 -> 1 - p|>]["SimulatedMeasurement", 3],
    Indeterminate,
    {},
    TestID -> "SymProb-SimulatedMeasurement-n"
]

(* the explicit DistributionInformation path with a sampling sub-property: the exact pair
   from the Entropy regression, Extract::psl1 + RandomVariate::unsdst *)
VerificationTest[
    QuantumMeasurement[<|0 -> p, 1 -> 1 - p|>]["DistributionInformation", "Entropy"],
    Indeterminate,
    {},
    TestID -> "SymProb-DistributionInformation-Entropy"
]

(* same path under N, the route the summary box used to take on display *)
VerificationTest[
    N @ QuantumMeasurement[<|0 -> p, 1 -> 1 - p|>]["DistributionInformation", "Entropy"],
    Indeterminate,
    {},
    TestID -> "SymProb-DistributionInformation-Entropy-N"
]

(* SimulatedStateMeasurement keys the state association by a simulated outcome *)
VerificationTest[
    QuantumMeasurement[<|0 -> p, 1 -> 1 - p|>]["SimulatedStateMeasurement"],
    Indeterminate,
    {},
    TestID -> "SymProb-SimulatedStateMeasurement"
]

(* list form was RandomVariate::array + Part::pkspec1 *)
VerificationTest[
    QuantumMeasurement[<|0 -> p, 1 -> 1 - p|>]["SimulatedStateMeasurement", 3],
    Indeterminate,
    {},
    TestID -> "SymProb-SimulatedStateMeasurement-n"
]

(* TopStateProbabilities depends on the (now guarded) TopProbabilities; was KeyMap::invak *)
VerificationTest[
    QuantumMeasurement[<|0 -> p, 1 -> 1 - p|>]["TopStateProbabilities"],
    Indeterminate,
    {},
    TestID -> "SymProb-TopStateProbabilities"
]

(* display path stays quiet on a symbolic measurement *)
VerificationTest[
    Module[{t}, ToBoxes @ QuantumCircuitOperator[{QuantumState[{1, 0}], "RY"[t], {1}}][]; True],
    True,
    {},
    TestID -> "SymProb-Display"
]

(* --- properties that are well defined symbolically must NOT degrade to Indeterminate --- *)

(* Categories is just the outcome support, independent of the probability values *)
VerificationTest[
    QuantumMeasurement[<|0 -> p, 1 -> 1 - p|>]["Categories"],
    QuantumMeasurement[<|0 -> p, 1 -> 1 - p|>]["Outcomes"],
    {},
    TestID -> "SymProb-Categories-Preserved"
]

(* ProbabilityArray is the symbolic weight vector, fully resolved (no leftover Information) *)
VerificationTest[
    With[{a = QuantumMeasurement[<|0 -> p, 1 -> 1 - p|>]["ProbabilityArray"]},
        Head[a] === List && FreeQ[a, _Information | _CategoricalDistribution]
    ],
    True,
    {},
    TestID -> "SymProb-ProbabilityArray-Preserved"
]

VerificationTest[
    Head @ QuantumMeasurement[<|0 -> p, 1 -> 1 - p|>]["ProbabilityTable"],
    Dataset,
    {},
    TestID -> "SymProb-ProbabilityTable-Preserved"
]

(* --- numeric measurements are unchanged and stay message-free --- *)

VerificationTest[
    Length @ QuantumMeasurementOperator["Z"][QuantumState[{1, Sqrt[3]} / 2]]["TopProbabilities"],
    2,
    {},
    TestID -> "SymProb-Numeric-TopProbabilities"
]

VerificationTest[
    Total @ QuantumMeasurementOperator["Z"][QuantumState[{1, Sqrt[3]} / 2]]["SimulatedCounts", 50],
    50,
    {},
    TestID -> "SymProb-Numeric-SimulatedCounts"
]

VerificationTest[
    With[{qm = QuantumMeasurementOperator["Z"][QuantumState[{1, Sqrt[3]} / 2]]},
        MemberQ[qm["Outcomes"], qm["SimulatedMeasurement"]]
    ],
    True,
    {},
    TestID -> "SymProb-Numeric-SimulatedMeasurement"
]

VerificationTest[
    NumericQ @ QuantumMeasurementOperator["Z"][QuantumState[{1, Sqrt[3]} / 2]]["DistributionInformation", "Entropy"],
    True,
    {},
    TestID -> "SymProb-Numeric-DistributionInformation"
]

EndTestSection[]


BeginTestSection["QuantumMeasurementOperator - degenerate observable eigenbasis"]

(* A rank-1 projector on a qutrit has a degenerate 0-eigenspace; at arbitrary precision
   Eigensystem returns non-orthonormal eigenvectors there, and projectors built from the
   skewed basis silently corrupted the outcome probabilities for some KCBS pentagon
   projectors (k = 2, 3) while others were fine. All five expectation values are exactly
   Cos[alpha]^2 = 1/Sqrt[5] by symmetry. *)
VerificationTest[
    Block[{sol, alpha, v, P, qs},
        sol = Solve[Cos[a]^2 + Sin[a]^2 Cos[4 Pi/5] == 0 && 0 < a < Pi/2, a];
        alpha = a /. First[sol];
        v[k_] := N[{Cos[alpha], Sin[alpha] Cos[4 Pi k/5], Sin[alpha] Sin[4 Pi k/5]}, 16];
        P[k_] := Outer[Times, v[k], v[k]];
        qs = QuantumState[{1, 0, 0}, {3}];
        Max @ Abs[N @ Table[QuantumMeasurementOperator[P[k]][qs]["Mean"], {k, 0, 4}] - N[1/Sqrt[5]]] < 1*^-10
    ],
    True,
    TestID -> "DegenerateEigenbasis-KCBS-arbitrary-precision"
]

(* same pentagon at machine precision *)
VerificationTest[
    Block[{sol, alpha, v, P, qs},
        sol = Solve[Cos[a]^2 + Sin[a]^2 Cos[4 Pi/5] == 0 && 0 < a < Pi/2, a];
        alpha = a /. First[sol];
        v[k_] := N @ {Cos[alpha], Sin[alpha] Cos[4 Pi k/5], Sin[alpha] Sin[4 Pi k/5]};
        P[k_] := Outer[Times, v[k], v[k]];
        qs = QuantumState[{1, 0, 0}, {3}];
        Max @ Abs[Table[QuantumMeasurementOperator[P[k]][qs]["Mean"], {k, 0, 4}] - N[1/Sqrt[5]]] < 1*^-8
    ],
    True,
    TestID -> "DegenerateEigenbasis-KCBS-machine-precision"
]

EndTestSection[]


BeginTestSection["QuantumMeasurementOperator - order re-seating on direct state application"]

(* A tensor-product observable whose factors were built with off-1 wire labels (output
   wire {3} each, renumbered to {3, 4} by the tensor product) failed with
   QuantumCircuitOperator::dim when applied directly to a two-qutrit state: phantom
   wires 1-2 were padded at qubit dimension 2. A directly applied measurement has no
   circuit context, so wire labels that do not fit the state are re-seated onto the
   state's qudits. *)
VerificationTest[
    Block[{P0, opAB, qmo, w, qs},
        P0 = Outer[Times, {1., 0., 0.}, {1., 0., 0.}];
        opAB = QuantumTensorProduct[QuantumOperator[P0, {3}, {1}], QuantumOperator[P0, {3}, {2}]];
        qmo = QuantumMeasurementOperator[opAB];
        w = Normalize[{1., 2., 0.}];
        qs = QuantumState[Flatten[Outer[Times, w, w]], {3, 3}];
        {Head[qmo[qs]], Abs[N @ qmo[qs]["Mean"] - 1/25] < 1*^-8}
    ],
    {QuantumMeasurement, True},
    TestID -> "OrderReseat-tensor-product-observable"
]

(* a measurement whose order fits inside a larger state is untouched by the re-seat *)
VerificationTest[
    Round[N @ Values @ QuantumMeasurementOperator["Z", {2}][
        QuantumCircuitOperator[{"X" -> 2}][QuantumState["000"]]]["Probabilities"], 0.001],
    {0., 1.},
    TestID -> "OrderReseat-preserves-targeted-measurement"
]

EndTestSection[]
