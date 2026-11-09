(* Applying one operator or channel to a small state asks the operator, its basis and the
   circuit it is wrapped in for their shapes and orders many thousands of times. These tests
   pin the answers of the routes that serve those questions: a basis's dimensions read straight
   off its representations, the circuit's option names collected once, a channel or measurement
   handing an operator property to its operator, and a network compiled with its basis built
   once. Each expected value is derived by hand, and the applications run with the property
   cache both off and on.

   Regimes covered:
     exact          X|0> = |1>; H then CNOT on |00> is (|00> + |11>)/Sqrt[2]
     symbolic       BitFlip[p] on |0> is diag(1 - p, p); H on (a, b) is (a + b, a - b)/Sqrt[2]
     numerical      an 8-qubit circuit of H, CNOT and R_Y layers against its unitary built by
                    hand from Kronecker products
     edge           the empty basis, a basis with a dimension-1 factor, a zero-dimensional
                    basis, operators with only input or only output wires, an unknown property *)

$cacheSettings = {False, True}

withCache[c_, expr_] := Block[{$QuantumFrameworkPropCache = c}, expr]
SetAttributes[withCache, HoldRest]


(* basis shapes *)

VerificationTest[
    keyGroupBy[<|{a, 1} -> 1, {b, 2} -> 2, {c, 1} -> 3|>, Last, Total],
    <|1 -> 4, 2 -> 2|>,
    TestID -> "keyGroupBy-reduces-each-group"
]

VerificationTest[
    keyGroupBy[<|{a, 1} -> 1, {b, 2} -> 2, {c, 1} -> 3|>, Last],
    <|1 -> <|{a, 1} -> 1, {c, 1} -> 3|>, 2 -> <|{b, 2} -> 2|>|>,
    TestID -> "keyGroupBy-default-keeps-groups"
]

VerificationTest[
    keyGroupBy[<||>, Last, Length],
    <||>,
    TestID -> "keyGroupBy-empty"
]

VerificationTest[
    Map[QuditBasis[{2, 3}], {"Shape", "Dimensions", "Dimension", "Qudits", "Size", "NameRank"}],
    {{2, 3}, {2, 3}, 6, 2, 6, 2},
    TestID -> "QuditBasis-shape-two-qudits"
]

VerificationTest[
    Map[QuditBasis[{2, 1, 3}], {"Shape", "Dimensions", "Dimension", "Qudits"}],
    {{2, 1, 3}, {2, 3}, 6, 2},
    TestID -> "QuditBasis-shape-drops-unit-factor"
]

VerificationTest[
    Map[QuditBasis[], {"Shape", "Dimensions", "Dimension", "Qudits"}],
    {{1}, {}, 1, 0},
    TestID -> "QuditBasis-shape-empty"
]

VerificationTest[
    Map[QuditBasis[0], {"Shape", "Dimensions", "Dimension", "Qudits"}],
    {{0}, {0}, 0, 1},
    TestID -> "QuditBasis-shape-zero-dimensional"
]

VerificationTest[
    Map[QuditBasis["Pauli"], {"Shape", "Dimensions", "Dimension", "Qudits"}],
    {{4}, {4}, 4, 1},
    TestID -> "QuditBasis-shape-Pauli"
]

VerificationTest[
    Map[QuantumBasis[{2, 3}, {3}], {
        "Dimension", "InputDimension", "OutputDimension", "Dimensions", "InputDimensions", "OutputDimensions",
        "Qudits", "InputQudits", "OutputQudits", "NameDimension", "NameDimensions"
    }],
    {18, 3, 6, {2, 3, 3}, {3}, {2, 3}, 3, 1, 2, 18, {2, 3, 3}},
    TestID -> "QuantumBasis-shape-read-off-output-and-input"
]


(* operator orders *)

VerificationTest[
    Map[QuantumOperator["CNOT", {3, 1}], {"FullOrder", "FullInputOrder", "FullOutputOrder"}],
    {{{3, 1}, {3, 1}}, {3, 1}, {3, 1}},
    TestID -> "QuantumOperator-full-order-unsorted"
]

VerificationTest[
    QuantumOperator["Cup", {-1, 0}]["FullOrder"],
    {{-1, 0}, {}},
    TestID -> "QuantumOperator-full-order-output-only"
]

VerificationTest[
    QuantumOperator[QuantumState["0"]["Dual"], {{}, {2}}]["FullOrder"],
    {{}, {2}},
    TestID -> "QuantumOperator-full-order-input-only"
]


(* applications, cache off and on *)

Scan[
    Function[cache,
        VerificationTest[
            withCache[cache, Normal @ QuantumOperator["X"][QuantumState["0"]]["StateVector"]],
            {0, 1},
            TestID -> "X-on-0-cache-" <> ToString[cache]
        ];

        VerificationTest[
            withCache[cache, Normal @ QuantumOperator["X"][QuantumState["0"], Method -> "Schrodinger"]["StateVector"]],
            {0, 1},
            TestID -> "X-on-0-Schrodinger-cache-" <> ToString[cache]
        ];

        VerificationTest[
            withCache[cache, Simplify[Normal @ QuantumOperator["H"][QuantumState[{a, b}]]["StateVector"] - {a + b, a - b} / Sqrt[2]]],
            {0, 0},
            TestID -> "H-on-symbolic-cache-" <> ToString[cache]
        ];

        VerificationTest[
            withCache[cache, Normal @ QuantumCircuitOperator[{"H", "CNOT"}][]["StateVector"]],
            {1, 0, 0, 1} / Sqrt[2],
            TestID -> "Bell-circuit-cache-" <> ToString[cache]
        ];

        VerificationTest[
            withCache[cache, Simplify[Normal @ QuantumChannel["BitFlip"[p]][QuantumState["0"]]["DensityMatrix"], 0 <= p <= 1]],
            {{1 - p, 0}, {0, p}},
            TestID -> "BitFlip-symbolic-cache-" <> ToString[cache]
        ];

        VerificationTest[
            withCache[cache, Normal @ QuantumChannel["BitFlip"[1 / 10]][QuantumState["0"]]["DensityMatrix"]],
            {{9 / 10, 0}, {0, 1 / 10}},
            TestID -> "BitFlip-exact-cache-" <> ToString[cache]
        ];

        VerificationTest[
            withCache[cache, Chop[Normal @ QuantumChannel["BitFlip"[0.1]][QuantumState["0"]]["DensityMatrix"] - {{0.9, 0}, {0, 0.1}}]],
            {{0, 0}, {0, 0}},
            TestID -> "BitFlip-machine-cache-" <> ToString[cache]
        ];

        VerificationTest[
            withCache[cache, Normal @ QuantumOperator["RZ"[\[Theta]]]["Dagger"]["Matrix"] ===
                ConjugateTranspose[Normal @ QuantumOperator["RZ"[\[Theta]]]["Matrix"]]],
            True,
            TestID -> "RZ-dagger-cache-" <> ToString[cache]
        ]
    ],
    $cacheSettings
]


(* an 8-qubit circuit against its unitary *)

hMatrix = {{1, 1}, {1, -1}} / Sqrt[2];

ryMatrix[t_] := {{Cos[t / 2], -Sin[t / 2]}, {Sin[t / 2], Cos[t / 2]}};

cnotMatrix = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}};

onWire[g_, q_, n_] := KroneckerProduct[IdentityMatrix[2 ^ (q - 1)], g, IdentityMatrix[2 ^ (n - q - Log2[Length[g]] + 1)]]

eightQubitSpec = Flatten[{Table["H" -> i, {i, 8}], Table["CNOT" -> {i, i + 1}, {i, 7}], Table["RY"[0.3 i] -> i, {i, 8}]}];

eightQubitUnitary = Dot @@ Reverse @ Join[
    Table[onWire[hMatrix, q, 8], {q, 8}],
    Table[onWire[cnotMatrix, q, 8], {q, 7}],
    Table[onWire[ryMatrix[0.3 q], q, 8], {q, 8}]
];

Scan[
    Function[cache,
        VerificationTest[
            withCache[cache, Chop[
                Normal @ QuantumCircuitOperator[eightQubitSpec][QuantumState["Register"[8]]]["StateVector"] -
                    eightQubitUnitary[[All, 1]],
                10 ^ -12
            ]],
            ConstantArray[0, 2 ^ 8],
            TestID -> "8-qubit-circuit-cache-" <> ToString[cache]
        ]
    ],
    $cacheSettings
]


(* properties a channel or measurement hands to its operator *)

VerificationTest[
    Map[QuantumChannel["BitFlip"[p]], {"Order", "OutputOrder", "InputOrder", "TraceOrder", "TraceQudits"}],
    {{{0, 1}, {1}}, {0, 1}, {1}, {0}, 1},
    TestID -> "QuantumChannel-orders"
]

VerificationTest[
    Map[QuantumMeasurementOperator["Z"], {"Order", "OutputOrder", "Eigenorder"}],
    {{{0, 1}, {1}}, {0, 1}, {0}},
    TestID -> "QuantumMeasurementOperator-orders"
]

VerificationTest[
    QuantumMeasurementOperator["Z"][QuantumState["+"]]["ProbabilitiesList"],
    {1 / 2, 1 / 2},
    TestID -> "QuantumMeasurement-delegated-probabilities"
]

VerificationTest[
    MatchQ[QuantumChannel["BitFlip"[p]]["Zzz"], _QuantumCircuitOperator["Zzz"]],
    True,
    {QuantumChannel::undefprop, QuantumCircuitOperator::undefprop},
    TestID -> "QuantumChannel-unknown-property-still-undefined"
]


(* circuit options and dimensions *)

VerificationTest[
    QuantumCircuitOperator[{"H", "CNOT"}, "WireLabels" -> None]["WireLabels"],
    None,
    TestID -> "circuit-option-stored"
]

VerificationTest[
    QuantumCircuitOperator[{"H", "CNOT"}]["ShowGateLabels"],
    None,
    TestID -> "circuit-option-absent-is-None"
]

VerificationTest[
    QuantumCircuitOperator[{"H", "CNOT"}, "Label" -> "Bell"]["Label"],
    "Bell",
    TestID -> "circuit-label-option"
]

VerificationTest[
    QuantumCircuitOperator[{"RX"[t]}, "Parameters" -> {t}]["Parameters"],
    {t},
    TestID -> "circuit-parameters-option"
]

VerificationTest[
    QuantumCircuitOperator[{"H"}, Parameters -> {x}, Picture -> "Heisenberg"]["Elements"][[1]]["Basis"]["Options"],
    {"Picture" -> "Heisenberg", "Label" -> "H", "ParameterSpec" -> {{x, 0, 1}}},
    TestID -> "circuit-own-options-as-symbols"
]

VerificationTest[
    Map[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 3}, "BitFlip"[p] -> 2}], {"InputDimensions", "OutputDimensions", "Order"}],
    {{2, 2, 2}, {2, 2, 2, 2}, {{0, 1, 2, 3}, {1, 2, 3}}},
    TestID -> "circuit-dimensions-with-channel"
]

VerificationTest[
    QuantumCircuitOperator[{"H"}, "Foo" -> 3]["Elements"][[1]]["Label"],
    "H",
    {OptionValue::nodef, OptionValue::nodef},
    TestID -> "circuit-unknown-option-warns"
]

VerificationTest[
    Internal`InheritedBlock[{QuantumCircuitOperator},
        SetOptions[QuantumCircuitOperator, "Parameters" -> {s}, "Picture" -> "Heisenberg"];
        QuantumCircuitOperator[{"H"}]["Elements"][[1]]["Basis"]["Options"]
    ],
    {"Picture" -> "Heisenberg", "Label" -> "H", "ParameterSpec" -> {{s, 0, 1}}},
    TestID -> "circuit-SetOptions-defaults-apply"
]

VerificationTest[
    QuantumCircuitOperator[{"H"}]["Elements"][[1]]["Basis"]["Options"],
    {"Picture" -> "Schrodinger", "Label" -> "H", "ParameterSpec" -> {}},
    TestID -> "circuit-defaults-after-SetOptions-block"
]

VerificationTest[
    Internal`InheritedBlock[{CircuitDraw, QuantumCircuitOperator},
        QuantumCircuitOperator[{"H"}];
        SetOptions[CircuitDraw, "WireLabels" -> None];
        Lookup[Options[QuantumCircuitOperator], "WireLabels"]
    ],
    None,
    TestID -> "circuit-options-follow-CircuitDraw-after-first-circuit"
]
