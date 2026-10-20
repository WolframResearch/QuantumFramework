(* Package loading is checked by RunTests.wls before any VerificationTest.
   These checks exercise the local shorthand cache installed during loading. *)
$InitializationShorthand = Wolfram`QuantumFramework`PackageScope`FromOperatorShorthand;

VerificationTest[
    With[{name = SymbolName[Unique["InitializationHit$"]], reference = QuantumOperator["H"]},
        Block[{initializationCalls = 0},
            Internal`InheritedBlock[{QuantumOperator},
                Unprotect[QuantumOperator];
                HoldPattern[QuantumOperator[name]] := (initializationCalls++; reference);
                {$InitializationShorthand[name] === reference,
                    $InitializationShorthand[name] === reference, initializationCalls}
            ]
        ]
    ],
    {True, True, 1},
    TestID -> "Initialization-shorthand-reuses-successful-results"
]

VerificationTest[
    With[{name = SymbolName[Unique["InitializationFailure$"]]},
        Block[{initializationCalls = 0},
            Internal`InheritedBlock[{QuantumOperator},
                Unprotect[QuantumOperator];
                QuantumOperator::initprobe = "temporary initialization test failure";
                HoldPattern[QuantumOperator[name]] := (initializationCalls++; Message[QuantumOperator::initprobe]; $Failed);
                {$InitializationShorthand[name], $InitializationShorthand[name], initializationCalls}
            ]
        ]
    ],
    {$Failed, $Failed, 2},
    {QuantumOperator::initprobe, QuantumOperator::initprobe},
    TestID -> "Initialization-shorthand-does-not-cache-failures-or-their-messages"
]

VerificationTest[
    With[{head = Unique["InitializationPattern$"]},
        Block[{initializationCalls = 0},
            Internal`InheritedBlock[{QuantumOperator},
                Unprotect[QuantumOperator];
                DownValues[QuantumOperator] = Prepend[DownValues[QuantumOperator],
                    HoldPattern[QuantumOperator[head[arg_]]] :> (initializationCalls++; HoldComplete[arg])];
                With[{fixtureValid = QuantumOperator[head[0]] === HoldComplete[0]},
                    initializationCalls = 0;
                    {fixtureValid, $InitializationShorthand[head[_]], $InitializationShorthand[head[1]],
                        $InitializationShorthand[head[2]], $InitializationShorthand[head[_]], initializationCalls}
                ]
            ]
        ]
    ],
    {True, HoldComplete[_], HoldComplete[1], HoldComplete[2], HoldComplete[_], 3},
    TestID -> "Initialization-pattern-keys-do-not-become-matching-rules"
]

VerificationTest[
    With[{name = SymbolName[Unique["InitializationRecursive$"]], reference = QuantumOperator["X"]},
        Block[{initializationCalls = 0},
            Internal`InheritedBlock[{QuantumOperator},
                Unprotect[QuantumOperator];
                HoldPattern[QuantumOperator[name]] := (initializationCalls++; reference);
                {$InitializationShorthand[{name, name}] === {reference, reference},
                    $InitializationShorthand[{name, name}] === {reference, reference}, initializationCalls}
            ]
        ]
    ],
    {True, True, 1},
    TestID -> "Initialization-recursive-shorthand-retains-caching"
]

VerificationTest[
    Head[$InitializationShorthand[]] === $InitializationShorthand &&
        Head[$InitializationShorthand[1, 2]] === $InitializationShorthand,
    True,
    TestID -> "Initialization-unanswered-shorthand-remains-unevaluated"
]

VerificationTest[
    With[{theta = Unique["initializationAngle$"]},
        FullSimplify[Normal[$InitializationShorthand["RX"[theta]]["MatrixRepresentation"]] -
            {{Cos[theta / 2], -I Sin[theta / 2]}, {-I Sin[theta / 2], Cos[theta / 2]}}] === ConstantArray[0, {2, 2}]
    ],
    True,
    TestID -> "Initialization-shorthand-preserves-exact-symbolic-operators"
]

VerificationTest[
    With[{name = SymbolName[Unique["InitializationWarning$"]], reference = QuantumOperator["H"]},
        Block[{initializationCalls = 0},
            Internal`InheritedBlock[{QuantumOperator},
                Unprotect[QuantumOperator];
                QuantumOperator::initwarningprobe = "temporary initialization test warning";
                HoldPattern[QuantumOperator[name]] :=
                    (initializationCalls++; Message[QuantumOperator::initwarningprobe]; reference);
                {$InitializationShorthand[{name}] === {reference},
                    $InitializationShorthand[{name}] === {reference}, initializationCalls}
            ]
        ]
    ],
    {True, True, 2},
    {QuantumOperator::initwarningprobe, QuantumOperator::initwarningprobe},
    TestID -> "Initialization-successful-warnings-are-not-cached"
]

VerificationTest[
    With[{name = SymbolName[Unique["InitializationAggregate$"]]},
        Block[{initializationCalls = 0},
            Internal`InheritedBlock[{QuantumOperator},
                Unprotect[QuantumOperator];
                QuantumOperator::initaggregateprobe = "temporary initialization aggregate failure";
                HoldPattern[QuantumOperator[name]] :=
                    (initializationCalls++; Message[QuantumOperator::initaggregateprobe]; $Failed);
                {$InitializationShorthand[{name}], $InitializationShorthand[{name}], initializationCalls}
            ]
        ]
    ],
    {{$Failed}, {$Failed}, 2},
    {QuantumOperator::initaggregateprobe, QuantumOperator::initaggregateprobe},
    TestID -> "Initialization-list-failures-repeat-their-diagnostics"
]

VerificationTest[
    With[{name = SymbolName[Unique["InitializationRetry$"]], reference = QuantumOperator["X"]},
        Block[{initializationCalls = 0},
            Internal`InheritedBlock[{QuantumOperator},
                Unprotect[QuantumOperator];
                HoldPattern[QuantumOperator[name]] :=
                    (initializationCalls++; If[initializationCalls === 1, Failure["InitializationRetry", <||>], reference]);
                With[{first = $InitializationShorthand[{name}], second = $InitializationShorthand[{name}],
                      third = $InitializationShorthand[{name}]},
                    {MatchQ[first, {Failure["InitializationRetry", _Association]}],
                        second === {reference}, third === {reference}, initializationCalls}
                ]
            ]
        ]
    ],
    {True, True, True, 2},
    TestID -> "Initialization-list-failure-retries-then-caches-success"
]

VerificationTest[
    With[{name = SymbolName[Unique["InitializationHeld$"]]},
        Block[{initializationCalls = 0, initializationHeldEffects = 0},
            Internal`InheritedBlock[{QuantumOperator},
                Unprotect[QuantumOperator];
                HoldPattern[QuantumOperator[name]] :=
                    (initializationCalls++; HoldComplete[initializationHeldEffects++]);
                With[{direct = QuantumOperator[name]},
                    initializationCalls = 0;
                    With[{first = $InitializationShorthand[name], second = $InitializationShorthand[name]},
                        {direct === HoldComplete[initializationHeldEffects++],
                            first === direct, second === direct, initializationCalls, initializationHeldEffects}
                    ]
                ]
            ]
        ]
    ],
    {True, True, True, 1, 0},
    TestID -> "Initialization-cache-preserves-held-results-without-side-effects"
]
