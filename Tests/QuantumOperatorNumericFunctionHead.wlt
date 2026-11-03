(* f[qo] for a head f with the NumericFunction attribute is the matrix function of the
   stored matrix; any other head leaves the operator alone. The rule that does this is
   an upvalue of QuantumOperator, so the kernel tries it on every expression with an
   operator among its arguments, including a List of operators. It must reject a head
   that is not a NumericFunction before it reads the arguments: a rule that matches the
   argument sequence first scans the ways to split it around one operator, and a list
   of N operators then costs about N^3 to evaluate, which is what building a circuit
   does.

   References are independent of the rule: for an involution P (P^2 = 1, as for a Pauli
   operator or a product of them), f(P) = (f(1) + f(-1))/2 + (f(1) - f(-1))/2 P, and
   f(theta P) the same with f(+-theta); for a symbolic rotation, MatrixFunction of the
   matrix written out by hand.

   Regimes covered:
     general symbolic     Sin, Cos of theta P for X, Y, Z, YY; Sin of R_X(theta)
     exactly solvable     Sin, Cos, Tanh, Sqrt, Log, Exp of X, Y, Z, YY, at the
                          principal branch for Sqrt and Log at -1
     arguments            Log[2, P] (a left argument) and P^(1/2) (a right one)
     edge                 a head that is not a NumericFunction, and a List, leave
                          the operators untouched
     scaling              a list of 20000 copies of one operator, and circuits of
                          4000 and 8000 gates, whose build time must grow linearly *)

involutionF[f_, p_] := (f[1] + f[-1]) / 2 IdentityMatrix[Length[p]] + (f[1] - f[-1]) / 2 p

involutionF[f_, p_, t_] := (f[t] + f[-t]) / 2 IdentityMatrix[Length[p]] + (f[t] - f[-t]) / 2 p

cm[qo_] := Normal[qo["Computational"]["Matrix"]]

zeroQ[m_] := AllTrue[Flatten[FullSimplify[m]], # === 0 &]

$involutions = <|
    "X" -> {QuantumOperator["X"], PauliMatrix[1]},
    "Y" -> {QuantumOperator["Y"], PauliMatrix[2]},
    "Z" -> {QuantumOperator["Z"], PauliMatrix[3]},
    "YY" -> {QuantumOperator["Y"] @ QuantumOperator["Y", {2}], KroneckerProduct[PauliMatrix[2], PauliMatrix[2]]}
|>;

KeyValueMap[
    Function[{name, op},
        Scan[
            Function[f,
                VerificationTest[
                    zeroQ[cm[f[op[[1]]]] - involutionF[f, op[[2]]]],
                    True,
                    TestID -> "NumericFunctionHead-" <> ToString[f] <> "-" <> name
                ]
            ],
            {Sin, Cos, Tanh, Sqrt, Log, Exp}
        ]
    ],
    $involutions
];

KeyValueMap[
    Function[{name, op},
        Scan[
            Function[f,
                VerificationTest[
                    zeroQ[cm[f[\[Theta] op[[1]]]] - involutionF[f, op[[2]], \[Theta]]],
                    True,
                    TestID -> "NumericFunctionHead-symbolic-" <> ToString[f] <> "-" <> name
                ]
            ],
            {Sin, Cos}
        ]
    ],
    $involutions
];

VerificationTest[
    zeroQ[
        cm[Sin[QuantumOperator["RX"[\[Theta]]]]] -
            MatrixFunction[Sin, {{Cos[\[Theta] / 2], -I Sin[\[Theta] / 2]}, {-I Sin[\[Theta] / 2], Cos[\[Theta] / 2]}}]
    ],
    True,
    TestID -> "NumericFunctionHead-symbolic-Sin-RX"
]

VerificationTest[
    zeroQ[cm[Log[2, QuantumOperator["Y"]]] - involutionF[Log, PauliMatrix[2]] / Log[2]],
    True,
    TestID -> "NumericFunctionHead-left-argument-Log2-Y"
]

VerificationTest[
    zeroQ[cm[QuantumOperator["Y"] ^ (1 / 2)] - involutionF[Sqrt, PauliMatrix[2]]],
    True,
    TestID -> "NumericFunctionHead-right-argument-Power-Y"
]

VerificationTest[
    With[{qo = QuantumOperator["Y", {2}]}, {Sin[qo]["Order"], Sin[qo]["Label"]}],
    {{{2}, {2}}, Sin["Y"]},
    TestID -> "NumericFunctionHead-keeps-order-and-label"
]

VerificationTest[
    With[{qo = QuantumOperator["Y"]}, {Head[notNumericF[qo]], First[notNumericF[qo]] === qo}],
    {notNumericF, True},
    TestID -> "NumericFunctionHead-other-head-untouched"
]

VerificationTest[
    With[{qo = QuantumOperator["Y"]}, AllTrue[{qo, qo, 2 qo}, QuantumOperatorQ] && First[{qo, qo}] === qo],
    True,
    TestID -> "NumericFunctionHead-list-untouched"
]

(* Tens of microseconds once a List is rejected at its head; hours when the rule scans
   the argument sequence. *)
VerificationTest[
    With[{qo = QuantumOperator["Y"]}, TimeConstrained[Length[ConstantArray[qo, 20000]], 10, $Aborted]],
    20000,
    TestID -> "NumericFunctionHead-long-operator-list"
]

(* Doubling the gate count doubles the build time when each gate costs the same; the
   N^3 scan multiplies it by 8. A build that hits the cap takes infinite time, so the
   ratio is Indeterminate and the test fails in seconds rather than hours. *)
VerificationTest[
    With[
        {
            build = Function[n,
                Min @ Table[
                    ClearSystemCache[];
                    Replace[
                        AbsoluteTiming[TimeConstrained[QuantumCircuitOperator[Table["H" -> {Mod[k, 10] + 1}, {k, n}]], 20, $Aborted]],
                        {{_, $Aborted} -> Infinity, {t_, _} :> t}
                    ],
                    3
                ]
            ]
        },
        TrueQ[build[8000] / build[4000] < 4]
    ],
    True,
    TestID -> "NumericFunctionHead-circuit-build-linear"
]
