(* f[qs], f[qb] and f[qmo] for a head f with the NumericFunction attribute act on the
   density matrix of a state, the matrix of a basis and the operator of a measurement;
   any other head leaves the object alone. Each rule is an upvalue, tried on every
   expression with such an object among its arguments, including a List of them. It
   must reject a head that is not a NumericFunction before it reads the arguments, or a
   list of N states or bases costs about N^3 to evaluate.

   References are built by hand: |psi><psi| from the amplitudes, the spectral form
   f(rho) = f((1 + r)/2) P+ + f((1 - r)/2) P- of rho = (1 + r X)/2 with P+- = (1 +- X)/2,
   and the conjugate of the amplitudes or basis vectors.

   Regimes covered:
     general symbolic     Sqrt and Log of (1 + r X)/2 for a symbol r
     exactly solvable     Sqrt of a pure state is the state, for |1>, |->, |-i> and a
                          Bell state; Conjugate and 2 qs of a vector state; Conjugate of
                          a basis; 2 qmo
     edge                 a head that is not a NumericFunction, and a List, leave the
                          objects untouched; the basis and target are kept
     scaling              lists of 20000 states, bases and measurement operators *)

dm[qs_] := Normal[qs["DensityMatrix"]]

zeroQ[m_] := AllTrue[Flatten[FullSimplify[m]], # === 0 &]

$pure = <|
    "1" -> {0, 1},
    "-" -> {1, -1} / Sqrt[2],
    "-i" -> {1, -I} / Sqrt[2],
    "Bell" -> {1, 0, 0, 1} / Sqrt[2]
|>;

KeyValueMap[
    Function[{name, v},
        VerificationTest[
            zeroQ[dm[Sqrt[QuantumState[v]]] - KroneckerProduct[v, Conjugate[v]]],
            True,
            TestID -> "StateNumericFunctionHead-Sqrt-pure-" <> name
        ]
    ],
    $pure
];

With[{
    rho = (IdentityMatrix[2] + r PauliMatrix[1]) / 2,
    pp = (IdentityMatrix[2] + PauliMatrix[1]) / 2,
    pm = (IdentityMatrix[2] - PauliMatrix[1]) / 2
},
    Scan[
        Function[f,
            VerificationTest[
                zeroQ[dm[f[QuantumState[rho]]] - (f[(1 + r) / 2] pp + f[(1 - r) / 2] pm)],
                True,
                TestID -> "StateNumericFunctionHead-symbolic-" <> ToString[f]
            ]
        ],
        {Sqrt, Log}
    ]
]

VerificationTest[
    With[{qs = QuantumState[$pure["-i"]]},
        {Normal[Conjugate[qs]["StateVector"]], Normal[(2 qs)["StateVector"]], Conjugate[qs]["Basis"] === qs["Basis"]}
    ],
    {{1, I} / Sqrt[2], {2, -2 I} / Sqrt[2], True},
    TestID -> "StateNumericFunctionHead-Conjugate-and-Times-vector"
]

VerificationTest[
    With[{m = {{1, I}, {1, -I}} / Sqrt[2]},
        With[{qb = QuditBasis[m]}, {Normal[Conjugate[qb]["Matrix"]], Conjugate[qb]["Names"] === qb["Names"]}]
    ],
    {{{1, 1}, {-I, I}} / Sqrt[2], True},
    TestID -> "StateNumericFunctionHead-Conjugate-basis"
]

VerificationTest[
    With[{qmo = QuantumMeasurementOperator["Z"]},
        {Normal[(2 qmo)["Operator"]["Matrix"]] === 2 Normal[qmo["Operator"]["Matrix"]], (2 qmo)["Target"]}
    ],
    {True, {1}},
    TestID -> "StateNumericFunctionHead-Times-measurement"
]

VerificationTest[
    With[{qs = QuantumState[$pure["-"]], qb = QuditBasis[2]},
        {Head[notNumericF[qs]], First[notNumericF[qs]] === qs, Head[notNumericF[qb]], First[{qs, qs}] === qs, First[{qb, qb}] === qb}
    ],
    {notNumericF, True, notNumericF, True, True},
    TestID -> "StateNumericFunctionHead-other-head-untouched"
]

(* Tens of microseconds once a List is rejected at its head; hours when the rule scans
   the argument sequence. *)
VerificationTest[
    TimeConstrained[Length /@ {
        ConstantArray[QuantumState[$pure["-"]], 20000],
        ConstantArray[QuditBasis[{{1, I}, {1, -I}} / Sqrt[2]], 20000],
        ConstantArray[QuantumMeasurementOperator["Z"], 20000]
    }, 10, $Aborted],
    {20000, 20000, 20000},
    TestID -> "StateNumericFunctionHead-long-lists"
]
