(* A measurement operator holds the operators of all its outcomes in one operator, with an
   extra wire for the outcome, so a function f of that operator is not f of each outcome's
   operator, and its matrix is not square. Multiplying by a number scales every outcome's
   operator, and Conjugate conjugates each one, so both are kept; every other head with the
   NumericFunction attribute gives one QuantumMeasurementOperator::numfunc message and a
   Failure. The square root of each effect, the one function of a measurement in use, is
   what the constructor from a list of effects already takes.

   References are built by hand: twice and a times the matrix of the Z measurement, the
   projectors (1 +- Y)/2, and the square roots of diagonal effects diag(p, q) and
   diag(1 - p, 1 - q).

   Regimes covered:
     refused              Sqrt, Sin, Exp, Log, a power and adding a number
     kept                 2 qmo and a qmo for a symbol a, against the matrix by hand; the
                          sum of two measurement operators; Conjugate of the measurement
                          with effects (1 +- Y)/2, which conjugate to (1 -+ Y)/2
     general symbolic     effects diag(p, q), diag(1 - p, 1 - q) for symbols p and q: the
                          stored operators are their square roots, and |0> gives the
                          probabilities {p, 1 - p} *)

$z = QuantumMeasurementOperator["Z"];

zMatrix = Normal[$z["Operator"]["Matrix"]];

Scan[
    Function[{fh},
        VerificationTest[
            First[fh][$z],
            Failure["UndefinedFunction", <|"MessageTemplate" :> QuantumMeasurementOperator::numfunc, "MessageParameters" -> {Last[fh]}|>],
            {QuantumMeasurementOperator::numfunc},
            TestID -> "MeasurementNumericFunction-refused-" <> ToString[Last[fh]]
        ]
    ],
    {{Sqrt, Sqrt}, {Sin, Sin}, {Exp, Exp}, {Log, Log}, {#^2 &, Power}, {1 + # &, Plus}}
]

VerificationTest[
    {
        Normal[(2 $z)["Operator"]["Matrix"]] === 2 zMatrix,
        Normal[(a $z)["Operator"]["Matrix"]] === a zMatrix,
        (a $z)["Target"],
        Head[$z + $z]
    },
    {True, True, {1}, QuantumMeasurementOperator},
    TestID -> "MeasurementNumericFunction-times-and-sum-kept"
]

(* conjugating the effects (1 +- Y)/2 swaps them, so |+i> goes from outcome 1 to outcome 2 *)
VerificationTest[
    With[{py = {(IdentityMatrix[2] + PauliMatrix[2]) / 2, (IdentityMatrix[2] - PauliMatrix[2]) / 2}, plusI = QuantumState[{1, I} / Sqrt[2]]},
        With[{qmo = QuantumMeasurementOperator[py]},
            With[{c = Conjugate[qmo]},
                {
                    Head[c],
                    Normal[c["Operator"]["Matrix"]] === Conjugate[Normal[qmo["Operator"]["Matrix"]]],
                    Normal /@ c["POVMElements"] === Reverse[py],
                    c["Target"],
                    {qmo[plusI]["ProbabilitiesList"], c[plusI]["ProbabilitiesList"]}
                }
            ]
        ]
    ],
    {QuantumMeasurementOperator, True, True, {1}, {{1, 0}, {0, 1}}},
    TestID -> "MeasurementNumericFunction-Conjugate-kept"
]

VerificationTest[
    With[{qmo = QuantumMeasurementOperator[{DiagonalMatrix[{p, q}], DiagonalMatrix[{1 - p, 1 - q}]}]},
        {
            Normal /@ qmo["POVMElements"],
            Values[Normal[#["Matrix"]] & /@ qmo["Operators"]],
            Simplify[qmo[QuantumState["0"]]["ProbabilitiesList"], 0 < p < 1]
        }
    ],
    {
        {DiagonalMatrix[{p, q}], DiagonalMatrix[{1 - p, 1 - q}]},
        {DiagonalMatrix[{Sqrt[p], Sqrt[q]}], DiagonalMatrix[{Sqrt[1 - p], Sqrt[1 - q]}]},
        {p, 1 - p}
    },
    TestID -> "MeasurementNumericFunction-effects-take-square-root"
]
