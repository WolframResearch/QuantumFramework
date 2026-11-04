(* A measurement operator holds the operators of all its outcomes in one operator, with an
   extra wire for the outcome, so a function f of that operator is not f of each outcome's
   operator, and its matrix is not square. Multiplying by a number scales every outcome's
   operator and is kept; every other head with the NumericFunction attribute gives one
   QuantumMeasurementOperator::numfunc message and a Failure. The square root of each
   effect, the one function of a measurement in use, is what the constructor from a list of
   effects already takes.

   References are built by hand: twice and a times the matrix of the Z measurement, and the
   square roots of diagonal effects diag(p, q) and diag(1 - p, 1 - q).

   Regimes covered:
     refused              Sqrt, Sin, Exp, Log, Conjugate, a power and adding a number
     kept                 2 qmo and a qmo for a symbol a, against the matrix by hand; the
                          sum of two measurement operators
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
    {{Sqrt, Sqrt}, {Sin, Sin}, {Exp, Exp}, {Log, Log}, {Conjugate, Conjugate}, {#^2 &, Power}, {1 + # &, Plus}}
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
