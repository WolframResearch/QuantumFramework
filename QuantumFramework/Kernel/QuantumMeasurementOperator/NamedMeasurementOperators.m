Package["Wolfram`QuantumFramework`"]


QuantumMeasurementOperator["RandomHermitian", args___, target : (_ ? orderQ) : {1}] := With[{
    basis = QuantumBasis[args]
},
    QuantumMeasurementOperator[
        QuantumOperator[
            With[{m = RandomComplex[1 + I, {basis["Dimension"], basis["Dimension"]}]}, (m + ConjugateTranspose[m]) / 2],
            basis
        ],
        target
]
]

