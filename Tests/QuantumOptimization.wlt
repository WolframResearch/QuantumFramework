(* Regression tests for the QuantumOptimization context, which had none. *)

Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`QuantumOptimization`"]

BeginTestSection["QuantumOptimization"]

(* The Fubini-Study metric of RY(theta)|0> is the textbook constant 1/4, and
   the state and vector forms must agree on it. *)
VerificationTest[
    With[{ps = QuantumCircuitOperator[{"RY"[\[FormalTheta]]}][QuantumState["0"]]},
        {
            Simplify @ Normal @ FubiniStudyMetricTensor[ps, "Matrix", "Parameters" -> {\[FormalTheta]}],
            Simplify @ Normal @ FubiniStudyMetricTensor[ps["StateVector"], "Matrix", "Parameters" -> {\[FormalTheta]}]
        }
    ],
    {{{1/4}}, {{1/4}}},
    TestID -> "Opt-Fubini-Study-metric-of-a-rotation"
]

VerificationTest[
    Chop[Last @ GradientDescent[Function[p, p^2], {1.5}, "MaxIterations" -> 40], 1*^-2],
    {0},
    TestID -> "Opt-gradient-descent-converges"
]

EndTestSection[]
