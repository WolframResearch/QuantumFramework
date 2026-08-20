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

(* SPSRGradientValues builds MatrixExp[I phi QuantumOperator[pauli]] for its shift gate; a bare
   gate-name string in that slot degrades to Power[E, ...] and exceeds $RecursionLimit in
   the circuit-shorthand parser. For a commuting generator G = theta X ox I the two split
   evolutions merge for every sample, so each random s draw yields the same difference and
   the estimator equals the exact derivative of <Z1> = Cos[2 theta], namely -2 Sin[2 theta]. *)
VerificationTest[
    Module[{gen},
        gen[t1_] = Normal @ QuantumOperator[t1*QuantumOperator[{"X" -> {1}, "I" -> {2}}], "Parameters" -> {t1}]["Matrix"];
        SPSRGradientValues[gen, "PauliX", "ParameterValues" -> {0.3}, "RandomNumberCount" -> 3]
    ],
    {{0.3, -2 Sin[0.6]}},
    SameTest -> (Norm[Flatten[#1 - #2]] < 1*^-6 &),
    TestID -> "Opt-SPSR-commuting-generator-exact-derivative"
]

(* The tutorial's own generator (non-commuting theta and drift parts) on a reduced sweep of
   exact parameter values, the same shape the default Subdivide[0, 2 Pi, 50] sweep produces:
   the estimator must come back as numeric-theta/real-value pairs with no messages. *)
VerificationTest[
    Module[{gen},
        gen[t1_] = Normal @ QuantumOperator[
            t1*QuantumOperator[{"X" -> {1}, "I" -> {2}}] -
              0.15*QuantumOperator[{"Z" -> {1}, "X" -> {2}}] +
              1.6*QuantumOperator[{"I" -> {1}, "X" -> {2}}],
            "Parameters" -> {t1}]["Matrix"];
        MatchQ[
            SPSRGradientValues[gen, "PauliX", "ParameterValues" -> Subdivide[0, 2 Pi, 2], "RandomNumberCount" -> 2],
            {{_?NumericQ, _Real} ..}
        ]
    ],
    True,
    TestID -> "Opt-SPSR-tutorial-generator-runs"
]

EndTestSection[]
