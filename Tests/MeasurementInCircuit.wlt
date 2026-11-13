BeginTestSection["MeasurementInCircuit"]

(* A measurement in a circuit follows Born's rule and keeps its record: what is done to the
   measured qubits afterwards cannot change the outcome probabilities, and inside a circuit a
   measurement gives what it gives applied directly to the state. *)

(* an X measurement of RY(t)|0> gives + with probability (1 + Sin[t]) / 2, whatever follows *)
VerificationTest[
    Chop[N @ {
        QuantumCircuitOperator[{"RY"[0.5] -> 1, QuantumMeasurementOperator["X", {1}]}][]["ProbabilitiesList"],
        QuantumCircuitOperator[{"RY"[0.5] -> 1, QuantumMeasurementOperator["X", {1}], "RY"[0.7] -> 1}][]["ProbabilitiesList"],
        QuantumCircuitOperator[{"RY"[0.5] -> 1, QuantumMeasurementOperator["X", {1}], "H" -> 1}][]["ProbabilitiesList"]
    } - ConstantArray[{1 + Sin[0.5], 1 - Sin[0.5]} / 2, 3], 10 ^ -10],
    ConstantArray[0, {3, 2}],
    TestID -> "later-gate-keeps-x-measurement-probabilities"
]

(* the same for a product-basis measurement on two wires out of order *)
VerificationTest[
    Chop[N[
        QuantumCircuitOperator[{"H" -> 1, "RY"[0.4] -> 3, "CNOT" -> {1, 2}, QuantumMeasurementOperator["XZ", {3, 1}], "H" -> 3, "RX"[0.3] -> 1}][]["ProbabilitiesList"] -
        QuantumCircuitOperator[{"H" -> 1, "RY"[0.4] -> 3, "CNOT" -> {1, 2}, QuantumMeasurementOperator["XZ", {3, 1}]}][]["ProbabilitiesList"]
    ], 10 ^ -10],
    {0, 0, 0, 0},
    TestID -> "later-gates-keep-product-basis-probabilities"
]

(* a second measurement on the same wire: the first outcome's marginal is unchanged and the
   second follows from the post-measurement eigenvector *)
VerificationTest[
    With[{p = ArrayReshape[N @ QuantumCircuitOperator[{"RY"[0.6] -> 1, QuantumMeasurementOperator["X", {1}], "RY"[0.3] -> 1, {1}}][]["ProbabilitiesList"], {2, 2}]},
        Chop[{Total[p, {2}][[1]] - (1 + Sin[0.6]) / 2, Total[p][[2]] - ((1 + Sin[0.6]) / 2 Sin[Pi / 4 + 0.15] ^ 2 + (1 - Sin[0.6]) / 2 Sin[-Pi / 4 + 0.15] ^ 2)}, 10 ^ -10]
    ],
    {0, 0},
    TestID -> "remeasured-wire-marginals"
]

(* an observable on wires {2, 1} measured in a circuit gives Born's probabilities |<v|psi>|^2,
   with v the eigenvectors of the observable written on wires {1, 2} (conjugated by the swap),
   in the order the direct measurement lists them *)
VerificationTest[
    With[{
        o = KroneckerProduct[PauliMatrix[1], PauliMatrix[1]] + KroneckerProduct[PauliMatrix[3], PauliMatrix[3]] / 2 + 3 / 10 KroneckerProduct[PauliMatrix[3], IdentityMatrix[2]],
        swap = {{1, 0, 0, 0}, {0, 0, 1, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}},
        v = Normalize[{1, 2 I, 3, -1 + I}]
    },
        With[{
            mo = QuantumMeasurementOperator[QuantumOperator[o, {2, 1}]],
            circuit = N @ QuantumCircuitOperator[{QuantumOperator[QuantumState[v], {{1, 2}, {}}], QuantumMeasurementOperator[QuantumOperator[o, {2, 1}]]}][]["ProbabilitiesList"]
        },
            Chop[N @ {
                Sort[circuit] - Sort[Abs[Conjugate[#] . v] ^ 2 & /@ Eigensystem[N[swap . o . swap]][[2]]],
                circuit - mo[QuantumState[v]]["ProbabilitiesList"]
            }, 10 ^ -10]
        ]
    ],
    {{0, 0, 0, 0}, {0, 0, 0, 0}},
    TestID -> "observable-in-circuit-follows-born"
]

VerificationTest[
    With[{mo = QuantumMeasurementOperator[QuantumOperator[-PauliMatrix[3], {1}]]},
        {N @ QuantumCircuitOperator[{QuantumOperator[QuantumState["0"], {{1}, {}}], mo}][]["ProbabilitiesList"], N @ mo[QuantumState["0"]]["ProbabilitiesList"]}
    ],
    {{0., 1.}, {0., 1.}},
    TestID -> "diagonal-observable-in-circuit-as-direct"
]

(* a state injected onto a wire in use replaces that wire's qubit: the reset of the wire *)
VerificationTest[
    With[{
        inj = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, QuantumOperator[QuantumState["+"], {{2}, {}}], "CNOT" -> {2, 1}}],
        res = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, "Reset"[QuantumState["+"]] -> {2}, "CNOT" -> {2, 1}}]
    },
        {inj[]["Dimensions"], Chop[Norm[N[Normal[inj[]["DensityMatrix"]] - Normal[res[]["DensityMatrix"]]]], 10 ^ -10]}
    ],
    {{2, 2}, 0},
    TestID -> "injection-on-wire-in-use-resets"
]

(* an injection onto a fresh wire is a new qubit, not a reset; this held before the reset rule
   too, and guards against the rule touching wires that carry nothing *)
VerificationTest[
    N @ Normal @ QuantumCircuitOperator[{QuantumOperator[QuantumState["+"], {{2}, {}}], "CNOT" -> {2, 1}}][]["StateVector"],
    N @ {1, 0, 0, 1} / Sqrt[2],
    SameTest -> (Chop[Norm[#1 - #2], 10 ^ -10] == 0 &),
    TestID -> "injection-on-fresh-wire-unchanged"
]

(* the same with symbolic angles: Born's (1 + Sin[t]) / 2, whatever rotation follows *)
VerificationTest[
    Simplify[
        QuantumCircuitOperator[{"RY"[t] -> 1, QuantumMeasurementOperator["X", {1}], "RY"[f] -> 1}][]["ProbabilitiesList"] - {1 + Sin[t], 1 - Sin[t]} / 2,
        Element[{t, f}, Reals]
    ],
    {0, 0},
    TestID -> "later-gate-keeps-symbolic-x-measurement-probabilities"
]

(* the Schrodinger route agrees *)
VerificationTest[
    Chop[N[QuantumCircuitOperator[{"RY"[0.5] -> 1, QuantumMeasurementOperator["X", {1}], "RY"[0.7] -> 1}][Method -> "Schrodinger"]["ProbabilitiesList"] - {1 + Sin[0.5], 1 - Sin[0.5]} / 2], 10 ^ -10],
    {0, 0},
    TestID -> "later-gate-keeps-x-measurement-probabilities-schrodinger"
]

(* a non-diagonal observable on the second wire of an entangled pair, against Born's rule
   computed by hand: |<v|rho|v>| for each eigenvector v of the reduced state of that wire *)
VerificationTest[
    With[{
        o = 2 PauliMatrix[1] + PauliMatrix[3],
        pre = QuantumCircuitOperator[{"RY"[0.5] -> 1, "RY"[1.2] -> 2, "CNOT" -> {1, 2}}]
    },
        With[{
            rho = Normal @ N @ QuantumPartialTrace[pre[], {1}]["DensityMatrix"],
            vecs = Eigensystem[N[o]][[2]]
        },
            Chop[
                Sort @ N @ QuantumCircuitOperator[Append[pre["Elements"], QuantumMeasurementOperator[QuantumOperator[o, {2}]]]][]["ProbabilitiesList"] -
                Sort @ Re[Conjugate[#] . rho . # & /@ vecs],
                10 ^ -10
            ]
        ]
    ],
    {0, 0},
    TestID -> "non-diagonal-observable-in-circuit-follows-born"
]

(* a degenerate observable on targets out of order: one outcome per eigenvector, as directly *)
VerificationTest[
    With[{
        mo = QuantumMeasurementOperator[QuantumOperator[KroneckerProduct[PauliMatrix[3], PauliMatrix[3]] + KroneckerProduct[PauliMatrix[1], PauliMatrix[1]], {2, 1}]],
        psi = QuantumState[Normalize[{1, 2 I, 3, -1 + I}]]
    },
        Chop[N[QuantumCircuitOperator[{QuantumOperator[psi, {{1, 2}, {}}], mo}][]["ProbabilitiesList"] - mo[psi]["ProbabilitiesList"]], 10 ^ -10]
    ],
    {0, 0, 0, 0},
    TestID -> "degenerate-observable-in-circuit-as-direct"
]

(* an injection onto one wire in use and one fresh wire: the wire in use is reset, the fresh
   one gains a qubit; compared with tracing wire 2 out of the Bell pair by hand *)
VerificationTest[
    With[{s = QuantumState[Normalize[{1, 2, 0, 3 I}]]}, With[{qc = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, QuantumOperator[s, {{2, 3}, {}}], "CNOT" -> {3, 1}}]},
        With[{ref = Normal @ N @ QuantumOperator["CNOT", {3, 1}][QuantumTensorProduct[QuantumPartialTrace[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}][], {2}], s]]["DensityMatrix"]},
            {
                Through[{qc[], qc[Method -> "Schrodinger"], qc["QuantumOperator"][]}["Dimensions"]],
                Chop[Norm[Normal @ N @ #["DensityMatrix"] - ref] & /@ {qc[], qc[Method -> "Schrodinger"], qc["QuantumOperator"][]}, 10 ^ -10]
            }
        ]
    ]],
    {ConstantArray[{2, 2, 2}, 3], {0, 0, 0}},
    TestID -> "injection-partly-on-wire-in-use-resets-that-wire"
]

(* with an input state, through every route: the injection resets the input's wire *)
VerificationTest[
    With[{
        inj = QuantumCircuitOperator[{"CNOT" -> {1, 2}, QuantumOperator[QuantumState["+"], {{2}, {}}], "CNOT" -> {2, 1}}],
        res = QuantumCircuitOperator[{"CNOT" -> {1, 2}, "Reset"[QuantumState["+"]] -> {2}, "CNOT" -> {2, 1}}],
        in = QuantumState[Normalize[{1, 2, 0, 1}]]
    },
        Chop[Norm[Normal @ N @ #["DensityMatrix"] - Normal @ N @ res[in]["DensityMatrix"]] & /@ {inj[in], inj[in, Method -> "Schrodinger"], inj["QuantumOperator"][in]}, 10 ^ -10]
    ],
    {0, 0, 0},
    TestID -> "injection-on-wire-in-use-resets-input-state"
]

(* a second injection on the same wire replaces the first: only the last state is left *)
VerificationTest[
    With[{qc = QuantumCircuitOperator[{"H" -> 1, QuantumOperator[QuantumState["+"], {{1}, {}}], QuantumOperator[QuantumState["1"], {{1}, {}}]}]},
        {qc[]["Dimensions"], Chop[N @ Normal @ qc[]["StateVector"]]}
    ],
    {{2}, {0, 1.}},
    TestID -> "second-injection-on-a-wire-replaces-the-first"
]

(* an injection onto a measured wire leaves the outcome probabilities as they were, and the
   wire holds the injected |1> whatever the outcome: one register qudit, one system qubit *)
VerificationTest[
    With[{r = QuantumCircuitOperator[{"RY"[0.6] -> 1, QuantumMeasurementOperator["X", {1}], QuantumOperator[QuantumState["1"], {{1}, {}}]}][]},
        {
            r["Dimensions"],
            Chop[N[r["ProbabilitiesList"] - {1 + Sin[0.6], 1 - Sin[0.6]} / 2], 10 ^ -10],
            Chop[N @ Normal @ QuantumPartialTrace[r["State"], {1}]["DensityMatrix"], 10 ^ -10]
        }
    ],
    {{2, 2}, {0, 0}, {{0, 0}, {0, 1.}}},
    TestID -> "injection-after-measurement-resets-the-measured-wire"
]

(* the fresh wire lies below the wire in use: |1> |0> on wires 1, 2, then CNOT gives |11> *)
VerificationTest[
    With[{qc = QuantumCircuitOperator[{"H" -> 2, QuantumOperator[QuantumState["10"], {{1, 2}, {}}], "CNOT" -> {1, 2}}]},
        Chop[N @ Normal @ #["StateVector"]] & /@ {qc[], qc["QuantumOperator"][]}
    ],
    {{0, 0, 0, 1.}, {0, 0, 0, 1.}},
    TestID -> "injection-with-fresh-wire-below-wire-in-use"
]

(* a wire ended by a post-selection carries nothing, so a later injection there is a new qubit:
   <0| on wire 2 of (|00> + |11>) / Sqrt[2] leaves |0> / Sqrt[2], then |+> and CNOT(2 -> 1) give
   (|00> + |11>) / 2. This held before the reset rule too, and guards against the rule treating
   an ended wire as one in use. *)
VerificationTest[
    With[{qc = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, QuantumOperator[QuantumState["0"]["Dual"], {{}, {2}}], QuantumOperator[QuantumState["+"], {{2}, {}}], "CNOT" -> {2, 1}}]},
        {qc[]["Dimensions"], Chop[N @ Normal @ qc[]["StateVector"]]}
    ],
    {{2, 2}, {0.5, 0, 0, 0.5}},
    TestID -> "injection-after-post-selection-is-a-new-qubit"
]

(* an injection of another dimension than the qudit it displaces: the qutrit is traced out and
   the qubit left on its wire, and the other way round (wire 1 of the Bell pair is left mixed) *)
VerificationTest[
    With[{
        qd = QuantumCircuitOperator[{QuantumOperator[{{0, 0, 1}, {1, 0, 0}, {0, 1, 0}}, {1}], QuantumOperator[QuantumState[{0, 1}], {{1}, {}}]}],
        qu = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, QuantumOperator[QuantumState[{0, 0, 1}], {{2}, {}}]}]
    },
        {
            qd[QuantumState[{1, 0, 0}]]["Dimensions"],
            Chop[N @ Normal @ qd[QuantumState[{1, 0, 0}]]["DensityMatrix"]],
            Chop[N @ Normal @ qd["QuantumOperator"][QuantumState[{1, 0, 0}]]["DensityMatrix"]],
            qu[]["Dimensions"],
            Chop[N @ Normal @ qu[]["DensityMatrix"]]
        }
    ],
    {
        {2}, {{0, 0}, {0, 1.}}, {{0, 0}, {0, 1.}},
        {2, 3}, N @ KroneckerProduct[IdentityMatrix[2] / 2, DiagonalMatrix[{0, 0, 1}]]
    },
    SameTest -> (Chop[N[Flatten[{#1}] - Flatten[{#2}]]] === ConstantArray[0, Length[Flatten[{#2}]]] &),
    TestID -> "injection-of-another-dimension-replaces-the-qudit"
]

EndTestSection[]
