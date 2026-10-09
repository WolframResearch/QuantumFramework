(* Regression tests for QuantumInterferometer: the Reck and Clements meshes must
   rebuild their unitary, and the photon statistics computed from the transfer
   matrix must agree with the Fock-space circuit of BeamSplitterOperator and
   PhaseShiftOperator. *)

Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"]

BeginTestSection["QuantumInterferometer"]

$exactMethods = {"Reck", "Clements"};

$methods = Append[$exactMethods, "CosineSine"];

rebuild[qi_] := QuantumInterferometer[qi["Elements"], qi["Modes"]]["Unitary"]

fockIndex[occupation_, d_] := FromDigits[occupation, d] + 1

(* W[[j, k]] = <1_j| U |1_k>, read off a Fock-space operator one photon at a time *)
transferMatrix[op_, m_, d_] := Transpose @ Map[
    Normal[op[FockState[#, d]]["StateVector"][[fockIndex[#, d] & /@ IdentityMatrix[m]]]] &,
    IdentityMatrix[m]
]


(* The single-photon blocks are those of the Fock-space operators. *)
VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[\[Theta], \[Phi]] -> {1, 2}, "PS"[\[Alpha]] -> 2}, 2]},
        FullSimplify[
            qi["Unitary"] == transferMatrix[qi["CircuitOperator", 2, Method -> "Recurrence"], 2, 2],
            (\[Theta] | \[Phi] | \[Alpha]) \[Element] Reals
        ]
    ],
    True,
    TestID -> "QI-element-blocks-match-Fock-operators"
]

VerificationTest[
    SeedRandom[11];
    Table[
        With[{w = RandomVariate[CircularUnitaryMatrixDistribution[m]]},
            Max[Abs[rebuild[QuantumInterferometer[w, Method -> #]] - w]] < 10^-10 & /@ $methods
        ],
        {m, 2, 8}
    ],
    ConstantArray[True, {7, 3}],
    TestID -> "QI-meshes-rebuild-Haar-unitaries"
]

VerificationTest[
    Table[
        With[{qi = QuantumInterferometer[FourierMatrix[m], Method -> method]},
            Precision[qi["Elements"]] === Infinity && AllTrue[Flatten[rebuild[qi] - FourierMatrix[m]], PossibleZeroQ]
        ],
        {m, 3, 4}, {method, $exactMethods}
    ],
    ConstantArray[True, {2, 2}],
    TestID -> "QI-exact-Fourier-decomposition"
]

(* The identities that move every phase of the mesh to its output, and keep theta in [0, Pi/2]. *)
VerificationTest[
    With[{
        bs = {{Cos[#1], - Exp[- I #2] Sin[#1]}, {Exp[I #2] Sin[#1], Cos[#1]}} &,
        phases = DiagonalMatrix[{Exp[I #1], Exp[I #2]}] &
    },
        FullSimplify[{
            bs[\[Theta], \[Gamma]] . phases[\[Alpha], \[Beta]] == phases[\[Alpha], \[Beta]] . bs[\[Theta], \[Gamma] + \[Alpha] - \[Beta]],
            bs[- \[Theta], \[Gamma]] == bs[\[Theta], \[Gamma] + Pi]
        }]
    ],
    {True, True},
    TestID -> "QI-phase-absorption-identities"
]

(* Every mesh element is a beam splitter with theta in [0, Pi/2], followed by one output phase column. *)
VerificationTest[
    SeedRandom[20];
    With[{w = RandomVariate[CircularUnitaryMatrixDistribution[6]]},
        Table[
            With[{elements = Join[qi["Elements"], qi["Dagger"]["Elements"]]},
                {
                    qi["BeamSplitterCount"] == 15,
                    qi["PhaseShifterCount"] <= 6,
                    AllTrue[Cases[elements, ("BS"[t_, _] -> _) :> t], 0 <= # <= Pi / 2 &],
                    MatchQ[qi["Elements"], {("BS"[__] -> _) .., ("PS"[_] -> _) ...}]
                }
            ],
            {qi, QuantumInterferometer[w, Method -> #] & /@ $methods}
        ]
    ],
    ConstantArray[True, {3, 4}],
    TestID -> "QI-mesh-is-beam-splitters-then-phases"
]

VerificationTest[
    SeedRandom[12];
    With[{w = RandomVariate[CircularUnitaryMatrixDistribution[6]]},
        {#["BeamSplitterCount"], #["Depth"]} & /@ (QuantumInterferometer[w, Method -> #] & /@ $exactMethods)
    ],
    {{15, 9}, {15, 6}},
    TestID -> "QI-mesh-sizes"
]

(* Degenerate nulling steps: zeros already in place, or zeros next to the entry. *)
VerificationTest[
    Table[
        With[{qi = QuantumInterferometer[w, Method -> method]}, AllTrue[Flatten[rebuild[qi] - w], PossibleZeroQ]],
        {w, {
            IdentityMatrix[4],
            PermutationMatrix[{3, 1, 4, 2}],
            DiagonalMatrix[Exp[I {1, 2, 3, 4}]],
            FourierMatrix[2]
        }},
        {method, $exactMethods}
    ],
    ConstantArray[True, {4, 2}],
    TestID -> "QI-degenerate-unitaries"
]

(* The cosine-sine mesh is numeric; these have every t_i at 0 or Pi/2, so the singular values are degenerate. *)
VerificationTest[
    Table[
        Max[Abs[rebuild[QuantumInterferometer[w, Method -> "CosineSine"]] - w]] < 10^-10,
        {w, {
            IdentityMatrix[5],
            PermutationMatrix[{3, 1, 4, 2}],
            PermutationMatrix[{3, 4, 1, 2}],
            PermutationMatrix[{4, 5, 1, 2, 3}],
            DiagonalMatrix[Exp[I {1, 2, 3, 4}]],
            FourierMatrix[2],
            FourierMatrix[7]
        }}
    ],
    ConstantArray[True, 7],
    TestID -> "QI-cosine-sine-degenerate-unitaries"
]

VerificationTest[
    QuantumInterferometer[IdentityMatrix[5]]["Elements"],
    {},
    TestID -> "QI-identity-has-no-elements"
]

VerificationTest[
    SeedRandom[13];
    With[{qi = QuantumInterferometer["Random"[3]]},
        With[{amplitudes = qi[{1, 1, 1}], out = qi["CircuitOperator", 4][FockState[{1, 1, 1}, 4]]["StateVector"]},
            Max[Abs[Lookup[amplitudes, Ket[#], 0] - out[[fockIndex[#, 4]]]] & /@ Tuples[Range[0, 3], 3]] < 10^-10
        ]
    ],
    True,
    TestID -> "QI-photon-sector-matches-Fock-circuit"
]

VerificationTest[
    SeedRandom[14];
    With[{qi = QuantumInterferometer["Random"[4]]},
        Chop[Normal[qi[FockState[{2, 0, 1, 0}, 4]]["StateVector"]] -
            Normal[qi["CircuitOperator", 4][FockState[{2, 0, 1, 0}, 4]]["StateVector"]]]
    ],
    ConstantArray[0, 4^4],
    TestID -> "QI-applies-to-Fock-states"
]

VerificationTest[
    SeedRandom[17];
    With[{qi = QuantumInterferometer["Random"[3]], in = (FockState[{2, 0, 0}, 3] + I FockState[{0, 1, 1}, 3]) / Sqrt[2]},
        Chop[Normal[qi[in]["StateVector"] - qi["CircuitOperator", 3][in]["StateVector"]]]
    ],
    ConstantArray[0, 3^3],
    TestID -> "QI-applies-to-superpositions"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[Pi / 4, 0] -> {1, 2}}, 2]},
        {
            qi["Amplitude", {1, 1} -> {1, 1}],
            qi["Probability", Ket[{1, 1}] -> Ket[{2, 0}]],
            qi["Probabilities", {1, 1}],
            qi[Ket[{1, 1}]]
        }
    ],
    {
        0,
        1/2,
        <|Ket[{2, 0}] -> 1/2, Ket[{1, 1}] -> 0, Ket[{0, 2}] -> 1/2|>,
        <|Ket[{2, 0}] -> - 1/Sqrt[2], Ket[{0, 2}] -> 1/Sqrt[2]|>
    },
    TestID -> "QI-Hong-Ou-Mandel"
]

(* One photon per port of the Fourier tritter: the six patterns whose mode labels
   do not sum to a multiple of 3 are dark. *)
VerificationTest[
    Keys @ Select[QuantumInterferometer[FourierMatrix[3]]["Probabilities", {1, 1, 1}], PossibleZeroQ],
    Ket /@ {{2, 1, 0}, {2, 0, 1}, {1, 2, 0}, {1, 0, 2}, {0, 2, 1}, {0, 1, 2}},
    TestID -> "QI-Fourier-tritter-suppression-law"
]

VerificationTest[
    SeedRandom[15];
    With[{qi = QuantumInterferometer["Random"[5]]},
        {
            Chop[Total[qi["Probabilities", {1, 1, 0, 1, 0}]] - 1],
            Chop[qi["Probability", {1, 1, 0, 1, 0} -> {0, 0, 2, 0, 1}] - qi["Probabilities", {1, 1, 0, 1, 0}][Ket[{0, 0, 2, 0, 1}]]]
        }
    ],
    {0, 0},
    TestID -> "QI-probabilities-normalized"
]

VerificationTest[
    QuantumInterferometer["Fourier"[3]]["Probabilities", {1, 1}],
    $Failed,
    {QuantumInterferometer::occ},
    TestID -> "QI-rejects-wrong-pattern-length"
]

VerificationTest[
    SeedRandom[16];
    With[{qi1 = QuantumInterferometer["Random"[4]], qi2 = QuantumInterferometer["Random"[4], Method -> "Reck"]},
        Max[Abs[{
            qi2[qi1]["Unitary"] - qi2["Unitary"] . qi1["Unitary"],
            rebuild[qi2[qi1]] - qi2["Unitary"] . qi1["Unitary"],
            rebuild[qi1["Dagger"]] - ConjugateTranspose[qi1["Unitary"]]
        }]] < 10^-10
    ],
    True,
    TestID -> "QI-composition-and-dagger"
]

VerificationTest[
    Head /@ {QuantumInterferometer["Random"[5]]["Diagram", "ShowParameters" -> True], QuantumInterferometer["Fourier"[3]]["Diagram"]},
    {Graphics, Graphics},
    TestID -> "QI-diagram"
]

VerificationTest[
    QuantumInterferometer[{{1, 1}, {1, 1}}],
    $Failed,
    {QuantumInterferometer::nonunitary},
    TestID -> "QI-rejects-nonunitary"
]

VerificationTest[
    QuantumInterferometer[{{a, b}, {c, d}}],
    $Failed,
    {QuantumInterferometer::symbolic},
    TestID -> "QI-rejects-symbolic-matrix"
]

(* with 2 levels the bunched outputs |2,0> and |0,2> do not fit and are dropped *)
VerificationTest[
    QuantumInterferometer[{"BS"[Pi / 4, 0] -> {1, 2}}, 2][FockState[{1, 1}, 2]]["Norm"],
    0,
    {QuantumInterferometer::cutoff},
    TestID -> "QI-warns-on-short-cutoff"
]

VerificationTest[
    QuantumInterferometer["Fourier"[3]]["CircuitOperator"],
    $Failed,
    {QuantumInterferometer::levels},
    TestID -> "QI-circuit-needs-levels"
]

(* the exact decomposition of FourierMatrix[5] stalls, so the named form is numeric *)
VerificationTest[
    With[{qi = TimeConstrained[QuantumInterferometer["Fourier"[8]], 10]},
        {Precision[qi["Unitary"]], Max[Abs[rebuild[qi] - FourierMatrix[8]]] < 10^-10}
    ],
    {MachinePrecision, True},
    TestID -> "QI-named-Fourier-is-numeric"
]

VerificationTest[
    QuantumInterferometer[{"BS"[\[Theta], \[Phi]] -> {1, 2}}, 2]["MatrixPlot"],
    $Failed,
    {QuantumInterferometer::plot},
    TestID -> "QI-matrix-plot-needs-numbers"
]

(* consecutive phases on one mode get their own positions, also as trailing phases and after "Dagger" *)
VerificationTest[
    With[{qi = QuantumInterferometer[{"PS"[a] -> 1, "PS"[b] -> 1, "BS"[t] -> {1, 2}, "PS"[c] -> 2, "PS"[e] -> 2}, 2]},
        {qi["Layers"], qi["Dagger"]["Layers"]}
    ],
    {{1/2, 3/2, 2, 5/2, 7/2}, {1/2, 3/2, 2, 5/2, 7/2}},
    TestID -> "QI-consecutive-phases-layout"
]

(* the coupler of a mesh cell shows its beam splitter's own (theta, phi); the cell phase is a separate box *)
VerificationTest[
    SeedRandom[19];
    With[{qi = QuantumInterferometer["Random"[4]]},
        With[{diagram = qi["Diagram", "ShowParameters" -> True]},
            {
                Count[diagram, _Rectangle, Infinity] == qi["PhaseShifterCount"],
                Cases[diagram, Tooltip[_, Row[{"(\[Theta], \[Phi]) = ", tf_, __}]] :> tf, Infinity] ==
                    Cases[qi["Elements"], ("BS"[t_, f_] -> _) :> {t, f}]
            }
        ]
    ],
    {True, True},
    TestID -> "QI-diagram-labels-beam-splitter-parameters"
]

VerificationTest[
    SeedRandom[18];
    With[{qi1 = QuantumInterferometer["Random"[4], Method -> "Reck"], qi2 = QuantumInterferometer["Random"[4]]},
        qi2[qi1]["Depth"] == qi1["Depth"] + qi2["Depth"]
    ],
    True,
    TestID -> "QI-composition-is-compact"
]

EndTestSection[]
