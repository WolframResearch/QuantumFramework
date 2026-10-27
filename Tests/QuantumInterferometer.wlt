(* Tests for QuantumInterferometer. Results are checked against independent references: the
   Fock-space circuit of BeamSplitterOperator, PhaseShiftOperator and SqueezeOperator, permanents,
   hafnians and closed forms. *)

Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"]

$exactMethods = {"Reck", "Clements"};
$methods = Append[$exactMethods, "CosineSine"];

rebuild[qi_] := QuantumInterferometer[qi["Elements"], qi["Modes"]]["Unitary"]

fockIndex[occupation_, d_] := FromDigits[occupation, d] + 1

(* W[[j, k]] = <1_j| U |1_k>, read off a Fock-space operator one photon at a time *)
transferMatrix[op_, m_, d_] := Transpose @ Map[
    Normal[op[FockState[#, d]]["StateVector"][[fockIndex[#, d] & /@ IdentityMatrix[m]]]] &,
    IdentityMatrix[m]
]

fockAmplitudes[circuit_, in_, d_, outs_] := Normal[circuit[FockState[in, d]]["StateVector"]][[FromDigits[#, d] + 1 & /@ outs]]

patterns[m_, nT_] := Select[Tuples[Range[0, nT], m], Total[#] == nT &];

(* substitute a concrete pattern into the formula for {n1, ..} -> {t1, ..} and Activate the sums *)
symbolicValue[f_, n_, t_] := Activate[f /. Thread[Array[Symbol["n" <> ToString[#]] &, Length[n]] -> n] /. Thread[Array[Symbol["t" <> ToString[#]] &, Length[t]] -> t]];

permanentAmplitudeReference[w_, s_, t_] := Permanent[w[[
    Catenate @ MapIndexed[ConstantArray[First[#2], #1] &, t], Catenate @ MapIndexed[ConstantArray[First[#2], #1] &, s]
]]] / Sqrt[Times @@ (s!) Times @@ (t!)]


BeginTestSection["QuantumInterferometer: meshes and construction"]

VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[\[Theta], \[Phi]] -> {1, 2}, "PS"[\[Alpha]] -> 2}, 2]},
        FullSimplify[qi["Unitary"] == transferMatrix[qi["CircuitOperator", 2], 2, 2], (\[Theta] | \[Phi] | \[Alpha]) \[Element] Reals]
    ],
    True,
    TestID -> "QI-element-blocks-match-Fock-operators"
]

VerificationTest[
    SeedRandom[11];
    {
        Table[With[{w = RandomVariate[CircularUnitaryMatrixDistribution[m]]},
            Max[Abs[rebuild[QuantumInterferometer[w, Method -> #]] - w]] < 10^-10 & /@ $methods], {m, 2, 8}],
        Table[With[{w = RandomVariate[CircularUnitaryMatrixDistribution[m]]},
            With[{u = rebuild[QuantumInterferometer[w, Method -> method]]},
                Max[Abs[u - w]] < 10^-10 && Max[Abs[u . ConjugateTranspose[u] - IdentityMatrix[m]]] < 10^-10]], {m, {16, 24, 32}}, {method, $exactMethods}]
    },
    {ConstantArray[True, {7, 3}], ConstantArray[True, {3, 2}]},
    TestID -> "QI-meshes-rebuild-Haar-unitaries"
]

VerificationTest[
    SeedRandom[404];
    With[{w = RandomVariate[CircularUnitaryMatrixDistribution[8]]},
        With[{before = Hash[w, "SHA256"]},
            And @@ Table[
                With[{a = QuantumInterferometer[w, Method -> method], b = QuantumInterferometer[w, Method -> method]},
                    a["Elements"] === b["Elements"] && a["Layers"] === b["Layers"]],
                {method, $methods}
            ] && Hash[w, "SHA256"] === before
        ]
    ],
    True,
    TestID -> "QI-decomposition-is-repeatable-and-leaves-its-input-untouched"
]

VerificationTest[
    TimeConstrained[
        Table[
            With[{qi = QuantumInterferometer[FourierMatrix[m], Method -> method]},
                Precision[qi["Elements"]] === Infinity &&
                    If[m <= 4,
                        AllTrue[Flatten[rebuild[qi] - FourierMatrix[m]], PossibleZeroQ],
                        Max[Abs[QuantumInterferometer[qi["Elements"] /. (h_[ps__] -> md_) :> (h[Sequence @@ N[{ps}]] -> md), m]["Unitary"] - N[FourierMatrix[m]]]] < 10^-10
                    ]
            ],
            {m, 3, 6}, {method, $exactMethods}
        ],
        120, $TimedOut
    ],
    ConstantArray[True, {4, 2}],
    TestID -> "QI-exact-Fourier-decompositions-stay-exact"
]

VerificationTest[
    {
        Table[With[{qi = QuantumInterferometer[w, Method -> method]}, AllTrue[Flatten[rebuild[qi] - w], PossibleZeroQ]],
            {w, {IdentityMatrix[4], PermutationMatrix[{3, 1, 4, 2}], DiagonalMatrix[Exp[I {1, 2, 3, 4}]], FourierMatrix[2]}}, {method, $exactMethods}],
        Table[Max[Abs[rebuild[QuantumInterferometer[w, Method -> "CosineSine"]] - w]] < 10^-10,
            {w, {IdentityMatrix[5], PermutationMatrix[{3, 1, 4, 2}], PermutationMatrix[{3, 4, 1, 2}], PermutationMatrix[{4, 5, 1, 2, 3}],
                DiagonalMatrix[Exp[I {1, 2, 3, 4}]], FourierMatrix[2], FourierMatrix[7]}}],
        QuantumInterferometer[IdentityMatrix[5]]["Elements"]
    },
    {ConstantArray[True, {4, 2}], ConstantArray[True, 7], {}},
    TestID -> "QI-degenerate-unitaries"
]

VerificationTest[
    With[{qi = TimeConstrained[QuantumInterferometer["Fourier"[8]], 10]},
        {Precision[qi["Unitary"]], Max[Abs[rebuild[qi] - FourierMatrix[8]]] < 10^-10}
    ],
    {MachinePrecision, True},
    TestID -> "QI-named-Fourier-is-numeric"
]

VerificationTest[
    FullSimplify[{
            QuantumInterferometer[{"PS"[\[Alpha]] -> 1, "PS"[\[Beta]] -> 2, "BS"[\[Theta], \[Gamma]] -> {1, 2}}, 2]["Unitary"] ==
                QuantumInterferometer[{"BS"[\[Theta], \[Gamma] + \[Alpha] - \[Beta]] -> {1, 2}, "PS"[\[Alpha]] -> 1, "PS"[\[Beta]] -> 2}, 2]["Unitary"],
            QuantumInterferometer[{"BS"[- \[Theta], \[Gamma]] -> {1, 2}}, 2]["Unitary"] ==
                QuantumInterferometer[{"BS"[\[Theta], \[Gamma] + Pi] -> {1, 2}}, 2]["Unitary"]
        }, Element[{\[Theta], \[Gamma], \[Alpha], \[Beta]}, Reals]],
    {True, True},
    TestID -> "QI-phase-absorption-identities"
]

(* every mesh element is a beam splitter with theta in [0, Pi/2], followed by one output phase column *)
VerificationTest[
    SeedRandom[20];
    With[{w = RandomVariate[CircularUnitaryMatrixDistribution[6]]},
        {
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
            ],
            #["Depth"] & /@ (QuantumInterferometer[w, Method -> #] & /@ $exactMethods)
        }
    ],
    {ConstantArray[True, {3, 4}], {9, 6}},
    TestID -> "QI-mesh-structure-and-depth"
]

VerificationTest[
    {
        QuantumInterferometer[{"BS"[0.5, 0.] -> {1, 2}, "BS"[0.5, 0.] -> {2, 3}, "BS"[0.5, 0.] -> {1, 2}}, 3]["Layers"],
        With[{qi = QuantumInterferometer[{"PS"[a] -> 1, "PS"[b] -> 1, "BS"[t] -> {1, 2}, "PS"[c] -> 2, "PS"[e] -> 2}, 2]},
            {qi["Layers"], qi["Dagger"]["Layers"]}]
    },
    {{1, 2, 3}, {{1/2, 3/2, 2, 5/2, 7/2}, {1/2, 3/2, 2, 5/2, 7/2}}},
    TestID -> "QI-layout-columns"
]

VerificationTest[
    SeedRandom[505];
    With[{w = RandomVariate[CircularUnitaryMatrixDistribution[6]]},
        With[{u = QuantumInterferometer[QuantumInterferometer[w]["Elements"], 6]["Unitary"]},
            Developer`PackedArrayQ[u] && Max[Abs[u - w]] < 10^-10
        ]
    ],
    True,
    TestID -> "QI-elements-unitary-packed-and-correct"
]

VerificationTest[
    SeedRandom[16];
    With[{qi1 = QuantumInterferometer["Random"[4]], qi2 = QuantumInterferometer["Random"[4], Method -> "Reck"]},
        Max[Abs[{
            qi2[qi1]["Unitary"] - qi2["Unitary"] . qi1["Unitary"],
            rebuild[qi2[qi1]] - qi2["Unitary"] . qi1["Unitary"],
            rebuild[qi1["Dagger"]] - ConjugateTranspose[qi1["Unitary"]]
        }]] < 10^-10 &&
        qi2[qi1]["Depth"] == qi1["Depth"] + qi2["Depth"] &&
        And @@ Map[
            With[{a = QuantumInterferometer[#[[1]], 4], b = QuantumInterferometer[#[[2]], 4]},
                b[a]["Layers"] === #[[3]] && Max[Abs[b[a]["Unitary"] - b["Unitary"] . a["Unitary"]]] < 10^-12] &,
            {
                {{"BS"[.4] -> {1, 2}}, {"BS"[.3] -> {3, 4}}, {1, 1}},
                {{"BS"[.4] -> {1, 2}}, {"BS"[.3] -> {2, 3}}, {1, 2}},
                {{"BS"[.4] -> {1, 4}}, {"BS"[.3] -> {2, 3}}, {1, 2}},
                {{"PS"[.2] -> 1}, {"PS"[.3] -> 1}, {1/2, 3/2}}
            }
        ]
    ],
    True,
    TestID -> "QI-composition-and-dagger"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[Pi / 4, 0] -> {1, 2}, "PS"[Pi / 3] -> 2}, 2]},
        {N[qi] === qi, N[qi]["Modes"], N[qi]["Elements"]}
    ],
    {True, 2, {"BS"[Pi / 4, 0] -> {1, 2}, "PS"[Pi / 3] -> {2}}},
    TestID -> "QI-N-keeps-the-object-valid"
]

(* exact tiny angles are elements, exact identities are dropped *)
VerificationTest[
    With[{eps = 1/10^12},
        With[{w = {{Sqrt[1 - eps^2], -eps}, {eps, Sqrt[1 - eps^2]}}, sq = QuantumInterferometer[{"S"[eps] -> 1}, 1],
              identities = QuantumInterferometer[{"BS"[0, Pi/7] -> {1, 2}, "PS"[2 Pi] -> 1, "S"[0] -> 2}, 2]},
            QuantumInterferometer[{"BS"[eps, 0] -> {1, 2}}, 2]["Amplitude", {1, 0} -> {0, 1}] === Sin[eps] &&
            QuantumInterferometer[{"PS"[eps] -> 1}, 1]["Amplitude", {1} -> {1}] === Exp[I eps] &&
            sq["SqueezerCount"] === 1 && Precision[sq["Amplitude", {0} -> {2}]] === Infinity &&
            QuantumInterferometer[w]["Unitary"] === w && QuantumInterferometer[w]["Amplitude", {1, 0} -> {0, 1}] === eps &&
            AllTrue[Flatten[FullSimplify[QuantumInterferometer[QuantumInterferometer[w]["Elements"], 2]["Unitary"] - w]], PossibleZeroQ] &&
            identities["Elements"] === {} && identities["Unitary"] === IdentityMatrix[2] && identities["Amplitude", {1, 2} -> {1, 2}] === 1
        ]
    ],
    True,
    TestID -> "QI-exact-small-elements-kept-exact-identities-dropped"
]

EndTestSection[]


BeginTestSection["QuantumInterferometer: elements without modes"]

VerificationTest[
    With[{xi = (2 + I)/7, explicit = QuantumInterferometer[{"BS"[theta, 0] -> {1, 2}}, 2],
          phased = QuantumInterferometer[{"BS"[theta, phi] -> {1, 2}}, 2]},
        QuantumInterferometer[{"BS"[theta]}] === explicit && QuantumInterferometer[{"BeamSplitter"[theta]}] === explicit &&
        (QuantumInterferometer /@ {{"BS"[{theta, phi}]}, {"BeamSplitter"[theta, phi]}, {"BeamSplitter"[{theta, phi}]}}) === ConstantArray[phased, 3] &&
        (#["Modes"] & /@ (QuantumInterferometer /@ {{"PS"[alpha]}, {"PhaseShift"[alpha]}, {"PS"[100]}})) === {1, 1, 1} &&
        QuantumInterferometer[{"S"[xi]}] === QuantumInterferometer[{"Squeeze"[xi]}] === QuantumInterferometer[{"S"[xi] -> 1}, 1] &&
        QuantumInterferometer[{"BS"[Pi/6], "PS"[Pi/7] -> 4, "PhaseShift"[Pi/11], "BeamSplitter"[Pi/5, Pi/13] -> {4, 2}}] ===
            QuantumInterferometer[{"BS"[Pi/6, 0] -> {1, 2}, "PS"[Pi/7] -> 4, "PS"[Pi/11] -> 1, "BS"[Pi/5, Pi/13] -> {4, 2}}, 4]
    ],
    True,
    TestID -> "QI-bare-elements-equal-their-explicit-forms"
]

VerificationTest[
    With[{repeated = QuantumInterferometer[{"BS"[Pi/6], "BeamSplitter"[Pi/5]}], wide = QuantumInterferometer[{"BS"[Pi/6], "PS"[Pi/7]}, 5]},
        {
            repeated["Modes"], repeated["Layers"],
            FullSimplify[repeated["Unitary"] - {{Cos[11 Pi/30], -Sin[11 Pi/30]}, {Sin[11 Pi/30], Cos[11 Pi/30]}}],
            wide["Modes"], wide["Unitary"][[3 ;;, 3 ;;]] === IdentityMatrix[3],
            #["Modes"] & /@ {QuantumInterferometer[{"BS"[0]}], QuantumInterferometer[{"S"[0]}], QuantumInterferometer[{"BS"[0], "PS"[0] -> 7}]},
            QuantumInterferometer[{{1}}]["Unitary"], QuantumInterferometer[{}, 4]["Unitary"] === IdentityMatrix[4]
        }
    ],
    {2, {1, 2}, ConstantArray[0, {2, 2}], 5, True, {2, 1, 7}, {{1}}, True},
    TestID -> "QI-bare-elements-mode-inference"
]

VerificationTest[
    {
        QuantumInterferometer[{"BS"[theta]}, 1],
        QuantumInterferometer[{"BS"[theta], "PS"[alpha] -> 4}, 3],
        QuantumInterferometer[{"BeamSplitter"[Pi/6, I]}],
        QuantumInterferometer[{"PhaseShift"[I/3]}]
    },
    ConstantArray[$Failed, 4],
    {QuantumInterferometer::modes, QuantumInterferometer::modes, QuantumInterferometer::angles, QuantumInterferometer::angles},
    TestID -> "QI-bare-elements-keep-mode-and-angle-checks"
]

VerificationTest[
    {
        QuantumInterferometer[{"BS"[theta] -> {1, 1}}],
        QuantumInterferometer[{"PS"[alpha] -> 0}],
        QuantumInterferometer[{"BS"[]}],
        QuantumInterferometer[{"UnknownElement"[theta]}]
    },
    ConstantArray[$Failed, 4],
    {QuantumInterferometer::elem, QuantumInterferometer::elem, QuantumInterferometer::elem, General::stop},
    TestID -> "QI-bare-elements-reject-malformed-elements"
]

(* the message names every element kind the constructor accepts *)
VerificationTest[
    StringContainsQ[QuantumInterferometer::elem, #] & /@ {"\"BS\"[", "\"PS\"[", "\"S\"[", "\"Loss\"["},
    {True, True, True, True},
    TestID -> "QI-element-message-lists-every-element-kind"
]

VerificationTest[
    And @@ Map[Function[p,
        With[{t = N[3/10, p], f = N[7/11, p]},
            With[{a = QuantumInterferometer[{"BS"[t], "PS"[f]}], b = QuantumInterferometer[{"BS"[t, 0] -> {1, 2}, "PS"[f] -> 1}, 2]},
                a === b && Precision[a["Amplitude", {2, 1} -> {1, 2}]] === Precision[b["Amplitude", {2, 1} -> {1, 2}]]
            ]
        ]], {MachinePrecision, 80}],
    True,
    TestID -> "QI-bare-elements-preserve-precision"
]

EndTestSection[]


BeginTestSection["QuantumInterferometer: validation"]

VerificationTest[
    {
        QuantumInterferometer[{{1, 1}, {1, 1}}],
        QuantumInterferometer[{{a, b}, {c, d}}],
        QuantumInterferometer[IdentityMatrix[2], Method -> "Unknown"],
        QuantumInterferometer[{"PS"[.3] -> 3}, 2],
        QuantumInterferometer[{"PS"[I] -> 1}, 1],
        QuantumInterferometer[{"BS"[I/5, 0] -> {1, 2}}, 2]
    },
    ConstantArray[$Failed, 6],
    {QuantumInterferometer::nonunitary, QuantumInterferometer::symbolic, QuantumInterferometer::method, QuantumInterferometer::modes,
        QuantumInterferometer::angles, QuantumInterferometer::angles},
    TestID -> "QI-constructor-errors"
]

VerificationTest[
    {
        QuantumInterferometer["Fourier"[3]]["Probabilities", {1, 1}],
        QuantumInterferometer["Fourier"[3]]["CircuitOperator"],
        QuantumInterferometer[{"BS"[\[Theta], \[Phi]] -> {1, 2}}, 2]["MatrixPlot"],
        QuantumInterferometer[{}, 2][QuantumInterferometer[{}, 3]],
        QuantumInterferometer[{}, 2][FockState[{1}, 3]]
    },
    ConstantArray[$Failed, 5],
    {QuantumInterferometer::occ, QuantumInterferometer::levels, QuantumInterferometer::notnumeric, QuantumInterferometer::compose, QuantumInterferometer::state},
    TestID -> "QI-property-errors"
]

(* with 2 levels the bunched outputs |2,0> and |0,2> do not fit and are dropped *)
VerificationTest[
    QuantumInterferometer[{"BS"[Pi / 4, 0] -> {1, 2}}, 2][FockState[{1, 1}, 2]]["Norm"],
    0,
    {QuantumInterferometer::truncation},
    TestID -> "QI-warns-on-short-truncation"
]

VerificationTest[
    {
        QuantumInterferometer[{"S"[.2] -> 1}, 1]["Probabilities", {0}, "Cutoff" -> "bad"],
        QuantumInterferometer[{}, 1]["Probabilities", "Coherent"[{1/3}], "Cutoff" -> "bad"],
        QuantumInterferometer[{"BS"[Pi/6, 0] -> {1, 2}}, 2]["Probabilities", {1, 0}, "Cutoff" -> -1],
        QuantumInterferometer[{"S"[.2] -> 1}, 1][{0}, "Cutoff" -> 1/2]
    },
    ConstantArray[$Failed, 4],
    {QuantumInterferometer::cutoff, QuantumInterferometer::cutoff, QuantumInterferometer::cutoff, General::stop},
    TestID -> "QI-cutoff-must-be-a-nonnegative-integer"
]

VerificationTest[
    With[{sq = QuantumInterferometer[{"S"[1/5] -> 1}, 1], passive = QuantumInterferometer[{}, 1]},
        {
            FullSimplify[Values[sq["Probabilities", {0}, "Cutoff" -> 0]]],
            passive["Probabilities", "Coherent"[{(1 + I)/3}], "Cutoff" -> 0],
            passive["Probabilities", {1}, "Cutoff" -> 0] === passive["Probabilities", {1}],
            sq[{1}, "Cutoff" -> 0]
        }
    ],
    {{Sech[1/5]}, <|QuditName[{0}] -> Exp[-2/9]|>, True, <||>},
    TestID -> "QI-zero-and-below-input-cutoffs"
]

(* a bad photon number next to a coherent amplitude is reported once per query, naming the query's input *)
VerificationTest[
    With[{qi = QuantumInterferometer["Fourier"[3]]},
        (Lookup[EvaluationData[#], {"Result", "MessagesExpressions"}] /. Hold[Message[name_, args___]] :> {Hold[name], {args}}) & /@ {
            Unevaluated[qi["Amplitude", {-1, "Coherent"[1/2], 0} -> {0, 1, 0}]],
            Unevaluated[qi["Probabilities", {-1, "Coherent"[1/2], 0}, "Cutoff" -> 2]]
        }
    ],
    ConstantArray[{$Failed, {{Hold[QuantumInterferometer::occ], {{-1, "Coherent"[1/2], 0}, 3}}}}, 2],
    {QuantumInterferometer::occ, QuantumInterferometer::occ},
    TestID -> "QI-coherent-bad-photon-numbers-warn-once"
]

(* an amplitude list in a "Coherent" slot is rejected, and lossy queries report the same way *)
VerificationTest[
    With[{qi = QuantumInterferometer["Fourier"[3]], lossy = QuantumInterferometer[{"BS"[Pi/6], "Loss"[2/5] -> 2}]},
        (Lookup[EvaluationData[#], {"Result", "MessagesExpressions"}] /. Hold[Message[name_, args___]] :> {Hold[name], {args}}) & /@ {
            Unevaluated[qi["Amplitude", {2, "Coherent"[{1, 2}], 0} -> {0, 1, 0}]],
            Unevaluated[lossy["Probability", {-1, "Coherent"[1/2]} -> {0, 1}]]
        }
    ],
    {
        {$Failed, {{Hold[QuantumInterferometer::occ], {{2, "Coherent"[{1, 2}], 0}, 3}}}},
        {$Failed, {{Hold[QuantumInterferometer::occ], {{-1, "Coherent"[1/2]}, 2}}}}
    },
    {QuantumInterferometer::occ, QuantumInterferometer::occ},
    TestID -> "QI-coherent-nested-and-lossy-bad-inputs-warn-once"
]

(* a "Coherent" slot with the wrong number of arguments is rejected the same way *)
VerificationTest[
    With[{qi = QuantumInterferometer["Fourier"[3]]},
        (Lookup[EvaluationData[#], {"Result", "MessagesExpressions"}] /. Hold[Message[name_, args___]] :> {Hold[name], {args}}) & /@ {
            Unevaluated[qi["Amplitude", {0, "Coherent"[1, 2], 0} -> {0, 1, 0}]],
            Unevaluated[qi["Probability", {0, "Coherent"[], 0} -> {0, 1, 0}]]
        }
    ],
    {
        {$Failed, {{Hold[QuantumInterferometer::occ], {{0, "Coherent"[1, 2], 0}, 3}}}},
        {$Failed, {{Hold[QuantumInterferometer::occ], {{0, "Coherent"[], 0}, 3}}}}
    },
    {QuantumInterferometer::occ, QuantumInterferometer::occ},
    TestID -> "QI-coherent-malformed-slots-warn-once"
]

(* a malformed "Coherent" output state is rejected the same way *)
VerificationTest[
    With[{qi = QuantumInterferometer["Fourier"[3]]},
        (Lookup[EvaluationData[#], {"Result", "MessagesExpressions"}] /. Hold[Message[name_, args___]] :> {Hold[name], {args}}) & /@ {
            Unevaluated[qi["Amplitude", {0, "Coherent"[1/2], 0} -> "Coherent"[{1, {2}, 0}]]],
            Unevaluated[qi["Probability", {0, "Coherent"[1/2], 0} -> "Coherent"[{1, {2}, 0}]]]
        }
    ],
    ConstantArray[{$Failed, {{Hold[QuantumInterferometer::occ], {"Coherent"[{1, {2}, 0}], 3}}}}, 2],
    {QuantumInterferometer::occ, QuantumInterferometer::occ},
    TestID -> "QI-coherent-malformed-output-warns-once"
]

EndTestSection[]


BeginTestSection["QuantumInterferometer: photons through passive networks"]

VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[Pi / 4, 0] -> {1, 2}}, 2]},
        {
            qi["Amplitude", {1, 1} -> {1, 1}],
            qi["Probability", Ket[{1, 1}] -> QuditName[{2, 0}]],
            qi["Probabilities", {1, 1}],
            qi[Ket[{1, 1}]]
        }
    ],
    {
        0,
        1/2,
        <|QuditName[{2, 0}] -> 1/2, QuditName[{1, 1}] -> 0, QuditName[{0, 2}] -> 1/2|>,
        <|QuditName[{2, 0}] -> - 1/Sqrt[2], QuditName[{0, 2}] -> 1/Sqrt[2]|>
    },
    TestID -> "QI-Hong-Ou-Mandel"
]

(* one photon per port of the Fourier tritter: the six patterns whose mode labels do not sum to a multiple of 3 are dark *)
VerificationTest[
    Keys @ Select[QuantumInterferometer[FourierMatrix[3]]["Probabilities", {1, 1, 1}], PossibleZeroQ],
    QuditName /@ {{2, 1, 0}, {2, 0, 1}, {1, 2, 0}, {1, 0, 2}, {0, 2, 1}, {0, 1, 2}},
    TestID -> "QI-Fourier-tritter-suppression-law"
]

VerificationTest[
    SeedRandom[15];
    With[{qi = QuantumInterferometer["Random"[5]], qi3 = QuantumInterferometer["Random"[3]]},
        {
            Chop[Total[qi["Probabilities", {1, 1, 0, 1, 0}]] - 1],
            Chop[qi["Probability", {1, 1, 0, 1, 0} -> {0, 0, 2, 0, 1}] - qi["Probabilities", {1, 1, 0, 1, 0}][QuditName[{0, 0, 2, 0, 1}]]],
            qi3["Probabilities", {1, 1, 0}, "Cutoff" -> 1] === qi3["Probabilities", {1, 1, 0}]
        }
    ],
    {0, 0, True},
    TestID -> "QI-probabilities-normalized-and-cutoff-free"
]

VerificationTest[
    SeedRandom[13];
    With[{qi3 = QuantumInterferometer["Random"[3]], qi4 = QuantumInterferometer["Random"[4]],
          superposition = (FockState[{2, 0, 0}, 3] + I FockState[{0, 1, 1}, 3]) / Sqrt[2]},
        {
            With[{amplitudes = qi3[{1, 1, 1}], out = qi3["CircuitOperator", 4][FockState[{1, 1, 1}, 4]]["StateVector"]},
                Max[Abs[Lookup[amplitudes, QuditName[#], 0] - out[[fockIndex[#, 4]]]] & /@ Tuples[Range[0, 3], 3]] < 10^-10],
            Max[Abs[Normal[qi4[FockState[{2, 0, 1, 0}, 4]]["StateVector"]] - Normal[qi4["CircuitOperator", 4][FockState[{2, 0, 1, 0}, 4]]["StateVector"]]]] < 10^-10,
            Max[Abs[Normal[qi3[superposition]["StateVector"] - qi3["CircuitOperator", 3][superposition]["StateVector"]]]] < 10^-10
        }
    ],
    {True, True, True},
    TestID -> "QI-states-match-Fock-circuit"
]

VerificationTest[
    With[{d = 4},
        With[{v = (FockState[{2, 0}, d] + I FockState[{0, 2}, d])/Sqrt[2], u = FockState[{1, 0}, d]},
            With[{rho = QuantumState[.7 Normal[v["DensityMatrix"]] + .3 Normal[u["DensityMatrix"]], {d, d}],
                  qi = QuantumInterferometer[{"BS"[.41, .23] -> {1, 2}, "PS"[.37] -> 2}, 2]},
                With[{matrix = Normal[qi[rho]["DensityMatrix"]], expected = PhaseShiftOperator[.37, d, {2}][
                        BeamSplitterOperator[{.41, .23}, d, {1, 2}, Method -> "Recurrence"][rho]]},
                    Max[Abs[matrix - Normal[expected["DensityMatrix"]]]] < 10^-12 &&
                        Max[Abs[matrix - ConjugateTranspose[matrix]]] < 10^-12 && Abs[Tr[matrix] - 1] < 10^-12
                ]
            ]
        ]
    ],
    True,
    TestID -> "QI-mixed-state-matches-Fock-operators"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[Pi/6, Pi/5] -> {1, 2}}, 2],
          in = QuantumState[SparseArray[{2 -> Exp[I phi] Sin[theta], 4 -> Cos[theta]}, 9], {3, 3}, "ParameterSpec" -> {{theta, -Pi, Pi}, {phi, -Pi, Pi}}]},
        With[{out = qi[in], values = <|theta -> Pi/7, phi -> Pi/11|>},
            out["ParameterSpec"] === in["ParameterSpec"] && out[<|theta -> Pi/7|>]["Parameters"] === {phi} &&
            AllTrue[Normal[FullSimplify[out[values]["StateVector"] - qi[in[values]]["StateVector"]]], PossibleZeroQ]
        ]
    ],
    True,
    TestID -> "QI-parametric-input-state-keeps-its-parameters"
]

(* bunched amplitudes against the permanent: every pattern with 1 to 4 photons in 4 modes, machine and 35 digits *)
VerificationTest[
    With[{w = DiagonalMatrix[{1, I, -1, -I}] . FourierMatrix[4] . DiagonalMatrix[{1, -I, -1, I}]},
        With[{machine = QuantumInterferometer[N[w]], digits = QuantumInterferometer[N[w, 35]]},
            And @@ Flatten @ Table[
                With[{ref = permanentAmplitudeReference[w, s, t]},
                    Abs[machine["Amplitude", s -> t] - N[ref]] < 10^-12 && Abs[digits["Amplitude", s -> t] - N[ref, 35]] < 10^-28],
                {n, 1, 4}, {s, FrobeniusSolve[ConstantArray[1, 4], n]}, {t, FrobeniusSolve[ConstantArray[1, 4], n]}
            ]
        ]
    ],
    True,
    TestID -> "QI-bunched-amplitudes-match-the-permanent"
]

VerificationTest[
    SeedRandom[4871];
    With[{w = RandomVariate[CircularUnitaryMatrixDistribution[16]], f = N[FourierMatrix[16]],
          s = ReplacePart[ConstantArray[0, 16], {2 -> 9, 11 -> 7}], t = ReplacePart[ConstantArray[0, 16], {3 -> 6, 8 -> 4, 14 -> 6}]},
        With[{ref = permanentAmplitudeReference[w, s, t], one = PadRight[ConstantArray[1, 12], 16], twelve = PadLeft[{12}, 16]},
            {
                Abs[QuantumInterferometer[w]["Amplitude", s -> t] - ref] < 10^-10 Abs[ref],
                Abs[QuantumInterferometer[ConjugateTranspose[w]]["Amplitude", t -> s] - Conjugate[ref]] < 10^-10 Abs[ref],
                Abs[QuantumInterferometer[f]["Amplitude", one -> twelve] - Sqrt[12!] Times @@ f[[-1, Range[12]]]] < 10^-12,
                Abs[QuantumInterferometer[N[DiagonalMatrix[{Exp[37 I / 100], Exp[- I / 5]}]]]["Amplitude", {24, 0} -> {24, 0}] - Exp[37 I / 100] ^ 24] < 10^-13
            }
        ]
    ],
    {True, True, True, True},
    TestID -> "QI-large-bunched-amplitudes-keep-accuracy"
]

VerificationTest[
    With[{u = N[FourierMatrix[16][[{5, 9}, {2, 13}]], 80],
          s = ReplacePart[ConstantArray[0, 16], {2 -> 8, 13 -> 8}], t = ReplacePart[ConstantArray[0, 16], {5 -> 10, 9 -> 6}]},
        With[{ref = Sqrt[10! 6! / (8! 8!)] Sum[Binomial[8, k] Binomial[8, 10 - k] u[[1, 1]] ^ k u[[2, 1]] ^ (8 - k) u[[1, 2]] ^ (10 - k) u[[2, 2]] ^ (k - 2), {k, 2, 8}]},
            Abs[QuantumInterferometer[N[FourierMatrix[16], 40]]["Amplitude", s -> t] - ref] < 10^-28 Abs[ref]
        ]
    ],
    True,
    TestID -> "QI-bunched-amplitudes-at-arbitrary-precision"
]

VerificationTest[
    {
        QuantumInterferometer[FourierMatrix[3]]["Amplitude", {4, 0, 0} -> {1, 2, 1}],
        QuantumInterferometer[{{1, -1}, {1, 1}} / Sqrt[2]]["Amplitude", {1, 1} -> {2, 0}],
        Precision[QuantumInterferometer[ReplacePart[IdentityMatrix[4], {4, 4} -> N[Exp[I/7]]]]["Amplitude", {1, 0, 0, 0} -> {1, 0, 0, 0}]]
    },
    {Sqrt[12] / 9, -1 / Sqrt[2], Infinity},
    TestID -> "QI-exact-entries-give-exact-amplitudes"
]

VerificationTest[
    SeedRandom[4873];
    With[{qi = QuantumInterferometer[RandomVariate[CircularUnitaryMatrixDistribution[3]]], s = {6, 6, 0}, outs = FrobeniusSolve[{1, 1, 1}, 12]},
        With[{amps = qi["Amplitude", s -> #] & /@ outs},
            {Abs[Total[Abs[amps] ^ 2] - 1] < 10^-11, Max[Abs[(qi["Dagger"]["Amplitude", # -> s] & /@ outs) - Conjugate[amps]]] < 10^-12}
        ]
    ],
    {True, True},
    TestID -> "QI-bunched-sector-normalization-and-adjoint"
]

VerificationTest[
    With[{exact = FourierMatrix[4]},
        With[{w = ReplacePart[N[exact, 50], {1, 1} -> N[exact[[1, 1]]]]},
            Abs[QuantumInterferometer[w]["Amplitude", {2, 1, 1, 0} -> {0, 1, 1, 2}] - N[permanentAmplitudeReference[exact, {2, 1, 1, 0}, {0, 1, 1, 2}], 50]] < 10^-12
        ]
    ],
    True,
    TestID -> "QI-mixed-precision-unitary"
]

EndTestSection[]


BeginTestSection["QuantumInterferometer: sampling"]

(* total variation distance of 5000 samples from "Probabilities", lossless and lossy; sampling noise alone is about 0.025 here *)
VerificationTest[
    SeedRandom[31];
    With[{
        qi = QuantumInterferometer["Random"[4]],
        lossy = QuantumInterferometer[{"BS"[0.7, 0.3] -> {1, 2}, "Loss"[0.6] -> 1, "BS"[0.4, 1.1] -> {2, 3}, "Loss"[0.85] -> 3, "BS"[1.0, -0.5] -> {1, 2}}, 3],
        k = 5000
    },
        Map[
            With[{probs = #[[1]]["Probabilities", #[[2]]], counts = #[[1]]["SampleCounts", #[[2]], k]},
                With[{keys = Union[Keys[probs], Keys[counts]]},
                    Total[counts] === k && Complement[Keys[counts], Keys[probs]] === {} &&
                        Total[Abs[Lookup[probs, keys, 0] - Lookup[counts, keys, 0] / k]] / 2 < 0.06
                ]
            ] &,
            {{qi, {1, 1, 1, 0}}, {qi, {2, 1, 0, 0}}, {lossy, {1, 1, 1}}, {lossy, {2, 1, 0}}}
        ]
    ],
    {True, True, True, True},
    TestID -> "QI-samples-match-probabilities"
]

VerificationTest[
    SeedRandom[32];
    With[{hom = QuantumInterferometer[{"BS"[Pi / 4, 0] -> {1, 2}}, 2], qi = QuantumInterferometer["Random"[3]]},
        {
            Keys[hom["SampleCounts", {1, 1}, 500]] === {QuditName[{2, 0}], QuditName[{0, 2}]},
            Union[Total[#["Name"]] & /@ qi["Samples", {2, 0, 1}, 50]],
            qi["Samples", {0, 0, 0}, 2],
            qi["Samples", {1, 0, 0}, 0],
            Total[qi["SampleCounts", {1, 1, 0}, 40]]
        }
    ],
    {True, {3}, {QuditName[{0, 0, 0}], QuditName[{0, 0, 0}]}, {}, 40},
    TestID -> "QI-samples-HOM-and-edge-cases"
]

VerificationTest[
    {
        QuantumInterferometer[{"S"[0.3] -> 1, "BS"[Pi / 4, 0] -> {1, 2}}, 2]["Samples", {1, 0}, 5],
        QuantumInterferometer[{"BS"[t, 0] -> {1, 2}}, 2]["Samples", {1, 0}, 5],
        QuantumInterferometer["Random"[2]]["Samples", {1, 0, 0}, 5]
    },
    {$Failed, $Failed, $Failed},
    {QuantumInterferometer::squeezed, QuantumInterferometer::notnumeric, QuantumInterferometer::occ},
    TestID -> "QI-samples-rejects-squeezers-symbolic-and-bad-patterns"
]

(* The exact law of the sampler: every order of the photons' columns and every sequence of output rows,
   weighted with the sampler's own step and binned to the retained modes. Several photons in one input
   mode, and loss through environment rows, must both reproduce "Probabilities" to round-off. *)
samplerLaw[qi_, s_] := With[{
    a = Wolfram`QuantumFramework`SecondQuantization`Interferometer`PackagePrivate`samplingColumns[qi, s],
    weights = Wolfram`QuantumFramework`SecondQuantization`Interferometer`PackagePrivate`photonWeights,
    m = qi["Modes"], n = Total[s]
},
    Merge[
        Map[
            Function[b, Merge[
                KeyValueMap[<|BinCounts[#1, {1, m + 1, 1}] -> #2|> &,
                    Nest[
                        Association @ Catenate @ KeyValueMap[
                            Function[{rows, p}, With[{w = weights[b, rows, Length[rows] + 1]}, MapIndexed[Append[rows, First[#2]] -> p #1 / Total[w] &, w]]],
                            #
                        ] &,
                        <|{} -> 1|>,
                        n
                    ]
                ],
                Total
            ]],
            a[[All, #]] & /@ Permutations[Range[n]]
        ] / n!,
        Total
    ]
]

VerificationTest[
    SeedRandom[33];
    With[{
        qi = QuantumInterferometer["Random"[3]],
        lossy = QuantumInterferometer[{"BS"[0.7, 0.3] -> {1, 2}, "Loss"[0.6] -> 1, "BS"[0.4, 1.1] -> {2, 3}, "Loss"[0.85] -> 3, "BS"[1.0, -0.5] -> {1, 2}}, 3]
    },
        Map[
            With[{exact = KeyMap[#["Name"] &, #[[1]]["Probabilities", #[[2]]]], law = samplerLaw @@ #},
                Max @ Abs[Lookup[exact, #, 0] - Lookup[law, #, 0]] &[Union[Keys[exact], Keys[law]]] < 10^-12
            ] &,
            {
                {QuantumInterferometer[{"BS"[Pi / 4, 0] -> {1, 2}}, 2], {1, 1}}, {QuantumInterferometer["Random"[2]], {2, 2}},
                {qi, {1, 1, 1}}, {qi, {2, 1, 0}}, {qi, {3, 0, 0}}, {lossy, {1, 1, 1}}, {lossy, {2, 1, 0}}
            }
        ]
    ],
    ConstantArray[True, 7],
    TestID -> "QI-sampler-law-is-exact-with-repeated-inputs-and-loss"
]

(* pure loss on one mode thins n photons binomially; a lossy mesh keeps at most the input photons *)
VerificationTest[
    SeedRandom[34];
    With[{
        attenuator = QuantumInterferometer[{"Loss"[1/4]}],
        lossy = QuantumInterferometer[{"BS"[0.7, 0.3] -> {1, 2}, "Loss"[0.6] -> 1, "BS"[0.4, 1.1] -> {2, 3}}, 3]
    },
        {
            With[{counts = attenuator["SampleCounts", {3}, 20000], outputs = QuditName /@ Range[0, 3]},
                Max @ Abs[Lookup[counts, outputs, 0] / 20000. - PDF[BinomialDistribution[3, 1/4], Range[0, 3]]] < 0.015
            ],
            SubsetQ[Range[0, 3], Union[Total[#["Name"]] & /@ lossy["Samples", {1, 1, 1}, 200]]],
            Total[lossy["SampleCounts", {2, 0, 1}, 50]]
        }
    ],
    {True, True, 50},
    TestID -> "QI-samples-lossy-binomial-thinning-and-photon-budget"
]

(* coherent light leaves a passive network, lossy or not, as |K a>: independent Poisson counts with means |(K a)_j|^2 *)
VerificationTest[
    SeedRandom[35];
    With[{
        lossy = QuantumInterferometer[{"BS"[0.7, 0.3] -> {1, 2}, "Loss"[0.6] -> 1, "BS"[0.4, 1.1] -> {2, 3}, "Loss"[0.85] -> 3}, 3],
        a = {0.8, 0.3 I, 0}, k = 20000
    },
        With[{counts = #["Name"] & /@ lossy["Samples", "Coherent"[a], k], means = Abs[lossy["TransferMatrix"] . a] ^ 2},
            {
                Max[Abs[Mean[counts] - means] / Sqrt[means / k]] < 5,
                Max[Abs[Variance[counts] - means] / Sqrt[(means + 2 means ^ 2) / k]] < 5,
                Max[Abs[Covariance[counts] - DiagonalMatrix[Variance[counts]]] / Sqrt[Outer[Times, means, means] / k]] < 5,
                Keys[lossy["SampleCounts", "Coherent"[{0, 0, 0}], 7]],
                Keys[lossy["SampleCounts", {"Coherent"[0], "Coherent"[0], "Coherent"[0]}, 7]],
                lossy["Samples", "Coherent"[a], 0]
            }
        ]
    ],
    {True, True, True, {QuditName[{0, 0, 0}]}, {QuditName[{0, 0, 0}]}, {}},
    TestID -> "QI-samples-coherent-poisson-counts"
]

VerificationTest[
    With[{qi = QuantumInterferometer["Random"[2]]},
        {
            qi["Samples", {1, "Coherent"[0.3]}, 3],
            qi["Samples", {1, 1}],
            qi["SampleCounts", {1, 1}, 3.],
            qi["Samples", "Coherent"[{x, 0}], 3],
            QuantumInterferometer[{"S"[0.3] -> 1, "BS"[Pi / 4, 0] -> {1, 2}}, 2]["Samples", "Coherent"[{1, 0}], 3]
        }
    ],
    ConstantArray[$Failed, 5],
    {
        QuantumInterferometer::samplecoherent, QuantumInterferometer::samplecount, QuantumInterferometer::samplecount,
        QuantumInterferometer::coherentnumeric, QuantumInterferometer::squeezed
    },
    TestID -> "QI-samples-messages-name-the-unsupported-input-or-count"
]

EndTestSection[]


BeginTestSection["QuantumInterferometer: symbolic photon numbers"]

VerificationTest[
    SeedRandom[41];
    With[{bs = QuantumInterferometer[{"BS"[Pi / 6, Pi / 3] -> {1, 2}}, 2], q3 = QuantumInterferometer["Random"[3]], q4 = QuantumInterferometer["Random"[4]]},
        With[{f2 = bs["Amplitude", {n1, n2} -> {t1, t2}], f3 = q3["Amplitude", {n1, n2, n3} -> {t1, t2, t3}], f4 = q4["Amplitude", {n1, n2, n3, n4} -> {t1, t2, t3, t4}]},
            {
                Union[FullSimplify[symbolicValue[f2, {2, 5}, #] - bs["Amplitude", {2, 5} -> #]] & /@ patterns[2, 7]],
                Max[Abs[N[symbolicValue[f3, {1, 2, 3}, #]] - q3["Amplitude", {1, 2, 3} -> #]] & /@ patterns[3, 6]] < 10^-12,
                Max[Abs[N[symbolicValue[f4, {1, 2, 0, 1}, #]] - q4["Amplitude", {1, 2, 0, 1} -> #]] & /@ patterns[4, 4]] < 10^-12,
                symbolicValue[f4, {1, 2, 0, 1}, {1, 1, 1, 2}],
                Length[Cases[f4, Inactive[Sum][_, iters__] :> {iters}, Infinity][[1]]]
            }
        ]
    ],
    {{0}, True, True, 0, 3},
    TestID -> "QI-symbolic-photon-numbers-match-permanents"
]

(* a mode no beam splitter touches keeps its photons *)
VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[Pi / 5, Pi / 7] -> {2, 3}, "PS"[Pi / 3] -> 1}, 3], single = QuantumInterferometer[{"PS"[Pi / 3] -> 1}, 1]},
        {
            Union[FullSimplify[symbolicValue[qi["Amplitude", {n1, n2, n3} -> {t1, t2, t3}], {1, 1, 0}, #] - qi["Amplitude", {1, 1, 0} -> #]] & /@ patterns[3, 2]],
            FullSimplify[single["Amplitude", {n1} -> {t1}] /. {n1 -> 2, t1 -> 2}] == single["Amplitude", {2} -> {2}]
        }
    ],
    {{0}, True},
    TestID -> "QI-symbolic-photon-numbers-untouched-modes"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[\[Theta], \[Phi]] -> {1, 2}}, 2], hom = QuantumInterferometer[{"BS"[Pi / 4, 0] -> {1, 2}}, 2]},
        {
            FullSimplify[qi["Amplitude", {nx, ny} -> {tx, nx + ny - tx}] == WignerD[{(nx + ny) / 2, (nx - ny) / 2, tx - (nx + ny) / 2}, \[Phi], 2 \[Theta], - \[Phi]],
                Element[{nx, ny, tx}, Integers] && nx >= 0 && ny >= 0 && 0 <= tx <= nx + ny],
            And @@ Flatten @ Table[FullSimplify[(hom["Amplitude", {nx, nx} -> {tx, 2 nx - tx}] /. {nx -> n, tx -> t}) - hom["Amplitude", {n, n} -> {t, 2 n - t}]] === 0, {n, 0, 5}, {t, 0, 2 n}],
            FullSimplify[With[{amplitudes = Table[qi["Amplitude", {nx, ny} -> {tx, nx + ny - tx}] /. {nx -> 2, ny -> 1, tx -> t}, {t, 0, 3}]},
                {Total[Abs[amplitudes] ^ 2], amplitudes /. \[Theta] -> 0, amplitudes /. \[Theta] -> Pi / 2}], Element[{\[Theta], \[Phi]}, Reals]]
        }
    ],
    {True, True, {1, {0, 0, 1, 0}, {0, - Exp[I \[Phi]], 0, 0}}},
    TestID -> "QI-symbolic-two-mode-WignerD-normalization-and-limits"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[Pi / 6, Pi / 7] -> {3, 1}, "BS"[Pi / 5, Pi / 11] -> {2, 3}, "BS"[Pi / 4, Pi / 13] -> {1, 2}}, 3]},
        With[{amplitude = qi["Amplitude", {nx, ny, nz} -> {tx, ty, tz}], vector = qi["CircuitOperator", 3][FockState[{1, 1, 0}, 3]]["StateVector"]},
            Max[Abs[N[Activate[amplitude /. Thread[{nx, ny, nz, tx, ty, tz} -> Join[{1, 1, 0}, #]]]] - N[vector[[fockIndex[#, 3]]]]] & /@ patterns[3, 2]] < 10^-12
        ]
    ],
    True,
    TestID -> "QI-symbolic-multimode-matches-Fock-circuit"
]

VerificationTest[
    {
        Map[Function[precision,
            With[{qi = QuantumInterferometer[{"BS"[N[3 / 10, precision], N[7 / 11, precision]] -> {1, 2}}, 2]},
                With[{value = qi["Amplitude", {nx, ny} -> {tx, nx + ny - tx}] /. {nx -> 2, ny -> 3, tx -> 1}, reference = qi["Amplitude", {2, 3} -> {1, 4}]},
                    If[precision === MachinePrecision,
                        Precision[value] === MachinePrecision && Abs[value - reference] < 10^-12,
                        Precision[value] > precision - 3 && Abs[value - reference] < 10 ^ (- precision + 5)]]]],
            {MachinePrecision, 60}],
        Map[Function[p,
            With[{phase = N[1/10^12, p]},
                With[{actual = QuantumInterferometer[{"PS"[phase] -> 1}, 1, "DropIdentities" -> False, Tolerance -> 0]["Amplitude", {nx} -> {nx}] /. nx -> 2, reference = Exp[2 I phase]},
                    SameQ[Precision[actual], Precision[reference]] && Abs[(actual - 1)/(reference - 1) - 1] < If[p === MachinePrecision, 10^-12, 10^-50]]]],
            {MachinePrecision, 60}]
    },
    {{True, True}, {True, True}},
    TestID -> "QI-symbolic-results-keep-numeric-precision"
]

(* summation indices never capture user symbols, even formal ones *)
VerificationTest[
    With[{parameter = Subscript[\[FormalK], 1]},
        With[{angle = QuantumInterferometer[{"BS"[parameter, 0] -> {1, 2}}, 2]["Amplitude", {1, 0} -> {tx, 1 - tx}],
              chain = QuantumInterferometer[{"BS"[Pi/6, Pi/7] -> {1, 2}, "BS"[Pi/5, Pi/11] -> {1, 2}}, 2]},
            FullSimplify[(angle /. {tx -> 1, parameter -> Pi/6}) - Sqrt[3]/2] === 0 &&
            And @@ Table[TrueQ[TimeConstrained[FullSimplify[
                Activate[chain["Amplitude", {parameter, 0} -> {tx, parameter - tx}] /. {parameter -> 2, tx -> t}] - chain["Amplitude", {2, 0} -> {t, 2 - t}]], 5, $TimedOut] === 0],
                {t, 0, 2}]
        ]
    ],
    True,
    TestID -> "QI-symbolic-sum-indices-do-not-capture-user-symbols"
]

(* on nonnegative integer photon numbers and real angles the results are bare formulas, with no
   condition restating that domain; photon conservation stays as a Kronecker delta *)
VerificationTest[
    With[{bs = QuantumInterferometer[{"BS"[\[Theta], \[Phi]] -> {1, 2}}, 2],
          hom = QuantumInterferometer[{"BS"[Pi/4, 0] -> {1, 2}}, 2]["Amplitude", {nx, nx} -> {kx, 2 nx - kx}],
          ps = QuantumInterferometer[{"PS"[\[Alpha]] -> 1}, 1],
          chain = QuantumInterferometer[{"BS"[Pi/6, Pi/7] -> {1, 2}, "BS"[Pi/5, Pi/11] -> {1, 2}}, 2]},
        With[{bs2 = bs["Amplitude", {nx, ny} -> {tx, nx + ny - tx}], general = bs["Amplitude", {nx, ny} -> {tx, ty}],
              sum = chain["Amplitude", {nx, 0} -> {tx, nx - tx}]},
            {
                MatchQ[bs2, _WignerD] && FullSimplify[bs2 - WignerD[{(nx + ny)/2, (nx - ny)/2, tx - (nx + ny)/2}, \[Phi], 2 \[Theta], - \[Phi]]] === 0,
                hom,
                {ps["Amplitude", {nx} -> {nx}], ps["Probability", {nx} -> {nx}], QuantumInterferometer[{}, 1]["Amplitude", {nx} -> {nx}]},
                MatchQ[general, KroneckerDelta[nx + ny, tx + ty] _WignerD] && (general /. {nx -> 1, ny -> 1, tx -> 1, ty -> 0}) === 0,
                MatchQ[sum, Inactive[Sum][_, {_K, 0, nx}]],
                FreeQ[{bs2, general, hom, sum}, If | Piecewise | ConditionalExpression],
                bs2 /. {nx -> 0, ny -> 0, tx -> 0},
                (* the photon distribution after the chain is normalized for every input *)
                Table[RootReduce[Total[Abs[Table[Activate[sum /. {nx -> n, tx -> t}], {t, 0, n}]] ^ 2]], {n, 0, 3}]
            }
        ]
    ],
    {True, WignerD[{nx, 0, kx - nx}, 0, Pi/2, 0], {E^(I \[Alpha] nx), 1, 1}, True, True, True, 1, {1, 1, 1, 1}},
    TestID -> "QI-symbolic-results-are-bare-formulas-on-the-physical-domain"
]

(* WignerD does not vanish outside its index range, so the formula keeps the one condition the
   photon numbers do not imply: photons in mode 3 reach mode 1 only through mode 2, after the
   first beam splitter has acted *)
VerificationTest[
    With[{exact = QuantumInterferometer[{"BS"[Pi/6, Pi/7] -> {1, 2}, "BS"[Pi/5, Pi/11] -> {2, 3}}, 3],
          machine = QuantumInterferometer[{"BS"[0.3, 0.7] -> {1, 2}, "BS"[1.1, 0.2] -> {2, 3}}, 3]},
        With[{a = exact["Amplitude", {n1, n2, n3} -> {t1, t2, t3}], b = machine["Amplitude", {n1, n2, n3} -> {t1, t2, t3}]},
            {
                MatchQ[a, KroneckerDelta[n1 + n2 + n3, t1 + t2 + t3] Piecewise[{{_, n1 + n2 >= t1}}, 0]],
                (a /. Thread[{n1, n2, n3, t1, t2, t3} -> {0, 0, 2, 1, 0, 1}]) === 0 && (b /. Thread[{n1, n2, n3, t1, t2, t3} -> {0, 0, 2, 1, 0, 1}]) === 0,
                And @@ Flatten @ Table[
                    FullSimplify[(a /. Thread[{n1, n2, n3, t1, t2, t3} -> Join[in, out]]) - exact["Amplitude", in -> out]] === 0 &&
                        Abs[(b /. Thread[{n1, n2, n3, t1, t2, t3} -> Join[in, out]]) - machine["Amplitude", in -> out]] < 10^-12,
                    {in, {{0, 1, 2}, {2, 0, 1}, {1, 1, 1}}}, {out, patterns[3, 3]}]
            }
        ]
    ],
    {True, True, True},
    TestID -> "QI-symbolic-photon-path-condition-is-kept"
]

VerificationTest[
    {
        QuantumInterferometer[{"S"[0.3] -> 1, "BS"[Pi / 4, 0] -> {1, 2}}, 2]["Amplitude", {nx, 0} -> {nx, 0}],
        QuantumInterferometer[{"BS"[Pi / 4, 0] -> {1, 2}}, 2]["Amplitude", {nx, -1} -> {nx, 0}]
    },
    {$Failed, $Failed},
    {QuantumInterferometer::squeezed, QuantumInterferometer::occ},
    TestID -> "QI-symbolic-photon-numbers-reject-squeezers-and-negatives"
]

(* symbolic amplitudes add no symbols to the caller's context *)
VerificationTest[
    With[{chain = QuantumInterferometer[{"BS"[Pi/6, Pi/7] -> {1, 2}, "BS"[Pi/5, Pi/11] -> {1, 2}}, 2], before = Names[$Context <> "*"]},
        chain["Amplitude", {nx, ny} -> {tx, ty}];
        chain["Probability", {nx, ny} -> {tx, ty}];
        Complement[Names[$Context <> "*"], before]
    ],
    {},
    TestID -> "QI-symbolic-amplitudes-create-no-symbols"
]

(* the summation indices K[j] skip past any K[j] already in the query *)
VerificationTest[
    With[{chain = QuantumInterferometer[{"BS"[Pi/6, Pi/7] -> {1, 2}, "BS"[Pi/5, Pi/11] -> {1, 2}}, 2]},
        With[{amplitude = chain["Amplitude", {K[1], 0} -> {K[2], K[1] - K[2]}]},
            {
                Cases[amplitude, Inactive[Sum][_, iters__] :> {iters}[[All, 1]], {0, Infinity}],
                Max[Abs[Table[N[Activate[amplitude /. {K[1] -> 2, K[2] -> t}] - chain["Amplitude", {2, 0} -> {t, 2 - t}]], {t, 0, 2}]]] < 10^-12,
                Cases[chain["Amplitude", {K[1][x], 0} -> {t1, t2}], Inactive[Sum][_, iters__] :> {iters}[[All, 1]], {0, Infinity}]
            }
        ]
    ],
    {{{K[3]}}, True, {{K[2]}}},
    TestID -> "QI-symbolic-sum-indices-skip-K-in-the-query",
    TimeConstraint -> 60
]

(* a K[j] used as an element parameter is skipped too, and only beam splitters take indices *)
VerificationTest[
    {
        Cases[QuantumInterferometer[{"BS"[K[5], 0] -> {1, 2}, "BS"[Pi/5, Pi/11] -> {1, 2}}, 2]["Amplitude", {n1, n2} -> {t1, t2}],
            Inactive[Sum][_, iters__] :> {iters}[[All, 1]], Infinity],
        Cases[QuantumInterferometer[{"PS"[al] -> 1, "BS"[th, ph] -> {1, 2}, "PS"[be] -> 2, "BS"[Pi/5, Pi/11] -> {1, 2}}, 2]["Amplitude", {n1, n2} -> {t1, t2}],
            Inactive[Sum][_, iters__] :> {iters}[[All, 1]], Infinity],
        With[{elements = {"BS"[K[1], 0] -> {1, 2}, "BS"[Pi/5, Pi/11] -> {1, 2}}},
            Max[Abs[Table[N[Activate[QuantumInterferometer[elements, 2]["Amplitude", {n1, n2} -> {t1, t2}] /. K[1] -> 3/10 /. {n1 -> 2, n2 -> 1, t1 -> t, t2 -> 3 - t}] -
                QuantumInterferometer[elements /. K[1] -> 3/10, 2]["Amplitude", {2, 1} -> {t, 3 - t}]], {t, 0, 3}]]] < 10^-12
        ]
    },
    {{{K[6]}}, {{K[1]}}, True},
    TestID -> "QI-symbolic-sum-indices-skip-K-in-the-elements",
    TimeConstraint -> 60
]

(* nested sums, whose inner bounds depend on the outer index, also skip a K[j] in the query *)
VerificationTest[
    With[{elements = {"BS"[th1, ph1] -> {1, 2}, "BS"[th2, ph2] -> {2, 3}, "BS"[th3, ph3] -> {1, 2}, "BS"[th4, ph4] -> {2, 3}},
          angles = {th1 -> 3/10, ph1 -> 7/10, th2 -> 11/10, ph2 -> 1/5, th3 -> 1/2, ph3 -> -2/5, th4 -> 9/10, ph4 -> 13/10}},
        With[{amplitude = QuantumInterferometer[elements, 3]["Amplitude", {K[1], K[2], 0} -> {t1, t2, t3}], numeric = QuantumInterferometer[elements /. angles, 3]},
            {
                Cases[amplitude, Inactive[Sum][_, iters__] :> {iters}[[All, 1]], Infinity],
                Max[Abs[Flatten @ Table[
                    N[Activate[amplitude /. angles /. {K[1] -> 2, K[2] -> 1, t1 -> a, t2 -> b, t3 -> 3 - a - b}] - numeric["Amplitude", {2, 1, 0} -> {a, b, 3 - a - b}]],
                    {a, 0, 3}, {b, 0, 3 - a}]]] < 10^-12
            }
        ]
    ],
    {{{K[3], K[4]}}, True},
    TestID -> "QI-symbolic-nested-sum-indices-skip-K-in-the-query"
]

(* a K[j] named only in the assumptions is skipped too, so the assumption does not constrain a sum *)
VerificationTest[
    With[{chain = QuantumInterferometer[{"BS"[Pi/6, Pi/7] -> {1, 2}, "BS"[Pi/5, Pi/11] -> {1, 2}}, 2],
          domain = Element[{n1, n2, t1, t2}, Integers] && n1 >= 0 && n2 >= 0 && t1 >= 0 && t2 >= 0},
        {
            Max[Abs[Table[N[Activate[chain["Amplitude", {n1, n2} -> {t1, t2}, Assumptions -> domain && K[1] < 0] /. {n1 -> 2, n2 -> 1, t1 -> t, t2 -> 3 - t}] - chain["Amplitude", {2, 1} -> {t, 3 - t}]], {t, 0, 3}]]] < 10^-12,
            Max[Abs[Table[N[Activate[Assuming[domain && K[1] < 0, chain["Probability", {n1, n2} -> {t1, t2}]] /. {n1 -> 2, n2 -> 1, t1 -> t, t2 -> 3 - t}] - chain["Probability", {2, 1} -> {t, 3 - t}]], {t, 0, 3}]]] < 10^-12
        }
    ],
    {True, True},
    TestID -> "QI-symbolic-sum-indices-skip-K-in-the-assumptions",
    TimeConstraint -> 60
]

EndTestSection[]


BeginTestSection["QuantumInterferometer: assumptions"]

VerificationTest[
    With[{hom = QuantumInterferometer[{"BS"[Pi/4, 0] -> {1, 2}}, 2], domain = Element[{nx, kx}, Integers] && nx >= 0 && 0 <= kx <= 2 nx, before = $Assumptions},
        {
            hom["Amplitude", {nx, nx} -> {kx, 2 nx - kx}, Assumptions -> domain],
            Assuming[domain, hom["Amplitude", {nx, nx} -> {kx, 2 nx - kx}]],
            $Assumptions === before,
            Assuming[nx == I, hom["Amplitude", {nx, nx} -> {kx, 2 nx - kx}, {Assumptions :> domain}]],
            Assuming[Element[{nx, kx}, Integers] && nx == 0 && kx == 0, hom["Amplitude", {nx, nx} -> {kx, 2 nx - kx}, Assumptions -> True]] ===
                Block[{$Assumptions = True}, hom["Amplitude", {nx, nx} -> {kx, 2 nx - kx}]],
            hom["Amplitude", {nx, nx} -> {kx, 2 nx - kx}]
        }
    ],
    {WignerD[{nx, 0, kx - nx}, 0, Pi/2, 0], WignerD[{nx, 0, kx - nx}, 0, Pi/2, 0], True, WignerD[{nx, 0, kx - nx}, 0, Pi/2, 0], True, WignerD[{nx, 0, kx - nx}, 0, Pi/2, 0]},
    TestID -> "QI-assumptions-give-the-compact-formula"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[theta, phi] -> {1, 2}}, 2]},
        With[{a = qi["Amplitude", {nx, ny} -> {kx, nx + ny - kx}, Assumptions -> Element[{nx, ny, kx}, Integers] && nx >= 0 && ny >= 0 && 0 <= kx <= nx + ny]},
            MatchQ[a, _WignerD] &&
            FullSimplify[a - WignerD[{(nx + ny)/2, (nx - ny)/2, kx - (nx + ny)/2}, phi, 2 theta, -phi]] === 0 &&
            (a /. {nx -> 0, ny -> 0, kx -> 0}) === 1 &&
            And @@ Flatten @ Table[FullSimplify[(a /. {nx -> n, ny -> 0, kx -> k, theta -> Pi/6, phi -> Pi/7}) -
                QuantumInterferometer[{"BS"[Pi/6, Pi/7] -> {1, 2}}, 2]["Amplitude", {n, 0} -> {k, n - k}]] === 0, {n, {1, 3}}, {k, 0, n}]
        ]
    ],
    True,
    TestID -> "QI-assumptions-general-two-mode-identity"
]

(* restating the physical domain changes nothing, a relation beyond it simplifies further, and
   assumptions that contradict it fail *)
VerificationTest[
    With[{hom = QuantumInterferometer[{"BS"[Pi/4, 0] -> {1, 2}}, 2], ps = QuantumInterferometer[{"PS"[alpha] -> 1}, 1],
          chain = QuantumInterferometer[{"BS"[Pi/6, Pi/7] -> {1, 2}, "BS"[Pi/5, Pi/11] -> {1, 2}}, 2]},
        {
            hom["Amplitude", {nx, nx} -> {kx, 2 nx - kx}, Assumptions -> Element[{nx, kx}, Integers] && nx >= 0 && kx >= 0] ===
                hom["Amplitude", {nx, nx} -> {kx, 2 nx - kx}],
            chain["Amplitude", {nx, 0} -> {kx, nx - kx}, Assumptions -> Element[nx, Integers] && nx >= 0] === chain["Amplitude", {nx, 0} -> {kx, nx - kx}],
            ps["Amplitude", {nx} -> {nx}, Assumptions -> alpha == Pi],
            hom["Amplitude", {nx, nx} -> {kx, 2 nx - kx}, Assumptions -> nx == 1/2],
            ps["Amplitude", {nx} -> {nx}, Assumptions -> alpha == I]
        }
    ],
    {True, True, (-1)^nx, $Failed, $Failed},
    {Refine::fas, Refine::fas},
    TestID -> "QI-assumptions-add-relations-beyond-the-physical-domain"
]

(* a list of assumptions reads as their conjunction, and an assumed endpoint transmissivity
   gives the endpoint law rather than 0^0 *)
VerificationTest[
    Module[{eta, n, k},
        With[{ps = QuantumInterferometer[{"PS"[alpha] -> 1}, 1], loss = QuantumInterferometer[{"Loss"[eta]}]},
            {
                ps["Amplitude", {nx} -> {nx}, Assumptions -> {alpha == Pi}],
                loss["Probability", {n} -> {k}, Assumptions -> eta == 0] === KroneckerDelta[k, 0],
                loss["Probability", {n} -> {k}, Assumptions -> {eta == 1}] === KroneckerDelta[n, k],
                Assuming[eta == 0, loss["Probability", {n} -> {k}]] === KroneckerDelta[k, 0]
            }
        ]
    ],
    {(-1)^nx, True, True, True},
    TestID -> "QI-assumptions-list-form-and-endpoint-transmissivities"
]

(* a non-numeric photon number that cannot be a nonnegative integer is rejected like a concrete one *)
VerificationTest[
    With[{hom = QuantumInterferometer[{"BS"[Pi/4, 0] -> {1, 2}}, 2]},
        {hom["Amplitude", {nx, Infinity} -> {tx, 0}], hom["Amplitude", {nx, "a"} -> {tx, 0}]}
    ],
    {$Failed, $Failed},
    {QuantumInterferometer::occ, QuantumInterferometer::occ},
    TestID -> "QI-symbolic-photon-numbers-reject-infinite-and-string-entries"
]

(* Null and None are rejected too, and ::occ names the pattern at fault *)
VerificationTest[
    With[{hom = QuantumInterferometer[{"BS"[Pi/4, 0] -> {1, 2}}, 2]},
        {hom["Amplitude", {nx, Null} -> {tx, 0}], QuantumInterferometer[{"Loss"[1/2]}]["Probability", {nx} -> {None}]}
    ],
    {$Failed, $Failed},
    {HoldForm[Message[QuantumInterferometer::occ, {nx, Null}, 2]], HoldForm[Message[QuantumInterferometer::occ, {None}, 1]]},
    TestID -> "QI-symbolic-photon-numbers-name-the-rejected-pattern"
]

VerificationTest[
    With[{hom = QuantumInterferometer[{"BS"[Pi/4, 0] -> {1, 2}}, 2]},
        hom["Amplitude", # -> {tx, 0}] & /@ {{nx, {1}}, {nx, Interval[{0, 1}]}}
    ],
    {$Failed, $Failed},
    {QuantumInterferometer::occ, QuantumInterferometer::occ},
    TestID -> "QI-symbolic-photon-numbers-reject-lists-and-placeholder-symbols"
]

(* the inverse, photon distributions and lossy states are simplified on the physical domain too *)
VerificationTest[
    Module[{eta},
        With[{bs = QuantumInterferometer[{"BS"[th, ph]}], loss = QuantumInterferometer[{"Loss"[eta]}]},
            {
                bs["Dagger"]["Unitary"],
                Values[bs["Probabilities", {1, 1}]],
                Normal[loss[QuantumState[{1, I}/Sqrt[2], 2]]["DensityMatrix"]] === {{1/2 + (1 - eta)/2, -I Sqrt[eta]/2}, {I Sqrt[eta]/2, eta/2}},
                FreeQ[Values[bs["Probabilities", "Coherent"[{a, 0}], "Cutoff" -> 1]], Im | Re],
                With[{q = QuantumInterferometer["Random"[3]]}, q["Dagger"]["Unitary"] === ConjugateTranspose[q["Unitary"]]]
            }
        ]
    ],
    {
        {{Cos[th], Sin[th]/E^(I ph)}, {-E^(I ph) Sin[th], Cos[th]}},
        {2 Cos[th]^2 Sin[th]^2, (Cos[th]^2 - Sin[th]^2)^2, 2 Cos[th]^2 Sin[th]^2},
        True, True, True
    },
    TestID -> "QI-dagger-probabilities-and-lossy-states-on-the-physical-domain"
]

(* the all-outputs form names symbolic photon numbers instead of staying unevaluated *)
VerificationTest[
    QuantumInterferometer[{"BS"[th, ph]}][{n1, n2}],
    $Failed,
    {QuantumInterferometer::occ},
    TestID -> "QI-symbolic-all-outputs-form-rejects-symbolic-photon-numbers"
]

(* the overlap phase of two coherent states has unit modulus, so probabilities are plain Gaussians *)
VerificationTest[
    Module[{eta},
        With[{loss = QuantumInterferometer[{"Loss"[eta]}], bs = QuantumInterferometer[{"BS"[th, ph]}]},
            {
                loss["Probability", "Coherent"[{a}] -> "Coherent"[{b}]] === E^(-Abs[b - a Sqrt[eta]]^2),
                FreeQ[bs["Probability", "Coherent"[{a, b}] -> "Coherent"[{c, d}]], Re | Im]
            }
        ]
    ],
    {True, True},
    TestID -> "QI-coherent-overlap-probabilities-carry-no-phase-term"
]

(* the Gaussian kernel blocks are (alpha^dag)^-1 and its products, so squeezed amplitudes come out compact *)
VerificationTest[
    {
        QuantumInterferometer[{"S"[xi]}]["Amplitude", {1} -> {3}],
        With[{elements = {"S"[x1] -> 1, "S"[x2] -> 2, "BS"[Pi/5, 1/3]}, values = {x1 -> 1/3 + I/4, x2 -> -1/2}},
            Abs[N[QuantumInterferometer[elements]["Amplitude", {1, 0} -> {1, 2}] /. values] -
                N[QuantumInterferometer[elements /. values]["Amplitude", {1, 0} -> {1, 2}]]] < 10^-12
        ]
    },
    {-Sqrt[3/2] E^(I Arg[xi]) Sech[Abs[xi]]^(5/2) Sinh[Abs[xi]], True},
    TestID -> "QI-squeezed-symbolic-amplitudes-use-the-closed-form-kernel"
]

(* states through the Fock-space circuit are simplified on the physical domain; numeric ones are untouched *)
VerificationTest[
    With[{bs = QuantumInterferometer[{"BS"[th, ph]}], rho = QuantumState[DiagonalMatrix[{0, 1/2, 1/2, 0}], {2, 2}]},
        {
            Normal[bs[rho]["DensityMatrix"]],
            FreeQ[Normal[QuantumInterferometer[{"S"[xi]}][QuantumState[DiagonalMatrix[{1/2, 1/2, 0}], 3]]["DensityMatrix"]], Conjugate],
            With[{q = QuantumInterferometer[{"BS"[1/3, 1/5]}]}, Normal[q[rho]["DensityMatrix"]] === Normal[q["CircuitOperator", 2][rho]["DensityMatrix"]]]
        }
    ],
    {
        {{0, 0, 0, 0}, {0, (Cos[th]^2 + Sin[th]^2)/2, 0, 0}, {0, 0, (Cos[th]^2 + Sin[th]^2)/2, 0}, {0, 0, 0, 0}},
        True, True
    },
    TestID -> "QI-mixed-and-squeezed-states-on-the-physical-domain"
]

(* a pattern that breaks a selection rule gives 0 with no message *)
VerificationTest[
    {
        QuantumInterferometer[{"BS"[Pi/4, 0] -> {1, 2}}, 2]["Amplitude", {nx, 0} -> {nx + 1, 0}],
        QuantumInterferometer[{"BS"[Pi/6, Pi/7] -> {2, 3}, "PS"[Pi/5] -> 1}, 3]["Amplitude", {nx, 1, 0} -> {nx + 1, uy, uz}],
        QuantumInterferometer[{"BS"[0.3, 0.7] -> {1, 2}}, 2]["Amplitude", {nx, 0} -> {nx + 1, 0}]
    },
    {0, 0, 0},
    TestID -> "QI-symbolic-off-shell-patterns-are-zero-without-messages"
]

(* inexact angles: a numeric WignerD outside its index range is Indeterminate, so the selection
   rules are held in the Piecewise condition and off-shell patterns give 0 *)
VerificationTest[
    With[{bs = QuantumInterferometer[{"BS"[0.3, 0.7] -> {1, 2}}, 2], side = QuantumInterferometer[{"BS"[0.3, 0.7] -> {2, 3}}, 3]},
        With[{f = bs["Amplitude", {n1, n2} -> {t1, t2}], g = side["Amplitude", {n1, n2, n3} -> {t1, t2, t3}]},
            {
                MatchQ[f, Piecewise[{{_WignerD, n1 + n2 == t1 + t2}}, 0]],
                {f /. {n1 -> 0, n2 -> 0, t1 -> 1, t2 -> 0}, f /. {n1 -> 1, n2 -> 0, t1 -> 2, t2 -> 0},
                    bs["Probability", {n1, n2} -> {t1, t2}] /. {n1 -> 1, n2 -> 0, t1 -> 2, t2 -> 0},
                    g /. Thread[{n1, n2, n3, t1, t2, t3} -> {1, 0, 0, 0, 1, 0}]},
                Max[Abs[Table[(f /. Thread[{n1, n2, t1, t2} -> Join[{2, 1}, p]]) - bs["Amplitude", {2, 1} -> p], {p, patterns[2, 3]}]]] < 10^-12
            }
        ]
    ],
    {True, {0, 0, 0, 0}, True},
    TestID -> "QI-symbolic-inexact-angles-hold-the-selection-rules"
]

(* the modulus squared of a selection rule is the rule, and self-contradictory assumptions fail *)
VerificationTest[
    Module[{eta, n, k},
        {
            QuantumInterferometer[{"PS"[alpha] -> 1}, 1]["Probability", {nx} -> {tx}],
            QuantumInterferometer[{"Loss"[eta]}]["Probability", {n} -> {k}, Assumptions -> {eta == 0, eta == 1}]
        }
    ],
    {KroneckerDelta[nx, tx], $Failed},
    {Refine::cas},
    TestID -> "QI-assumptions-contradictory-and-squared-selection-rules"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[Pi/6, Pi/7] -> {2, 3}, "PS"[Pi/5] -> 1}, 3]},
        With[{a = qi["Amplitude", {nx, 1, 0} -> {tx, uy, uz}, Assumptions -> Element[{nx, tx, uy, uz}, Integers] && Min[nx, tx, uy, uz] >= 0]},
            {a /. {nx -> 2, tx -> 1, uy -> 1, uz -> 1}, a /. {nx -> 2, tx -> 2, uy -> 0, uz -> 0},
                FullSimplify[(a /. {nx -> 2, tx -> 2, uy -> 0, uz -> 1}) - qi["Amplitude", {2, 1, 0} -> {2, 0, 1}]]}
        ]
    ],
    {0, 0, 0},
    TestID -> "QI-assumptions-preserve-photon-conservation"
]

VerificationTest[
    Map[Function[p,
        With[{qi = QuantumInterferometer[{"BS"[N[3/10, p], N[7/11, p]] -> {1, 2}}, 2]},
            With[{v = qi["Amplitude", {nx, ny} -> {kx, nx + ny - kx}, Assumptions -> Element[{nx, ny, kx}, Integers] && nx >= 0 && ny >= 0 && 0 <= kx <= nx + ny] /. {nx -> 2, ny -> 3, kx -> 1},
                  old = qi["Amplitude", {nx, ny} -> {kx, nx + ny - kx}] /. {nx -> 2, ny -> 3, kx -> 1}, direct = qi["Amplitude", {2, 3} -> {1, 4}]},
                v === old && Precision[v] === Precision[old] && Abs[v - direct] < If[p === MachinePrecision, 10^-12, 10^(-p + 5)] &&
                qi["Amplitude", {2, 3} -> {1, 4}, Assumptions -> nx > 0] === direct
            ]
        ]], {MachinePrecision, 80}],
    {True, True},
    TestID -> "QI-assumptions-preserve-precision"
]

VerificationTest[
    With[{ps = QuantumInterferometer[{"PS"[alpha] -> 1}, 1], identity = QuantumInterferometer[{}, 1]},
        {
            ps["Probability", {2} -> {2}, Assumptions -> Element[alpha, Reals]],
            FullSimplify[identity["Amplitude", "Coherent"[{drive}] -> {1}, Assumptions -> Im[drive] > 0] /. drive -> (1 + I)/3],
            FullSimplify[identity["Probability", "Coherent"[{drive}] -> {1}, Assumptions -> Element[theta, Reals]] /. drive -> (1 + I)/3]
        }
    ],
    {1, (1 + I) Exp[-1/9]/3, 2 Exp[-2/9]/9},
    TestID -> "QI-assumptions-on-angles-and-coherent-amplitudes"
]

VerificationTest[
    {
        QuantumInterferometer[{"BS"[Pi/4, 0] -> {1, 2}}, 2]["Amplitude", {1/2, 0} -> {1/2, 0}, Assumptions -> True],
        QuantumInterferometer[{"S"[1/5] -> 1}, 1]["Amplitude", {nx} -> {nx}, Assumptions -> Element[nx, Integers] && nx >= 0],
        QuantumInterferometer[{"BS"[Pi/4, 0] -> {1, 2}}, 2]["Amplitude", {nx, nx} -> {kx, 2 nx - kx}, Assumptions -> False]
    },
    {$Failed, $Failed, $Failed},
    {QuantumInterferometer::occ, QuantumInterferometer::squeezed, Refine::fas},
    TestID -> "QI-assumptions-cannot-bypass-checks"
]

EndTestSection[]


BeginTestSection["QuantumInterferometer: squeezers"]

hafnian = Wolfram`QuantumFramework`SecondQuantization`Interferometer`PackagePrivate`hafnian;

(* all-ones matrices count perfect matchings, (n - 1)!!; Haf([[0, G], [G^T, 0]]) = Perm(G) on the pairing and power-trace paths *)
VerificationTest[
    SeedRandom[21];
    With[{g = RandomComplex[{-1 - I, 1 + I}, {6, 6}], h = RandomInteger[{-3, 3}, {3, 3}]},
        {
            hafnian[ConstantArray[1, {#, #}]] & /@ {2, 4, 6, 8},
            Chop[hafnian[ConstantArray[1., {10, 10}]] - 945],
            Chop[hafnian[ConstantArray[N[1, 30], {10, 10}]] - 945],
            hafnian[ArrayFlatten[{{0, h}, {Transpose[h], 0}}]] == Permanent[h],
            Abs[hafnian[ArrayFlatten[{{0, g}, {Transpose[g], 0}}]] - Permanent[g]] < 10^-10
        }
    ],
    {{1, 3, 15, 105}, 0, 0, True, True},
    TestID -> "QI-hafnian"
]

(* the mesh conserves photon number, so the Fock circuit is exact below d photons *)
VerificationTest[
    SeedRandom[22];
    With[{qi = QuantumInterferometer[Join[{"S"[0.5] -> 1, "S"[0.3 Exp[I]] -> 2, "S"[0.4] -> 3}, QuantumInterferometer["Random"[3]]["Elements"]], 3],
          outs = Select[Tuples[Range[0, 6], 3], Total[#] <= 6 &]},
        Max[Abs[(qi["Amplitude", {0, 0, 0} -> #] & /@ outs) - fockAmplitudes[qi["CircuitOperator", 7], {0, 0, 0}, 7, outs]]] < 10^-12
    ],
    True,
    TestID -> "QI-squeezers-then-mesh-match-Fock-circuit"
]

(* squeezers between beam splitters, with photons in: the truncated circuit converges to the hafnians *)
VerificationTest[
    With[{qi = QuantumInterferometer[{"S"[0.4 Exp[0.5 I]] -> 1, "BS"[0.7, 0.3] -> {1, 2}, "S"[0.3 Exp[-1.2 I]] -> 2, "BS"[0.4, 1.] -> {1, 2}}, 2],
          outs = Tuples[Range[0, 4], 2]},
        Max @ Table[Abs[(qi["Amplitude", in -> #] & /@ outs) - fockAmplitudes[qi["CircuitOperator", 40], in, 40, outs]], {in, {{0, 0}, {1, 0}, {0, 2}}}] < 10^-9
    ],
    True,
    TestID -> "QI-interleaved-squeezers-match-Fock-circuit"
]

(* S(r) and S(-r) on a 50:50 splitter give the two-mode squeezed vacuum, Sum_n Tanh[r]^n |n, n> / Cosh[r] *)
VerificationTest[
    With[{qi = QuantumInterferometer[{"S"[0.5] -> 1, "S"[-0.5] -> 2, "BS"[Pi / 4, 0] -> {1, 2}}, 2]},
        {
            With[{outs = Select[Tuples[Range[0, 6], 2], Total[#] <= 6 && EvenQ[Total[#]] &], probs = qi["Probabilities", {0, 0}, "Cutoff" -> 6]},
                Sort[Keys[probs]] === Sort[QuditName /@ outs] &&
                    Max[Abs[Lookup[probs, QuditName /@ outs] - (If[Equal @@ #, Tanh[.5] ^ (2 First[#]) / Cosh[.5] ^ 2, 0] & /@ outs)]] < 10^-12 &&
                    Abs[Total[probs] - (1 - Tanh[.5] ^ 8)] < 10^-12
            ],
            qi["Unitary"],
            qi["SqueezerCount"]
        }
    ],
    {True, $Failed, 2},
    {QuantumInterferometer::squeezed},
    TestID -> "QI-two-mode-squeezed-vacuum"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"S"[0.4 Exp[0.5 I]] -> 1, "BS"[0.7, 0.3] -> {1, 2}, "S"[0.3 Exp[-1.2 I]] -> 2, "PS"[0.9] -> 2}, 2]},
        {
            Max[Abs[Table[qi["Dagger"][qi]["Amplitude", in -> out], {in, {{0, 0}, {1, 0}, {1, 1}}}, {out, {{0, 0}, {1, 0}, {1, 1}}}] - IdentityMatrix[3]]] < 10^-12,
            With[{probs = qi["Probabilities", {1, 0}]},
                1 - 10^-4 < Total[probs] <= 1 + 10^-12 && AnyTrue[Keys[probs], Total[#["Name"]] > 1 &] &&
                    Max[Abs[Values[probs] - (qi["Probability", {1, 0} -> #["Name"]] & /@ Keys[probs])]] < 10^-12
            ]
        }
    ],
    {True, True},
    TestID -> "QI-squeezers-dagger-and-auto-cutoff"
]

VerificationTest[
    {
        With[{probs = QuantumInterferometer[{"S"[.6] -> 1}, 1]["Probabilities", {0}]},
            With[{top = Max[First[Flatten[{#["Name"]}]] & /@ Keys[probs]]},
                With[{outs = List /@ Range[0, top, 2], ref = Table[Binomial[2 k, k] Tanh[.6] ^ (2 k) / (4 ^ k Cosh[.6]), {k, 0, top / 2}]},
                    top > 0 && Sort[Keys[probs]] === Sort[QuditName /@ outs] && Max[Abs[Lookup[probs, QuditName /@ outs] - ref]] < 10^-12 &&
                        0 < 1 - Total[probs] < 10^-4 && Total[Most[ref]] < 1 - 10^-4]]],
        With[{probs = QuantumInterferometer[{"S"[1.2] -> 1}, 1]["Probabilities", {0}], ref = Table[Binomial[2 k, k] Tanh[1.2] ^ (2 k) / (4 ^ k Cosh[1.2]), {k, 0, 8}]},
            Sort[Keys[probs]] === Sort[QuditName /@ List /@ Range[0, 16, 2]] && Max[Abs[Lookup[probs, QuditName /@ List /@ Range[0, 16, 2]] - ref]] < 10^-12 &&
                0 < Total[probs] < 1 - 10^-4]
    },
    {True, True},
    {QuantumInterferometer::autocutoff},
    TestID -> "QI-squeezed-auto-cutoff-keeps-absolute-weights"
]

(* S(xi)|0> and S(xi)|1> at high photon number, across the generic and repeated-row hafnian routes *)
VerificationTest[
    And @@ Flatten @ Table[
        With[{r = If[p === MachinePrecision, N[2/5], N[2/5, p]], phi = If[p === MachinePrecision, N[7/13], N[7/13, p]]},
            With[{qi = QuantumInterferometer[{"S"[r Exp[I phi]] -> 1}, 1]},
                Table[
                    With[{n = 16 + in + 2 j},
                        With[{ref = Sqrt[n!] (- Exp[I phi] Tanh[r]) ^ ((n - in) / 2) / (2 ^ ((n - in) / 2) ((n - in) / 2)! Cosh[r] ^ (in + 1/2)),
                              actual = qi["Amplitude", {in} -> {n}]},
                            Abs[ref] > 10^-8 && Abs[actual - ref] < If[p === MachinePrecision, 10^-11, 10^-30] Abs[ref]
                        ]
                    ],
                    {in, {0, 1}}, {j, 0, 2}
                ]
            ]
        ],
        {p, {MachinePrecision, 40}}
    ],
    True,
    TestID -> "QI-high-photon-squeezing-precision"
]

VerificationTest[
    With[{d = 8, xi = .25 Exp[.4 I]},
        With[{v = (FockState[{0, 0}, d] + I FockState[{1, 1}, d])/Sqrt[2], u = FockState[{1, 0}, d],
              qi = QuantumInterferometer[{"S"[xi] -> 1, "BS"[.37, .29] -> {1, 2}, "PS"[.3] -> 2}, 2],
              circuit = PhaseShiftOperator[.3, d, {2}] @ BeamSplitterOperator[{.37, .29}, d, {1, 2}, Method -> "Recurrence"] @ SqueezeOperator[xi, d, {1}]},
            With[{rho = QuantumState[.6 Normal[v["DensityMatrix"]] + .4 Normal[u["DensityMatrix"]], {d, d}]},
                {
                    Max[Abs[Normal[qi[v]["StateVector"]] - Normal[circuit[v]["StateVector"]]]] < 10^-12,
                    Max[Abs[Normal[qi[rho]["DensityMatrix"]] - Normal[circuit[rho]["DensityMatrix"]]]] < 10^-12
                }
            ]
        ]
    ],
    {True, True},
    TestID -> "QI-active-states-match-Fock-operators"
]

EndTestSection[]


BeginTestSection["QuantumInterferometer: coherent light"]

$coherentW = DiagonalMatrix[{1, I}] . {{3/5, -4 Exp[-I Pi/3]/5}, {4 Exp[I Pi/3]/5, 3/5}};

VerificationTest[
    With[{out = QuantumInterferometer[$coherentW]["Coherent"[{a, b}]]},
        MatchQ[out, "Coherent"[_List]] && FullSimplify[First[out] - {3 a/5 - 4 Exp[-I Pi/3] b/5, 4 I Exp[I Pi/3] a/5 + 3 I b/5}] === {0, 0}
    ],
    True,
    TestID -> "QI-coherent-symbolic-transfer"
]

VerificationTest[
    With[{w = N[$coherentW], a = N[{(2 + I)/7, (-1 + 2 I)/6}], outs = Select[Tuples[Range[0, 6], 2], Total[#] <= 6 &]},
        With[{qi = QuantumInterferometer[w], beta = w . a},
            With[{probs = qi["Probabilities", "Coherent"[a], "Cutoff" -> 6], ref = Exp[-Total[Abs[a] ^ 2]/2] (Times @@ (beta ^ # / Sqrt[#!]) & /@ outs)},
                Sort[Keys[probs]] === Sort[QuditName /@ outs] &&
                    Max[Abs[(qi["Amplitude", "Coherent"[a] -> #] & /@ outs) - ref]] < 10^-12 &&
                    Max[Abs[Lookup[probs, QuditName /@ outs] - Abs[ref] ^ 2]] < 10^-12 &&
                    Abs[Total[probs] - Exp[-Total[Abs[a] ^ 2]] Sum[Total[Abs[a] ^ 2] ^ k / k!, {k, 0, 6}]] < 10^-12
            ]
        ]
    ],
    True,
    TestID -> "QI-coherent-Fock-amplitudes-match-Poisson"
]

VerificationTest[
    With[{w = N[$coherentW], a = N[{(2 + I)/7, (-1 + 2 I)/6}], g = N[{(1 + 2 I)/10, (-2 + I)/10}], alpha = .4 + .2 I},
        With[{qi = QuantumInterferometer[w], inputs = {"Coherent"[a], {1, 1}, {2, "Coherent"[alpha]}},
              refs = {
                Exp[-(Total[Abs[a] ^ 2] + Total[Abs[g] ^ 2])/2 + Conjugate[g] . w . a],
                Exp[-Total[Abs[g] ^ 2]/2] Times @@ (Conjugate[g] . w),
                Exp[-(Abs[alpha] ^ 2 + Total[Abs[g] ^ 2])/2 + alpha Conjugate[g] . w[[All, 2]]] (Conjugate[g] . w[[All, 1]]) ^ 2 / Sqrt[2]
              }},
            Max[Abs[MapThread[qi["Amplitude", #1 -> "Coherent"[g]] - #2 &, {inputs, refs}]]] < 10^-12 &&
                Max[Abs[MapThread[qi["Probability", #1 -> "Coherent"[g]] - Abs[#2] ^ 2 &, {inputs, refs}]]] < 10^-12
        ]
    ],
    True,
    TestID -> "QI-coherent-output-overlaps"
]

(* |2,alpha> maps to (W_1.a^dag)^2 Exp[alpha W_2.a^dag]|0>/Sqrt[2]; a degree-four polynomial gives every coefficient through six photons *)
VerificationTest[
    With[{w = N[$coherentW], alpha = .4 + .2 I, outs = Select[Tuples[Range[0, 6], 2], 2 <= Total[#] <= 6 &]},
        Module[{x, y},
            With[{qi = QuantumInterferometer[w], polynomial = Expand[({x, y} . w[[All, 1]]) ^ 2 Sum[(alpha {x, y} . w[[All, 2]]) ^ k/k!, {k, 0, 4}]/Sqrt[2]]},
                With[{ref = (Exp[-Abs[alpha] ^ 2/2] Sqrt[Times @@ (#!)] Coefficient[Coefficient[polynomial, x, #[[1]]], y, #[[2]]] & /@ outs),
                      probs = qi["Probabilities", {2, "Coherent"[alpha]}, "Cutoff" -> 6]},
                    Sort[Keys[probs]] === Sort[QuditName /@ outs] &&
                        Max[Abs[(qi["Amplitude", {2, "Coherent"[alpha]} -> #] & /@ outs) - ref]] < 10^-12 &&
                        Max[Abs[Lookup[probs, QuditName /@ outs] - Abs[ref] ^ 2]] < 10^-12 &&
                        Max[Abs[qi["Amplitude", {2, "Coherent"[alpha]} -> #] & /@ {{0, 0}, {1, 0}, {0, 1}}]] < 10^-12
                ]
            ]
        ]
    ],
    True,
    TestID -> "QI-mixed-Fock-coherent-input"
]

VerificationTest[
    {
        With[{probs = QuantumInterferometer[{}, 1]["Probabilities", "Coherent"[{(3 + 2 I)/5}]], lambda = Abs[(3 + 2 I)/5] ^ 2},
            With[{top = Max[First[Flatten[{#["Name"]}]] & /@ Keys[probs]]},
                With[{ref = Table[Exp[-lambda] lambda ^ n/n!, {n, 0, top}]},
                    top > 0 && Max[Abs[Lookup[probs, QuditName /@ List /@ Range[0, top]] - ref]] < 10^-12 && 1 - 10^-4 <= Total[probs] < 1 && Total[Most[ref]] < 1 - 10^-4]]],
        FullSimplify[Values[QuantumInterferometer[{}, 1]["Probabilities", "Coherent"[{alpha}], "Cutoff" -> 3]] - Table[Exp[-Abs[alpha] ^ 2] Abs[alpha] ^ (2 n)/n!, {n, 0, 3}]]
    },
    {True, ConstantArray[0, 4]},
    TestID -> "QI-coherent-cutoffs"
]

VerificationTest[
    {
        QuantumInterferometer[{"BS"[.4, .2] -> {1, 2}}, 2]["Coherent"[{.3}]],
        QuantumInterferometer[{"BS"[.4, .2] -> {1, 2}}, 2]["Amplitude", "Coherent"[{.3}] -> {0, 0}],
        QuantumInterferometer[{"S"[.3] -> 1}, 1]["Coherent"[{.2}]],
        QuantumInterferometer[{"S"[.3] -> 1}, 1]["Probabilities", "Coherent"[{.2}], "Cutoff" -> 4],
        QuantumInterferometer[{}, 1]["Probabilities", "Coherent"[{alpha}]],
        QuantumInterferometer[{"S"[r] -> 1}, 1]["Probabilities", {0}]
    },
    ConstantArray[$Failed, 6],
    {QuantumInterferometer::occ, QuantumInterferometer::occ, QuantumInterferometer::squeezed, QuantumInterferometer::squeezed,
        QuantumInterferometer::symcutoff, QuantumInterferometer::symcutoff},
    TestID -> "QI-coherent-errors"
]

EndTestSection[]


BeginTestSection["QuantumInterferometer: parameters"]

VerificationTest[
    With[{tm = QuantumInterferometer[{"S"[r] -> 1, "S"[-r] -> 2, "BS"[Pi/4, 0] -> {1, 2}}, 2, "Parameters" -> {r}],
          bs = QuantumInterferometer[{"BS"[\[Theta], \[Phi]] -> {1, 2}}, 2, "Parameters" -> {\[Theta], \[Phi]}]},
        {
            tm["Parameters"],
            tm[0.5] === QuantumInterferometer[{"S"[0.5] -> 1, "S"[-0.5] -> 2, "BS"[Pi/4, 0] -> {1, 2}}, 2] === tm[<|r -> 0.5|>],
            tm[0]["Elements"],
            {bs[Pi/4]["Parameters"], bs[<|\[Phi] -> 0|>]["Elements"]},
            {bs["Dagger"]["Parameters"], bs[bs]["Parameters"]},
            QuantumInterferometer[{"S"[r] -> 1}, 1, "Parameters" -> {r}, "DropIdentities" -> False][0]["Elements"],
            MatchQ[QuantumInterferometer[{"BS"[Pi/4, 0] -> {1, 2}}, 2][0.8], _QuantumInterferometer[0.8]]
        }
    ],
    {{r}, True, {"BS"[Pi/4, 0] -> {1, 2}}, {{\[Phi]}, {"BS"[\[Theta], 0] -> {1, 2}}}, {{\[Theta], \[Phi]}, {\[Theta], \[Phi]}}, {"S"[0] -> {1}}, True},
    TestID -> "QI-parameters-substitute-values"
]

EndTestSection[]


BeginTestSection["QuantumInterferometer: Bloch-Messiah and GBSForm"]

VerificationTest[
    With[{nets = {
        QuantumInterferometer[{"S"[0.4] -> 1, "BS"[0.7, 0.3] -> {1, 2}, "S"[0.3 I] -> 2, "BS"[1.1, -0.5] -> {2, 3}, "S"[-0.5] -> 3, "PS"[0.9] -> 1, "BS"[0.4, 1.2] -> {1, 3}}, 3],
        (SeedRandom[7]; QuantumInterferometer[Join[Table["S"[0.5] -> k, {k, 3}], QuantumInterferometer["Random"[3]]["Elements"]], 3])
    }, pats = Catenate[FrobeniusSolve[{1, 1, 1}, #] & /@ Range[0, 3]]},
        {Max[Abs[Flatten @ Table[qi["Probability", s -> t] - qi["BlochMessiah"]["Probability", s -> t], {qi, nets}, {s, Take[pats, 4]}, {t, pats}]]] < 10^-10,
         Chop[Cases[nets[[2]]["BlochMessiah"]["Elements"], ("S"[x_] -> _) :> x] - 0.5, 10^-10]}
    ],
    {True, {0, 0, 0}},
    TestID -> "QI-BlochMessiah-reproduces-the-network-including-equal-squeezing"
]

VerificationTest[
    With[{nets = {
        QuantumInterferometer[{"S"[0.5] -> 1, "S"[-0.5] -> 2, "BS"[Pi / 4, 0] -> {1, 2}}, 2],
        (SeedRandom[3]; QuantumInterferometer[Join[{"S"[0.5] -> 1, "S"[0.5] -> 2, "S"[0.3] -> 3}, QuantumInterferometer["Random"[4]]["Elements"]], 4])
    }},
        {Table[With[{m = qi["Modes"]},
            Max[Abs[(qi["Probability", ConstantArray[0, m] -> #] - qi["GBSForm"]["Probability", ConstantArray[0, m] -> #]) & /@ Catenate[FrobeniusSolve[ConstantArray[1, m], #] & /@ Range[0, 4]]]] < 10^-10], {qi, nets}],
         With[{passive = QuantumInterferometer[{"BS"[Pi / 4, 0] -> {1, 2}}, 2]}, passive["GBSForm"] === passive["BlochMessiah"] === passive],
         QuantumInterferometer[{"S"[r] -> 1}, 1]["GBSForm"]}
    ],
    {{True, True}, True, $Failed},
    {QuantumInterferometer::notnumeric},
    TestID -> "QI-GBSForm-puts-squeezers-first-for-vacuum-input"
]

EndTestSection[]


BeginTestSection["QuantumInterferometer: drawing"]

VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[.4, .2] -> {1, 3}, "PS"[.3] -> 2}, 3]},
        With[{diagram = qi["Diagram"], unlabelled = qi["Diagram", "ModeLabels" -> False]},
            Head[diagram] === Graphics && Count[diagram, _Line, Infinity] >= 4 && Count[diagram, _Rectangle, Infinity] === 1 &&
                Cases[diagram, Text[Style[p_Integer, 10, ___], _, _] :> p, Infinity] === Range[3] &&
                Cases[unlabelled, Text[Style[p_Integer, 10, ___], _, _] :> p, Infinity] === {} &&
                ! FreeQ[diagram, Line[{{{.75, -1}, {1.25, -3}}, {{.75, -3}, {1.25, -1}}}]]
        ]
    ],
    True,
    TestID -> "QI-diagram"
]

(* the coupler of a mesh cell shows its beam splitter's own (theta, phi); the cell phase is a separate box *)
VerificationTest[
    SeedRandom[19];
    With[{qi = QuantumInterferometer["Random"[4]]},
        With[{diagram = qi["Diagram", "ShowParameters" -> True], hidden = qi["Diagram", "ShowParameters" -> False]},
            {
                Count[diagram, _Rectangle, Infinity] == qi["PhaseShifterCount"],
                Cases[diagram, Tooltip[_, Row[{"(\[Theta], \[Phi]) = ", tf_, __}]] :> tf, Infinity] == Cases[qi["Elements"], ("BS"[t_, f_] -> _) :> {t, f}],
                Count[diagram, Text[Style[_Row | _NumberForm, 8], _], Infinity] === qi["BeamSplitterCount"] + qi["PhaseShifterCount"],
                Count[hidden, Text[Style[_Row | _NumberForm, 8], _], Infinity] === 0
            }
        ]
    ],
    {True, True, True, True},
    TestID -> "QI-diagram-labels-beam-splitter-parameters"
]

VerificationTest[
    With[{plot = QuantumInterferometer[{{3, -4}, {4, 3}}/5]["MatrixPlot", ImageSize -> 111]},
        With[{rasters = Cases[plot, Raster[data_, ___] :> data, Infinity],
              intensity = Map[{1 - .95 #, 1 - .75 #, 1 - .4 #} &, {{16/25, 9/25}, {9/25, 16/25}}, {2}],
              phase = {{{.95, .285, .285}, {.95, .285, .285}}, {{.95, .285, .285}, {.285, .95, .95}}}},
            Length[rasters] === 2 && Max[Abs[rasters - {intensity, phase}]] < 10^-12 &&
                Cases[plot, g_Graphics :> (ImageSize /. Options[g, ImageSize]), Infinity] === {111, 111}
        ]
    ],
    True,
    TestID -> "QI-matrix-plot"
]

EndTestSection[]


BeginTestSection["QuantumInterferometer: native pure loss"]

(* Native vacuum loss: independent count laws, preserved coherences and unitary
   dilations test a channel, rather than treating attenuation as a state vector. *)
VerificationTest[
    With[{a = QuantumInterferometer[{"Loss"[2/3]}],
          b = QuantumInterferometer[{"Loss"[2/3] -> {1}}],
          c = QuantumInterferometer[{"Loss"[1] -> 5}]},
        a === b && a["Modes"] === 1 && a["TransferMatrix"] === {{Sqrt[2/3]}} &&
            a["LossCount"] === 1 && c["Modes"] === 5 && c["Unitary"] === IdentityMatrix[5]
    ],
    True,
    TestID -> "QI-loss-normalization-default-mode-and-inference-before-pruning"
]

VerificationTest[
    QuantumInterferometer[{"Loss"[#]}] & /@ {-1/3, 4/3},
    {$Failed, $Failed}, {QuantumInterferometer::loss, QuantumInterferometer::loss},
    TestID -> "QI-loss-rejects-out-of-range-transmissivities"
]

VerificationTest[
    QuantumInterferometer[{"Loss"[#]}] & /@ {I, 1 + I},
    {$Failed, $Failed}, {QuantumInterferometer::loss, QuantumInterferometer::loss},
    TestID -> "QI-loss-rejects-nonreal-transmissivities"
]

VerificationTest[
    {QuantumInterferometer[{"Loss"[]}], QuantumInterferometer[{"Loss"[1/2, 0]}]},
    {$Failed, $Failed}, {QuantumInterferometer::elem, QuantumInterferometer::elem},
    TestID -> "QI-loss-preserves-element-arity-validation"
]

VerificationTest[
    {QuantumInterferometer[{"Loss"[1/2] -> 0}], QuantumInterferometer[{"Loss"[1/2] -> 3}, 2]},
    {$Failed, $Failed}, {QuantumInterferometer::elem, QuantumInterferometer::modes},
    TestID -> "QI-loss-preserves-mode-validation"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"Loss"[2/5]}]},
        Table[qi["Probability", {3} -> {k}], {k, 0, 4}] === {27/125, 54/125, 36/125, 8/125, 0} &&
            Total[qi["Probabilities", {3}]] === 1 &&
            qi["Probability", {0} -> {0}] === 1 && qi["Probability", {0} -> {1}] === 0
    ],
    True,
    TestID -> "QI-loss-three-photon-binomial-law-vacuum-and-no-photon-creation"
]

VerificationTest[
    TimeConstrained[
        QuantumInterferometer[{"Loss"[3/4]}]["Probability", {80} -> {79}] === 20 (3/4)^79,
        5, $TimedOut
    ],
    True,
    TestID -> "QI-loss-large-occupation-uses-direct-count-law-without-cutoff"
]

(* the count law is the bare binomial: normalized, with mean n eta, and no photon is created *)
VerificationTest[
    Module[{eta, n, k},
        With[{probability = QuantumInterferometer[{"Loss"[eta]}]["Probability", {n} -> {k}]},
            {
                probability === Binomial[n, k] eta^k (1 - eta)^(n - k),
                QuantumInterferometer[{"Loss"[eta]}]["Probability", {n} -> {k},
                    Assumptions -> Element[{n, k}, Integers] && 0 <= k <= n && 0 < eta < 1] === probability,
                Sum[probability, {k, 0, n}],
                Simplify[Sum[k probability, {k, 0, n}]] === n eta,
                Table[probability /. {n -> input, k -> output}, {input, 0, 2}, {output, input + 1, 3}],
                (probability /. {n -> 2, k -> 1, eta -> 1/2}) === 1/2
            }
        ]
    ],
    {True, True, 1, True, {{0, 0, 0}, {0, 0}, {0}}, True},
    TestID -> "QI-loss-symbolic-count-law-is-the-bare-binomial"
]

VerificationTest[
    Module[{eta, n, k},
        With[{probability = QuantumInterferometer[{"Loss"[eta]}]["Probability", {n} -> {k}]},
            FreeQ[probability, _QuantumInterferometer] &&
                FreeQ[probability, s_Symbol /; StringStartsQ[Context[s], "Wolfram`QuantumFramework`"]] &&
                FreeQ[probability, ConditionalExpression | Piecewise | If] &&
                FullSimplify[probability - Binomial[n, k] eta^k (1 - eta)^(n - k),
                    Element[{n, k}, Integers] && 0 <= k <= n && 0 < eta < 1] === 0
        ]
    ],
    True,
    TestID -> "QI-loss-symbolic-count-output-contains-only-native-expressions"
]

VerificationTest[
    Module[{eta, mu, n, r, k, l, q},
        With[{probability = QuantumInterferometer[{"Loss"[eta] -> 1, "PS"[Pi/7] -> 1,
                "Loss"[2/3] -> 1, "Loss"[mu] -> 2}, 3]["Probability", {n, r, 1} -> {k, l, q}]},
            FreeQ[probability, _QuantumInterferometer] &&
                FreeQ[probability, s_Symbol /; StringStartsQ[Context[s], "Wolfram`QuantumFramework`"]] &&
                FullSimplify[probability - Binomial[n, k] (2 eta/3)^k (1 - 2 eta/3)^(n - k) *
                    Binomial[r, l] mu^l (1 - mu)^(r - l),
                    Element[{n, r, k, l}, Integers] && 0 <= k <= n && 0 <= l <= r &&
                        q == 1 && 0 < eta < 1 && 0 < mu < 1] === 0 &&
                (probability /. {n -> 3, r -> 2, k -> 1, l -> 1, q -> 0, eta -> 3/5, mu -> 1/4}) === 0 &&
                (probability /. {n -> 3, r -> 2, k -> 1, l -> 1, q -> 1, eta -> 3/5, mu -> 1/4}) === 81/500
        ]
    ],
    True,
    TestID -> "QI-loss-symbolic-count-output-composes-independent-attenuators-and-untouched-modes"
]

(* the symbolic law reaches the endpoints by continuity; an endpoint rate known in the evaluator
   gives its limit directly, and an untouched mode conserves its photons *)
VerificationTest[
    Module[{eta, n, k, r, l},
        With[{probability = QuantumInterferometer[{"Loss"[eta]}]["Probability", {n} -> {k}]},
            {
                And @@ Flatten[Table[
                    Limit[probability /. {n -> input, k -> output}, eta -> transmissivity] ===
                        If[transmissivity == 0, KroneckerDelta[output, 0], KroneckerDelta[input, output]],
                    {transmissivity, {0, 1}}, {input, {0, 1, 3}}, {output, Range[0, 4]}]],
                QuantumInterferometer[{"Loss"[0]}]["Probability", {n} -> {k}] === KroneckerDelta[k, 0],
                QuantumInterferometer[{"Loss"[eta] -> 1}, 2]["Probability", {n, r} -> {k, l}] === probability KroneckerDelta[r, l],
                Table[QuantumInterferometer[{"Loss"[0]}]["Probability", {3} -> {output}], {output, 0, 4}]
            }
        ]
    ],
    {True, True, True, {1, 0, 0, 0, 0}},
    TestID -> "QI-loss-symbolic-count-endpoints-by-continuity-and-in-the-evaluator"
]

VerificationTest[
    Module[{eta, n, k},
        With[{symbolic = QuantumInterferometer[{"Loss"[eta]}]["Probability", {n} -> {k}]},
            And @@ Map[Function[transmissivity,
                With[{qi = QuantumInterferometer[{"Loss"[transmissivity]}]},
                    {probability = qi["Probability", {n} -> {k}]},
                    {lost = probability /. {n -> 1, k -> 0}, value = probability /. {n -> 3, k -> 1}},
                    lost === qi["Probability", {1} -> {0}] && lost === 1 - transmissivity && lost > 0 &&
                        (symbolic /. {n -> 1, k -> 0, eta -> transmissivity}) === lost &&
                        Precision[lost] === Precision[1 - transmissivity] &&
                        If[Precision[transmissivity] === Infinity,
                            value === 3 transmissivity (1 - transmissivity)^2 && Precision[value] === Infinity,
                            Abs[value/(3 transmissivity (1 - transmissivity)^2) - 1] < 10^-13 &&
                                Precision[value] >= Min[Precision[transmissivity], Precision[1 - transmissivity]] - 6
                        ]
                ]], {2/5, 1 - 10^-40, 1. - 2^-50, N[2/5, 60], N[1 - 10^-40, 80]}]
        ]
    ],
    True,
    TestID -> "QI-loss-symbolic-count-evaluation-preserves-exact-machine-and-arbitrary-precision"
]

(* inexact endpoint transmissivities: concrete counts keep the native binomial law and its
   precision, and the symbolic-count law built at that endpoint agrees with them *)
VerificationTest[
    Module[{n, k, parameter},
        And @@ Map[Function[transmissivity,
            With[{qi = QuantumInterferometer[{"Loss"[transmissivity]}]},
                {bound = qi["Probability", {n} -> {k}]},
                FreeQ[bound, Piecewise | If] &&
                And @@ Flatten[Table[
                    With[{reference = PDF[BinomialDistribution[input, parameter], output] /. parameter -> transmissivity},
                        qi["Probability", {input} -> {output}] === reference && (bound /. {n -> input, k -> output}) == reference
                    ],
                    {input, {0, 1, 3}}, {output, Range[0, 4]}]]
            ]], {0., N[0, 60]}]
    ],
    True,
    TestID -> "QI-loss-inexact-endpoint-counts-keep-native-precision"
]

VerificationTest[
    Module[{eta},
        With[{qi = QuantumInterferometer[{"Loss"[eta]}]},
            qi["Probability", {1} -> {0}] === 1 - eta &&
            qi["Probability", {1} -> {0}, Assumptions -> 0 <= eta <= 1] === 1 - eta &&
            Assuming[0 <= eta <= 1, qi["Probability", {1} -> {0}]] === 1 - eta
        ]
    ],
    True,
    TestID -> "QI-loss-transmissivity-domain-and-native-assumptions"
]

VerificationTest[
    QuantumInterferometer[{"Loss"[1/2]}]["Probability", {1} -> {0}, Assumptions -> False],
    $Failed,
    {Refine::fas},
    TestID -> "QI-loss-probability-rejects-false-assumptions"
]

VerificationTest[
    With[{a = QuantumInterferometer[{"Loss"[2/3]}], b = QuantumInterferometer[{"Loss"[3/5]}]},
        b[a]["TransferMatrix"] === {{Sqrt[2/5]}} &&
            Table[b[a]["Probability", {3} -> {k}], {k, 0, 3}] === {27/125, 54/125, 36/125, 8/125}
    ],
    True,
    TestID -> "QI-loss-composition-multiplies-transmissivities"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"Loss"[1/3]}], alpha = (1 + 2 I)/3},
        qi["Coherent"[{alpha}]] === "Coherent"[{alpha/Sqrt[3]}] &&
            FullSimplify[qi["Probability", "Coherent"[{alpha}] -> {4}] - Exp[-5/27] (5/27)^4/4!] === 0 &&
            qi["Amplitude", "Coherent"[{alpha}] -> "Coherent"[{alpha/Sqrt[3]}]] === 1 &&
            FullSimplify[Abs[qi["Amplitude", "Coherent"[{alpha}] -> {1}]]^2 - Exp[-5/27] 5/27] === 0
    ],
    True,
    TestID -> "QI-loss-complex-coherent-input-stays-pure-with-Poisson-counts"
]

VerificationTest[
    Module[{eta},
        FullSimplify[Normal[QuantumInterferometer[{"Loss"[eta]}][QuantumState[{1, I}/Sqrt[2], 2]]["DensityMatrix"]] ==
            {{1 - eta/2, -I Sqrt[eta]/2}, {I Sqrt[eta]/2, eta/2}}, 0 < eta < 1]
    ],
    True,
    TestID -> "QI-loss-symbolic-superposition-coherences-scale-by-square-root"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"Loss"[2/5]}], rho = {{3/5, I/10}, {-I/10, 2/5}}},
        With[{out = Normal[qi[QuantumState[rho, 2]]["DensityMatrix"]]},
            out === {{21/25, I Sqrt[2/5]/10}, {-I Sqrt[2/5]/10, 4/25}} &&
                Tr[out] === 1 && out === ConjugateTranspose[out] && Min[Eigenvalues[out]] > 0
        ]
    ],
    True,
    TestID -> "QI-loss-mixed-state-trace-Hermiticity-positivity-and-coherences"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"Loss"[1/3] -> 2}, 2],
          input = QuantumState[{0, I, 1, 0}/Sqrt[2], {2, 2}],
          survivor = {0, I/Sqrt[3], 1, 0}/Sqrt[2]},
        FullSimplify[Normal[qi[input]["DensityMatrix"]] -
            (Outer[Times, survivor, Conjugate[survivor]] + DiagonalMatrix[{1/3, 0, 0, 0}])] === ConstantArray[0, {4, 4}]
    ],
    True,
    TestID -> "QI-loss-target-mode-and-entangled-spectator-coherences"
]

VerificationTest[
    Module[{a, b},
        With[{input = QuantumState[{a, b}, 2, "Parameters" -> {a, b}], qi = QuantumInterferometer[{"Loss"[1/2]}]},
            With[{out = qi[input]},
                out["ParameterSpec"] === input["ParameterSpec"] &&
                    FullSimplify[Normal[out[<|a -> 1/Sqrt[2], b -> I/Sqrt[2]|>]["DensityMatrix"]] -
                        {{3/4, -I/(2 Sqrt[2])}, {I/(2 Sqrt[2]), 1/4}}] === ConstantArray[0, {2, 2}]
            ]
        ]
    ],
    True,
    TestID -> "QI-loss-state-parameters-survive-and-bind-after-propagation"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"Loss"[2/5]}], input = QuantumState[{1, I, -1}/Sqrt[3], 3]},
        FullSimplify[Normal[qi["Channel", 3][input]["DensityMatrix"]] - Normal[qi[input]["DensityMatrix"]]] === ConstantArray[0, {3, 3}]
    ],
    True,
    TestID -> "QI-loss-native-QuantumChannel-interop-preserves-bosonic-dimension"
]

VerificationTest[
    With[{before = QuantumInterferometer[{"Loss"[1/2] -> 1, "BS"[Pi/4]}],
          after = QuantumInterferometer[{"BS"[Pi/4], "Loss"[1/2] -> 1}]},
        {before["Probability", {1, 0} -> {0, 0}], after["Probability", {1, 0} -> {0, 0}],
         before["Probability", {1, 0} -> {0, 1}], after["Probability", {1, 0} -> {0, 1}]} === {1/2, 1/4, 1/4, 1/2}
    ],
    True,
    TestID -> "QI-loss-cannot-commute-mode-dependent-attenuation-through-mixing"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[Pi/6, Pi/2], "Loss"[2/3] -> 2,
            "BS"[Pi/4, Pi/2], "Loss"[3/5] -> 1}],
          dilation = QuantumInterferometer[{"BS"[Pi/6, Pi/2] -> {1, 2},
            "BS"[ArcCos[Sqrt[2/3]], 0] -> {2, 3}, "BS"[Pi/4, Pi/2] -> {1, 2},
            "BS"[ArcCos[Sqrt[3/5]], 0] -> {1, 4}}, 4]},
        And @@ Map[Function[in,
            With[{reference = Merge[KeyValueMap[<|Take[#1["Name"], 2] -> #2|> &,
                    dilation["Probabilities", PadRight[in, 4]]], Total],
                  direct = qi["Probabilities", in]},
                FullSimplify[Total[direct]] === 1 &&
                    AllTrue[KeyValueMap[FullSimplify[#2 - Lookup[reference, Key[#1["Name"]], 0]] === 0 &, direct], TrueQ] &&
                    AllTrue[KeyValueMap[FullSimplify[qi["Probability", in -> #1["Name"]] - #2] === 0 &, direct], TrueQ]
            ]], {{2, 0}, {1, 2}}]
    ],
    True,
    TestID -> "QI-loss-colliding-multimode-probabilities-match-independent-vacuum-dilation"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[Pi/6, Pi/2], "Loss"[2/3] -> 2, "BS"[Pi/4, 0], "Loss"[3/5] -> 1}],
          input = FockState[{1, 1}, 3]},
        With[{rho = Normal[qi[input]["DensityMatrix"]], channelRho = Normal[qi["Channel", 3][input]["DensityMatrix"]]},
            FullSimplify[rho - channelRho] === ConstantArray[0, {9, 9}] &&
                AllTrue[KeyValueMap[FullSimplify[rho[[fockIndex[#1["Name"], 3], fockIndex[#1["Name"], 3]]] - #2] === 0 &,
                    qi["Probabilities", {1, 1}]], TrueQ]
        ]
    ],
    True,
    TestID -> "QI-loss-interleaved-state-channel-and-count-routes-agree"
]

VerificationTest[
    With[{zero = QuantumInterferometer[{"Loss"[0]}], one = QuantumInterferometer[{"Loss"[1]}, "DropIdentities" -> False]},
        zero["Probability", {3} -> {0}] === 1 && zero["Probability", {3} -> {1}] === 0 &&
            Normal[zero[QuantumState[{1, I, 1}/Sqrt[3], 3]]["DensityMatrix"]] === DiagonalMatrix[{1, 0, 0}] &&
            one["Amplitude", {3} -> {3}] === 1 && one["Dagger"]["Unitary"] === {{1}} &&
            zero["Probability", {0} -> {0}] === 1
    ],
    True,
    TestID -> "QI-loss-endpoints-reset-or-identity-without-zero-powers"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"Loss"[1/2]}]}, {qi["Amplitude", {1} -> {0}], qi[{1}]}],
    {$Failed, $Failed}, {QuantumInterferometer::lossy, QuantumInterferometer::lossy},
    TestID -> "QI-loss-does-not-invent-unconditional-Fock-amplitudes"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"Loss"[1/2]}]}, {qi["Unitary"], qi["Dagger"]}],
    {$Failed, $Failed}, {QuantumInterferometer::lossy, QuantumInterferometer::lossy},
    TestID -> "QI-loss-does-not-invent-unitarity-or-an-inverse"
]

VerificationTest[
    QuantumInterferometer[QuantumInterferometer[{"Loss"[1/2]}]],
    $Failed, {QuantumInterferometer::lossy},
    TestID -> "QI-loss-does-not-pretend-to-support-unitary-remeshing"
]

VerificationTest[
    {QuantumInterferometer[{"S"[1/5], "Loss"[1/2]}], QuantumInterferometer[{"Loss"[1/2], "S"[1/5]}]},
    {$Failed, $Failed}, {QuantumInterferometer::lossactive, QuantumInterferometer::lossactive},
    TestID -> "QI-loss-with-squeezing-is-rejected-in-both-element-orders"
]

VerificationTest[
    {QuantumInterferometer[{"S"[1/5]}][QuantumInterferometer[{"Loss"[1/2]}]],
     QuantumInterferometer[{"Loss"[1/2]}][QuantumInterferometer[{"S"[1/5]}]]},
    {$Failed, $Failed}, {QuantumInterferometer::lossactive, QuantumInterferometer::lossactive},
    TestID -> "QI-loss-with-squeezing-is-rejected-in-both-composition-orders"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"Loss"[1/2]}]},
        {qi["Probability", {-1} -> {0}], qi["Probability", {1/2} -> {0}]}],
    {$Failed, $Failed}, {QuantumInterferometer::occ, QuantumInterferometer::occ},
    TestID -> "QI-loss-preserves-nonnegative-integer-occupation-guards"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"Loss"[1/2]}]},
        {qi["Probability", {1} -> {0, 0}], qi["Probabilities", {1}, "Cutoff" -> -1],
         qi["Probability", {1, "Coherent"[1]} -> {0}]}],
    {$Failed, $Failed, $Failed}, {QuantumInterferometer::occ, QuantumInterferometer::cutoff, QuantumInterferometer::occ},
    TestID -> "QI-loss-preserves-occupation-length-and-cutoff-guards"
]

VerificationTest[
    QuantumInterferometer[{"BS"[Pi/4], "Loss"[1/2] -> 1}]["Probability", {1, "Coherent"[1]} -> {1, 0}],
    $Failed,
    {QuantumInterferometer::losscoherent},
    TestID -> "QI-loss-does-not-treat-mixed-Fock-coherent-input-as-pure"
]

VerificationTest[
    Module[{eta, theta},
        With[{qi = QuantumInterferometer[{"BS"[theta], "Loss"[eta] -> 2},
                "Parameters" -> {theta, eta}, "DropIdentities" -> False]},
            qi[<|theta -> 0, eta -> 1|>]["Unitary"] === IdentityMatrix[2] &&
                qi[Pi/4, 1/2]["Probability", {1, 0} -> {0, 0}] === 1/4 &&
                qi[Pi/4, 1/2]["Parameters"] === {}
        ]
    ],
    True,
    TestID -> "QI-loss-parameter-binding-revalidates-and-restores-lossless-limit"
]

VerificationTest[
    With[{delta = 10^-40},
        QuantumInterferometer[{"Loss"[1 - delta]}]["Probability", {1} -> {0}] === delta &&
            QuantumInterferometer[{"BS"[Pi/4], "Loss"[1 - delta] -> 1}]["Probability", {1, 0} -> {0, 0}] === delta/2 &&
            Precision[QuantumInterferometer[{"Loss"[1 - delta]}]["TransferMatrix"]] === Infinity
    ],
    True,
    TestID -> "QI-loss-exact-weak-attenuation-is-never-pruned-by-tolerance"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[.31, .27], "Loss"[.8] -> 1,
            "BS"[.49, -.38], "Loss"[.6] -> 2}]},
        Max[Abs[ConjugateTranspose[qi["TransferMatrix"]] . qi["TransferMatrix"] +
            First[qi]["LossGram"] - IdentityMatrix[2]]] < 10^-14 &&
            Abs[Total[qi["Probabilities", {2, 1}]] - 1] < 10^-13 &&
            Min[Values[qi["Probabilities", {2, 1}]]] >= -10^-14
    ],
    True,
    TestID -> "QI-loss-compact-Gram-completeness-and-machine-probability-normalization"
]

VerificationTest[
    And @@ Table[
        With[{eta = N[2/5, p], theta = N[Pi/6, p]},
            With[{qi = QuantumInterferometer[{"BS"[theta, 0], "Loss"[eta] -> 1}],
                  reference = N[3 (2/5)^2 (3/4)^2 (1/4), p]},
                With[{prob = qi["Probability", {3, 0} -> {2, 1}]},
                    Abs[prob - reference] < If[p === MachinePrecision, 10^-13, 10^(-p + 8)] &&
                        If[p === MachinePrecision, MachineNumberQ[prob], Precision[prob] >= p - 6]
                ]
            ]
        ], {p, {MachinePrecision, 60, 80}}],
    True,
    TestID -> "QI-loss-machine-and-arbitrary-precision-agree-with-analytic-nonzero-reference"
]

VerificationTest[
    With[{eta = N[1 - 10^-40, 80]},
        With[{prob = QuantumInterferometer[{"BS"[N[Pi/4, 80]], "Loss"[eta] -> 1}]["Probability", {1, 0} -> {0, 0}]},
            prob > 0 && Abs[prob/(10^-40/2) - 1] < 10^-35 && Accuracy[prob] > 70
        ]
    ],
    True,
    TestID -> "QI-loss-Gram-accumulation-preserves-weak-arbitrary-precision-loss"
]

VerificationTest[
    Module[{theta},
        With[{qi = QuantumInterferometer[{"BS"[theta[1], 0], "Loss"[.5] -> 1}]},
            With[{p = qi["Probability", {1, 0} -> {0, 0}]},
                FreeQ[p, theta[1.]] && ! FreeQ[p, theta[1]] &&
                    Abs[(p /. theta[1] -> Pi/6) - 3/8] < 10^-14
            ]
        ]
    ],
    True,
    TestID -> "QI-loss-numerical-precision-work-never-rewrites-indexed-symbolic-parameters"
]

VerificationTest[
    QuantumInterferometer[{"Loss"[#]}] & /@ {{1/2}, "bad"},
    {$Failed, $Failed}, {QuantumInterferometer::loss, QuantumInterferometer::loss},
    TestID -> "QI-loss-rejects-nonscalar-and-nonnumeric-parameter-data"
]

VerificationTest[
    QuantumInterferometer[{"Loss"[#]}] & /@ {Infinity, Indeterminate},
    {$Failed, $Failed}, {QuantumInterferometer::loss, QuantumInterferometer::loss},
    TestID -> "QI-loss-rejects-infinite-and-indeterminate-transmissivity"
]

VerificationTest[
    With[{eta = 1. - 2^-50},
        With[{qi = QuantumInterferometer[{"Loss"[eta]}]},
            qi["LossCount"] === 1 && qi["Probability", {1} -> {0}] > 0 &&
                qi["Probability", {1} -> {0}] === 1. - eta
        ]
    ],
    True,
    TestID -> "QI-loss-representable-machine-weak-loss-is-not-an-approximate-identity"
]

VerificationTest[
    Module[{eta, theta},
        With[{reference = Block[{$Assumptions = True}, QuantumInterferometer[{"BS"[theta], "Loss"[eta] -> 1}]]},
            Block[{$Assumptions = eta < 0 && theta == 0},
                QuantumInterferometer[{"BS"[theta], "Loss"[eta] -> 1}] === reference &&
                    FullSimplify[reference["Probability", {1, 0} -> {0, 0},
                        Assumptions -> 0 < eta < 1 && Element[theta, Reals]] - (1 - eta) Cos[theta]^2,
                        Assumptions -> 0 < eta < 1 && Element[theta, Reals]] === 0
            ]
        ]
    ],
    True,
    TestID -> "QI-loss-construction-and-explicit-query-domain-ignore-conflicting-ambient-assumptions"
]

VerificationTest[
    QuantumInterferometer[{"Loss"[#]}] & /@ {SparseArray[{1/2, 0}], SparseArray[{{1/2}}]},
    {$Failed, $Failed}, {QuantumInterferometer::loss, QuantumInterferometer::loss},
    TestID -> "QI-loss-rejects-sparse-vector-and-matrix-transmissivities"
]

VerificationTest[
    QuantumInterferometer[{"Loss"[#]}] & /@ {NumericArray[{.5}, "Real64"], Quantity[1/2, "Meters"]},
    {$Failed, $Failed}, {QuantumInterferometer::loss, QuantumInterferometer::loss},
    TestID -> "QI-loss-rejects-numeric-arrays-and-dimensional-quantities"
]

VerificationTest[
    Module[{eta},
        With[{qi = QuantumInterferometer[{"Loss"[eta]}]},
            With[{probabilities = qi["Probabilities", "Coherent"[{1}], "Cutoff" -> 2]},
                And @@ Table[Lookup[probabilities, QuditName[{k}]] ===
                    qi["Probability", "Coherent"[{1}] -> {k}], {k, 0, 2}] &&
                    Values[probabilities] === {E^(-eta), eta E^(-eta), eta^2 E^(-eta)/2} &&
                    Values[probabilities /. eta -> 0] === {1, 0, 0}
            ]
        ]
    ],
    True,
    TestID -> "QI-loss-coherent-probability-table-is-the-bare-Poisson-law"
]

VerificationTest[
    Module[{eta},
        FullSimplify[QuantumInterferometer[{"Loss"[eta[1]]}]["Probability", {2} -> {1}] -
            2 eta[1] (1 - eta[1]), Assumptions -> 0 < eta[1] < 1] === 0 &&
            FullSimplify[QuantumInterferometer[{"Loss"[Sin[eta]]}]["Probability", {2} -> {1}] -
                2 Sin[eta] (1 - Sin[eta]), Assumptions -> 0 < eta < Pi/2] === 0
    ],
    True,
    TestID -> "QI-loss-scalar-validation-preserves-indexed-and-native-symbolic-expressions"
]

VerificationTest[
    QuantumInterferometer[{"Loss"[#]}] & /@ {VectorSymbol["v", 2], MatrixSymbol["M", {2, 3}]},
    {$Failed, $Failed}, {QuantumInterferometer::loss, QuantumInterferometer::loss},
    TestID -> "QI-loss-rejects-native-symbolic-vector-and-matrix-transmissivities"
]

VerificationTest[
    QuantumInterferometer[{"Loss"[#]}] & /@ {
        ArraySymbol["A", {2, 3, 4}], Inactive[TensorProduct][{{1/2, 0}, {0, 1/2}}, {{1/3, 0}, {0, 1/3}}]},
    {$Failed, $Failed}, {QuantumInterferometer::loss, QuantumInterferometer::loss},
    TestID -> "QI-loss-rejects-symbolic-tensors-and-array-valued-structural-expressions"
]

VerificationTest[
    QuantumInterferometer[{"Loss"[#]}] & /@ {
        ByteArray[{1, 2}], Wolfram`Arrays`ArrayObject[SparseArray[{1/2, 0}]]},
    {$Failed, $Failed}, {QuantumInterferometer::loss, QuantumInterferometer::loss},
    TestID -> "QI-loss-rejects-wrapped-array-transmissivities"
]

VerificationTest[
    Module[{eta},
        With[{reference = QuantumInterferometer[{"Loss"[eta]}]},
            Block[{$Assumptions = Element[eta, Matrices[{2, 2}]]},
                QuantumInterferometer[{"Loss"[eta]}] === reference &&
                    reference["Probability", {1} -> {0}, Assumptions -> 0 < eta < 1] === 1 - eta
            ]
        ]
    ],
    True,
    TestID -> "QI-loss-scalar-construction-is-independent-of-ambient-array-assumptions"
]

VerificationTest[
    Module[{eta},
        With[{qi = QuantumInterferometer[{"Loss"[eta]}]},
            {qi["Probability", {1} -> {0}, Assumptions -> eta == I],
             qi["Probability", {1} -> {0}, Assumptions -> eta > 1],
             Values[qi["Probabilities", {1}]] === {1 - eta, eta}}
        ]
    ],
    {$Failed, $Failed, True}, {Refine::fas, Refine::fas},
    TestID -> "QI-loss-assumptions-outside-the-physical-domain-fail"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"Loss"[2/5], "BS"[Pi/6]}]},
        ContainsAll[qi["Properties"], {"Channel", "LossCount", "BlochMessiah", "GBSForm"}] &&
            DuplicateFreeQ[qi["Properties"]]
    ],
    True,
    TestID -> "QI-loss-integration-preserves-upstream-and-channel-property-discovery"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[Pi/6], "Loss"[2/5] -> 2}]},
        {qi["BlochMessiah"], qi["GBSForm"]}
    ],
    {$Failed, $Failed}, {QuantumInterferometer::lossy, QuantumInterferometer::lossy},
    TestID -> "QI-loss-unitary-normal-forms-reject-numeric-attenuation"
]

VerificationTest[
    Module[{eta},
        With[{qi = QuantumInterferometer[{"Loss"[eta], "BS"[Pi/6]},
                "Parameters" -> {eta}, "DropIdentities" -> False]},
            {qi["BlochMessiah"], qi["GBSForm"],
             With[{lossless = qi[1]},
                lossless["BlochMessiah"] === lossless && lossless["GBSForm"] === lossless &&
                    lossless["LossCount"] === 1 && Precision[lossless["Unitary"]] === Infinity &&
                    lossless["Unitary"] === QuantumInterferometer[{"BS"[Pi/6]}]["Unitary"]
             ]}
        ]
    ],
    {$Failed, $Failed, True}, {QuantumInterferometer::lossy, QuantumInterferometer::lossy},
    TestID -> "QI-loss-unitary-normal-forms-guard-symbolic-loss-and-recover-identity-limit"
]

EndTestSection[]
