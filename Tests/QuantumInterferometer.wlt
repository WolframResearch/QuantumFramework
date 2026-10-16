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
            qi["Unitary"] == transferMatrix[qi["CircuitOperator", 2], 2, 2],
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

(* The nulling sweep owns its working matrix and mutates it in place; it must leave
   the caller's input matrix untouched (a regression that dropped the copy-on-first-
   write and aliased the input would fail this). *)
VerificationTest[
    SeedRandom[101];
    With[{w = RandomVariate[CircularUnitaryMatrixDistribution[12]]},
        With[{before = w + 0},
            QuantumInterferometer[w, Method -> "Reck"];
            QuantumInterferometer[w, Method -> "Clements"];
            w === before
        ]
    ],
    True,
    TestID -> "QI-decomposition-does-not-mutate-input"
]

(* The Reck and Clements meshes rebuild Haar unitaries and stay unitary at larger m,
   driving hundreds of in-place elimination steps past the sizes above. *)
VerificationTest[
    SeedRandom[102];
    Table[
        With[{w = RandomVariate[CircularUnitaryMatrixDistribution[m]]},
            With[{u = rebuild[QuantumInterferometer[w, Method -> method]]},
                Max[Abs[u - w]] < 10^-10 && Max[Abs[u . ConjugateTranspose[u] - IdentityMatrix[m]]] < 10^-10
            ]
        ],
        {m, {16, 24, 32}}, {method, $exactMethods}
    ],
    ConstantArray[True, {3, 2}],
    TestID -> "QI-rebuilds-and-stays-unitary-at-larger-m"
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

(* Fix 2 (closed-form Givens cell + RootReduce on the exact path): the exact decomposition
   of FourierMatrix at a cyclotomic degree that previously did not terminate now completes
   and reconstructs the matrix exactly, with exact (Precision Infinity) parameters. The
   TimeConstrained is a regression guard against the trig-of-arctan swell returning. *)
VerificationTest[
    TimeConstrained[
        Table[
            With[{qi = QuantumInterferometer[FourierMatrix[m], Method -> method]},
                Precision[qi["Elements"]] === Infinity &&
                    Max[Abs[QuantumInterferometer[
                            qi["Elements"] /. (h_[ps__] -> md_) :> (h[Sequence @@ N[{ps}]] -> md), m
                        ]["Unitary"] - N[FourierMatrix[m]]]] < 10^-10
            ],
            {m, {5, 6}}, {method, $exactMethods}
        ],
        120, $TimedOut
    ],
    ConstantArray[True, {2, 2}],
    TestID -> "QI-exact-high-degree-decomposition-terminates"
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

(* Fix 5 (in-place single-owner phase-absorption and layout folds): placing a fixed list
   of beam splitters is deterministic and each one lands one column past the last column
   used on the modes it spans. A mis-threaded in-place position vector would shift columns. *)
VerificationTest[
    QuantumInterferometer[{"BS"[0.5, 0.] -> {1, 2}, "BS"[0.5, 0.] -> {2, 3}, "BS"[0.5, 0.] -> {1, 2}}, 3]["Layers"],
    {1, 2, 3},
    TestID -> "QI-layout-columns-are-deterministic"
]

(* Fix 5: the in-place mutable vectors in meshItems and layoutItems are owned per call, so
   decomposing the same matrix twice gives identical Elements and Layers; a leaked or aliased
   accumulator would make a second call differ from the first. *)
VerificationTest[
    SeedRandom[404];
    With[{w = RandomVariate[CircularUnitaryMatrixDistribution[8]]},
        Table[
            With[{a = QuantumInterferometer[w, Method -> method], b = QuantumInterferometer[w, Method -> method]},
                a["Elements"] === b["Elements"] && a["Layers"] === b["Layers"]
            ],
            {method, $methods}
        ]
    ],
    {True, True, True},
    TestID -> "QI-decomposition-is-repeatable"
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
            Max[Abs[Lookup[amplitudes, QuditName[#], 0] - out[[fockIndex[#, 4]]]] & /@ Tuples[Range[0, 3], 3]] < 10^-10
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

(* One photon per port of the Fourier tritter: the six patterns whose mode labels
   do not sum to a multiple of 3 are dark. *)
VerificationTest[
    Keys @ Select[QuantumInterferometer[FourierMatrix[3]]["Probabilities", {1, 1, 1}], PossibleZeroQ],
    QuditName /@ {{2, 1, 0}, {2, 0, 1}, {1, 2, 0}, {1, 0, 2}, {0, 2, 1}, {0, 1, 2}},
    TestID -> "QI-Fourier-tritter-suppression-law"
]

VerificationTest[
    SeedRandom[15];
    With[{qi = QuantumInterferometer["Random"[5]]},
        {
            Chop[Total[qi["Probabilities", {1, 1, 0, 1, 0}]] - 1],
            Chop[qi["Probability", {1, 1, 0, 1, 0} -> {0, 0, 2, 0, 1}] - qi["Probabilities", {1, 1, 0, 1, 0}][QuditName[{0, 0, 2, 0, 1}]]]
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

(* ---- symbolic photon numbers: nested sums of Wigner d-functions, evaluated back to the permanents ---- *)

(* substitute a concrete pattern into the formula for {n1, ..} -> {t1, ..} and Activate the sums *)
symbolicValue[f_, n_, t_] := Activate[f /. Thread[Array[Symbol["n" <> ToString[#]] &, Length[n]] -> n] /. Thread[Array[Symbol["t" <> ToString[#]] &, Length[t]] -> t]];
patterns[m_, nT_] := Select[Tuples[Range[0, nT], m], Total[#] == nT &];

VerificationTest[
    SeedRandom[41];
    With[{
        bs = QuantumInterferometer[{"BS"[Pi / 6, Pi / 3] -> {1, 2}}, 2],
        q3 = QuantumInterferometer["Random"[3]],
        q4 = QuantumInterferometer["Random"[4]]
    },
        With[{
            f2 = bs["Amplitude", {n1, n2} -> {t1, t2}],
            f3 = q3["Amplitude", {n1, n2, n3} -> {t1, t2, t3}],
            f4 = q4["Amplitude", {n1, n2, n3, n4} -> {t1, t2, t3, t4}]
        },
            {
                Union[FullSimplify[symbolicValue[f2, {2, 5}, #] - bs["Amplitude", {2, 5} -> #]] & /@ patterns[2, 7]],
                Max[Abs[N[symbolicValue[f3, {1, 2, 3}, #]] - q3["Amplitude", {1, 2, 3} -> #]] & /@ patterns[3, 6]] < 10^-12,
                Max[Abs[N[symbolicValue[f4, {1, 2, 0, 1}, #]] - q4["Amplitude", {1, 2, 0, 1} -> #]] & /@ patterns[4, 4]] < 10^-12,
                symbolicValue[f4, {1, 2, 0, 1}, {1, 1, 1, 2}]
            }
        ]
    ],
    {{0}, True, True, 0},
    TestID -> "QI-symbolic-photon-numbers-match-permanents"
]

(* a full m-mode mesh leaves (m - 1)(m - 2)/2 free indices; two modes give a single d-function *)
VerificationTest[
    SeedRandom[42];
    With[{bs = QuantumInterferometer[{"BS"[\[Theta], \[Phi]] -> {1, 2}}, 2], q4 = QuantumInterferometer["Random"[4]]},
        {
            FreeQ[bs["Amplitude", {nx, 0} -> {nx / 2, nx / 2}], Inactive[Sum]] && ! FreeQ[bs["Amplitude", {nx, 0} -> {nx / 2, nx / 2}], nx],
            Length[Cases[q4["Amplitude", {n1, n2, n3, n4} -> {t1, t2, t3, t4}], Inactive[Sum][_, iters__] :> {iters}, Infinity][[1]]],
            Chop[N[symbolicValue[q4["Probability", {n1, n2, n3, n4} -> {t1, t2, t3, t4}], {2, 1, 0, 1}, {0, 1, 2, 1}] - q4["Probability", {2, 1, 0, 1} -> {0, 1, 2, 1}]]]
        }
    ],
    {True, 3, 0},
    TestID -> "QI-symbolic-photon-numbers-structure"
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

(* ---- sampling: Clifford-Clifford draws agree with the exact distribution ---- *)

(* total variation distance of 5000 samples from "Probabilities", for one photon per mode and for a
   doubly occupied input; sampling noise alone is about 0.025 here *)
VerificationTest[
    SeedRandom[31];
    With[{qi = QuantumInterferometer["Random"[4]], k = 5000},
        Map[
            Total[Abs[Values[qi["Probabilities", #]] - Lookup[qi["SampleCounts", #, k], Keys[qi["Probabilities", #]], 0] / k]] / 2 < 0.06 &,
            {{1, 1, 1, 0}, {2, 1, 0, 0}}
        ]
    ],
    {True, True},
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
    {QuantumInterferometer::squeezed, QuantumInterferometer::sample, QuantumInterferometer::occ},
    TestID -> "QI-samples-rejects-squeezers-symbolic-and-bad-patterns"
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
    {QuantumInterferometer::truncation},
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


(* ---- squeezers: <t| U |s> = <0|U|0> Haf(K_{t+s}) / Sqrt[t! s!] ---- *)

hafnian = Wolfram`QuantumFramework`SecondQuantization`Interferometer`PackagePrivate`hafnian;

(* all-ones matrices count perfect matchings, (n - 1)!!; Haf([[0, G], [G^T, 0]]) = Perm(G) on both the pairing and the power-trace paths *)
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

fockAmplitudes[circuit_, in_, d_, outs_] := Normal[circuit[FockState[in, d]]["StateVector"]][[FromDigits[#, d] + 1 & /@ outs]]

(* Gaussian boson sampling: the mesh conserves photon number, so the Fock circuit is exact below d photons *)
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
            Max @ Abs @ KeyValueMap[#2 - If[Equal @@ #1["Name"], Tanh[0.5] ^ (2 #1["Name"][[1]]) / Cosh[0.5] ^ 2, 0] &, qi["Probabilities", {0, 0}, "Cutoff" -> 6]] < 10^-12,
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
            Total[qi["Probabilities", {1, 0}]] > 1 - 10^-4
        }
    ],
    {True, True},
    TestID -> "QI-squeezers-dagger-and-auto-cutoff"
]

(* passive networks ignore "Cutoff": photon number is fixed *)
VerificationTest[
    SeedRandom[23];
    With[{qi = QuantumInterferometer["Random"[3]]},
        qi["Probabilities", {1, 1, 0}, "Cutoff" -> 1] === qi["Probabilities", {1, 1, 0}] && qi["SqueezerCount"] == 0
    ],
    True,
    TestID -> "QI-passive-ignores-cutoff"
]

EndTestSection[]
