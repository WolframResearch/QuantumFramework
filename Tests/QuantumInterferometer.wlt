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
        With[{before = Hash[w, "SHA256"]},
            QuantumInterferometer[w, Method -> "Reck"];
            QuantumInterferometer[w, Method -> "Clements"];
            Hash[w, "SHA256"] === before
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
   and gives exact (Precision Infinity) parameters with a numerical reconstruction check. The
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
    FullSimplify[{
            QuantumInterferometer[{"PS"[\[Alpha]] -> 1, "PS"[\[Beta]] -> 2, "BS"[\[Theta], \[Gamma]] -> {1, 2}}, 2]["Unitary"] ==
                QuantumInterferometer[{"BS"[\[Theta], \[Gamma] + \[Alpha] - \[Beta]] -> {1, 2}, "PS"[\[Alpha]] -> 1, "PS"[\[Beta]] -> 2}, 2]["Unitary"],
            QuantumInterferometer[{"BS"[- \[Theta], \[Gamma]] -> {1, 2}}, 2]["Unitary"] ==
                QuantumInterferometer[{"BS"[\[Theta], \[Gamma] + Pi] -> {1, 2}}, 2]["Unitary"]
        }, Element[{\[Theta], \[Gamma], \[Alpha], \[Beta]}, Reals]],
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

(* Fix 3 (per-gate two-row update in elementsUnitary): building the transfer matrix from an
   element list reconstructs the unitary and returns a packed array. A regression dropping
   Developer`ToPackedArray would leave it unpacked; a wrong row convention would miss the unitary. *)
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

(* N leaves the mode numbers and layers exact, so the object stays valid *)
VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[Pi / 4, 0] -> {1, 2}, "PS"[Pi / 3] -> 2}, 2]},
        {N[qi] === qi, N[qi]["Modes"], N[qi]["Elements"]}
    ],
    {True, 2, {"BS"[Pi / 4, 0] -> {1, 2}, "PS"[Pi / 3] -> {2}}},
    TestID -> "QI-N-keeps-the-object-valid"
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

(* a mode no beam splitter touches keeps its photons: its equation survives as a KroneckerDelta *)
VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[Pi / 5, Pi / 7] -> {2, 3}, "PS"[Pi / 3] -> 1}, 3], single = QuantumInterferometer[{"PS"[Pi / 3] -> 1}, 1]},
        With[{f = qi["Amplitude", {n1, n2, n3} -> {t1, t2, t3}]},
            {
                Union[FullSimplify[symbolicValue[f, {1, 1, 0}, #] - qi["Amplitude", {1, 1, 0} -> #]] & /@ patterns[3, 2]],
                FullSimplify[single["Amplitude", {n1} -> {t1}] /. {n1 -> 2, t1 -> 2}] == single["Amplitude", {2} -> {2}]
            }
        ]
    ],
    {{0}, True},
    TestID -> "QI-symbolic-photon-numbers-untouched-modes"
]

(* Native WignerD exposes the general two-mode rotation; occupations remain discrete. *)
VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[\[Theta], \[Phi]] -> {1, 2}}, 2]},
        With[{amplitude = qi["Amplitude", {nx, ny} -> {tx, nx + ny - tx}]},
            FullSimplify[
                amplitude == WignerD[{(nx + ny) / 2, (nx - ny) / 2, tx - (nx + ny) / 2}, \[Phi], 2 \[Theta], - \[Phi]],
                Element[{nx, ny, tx}, Integers] && nx >= 0 && ny >= 0 && 0 <= tx <= nx + ny
            ]
        ]
    ],
    True,
    TestID -> "QI-symbolic-native-WignerD-general-two-mode-identity"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[Pi / 4, 0] -> {1, 2}}, 2]},
        With[{amplitude = qi["Amplitude", {nx, nx} -> {tx, 2 nx - tx}]},
            And @@ Flatten @ Table[
                FullSimplify[(amplitude /. {nx -> n, tx -> t}) - qi["Amplitude", {n, n} -> {t, 2 n - t}]] === 0,
                {n, 0, 5}, {t, 0, 2 n}
            ]
        ]
    ],
    True,
    TestID -> "QI-symbolic-native-WignerD-balanced-occupations"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[\[Theta], \[Phi]] -> {1, 2}}, 2]},
        With[{amplitude = qi["Amplitude", {nx, ny} -> {tx, nx + ny - tx}]},
            {
                amplitude /. {nx -> 0, ny -> 0, tx -> 0},
                (amplitude /. Thread[{nx, ny, tx} -> #]) & /@
                    {{0, 0, -1}, {0, 0, 1}, {-1, -1, 0}, {-1, 2, 0}, {2, -1, 0}, {1, 1, -1}, {1, 1, 3}},
                qi["Amplitude", {nx, ny} -> {tx, ty}] /. {nx -> 1, ny -> 1, tx -> 1, ty -> 0}
            }
        ]
    ],
    {1, ConstantArray[0, 7], 0},
    TestID -> "QI-symbolic-native-WignerD-range-vacuum-and-number-guards"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[\[Theta], \[Phi]] -> {1, 2}}, 2]},
        With[{amplitude = qi["Amplitude", {nx, ny} -> {tx, nx + ny - tx}]},
            With[{amplitudes = Table[amplitude /. {nx -> 2, ny -> 1, tx -> t}, {t, 0, 3}]},
                FullSimplify[
                    {
                        Total[Abs[amplitudes] ^ 2],
                        amplitudes /. \[Theta] -> 0,
                        amplitudes /. \[Theta] -> Pi / 2
                    },
                    Element[{\[Theta], \[Phi]}, Reals]
                ]
            ]
        ]
    ],
    {1, {0, 0, 1, 0}, {0, - Exp[I \[Phi]], 0, 0}},
    TestID -> "QI-symbolic-native-WignerD-normalization-and-splitter-limits"
]

VerificationTest[
    Map[
        Function[precision,
            With[{qi = QuantumInterferometer[{"BS"[N[3 / 10, precision], N[7 / 11, precision]] -> {1, 2}}, 2]},
                With[{amplitude = qi["Amplitude", {nx, ny} -> {tx, nx + ny - tx}]},
                    With[{value = amplitude /. {nx -> 2, ny -> 3, tx -> 1}, reference = qi["Amplitude", {2, 3} -> {1, 4}]},
                        If[precision === MachinePrecision,
                            Precision[value] === MachinePrecision && Abs[value - reference] < 10^-12,
                            Precision[value] > precision - 3 && Abs[value - reference] < 10 ^ (- precision + 5)
                        ]
                    ]
                ]
            ]
        ],
        {MachinePrecision, 60}
    ],
    {True, True},
    TestID -> "QI-symbolic-native-WignerD-preserves-numeric-precision"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{
        "BS"[Pi / 6, Pi / 7] -> {3, 1}, "BS"[Pi / 5, Pi / 11] -> {2, 3}, "BS"[Pi / 4, Pi / 13] -> {1, 2}
    }, 3]},
        With[{amplitude = qi["Amplitude", {nx, ny, nz} -> {tx, ty, tz}],
              vector = qi["CircuitOperator", 3][FockState[{1, 1, 0}, 3]]["StateVector"]},
            Max[Abs[
                N[Activate[amplitude /. Thread[{nx, ny, nz, tx, ty, tz} -> Join[{1, 1, 0}, #]]]] -
                N[vector[[fockIndex[#, 3]]]]
            ] & /@ patterns[3, 2]] < 10^-12
        ]
    ],
    True,
    TestID -> "QI-symbolic-native-WignerD-multimode-matches-Fock-circuit"
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
            With[{probs = qi["Probabilities", #], counts = qi["SampleCounts", #, k]},
                With[{keys = Union[Keys[probs], Keys[counts]]},
                    Total[counts] === k && Complement[Keys[counts], Keys[probs]] === {} &&
                        Total[Abs[Lookup[probs, keys, 0] - Lookup[counts, keys, 0] / k]] / 2 < 0.06
                ]
            ] &,
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
    {QuantumInterferometer::squeezed, QuantumInterferometer::notnumeric, QuantumInterferometer::occ},
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
    With[{qi = QuantumInterferometer[{"BS"[.4, .2] -> {1, 3}, "PS"[.3] -> 2}, 3]},
        With[{diagram = qi["Diagram"], unlabelled = qi["Diagram", "ModeLabels" -> False]},
            Head[diagram] === Graphics && Count[diagram, _Line, Infinity] >= 4 &&
                Count[diagram, _Rectangle, Infinity] === 1 &&
                Cases[diagram, Text[Style[p_Integer, 10, ___], _, _] :> p, Infinity] === Range[3] &&
                Cases[unlabelled, Text[Style[p_Integer, 10, ___], _, _] :> p, Infinity] === {} &&
                ! FreeQ[diagram, Line[{{{.75, -1}, {1.25, -3}}, {{.75, -3}, {1.25, -1}}}]]
        ]
    ],
    True,
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

(* The named Fourier constructor deliberately uses machine numbers, including at larger sizes. *)
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
    {QuantumInterferometer::notnumeric},
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
        With[{diagram = qi["Diagram", "ShowParameters" -> True], hidden = qi["Diagram", "ShowParameters" -> False]},
            {
                Count[diagram, _Rectangle, Infinity] == qi["PhaseShifterCount"],
                Cases[diagram, Tooltip[_, Row[{"(\[Theta], \[Phi]) = ", tf_, __}]] :> tf, Infinity] ==
                    Cases[qi["Elements"], ("BS"[t_, f_] -> _) :> {t, f}],
                Count[diagram, Text[Style[_Row | _NumberForm, 8], _], Infinity] ===
                    qi["BeamSplitterCount"] + qi["PhaseShifterCount"],
                Count[hidden, Text[Style[_Row | _NumberForm, 8], _], Infinity] === 0
            }
        ]
    ],
    {True, True, True, True},
    TestID -> "QI-diagram-labels-beam-splitter-parameters"
]

VerificationTest[
    SeedRandom[18];
    With[{qi1 = QuantumInterferometer["Random"[4], Method -> "Reck"], qi2 = QuantumInterferometer["Random"[4]]},
        qi2[qi1]["Depth"] == qi1["Depth"] + qi2["Depth"] &&
            And @@ Map[
                With[{a = QuantumInterferometer[#[[1]], 4], b = QuantumInterferometer[#[[2]], 4], layers = #[[3]]},
                    b[a]["Layers"] === layers &&
                        Max[Abs[b[a]["Unitary"] - b["Unitary"] . a["Unitary"]]] < 10^-12
                ] &,
                {
                    {{"BS"[.4] -> {1, 2}}, {"BS"[.3] -> {3, 4}}, {1, 1}},
                    {{"BS"[.4] -> {1, 2}}, {"BS"[.3] -> {2, 3}}, {1, 2}},
                    {{"BS"[.4] -> {1, 4}}, {"BS"[.3] -> {2, 3}}, {1, 2}},
                    {{"PS"[.2] -> 1}, {"PS"[.3] -> 1}, {1/2, 3/2}}
                }
            ]
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
            With[{outs = Select[Tuples[Range[0, 6], 2], Total[#] <= 6 && EvenQ[Total[#]] &],
                  probs = qi["Probabilities", {0, 0}, "Cutoff" -> 6]},
                AssociationQ[probs] && Sort[Keys[probs]] === Sort[QuditName /@ outs] &&
                    Max[Abs[Lookup[probs, QuditName /@ outs] -
                        (If[Equal @@ #, Tanh[.5] ^ (2 First[#]) / Cosh[.5] ^ 2, 0] & /@ outs)]] < 10^-12 &&
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
                1 - 10^-4 < Total[probs] <= 1 + 10^-12 &&
                    AnyTrue[Keys[probs], Total[#["Name"]] > 1 &] &&
                    Max[Abs[Values[probs] - (qi["Probability", {1, 0} -> #["Name"]] & /@ Keys[probs])]] < 10^-12
            ]
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

(* ---- permanent dispatch: calibrated costs, exact results and numerical invariants ---- *)

$qiPermanentAmplitude = Wolfram`QuantumFramework`SecondQuantization`Interferometer`PackagePrivate`permanentAmplitude;
$qiUseRepeatedPermanentQ = Wolfram`QuantumFramework`SecondQuantization`Interferometer`PackagePrivate`useRepeatedPermanentQ;

VerificationTest[
    With[{w = FourierMatrix[3], symbolic = {{a, b}, {c, d}}},
        {
            $qiPermanentAmplitude[w, {2, 1, 0}, {1, 2, 0}] === Permanent[w[[{1, 2, 2}, {1, 1, 2}]]] / 2,
            $qiUseRepeatedPermanentQ[w, {2, 1, 0}, {1, 2, 0}],
            $qiPermanentAmplitude[symbolic, {1, 1}, {1, 1}] === a d + b c,
            $qiUseRepeatedPermanentQ[symbolic, {1, 1}, {1, 1}]
        }
    ],
    {True, False, True, False},
    TestID -> "QI-permanent-dispatch-preserves-exact-and-symbolic-inputs"
]

VerificationTest[
    {
        $qiPermanentAmplitude[IdentityMatrix[2], {0, 0}, {0, 0}],
        $qiPermanentAmplitude[IdentityMatrix[2], {1, 0}, {0, 0}]
    },
    {1, 0},
    TestID -> "QI-permanent-dispatch-vacuum-and-number-conservation"
]

VerificationTest[
    With[{w = {{1, -1}, {1, 1}} / Sqrt[2]},
        {
            $qiPermanentAmplitude[w, {1, 1}, {1, 1}],
            Abs[$qiPermanentAmplitude[w, {1, 1}, {2, 0}]] ^ 2,
            Abs[$qiPermanentAmplitude[w, {1, 1}, {0, 2}]] ^ 2
        }
    ],
    {0, 1/2, 1/2},
    TestID -> "QI-permanent-dispatch-exact-HOM"
]

VerificationTest[
    With[{s = PadRight[Join[{2}, ConstantArray[1, 10]], 16]},
        $qiUseRepeatedPermanentQ[N[FourierMatrix[16]], s, Reverse[s]]
    ],
    False,
    TestID -> "QI-permanent-dispatch-avoids-mild-collision-slowdown"
]

VerificationTest[
    With[{s = PadRight[{8, 8}, 16]},
        $qiUseRepeatedPermanentQ[N[FourierMatrix[16]], s, Reverse[s]]
    ],
    True,
    TestID -> "QI-permanent-dispatch-retains-strong-collision-route"
]

VerificationTest[
    With[{s = PadRight[{6, 6}, 128]},
        $qiUseRepeatedPermanentQ[N[FourierMatrix[128]], s, Reverse[s]]
    ],
    True,
    TestID -> "QI-permanent-dispatch-accounts-for-mode-conversion-cost"
]

VerificationTest[
    With[{s = PadRight[{4, 4}, 16]},
        {
            $qiUseRepeatedPermanentQ[Developer`ToPackedArray[N[FourierMatrix[16]], Complex], s, Reverse[s]],
            $qiUseRepeatedPermanentQ[N[FourierMatrix[16], 35], s, Reverse[s]]
        }
    ],
    {False, True},
    TestID -> "QI-permanent-dispatch-distinguishes-numerical-precisions"
]

VerificationTest[
    With[{w = N[DiagonalMatrix[{Exp[37 I / 100], Exp[- I / 5]}]], s = {24, 0}},
        {
            $qiUseRepeatedPermanentQ[w, s, s],
            Abs[$qiPermanentAmplitude[w, s, s] - w[[1, 1]] ^ 24] < 10^-13
        }
    ],
    {True, True},
    TestID -> "QI-permanent-dispatch-preserves-large-bunching-accuracy"
]

VerificationTest[
    With[{m = 16, n = 12},
        With[{w = N[FourierMatrix[m]], s = PadRight[ConstantArray[1, n], m], t = PadLeft[{n}, m]},
            Abs[$qiPermanentAmplitude[w, s, t] - Sqrt[n!] Times @@ w[[-1, Range[n]]]] < 10^-12
        ]
    ],
    True,
    TestID -> "QI-permanent-dispatch-output-collision-normalization"
]

(* Exact four-mode references include every input/output pair with one through four photons:
   16 + 100 + 400 + 1225 = 1741 pairs, checked at machine and 35-digit precision. *)
VerificationTest[
    With[{w = DiagonalMatrix[{1, I, -1, -I}] . FourierMatrix[4] . DiagonalMatrix[{1, -I, -1, I}]},
        And @@ Flatten @ Table[
            With[{patterns = FrobeniusSolve[ConstantArray[1, 4], n]},
                And @@ Flatten @ Table[
                    With[{
                        rows = Catenate @ MapIndexed[ConstantArray[First[#2], #1] &, t],
                        cols = Catenate @ MapIndexed[ConstantArray[First[#2], #1] &, s]
                    },
                        With[{a = Permanent[w[[rows, cols]]] / Sqrt[Times @@ (s!) Times @@ (t!)]},
                            Abs[$qiPermanentAmplitude[N[w], s, t] - N[a]] < 10^-12 &&
                            Abs[$qiPermanentAmplitude[N[w, 35], s, t] - N[a, 35]] < 10^-28
                        ]
                    ],
                    {s, patterns}, {t, patterns}
                ]
            ],
            {n, 1, 4}
        ]
    ],
    True,
    TestID -> "QI-permanent-dispatch-all-small-occupation-pairs"
]

VerificationTest[
    With[{qi = QuantumInterferometer["Fourier"[3]], s = {2, 1, 1}, outs = FrobeniusSolve[{1, 1, 1}, 4]},
        Abs[Total[Abs[qi["Amplitude", s -> #]] ^ 2 & /@ outs] - 1] < 10^-12
    ],
    True,
    TestID -> "QI-permanent-dispatch-scalar-amplitudes-normalize-sector"
]

VerificationTest[
    With[{w = SparseArray[N[DiagonalMatrix[{Exp[I / 3], Exp[I / 5], Exp[I / 7], Exp[I / 11]}]]], s = {4, 4, 0, 0}},
        Abs[$qiPermanentAmplitude[w, s, s] - w[[1, 1]] ^ 4 w[[2, 2]] ^ 4] < 10^-12
    ],
    True,
    TestID -> "QI-permanent-dispatch-sparse-unitary"
]

VerificationTest[
    With[{exact = FourierMatrix[4]},
        With[{w = ReplacePart[N[exact, 50], {1, 1} -> N[exact[[1, 1]]]], s = {2, 1, 1, 0}, t = {0, 1, 1, 2}},
            Abs[$qiPermanentAmplitude[w, s, t] - N[Permanent[exact[[{2, 3, 4, 4}, {1, 1, 2, 3}]]] / 2, 50]] < 10^-12
        ]
    ],
    True,
    TestID -> "QI-permanent-dispatch-mixed-precision-unitary"
]

(* ---- occupied-mode reduction: rectangular blocks, precision and invariants ---- *)

$qiRepeatedPermanent = Wolfram`QuantumFramework`SecondQuantization`Interferometer`PackagePrivate`repeatedPermanent;

VerificationTest[
    With[{w = Array[z, {3, 3}]},
        Expand[$qiRepeatedPermanent[w, {2, 0, 1}, {1, 1, 1}] -
            $qiRepeatedPermanent[w[[All, {1, 3}]], {2, 1}, {1, 1, 1}]] === 0
    ],
    True,
    TestID -> "QI-occupied-modes-symbolic-rectangular-identity"
]

VerificationTest[
    BlockRandom[SeedRandom[4871];
        With[{w = RandomVariate[CircularUnitaryMatrixDistribution[16]],
              s = ReplacePart[ConstantArray[0, 16], {2 -> 9, 11 -> 7}],
              t = ReplacePart[ConstantArray[0, 16], {3 -> 6, 8 -> 4, 14 -> 6}]},
            With[{ref = Permanent[w[[Join[ConstantArray[3, 6], ConstantArray[8, 4], ConstantArray[14, 6]],
                    Join[ConstantArray[2, 9], ConstantArray[11, 7]]]]] / Sqrt[9! 7! 6! 4! 6!]},
                Abs[ref] > 10^-14 &&
                Abs[$qiPermanentAmplitude[w, s, t] - ref] < 10^-10 Abs[ref] &&
                Abs[$qiPermanentAmplitude[ConjugateTranspose[w], t, s] - Conjugate[ref]] < 10^-10 Abs[ref]
            ]
        ]
    ],
    True,
    TestID -> "QI-occupied-modes-noncontiguous-rectangular-both-orientations"
]

VerificationTest[
    With[{w = N[FourierMatrix[16]],
          s = ReplacePart[ConstantArray[0, 16], {2 -> 9, 11 -> 7}],
          t = ReplacePart[ConstantArray[0, 16], {3 -> 6, 8 -> 4, 14 -> 6}]},
        Module[{amplitudeProbe, precisionProbe},
            (* SetPrecision is invisible to Trace; instrument a temporary copy of the live definition. *)
            precisionProbe[a_, p_] := (Sow[Dimensions[a], "occupied-promotion"]; SetPrecision[a, p]);
            With[{f = $qiPermanentAmplitude},
                DownValues[amplitudeProbe] = DownValues[f] /. f -> amplitudeProbe /. SetPrecision -> precisionProbe
            ];
            Reap[amplitudeProbe[w, s, t], "occupied-promotion"][[2]] === {{{3, 2}}}
        ]
    ],
    True,
    TestID -> "QI-occupied-modes-promotes-only-occupied-block"
]

VerificationTest[
    With[{w = ReplacePart[IdentityMatrix[16], {16, 16} -> N[Exp[I/7]]],
          s = ReplacePart[ConstantArray[0, 16], {2 -> 8, 6 -> 8}]},
        With[{a = $qiPermanentAmplitude[w, s, s]},
            Precision[a] === MachinePrecision && Abs[a - 1] < 10^-12
        ]
    ],
    True,
    TestID -> "QI-occupied-modes-retains-machine-precision-outside-block"
]

VerificationTest[
    With[{w = ReplacePart[N[DiagonalMatrix[Exp[I Range[16] / 17]], 60],
                {16, 16} -> N[Exp[16 I/17], 30]],
          s = ReplacePart[ConstantArray[0, 16], {2 -> 8, 6 -> 8}]},
        With[{a = $qiPermanentAmplitude[w, s, s]},
            Abs[Precision[a] - Precision[w]] < 10^-8 && Abs[a - N[Exp[64 I/17], 80]] < 10^-27
        ]
    ],
    True,
    TestID -> "QI-occupied-modes-retains-lower-arbitrary-precision-outside-block"
]

VerificationTest[
    With[{w = SparseArray[N[DiagonalMatrix[Exp[I Range[16] / 17]]]],
          s = ReplacePart[ConstantArray[0, 16], {2 -> 8, 6 -> 8}],
          t = ReplacePart[ConstantArray[0, 16], {2 -> 8, 7 -> 8}]},
        Abs[$qiPermanentAmplitude[w, s, s] - Exp[64 I/17]] < 10^-12 &&
        Abs[$qiPermanentAmplitude[w, s, t]] < 10^-12
    ],
    True,
    TestID -> "QI-occupied-modes-sparse-complex-unitary-and-zero-amplitude"
]

VerificationTest[
    BlockRandom[SeedRandom[4872];
        With[{v = RandomVariate[CircularUnitaryMatrixDistribution[2]]},
            With[{w = ReplacePart[IdentityMatrix[32], Thread[Tuples[{2, 17}, 2] -> Flatten[v]]],
                  s = ReplacePart[ConstantArray[0, 32], {2 -> 8, 17 -> 8}],
                  t = ReplacePart[ConstantArray[0, 32], {2 -> 6, 17 -> 10}]},
                Abs[$qiPermanentAmplitude[w, s, t] - $qiPermanentAmplitude[v, {8, 8}, {6, 10}]] < 10^-12
            ]
        ]
    ],
    True,
    TestID -> "QI-occupied-modes-vacuum-padding-and-port-embedding"
]

VerificationTest[
    BlockRandom[SeedRandom[4873];
        With[{qi = QuantumInterferometer[RandomVariate[CircularUnitaryMatrixDistribution[3]]],
              s = {6, 6, 0}, outs = FrobeniusSolve[{1, 1, 1}, 12]},
            With[{amps = qi["Amplitude", s -> #] & /@ outs, adj = qi["Dagger"]},
                {
                    Abs[Total[Abs[amps] ^ 2] - 1] < 10^-11,
                    Max[Abs[(adj["Amplitude", # -> s] & /@ outs) - Conjugate[amps]]] < 10^-12,
                    AnyTrue[outs, $qiUseRepeatedPermanentQ[qi["Unitary"], s, #] &]
                }
            ]
        ]
    ],
    {True, True, True},
    TestID -> "QI-occupied-modes-full-sector-normalization-and-adjoint"
]

VerificationTest[
    With[{w = FourierMatrix[3]},
        {
            $qiPermanentAmplitude[w, {4, 0, 0}, {1, 2, 1}] === Sqrt[12] / 9,
            $qiUseRepeatedPermanentQ[w, {4, 0, 0}, {1, 2, 1}],
            Module[{selectorProbe, planProbe, preprocessingProbe},
                SetAttributes[preprocessingProbe, HoldAll];
                preprocessingProbe[op_, args___] := (Sow[HoldForm[op], "exact-preprocessing"]; op[args]);
                With[{f = $qiUseRepeatedPermanentQ,
                      plan = Wolfram`QuantumFramework`SecondQuantization`Interferometer`PackagePrivate`repeatedPermanentPlan},
                    DownValues[selectorProbe] = DownValues[f] /. {f -> selectorProbe, plan -> planProbe};
                    DownValues[planProbe] = DownValues[plan] /. plan -> planProbe /.
                        HoldPattern[(h : Count | Position | Pick | Part | Developer`PackedArrayQ | Developer`ToPackedArray | SetPrecision)[args___]] :>
                            preprocessingProbe[h, args]
                ];
                Reap[selectorProbe[w, {4, 0, 0}, {1, 2, 1}], "exact-preprocessing"] === {False, {}}
            ]
        }
    ],
    {True, False, True},
    TestID -> "QI-occupied-modes-exact-amplitude-and-selector-guard"
]

VerificationTest[
    With[{w = N[FourierMatrix[16], 40], u = N[FourierMatrix[16][[{5, 9}, {2, 13}]], 80],
          s = ReplacePart[ConstantArray[0, 16], {2 -> 8, 13 -> 8}],
          t = ReplacePart[ConstantArray[0, 16], {5 -> 10, 9 -> 6}]},
        With[{ref = Sqrt[10! 6! / (8! 8!)] Sum[
            Binomial[8, k] Binomial[8, 10 - k] u[[1, 1]] ^ k u[[2, 1]] ^ (8 - k)
                u[[1, 2]] ^ (10 - k) u[[2, 2]] ^ (k - 2), {k, 2, 8}]},
            Abs[ref] > 10^-20 && Abs[$qiPermanentAmplitude[w, s, t] - ref] < 10^-28 Abs[ref]
        ]
    ],
    True,
    TestID -> "QI-occupied-modes-arbitrary-precision-independent-coefficient-reference"
]

VerificationTest[
    With[{w = N[N[FourierMatrix[16], 80]], s = PadRight[{4, 4}, 16]},
        With[{packed = Developer`ToPackedArray[w, Complex]},
            {
                $qiUseRepeatedPermanentQ[w, s, Reverse[s]],
                $qiUseRepeatedPermanentQ[packed, s, Reverse[s]],
                $qiUseRepeatedPermanentQ[Developer`FromPackedArray[packed], s, Reverse[s]],
                Max[Abs[($qiPermanentAmplitude[#, s, Reverse[s]] & /@ {w, packed, Developer`FromPackedArray[packed]}) -
                    $qiPermanentAmplitude[packed, s, Reverse[s]]]] < 10^-12
            }
        ]
    ],
    {True, False, False, True},
    TestID -> "QI-occupied-modes-dispatch-distinguishes-mixed-and-uniform-storage"
]

VerificationTest[
    With[{s = PadRight[Join[{2}, ConstantArray[1, 10]], 16], free = PadRight[ConstantArray[1, 12], 16]},
        {
            $qiUseRepeatedPermanentQ[N[FourierMatrix[16], 35], s, Reverse[s]],
            $qiUseRepeatedPermanentQ[N[FourierMatrix[16], 35], free, Reverse[free]],
            $qiUseRepeatedPermanentQ[FourierMatrix[16], s, Reverse[s]]
        }
    ],
    {True, False, False},
    TestID -> "QI-occupied-modes-arbitrary-precision-collision-and-exact-guards"
]

(* ---- coverage added after the test-meaning audit: independent weights and missing routes ---- *)

$qiRepeatedHafnian = Wolfram`QuantumFramework`SecondQuantization`Interferometer`PackagePrivate`repeatedHafnian;

VerificationTest[
    With[{qi = QuantumInterferometer[{"S"[.6] -> 1}, 1]},
        With[{probs = qi["Probabilities", {0}]},
            AssociationQ[probs] && Length[probs] > 0 && With[{top = Max[First[Flatten[{#["Name"]}]] & /@ Keys[probs]]},
                With[{outs = List /@ Range[0, top, 2],
                      ref = Table[Binomial[2 k, k] Tanh[.6] ^ (2 k) / (4 ^ k Cosh[.6]), {k, 0, top / 2}]},
                    top > 0 && Sort[Keys[probs]] === Sort[QuditName /@ outs] &&
                        Max[Abs[Lookup[probs, QuditName /@ outs] - ref]] < 10^-12 &&
                        0 < 1 - Total[probs] < 10^-4 && Total[Most[ref]] < 1 - 10^-4 &&
                        Max[Abs[Lookup[qi["Probabilities", {0}, "Cutoff" -> top + 4], QuditName /@ outs] - ref]] < 10^-12
                ]
            ]
        ]
    ],
    True,
    TestID -> "QI-auto-squeezed-cutoff-keeps-absolute-weights-and-bounded-tail"
]

VerificationTest[
    With[{probs = QuantumInterferometer[{"S"[1.2] -> 1}, 1]["Probabilities", {0}],
          outs = List /@ Range[0, 16, 2],
          ref = Table[Binomial[2 k, k] Tanh[1.2] ^ (2 k) / (4 ^ k Cosh[1.2]), {k, 0, 8}]},
        AssociationQ[probs] && Sort[Keys[probs]] === Sort[QuditName /@ outs] &&
            Max[Abs[Lookup[probs, QuditName /@ outs] - ref]] < 10^-12 &&
            Abs[Total[probs] - Total[ref]] < 10^-12 && 0 < Total[probs] < 1 - 10^-4
    ],
    True,
    {QuantumInterferometer::autocutoff},
    TestID -> "QI-auto-squeezed-cutoff-cap-warns-without-renormalizing"
]

VerificationTest[
    With[{w = {{a, b}, {b, c}}},
        {
            $qiRepeatedHafnian[w, {0, 0}],
            $qiRepeatedHafnian[w, {1, 2}],
            Expand[$qiRepeatedHafnian[w, {2, 2}] - a c - 2 b ^ 2]
        }
    ],
    {1, 0, 0},
    TestID -> "QI-repeated-hafnian-symbolic-matchings-vacuum-and-odd-parity"
]

(* S(xi)|0> and S(xi)|1> give independent factorial coefficients. Total input plus
   output photon numbers straddle the generic/repeated hafnian boundary. *)
VerificationTest[
    And @@ Flatten @ Table[
        With[{r = If[p === MachinePrecision, N[2/5], N[2/5, p]],
              phi = If[p === MachinePrecision, N[7/13], N[7/13, p]]},
            With[{qi = QuantumInterferometer[{"S"[r Exp[I phi]] -> 1}, 1]},
                Table[
                    With[{n = 16 + in + 2 j},
                        With[{ref = Sqrt[n!] (- Exp[I phi] Tanh[r]) ^ ((n - in) / 2) /
                                (2 ^ ((n - in) / 2) ((n - in) / 2)! Cosh[r] ^ (in + 1/2)),
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
    TestID -> "QI-high-photon-squeezing-phase-precision-and-nonvacuum-reference"
]

VerificationTest[
    With[{w = N[FourierMatrix[16]], s = PadRight[ConstantArray[1, 12], 16], t = PadLeft[{12}, 16]},
        Module[{amplitudeProbe, repeatedProbe},
            repeatedProbe[v_, in_, out_] := (Sow[{Dimensions[v], Times @@ (in + 1)}, "permanent-work"];
                $qiRepeatedPermanent[v, in, out]);
            With[{f = $qiPermanentAmplitude, r = $qiRepeatedPermanent},
                DownValues[amplitudeProbe] = DownValues[f] /. f -> amplitudeProbe /. r -> repeatedProbe
            ];
            With[{result = Reap[{amplitudeProbe[w, s, t], amplitudeProbe[ConjugateTranspose[w], t, s]}, "permanent-work"],
                  ref = Sqrt[12!] Times @@ w[[-1, Range[12]]]},
                result[[2]] === {{{{12, 1}, 13}, {{12, 1}, 13}}} && Abs[ref] > 10^-6 &&
                    Max[Abs[result[[1]] - {ref, Conjugate[ref]}]] < 10^-12 Abs[ref]
            ]
        ]
    ],
    True,
    TestID -> "QI-permanent-amplitude-executes-cheaper-orientation-in-both-directions"
]

VerificationTest[
    With[{w = DiagonalMatrix[{1, I}] . {{3/5, -4 Exp[-I Pi/3]/5}, {4 Exp[I Pi/3]/5, 3/5}}},
        With[{out = QuantumInterferometer[w]["Coherent"[{a, b}]]},
            MatchQ[out, "Coherent"[_List]] &&
                FullSimplify[First[out] - {3 a/5 - 4 Exp[-I Pi/3] b/5, 4 I Exp[I Pi/3] a/5 + 3 I b/5}] === {0, 0}
        ]
    ],
    True,
    TestID -> "QI-coherent-symbolic-transfer-is-phase-sensitive"
]

VerificationTest[
    With[{w = N[DiagonalMatrix[{1, I}] . {{3/5, -4 Exp[-I Pi/3]/5}, {4 Exp[I Pi/3]/5, 3/5}}],
          a = N[{(2 + I)/7, (-1 + 2 I)/6}], outs = Select[Tuples[Range[0, 6], 2], Total[#] <= 6 &]},
        With[{qi = QuantumInterferometer[w], beta = w . a},
            With[{probs = qi["Probabilities", "Coherent"[a], "Cutoff" -> 6],
                  ref = Exp[-Total[Abs[a] ^ 2]/2] (Times @@ (beta ^ # / Sqrt[#!]) & /@ outs)},
                AssociationQ[probs] && Sort[Keys[probs]] === Sort[QuditName /@ outs] &&
                    Max[Abs[(qi["Amplitude", "Coherent"[a] -> #] & /@ outs) - ref]] < 10^-12 &&
                    Max[Abs[Lookup[probs, QuditName /@ outs] - Abs[ref] ^ 2]] < 10^-12 &&
                    Abs[Total[probs] - Exp[-Total[Abs[a] ^ 2]] Sum[Total[Abs[a] ^ 2] ^ k / k!, {k, 0, 6}]] < 10^-12
            ]
        ]
    ],
    True,
    TestID -> "QI-coherent-Fock-amplitudes-and-complete-probabilities-match-Poisson"
]

VerificationTest[
    With[{w = N[DiagonalMatrix[{1, I}] . {{3/5, -4 Exp[-I Pi/3]/5}, {4 Exp[I Pi/3]/5, 3/5}}],
          a = N[{(2 + I)/7, (-1 + 2 I)/6}], g = N[{(1 + 2 I)/10, (-2 + I)/10}], alpha = .4 + .2 I},
        With[{qi = QuantumInterferometer[w]},
            With[{inputs = {"Coherent"[a], {1, 1}, {2, "Coherent"[alpha]}},
                  refs = {
                    Exp[-(Total[Abs[a] ^ 2] + Total[Abs[g] ^ 2])/2 + Conjugate[g] . w . a],
                    Exp[-Total[Abs[g] ^ 2]/2] Times @@ (Conjugate[g] . w),
                    Exp[-(Abs[alpha] ^ 2 + Total[Abs[g] ^ 2])/2 + alpha Conjugate[g] . w[[All, 2]]] *
                        (Conjugate[g] . w[[All, 1]]) ^ 2 / Sqrt[2]
                  }},
                Abs[Im[N[First[refs]]]] > 10^-4 &&
                    Max[Abs[MapThread[qi["Amplitude", #1 -> "Coherent"[g]] - #2 &, {inputs, refs}]]] < 10^-12 &&
                    Max[Abs[MapThread[qi["Probability", #1 -> "Coherent"[g]] - Abs[#2] ^ 2 &, {inputs, refs}]]] < 10^-12
            ]
        ]
    ],
    True,
    TestID -> "QI-coherent-output-overlaps-for-coherent-Fock-and-mixed-inputs"
]

(* |2,alpha> maps to (W_1.a^dag)^2 Exp[alpha W_2.a^dag]|0>/Sqrt[2].
   A degree-four exponential polynomial gives every coefficient through six photons exactly. *)
VerificationTest[
    With[{w = N[DiagonalMatrix[{1, I}] . {{3/5, -4 Exp[-I Pi/3]/5}, {4 Exp[I Pi/3]/5, 3/5}}],
          alpha = .4 + .2 I, outs = Select[Tuples[Range[0, 6], 2], 2 <= Total[#] <= 6 &]},
        Module[{x, y},
            With[{qi = QuantumInterferometer[w],
                  polynomial = Expand[({x, y} . w[[All, 1]]) ^ 2 Sum[(alpha {x, y} . w[[All, 2]]) ^ k/k!, {k, 0, 4}]/Sqrt[2]]},
                With[{ref = (Exp[-Abs[alpha] ^ 2/2] Sqrt[Times @@ (#!)] Coefficient[Coefficient[polynomial, x, #[[1]]], y, #[[2]]] & /@ outs),
                      probs = qi["Probabilities", {2, "Coherent"[alpha]}, "Cutoff" -> 6]},
                    AssociationQ[probs] && Sort[Keys[probs]] === Sort[QuditName /@ outs] &&
                        Max[Abs[(qi["Amplitude", {2, "Coherent"[alpha]} -> #] & /@ outs) - ref]] < 10^-12 &&
                        Max[Abs[Lookup[probs, QuditName /@ outs] - Abs[ref] ^ 2]] < 10^-12 &&
                        Max[Abs[qi["Amplitude", {2, "Coherent"[alpha]} -> #] & /@ {{0, 0}, {1, 0}, {0, 1}}]] < 10^-12 &&
                        Abs[Total[probs] - Exp[-Abs[alpha] ^ 2] Sum[Abs[alpha] ^ (2 k)/k!, {k, 0, 4}]] < 10^-12
                ]
            ]
        ]
    ],
    True,
    TestID -> "QI-mixed-Fock-coherent-amplitudes-match-independent-polynomial"
]

VerificationTest[
    With[{alpha = (3 + 2 I)/5, qi = QuantumInterferometer[{}, 1]},
        With[{probs = qi["Probabilities", "Coherent"[{alpha}]], lambda = Abs[alpha] ^ 2},
            AssociationQ[probs] && Length[probs] > 0 && With[{top = Max[First[Flatten[{#["Name"]}]] & /@ Keys[probs]]},
                With[{outs = List /@ Range[0, top], ref = Table[Exp[-lambda] lambda ^ n/n!, {n, 0, top}]},
                    top > 0 && Sort[Keys[probs]] === Sort[QuditName /@ outs] &&
                        Max[Abs[Lookup[probs, QuditName /@ outs] - ref]] < 10^-12 &&
                        1 - 10^-4 <= Total[probs] < 1 && Total[Most[ref]] < 1 - 10^-4
                ]
            ]
        ]
    ],
    True,
    TestID -> "QI-coherent-auto-cutoff-preserves-absolute-Poisson-tail"
]

VerificationTest[
    With[{probs = QuantumInterferometer[{}, 1]["Probabilities", "Coherent"[{alpha}], "Cutoff" -> 3]},
        AssociationQ[probs] && Sort[Keys[probs]] === Sort[QuditName /@ List /@ Range[0, 3]] &&
            FullSimplify[Values[probs] - Table[Exp[-Abs[alpha] ^ 2] Abs[alpha] ^ (2 n)/n!, {n, 0, 3}]] === ConstantArray[0, 4]
    ],
    True,
    TestID -> "QI-coherent-symbolic-amplitude-accepts-explicit-cutoff"
]

VerificationTest[
    QuantumInterferometer[{"BS"[.4, .2] -> {1, 2}}, 2]["Coherent"[{.3}]],
    $Failed,
    {QuantumInterferometer::occ},
    TestID -> "QI-coherent-direct-input-rejects-wrong-length"
]

VerificationTest[
    QuantumInterferometer[{"BS"[.4, .2] -> {1, 2}}, 2]["Amplitude", "Coherent"[{.3}] -> {0, 0}],
    $Failed,
    {QuantumInterferometer::occ},
    TestID -> "QI-coherent-amplitude-rejects-wrong-input-length"
]

VerificationTest[
    QuantumInterferometer[{"BS"[.4, .2] -> {1, 2}}, 2]["Amplitude", "Coherent"[{.3, .2}] -> "Coherent"[{.1}]],
    $Failed,
    {QuantumInterferometer::occ},
    TestID -> "QI-coherent-amplitude-rejects-wrong-output-length"
]

VerificationTest[
    QuantumInterferometer[{"S"[.3] -> 1}, 1]["Coherent"[{.2}]],
    $Failed,
    {QuantumInterferometer::squeezed},
    TestID -> "QI-coherent-direct-input-rejects-active-network"
]

VerificationTest[
    QuantumInterferometer[{"S"[.3] -> 1}, 1]["Amplitude", "Coherent"[{.2}] -> {0}],
    $Failed,
    {QuantumInterferometer::squeezed},
    TestID -> "QI-coherent-amplitude-rejects-active-network"
]

VerificationTest[
    QuantumInterferometer[{"S"[.3] -> 1}, 1]["Probabilities", "Coherent"[{.2}], "Cutoff" -> 4],
    $Failed,
    {QuantumInterferometer::squeezed},
    TestID -> "QI-coherent-probabilities-reject-active-network"
]

VerificationTest[
    {QuantumInterferometer[{}, 1]["Probabilities", "Coherent"[{alpha}]],
        QuantumInterferometer[{"S"[r] -> 1}, 1]["Probabilities", {0}]},
    {$Failed, $Failed},
    {QuantumInterferometer::symcutoff, QuantumInterferometer::symcutoff},
    TestID -> "QI-symbolic-coherent-and-squeezed-auto-cutoff-needs-numbers"
]

VerificationTest[
    With[{d = 4},
        With[{v = (FockState[{2, 0}, d] + I FockState[{0, 2}, d])/Sqrt[2], u = FockState[{1, 0}, d]},
            With[{rho = QuantumState[.7 Normal[v["DensityMatrix"]] + .3 Normal[u["DensityMatrix"]], {d, d}],
                  qi = QuantumInterferometer[{"BS"[.41, .23] -> {1, 2}, "PS"[.37] -> 2}, 2]},
                With[{actual = qi[rho], expected = PhaseShiftOperator[.37, d, {2}][
                        BeamSplitterOperator[{.41, .23}, d, {1, 2}, Method -> "Recurrence"][rho]]},
                    With[{matrix = Normal[actual["DensityMatrix"]]},
                        rho["StateType"] === "Matrix" && actual["Dimensions"] === {d, d} &&
                            Max[Abs[Im[Normal[rho["DensityMatrix"]]]]] > .1 &&
                            Max[Abs[matrix - Normal[expected["DensityMatrix"]]]] < 10^-12 &&
                            Max[Abs[matrix - ConjugateTranspose[matrix]]] < 10^-12 &&
                            Abs[Tr[matrix] - 1] < 10^-12 && Min[Re[Eigenvalues[matrix]]] > -10^-12
                    ]
                ]
            ]
        ]
    ],
    True,
    TestID -> "QI-passive-mixed-state-matches-explicit-Fock-operators"
]

VerificationTest[
    With[{d = 8, xi = .25 Exp[.4 I]},
        With[{state = (FockState[{0, 0}, d] + I FockState[{1, 1}, d])/Sqrt[2],
              qi = QuantumInterferometer[{"S"[xi] -> 1, "BS"[.37, .29] -> {1, 2}, "PS"[.3] -> 2}, 2]},
            With[{actual = qi[state], expected = PhaseShiftOperator[.3, d, {2}][
                    BeamSplitterOperator[{.37, .29}, d, {1, 2}, Method -> "Recurrence"][SqueezeOperator[xi, d, {1}][state]]]},
                actual["Dimensions"] === {d, d} &&
                    Max[Abs[Normal[actual["StateVector"]] - Normal[expected["StateVector"]]]] < 10^-12 &&
                    .999 < actual["Norm"] <= 1 + 10^-12
            ]
        ]
    ],
    True,
    TestID -> "QI-active-superposition-matches-explicit-Fock-operators"
]

VerificationTest[
    With[{d = 8, xi = .25 Exp[.4 I]},
        With[{v = (FockState[{0, 0}, d] + I FockState[{1, 1}, d])/Sqrt[2], u = FockState[{1, 0}, d]},
            With[{rho = QuantumState[.6 Normal[v["DensityMatrix"]] + .4 Normal[u["DensityMatrix"]], {d, d}],
                  qi = QuantumInterferometer[{"S"[xi] -> 1, "BS"[.37, .29] -> {1, 2}, "PS"[.3] -> 2}, 2]},
                With[{actual = qi[rho], expected = PhaseShiftOperator[.3, d, {2}][
                        BeamSplitterOperator[{.37, .29}, d, {1, 2}, Method -> "Recurrence"][SqueezeOperator[xi, d, {1}][rho]]]},
                    With[{matrix = Normal[actual["DensityMatrix"]]},
                        actual["StateType"] === "Matrix" && actual["Dimensions"] === {d, d} &&
                            Max[Abs[matrix - Normal[expected["DensityMatrix"]]]] < 10^-12 &&
                            Abs[Tr[matrix] - Tr[Normal[expected["DensityMatrix"]]]] < 10^-12 &&
                            .999 < Re[Tr[matrix]] <= 1 + 10^-12 && Min[Re[Eigenvalues[matrix]]] > -10^-12
                    ]
                ]
            ]
        ]
    ],
    True,
    TestID -> "QI-active-mixed-state-preserves-coherences-and-truncated-trace"
]

VerificationTest[
    With[{plot = QuantumInterferometer[{{3, -4}, {4, 3}}/5]["MatrixPlot", ImageSize -> 111]},
        With[{rasters = Cases[plot, Raster[data_, ___] :> data, Infinity],
              intensity = Map[{1 - .95 #, 1 - .75 #, 1 - .4 #} &, {{16/25, 9/25}, {9/25, 16/25}}, {2}],
              phase = {{{.95, .285, .285}, {.95, .285, .285}}, {{.95, .285, .285}, {.285, .95, .95}}}},
            Length[rasters] === 2 && Dimensions[rasters] === {2, 2, 2, 3} &&
                Max[Abs[rasters - {intensity, phase}]] < 10^-12 &&
                Cases[plot, g_Graphics :> (ImageSize /. Options[g, ImageSize]), Infinity] === {111, 111}
        ]
    ],
    True,
    TestID -> "QI-numeric-matrix-plot-renders-intensity-phase-and-options"
]

VerificationTest[
    {
        QuantumInterferometer[IdentityMatrix[2], Method -> "Unknown"],
        QuantumInterferometer[{"BS"[.4] -> {1, 1}}, 2],
        QuantumInterferometer[{"PS"[.3] -> 3}, 2],
        QuantumInterferometer[{}, 2][QuantumInterferometer[{}, 3]],
        QuantumInterferometer[{}, 2][FockState[{1}, 3]]
    },
    ConstantArray[$Failed, 5],
    {QuantumInterferometer::method, QuantumInterferometer::elem, QuantumInterferometer::modes,
        QuantumInterferometer::compose, QuantumInterferometer::state},
    TestID -> "QI-rejects-invalid-method-elements-composition-and-state-dimensions"
]

(* The fallback keeps exact occupied entries exact even beside unused approximate modes. *)
VerificationTest[
    With[{w = ReplacePart[IdentityMatrix[4], {4, 4} -> N[Exp[I/7]]], in = {1, 0, 0, 0}},
        With[{a = $qiPermanentAmplitude[w, in, in]},
            {Precision[w] === MachinePrecision, $qiUseRepeatedPermanentQ[w, in, in], a, Precision[a]}
        ]
    ],
    {True, False, 1, Infinity},
    TestID -> "QI-permanent-fallback-preserves-exact-occupied-block"
]

(* Enumerating the six assignments gives 2 b^2 f + 4 b c e; the zero-output row contributes nothing. *)
VerificationTest[
    Expand[$qiRepeatedPermanent[{{a, d}, {b, e}, {c, f}}, {2, 1}, {0, 2, 1}] - 2 b^2 f - 4 b c e],
    0,
    TestID -> "QI-repeated-permanent-hoisted-rows-exact-coefficient"
]

(* Each regression below targets a symbolic correctness or validation failure. *)
VerificationTest[
    With[{eps = 1/10^12},
        With[{bs = QuantumInterferometer[{"BS"[eps, 0] -> {1, 2}}, 2],
              ps = QuantumInterferometer[{"PS"[eps] -> 1}, 1],
              sq = QuantumInterferometer[{"S"[eps] -> 1}, 1]},
            bs["Amplitude", {1, 0} -> {0, 1}] === Sin[eps] &&
            ps["Amplitude", {1} -> {1}] === Exp[I eps] &&
            sq["SqueezerCount"] === 1 && Precision[sq["Amplitude", {0} -> {2}]] === Infinity &&
            Abs[N[sq["Amplitude", {0} -> {2}]/(-Sqrt[Sech[eps]] Tanh[eps]/Sqrt[2]), 40] - 1] < 10^-30
        ]
    ],
    True,
    TestID -> "QI-fix-1-exact-small-elements-are-retained"
]

VerificationTest[
    With[{eps = 1/10^12},
        With[{w = {{Sqrt[1 - eps^2], -eps}, {eps, Sqrt[1 - eps^2]}}},
            With[{qi = QuantumInterferometer[w]},
                qi["Unitary"] === w && qi["Amplitude", {1, 0} -> {0, 1}] === eps &&
                FullSimplify[(qi["Amplitude", {nx, 0} -> {0, nx}] /. nx -> 1) - eps] === 0 &&
                AllTrue[Flatten[FullSimplify[QuantumInterferometer[qi["Elements"], 2]["Unitary"] - w]], PossibleZeroQ]
            ]
        ]
    ],
    True,
    TestID -> "QI-fix-1-exact-small-matrix-mesh-and-symbolic-route-agree"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[0, Pi/7] -> {1, 2}, "PS"[2 Pi] -> 1, "S"[0] -> 2}, 2]},
        qi["Elements"] === {} && qi["Unitary"] === IdentityMatrix[2] &&
        qi["Amplitude", {1, 2} -> {1, 2}] === 1
    ],
    True,
    TestID -> "QI-fix-1-exact-identities-still-drop"
]

VerificationTest[
    Map[Function[p,
        With[{phase = N[1/10^12, p]},
            With[{qi = QuantumInterferometer[{"PS"[phase] -> 1}, 1, "DropIdentities" -> False, Tolerance -> 0]},
                With[{actual = qi["Amplitude", {nx} -> {nx}] /. nx -> 2, reference = Exp[2 I phase]},
                    SameQ[Precision[actual], Precision[reference]] &&
                    Abs[(actual - 1)/(reference - 1) - 1] < If[p === MachinePrecision, 10^-12, 10^-50]
                ]
            ]
        ]], {MachinePrecision, 60}],
    {True, True},
    TestID -> "QI-fix-2-symbolic-small-phase-and-precision-survive"
]

VerificationTest[
    With[{parameter = Subscript[\[FormalK], 1]},
        With[{qi = QuantumInterferometer[{"BS"[parameter, 0] -> {1, 2}}, 2]},
            With[{formula = qi["Amplitude", {1, 0} -> {tx, 1 - tx}]},
                FullSimplify[(formula /. {tx -> 1, parameter -> Pi/6}) - Sqrt[3]/2] === 0 &&
                FullSimplify[(formula /. {tx -> 0, parameter -> Pi/6}) - 1/2] === 0
            ]
        ]
    ],
    True,
    TestID -> "QI-fix-3-user-angle-is-not-captured-by-sum-indices"
]

VerificationTest[
    With[{parameter = Subscript[\[FormalK], 1],
          qi = QuantumInterferometer[{"BS"[Pi/6, Pi/7] -> {1, 2}, "BS"[Pi/5, Pi/11] -> {1, 2}}, 2]},
        With[{formula = qi["Amplitude", {parameter, 0} -> {tx, parameter - tx}]},
            And @@ Table[
                TrueQ[TimeConstrained[FullSimplify[
                    Activate[formula /. {parameter -> 2, tx -> t}] - qi["Amplitude", {2, 0} -> {t, 2 - t}]
                ], 5, $TimedOut] === 0],
                {t, 0, 2}
            ]
        ]
    ],
    True,
    TestID -> "QI-fix-3-user-count-is-not-captured-by-sum-indices"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[Pi/6, Pi/5] -> {1, 2}}, 2],
          in = QuantumState[SparseArray[{2 -> Exp[I phi] Sin[theta], 4 -> Cos[theta]}, 9], {3, 3},
              "ParameterSpec" -> {{theta, -Pi, Pi}, {phi, -Pi, Pi}}]},
        With[{out = qi[in], values = <|theta -> Pi/7, phi -> Pi/11|>},
            out["ParameterSpec"] === in["ParameterSpec"] &&
            out[<|theta -> Pi/7|>]["Parameters"] === {phi} &&
            AllTrue[Normal[FullSimplify[out[values]["StateVector"] - qi[in[values]]["StateVector"]]], PossibleZeroQ] &&
            AllTrue[Normal[FullSimplify[out[Pi/7, Pi/11]["StateVector"] - out[values]["StateVector"]]], PossibleZeroQ]
        ]
    ],
    True,
    TestID -> "QI-fix-4-pure-state-retains-parameter-spec-and-binding"
]

VerificationTest[
    {
        QuantumInterferometer[{"PS"[I] -> 1}, 1],
        QuantumInterferometer[{"BS"[I/5, 0] -> {1, 2}}, 2],
        QuantumInterferometer[{"BeamSplitter"[{Pi/7, I/11}] -> {1, 2}}, 2]
    },
    ConstantArray[$Failed, 3],
    {QuantumInterferometer::angles, QuantumInterferometer::angles, QuantumInterferometer::angles, General::stop},
    TestID -> "QI-fix-5-rejects-nonreal-passive-angles"
]

VerificationTest[
    With[{xi = (2 + I)/7, alpha = {(1 + I)/3, (2 - I)/5},
          qi = QuantumInterferometer[{"BS"[Pi/6, Pi/7] -> {1, 2}}, 2]},
        With[{sq = QuantumInterferometer[{"S"[xi] -> 1}, 1]},
            FullSimplify[First[qi["Coherent"[alpha]]] - qi["Unitary"] . alpha] === {0, 0} &&
            Abs[N[sq["Amplitude", {0} -> {2}] + Exp[I Arg[xi]] Sqrt[Sech[Abs[xi]]] Tanh[Abs[xi]]/Sqrt[2], 40]] < 10^-30 &&
            FullSimplify[QuantumInterferometer[{"BS"[th, ph] -> {1, 2}, "PS"[al] -> 1}, 2]["Unitary"] .
                ConjugateTranspose[QuantumInterferometer[{"BS"[th, ph] -> {1, 2}, "PS"[al] -> 1}, 2]["Unitary"]],
                Element[{th, ph, al}, Reals]] === IdentityMatrix[2]
        ]
    ],
    True,
    TestID -> "QI-fix-5-real-angle-validation-keeps-complex-light-parameters"
]

VerificationTest[
    With[{identity = QuantumInterferometer[{}, 1]["Amplitude", {nx} -> {nx}],
          bs = QuantumInterferometer[{"BS"[Pi/6, Pi/7] -> {1, 2}}, 2]["Amplitude", {nx, nx} -> {tx, 2 nx - tx}]},
        {
            (identity /. nx -> #) & /@ {-1, 1/2, I},
            (bs /. Thread[{nx, tx} -> #]) & /@ {{1/2, 1/2}, {I, 0}, {1, I}},
            identity /. nx -> 3,
            bs /. {nx -> 1, tx -> 1}
        }
    ],
    {ConstantArray[0, 3], ConstantArray[0, 3], 1, 1/2},
    TestID -> "QI-fix-6-symbolic-occupations-have-a-discrete-domain"
]

VerificationTest[
    With[{qi = QuantumInterferometer[{"BS"[Pi/6, Pi/7] -> {1, 2}, "BS"[Pi/5, Pi/11] -> {1, 2}}, 2]},
        With[{formula = qi["Amplitude", {nx, 0} -> {tx, nx - tx}]},
            (Activate[formula /. Thread[{nx, tx} -> #]] &) /@
                {{I, 0}, {2, I}, {1/2, 0}, {-1, 0}, {2, -1}, {2, 3}}
        ]
    ],
    ConstantArray[0, 6],
    TestID -> "QI-fix-6-invalid-substitutions-do-not-evaluate-sum-bounds"
]

VerificationTest[
    {
        QuantumInterferometer[{"S"[.2] -> 1}, 1]["Probabilities", {0}, "Cutoff" -> "bad"],
        QuantumInterferometer[{}, 1]["Probabilities", "Coherent"[{1/3}], "Cutoff" -> "bad"],
        QuantumInterferometer[{"BS"[Pi/6, 0] -> {1, 2}}, 2]["Probabilities", {1, 0}, "Cutoff" -> -1]
    },
    ConstantArray[$Failed, 3],
    {QuantumInterferometer::cutoff, QuantumInterferometer::cutoff, QuantumInterferometer::cutoff, General::stop},
    TestID -> "QI-fix-7-probability-distributions-validate-cutoffs"
]

VerificationTest[
    {
        QuantumInterferometer[{"S"[.2] -> 1}, 1][{0}, "Cutoff" -> -1],
        QuantumInterferometer[{"S"[.2] -> 1}, 1][{0}, "Cutoff" -> 1/2],
        QuantumInterferometer[{}, 1][{1}, "Cutoff" -> "bad"]
    },
    ConstantArray[$Failed, 3],
    {QuantumInterferometer::cutoff, QuantumInterferometer::cutoff, QuantumInterferometer::cutoff, General::stop},
    TestID -> "QI-fix-7-amplitude-distributions-validate-cutoffs"
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
    TestID -> "QI-fix-7-zero-and-below-input-cutoffs-remain-valid"
]


EndTestSection[]
