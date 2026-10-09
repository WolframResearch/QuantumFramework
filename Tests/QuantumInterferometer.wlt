(* Regression tests for QuantumInterferometer: the Reck and Clements meshes must
   rebuild their unitary, and the photon statistics computed from the transfer
   matrix must agree with the Fock-space circuit of BeamSplitterOperator and
   PhaseShiftOperator. *)

Needs["Wolfram`QuantumFramework`"]
Needs["Wolfram`QuantumFramework`SecondQuantization`"]

BeginTestSection["QuantumInterferometer"]

$methods = {"Reck", "Clements", "ClementsPhaseEnd"};

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
        {m, 3, 4}, {method, $methods}
    ],
    ConstantArray[True, {2, 3}],
    TestID -> "QI-exact-Fourier-decomposition"
]

(* The identity that moves the middle phase column to the output (Clements Eq. 5). *)
VerificationTest[
    With[{
        cell = {{Exp[I #2] Cos[#1], - Sin[#1]}, {Exp[I #2] Sin[#1], Cos[#1]}} &,
        phases = DiagonalMatrix[{Exp[I #1], Exp[I #2]}] &
    },
        FullSimplify[
            Inverse[cell[\[Theta], \[Phi]]] . phases[\[Alpha], \[Beta]] ==
                phases[\[Beta] - \[Phi] + Pi, \[Beta]] . cell[\[Theta], \[Alpha] - \[Beta] + Pi]
        ]
    ],
    True,
    TestID -> "QI-phase-end-commutation-identity"
]

VerificationTest[
    SeedRandom[12];
    With[{w = RandomVariate[CircularUnitaryMatrixDistribution[6]]},
        {#["BeamSplitterCount"], #["Depth"]} & /@ (QuantumInterferometer[w, Method -> #] & /@ $methods)
    ],
    {{15, 9}, {15, 6}, {15, 6}},
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
        {method, $methods}
    ],
    ConstantArray[True, {4, 3}],
    TestID -> "QI-degenerate-unitaries"
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
    Keys @ Select[QuantumInterferometer["Fourier"[3]]["Probabilities", {1, 1, 1}], PossibleZeroQ],
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

EndTestSection[]
