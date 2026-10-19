(* QuantumEvolve stores the solver's array-valued InterpolatingFunction whole,
   as a lazy array container, instead of splitting it into one scalar
   interpolation per amplitude.  These tests pin that the container is what is
   stored, that reading it agrees with the split form it replaces, and that
   "Expand" -> True still produces the split form. *)

(* ArrayLazyQ and ArrayDimensions below are Wolfram`Arrays` functions named by
   their short names, so that context has to be on $ContextPath when this file
   is read for those names to resolve rather than intern as inert Global
   symbols.  The runner puts it there for every file; naming it here lets this
   one also run under a bare TestReport. *)
Needs["Wolfram`Arrays`"]

BeginTestSection["QuantumEvolve - lazy container"]

$evolved := QuantumEvolve[
    QuantumOperator["PauliX"], {} -> {},
    QuantumState["0"],
    {\[FormalT], 0, 2}
]

$expanded := QuantumEvolve[
    QuantumOperator["PauliX"], {} -> {},
    QuantumState["0"],
    {\[FormalT], 0, 2},
    "Expand" -> True
]

(* The stored amplitudes are one applied InterpolatingFunction, not an array of
   scalar ones: the container's head-of-head is the whole solver result. *)
VerificationTest[
    Head @ Head @ $evolved["State"],
    InterpolatingFunction,
    TestID -> "Evolve-stores-InterpolatingFunction-whole"
]

VerificationTest[
    ArrayLazyQ @ $evolved["State"],
    True,
    TestID -> "Evolve-container-is-lazy"
]

(* The container reports the state's shape without being materialized, which is
   what lets validity and every shape property answer for it. *)
VerificationTest[
    ArrayDimensions @ $evolved["State"],
    {2},
    TestID -> "Evolve-lazy-container-reports-shape"
]

VerificationTest[
    {$evolved["Dimension"], $evolved["StateType"], $evolved["Qudits"]},
    {2, "Vector", 1},
    TestID -> "Evolve-lazy-state-answers-shape-properties"
]

(* The evolution parameter survives into the state, so the state is callable. *)
VerificationTest[
    {$evolved["ParameterArity"], Length @ $evolved["Parameters"]},
    {1, 1},
    TestID -> "Evolve-lazy-state-keeps-its-parameter"
]

EndTestSection[]


BeginTestSection["QuantumEvolve - lazy and expanded agree"]

(* Binding the parameter evaluates the array-valued function once and gives an
   ordinary explicit state, which is the dominant usage. *)
VerificationTest[
    Head @ $evolved[1.]["State"],
    SparseArray,
    TestID -> "Evolve-bound-state-is-explicit"
]

(* Under a PauliX Hamiltonian from |0>, the populations are cos^2 t and sin^2 t.
   The tolerance is the solver's, not the container's: both routes interpolate
   the same NDSolve grid. *)
VerificationTest[
    Max @ Abs[$evolved[1.]["ProbabilitiesList"] - {Cos[1.] ^ 2, Sin[1.] ^ 2}] < 10 ^ -5,
    True,
    TestID -> "Evolve-lazy-matches-closed-form"
]

VerificationTest[
    Max @ Abs[$expanded[1.]["ProbabilitiesList"] - {Cos[1.] ^ 2, Sin[1.] ^ 2}] < 10 ^ -5,
    True,
    TestID -> "Evolve-expanded-matches-closed-form"
]

(* The two routes agree with each other far more tightly than either agrees
   with the closed form, since they share the solver's grid and differ only in
   how the interpolation is arranged. *)
VerificationTest[
    Max @ Abs[$evolved[1.]["ProbabilitiesList"] - $expanded[1.]["ProbabilitiesList"]] < 10 ^ -6,
    True,
    TestID -> "Evolve-lazy-equals-expanded"
]

VerificationTest[
    Max @ Abs[$evolved[0.3]["ProbabilitiesList"] - $expanded[0.3]["ProbabilitiesList"]] < 10 ^ -6,
    True,
    TestID -> "Evolve-lazy-equals-expanded-off-grid"
]

(* "Expand" -> True restores the pre-container form: an array of scalar
   interpolations rather than one applied InterpolatingFunction. *)
VerificationTest[
    Head @ Head @ $expanded["State"],
    Symbol,
    TestID -> "Evolve-Expand-splits-per-amplitude"
]

VerificationTest[
    ArrayLazyQ @ $expanded["State"],
    False,
    TestID -> "Evolve-Expand-container-is-not-lazy"
]

(* Norm is preserved by unitary evolution on both routes. *)
VerificationTest[
    Abs[$evolved[1.7]["Norm"] - 1] < 10 ^ -5,
    True,
    TestID -> "Evolve-lazy-preserves-norm"
]

EndTestSection[]


BeginTestSection["QuantumEvolve - open-system and operator routes"]

(* A Lindblad evolution produces a density-matrix state; the container has to
   report a square shape for the state to be a valid matrix state. *)
VerificationTest[
    Block[{r = QuantumEvolve[
        QuantumOperator["PauliZ"], {QuantumOperator["PauliX"]} -> {0.4},
        QuantumState["0"],
        {\[FormalT], 0, 1}
    ]},
        {r["StateType"], ArrayDimensions[r["State"]]}
    ],
    {"Matrix", {2, 2}},
    TestID -> "Evolve-Lindblad-lazy-matrix-state"
]

(* A decohering evolution drives the state toward the maximally mixed one, so
   purity falls below 1 and trace is preserved. *)
VerificationTest[
    Block[{r = QuantumEvolve[
        QuantumOperator["PauliZ"], {QuantumOperator["PauliX"]} -> {0.4},
        QuantumState["0"],
        {\[FormalT], 0, 1}
    ][1.]},
        {r["Purity"] < 1, Abs[Total[r["ProbabilitiesList"]] - 1] < 10 ^ -5}
    ],
    {True, True},
    TestID -> "Evolve-Lindblad-decoheres"
]

(* With no state given, the solve produces the propagator as an operator; the
   square-shape branch has to recognize the lazy container as square. *)
VerificationTest[
    Block[{u = QuantumEvolve[QuantumOperator["PauliX"], {} -> {}, None, {\[FormalT], 0, 2}]},
        {Head[u], u["OutputDimension"], u["InputDimension"]}
    ],
    {QuantumOperator, 2, 2},
    TestID -> "Evolve-propagator-is-an-operator"
]

(* The propagator at t reproduces MatrixExp[-I X t] applied to |0>. *)
VerificationTest[
    Block[{u = QuantumEvolve[QuantumOperator["PauliX"], {} -> {}, None, {\[FormalT], 0, 2}]},
        Max @ Abs[
            Normal[(u[1.] @ QuantumState["0"])["ProbabilitiesList"]] - {Cos[1.] ^ 2, Sin[1.] ^ 2}
        ] < 10 ^ -5
    ],
    True,
    TestID -> "Evolve-propagator-matches-closed-form"
]

EndTestSection[]


BeginTestSection["QuantumEvolve - phase-space picture"]

(* Evolving a phase-space state with a symbolic time returns a phase-space state
   carrying the time, not an unevaluated DSolveValue with QuantumEvolve::error.
   The evolution is the linear system W'[t] = R.W[t] with a time-independent rate
   matrix R, so its closed form is the matrix exponential. *)

VerificationTest[
    With[{ev = QuantumEvolve[QuantumOperator["Z"], QuantumPhaseSpaceTransform[QuantumState[{1, 1}/Sqrt[2]]], \[FormalT]]},
        {Head[ev], ev["Picture"]}
    ],
    {QuantumState, "PhaseSpace"},
    TestID -> "Evolve-phase-space-symbolic-time"
]

(* The evolved quasiprobability vector equals exp(R t) W(0) with R the transition
   rate matrix, for real time. *)
VerificationTest[
    With[{
        ev = QuantumEvolve[QuantumOperator["Z"], QuantumPhaseSpaceTransform[QuantumState[{1, 1}/Sqrt[2]]], \[FormalT]],
        rate = Normal @ HamiltonianTransitionRate[QuantumOperator["Z"]],
        w0 = Normal @ QuantumPhaseSpaceTransform[QuantumState[{1, 1}/Sqrt[2]]]["StateVector"]
    },
        Simplify[Normal[ev["StateVector"]] - MatrixExp[\[FormalT] rate] . w0, Element[\[FormalT], Reals]]
    ],
    {0, 0, 0, 0},
    TestID -> "Evolve-phase-space-matches-matrix-exponential"
]

(* A dense state with a symbolic time still evolves, on the unchanged solver path. *)
VerificationTest[
    Head @ QuantumEvolve[QuantumOperator["Z"], QuantumState[{1, 1}/Sqrt[2]], \[FormalT]],
    QuantumState,
    TestID -> "Evolve-dense-symbolic-time-still-works"
]

(* A phase-space state with a time interval still integrates numerically. *)
VerificationTest[
    With[{ev = QuantumEvolve[QuantumOperator["Z"], QuantumPhaseSpaceTransform[QuantumState[{1, 1}/Sqrt[2]]], {\[FormalT], 0, 1}]},
        {Head[ev], ev["Picture"]}
    ],
    {QuantumState, "PhaseSpace"},
    TestID -> "Evolve-phase-space-time-interval-still-works"
]

(* A numerically integrated density matrix carries roundoff that can fail
   PhysicalQ, and so can the state "Physical" repairs it into. Reading its
   probabilities must not ask the repaired state for "Weights" again, which
   recursed until $RecursionLimit on this driven, damped five-level ladder. *)
VerificationTest[
    Module[{b = QuantumOperator[DiagonalMatrix[Sqrt[Range[4]], 1], {1}, 5], h, rho, p},
        h = (-Pi/5) (b["Dagger"] @ b["Dagger"] @ b @ b) + ((Pi/12) (1 - Cos[2 Pi \[FormalT]/12])/2) (b + b["Dagger"]);
        rho = QuantumEvolve[h, {b} -> {1/100000}, QuantumState[UnitVector[5, 1], 5], {\[FormalT], 0, 12}];
        p = rho[12]["ProbabilitiesList"];
        {Length[p], Chop[Total[p] - 1, 1*^-8], 0.96 < p[[2]] < 0.97}
    ],
    {5, 0, True},
    TestID -> "Evolve-Lindblad-probabilities-of-a-roundoff-state-do-not-recurse"
]

(* A jump operator on a subsystem acts as the identity on the rest. Damped Jaynes-Cummings
   with cavity loss: in the one-excitation sector the population of |e,0> is
   e^(-k t/2) (cos w t + k/(4 w) sin w t)^2 with w = sqrt(g^2 - k^2/16). *)
VerificationTest[
    Module[{n = 12, g = 1., k = 0.05, w, a, sm, h, psi0, rho},
        w = Sqrt[g ^ 2 - k ^ 2 / 16];
        a = QuantumOperator[DiagonalMatrix[Sqrt[Range[n - 1]], 1], {2}, n];
        sm = QuantumOperator[{{0, 1}, {0, 0}}, {1}];
        h = g (sm["Dagger"] @ a + sm @ a["Dagger"]);
        psi0 = QuantumState[UnitVector[2 n, n + 1], {2, n}];
        rho = QuantumEvolve[h, {a} -> {k}, psi0, {\[FormalT], 0, 2}];
        Abs[rho[2]["ProbabilitiesList"][[n + 1]] - Exp[-k] (Cos[2 w] + k / (4 w) Sin[2 w]) ^ 2] < 10 ^ -5
    ],
    True,
    TestID -> "Evolve-Lindblad-jump-on-a-subsystem"
]

(* Loss on both modes of a NOON state, with jump operators built on the full space through
   a nine-dimensional identity: the coherence between |n,0> and |0,n> decays as
   e^(-n t) / 2 at unit rates. *)
VerificationTest[
    Module[{d = 9, n = 8, t = -Log[0.9], a, a1, a2, psi, rho},
        a = DiagonalMatrix[Sqrt[Range[d - 1]], 1];
        a1 = QuantumTensorProduct[QuantumOperator[a, {1}, d], QuantumOperator[IdentityMatrix[d], {2}]];
        a2 = QuantumTensorProduct[QuantumOperator[IdentityMatrix[d], {1}], QuantumOperator[a, {2}, d]];
        psi = QuantumState[(UnitVector[d ^ 2, n d + 1] + UnitVector[d ^ 2, n + 1]) / Sqrt[2.], {d, d}];
        rho = QuantumEvolve[0. (a1["Dagger"] @ a1), {a1, a2} -> {1., 1.}, psi, {\[FormalT], 0, t}];
        {a1["Dimensions"], Abs[Abs[Normal[rho[t]["DensityMatrix"]][[n d + 1, n + 1]]] - Exp[-n t] / 2] < 10 ^ -6}
    ],
    {{9, 9, 9, 9}, True},
    TestID -> "Evolve-Lindblad-two-mode-loss-at-a-composite-cutoff"
]

(* Above $QuantumEvolveSparseThreshold unknowns the numeric solve keeps its operators
   sparse, splitting a time-dependent one into numeric arrays with time-dependent
   coefficients; lowering the threshold puts small systems on that route. Under
   H = cos(t) X from |0>, the population of |1> is sin(sin t)^2. *)
VerificationTest[
    Block[{$QuantumEvolveSparseThreshold = 0},
        Abs[QuantumEvolve[QuantumOperator[Cos[\[FormalT]] PauliMatrix[1]], QuantumState["0"], {\[FormalT], 0, 2}][2]["ProbabilitiesList"][[2]] - Sin[Sin[2.]] ^ 2] < 10 ^ -6
    ],
    True,
    TestID -> "Evolve-time-dependent-sparse-matches-closed-form"
]

(* An entry can hold an exact and an inexact multiple of the same time factor, which Plus
   keeps apart; both reach the sparse route. A DRAG-corrected pi pulse on a five-level
   transmon gives the same final state on both routes. *)
VerificationTest[
    Module[{b = QuantumOperator[DiagonalMatrix[Sqrt[Range[4]], 1], {1}, 5], env, h, final},
        env = (Pi / 12) (1 - Cos[2 Pi \[FormalT] / 12]);
        h = (-Pi / 5) (b["Dagger"] @ b["Dagger"] @ b @ b) + 0.03 (b["Dagger"] @ b) + (env / 2) (b + b["Dagger"]) -
            (0.7 D[env, \[FormalT]] / (-4 Pi / 5)) (I (b["Dagger"] - b));
        final[] := Normal @ QuantumEvolve[h, QuantumState[UnitVector[5, 1], 5], {\[FormalT], 0, 12}][12]["StateVector"];
        Max @ Abs[Block[{$QuantumEvolveSparseThreshold = 0}, final[]] - final[]] < 10 ^ -7
    ],
    True,
    TestID -> "Evolve-time-dependent-sparse-keeps-exact-and-inexact-terms"
]

VerificationTest[
    Module[{a = QuantumOperator[DiagonalMatrix[Sqrt[Range[5]], 1], {1}, 6], h, rho, ref},
        h = 0.25 (a["Dagger"] @ a["Dagger"] @ a @ a) + (0.7 Cos[\[FormalT]]) (a + a["Dagger"]);
        rho = Block[{$QuantumEvolveSparseThreshold = 0}, QuantumEvolve[h, {a} -> {0.2}, QuantumState[UnitVector[6, 1], 6], {\[FormalT], 0, 3}]];
        ref = QuantumEvolve[h, {a} -> {0.2}, QuantumState[UnitVector[6, 1], 6], {\[FormalT], 0, 3}, "MergeInterpolatingFunctions" -> False];
        Max @ Abs[Normal[rho[3]["DensityMatrix"]] - Normal[ref[3]["DensityMatrix"]]] < 10 ^ -6
    ],
    True,
    TestID -> "Evolve-Lindblad-time-dependent-sparse-matches-dense-route"
]

EndTestSection[]


BeginTestSection["QuantumEvolve - short pulses"]

(* A pi pulse of width 0.01 at t = 50 in [0, 100], under H = Omega(t) X / 2: step control
   steps over it unless the integration restarts at times read off the time
   dependence. Captured, the population of |1> at t = 100 is 1. *)
$piPulse[omega_] := QuantumEvolve[(omega / 2) QuantumOperator["X"], QuantumState["0"], {\[FormalT], 0, 100}][100]["ProbabilitiesList"][[2]]

VerificationTest[
    Quiet[$piPulse[(Pi / (Sqrt[2 Pi] 0.0025)) Exp[-(\[FormalT] - 50) ^ 2 / (2 0.0025 ^ 2)]], General::munfl] > 1 - 10 ^ -6,
    True,
    TestID -> "Evolve-short-Gaussian-pulse-captured"
]

VerificationTest[
    $piPulse[Piecewise[{{Pi / 0.01, Abs[\[FormalT] - 50] < 0.005}}, 0]] > 1 - 10 ^ -6,
    True,
    TestID -> "Evolve-short-square-pulse-captured"
]

(* Each discontinuity is crossed at its own time, and restarts sit between discontinuities and
   across exponential profiles: nothing for a smooth drive, a crossing for the single edge of a
   step, crossings at both edges of a window with a restart halfway across it, and restarts at
   the centre of a Gaussian and two and four standard deviations either side. *)
VerificationTest[
    Map[
        With[{events = Cases[First @ QuantumEvolve[QuantumOperator[# PauliMatrix[1]], QuantumState["0"], {\[FormalT], 0, 10}, "ReturnEquations" -> True], _WhenEvent]},
            Cases[events, WhenEvent[_ == x_, #] :> Round[x, 10 ^ -6]] & /@ {"RestartIntegration", "CrossDiscontinuity"}
        ] &,
        {Cos[\[FormalT]], UnitStep[\[FormalT] - 3], Piecewise[{{1, 3 < \[FormalT] < 4}}, 0], Exp[-(\[FormalT] - 5) ^ 2 / 2 / 0.01]}
    ],
    {{{}, {}}, {{}, {3}}, {{7 / 2}, {3, 4}}, {{23 / 5, 24 / 5, 5, 26 / 5, 27 / 5}, {}}},
    TestID -> "Evolve-restart-times"
]

(* A solution whose grid repeats the time of a discontinuity, with a different value on
   each side, expands to the same state on both sides of it: just past the pulse edge and
   at the end. *)
VerificationTest[
    Module[{h = (Piecewise[{{Pi / 0.01, Abs[\[FormalT] - 50] < 0.005}}, 0] / 2) QuantumOperator["X"], lazy, expanded},
        lazy = QuantumEvolve[h, QuantumState["0"], {\[FormalT], 0, 100}];
        expanded = QuantumEvolve[h, QuantumState["0"], {\[FormalT], 0, 100}, "Expand" -> True];
        Max @ Abs[Flatten[Table[Normal[expanded[x]["StateVector"]] - Normal[lazy[x]["StateVector"]], {x, {50.006, 50.01, 60., 100.}}]]] < 10 ^ -7
    ],
    True,
    TestID -> "Evolve-Expand-across-a-discontinuity"
]

(* Through the same pulse, a time tau after its first edge and up to its last, the state is
   {Cos[Omega tau / 2], -I Sin[Omega tau / 2]} with Omega = Pi / 0.01, and the expanded state
   follows it on both sides of each edge: NDSolve brackets each edge of this pulse by two grid
   times a rounding apart, and the interpolation runs on either side are kept apart. *)
VerificationTest[
    With[{expanded = QuantumEvolve[(Piecewise[{{Pi / 0.01, Abs[\[FormalT] - 50] < 0.005}}, 0] / 2) QuantumOperator["X"], QuantumState["0"], {\[FormalT], 0, 100}, "Expand" -> True]},
        Max @ Abs[Flatten[Table[
            Normal[expanded[x]["StateVector"]] - With[{tau = Min[x, 50.005] - 49.995}, {Cos[Pi / 0.01 tau / 2], - I Sin[Pi / 0.01 tau / 2]}],
            {x, {49.9950001, 49.9951, 50.0049, 50.0051}}
        ]]] < 10 ^ -6
    ],
    True,
    TestID -> "Evolve-Expand-inside-a-pulse-matches-closed-form"
]

(* With "MergeInterpolatingFunctions" -> False the initial state sits at the start of the
   time range: a Rabi flop over [-1, 1] leaves cos^2 2 in |0>. *)
VerificationTest[
    Abs[QuantumEvolve[QuantumOperator["X"], QuantumState["0"], {\[FormalT], -1, 1}, "MergeInterpolatingFunctions" -> False][1]["ProbabilitiesList"][[1]] - Cos[2.] ^ 2] < 10 ^ -6,
    True,
    TestID -> "Evolve-unmerged-route-starts-at-the-range-start"
]

EndTestSection[]


BeginTestSection["QuantumEvolve - restart times"]

(* Under H = Omega(t) X / 2 from |0>, the drive commutes with itself at all times, so the
   population of |1> after it is Sin[A / 2] ^ 2, with A the area under Omega. A Gaussian narrow
   enough to fit inside the steps NDSolve takes before it is captured only through its restart
   times. At machine precision such a Gaussian underflows far from its centre, which NDSolve
   reports with General::munfl. *)
$excited[omega_, {t0_, t1_}, opts___] := QuantumEvolve[(omega / 2) QuantumOperator["X"], QuantumState["0"], {\[FormalT], t0, t1}, opts][t1]["ProbabilitiesList"][[2]]

$gaussian[area_, centre_, sigma_] := area / (Sqrt[2 Pi] sigma) Exp[-(\[FormalT] - centre) ^ 2 / (2 sigma ^ 2)]

(* A pi/2 square pulse from the start of the range, written as a product of step functions of
   which Boole[t >= 0] has no edge inside the range, then a Gaussian pi/2 pulse half the range
   later: A = Pi. The dense, sparse and unmerged routes solve the same equations; NDSolve cannot
   see into the right-hand side of the unmerged one. *)
$squareThenGaussian = (Pi / 2) Boole[\[FormalT] >= 0] Boole[\[FormalT] <= 1] + $gaussian[Pi / 2, 50, 0.025];

VerificationTest[
    Abs[$excited[$squareThenGaussian, {0, 100}] - 1] < 10 ^ -6,
    True,
    {General::munfl},
    TestID -> "Evolve-square-from-range-start-then-Gaussian"
]

VerificationTest[
    Abs[Block[{$QuantumEvolveSparseThreshold = 0}, $excited[$squareThenGaussian, {0, 100}]] - 1] < 10 ^ -6,
    True,
    {General::munfl},
    TestID -> "Evolve-square-from-range-start-then-Gaussian-sparse"
]

VerificationTest[
    Abs[$excited[$squareThenGaussian, {0, 100}, "MergeInterpolatingFunctions" -> False] - 1] < 10 ^ -6,
    True,
    {General::munfl},
    TestID -> "Evolve-square-from-range-start-then-Gaussian-unmerged"
]

VerificationTest[
    Abs[$excited[(Pi / 2) UnitStep[\[FormalT]] UnitStep[1 - \[FormalT]] + $gaussian[Pi / 2, 50, 0.025], {0, 100}] - 1] < 10 ^ -6,
    True,
    {General::munfl},
    TestID -> "Evolve-UnitStep-square-from-range-start-then-Gaussian"
]

(* A Gaussian pi/2 pulse, then a pi/2 square pulse that ends where the range ends: A = Pi. *)
VerificationTest[
    Abs[$excited[$gaussian[Pi / 2, 50, 0.025] + (Pi / 2) Boole[\[FormalT] >= 99] Boole[\[FormalT] <= 100], {0, 100}] - 1] < 10 ^ -6,
    True,
    {General::munfl},
    TestID -> "Evolve-square-ending-at-range-end"
]

(* A square pulse whose edges are the two ends of the range, pi/2 spread over all of it, under a
   Gaussian pi/2 pulse at its middle: A = Pi. *)
VerificationTest[
    Abs[$excited[(Pi / 200) Boole[\[FormalT] >= 0] Boole[\[FormalT] <= 100] + $gaussian[Pi / 2, 50, 0.025], {0, 100}] - 1] < 10 ^ -6,
    True,
    {General::munfl},
    TestID -> "Evolve-square-over-the-whole-range"
]

(* The two raised-cosine pi/2 windows of width 20 of QuTiP issue #1265, each a product of two
   Boole factors, the first from the start of the range, then a Gaussian pi/4 pulse after a long
   wait, evolved as a density matrix: A = 5 Pi / 4. *)
VerificationTest[
    With[{window = Function[centre, With[{x = \[FormalT] - centre}, (Pi / 2) / 20 (1 + Cos[2 Pi x / 20]) Boole[x + 20 / 2. > 0] Boole[- x + 20 / 2. > 0]]]},
        Abs[
            QuantumEvolve[((window[10.] + window[70.] + $gaussian[Pi / 4, 300, 0.05]) / 2) QuantumOperator["X"], QuantumState["0"]["MatrixState"], {\[FormalT], 0, 400}][400]["ProbabilitiesList"][[2]] -
                Sin[5 Pi / 8] ^ 2
        ] < 10 ^ -6
    ],
    True,
    {General::munfl},
    TestID -> "Evolve-raised-cosine-windows-then-Gaussian"
]

(* A drive Omega E^(I delta t) detuned by delta, on a narrow window written as one condition,
   whose two edges are the two roots of one switching function. In the frame of the drive the
   window is a constant Rabi problem: P1 = Omega^2 / (Omega^2 + delta^2) Sin[Sqrt[Omega^2 + delta^2] w / 2] ^ 2. *)
VerificationTest[
    With[{w = 0.001, omega = Pi / 0.001, delta = Pi / 0.001},
        With[{window = Abs[\[FormalT] - 50] < w / 2},
            Abs[
                QuantumEvolve[
                    QuantumOperator[{{0, Piecewise[{{omega Exp[I delta \[FormalT]], window}}, 0]}, {Piecewise[{{omega Exp[- I delta \[FormalT]], window}}, 0], 0}} / 2],
                    QuantumState["0"],
                    {\[FormalT], 0, 100}
                ][100]["ProbabilitiesList"][[2]] -
                    omega ^ 2 / (omega ^ 2 + delta ^ 2) Sin[Sqrt[omega ^ 2 + delta ^ 2] w / 2] ^ 2
            ] < 10 ^ -6
        ]
    ],
    True,
    TestID -> "Evolve-detuned-narrow-window"
]

(* A pi pulse of 10 ns with the time in seconds, and a pi pulse of width 0.01 a million time
   units from the time origin: each window is short against its distance from the origin, and
   P1 = 1 after it. *)
VerificationTest[
    Abs[$excited[(Pi / 1*^-8) Boole[Abs[\[FormalT] - 5*^-5] < 0.5*^-8], {0, 1*^-4}] - 1] < 10 ^ -6,
    True,
    TestID -> "Evolve-nanosecond-pulse-in-seconds"
]

VerificationTest[
    Abs[$excited[(Pi / 0.01) Boole[Abs[\[FormalT] - (10 ^ 6 + 50)] < 0.005], {10 ^ 6, 10 ^ 6 + 100}] - 1] < 10 ^ -6,
    True,
    TestID -> "Evolve-pulse-far-from-the-time-origin"
]

(* A three-level transmon in the frame rotating at its qubit frequency, anharmonicity
   alpha / 2 Pi = -200 MHz with time in ns, from its ground state: a pi pulse windowed on [0, 16]
   as a product of Boole factors, a 1000 ns idle, then a Gaussian pi pulse of sigma 4 ns. The
   reference integrates each pulse alone, at WorkingPrecision 40 and 30, and carries the state
   across the idle by its exact phases; it agrees with a WorkingPrecision 50 solve to 1e-15. The
   leaked |2> amplitude keeps turning through the idle, and NDSolve's default tolerances leave an
   error near 1e-6 over the 1100 ns, against an error of order one if the second pulse is
   stepped over. *)
VerificationTest[
    Module[{alpha = -2 Pi 2 / 10, sigma = 4, window = 16, idle = 1000, centre, end, a, h0, hd, first, second, pulse, reference},
        a = {{0, 1, 0}, {0, 0, Sqrt[2]}, {0, 0, 0}};
        h0 = alpha / 2 ConjugateTranspose[a] . ConjugateTranspose[a] . a . a;
        hd = (a + ConjugateTranspose[a]) / 2;
        centre = window + idle + 48;
        end = centre + 48;
        first = Pi / (sigma Sqrt[2 Pi] Erf[window / (2 Sqrt[2] sigma)]) Exp[-(\[FormalT] - window / 2) ^ 2 / (2 sigma ^ 2)];
        second = Pi / (Sqrt[2 Pi] sigma) Exp[-(\[FormalT] - centre) ^ 2 / (2 sigma ^ 2)];
        pulse[drive_, psi0_, {ta_, tb_}, wp_] := Block[{psi},
            NDSolveValue[{psi'[\[FormalT]] == -I (h0 + drive hd) . psi[\[FormalT]], psi[ta] == psi0}, psi[tb], {\[FormalT], ta, tb},
                WorkingPrecision -> wp, PrecisionGoal -> wp / 2, AccuracyGoal -> wp / 2, MaxStepSize -> 1 / 20, MaxSteps -> Infinity
            ]
        ];
        reference = Abs[pulse[second, Exp[-I Diagonal[h0] idle] pulse[first, {1, 0, 0}, {0, window}, 40], {window + idle, end}, 30]] ^ 2;
        Max @ Abs[
            QuantumEvolve[QuantumOperator[h0 + (first Boole[\[FormalT] >= 0] Boole[\[FormalT] <= window] + second) hd], QuantumState[{1, 0, 0}], {\[FormalT], 0, end}][end]["ProbabilitiesList"] -
                reference
        ] < 10 ^ -6
    ],
    True,
    TestID -> "Evolve-transmon-window-idle-Gaussian"
]

(* A Gaussian pi pulse of sigma 0.01 given as samples 0.0025 apart, beside a constant energy
   offset that only shifts the phase. No restart times are read off data, and its samples lie
   closer together than the steps NDSolve may take, which QuantumEvolve says. With MaxStepSize at
   the sample spacing the pulse is captured: P1 = Sin[A / 2] ^ 2, with A the trapezoidal area of
   the samples, which their linear interpolation integrates exactly. *)
VerificationTest[
    With[{samples = Join[{{0., 0.}}, Table[{x, Pi / (Sqrt[2 Pi] 0.01) Exp[-(x - 10) ^ 2 / (2 0.01 ^ 2)]}, {x, 9.95, 10.05, 0.0025}], {{20., 0.}}]},
        With[{
            h = (Interpolation[samples, InterpolationOrder -> 1][\[FormalT]] / 2) QuantumOperator["X"] + QuantumOperator["I"],
            area = Total[Differences[samples[[All, 1]]] MovingAverage[samples[[All, 2]], 2]]
        },
            QuantumEvolve[h, QuantumState["0"], {\[FormalT], 0, 20}, "ReturnEquations" -> True];
            Abs[QuantumEvolve[h, QuantumState["0"], {\[FormalT], 0, 20}, MaxStepSize -> 0.0025][20]["ProbabilitiesList"][[2]] - Sin[area / 2] ^ 2] < 10 ^ -6
        ]
    ],
    True,
    {QuantumEvolve::pulsedata},
    TestID -> "Evolve-sampled-pulse-names-MaxStepSize"
]

(* Floor[1 / (t - 50)] jumps at 50 + 1 / k for every integer k, times that accumulate at 50, so
   its discontinuities are not a list of times; QuantumEvolve says so, and crosses only t = 50,
   the one time Reduce lists. *)
VerificationTest[
    Cases[
        First @ QuantumEvolve[(Floor[1 / (\[FormalT] - 50)] / 2) QuantumOperator["X"], QuantumState["0"], {\[FormalT], 0, 100}, "ReturnEquations" -> True],
        WhenEvent[_ == x_, action_] :> {Round[x, 10 ^ -6], action}
    ],
    {{50, "CrossDiscontinuity"}},
    {QuantumEvolve::pulseedges},
    TestID -> "Evolve-accumulating-discontinuities-are-named"
]

(* A smooth drive has no discontinuity and no exponential profile, so no restart times: under
   Omega = 2 Cos[3 t] the population of |1> at t = 10 is Sin[Sin[30] / 3] ^ 2. *)
VerificationTest[
    {
        Abs[$excited[2 Cos[3 \[FormalT]], {0, 10}] - Sin[Sin[30.] / 3] ^ 2] < 10 ^ -6,
        Cases[First @ QuantumEvolve[Cos[3 \[FormalT]] QuantumOperator["X"], QuantumState["0"], {\[FormalT], 0, 10}, "ReturnEquations" -> True], _WhenEvent]
    },
    {True, {}},
    TestID -> "Evolve-smooth-drive-has-no-restart-times"
]

EndTestSection[]


BeginTestSection["QuantumEvolve - malformed input"]

(* A call that fits none of QuantumEvolve's forms returns a Failure at once, with a
   message naming the argument at fault. TimeConstrained turns a call that rewrites
   itself without end into a failed test instead of a stalled run. First of a Failure
   is its tag. *)

$sigmaMinus = QuantumOperator[{{0, 1}, {0, 0}}];

VerificationTest[
    With[{f = TimeConstrained[QuantumEvolve[QuantumOperator["Z"], {foo[2]} -> {1}, QuantumState["1"], \[FormalT]], 10, {"TIMEOUT"}]},
        {First[f], f["MessageParameters"]}
    ],
    {"InvalidJumpOperators", {{foo[2]}}},
    {QuantumEvolve::jump},
    TestID -> "Evolve-nonoperator-jump-fails"
]

VerificationTest[
    First @ TimeConstrained[QuantumEvolve[QuantumOperator["Z"], foo[2] -> 1, QuantumState["1"], \[FormalT]], 10, {"TIMEOUT"}],
    "InvalidJumpOperators",
    {QuantumEvolve::jump},
    TestID -> "Evolve-nonoperator-jump-rule-fails"
]

VerificationTest[
    First @ TimeConstrained[QuantumEvolve[QuantumOperator["Z"], {$sigmaMinus, foo} -> {1, 1}, QuantumState["1"], \[FormalT]], 10, {"TIMEOUT"}],
    "InvalidJumpOperators",
    {QuantumEvolve::jump},
    TestID -> "Evolve-mixed-jump-list-fails"
]

(* Through the observable form the Failure comes back as it is. *)
VerificationTest[
    First @ TimeConstrained[QuantumEvolve[QuantumOperator["Z"], {foo} -> {1}, QuantumOperator["X"], \[FormalT]], 10, {"TIMEOUT"}],
    "InvalidJumpOperators",
    {QuantumEvolve::jump},
    TestID -> "Evolve-observable-form-returns-the-failure"
]

VerificationTest[
    First @ TimeConstrained[QuantumEvolve[QuantumOperator["Z"], QuantumState["1"], {\[FormalT], 0, T}], 10, {"TIMEOUT"}],
    "InvalidArguments",
    {QuantumEvolve::args},
    TestID -> "Evolve-symbolic-time-range-fails"
]

VerificationTest[
    First @ TimeConstrained[QuantumEvolve[QuantumOperator["Z"], {$sigmaMinus} -> {1}, "1", \[FormalT]], 10, {"TIMEOUT"}],
    "InvalidArguments",
    {QuantumEvolve::args},
    TestID -> "Evolve-string-state-fails"
]

(* Rates go in a rule with the jump operators, never as a separate argument. *)
VerificationTest[
    First @ TimeConstrained[QuantumEvolve[QuantumOperator["Z"], {$sigmaMinus}, {1}, QuantumState["1"], \[FormalT]], 10, {"TIMEOUT"}],
    "InvalidArguments",
    {QuantumEvolve::args},
    TestID -> "Evolve-positional-rates-fail"
]

VerificationTest[
    First @ TimeConstrained[QuantumEvolve[QuantumOperator["Z"], {$sigmaMinus} -> {1}, QuantumState["1"], \[FormalT], 5], 10, {"TIMEOUT"}],
    "InvalidArguments",
    {QuantumEvolve::args},
    TestID -> "Evolve-extra-argument-fails"
]

(* A rate list of the wrong length is refused, neither padded with unit rates nor cut. *)
VerificationTest[
    With[{f = QuantumEvolve[QuantumOperator["Z"], {$sigmaMinus, $sigmaMinus} -> {\[FormalG]}, QuantumState["1"], \[FormalT]]},
        {First[f], f["MessageParameters"]}
    ],
    {"InvalidRates", {{\[FormalG]}, 2}},
    {QuantumEvolve::rates},
    TestID -> "Evolve-too-few-rates-fail"
]

VerificationTest[
    First @ QuantumEvolve[QuantumOperator["Z"], {$sigmaMinus} -> {1, 2}, QuantumState["1"], \[FormalT]],
    "InvalidRates",
    {QuantumEvolve::rates},
    TestID -> "Evolve-too-many-rates-fail"
]

VerificationTest[
    First @ QuantumEvolve[QuantumOperator["Z"], {$sigmaMinus} -> {{1, 0}, {0, 1}}, QuantumState["1"], \[FormalT]],
    "InvalidRates",
    {QuantumEvolve::rates},
    TestID -> "Evolve-rate-matrix-of-the-wrong-size-fails"
]

(* An option in the second argument is an option, not a jump operator. *)
VerificationTest[
    Length @ TimeConstrained[QuantumEvolve[QuantumOperator["Z"], "ReturnEquations" -> True], 10, {"TIMEOUT"}],
    3,
    TestID -> "Evolve-options-without-a-state"
]

EndTestSection[]


BeginTestSection["QuantumEvolve - Lindblad forms"]

(* Amplitude damping of |1> by sigma-minus at rate g: the excited population decays as
   Exp[-g t] whether the rate is a rule, folded into the jump operator, or given as a
   1 by 1 rate matrix. *)
$dampingClosedForm = {{1 - Exp[-\[FormalG] \[FormalT]], 0}, {0, Exp[-\[FormalG] \[FormalT]]}};

VerificationTest[
    Simplify[
        Normal[QuantumEvolve[QuantumOperator["Z"], {$sigmaMinus} -> {\[FormalG]}, QuantumState["1"], \[FormalT]]["DensityMatrix"]] - $dampingClosedForm,
        \[FormalG] > 0
    ],
    {{0, 0}, {0, 0}},
    TestID -> "Evolve-damping-closed-form"
]

VerificationTest[
    Simplify[
        {
            Normal[QuantumEvolve[QuantumOperator["Z"], {Sqrt[\[FormalG]] $sigmaMinus}, QuantumState["1"], \[FormalT]]["DensityMatrix"]],
            Normal[QuantumEvolve[QuantumOperator["Z"], $sigmaMinus -> \[FormalG], QuantumState["1"], \[FormalT]]["DensityMatrix"]],
            Normal[QuantumEvolve[QuantumOperator["Z"], {$sigmaMinus} -> {{\[FormalG]}}, QuantumState["1"], \[FormalT]]["DensityMatrix"]]
        } - {$dampingClosedForm, $dampingClosedForm, $dampingClosedForm},
        \[FormalG] > 0
    ],
    ConstantArray[0, {3, 2, 2}],
    TestID -> "Evolve-damping-rate-forms-agree"
]

(* No rates means a unit rate for each jump operator. *)
VerificationTest[
    Normal[QuantumEvolve[QuantumOperator["Z"], {$sigmaMinus} -> {}, QuantumState["1"], \[FormalT]]["DensityMatrix"]] ===
        Normal[QuantumEvolve[QuantumOperator["Z"], {$sigmaMinus}, QuantumState["1"], \[FormalT]]["DensityMatrix"]] ===
        Normal[QuantumEvolve[QuantumOperator["Z"], {$sigmaMinus} -> {1}, QuantumState["1"], \[FormalT]]["DensityMatrix"]],
    True,
    TestID -> "Evolve-no-rates-are-unit-rates"
]

(* The damped state relaxes to the ground state |0>. *)
VerificationTest[
    Limit[
        Simplify[Normal[QuantumEvolve[QuantumOperator["Z"], {$sigmaMinus} -> {\[FormalG]}, QuantumState["1"], \[FormalT]]["DensityMatrix"]], \[FormalG] > 0],
        \[FormalT] -> Infinity,
        Assumptions -> \[FormalG] > 0
    ],
    {{1, 0}, {0, 0}},
    TestID -> "Evolve-damping-relaxes-to-the-ground-state"
]

(* A numeric time range reproduces the closed form. *)
VerificationTest[
    Max @ Abs[
        Normal[QuantumEvolve[QuantumOperator["Z"], {$sigmaMinus} -> {1 / 2}, QuantumState["1"], {\[FormalT], 0, 3}][3]["DensityMatrix"]] -
            N[$dampingClosedForm /. {\[FormalG] -> 1 / 2, \[FormalT] -> 3}]
    ] < 10 ^ -6,
    True,
    TestID -> "Evolve-damping-numeric-time-matches-closed-form"
]

(* Without jump operators a symbolic time gives the unitary closed form. *)
VerificationTest[
    Simplify[Normal[QuantumEvolve[QuantumOperator["X"], QuantumState["0"], \[FormalT]]["StateVector"]] - {Cos[\[FormalT]], -I Sin[\[FormalT]]}],
    {0, 0},
    TestID -> "Evolve-symbolic-time-unitary"
]

EndTestSection[]
