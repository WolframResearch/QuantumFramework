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
   steps over it unless the integration restarts at the landmarks of the time
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

(* The landmarks are restarts at each discontinuity and across each Gaussian profile; a
   smooth drive has none. *)
VerificationTest[
    {
        Count[First @ QuantumEvolve[QuantumOperator[Cos[\[FormalT]] PauliMatrix[1]], QuantumState["0"], {\[FormalT], 0, 10}, "ReturnEquations" -> True], _WhenEvent],
        Length @ Cases[
            First @ QuantumEvolve[QuantumOperator[Exp[-(\[FormalT] - 5) ^ 2 / 2 / 0.01] PauliMatrix[1]], QuantumState["0"], {\[FormalT], 0, 10}, "ReturnEquations" -> True],
            WhenEvent[_ == x_, "RestartIntegration"] :> Round[x, 10 ^ -6]
        ],
        Cases[
            First @ QuantumEvolve[QuantumOperator[UnitStep[\[FormalT] - 3] PauliMatrix[1]], QuantumState["0"], {\[FormalT], 0, 10}, "ReturnEquations" -> True],
            WhenEvent[_ == x_, "RestartIntegration"] :> Round[x, 10 ^ -6]
        ]
    },
    {0, 5, {3}},
    TestID -> "Evolve-restart-landmarks"
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

(* With "MergeInterpolatingFunctions" -> False the initial state sits at the start of the
   time range: a Rabi flop over [-1, 1] leaves cos^2 2 in |0>. *)
VerificationTest[
    Abs[QuantumEvolve[QuantumOperator["X"], QuantumState["0"], {\[FormalT], -1, 1}, "MergeInterpolatingFunctions" -> False][1]["ProbabilitiesList"][[1]] - Cos[2.] ^ 2] < 10 ^ -6,
    True,
    TestID -> "Evolve-unmerged-route-starts-at-the-range-start"
]

EndTestSection[]
