(* QuantumEvolve stores the solver's array-valued InterpolatingFunction whole,
   as a lazy array container, instead of splitting it into one scalar
   interpolation per amplitude.  These tests pin that the container is what is
   stored, that reading it agrees with the split form it replaces, and that
   "Expand" -> True still produces the split form. *)

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
