(* ==========================================================================
   Tests/HybridInterop.wlt -- Phase 7.1 hybrid interop UpValues.

   Cross-head dispatch: QuantumMeasurementOperator / QuantumChannel applied to
   PauliStabilizer / StabilizerFrame. Pauli-basis QMOs stay in the tableau
   (route to ps["M", pauli]); non-Pauli bases emit ::nonpaulibasis info and
   fall back to legacy/dense paths.

   Phase 7.1 covers:
     - Pauli-string QMO on a PauliStabilizer (Pauli, Bell, 5Q-code stabilizers).
     - Non-Pauli QMO -> info message + fallback (smoke test).
     - QuantumChannel on PauliStabilizer / StabilizerFrame -> info + fallback.
   ========================================================================== *)

Needs["Wolfram`QuantumFramework`"];

(* Local validity helper: re-export the package-scoped predicate so tests can
   assert on the receiver's structure. *)
psValidQ = Wolfram`QuantumFramework`PackageScope`PauliStabilizerQ;


(* ============================================================================ *)
(* TIER A -- Pauli-basis fast path: qmo[ps] -> ps["M", pauliString]            *)
(* ============================================================================ *)

(* Two-qubit Z-on-q1 measurement on |00>: deterministic outcome 0.
   Note: a 1-qubit Z measurement currently exposes ROADMAP A.11 (the single-
   element KroneckerProduct in pauliStringMatrix), so we test on n=2. *)
VerificationTest[
    With[{ps = PauliStabilizer[2], qmo = QuantumMeasurementOperator["ZI", {1, 2}]},
        Sort @ Keys @ qmo[ps]
    ],
    {0},
    {},
    TestID -> "Phase7-QMO-ZI-on-zero-Deterministic"
]

(* Single-qubit Z-measurement on |+> is non-deterministic (both outcomes). *)
VerificationTest[
    With[{ps = PauliStabilizer[1]["H", 1], qmo = QuantumMeasurementOperator["Z", {1}]},
        Sort @ Keys @ qmo[ps]
    ],
    {0, 1},
    {},
    TestID -> "Phase7-QMO-Z-on-Plus-NonDeterministic"
]

(* Two-qubit Pauli-string measurement: Bell ZZ gives deterministic outcome 0. *)
VerificationTest[
    With[{psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]],
          qmo = QuantumMeasurementOperator["ZZ", {1, 2}]},
        Sort @ Keys @ qmo[psBell]
    ],
    {0},
    {},
    TestID -> "Phase7-QMO-ZZ-on-Bell-Deterministic"
]

(* 5Q-code stabilizer measurement gives deterministic outcome 0. *)
VerificationTest[
    With[{ps5 = PauliStabilizer["5QubitCode"],
          qmo = QuantumMeasurementOperator["XZZXI", Range[5]]},
        Sort @ Keys @ qmo[ps5]
    ],
    {0},
    {},
    TestID -> "Phase7-QMO-XZZXI-on-5Q-Deterministic"
]

(* Pauli-basis QMO returns the same Association as ps["M", pauliString] -- this
   is the direct equivalence test that distinguishes the fast path from the
   fallback. *)
VerificationTest[
    With[{psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]],
          qmo = QuantumMeasurementOperator["XX", {1, 2}]},
        Keys @ qmo[psBell] === Keys @ psBell["M", "XX"]
    ],
    True,
    {},
    TestID -> "Phase7-QMO-XX-on-Bell-MatchesNativeM"
]


(* ============================================================================ *)
(* TIER B -- Non-Pauli basis: info message + dense fallback                    *)
(* ============================================================================ *)

(* A computational-basis projector measurement (label not Pauli) triggers the
   ::nonpaulibasis fallback. The Phase 7.1 fallback routes through the legacy
   PauliStabilizerApply path which keeps the result as a PauliStabilizer (and
   may itself emit nonclifford for unrecognized gates). Smoke check: message is
   emitted, the call doesn't crash. *)
VerificationTest[
    With[{ps = PauliStabilizer[1],
          qmo = QuantumMeasurementOperator[QuantumBasis["Computational"], {1}]},
        Head @ qmo[ps]
    ],
    PauliStabilizer,
    {PauliStabilizer::nonpaulibasis, PauliStabilizer::nonclifford},
    TestID -> "Phase7-QMO-NonPauliBasis-Fallback-EmitsMessage"
]


(* ============================================================================ *)
(* TIER C -- QuantumChannel on stabilizer inputs (Phase 7.2 routes named       *)
(*           Pauli channels through tableau)                                   *)
(* ============================================================================ *)

(* Phase 7.2: BitFlip[p] on a stabilizer state is a named Pauli channel and   *)
(* now stays in the tableau, returning a list {{prob, ps_branch}, ...}.       *)

VerificationTest[
    With[{ps = PauliStabilizer[1],
          qc = QuantumChannel["BitFlip"[1/3], {1}]},
        (* Two branches: (2/3, no-op identity), (1/3, X applied) *)
        Length @ qc[ps]
    ],
    2,
    {},
    TestID -> "Phase7.2-QC-BitFlip-NumBranches"
]

VerificationTest[
    With[{ps = PauliStabilizer[1],
          qc = QuantumChannel["BitFlip"[1/3], {1}]},
        (* The probabilities sum to 1 *)
        Total[First /@ qc[ps]]
    ],
    1,
    {},
    TestID -> "Phase7.2-QC-BitFlip-ProbabilitiesSumToOne"
]

VerificationTest[
    With[{ps = PauliStabilizer[1],
          qc = QuantumChannel["BitFlip"[1/3], {1}]},
        (* Identity branch keeps ps; X branch applies X to |0>, giving |1> stab Z->-Z *)
        SameQ[
            qc[ps][[1, 2]]["Stabilizers"],
            ps["Stabilizers"]
        ]
    ],
    True,
    {},
    TestID -> "Phase7.2-QC-BitFlip-IdentityBranchUnchanged"
]

VerificationTest[
    With[{ps = PauliStabilizer[1],
          qc = QuantumChannel["PhaseFlip"[1/4], {1}]},
        (* PhaseFlip applies Z; Z on |0> flips no sign (Z|0> = |0>). Identity branch
           and Z branch produce the same state on |0> -- both leave stabilizer Z. *)
        Length @ qc[ps]
    ],
    2,
    {},
    TestID -> "Phase7.2-QC-PhaseFlip-NumBranches"
]

VerificationTest[
    With[{ps = PauliStabilizer[1]["H", 1],
          qc = QuantumChannel["BitPhaseFlip"[1/5], {1}]},
        (* BitPhaseFlip applies Y. Two branches. *)
        Length @ qc[ps]
    ],
    2,
    {},
    TestID -> "Phase7.2-QC-BitPhaseFlip-NumBranches"
]

VerificationTest[
    With[{ps = PauliStabilizer[1],
          qc = QuantumChannel["Depolarizing"[1/2], {1}]},
        (* Depolarizing[p] has 4 branches: I, X, Y, Z. *)
        Length @ qc[ps]
    ],
    4,
    {},
    TestID -> "Phase7.2-QC-Depolarizing-FourBranches"
]

VerificationTest[
    With[{ps = PauliStabilizer[1],
          qc = QuantumChannel["Depolarizing"[1/2], {1}]},
        Total[First /@ qc[ps]]
    ],
    1,
    {},
    TestID -> "Phase7.2-QC-Depolarizing-ProbabilitiesSumToOne"
]

(* Non-Clifford channel (AmplitudeDamping) still falls back. *)
VerificationTest[
    With[{ps = PauliStabilizer[1],
          qc = QuantumChannel["AmplitudeDamping"[1/2], {1}]},
        Head @ qc[ps]
    ],
    QuantumState,
    {PauliStabilizer::nonpaulibasis},
    TestID -> "Phase7.2-QC-AmplitudeDamping-FallbackMaterializes"
]


(* ============================================================================ *)
(* TIER D -- ConcretePauliStabilizerQ guard: symbolic-phase ps falls through   *)
(* ============================================================================ *)

(* After ps["SymbolicMeasure", q] the receiver has symbolic phases. The fast
   path is gated on ConcretePauliStabilizerQ (numeric +-1 signs only); a
   subsequent qmo[ps] on the symbolic ps must NOT fire the Pauli fast path. *)
VerificationTest[
    Module[{psSym, qmo, h},
        psSym = PauliStabilizer[1]["H", 1]["SymbolicMeasure", 1];
        qmo = QuantumMeasurementOperator["Z", {1}];
        h = Head @ qmo[psSym];
        (* The dispatch should NOT match the Pauli fast path (which requires
           ConcretePauliStabilizerQ). It falls to whatever generic dispatch
           does -- we only require it doesn't crash and stays evaluable. *)
        FreeQ[h, $Failed]
    ],
    True,
    TestID -> "Phase7-QMO-on-SymbolicPS-FallsThrough"
]


(* ============================================================================ *)
(* TIER E -- Pauli-string QMO on named QEC code states                         *)
(* ============================================================================ *)

(* Steane-7Q stabilizer measurement: deterministic outcome 0 on |0_L>.
   Use the actual stabilizers from the kernel (avoids hand-coding the X/Z
   layout). Pick a known stabilizer row (e.g. the last "XXXXXXX"). *)
VerificationTest[
    With[{ps = PauliStabilizer["SteaneCode"]},
        With[{stab = ps["Stabilizers"][[7]],   (* "XXXXXXX" by API.md *)
              qmo = QuantumMeasurementOperator[ps["Stabilizers"][[7]], Range[7]]},
            Sort @ Keys @ qmo[ps]
        ]
    ],
    {0},
    {},
    TestID -> "Phase7-QMO-Steane-Stabilizer-Deterministic"
]

(* GHZ-3 measured in the ZZI Pauli basis: ZZI is in the stabilizer group <ZZI, IZZ, XXX>?
   For GHZ-3 with stabilizer {XXX, ZZI, IZZ}, the operator ZZI is one of the generators
   so the measurement is deterministic with outcome 0. *)
VerificationTest[
    With[{psGHZ = PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "CNOT" -> {1, 2}, "CNOT" -> {2, 3}},
          qmo = QuantumMeasurementOperator["ZZI", {1, 2, 3}]},
        Sort @ Keys @ qmo[psGHZ]
    ],
    {0},
    {},
    TestID -> "Phase7-QMO-GHZ3-Stab-ZZI-Deterministic"
]

(* GHZ-3 measured against XXX: deterministic outcome 0 (XXX is a stabilizer). *)
VerificationTest[
    With[{psGHZ = PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "CNOT" -> {1, 2}, "CNOT" -> {2, 3}},
          qmo = QuantumMeasurementOperator["XXX", {1, 2, 3}]},
        Sort @ Keys @ qmo[psGHZ]
    ],
    {0},
    {},
    TestID -> "Phase7-QMO-GHZ3-Stab-XXX-Deterministic"
]

(* GHZ-3 measured against an anticommuting Pauli (XII anticommutes with ZZI):
   outcome must be non-deterministic (both 0 and 1 keys present). *)
VerificationTest[
    With[{psGHZ = PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "CNOT" -> {1, 2}, "CNOT" -> {2, 3}},
          qmo = QuantumMeasurementOperator["XII", {1, 2, 3}]},
        Sort @ Keys @ qmo[psGHZ]
    ],
    {0, 1},
    {},
    TestID -> "Phase7-QMO-GHZ3-Anticommuting-NonDeterministic"
]


(* ============================================================================ *)
(* TIER F -- Channel edge cases                                                 *)
(* ============================================================================ *)

(* BitFlip[0] is a no-op identity channel: all probability on the I branch.    *)
VerificationTest[
    With[{ps = PauliStabilizer[1], qc = QuantumChannel["BitFlip"[0], {1}]},
        First /@ qc[ps]                  (* probabilities *)
    ],
    {1, 0},
    {},
    TestID -> "Phase7.2-QC-BitFlip-zero-IdentityNoOp"
]

(* BitFlip[1] applies X with probability 1: I branch has prob 0. *)
VerificationTest[
    With[{ps = PauliStabilizer[1], qc = QuantumChannel["BitFlip"[1], {1}]},
        First /@ qc[ps]
    ],
    {0, 1},
    {},
    TestID -> "Phase7.2-QC-BitFlip-one-FullX"
]

(* Symbolic probability: BitFlip[p] returns symbolic branches. *)
VerificationTest[
    With[{ps = PauliStabilizer[1], qc = QuantumChannel["BitFlip"[\[FormalP]], {1}]},
        Total[First /@ qc[ps]] // Simplify
    ],
    1,
    {},
    TestID -> "Phase7.2-QC-BitFlip-Symbolic-ProbSumToOne"
]

(* PhaseFlip[p] on |0>: Z branch's post-state is Z|0> stabilizer = "Z" (sign-
   flipped is "-Z" only if Z was a sign-anticommuting stab; for |0>, Z|0> = |0>
   so the stabilizer string stays "Z"). The branch shape is {prob, ps_after}. *)
VerificationTest[
    With[{ps = PauliStabilizer[1], qc = QuantumChannel["PhaseFlip"[1/3], {1}]},
        First @ qc[ps][[2, 2]]["Stabilizers"]   (* second branch -> ps -> Stabilizers -> first elt *)
    ],
    "Z",
    {},
    TestID -> "Phase7.2-QC-PhaseFlip-Z-on-zero-StaysZ"
]

(* Channel on a multi-qubit state targeting one of the qubits.                 *)
VerificationTest[
    With[{psBell = PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "CNOT" -> {1, 2}},
          qc = QuantumChannel["BitFlip"[1/4], {1}]},
        Length @ qc[psBell]
    ],
    2,
    {},
    TestID -> "Phase7.2-QC-BitFlip-on-Bell-Q1-NumBranches"
]

(* Depolarizing branches: identity branch state matches input. *)
VerificationTest[
    With[{ps = PauliStabilizer[1], qc = QuantumChannel["Depolarizing"[1/3], {1}]},
        SameQ[qc[ps][[1, 2]]["Stabilizers"], ps["Stabilizers"]]
    ],
    True,
    {},
    TestID -> "Phase7.2-QC-Depolarizing-IdentityBranchUnchanged"
]


(* ============================================================================ *)
(* TIER G -- QMO direct equivalence to ps["M", ...] across multiple states     *)
(* ============================================================================ *)

(* Build a 2-qubit state via H1 + S1 + H2 then measure ZZ via QMO and via ps["M"];
   compare keys and post-state stabilizer sets. *)
VerificationTest[
    Module[{ps, qmo},
        ps = PauliStabilizer[2]["H", 1]["S", 1]["H", 2];
        qmo = QuantumMeasurementOperator["ZZ", {1, 2}];
        Sort[Keys[qmo[ps]]] === Sort[Keys[ps["M", "ZZ"]]]
    ],
    True,
    {},
    TestID -> "Phase7-QMO-vs-PS-M-Equivalence-ZZ"
]

(* Iterate over a few random Cliffords and check qmo vs ps["M"] keys match. *)
VerificationTest[
    Block[{},
        SeedRandom[20260507];
        AllTrue[
            Table[
                Module[{ps = PauliStabilizer["Random", 3], qmo = QuantumMeasurementOperator["XYZ", {1, 2, 3}]},
                    Sort[Keys[qmo[ps]]] === Sort[Keys[ps["M", "XYZ"]]]
                ],
                {12}
            ],
            TrueQ
        ]
    ],
    True,
    {},
    TestID -> "Phase7-QMO-vs-PS-M-Equivalence-Random3Q-12reps"
]


(* ============================================================================ *)
(* TIER H -- QMO + StabilizerFrame fallback path                                *)
(* ============================================================================ *)

(* Measurement on a StabilizerFrame produced by a T gate: must emit nonpaulibasis
   (frames currently always materialize) and produce a non-failure result. *)
VerificationTest[
    Module[{frame, qmo, result},
        frame = PauliStabilizer[1]["H", 1]["T", 1];   (* Frame from T gate *)
        qmo = QuantumMeasurementOperator["Z", {1}];
        result = qmo[frame];
        FreeQ[result, $Failed]
    ],
    True,
    {PauliStabilizer::nonpaulibasis},
    TestID -> "Phase7-QMO-on-StabilizerFrame-FallbackEvaluates"
]


(* ============================================================================ *)
(* TIER I -- Pauli sign prefix in QMO operator label                            *)
(* ============================================================================ *)

(* Document a sharp edge of Phase 7.1's label gate: a QMO built from           *)
(* QuantumOperator[-"XX"] does NOT carry a string-typed label "-XX" -- the      *)
(* unary minus produces a Times[-1, Superscript[X, CircleTimes[2]]] expression.*)
(* Phase 7.3 (2026-05-06) extended the detector to recognize this expression  *)
(* and route through the Pauli fast path.                                      *)
VerificationTest[
    Module[{psBell, label},
        psBell = PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "CNOT" -> {1, 2}};
        label = QuantumMeasurementOperator[QuantumOperator[-"XX"], {1, 2}]["Operator"]["Label"];
        StringQ[label]   (* still: False -- it is a Times expression *)
    ],
    False,
    {},
    TestID -> "Phase7-QMO-NegativePauli-Label-NotStringForm"
]


(* ============================================================================ *)
(* TIER J -- Phase 7.3 extended Pauli-label detection                          *)
(*                                                                              *)
(* QuantumOperator[-"XX"] has label Times[-1, Superscript["X", CircleTimes 2]].*)
(* Phase 7.3 maps that to "-XX" and routes via the AG fast path.                *)
(* ============================================================================ *)

(* Detector recognizes -Superscript[X, CircleTimes[2]] as "-XX". *)
VerificationTest[
    Wolfram`QuantumFramework`PackageScope`stabilizerPauliLabelFromQMO @
        QuantumMeasurementOperator[QuantumOperator[-"XX"], {1, 2}],
    "-XX",
    {},
    TestID -> "Phase7.3-Detector-NegXX-Recognized"
]

(* Detector recognizes Superscript[Z, CircleTimes[3]] as "ZZZ". *)
VerificationTest[
    Wolfram`QuantumFramework`PackageScope`stabilizerPauliLabelFromQMO @
        QuantumMeasurementOperator[QuantumOperator["ZZZ"], {1, 2, 3}],
    "ZZZ",
    {},
    TestID -> "Phase7.3-Detector-ZZZ-Recognized"
]


(* qmo[-XX measurement] on Bell state stays in tableau (no fallback message). *)
VerificationTest[
    Module[{psBell = PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "CNOT" -> {1, 2}},
            qmo = QuantumMeasurementOperator[QuantumOperator[-"XX"], {1, 2}]},
        Sort @ Keys @ qmo[psBell]
    ],
    {1},   (* -XX is the negative of the XX stabilizer; outcome is 1 *)
    {},
    TestID -> "Phase7.3-NegXX-on-Bell-DeterministicOutcomeOne"
]

(* qmo[XX measurement on Bell] still stays in tableau (positive Pauli). *)
VerificationTest[
    Module[{psBell = PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "CNOT" -> {1, 2}},
            qmo = QuantumMeasurementOperator[QuantumOperator["XX"], {1, 2}]},
        Sort @ Keys @ qmo[psBell]
    ],
    {0},
    {},
    TestID -> "Phase7.3-XX-on-Bell-DeterministicOutcomeZero"
]


(* ============================================================================ *)
(* TIER K -- Phase 7.3 detector smoke-checks                                    *)
(* ============================================================================ *)

(* PauliX-basis QMO does NOT match Phase 7.3 detection (label is Symbol     *)
(* "PauliX", not a Pauli string or Superscript form). It falls through to the *)
(* legacy fallback path with the nonpaulibasis info message. Document that.   *)
VerificationTest[
    Wolfram`QuantumFramework`PackageScope`stabilizerPauliLabelFromQMO @
        QuantumMeasurementOperator[QuantumBasis["PauliX"], {1}],
    Missing["NonPauliBasis"],
    {},
    TestID -> "Phase7.3-Detector-PauliXBasis-FallsThrough"
]


(* ============================================================================ *)
(* TIER L -- Phase 7.4 matrix-iteration Pauli detector                          *)
(*                                                                              *)
(* For QMOs constructed from explicit matrices (where the symbolic label is    *)
(* None or a non-Pauli expression), the detector iterates over 4^n * {+1,-1}  *)
(* Pauli candidates and checks against the QMO's MatrixRepresentation. For    *)
(* n <= 4 qubits the search is bounded.                                       *)
(* ============================================================================ *)

(* Detector recognizes single-qubit X matrix as "X". *)
VerificationTest[
    Wolfram`QuantumFramework`PackageScope`stabilizerPauliLabelFromQMO @
        QuantumMeasurementOperator[QuantumOperator[PauliMatrix[1]], {1}],
    "X",
    {},
    TestID -> "Phase7.4-Detector-X-Matrix-Recognized"
]

(* Detector recognizes single-qubit Y matrix as "Y" (complex matrix). *)
VerificationTest[
    Wolfram`QuantumFramework`PackageScope`stabilizerPauliLabelFromQMO @
        QuantumMeasurementOperator[QuantumOperator[PauliMatrix[2]], {1}],
    "Y",
    {},
    TestID -> "Phase7.4-Detector-Y-Matrix-Recognized"
]

(* Detector recognizes -Z matrix with sign. *)
VerificationTest[
    Wolfram`QuantumFramework`PackageScope`stabilizerPauliLabelFromQMO @
        QuantumMeasurementOperator[QuantumOperator[-PauliMatrix[3]], {1}],
    "-Z",
    {},
    TestID -> "Phase7.4-Detector-NegZ-Matrix-Recognized"
]

(* Detector recognizes a 2-qubit Pauli tensor product matrix. *)
VerificationTest[
    Wolfram`QuantumFramework`PackageScope`stabilizerPauliLabelFromQMO @
        QuantumMeasurementOperator[QuantumOperator[KroneckerProduct[PauliMatrix[1], PauliMatrix[3]]], {1, 2}],
    "XZ",
    {},
    TestID -> "Phase7.4-Detector-XZ-Matrix-Recognized"
]

(* 2-qubit -YY: complex multi-qubit Pauli with sign. *)
VerificationTest[
    Wolfram`QuantumFramework`PackageScope`stabilizerPauliLabelFromQMO @
        QuantumMeasurementOperator[QuantumOperator[-KroneckerProduct[PauliMatrix[2], PauliMatrix[2]]], {1, 2}],
    "-YY",
    {},
    TestID -> "Phase7.4-Detector-NegYY-Matrix-Recognized"
]

(* 3-qubit XYZ. *)
VerificationTest[
    Wolfram`QuantumFramework`PackageScope`stabilizerPauliLabelFromQMO @
        QuantumMeasurementOperator[
            QuantumOperator[KroneckerProduct[PauliMatrix[1], PauliMatrix[2], PauliMatrix[3]]],
            {1, 2, 3}
        ],
    "XYZ",
    {},
    TestID -> "Phase7.4-Detector-XYZ-Matrix-Recognized"
]

(* Identity matrix recognized as I^n. *)
VerificationTest[
    Wolfram`QuantumFramework`PackageScope`stabilizerPauliLabelFromQMO @
        QuantumMeasurementOperator[QuantumOperator[IdentityMatrix[4]], {1, 2}],
    "II",
    {},
    TestID -> "Phase7.4-Detector-Identity2Q-Matrix-Recognized"
]

(* Non-Pauli matrix returns Missing. The Hadamard matrix (1/Sqrt[2]){1,1;1,-1}
   is NOT a Pauli matrix (it's Clifford but acts as a basis change, not as a
   Pauli generator). *)
VerificationTest[
    Wolfram`QuantumFramework`PackageScope`stabilizerPauliLabelFromQMO @
        QuantumMeasurementOperator[QuantumOperator[(1/Sqrt[2]) {{1, 1}, {1, -1}}], {1}],
    Missing["NonPauliBasis"],
    {},
    TestID -> "Phase7.4-Detector-Hadamard-Matrix-Missing"
]

(* Off-diagonal non-Hermitian matrix (rank-1 projector |0><1|) returns Missing. *)
VerificationTest[
    Wolfram`QuantumFramework`PackageScope`stabilizerPauliLabelFromQMO @
        QuantumMeasurementOperator[QuantumOperator[{{0, 1}, {0, 0}}], {1}],
    Missing["NonPauliBasis"],
    {},
    TestID -> "Phase7.4-Detector-RankOneProjector-Missing"
]

(* Direct test that the fast path is taken with no nonpaulibasis message:    *)
(* matrix-form XX-measurement on Bell |Phi+> = (|00>+|11>)/Sqrt[2] is        *)
(* deterministic outcome 0 (XX is a Bell stabilizer). With Phase 7.4, this   *)
(* routes through the AG fast path. (n=2 chosen to dodge ROADMAP A.11.)      *)
VerificationTest[
    With[{psBell = PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "CNOT" -> {1, 2}},
          qmo = QuantumMeasurementOperator[
              QuantumOperator[KroneckerProduct[PauliMatrix[1], PauliMatrix[1]]], {1, 2}
          ]},
        Sort @ Keys @ qmo[psBell]
    ],
    {0},
    {},
    TestID -> "Phase7.4-MatrixXX-on-Bell-Deterministic-FastPath"
]

(* Same for matrix-form -XX measurement on Bell: -XX has eigenvalue -1 on    *)
(* |Phi+> (Bell is +1 eigenstate of XX, so -XX is -1 -> outcome 1).         *)
VerificationTest[
    With[{psBell = PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "CNOT" -> {1, 2}},
          qmo = QuantumMeasurementOperator[
              QuantumOperator[-KroneckerProduct[PauliMatrix[1], PauliMatrix[1]]], {1, 2}
          ]},
        Sort @ Keys @ qmo[psBell]
    ],
    {1},
    {},
    TestID -> "Phase7.4-MatrixNegXX-on-Bell-Deterministic-FastPath"
]

(* Matrix-form ZZ-measurement on Bell: Bell state has stabilizers {XX, ZZ} so
   ZZ measurement is deterministic outcome 0. Tests the matrix path matches the
   string path. *)
VerificationTest[
    With[{psBell = PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "CNOT" -> {1, 2}},
          qmoMat = QuantumMeasurementOperator[
              QuantumOperator[KroneckerProduct[PauliMatrix[3], PauliMatrix[3]]], {1, 2}
          ],
          qmoStr = QuantumMeasurementOperator["ZZ", {1, 2}]},
        Sort[Keys[qmoMat[psBell]]] === Sort[Keys[qmoStr[psBell]]]
    ],
    True,
    {},
    TestID -> "Phase7.4-MatrixZZ-vs-StringZZ-Equivalence"
]

(* Matrix-form anticommuting measurement on GHZ-3 (XII via matrix): XII
   anticommutes with ZZI stabilizer of GHZ, so non-deterministic. *)
VerificationTest[
    With[{psGHZ = PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "CNOT" -> {1, 2}, "CNOT" -> {2, 3}},
          qmoMat = QuantumMeasurementOperator[
              QuantumOperator[KroneckerProduct[PauliMatrix[1], PauliMatrix[0], PauliMatrix[0]]],
              {1, 2, 3}
          ]},
        Sort @ Keys @ qmoMat[psGHZ]
    ],
    {0, 1},
    {},
    TestID -> "Phase7.4-MatrixXII-on-GHZ-NonDeterministic"
]

(* Cap test: bumping the search cap allows n = 5 detection.                  *)
(* Wrap in Block to scope the override. *)
VerificationTest[
    Block[{Wolfram`QuantumFramework`PackageScope`$stabilizerPauliMatrixSearchMaxQubits = 5},
        Wolfram`QuantumFramework`PackageScope`stabilizerPauliFromMatrix[
            KroneckerProduct[PauliMatrix[1], PauliMatrix[2], PauliMatrix[3], PauliMatrix[1], PauliMatrix[2]],
            5
        ]
    ],
    "XYZXY",
    {},
    TestID -> "Phase7.4-Detector-MaxQubitCap-OverrideBlock"
]

(* Default cap rejects n = 5 with TooManyQubits. *)
VerificationTest[
    Wolfram`QuantumFramework`PackageScope`stabilizerPauliFromMatrix[
        KroneckerProduct @@ ConstantArray[PauliMatrix[1], 5],
        5
    ],
    Missing["TooManyQubits"],
    {},
    TestID -> "Phase7.4-Detector-MaxQubitCap-DefaultRejects"
]

(* Cross-check: the matrix-iteration detector agrees with the string-form fast
   path on a battery of seeded random Cliffords. *)
VerificationTest[
    Block[{},
        SeedRandom[20260507];
        AllTrue[
            Table[
                Module[{ps, qmoStr, qmoMat, sStr, sMat},
                    ps = PauliStabilizer["Random", 3];
                    qmoStr = QuantumMeasurementOperator["XYZ", {1, 2, 3}];
                    qmoMat = QuantumMeasurementOperator[
                        QuantumOperator[KroneckerProduct[PauliMatrix[1], PauliMatrix[2], PauliMatrix[3]]],
                        {1, 2, 3}
                    ];
                    sStr = Sort[Keys[qmoStr[ps]]];
                    sMat = Sort[Keys[qmoMat[ps]]];
                    sStr === sMat
                ],
                {12}
            ],
            TrueQ
        ]
    ],
    True,
    {},
    TestID -> "Phase7.4-MatrixForm-vs-StringForm-Random3Q-12reps"
]

(* Phase 7.4 detector handles single-qubit Pauli matrices correctly. ROADMAP A.11
   bug (KroneckerProduct on single element) is worked around in the helper. *)
VerificationTest[
    Wolfram`QuantumFramework`PackageScope`stabilizerPauliFromMatrix[PauliMatrix[1], 1],
    "X",
    {},
    TestID -> "Phase7.4-Detector-SingleQubit-Workaround"
]

(* Detector returns DimMismatch on a non-square matrix. *)
VerificationTest[
    Wolfram`QuantumFramework`PackageScope`stabilizerPauliFromMatrix[
        {{1, 0}, {0, 0}, {0, 0}, {0, 1}}, 1
    ],
    Missing["DimMismatch"],
    {},
    TestID -> "Phase7.4-Detector-NonSquareMatrix-DimMismatch"
]
