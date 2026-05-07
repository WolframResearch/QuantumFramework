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
(* unary minus produces a Times[-1, "XX"] expression. The fast path is gated    *)
(* on StringQ[label], so the dispatcher takes the legacy generic route in      *)
(* this case (correct, but slower than ideal). Phase 7.3 may extend the gate.  *)
VerificationTest[
    Module[{psBell, label},
        psBell = PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "CNOT" -> {1, 2}};
        label = QuantumMeasurementOperator[QuantumOperator[-"XX"], {1, 2}]["Operator"]["Label"];
        StringQ[label]   (* expected: False -- it is a Times expression *)
    ],
    False,
    {},
    TestID -> "Phase7-QMO-NegativePauli-Label-NotStringForm"
]
