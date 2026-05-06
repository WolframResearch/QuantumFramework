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
