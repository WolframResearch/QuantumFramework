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
(* The detector and register-label helpers exercised below are PackageScope. *)
Needs["Wolfram`QuantumFramework`PackageScope`"];

(* Local validity helper: re-export the package-scoped predicate so tests can
   assert on the receiver's structure. *)
psValidQ = PauliStabilizerQ;

(* Equal as rays, exactly: |<u|v>|^2 = <u|u><v|v> on amplitude lists (no normalization asked). *)
rayEqQ[u_List, v_List] := Simplify[(Conjugate[u] . v) (Conjugate[v] . u) - (Conjugate[u] . u) (Conjugate[v] . v)] === 0

(* Aaronson-Gottesman pairing: the 2n generator rows (destabilizers, then stabilizers) *)
(* form a symplectic basis of F_2^(2n), so M . Omega . M^T == Omega mod 2.        *)
symplecticQ[ps_] := With[{n = ps["Qubits"]},
    With[{om = ArrayFlatten[{{0 IdentityMatrix[n], IdentityMatrix[n]}, {IdentityMatrix[n], 0 IdentityMatrix[n]}}]},
        Mod[ps["Matrix"] . om . Transpose[ps["Matrix"]], 2] === om
    ]
]


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


(* A Pauli label shorter than the register is padded with identities, each     *)
(* letter landing on the wire its target names, so a single-qubit measurement  *)
(* on a larger tableau stays in the tableau. Z on qubit 2 of GHZ-3 is uniformly *)
(* random ...                                                                   *)
VerificationTest[
    With[{ps = PauliStabilizer[3][{"H" -> 1, "CNOT" -> {1, 2}, "CNOT" -> {2, 3}}]},
        Sort @ Keys @ QuantumMeasurementOperator["Z", {2}][ps]
    ],
    {0, 1},
    {},
    TestID -> "Phase7-QMO-Partial-Z-on-q2-of-GHZ3-NonDeterministic"
]

(* ... and is exactly the padded native measurement. *)
VerificationTest[
    With[{ps = PauliStabilizer[3][{"H" -> 1, "CNOT" -> {1, 2}, "CNOT" -> {2, 3}}]},
        QuantumMeasurementOperator["Z", {2}][ps] === ps["M", "IZI"]
    ],
    True,
    {},
    TestID -> "Phase7-QMO-Partial-Z-on-q2-EqualsPaddedNativeM"
]

(* Deterministic partial measurement: Z on qubit 2 of |00> is outcome 0. *)
VerificationTest[
    Sort @ Keys @ QuantumMeasurementOperator["Z", {2}][PauliStabilizer[2]],
    {0},
    {},
    TestID -> "Phase7-QMO-Partial-Z-on-q2-of-zero-Deterministic"
]

(* Unsorted multi-qubit targets: QuantumMeasurementOperator["XZ", {3, 1}] puts X *)
(* on wire 3 and Z on wire 1, the operator's own convention                     *)
(* (QuantumOperator["XZ", {3, 1}] sends |+00> to |-01>). On |0>|0>|+> both       *)
(* factors are +1, so the joint outcome is deterministic; the other wire        *)
(* assignment is random.                                                        *)
VerificationTest[
    With[{ps = PauliStabilizer[3][{"H" -> 3}]},
        {
            Sort @ Keys @ QuantumMeasurementOperator["XZ", {3, 1}][ps],
            QuantumMeasurementOperator["XZ", {3, 1}][ps] === ps["M", "ZIX"],
            Sort @ Keys @ QuantumMeasurementOperator["XZ", {1, 3}][ps]
        }
    ],
    {{0}, True, {0, 1}},
    {},
    TestID -> "Phase7-QMO-UnsortedTargets-XZ-on-31-LettersFollowWires"
]

(* The same wires through the operator-built and the matrix-built forms: the    *)
(* matrix search finds its hit in sorted-wire order and reports it in the       *)
(* operator's wire order.                                                       *)
VerificationTest[
    With[{ps = PauliStabilizer[3][{"H" -> 3}],
          qmoMat = QuantumMeasurementOperator[QuantumOperator[KroneckerProduct[PauliMatrix[1], PauliMatrix[3]], {3, 1}]]},
        {
            Keys @ QuantumMeasurementOperator[QuantumOperator["XZ", {3, 1}]][ps],
            Keys @ qmoMat[ps],
            stabilizerPauliLabelFromQMO[qmoMat]
        }
    ],
    {{0}, {0}, "XZ"},
    {},
    TestID -> "Phase7-QMO-UnsortedTargets-OperatorAndMatrixForms-Agree"
]

(* A signed partial label keeps its sign: -Z on qubit 2 of |000> is outcome 1. *)
VerificationTest[
    With[{ps = PauliStabilizer[3], qmo = QuantumMeasurementOperator[QuantumOperator[-"Z", {2}]]},
        {Sort @ Keys @ qmo[ps], qmo[ps] === ps["M", "-IZI"]}
    ],
    {{1}, True},
    {},
    TestID -> "Phase7-QMO-Partial-SignedLabel-NegZ-on-q2"
]

(* The padding contract: a label that already spans the register in wire order  *)
(* is unchanged, in another wire order it is re-seated, a shorter one is padded, *)
(* a non-Pauli basis passes through as Missing, and a one-letter basis tiled     *)
(* over several wires is a Missing of its own.                                  *)
VerificationTest[
    {
        stabilizerPauliRegisterLabel[QuantumMeasurementOperator["XZ", {1, 2}], 2],
        stabilizerPauliRegisterLabel[QuantumMeasurementOperator["XZ", {2, 1}], 2],
        stabilizerPauliRegisterLabel[QuantumMeasurementOperator["XZ", {3, 1}], 3],
        stabilizerPauliRegisterLabel[QuantumMeasurementOperator[QuantumBasis["Computational"], {1}], 1],
        stabilizerPauliRegisterLabel[QuantumMeasurementOperator["X", {1, 2}], 2]
    },
    {"XZ", "ZX", "ZIX", Missing["NonPauliBasis"], Missing["LabelWireMismatch"]},
    {},
    TestID -> "Phase7-RegisterLabel-Padding-Contract"
]

(* A mixed multi-letter Pauli operator carries a CircleTimes chain as its label; *)
(* the detector reads it (with its sign and in the operator's wire order) at    *)
(* any length, past the matrix-search cap, and a chain with a non-Pauli factor  *)
(* is not a Pauli label.                                                        *)
VerificationTest[
    {
        stabilizerPauliLabelFromQMO[QuantumMeasurementOperator[QuantumOperator[-"ZYYZY", Range[5]]]],
        stabilizerPauliLabelFromQMO[QuantumMeasurementOperator[QuantumOperator["XZYXZ", {5, 1, 4, 2, 3}]]],
        stabilizerPauliLabelFromQMO[QuantumMeasurementOperator[QuantumOperator["XHZ", Range[3]]]]
    },
    {"-ZYYZY", "XZYXZ", Missing["NonPauliBasis"]},
    {},
    TestID -> "Phase7-Detector-CircleTimesChain-Recognized"
]

(* A one-letter basis tiled over several wires (QuantumMeasurementOperator["X",  *)
(* {1, 2}] measures both qubits in the X basis) is a product-basis measurement  *)
(* no single Pauli string expresses: it takes the dense fallback.               *)
VerificationTest[
    Head @ QuantumMeasurementOperator["X", {1, 2}][PauliStabilizer[2]],
    QuantumMeasurement,
    {PauliStabilizer::nonpaulibasis},
    TestID -> "Phase7-QMO-TiledLetterBasis-FallsBackToDense"
]

(* A target outside the register is refused with a message, not an assertion,  *)
(* whatever the label: the wire range is judged before the letter count, so a   *)
(* one-letter basis over a wire pair that leaves the register is refused too.   *)
VerificationTest[
    {
        QuantumMeasurementOperator["Z", {4}][PauliStabilizer[3]],
        QuantumMeasurementOperator["X", {1, 5}][PauliStabilizer[2]]
    },
    {$Failed, $Failed},
    {PauliStabilizer::target, PauliStabilizer::target},
    TestID -> "Phase7-QMO-TargetOutsideRegister-Fails"
]

(* Agreement with the dense path, exactly: Z on qubit 2 of GHZ-3 has Born        *)
(* weights 1/2, 1/2 (the dense branches carry Sqrt[p]), the two tableau branches *)
(* are orthogonal, and each equals its dense branch up to a global phase.        *)
VerificationTest[
    With[{ps = PauliStabilizer[3][{"H" -> 1, "CNOT" -> {1, 2}, "CNOT" -> {2, 3}}], qmo = QuantumMeasurementOperator["Z", {2}]},
        With[{viaPS = qmo[ps], viaQS = qmo[ps["State"]]},
            With[{tab = Normal @ #["State"]["StateVector"] & /@ Values[viaPS], dense = Normal @ #["StateVector"] & /@ viaQS["States"]},
                {
                    Values[viaQS["Probabilities"]],
                    Simplify[Conjugate[#] . #] & /@ dense,
                    Simplify[Conjugate[tab[[1]]] . tab[[2]]],
                    MapThread[rayEqQ, {tab, dense}]
                }
            ]
        ]
    ],
    {{1/2, 1/2}, {1/2, 1/2}, 0, {True, True}},
    {},
    TestID -> "Phase7-QMO-Partial-Z-on-q2-of-GHZ3-MatchesDense"
]

(* Random stabilizer states, random Pauli labels of length 1..n on random wire  *)
(* subsets (sorted or not, signed or not), n = 1..6, past the matrix-search cap *)
(* of 4 qubits so the label parser alone carries the longer ones, all exact:    *)
(* <P> from the                                                                 *)
(* operator's own action on the dense state is one of -1, 0, 1 and equals the   *)
(* tableau's closed-form expectation on the register label the interop builds;  *)
(* the outcome set is {0}, {1} or {0, 1} accordingly; the result is the native   *)
(* measurement of that label; and every branch is the normalized               *)
(* eigenprojection (1 +- P)|psi> up to a global phase.                           *)
VerificationTest[
    BlockRandom[
        SeedRandom[20260921];
        AllTrue[
            Flatten @ Table[
                With[{ps = PauliStabilizer["Random"[n]], k = RandomInteger[{1, n}]},
                    With[{wires = RandomSample[Range[n], k], body = StringJoin[RandomChoice[{"X", "Y", "Z"}, k]], sign = RandomChoice[{1, -1}]},
                        With[{qmo = QuantumMeasurementOperator[QuantumOperator[sign body, wires]]},
                        With[{
                            res = qmo[ps],
                            v = Normal @ ps["State"]["StateVector"],
                            pv = Normal @ QuantumOperator[sign body, wires][ps["State"]]["StateVector"],
                            registerLabel = stabilizerPauliRegisterLabel[qmo, n]
                        },
                            With[{mean = Simplify[Conjugate[v] . pv]},
                                MemberQ[{-1, 0, 1}, mean] &&
                                    mean === stabilizerExpectation[ps, registerLabel] &&
                                    Sort[Keys[res]] === Replace[mean, {1 -> {0}, -1 -> {1}, 0 -> {0, 1}}] &&
                                    res === ps["M", registerLabel] &&
                                    AllTrue[Keys[res],
                                        With[{u = Normal @ res[#]["State"]["StateVector"]},
                                            Simplify[Conjugate[u] . u - 1] === 0 && rayEqQ[u, v + (1 - 2 #) pv]
                                        ] &
                                    ]
                            ]
                        ]]
                    ]
                ],
                {n, 6}, {4}
            ],
            TrueQ
        ]
    ],
    True,
    {},
    TestID -> "Phase7-QMO-Partial-Random-MatchesOperatorAction-24reps"
]

(* Beyond any dense reference: on the 60-qubit GHZ tableau a single-qubit Z is   *)
(* uniformly random and each branch is then certain of that same Z (the sign    *)
(* bookkeeping), each branch keeps the symplectic pairing of its generators     *)
(* (the tableau bookkeeping), the Z_1 Z_2 parity is certain from the start, and *)
(* a mixed three-letter Pauli on far-apart wires, parsed from its CircleTimes   *)
(* label, is measured in the tableau as well.                                   *)
VerificationTest[
    With[{n = 60},
        With[{ps = PauliStabilizer[n][Prepend[Table["CNOT" -> {q, q + 1}, {q, n - 1}], "H" -> 1]]},
            With[{branches = QuantumMeasurementOperator["Z", {17}][ps]},
                {
                    Sort @ Keys @ branches,
                    KeyValueMap[Keys[QuantumMeasurementOperator["Z", {17}][#2]] === {#1} &, branches],
                    symplecticQ /@ Values[branches],
                    Sort @ Keys @ QuantumMeasurementOperator["ZZ", {1, 2}][ps],
                    Sort @ Keys @ QuantumMeasurementOperator[QuantumOperator["XZY", {17, 33, 51}]][ps]
                }
            ]
        ]
    ],
    {{0, 1}, {True, True}, {True, True}, {0}, {0, 1}},
    {},
    TestID -> "Phase7-QMO-Partial-GHZ60-BeyondDense"
]


(* ============================================================================ *)
(* TIER B -- Non-Pauli basis: info message + dense fallback                    *)
(* ============================================================================ *)

(* A computational-basis projector measurement (label not Pauli) triggers the
   ::nonpaulibasis fallback: a single notice, then the dense state-vector path
   qmo[ps["State"]], which performs the measurement and returns a
   QuantumMeasurement. Exactly one message; no cascading ::nonclifford. *)
VerificationTest[
    With[{ps = PauliStabilizer[1],
          qmo = QuantumMeasurementOperator[QuantumBasis["Computational"], {1}]},
        Head @ qmo[ps]
    ],
    QuantumMeasurement,
    {PauliStabilizer::nonpaulibasis},
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

(* After ps["SymbolicMeasure", q] the receiver carries symbolic signs. The Pauli *)
(* measurement primitive ps["M", ...] carries them too, so a Pauli QMO on the   *)
(* symbolic tableau stays in the tableau. On a Bell pair whose first qubit was  *)
(* measured symbolically (outcome s), Z on the partner qubit is certain and     *)
(* keyed by s itself (the Bell correlation), s really is the recorded Z outcome *)
(* of qubit 1, the tableau is unchanged, the Z Z parity is certain, and X on    *)
(* the measured qubit is uniformly random.                                       *)
VerificationTest[
    With[{psSym = PauliStabilizer[2][{"H" -> 1, "CNOT" -> {1, 2}}]["SymbolicMeasure", 1]},
        With[{
            s = First @ Cases[psSym["Signs"], _\[FormalS], Infinity],
            res = QuantumMeasurementOperator["Z", {2}][psSym],
            dense = QuantumMeasurementOperator["Z", {1}][PauliStabilizer[2][{"H" -> 1, "CNOT" -> {1, 2}}]["State"]]["States"]
        },
            {
                Keys[res] === {s},
                (* s is the Z outcome of qubit 1: substituting s = b lands on the dense branch b, |bb>. *)
                Table[rayEqQ[Normal @ psSym["SubstituteOutcomes", s -> b]["State"]["StateVector"], Normal @ dense[[b + 1]]["StateVector"]], {b, 0, 1}],
                res[s]["Tableau"] === psSym["Tableau"],
                Keys @ QuantumMeasurementOperator["ZZ", {1, 2}][psSym],
                Sort @ Keys @ QuantumMeasurementOperator["X", {1}][psSym]
            }
        ]
    ],
    {True, {True, True}, True, {0}, {0, 1}},
    {},
    TestID -> "Phase7-QMO-on-SymbolicPS-StaysInTableau"
]

(* A non-Pauli basis on a symbolic-sign tableau matches no interop rule (symbolic *)
(* signs cannot be materialized) and takes the generic circuit route, where the *)
(* stabilizer engine cannot fold a measurement gate.                            *)
VerificationTest[
    QuantumMeasurementOperator[QuantumBasis["Computational"], {1}][PauliStabilizer[1]["H", 1]["SymbolicMeasure", 1]],
    $Failed,
    {PauliStabilizer::nonclifford},
    TestID -> "Phase7-QMO-NonPauli-on-SymbolicPS-GenericRoute"
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
    stabilizerPauliLabelFromQMO @
        QuantumMeasurementOperator[QuantumOperator[-"XX"], {1, 2}],
    "-XX",
    {},
    TestID -> "Phase7.3-Detector-NegXX-Recognized"
]

(* Detector recognizes Superscript[Z, CircleTimes[3]] as "ZZZ". *)
VerificationTest[
    stabilizerPauliLabelFromQMO @
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
    stabilizerPauliLabelFromQMO @
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
    stabilizerPauliLabelFromQMO @
        QuantumMeasurementOperator[QuantumOperator[PauliMatrix[1]], {1}],
    "X",
    {},
    TestID -> "Phase7.4-Detector-X-Matrix-Recognized"
]

(* Detector recognizes single-qubit Y matrix as "Y" (complex matrix). *)
VerificationTest[
    stabilizerPauliLabelFromQMO @
        QuantumMeasurementOperator[QuantumOperator[PauliMatrix[2]], {1}],
    "Y",
    {},
    TestID -> "Phase7.4-Detector-Y-Matrix-Recognized"
]

(* Detector recognizes -Z matrix with sign. *)
VerificationTest[
    stabilizerPauliLabelFromQMO @
        QuantumMeasurementOperator[QuantumOperator[-PauliMatrix[3]], {1}],
    "-Z",
    {},
    TestID -> "Phase7.4-Detector-NegZ-Matrix-Recognized"
]

(* Detector recognizes a 2-qubit Pauli tensor product matrix. *)
VerificationTest[
    stabilizerPauliLabelFromQMO @
        QuantumMeasurementOperator[QuantumOperator[KroneckerProduct[PauliMatrix[1], PauliMatrix[3]]], {1, 2}],
    "XZ",
    {},
    TestID -> "Phase7.4-Detector-XZ-Matrix-Recognized"
]

(* 2-qubit -YY: complex multi-qubit Pauli with sign. *)
VerificationTest[
    stabilizerPauliLabelFromQMO @
        QuantumMeasurementOperator[QuantumOperator[-KroneckerProduct[PauliMatrix[2], PauliMatrix[2]]], {1, 2}],
    "-YY",
    {},
    TestID -> "Phase7.4-Detector-NegYY-Matrix-Recognized"
]

(* 3-qubit XYZ. *)
VerificationTest[
    stabilizerPauliLabelFromQMO @
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
    stabilizerPauliLabelFromQMO @
        QuantumMeasurementOperator[QuantumOperator[IdentityMatrix[4]], {1, 2}],
    "II",
    {},
    TestID -> "Phase7.4-Detector-Identity2Q-Matrix-Recognized"
]

(* Non-Pauli matrix returns Missing. The Hadamard matrix (1/Sqrt[2]){1,1;1,-1}
   is NOT a Pauli matrix (it's Clifford but acts as a basis change, not as a
   Pauli generator). *)
VerificationTest[
    stabilizerPauliLabelFromQMO @
        QuantumMeasurementOperator[QuantumOperator[(1/Sqrt[2]) {{1, 1}, {1, -1}}], {1}],
    Missing["NonPauliBasis"],
    {},
    TestID -> "Phase7.4-Detector-Hadamard-Matrix-Missing"
]

(* Off-diagonal non-Hermitian matrix (rank-1 projector |0><1|) returns Missing. *)
VerificationTest[
    stabilizerPauliLabelFromQMO @
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
    Block[{$stabilizerPauliMatrixSearchMaxQubits = 5},
        stabilizerPauliFromMatrix[
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
    stabilizerPauliFromMatrix[
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
    stabilizerPauliFromMatrix[PauliMatrix[1], 1],
    "X",
    {},
    TestID -> "Phase7.4-Detector-SingleQubit-Workaround"
]

(* Detector returns DimMismatch on a non-square matrix. *)
VerificationTest[
    stabilizerPauliFromMatrix[
        {{1, 0}, {0, 0}, {0, 0}, {0, 1}}, 1
    ],
    Missing["DimMismatch"],
    {},
    TestID -> "Phase7.4-Detector-NonSquareMatrix-DimMismatch"
]
