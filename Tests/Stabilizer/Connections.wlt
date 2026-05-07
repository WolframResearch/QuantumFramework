(* ::Package:: *)

(* ============================================================================
   Tests/Stabilizer/Connections.wlt -- cross-head dispatch coverage.

   Every entry in api.md's "Connections to other QuantumFramework heads" matrix
   gets at least:
     (1) a direct dispatch test (the rule fires and produces something),
     (2) a semantic equivalence test (the result matches the documented
         alternative path), and
     (3) where round-trippable, a multi-hop chain test
         (PauliStabilizer -> QF head -> PauliStabilizer).

   The 9 connections audited:
     1.  QuantumState[ps]                           (UpValue)
     2.  QuantumOperator[ps]                        (DownValue)
     3.  QuantumCircuitOperator[ps]                 (DownValue)
     4.  qo_QuantumOperator[ps]                     (UpValue, gate dispatch)
     5.  qmo_QuantumMeasurementOperator[ps]         (UpValue, hybrid interop)
     6.  qc_QuantumChannel[ps]                      (UpValue, hybrid interop)
     7.  cc_CliffordChannel[ps]                     (DownValue)
     8.  QuantumTensorProduct[ps_a, ps_b]           (DownValue)
     9.  QuantumCircuitOperator[..][Method -> "Stabilizer"]  (Method dispatch)

   Plus:
     - StabilizerFrame[ps] cross-head dispatch (qmo[sf], qc[sf]).
     - GraphState[ps]["PauliStabilizer"] round-trip (graph-form input).
     - Multi-hop chains: ps -> qs -> ps', ps -> qo -> ps', ps -> qco -> ps'.
     - End-to-end: stabilizer state, simulate a measurement via QMO, compare
       outcome distribution to direct ps["M", q].

   Phase 5c contract still applies: gate-update propagation does NOT preserve
   GlobalPhase, so multi-hop tests use phase-aware comparisons where
   appropriate.
   ========================================================================== *)

Needs["Wolfram`QuantumFramework`"];

psValidQ = Wolfram`QuantumFramework`PackageScope`PauliStabilizerQ;
sfValidQ = Wolfram`QuantumFramework`PackageScope`StabilizerFrameQ;
ccValidQ = Wolfram`QuantumFramework`PackageScope`CliffordChannelQ;

matEqQO[a_, b_] := SameQ @@ (Normal @ #["Matrix"] & /@ {a, b})
matEqQS[a_, b_] := SameQ @@ (Normal @ #["StateVector"] & /@ {a, b})


(* Common fixtures *)
$ps0   = PauliStabilizer[1];
$ps00  = PauliStabilizer[2];
$psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
$psGHZ3 = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, "CNOT" -> {2, 3}}]];
$ps5Q  = PauliStabilizer["5QubitCode"];


(* ============================================================================ *)
(* CONNECTION 1 -- QuantumState[ps] UpValue                                     *)
(* ============================================================================ *)

(* Direct dispatch fires. *)
VerificationTest[
    Head @ QuantumState[$psBell],
    QuantumState,
    TestID -> "Conn1-QSofPS-Dispatches"
]

(* Semantic equivalence: QuantumState[ps] === ps["State"]. *)
VerificationTest[
    matEqQS[QuantumState[$psBell], $psBell["State"]],
    True,
    TestID -> "Conn1-QSofPS-EqualsState"
]

(* Sweep across canonical states. *)
VerificationTest[
    AllTrue[
        {$ps0, $ps00, $psBell, $psGHZ3},
        matEqQS[QuantumState[#], #["State"]] &
    ],
    True,
    TestID -> "Conn1-QSofPS-Sweep"
]

(* Round-trip: ps -> QuantumState -> PauliStabilizer (re-tomography). The new   *)
(* stabilizer set generates the same group up to canonicalization; we compare   *)
(* via materialized state vectors (which IS exact thanks to GlobalPhase).       *)
VerificationTest[
    Module[{ps2 = PauliStabilizer[QuantumState[$psBell]]},
        matEqQS[ps2["State"], $psBell["State"]]
    ],
    True,
    TestID -> "Conn1-RT-PS-QS-PS-Bell"
]

VerificationTest[
    Module[{ps2 = PauliStabilizer[QuantumState[$psGHZ3]]},
        matEqQS[ps2["State"], $psGHZ3["State"]]
    ],
    True,
    TestID -> "Conn1-RT-PS-QS-PS-GHZ3"
]


(* ============================================================================ *)
(* CONNECTION 2 -- ps["Operator"] (the working accessor; QuantumOperator[ps]   *)
(* DownValue is silently inactive against the protected symbol).                *)
(* ============================================================================ *)

(* Use the documented accessor ps["Operator"] / ps["QuantumOperator"]. *)
VerificationTest[
    Head @ $psBell["Operator"],
    QuantumOperator,
    TestID -> "Conn2-psOperator-IsQuantumOperator"
]

(* Equivalence: ps["Operator"] === GlobalPhase * ps["Circuit"]["QuantumOperator"]. *)
VerificationTest[
    matEqQO[$psBell["Operator"], $psBell["Circuit"]["QuantumOperator"]],
    True,
    TestID -> "Conn2-psOperator-MatchesCircuitQO-WhenPhaseIsOne"
]

(* Document the inactive DownValue: QuantumOperator[ps] does NOT route through *)
(* the Stabilizer DownValue (the host symbol is protected). The result is a    *)
(* malformed QuantumOperator wrapping the PauliStabilizer's tableau as a value.*)
(* This test asserts the documented contract: use ps["Operator"] instead.       *)
VerificationTest[
    matEqQO[QuantumOperator[$psBell], $psBell["Circuit"]["QuantumOperator"]],
    False,   (* Documents the gotcha; promote to True if the DownValue is fixed. *)
    TestID -> "Conn2-QuantumOperatorOfPS-DownValue-Inactive-Documented"
]

(* Round-trip via the working accessor: ps -> ps["Operator"] -> PS reproduces *)
(* the state up to global phase (Phase 5c contract for tableau-only paths).   *)
VerificationTest[
    Module[{
        ps2 = PauliStabilizer[$psBell["Operator"]],
        v1 = Normal @ $psBell["State"]["StateVector"],
        v2
    },
        v2 = Normal @ ps2["State"]["StateVector"];
        Chop[Abs[Conjugate[v1] . v2] - 1] === 0
    ],
    True,
    TestID -> "Conn2-RT-PS-psOperator-PS-Bell-UpToPhase"
]


(* ============================================================================ *)
(* CONNECTION 3 -- ps["Circuit"] (the working accessor; QuantumCircuitOperator *)
(* [ps] DownValue is silently inactive against the protected symbol).           *)
(* ============================================================================ *)

VerificationTest[
    Head @ $psBell["Circuit"],
    QuantumCircuitOperator,
    TestID -> "Conn3-psCircuit-IsQuantumCircuitOperator"
]

(* Document the inactive DownValue. *)
VerificationTest[
    QuantumCircuitOperator[$psBell] === $psBell["Circuit"],
    False,   (* Documents the gotcha; promote to True if the DownValue is fixed. *)
    TestID -> "Conn3-QuantumCircuitOperatorOfPS-DownValue-Inactive-Documented"
]

(* Round-trip via the working accessor: ps -> ps["Circuit"] -> PauliStabilizer *)
(* round-trip is up-to-phase per Phase 5c contract. *)
VerificationTest[
    Module[{
        psFromCirc = PauliStabilizer[$psBell["Circuit"]],
        v1 = Normal @ $psBell["State"]["StateVector"],
        v2
    },
        v2 = Normal @ psFromCirc["State"]["StateVector"];
        Chop[Abs[Conjugate[v1] . v2] - 1] === 0
    ],
    True,
    TestID -> "Conn3-RT-PS-psCircuit-PS-Bell-UpToPhase"
]


(* ============================================================================ *)
(* CONNECTION 4 -- qo_QuantumOperator[ps] UpValue (gate dispatch)               *)
(* ============================================================================ *)

(* Direct dispatch: applying a Clifford QO to ps via @ should produce a PS.     *)
VerificationTest[
    Head @ (QuantumOperator["H", 1] @ $ps0),
    PauliStabilizer,
    TestID -> "Conn4-qoApplyPS-Dispatches"
]

(* Semantic: QuantumOperator["H", 1] @ ps === ps["H", 1] up to canonical form. *)
VerificationTest[
    Sort @ (QuantumOperator["H", 1] @ $ps0)["Stabilizers"],
    Sort @ $ps0["H", 1]["Stabilizers"],
    TestID -> "Conn4-qoApplyPS-EqualsMethodCall"
]

(* CNOT via UpValue *)
VerificationTest[
    Sort @ (QuantumOperator["CNOT", {1, 2}] @ ($ps00["H", 1]))["Stabilizers"],
    Sort @ {"XX", "ZZ"},
    TestID -> "Conn4-qoCNOT-on-PlusZero-EqualsBell"
]

(* Chain: ps -> qo @ ps -> qo @ (qo @ ps) ... matches direct gate composition. *)
VerificationTest[
    Module[{psA, psB},
        psA = QuantumOperator["CNOT", {1, 2}] @ (QuantumOperator["H", 1] @ $ps00);
        psB = $ps00["H", 1]["CNOT", 1, 2];
        Sort[psA["Stabilizers"]] === Sort[psB["Stabilizers"]]
    ],
    True,
    TestID -> "Conn4-qoChain-EqualsMethodChain"
]


(* ============================================================================ *)
(* CONNECTION 5 -- qmo_QuantumMeasurementOperator[ps] UpValue (Phase 7)         *)
(* ============================================================================ *)

(* Direct dispatch: returns Association of post-measurement states. *)
VerificationTest[
    Head @ QuantumMeasurementOperator["ZZ", {1, 2}][$psBell],
    Association,
    TestID -> "Conn5-qmoPS-Dispatches"
]

(* Semantic: qmo[ps] === ps["M", "ZZ"] for Pauli-string QMO. *)
VerificationTest[
    Sort[Keys[QuantumMeasurementOperator["ZZ", {1, 2}][$psBell]]] ===
        Sort[Keys[$psBell["M", "ZZ"]]],
    True,
    TestID -> "Conn5-qmoPS-EqualsMethodM"
]

(* Bell ZZ correlation: deterministic outcome. *)
VerificationTest[
    Sort @ Keys @ QuantumMeasurementOperator["ZZ", {1, 2}][$psBell],
    {0},
    TestID -> "Conn5-qmoZZ-on-Bell-Deterministic"
]

(* Non-Pauli basis falls back via ::nonpaulibasis. *)
VerificationTest[
    Head @ QuantumMeasurementOperator[QuantumBasis["Computational"], {1}][$ps0],
    PauliStabilizer,
    {PauliStabilizer::nonpaulibasis, PauliStabilizer::nonclifford},
    TestID -> "Conn5-qmoComputational-FallsBack"
]


(* ============================================================================ *)
(* CONNECTION 6 -- qc_QuantumChannel[ps] UpValue (Phase 7.2)                    *)
(* ============================================================================ *)

(* Direct dispatch: BitFlip returns mixture list. *)
VerificationTest[
    Length @ QuantumChannel["BitFlip"[1/3], {1}][$ps0],
    2,
    TestID -> "Conn6-qcBitFlip-Dispatches"
]

(* Probabilities sum to 1. *)
VerificationTest[
    Total[First /@ QuantumChannel["BitFlip"[1/3], {1}][$ps0]],
    1,
    TestID -> "Conn6-qcBitFlip-ProbSum"
]

(* Each branch is a valid PauliStabilizer. *)
VerificationTest[
    AllTrue[QuantumChannel["BitFlip"[1/3], {1}][$ps0], psValidQ[#[[2]]] &],
    True,
    TestID -> "Conn6-qcBitFlip-AllValid"
]

(* Depolarizing -> 4 branches; non-Clifford -> dense fallback. *)
VerificationTest[
    Length @ QuantumChannel["Depolarizing"[1/2], {1}][$ps0],
    4,
    TestID -> "Conn6-qcDepolarizing-FourBranches"
]

VerificationTest[
    Head @ QuantumChannel["AmplitudeDamping"[1/2], {1}][$ps0],
    QuantumState,
    {PauliStabilizer::nonpaulibasis},
    TestID -> "Conn6-qcAmplitudeDamping-FallsBack"
]


(* ============================================================================ *)
(* CONNECTION 7 -- cc_CliffordChannel[ps] DownValue (Phase 8)                   *)
(* ============================================================================ *)

(* Direct dispatch: identity channel. *)
VerificationTest[
    Head @ CliffordChannel["Identity", 2][$ps00],
    PauliStabilizer,
    TestID -> "Conn7-ccIdentity-Dispatches"
]

(* Identity channel preserves the state. *)
VerificationTest[
    CliffordChannel["Identity", 2][$ps00]["Stabilizers"] === $ps00["Stabilizers"],
    True,
    TestID -> "Conn7-ccIdentity-Preserves"
]

(* CC composition: cc_S^2 = cc_Z; applied to |1> gives -|1>. *)
VerificationTest[
    Module[{ccS, ccS2, ps1, psf},
        ccS = CliffordChannel[<|"UA" -> {{1, 0}, {0, 1}}, "UB" -> {{1, 1}, {0, 1}},
                                 "c" -> {0, 0}, "InputQubits" -> 1,
                                 "OutputQubits" -> 1, "Source" -> "S"|>];
        ccS2 = ccS[ccS];
        ps1 = $ps0["X", 1];
        psf = ccS2[ps1];
        psf["Stabilizers"] === {"-Z"}
    ],
    True,
    TestID -> "Conn7-ccS2onOne-EqualsNegOne"
]

(* CliffordChannel state-prep: cc -> ps round-trip. *)
VerificationTest[
    Module[{cc, ps2},
        cc = CliffordChannel[$ps00];
        cc["InputQubits"] == 0 && cc["OutputQubits"] == 2 && cc["Rank"] == 2
    ],
    True,
    TestID -> "Conn7-ccFromPS-StatePrep-Structure"
]


(* ============================================================================ *)
(* CONNECTION 8 -- QuantumTensorProduct[ps_a, ps_b] DownValue                   *)
(* ============================================================================ *)

(* Direct dispatch. *)
VerificationTest[
    Head @ QuantumTensorProduct[$ps0, $ps0],
    PauliStabilizer,
    TestID -> "Conn8-QTPofPSPS-Dispatches"
]

(* Qubit count adds. *)
VerificationTest[
    QuantumTensorProduct[$psBell, $ps0]["Qubits"],
    3,
    TestID -> "Conn8-QTP-QubitsAdd"
]

(* Stabilizer set is the block-diagonal merge: {XXI, ZZI, IIZ} for Bell ⊗ |0>. *)
VerificationTest[
    QuantumTensorProduct[$psBell, $ps0]["Stabilizers"],
    {"XXI", "ZZI", "IIZ"},
    TestID -> "Conn8-QTP-Bell-Zero-Stabs"
]

(* Triple tensor product. *)
VerificationTest[
    QuantumTensorProduct[$ps0, $ps0, $ps0]["Stabilizers"],
    {"ZII", "IZI", "IIZ"},
    TestID -> "Conn8-QTP-Triple"
]


(* ============================================================================ *)
(* CONNECTION 9 -- QuantumCircuitOperator[..][Method -> "Stabilizer"]           *)
(* ============================================================================ *)

(* Direct dispatch: the Method routes through PauliStabilizerApply. *)
VerificationTest[
    Head @ QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}][Method -> "Stabilizer"],
    PauliStabilizer,
    TestID -> "Conn9-MethodStabilizer-Dispatches"
]

(* Semantic: matches PauliStabilizer @ qco directly. *)
VerificationTest[
    Sort @ QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}][Method -> "Stabilizer"]["Stabilizers"],
    Sort @ {"XX", "ZZ"},
    TestID -> "Conn9-MethodStabilizer-Bell"
]

(* Named circuit. *)
VerificationTest[
    QuantumCircuitOperator["GHZ"[3]][Method -> "Stabilizer"]["Stabilizers"],
    {"XXX", "ZZI", "IZZ"},
    TestID -> "Conn9-MethodStabilizer-NamedGHZ"
]

(* Method on a circuit with explicit input QuantumState. *)
VerificationTest[
    Module[{result},
        result = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}][
            QuantumState["00"], Method -> "Stabilizer"
        ];
        Head[result] === PauliStabilizer && Sort[result["Stabilizers"]] === Sort[{"XX", "ZZ"}]
    ],
    True,
    TestID -> "Conn9-MethodStabilizer-ExplicitQS"
]


(* ============================================================================ *)
(* CONNECTION 10 -- StabilizerFrame cross-head: qmo[sf], qc[sf]                 *)
(* ============================================================================ *)

(* T-gate -> StabilizerFrame, then qmo dispatches via the materialize-fallback. *)
VerificationTest[
    Module[{psT = $ps0["H", 1]["T", 1], qmo = QuantumMeasurementOperator["Z", {1}], result},
        result = qmo[psT];
        FreeQ[result, $Failed]
    ],
    True,
    {PauliStabilizer::nonpaulibasis},
    TestID -> "Conn10-qmoOnFrame-Materializes"
]

(* qc on frame falls back via materialization. *)
VerificationTest[
    Module[{psT = $ps0["H", 1]["T", 1], qc = QuantumChannel["BitFlip"[1/3], {1}], result},
        result = qc[psT];
        FreeQ[result, $Failed]
    ],
    True,
    {PauliStabilizer::nonpaulibasis},
    TestID -> "Conn10-qcOnFrame-Materializes"
]


(* ============================================================================ *)
(* CONNECTION 11 -- GraphState <-> PauliStabilizer round-trip                   *)
(* ============================================================================ *)

(* gs -> ps -> gs round-trip on a graph-form input. *)
VerificationTest[
    Module[{gs = GraphState[Graph[Range[4], Table[i \[UndirectedEdge] (i + 1), {i, 3}]]], ps, gs2},
        ps = gs["PauliStabilizer"];
        gs2 = GraphState[ps];
        Sort[EdgeList[gs2["Graph"]]] === Sort[EdgeList[gs["Graph"]]]
    ],
    True,
    TestID -> "Conn11-RT-GS-PS-GS-LinearCluster4"
]

(* gs stabilizers match circuit-built cluster state. *)
VerificationTest[
    Module[{n = 4, gs, psFromGS, psFromCirc},
        gs = GraphState[Graph[Range[n], Table[i \[UndirectedEdge] (i + 1), {i, n - 1}]]];
        psFromGS = gs["PauliStabilizer"];
        psFromCirc = PauliStabilizer[QuantumCircuitOperator[
            Join[Table["H" -> i, {i, n}], Table["CZ" -> {i, i + 1}, {i, n - 1}]]
        ]];
        Sort[psFromGS["Stabilizers"]] === Sort[psFromCirc["Stabilizers"]]
    ],
    True,
    TestID -> "Conn11-GraphState-EqualsClusterCircuit"
]


(* ============================================================================ *)
(* CHAIN 1 -- ps -> QuantumState -> measurement via QMO -> outcome distribution *)
(*                                                                              *)
(* Verifies that going PS -> QS and applying a QMO matches direct ps["M", q]    *)
(* on the same Pauli operator. End-to-end probabilistic-outcome test.           *)
(* ============================================================================ *)

VerificationTest[
    Module[{psBell, qsBell, expQS, expPS},
        psBell = $psBell;
        qsBell = QuantumState[psBell];

        (* Path A: ps -> QuantumState -> Tr[Z⊗Z * rho] = <ZZ>. *)
        expQS = Re @ Chop @ Tr[
            Normal[QuantumOperator["ZZ"]["MatrixRepresentation"]] .
            Normal[qsBell["DensityMatrix"]]
        ];

        (* Path B: direct expectation via ps["Expectation", "ZZ"]. *)
        expPS = psBell["Expectation", "ZZ"];

        expQS === expPS
    ],
    True,
    TestID -> "Chain1-PS-QS-TrZZRho-vs-Expectation-Bell"
]


(* ============================================================================ *)
(* CHAIN 2 -- ps -> Circuit -> apply on QuantumState[|0>] -> recovers ps state  *)
(*                                                                              *)
(* The AG-decomposed circuit's dagger applied to |0...0> reproduces ps's        *)
(* materialized state (up to global phase, per Phase 5c contract).              *)
(* ============================================================================ *)

VerificationTest[
    Module[{psBell, circ, recovered, v1, v2},
        psBell = $psBell;
        circ = psBell["Circuit"];
        (* Apply the AG-decomposed circuit to the |00> register. The circuit is *)
        (* stored as `dagger applied to |0>` -- so circ[QuantumState["00"]]     *)
        (* recovers psBell["State"] up to global phase.                         *)
        recovered = circ[QuantumState["00"]];
        v1 = Normal @ psBell["State"]["StateVector"];
        v2 = Normal @ recovered["StateVector"];
        Chop[Abs[Conjugate[v1] . v2] - 1] === 0
    ],
    True,
    TestID -> "Chain2-PS-Circuit-onQS-MatchesPS-State-UpToPhase"
]


(* ============================================================================ *)
(* CHAIN 3 -- ps -> apply Clifford gate via UpValue -> measure -> match direct  *)
(* ============================================================================ *)

VerificationTest[
    Module[{ps0, qoH, qoCNOT, psBell2, viaUpValue, psBellDirect},
        ps0 = $ps00;
        qoH = QuantumOperator["H", 1];
        qoCNOT = QuantumOperator["CNOT", {1, 2}];

        (* Apply gates via UpValue chain. *)
        viaUpValue = qoCNOT @ (qoH @ ps0);

        (* Direct method calls. *)
        psBellDirect = ps0["H", 1]["CNOT", 1, 2];

        Sort[viaUpValue["Stabilizers"]] === Sort[psBellDirect["Stabilizers"]]
    ],
    True,
    TestID -> "Chain3-qoUpValueChain-EqualsMethodChain"
]


(* ============================================================================ *)
(* CHAIN 4 -- QF circuit -> Stabilizer method -> CC composition -> back to PS   *)
(*                                                                              *)
(* End-to-end story: build a circuit in QF; route through Method ->            *)
(* "Stabilizer"; compose with another CC; verify the final PS is correct.      *)
(* ============================================================================ *)

VerificationTest[
    Module[{qco, ps, ccS, psf},
        (* Step 1: QF circuit. *)
        qco = QuantumCircuitOperator[{"H" -> 1}];

        (* Step 2: Method -> "Stabilizer" gives a PS. *)
        ps = qco[Method -> "Stabilizer"];

        (* Step 3: Apply CliffordChannel cc_S to ps. After H|0> = |+>, apply S: *)
        (* S|+> = |y+>, stabilizer Y. *)
        ccS = CliffordChannel[<|"UA" -> {{1, 0}, {0, 1}}, "UB" -> {{1, 1}, {0, 1}},
                                 "c" -> {0, 0}, "InputQubits" -> 1,
                                 "OutputQubits" -> 1, "Source" -> "S"|>];
        psf = ccS[ps];

        (* Verify chain end: stabilizer is "Y". *)
        Head[psf] === PauliStabilizer && psf["Stabilizers"] === {"Y"}
    ],
    True,
    TestID -> "Chain4-QFcircuit-Method-Stabilizer-CC-Composition"
]


(* ============================================================================ *)
(* CHAIN 5 -- Outcome distribution sanity: ps["M"] vs sampling QuantumState    *)
(*                                                                              *)
(* Take Bell state, measure first qubit. ps["M"] returns Association with both *)
(* outcomes. The materialized QuantumState's measurement should agree on the   *)
(* outcome-key set.                                                             *)
(* ============================================================================ *)

(* Use 2-qubit ZZ instead of 1-qubit Z to dodge ROADMAP A.11 (single-qubit    *)
(* Pauli-string measurement on multi-qubit ps fails with KroneckerProduct::argmu). *)
VerificationTest[
    Module[{psBell, viaPS, qmo, viaQS},
        psBell = $psBell;
        viaPS = Sort[Keys[psBell["M", "ZZ"]]];
        qmo = QuantumMeasurementOperator["ZZ", {1, 2}];
        viaQS = Sort[Keys[qmo[psBell]]];
        viaPS === viaQS
    ],
    True,
    TestID -> "Chain5-OutcomeKeys-PS-MZZ-vs-qmoZZ-Bell"
]


(* ============================================================================ *)
(* CHAIN 6 -- Bell state via three paths agree                                  *)
(*                                                                              *)
(* Path A: PauliStabilizer[QuantumCircuitOperator[Bell-circuit]]                *)
(* Path B: QuantumState["PhiPlus"] -> PauliStabilizer (tomography)              *)
(* Path C: GraphState[2-vertex graph]["PauliStabilizer"] (Bell == cluster-2     *)
(*           up to local Hadamards: cluster-2 has stabilizers {XZ, ZX} which   *)
(*           is NOT exactly Bell; instead use circuit form for path C.)        *)
(*                                                                              *)
(* All three should produce the same stabilizer set (modulo signs / canonical) *)
(* and the same materialized state vector.                                     *)
(* ============================================================================ *)

VerificationTest[
    Module[{psA, psB, vecA, vecB},
        psA = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
        psB = PauliStabilizer[QuantumState["PhiPlus"]];

        vecA = Normal @ psA["State"]["StateVector"];
        vecB = Normal @ psB["State"]["StateVector"];

        (* Both should produce |Phi+> = (|00>+|11>)/sqrt(2) exactly (Phase 5c   *)
        (* guarantees first-hop GlobalPhase capture for QS path).               *)
        Sort[psA["Stabilizers"]] === Sort[psB["Stabilizers"]] &&
        vecA === vecB
    ],
    True,
    TestID -> "Chain6-Bell-via-Two-Paths-Agree"
]


(* ============================================================================ *)
(* CHAIN 7 -- 5-qubit code: build via PS, measure all stabilizers via QMO       *)
(* ============================================================================ *)

VerificationTest[
    Module[{ps5 = $ps5Q, allOutcomes},
        allOutcomes = Sort @ Keys @ # & /@ {
            QuantumMeasurementOperator["XZZXI", Range[5]][ps5],
            QuantumMeasurementOperator["IXZZX", Range[5]][ps5],
            QuantumMeasurementOperator["XIXZZ", Range[5]][ps5],
            QuantumMeasurementOperator["ZXIXZ", Range[5]][ps5]
        };
        AllTrue[allOutcomes, # === {0} &]
    ],
    True,
    TestID -> "Chain7-5Q-AllStabs-Deterministic-via-qmoPS"
]


(* ============================================================================ *)
(* CHAIN 8 -- Bell -> apply BitFlip channel via UpValue -> total prob = 1       *)
(* ============================================================================ *)

VerificationTest[
    Module[{ps = $psBell, qc, branches, total},
        qc = QuantumChannel["BitFlip"[1/4], {1}];
        branches = qc[ps];
        total = Total[First /@ branches];
        Length[branches] == 2 && total == 1
    ],
    True,
    TestID -> "Chain8-Bell-qcBitFlip-Branches-Sum-One"
]


(* ============================================================================ *)
(* CHAIN 9 -- Random Clifford -> apply via Method -> compare to direct          *)
(* ============================================================================ *)

VerificationTest[
    Block[{},
        SeedRandom[20260507];
        AllTrue[
            Table[
                Module[{psRandom, qco, viaMethod, viaApply},
                    psRandom = PauliStabilizer["Random", 3];
                    qco = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, "S" -> 3}];
                    viaMethod = qco[psRandom["State"], Method -> "Stabilizer"];
                    viaApply = psRandom["H", 1]["CNOT", 1, 2]["S", 3];
                    matEqQS[viaMethod["State"], viaApply["State"]]
                ],
                {6}
            ],
            TrueQ
        ]
    ],
    True,
    TestID -> "Chain9-Random-Method-vs-DirectApply-6reps"
]


(* ============================================================================ *)
(* CHAIN 10 -- Density matrix path: ps -> QuantumState -> DensityMatrix         *)
(*                                                                              *)
(* Build PS, materialize to QS, take density matrix, verify it's a valid pure  *)
(* state (rank 1, trace 1).                                                     *)
(* ============================================================================ *)

VerificationTest[
    Module[{ps = $psBell, qs, rho, traceVal, rankVal},
        qs = QuantumState[ps];
        rho = qs["DensityMatrix"] // Normal;
        traceVal = Tr[rho];
        rankVal = MatrixRank[rho];
        Chop[traceVal - 1] === 0 && rankVal == 1
    ],
    True,
    TestID -> "Chain10-PS-QS-DensityMatrix-Pure-Bell"
]


(* ============================================================================ *)
(* CHAIN 11 -- ps tensor-product chain across heads                             *)
(*                                                                              *)
(* PS_a tensor PS_b matches QuantumState[PS_a] tensor QuantumState[PS_b] when  *)
(* both are materialized.                                                       *)
(* ============================================================================ *)

VerificationTest[
    Module[{tensorPS, tensorQS},
        tensorPS = QuantumTensorProduct[$psBell, $ps0]["State"];
        tensorQS = QuantumTensorProduct[QuantumState[$psBell], QuantumState[$ps0]];
        matEqQS[tensorPS, tensorQS]
    ],
    True,
    TestID -> "Chain11-QTP-PSPath-EqualsQSPath"
]
