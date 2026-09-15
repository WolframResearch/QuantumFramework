(* ::Package:: *)

(* ============================================================================
   Tests/QEC/ErrorCorrection.wlt

   Steane error correction (Got26 sec. 12.3), and the slot it fills in the
   fault-tolerant Pauli measurement of Measurement.wl.

   The load-bearing test is QEC-EC-additivity-one-fault-one-data-error, which is
   sec. 12.3.3's claim: the classical error is e + f + g, so correcting by it
   leaves f + g, and one ancilla fault can leave only one data error.  That is
   what makes Steane EC need no repetition, and it is why this gadget and not
   Shor EC fills the measurement's slot.

   Its companion is QEC-EC-preparation-is-where-fault-tolerance-is-missing.  The
   ancilla is prepared by the code's own non-fault-tolerant encoder, so a fault
   there does reach several data qubits.  Got26 says so plainly and defers the
   fix to chapter 13; the two tests together turn that from a caveat into a
   measured split.

   The ancilla states are checked against the engine, not asserted: |0_L> and
   |+_L> have to stabilise every generator and be +1 eigenstates of Z-bar and
   X-bar respectively.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

QECClearCache[];

qecScope = "Wolfram`QuantumFramework`QEC`PackageScope`";

qecFrame      = Symbol[qecScope <> "framePropagate"];
qecEngGates   = Symbol[qecScope <> "instructionEngineGates"];
qecEngInstr   = Symbol[qecScope <> "engineGateInstructions"];
qecApply      = Symbol[qecScope <> "applyGates"];
qecCodeData   = Symbol[qecScope <> "codeData"];
qecWires      = Symbol[qecScope <> "codeLogicalWires"];
qecZeroInstr  = Symbol[qecScope <> "codeEncodedZeroInstructions"];
qecPlusInstr  = Symbol[qecScope <> "codeEncodedPlusInstructions"];
qecRegions    = Symbol[qecScope <> "steaneRegions"];
qecOneQubit   = Symbol[qecScope <> "$circuitOneQubitOps"];
qecTwoQubit   = Symbol[qecScope <> "$circuitTwoQubitOps"];

steane = QECCode["SteaneCode"];
bitFlip = QECCode["BitFlipCode"];
five = QECCode["5QubitCode"];

ec = QECErrorCorrection[steane];

(* Run a preparation as an actual stabilizer state, to check what it prepares. *)
qecPrepared[instrs_, n_] := qecApply[PauliStabilizer[n], qecEngGates[instrs]]


(* ============================================================================
   The new instruction vocabulary
   ============================================================================ *)

(* The encoder speaks SWAP and Pauli gates, so the instruction language had to
   grow them.  SWAP exchanges both halves of the frame. *)
VerificationTest[
    qecFrame[{{"SWAP", 1, 2}}, 2, {{0, 1, {1, 0}}, {0, 2, {0, 1}}}]["Frame"],
    {{0, 1}, {1, 0}},
    TestID -> "QEC-EC-swap-exchanges-the-frame"
]

(* A deterministic Pauli conjugates every Pauli to itself up to a sign, and signs
   are dropped, so X and Z move no frame bit -- but they are still instructions,
   because a gate the noise model cannot see is a fault location silently lost. *)
VerificationTest[
    {qecFrame[{{"X", 1}, {"Z", 1}}, 1, {{0, 1, {1, 1}}}]["Frame"],
     MemberQ[qecOneQubit, "X"], MemberQ[qecOneQubit, "Z"],
     MemberQ[qecTwoQubit, "SWAP"]},
    {{{1}, {1}}, True, True, True},
    TestID -> "QEC-EC-pauli-gates-are-locations-but-not-frame-moves"
]

(* The two converters are inverses, which is what keeps the encoder and the
   instruction list from drifting apart. *)
VerificationTest[
    With[{instr = {{"H", 1}, {"CNOT", 1, 2}, {"CZ", 2, 3}, {"SWAP", 1, 3}, {"X", 2}, {"Z", 3}}},
        qecEngInstr[qecEngGates[instr], 0] === instr
    ],
    True,
    TestID -> "QEC-EC-engine-and-instruction-forms-round-trip"
]


(* ============================================================================
   The encoded ancilla states
   ============================================================================ *)

(* Derived, not assumed: the wire whose |+> input makes the encoded state a +1
   eigenstate of logical X.  It lands on the last k, which is the standard form's
   convention -- recorded here rather than relied on. *)
VerificationTest[
    {qecWires[qecCodeData[steane]], qecWires[qecCodeData[five]], qecWires[qecCodeData[bitFlip]]},
    {{7}, {5}, {3}},
    TestID -> "QEC-EC-logical-wires-are-the-last-k"
]

(* |0_L>: stabilised by every generator, and a +1 eigenstate of Z-bar. *)
VerificationTest[
    With[{s = qecPrepared[qecZeroInstr[qecCodeData[steane], 0], 7]},
        {AllTrue[steane["Generators"], s["Expectation", #] === 1 &],
         s["Expectation", First[steane["LogicalZ"]]]}
    ],
    {True, 1},
    TestID -> "QEC-EC-prepares-encoded-zero"
]

(* |+_L>: the same generators, and a +1 eigenstate of X-bar instead. *)
VerificationTest[
    With[{s = qecPrepared[qecPlusInstr[qecCodeData[steane], 0], 7]},
        {AllTrue[steane["Generators"], s["Expectation", #] === 1 &],
         s["Expectation", First[steane["LogicalX"]]]}
    ],
    {True, 1},
    TestID -> "QEC-EC-prepares-encoded-plus"
]


(* ============================================================================
   Shape
   ============================================================================ *)

(* Two halves, each an encoded block of n ancillas and n measurements. *)
VerificationTest[
    {ec["DataQubits"], ec["AncillaQubits"], ec["Qubits"], ec["Measurements"]},
    {7, Range[8, 14], 14, 14},
    TestID -> "QEC-EC-two-halves-on-one-reused-ancilla-block"
]

(* Bit-flip half: data is the control, so bit flips copy into the ancilla
   (fig. 12.8a).  Phase half: ancilla is the control and a transversal Hadamard
   precedes the readout (fig. 12.8b).

   Only the tail of each half is the interaction -- everything before it is the
   encoder, which has CNOTs and Hadamards of its own on the ancilla block.  Taking
   the last 2n and 3n instructions is how the two are told apart here; the gadget
   itself does it by region. *)
VerificationTest[
    {Take[ec["BitFlipInstructions"], -14],
     Take[ec["PhaseInstructions"], -21]},
    {Join[Table[{"CNOT", q, q + 7}, {q, 7}], Table[{"M", q + 7}, {q, 7}]],
     Join[Table[{"CNOT", q + 7, q}, {q, 7}], Table[{"H", q + 7}, {q, 7}],
          Table[{"M", q + 7}, {q, 7}]]},
    TestID -> "QEC-EC-the-two-transversal-CNOTs-run-opposite-ways"
]

(* The data-ancilla interaction is one CNOT per position and nothing else. *)
VerificationTest[
    {ec["TransversalQ"], QECErrorCorrection[bitFlip]["TransversalQ"],
     QECErrorCorrection[steane, 7, "Order" -> "PhaseFirst"]["TransversalQ"]},
    {True, True, True},
    TestID -> "QEC-EC-data-ancilla-interaction-is-transversal"
]

(* Section 12.3.3: no repetition at all, unlike Shor EC and unlike the cat
   measurement.  That is the property that makes this the right sub-gadget. *)
VerificationTest[
    ec["RepetitionsNeeded"],
    1,
    TestID -> "QEC-EC-needs-no-repetition"
]

(* The regions are a partition: every fault index is preparation or interaction,
   never both and never neither. *)
VerificationTest[
    With[{r = ec["Regions"]},
        {Sort[Join[r["Preparation"], r["Interaction"]]] === Range[0, ec["InstructionCount"]],
         Intersection[r["Preparation"], r["Interaction"]]}
    ],
    {True, {}},
    TestID -> "QEC-EC-regions-partition-the-gadget"
]


(* ============================================================================
   The property the gadget exists for, and the one it still lacks
   ============================================================================ *)

(* THE test.  Additivity: a single fault in the interaction leaves at most one
   data error, for either order and either CSS code. *)
VerificationTest[
    {ec["DataWeights"],
     QECErrorCorrection[steane, 7, "Order" -> "PhaseFirst"]["DataWeights"],
     QECErrorCorrection[bitFlip]["DataWeights"]},
    {{0, 1}, {0, 1}, {0, 1}},
    TestID -> "QEC-EC-additivity-one-fault-one-data-error"
]

(* And its companion.  The ancilla preparation is the code's non-fault-tolerant
   encoder, so a fault there does reach several data qubits.  This is the open
   assumption, measured rather than caveated. *)
VerificationTest[
    {Max[ec["PreparationDataWeights"]], Max[QECErrorCorrection[bitFlip]["PreparationDataWeights"]]},
    {4, 3},
    TestID -> "QEC-EC-preparation-is-where-fault-tolerance-is-missing"
]

VerificationTest[
    Length[ec["OpenAssumptions"]] === 1 &&
        StringContainsQ[First[ec["OpenAssumptions"]], "not fault"],
    True,
    TestID -> "QEC-EC-names-its-open-assumption"
]


(* ============================================================================
   Filling the measurement gadget's slot
   ============================================================================ *)

(* Without the sub-gadget the measurement handles ancilla faults only; with it,
   Theorem 12.1's structure is complete. *)
VerificationTest[
    With[{
        m0 = QECPauliMeasurement[steane, First[steane["LogicalZ"]]],
        m1 = QECPauliMeasurement[steane, First[steane["LogicalZ"]], "ErrorCorrection" -> "Steane"]
    },
        {m0["MeasurementCorrectQ"], m1["MeasurementCorrectQ"],
         m0["Qubits"], m1["Qubits"], m1["DataWeights"]}
    ],
    {False, True, 11, 18, {0, 1}},
    TestID -> "QEC-EC-fills-the-measurement-slot"
]

(* The measurement inherits the sub-gadget's open assumption rather than hiding it
   behind that True. *)
VerificationTest[
    QECPauliMeasurement[steane, First[steane["LogicalZ"]], "ErrorCorrection" -> "Steane"]["OpenAssumptions"]
        === ec["OpenAssumptions"],
    True,
    TestID -> "QEC-EC-open-assumptions-are-inherited"
]


(* ============================================================================
   Refusals
   ============================================================================ *)

(* Steane EC runs on transversal CNOT being the logical CNOT, which is a CSS
   property.  The five-qubit code is not CSS and is refused rather than silently
   producing a gadget that does not measure what it claims. *)
VerificationTest[
    QECErrorCorrection[five],
    $Failed,
    {QECErrorCorrection::css},
    TestID -> "QEC-EC-needs-a-CSS-code"
]

VerificationTest[
    QECPauliMeasurement[five, "ZZZZZ", "ErrorCorrection" -> "Steane"],
    $Failed,
    {QECPauliMeasurement::eccss},
    TestID -> "QEC-EC-the-measurement-refuses-a-non-CSS-code-too"
]

VerificationTest[
    QECErrorCorrection[steane, 7, "Order" -> "Sideways"],
    $Failed,
    {QECErrorCorrection::order},
    TestID -> "QEC-EC-order-must-name-a-half"
]

VerificationTest[
    ec["Nonsense"],
    Missing["NotFound", "Nonsense"],
    {QECErrorCorrection::noprop},
    TestID -> "QEC-EC-unknown-property-is-refused"
]
