(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Measurement.wlt

   Fault-tolerant measurement of a Pauli through cat states (Got26 sec. 12.1.2,
   12.1.4, Theorem 12.1).

   The load-bearing test in this file is QEC-Measure-one-fault-one-data-error.
   Everything else is structure; that one is the property the construction exists
   for, and it is stated against the bare-ancilla circuit measuring the SAME
   operator so the comparison is like for like.

   Two corrections make that number mean something, and both are tested on their
   own so neither can be mistaken for bookkeeping:

     - the residual is read modulo P, because X^(tensor m) stabilises the cat and
       a pattern and its complement are one physical error.  The accepted weights
       are complement-symmetric, which is the fingerprint of that quotient being
       the right one rather than a convenient one.
     - only accepted runs count.  Drop the herald filter and weight two comes
       back, so the cat checks are visibly load-bearing.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

QECClearCache[];

qecScope = "Wolfram`QuantumFramework`QEC`PackageScope`";

qecFrame    = Symbol[qecScope <> "framePropagate"];
qecWeights  = Symbol[qecScope <> "pauliMeasureDataWeights"];
qecOutcome  = Symbol[qecScope <> "pauliMeasureOutcome"];
qecGenInstr = Symbol[qecScope <> "generatorInstructions"];
qecCatMeas  = Symbol[qecScope <> "catMeasureInstructions"];
qecReps     = Symbol[qecScope <> "pauliMeasureRepetitions"];

five  = QECCode["5QubitCode"];
steane = QECCode["SteaneCode"];

(* The five-qubit logical Z, weight five, three repetitions. *)
mz = QECPauliMeasurement[five, "ZZZZZ"];

(* One repetition of the first generator, to compare against the bare-ancilla
   circuit that measures exactly the same operator. *)
mg = QECPauliMeasurement[five, "XZZXI", 1];

qecParities[m_, faults_] := Mod[
    Total /@ Partition[qecFrame[m["Instructions"], m["Qubits"], faults]["Record"], m["Weight"]],
    2
]


(* ============================================================================
   Shape
   ============================================================================ *)

(* One ancilla qubit per single-qubit Pauli making up P, which is what makes the
   controlled-P transversal in the first place. *)
VerificationTest[
    {mz["Weight"], mz["CatQubits"], mz["CheckQubit"], mz["Qubits"]},
    {5, {6, 7, 8, 9, 10}, 11, 11},
    TestID -> "QEC-Measure-one-cat-qubit-per-letter"
]

(* A code correcting t errors gets 2t+1 repetitions (Theorem 12.1).  Distance one
   gets one, which is the honest answer rather than a special case. *)
VerificationTest[
    {qecReps[1], qecReps[2], qecReps[3], qecReps[5], qecReps[7]},
    {1, 1, 3, 5, 7},
    TestID -> "QEC-Measure-repetitions-are-two-t-plus-one"
]

VerificationTest[
    {QECPauliMeasurement[five, "ZZZZZ"]["Repetitions"],
     QECPauliMeasurement[QECCode["BitFlipCode"], "ZZI"]["Repetitions"]},
    {3, 1},
    TestID -> "QEC-Measure-repetitions-follow-the-code-distance"
]

(* Each repetition prepares a fresh cat, so the readout count scales with both. *)
VerificationTest[
    {mz["Measurements"], mz["Heralds"], mg["Measurements"], mg["Heralds"]},
    {15, 6, 4, 1},
    TestID -> "QEC-Measure-a-fresh-cat-per-repetition"
]

(* The readout is the Hadamard transform of eqs. 12.1-12.3, not a disentangling
   circuit: one H on each cat qubit, then measure them all. *)
VerificationTest[
    Take[qecCatMeas[QECPauliVector["ZZZZZ"], 5, {6, 7, 8, 9, 10}], -10],
    Join[Table[{"H", q}, {q, 6, 10}], Table[{"M", q}, {q, 6, 10}]],
    TestID -> "QEC-Measure-readout-is-a-hadamard-transform"
]

(* The controlled-P is the same U, CZ, U^-1 sandwich generatorInstructions uses,
   with CZ where that has CNOT -- shared rather than re-derived, so the two
   constructions cannot drift apart on what a letter's rotation is. *)
VerificationTest[
    Cases[qecCatMeas[QECPauliVector["XZZXI"], 5, {6, 7, 8, 9}], {"H" | "V" | "Vdg", q_} /; q <= 5],
    Cases[qecGenInstr[First[five["CheckMatrix"]], 5, 6], {"H" | "V" | "Vdg", q_} /; q <= 5],
    TestID -> "QEC-Measure-rotations-are-shared-with-the-extraction-circuit"
]

(* Structural transversality: no cat qubit shares a two-qubit gate with more than
   one data qubit. *)
VerificationTest[
    {mz["TransversalQ"], mg["TransversalQ"],
     QECPauliMeasurement[steane, First[steane["LogicalZ"]]]["TransversalQ"]},
    {True, True, True},
    TestID -> "QEC-Measure-is-transversal"
]


(* ============================================================================
   The property the construction exists for
   ============================================================================ *)

(* THE test.  Conditioned on the cat checks accepting, and modulo P, every single
   fault the gadget admits leaves at most one data error. *)
VerificationTest[
    {mz["DataWeights"], mg["DataWeights"],
     QECPauliMeasurement[steane, First[steane["LogicalZ"]]]["DataWeights"]},
    {{0, 1}, {0, 1}, {0, 1}},
    TestID -> "QEC-Measure-one-fault-one-data-error"
]

(* The same operator measured the old way.  XZZXI with a bare shared ancilla lets
   one fault reach four data qubits -- the hook error, priced in the same units. *)
VerificationTest[
    Max @ qecWeights[
        qecGenInstr[First[five["CheckMatrix"]], 5, 6], 6, 5, ConstantArray[0, 10]
    ],
    4,
    TestID -> "QEC-Measure-bare-ancilla-reaches-four-data-qubits"
]

(* The quotient by P is forced, not chosen: X^(tensor m) stabilises the cat, so the
   accepted weights must come in complementary pairs w and m - w.  They do, which
   is what says the quotient is the right one. *)
VerificationTest[
    With[{w = Sort @ qecWeights[mg["Instructions"], mg["Qubits"], 5, ConstantArray[0, 10]]},
        {w, Sort[mg["Weight"] - w] === w}
    ],
    {{0, 1, 3, 4}, True},
    TestID -> "QEC-Measure-accepted-weights-are-complement-symmetric"
]

(* And the cat checks are load-bearing rather than decorative: stop discarding the
   runs they reject and weight two reappears. *)
VerificationTest[
    Max @ DeleteDuplicates @ Flatten @ Table[
        With[{
            fr = qecFrame[mg["Instructions"], mg["Qubits"], {{i, q, pauli}}]["Frame"],
            pv = Take[QECPauliVector["XZZXI"], 10]
        },
            With[{row = Join[Take[fr[[1]], 5], Take[fr[[2]], 5]]},
                Min[
                    Total[Max /@ Transpose[{Take[row, 5], Take[row, -5]}]],
                    Total[Max /@ Transpose[{Take[BitXor[row, pv], 5], Take[BitXor[row, pv], -5]}]]
                ]
            ]
        ],
        {i, 0, mg["InstructionCount"]}, {q, mg["Qubits"]}, {pauli, {{1, 0}, {0, 1}, {1, 1}}}
    ],
    2,
    TestID -> "QEC-Measure-the-cat-checks-are-load-bearing"
]


(* ============================================================================
   Reading the answer
   ============================================================================ *)

(* Each repetition reports the parity of its m record bits; the gadget reports the
   majority.  An odd repetition count means a majority always exists. *)
VerificationTest[
    {qecOutcome[Flatten[{{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}}], 5, 3],
     qecOutcome[Flatten[{{1, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}}], 5, 3],
     qecOutcome[Flatten[{{1, 0, 0, 0, 0}, {1, 0, 0, 0, 0}, {0, 0, 0, 0, 0}}], 5, 3],
     qecOutcome[Flatten[{{1, 1, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}}], 5, 3]},
    {0, 0, 1, 0},
    TestID -> "QEC-Measure-majority-over-repetition-parities"
]

VerificationTest[
    {qecParities[mz, {}], qecOutcome[qecFrame[mz["Instructions"], mz["Qubits"], {}]["Record"], 5, 3]},
    {{0, 0, 0}, 0},
    TestID -> "QEC-Measure-no-faults-reads-the-plus-one-eigenvalue"
]

(* A fault inside one repetition corrupts that repetition and no other, which is
   the whole point of a fresh cat each time -- so the majority survives it. *)
VerificationTest[
    With[{i = FirstPosition[mz["Instructions"], {"M", 6}][[1]]},
        {qecParities[mz, {{i - 1, 6, {1, 0}}}],
         qecOutcome[qecFrame[mz["Instructions"], mz["Qubits"], {{i - 1, 6, {1, 0}}}]["Record"], 5, 3]}
    ],
    {{1, 0, 0}, 0},
    TestID -> "QEC-Measure-one-ancilla-fault-loses-one-repetition-only"
]


(* ============================================================================
   The failure mode repetition cannot fix (sec. 12.1.4)
   ============================================================================ *)

(* A data error that ANTICOMMUTES with P flips every repetition alike, so the
   majority is confidently wrong -- "no matter how many times we repeat it".  This
   is why Theorem 12.1 needs FTEC interspersed between repetitions and not just
   more of them. *)
VerificationTest[
    {qecParities[mz, {{0, 1, {1, 0}}}],
     qecOutcome[qecFrame[mz["Instructions"], mz["Qubits"], {{0, 1, {1, 0}}}]["Record"], 5, 3]},
    {{1, 1, 1}, 1},
    TestID -> "QEC-Measure-anticommuting-data-error-defeats-every-repetition"
]

(* A data error that commutes with P is invisible to it, as it should be. *)
VerificationTest[
    qecParities[mz, {{0, 1, {0, 1}}}],
    {0, 0, 0},
    TestID -> "QEC-Measure-commuting-data-error-is-invisible"
]

(* So the gadget does not yet satisfy the measurement correctness property, and
   says so rather than claiming it.  The slot is there; Steane EC fills it. *)
VerificationTest[
    {mz["MeasurementCorrectQ"], mz["ErrorCorrection"]},
    {False, None},
    TestID -> "QEC-Measure-MCP-is-false-while-the-EC-slot-is-empty"
]


(* ============================================================================
   Refusals
   ============================================================================ *)

VerificationTest[
    QECPauliMeasurement[five, "IIIII"],
    $Failed,
    {QECPauliMeasurement::identity},
    TestID -> "QEC-Measure-the-identity-has-no-eigenvalue"
]

VerificationTest[
    QECPauliMeasurement[five, "ZZZZZ", 0],
    $Failed,
    {QECPauliMeasurement::reps},
    TestID -> "QEC-Measure-repetitions-must-be-positive"
]

(* Steane EC exists now (ErrorCorrection.wl), but it is CSS-only, and the 5-qubit
   code is not CSS -- so the refusal here is the CSS one, not "not built yet". *)
VerificationTest[
    QECPauliMeasurement[five, "ZZZZZ", "ErrorCorrection" -> "Steane"],
    $Failed,
    {QECPauliMeasurement::eccss},
    TestID -> "QEC-Measure-Steane-error-correction-needs-a-CSS-code"
]

(* And a gadget that is not built at all is refused by name.  Shor EC is the one
   with no CSS restriction, which is exactly why its absence is worth a message. *)
VerificationTest[
    QECPauliMeasurement[five, "ZZZZZ", "ErrorCorrection" -> "Shor"],
    $Failed,
    {QECPauliMeasurement::ec},
    TestID -> "QEC-Measure-an-unbuilt-error-correction-gadget-is-refused"
]

VerificationTest[
    QECPauliMeasurement[five, "ZZZZZ"]["Nonsense"],
    Missing["NotFound", "Nonsense"],
    {QECPauliMeasurement::noprop},
    TestID -> "QEC-Measure-unknown-property-is-refused"
]
