(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Idle.wlt

   The layered schedule, and idle noise on top of it.

   A qubit that no instruction of a time step touches is waiting, and waiting is a
   fault location of its own (book sec. 10.1.1, Definition 10.1).  The instruction
   list stays flat; the schedule is derived from it by ASAP list scheduling, and
   idle mechanisms hang off the layers.

   The load-bearing test in this file is the regression one: with Idle -> 0, the
   default, idleMechanisms returns {} and every number the package produced before
   the schedule existed is unchanged.  Everything else here is about the schedule
   being a valid schedule.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

QECClearCache[];

qecScope = "Wolfram`QuantumFramework`QEC`PackageScope`";

qecSchedule   = Symbol[qecScope <> "circuitSchedule"];
qecIdleQubits = Symbol[qecScope <> "circuitIdleQubits"];
qecIdleMech   = Symbol[qecScope <> "idleMechanisms"];
qecMechanisms = Symbol[qecScope <> "circuitFaultMechanisms"];
qecCodeData   = Symbol[qecScope <> "codeData"];
qecCircInstr  = Symbol[qecScope <> "codeCircuitInstructions"];
qecFrame      = Symbol[qecScope <> "framePropagate"];
qecFaultEffect = Symbol[qecScope <> "faultEffect"];
qecStabCount  = Symbol[qecScope <> "codeStabilizerCount"];
qecInstrLayers = Symbol[qecScope <> "circuitInstructionLayers"];
qecAnchor     = Symbol[qecScope <> "circuitIdleAnchor"];
qecSlots      = Symbol[qecScope <> "circuitIdleSlots"];

(* the circuit under test, and its register size *)
qecInstr[code_, r_ : 1] := qecCircInstr[qecCodeData[code], r];
qecNq[code_] := code["Qubits"] + qecStabCount[qecCodeData[code]];
qecLayers[code_, r_ : 1] := qecSchedule[qecInstr[code, r], qecNq[code]];

(* which layer each instruction landed in *)
qecLayerOf[layers_, i_] := First[FirstPosition[layers, i]];


(* ============================================================================
   The schedule is a schedule
   ============================================================================ *)

(* Every instruction appears exactly once. *)
VerificationTest[
    With[{c = QECCode["5QubitCode"]},
        Sort[Catenate[qecLayers[c]]] === Range[Length[qecInstr[c]]]
    ],
    True,
    TestID -> "QEC-Idle-schedule-is-a-partition"
]

(* No qubit is used twice inside one layer: that is what makes a layer a parallel
   time step rather than just a bucket. *)
VerificationTest[
    With[{c = QECCode["5QubitCode"], instr = qecInstr[QECCode["5QubitCode"]]},
        AllTrue[
            qecLayers[c],
            With[{qs = Catenate[Rest[instr[[#]]] & /@ #]}, DuplicateFreeQ[qs]] &
        ]
    ],
    True,
    TestID -> "QEC-Idle-no-qubit-twice-in-a-layer"
]

(* Per-qubit order is preserved: for each qubit, the layers of the instructions
   touching it are strictly increasing.  This is the correctness condition -- a
   schedule that reordered two gates on one qubit would be a different circuit. *)
VerificationTest[
    With[
        {c = QECCode["5QubitCode"]},
        With[
            {instr = qecInstr[c], layers = qecLayers[c], nq = qecNq[c]},
            AllTrue[
                Range[nq],
                Function[q,
                    OrderedQ[
                        qecLayerOf[layers, #] & /@
                            Select[Range[Length[instr]], MemberQ[Rest[instr[[#]]], q] &],
                        Less
                    ]
                ]
            ]
        ]
    ],
    True,
    TestID -> "QEC-Idle-preserves-per-qubit-order"
]

(* Layering shortens the circuit.  It has to: two instructions on disjoint qubits
   share a time step, and the bit-flip circuit has several such pairs. *)
VerificationTest[
    With[{c = QECCode["BitFlipCode"]},
        Length[qecLayers[c]] < Length[qecInstr[c]]
    ],
    True,
    TestID -> "QEC-Idle-layering-shortens-the-circuit"
]

(* The concrete case, pinned: the bit-flip circuit is 8 instructions in 6 layers.
   {R,4} and {R,5} share the first; {M,4} and {CNOT,2,5} share the fourth. *)
VerificationTest[
    With[{c = QECCode["BitFlipCode"]},
        {Length[qecInstr[c]], Length[qecLayers[c]]}
    ],
    {8, 6},
    TestID -> "QEC-Idle-bitflip-schedule-is-six-layers"
]

(* An empty instruction list has an empty schedule, not an error: rounds -> 0 is a
   legitimate circuit (the code-capacity limit -- an error and nothing to measure
   it with). *)
VerificationTest[
    qecSchedule[{}, 5],
    {},
    TestID -> "QEC-Idle-empty-circuit-empty-schedule"
]

(* Busy and idle partition the register in every layer. *)
VerificationTest[
    With[
        {c = QECCode["BitFlipCode"]},
        With[
            {instr = qecInstr[c], nq = qecNq[c]},
            AllTrue[
                qecLayers[c],
                With[{busy = Union @@ (Rest[instr[[#]]] & /@ #)},
                    Sort[Join[busy, qecIdleQubits[instr, nq, #]]] === Range[nq]
                ] &
            ]
        ]
    ],
    True,
    TestID -> "QEC-Idle-busy-and-idle-partition-the-register"
]

(* A qubit the circuit never touches is idle in every layer.  This is why the
   register size is a parameter and not derived from the instructions: spectator
   qubits decohere, and nothing in the instruction list knows they exist. *)
VerificationTest[
    With[
        {c = QECCode["BitFlipCode"]},
        With[
            {instr = qecInstr[c], nq = qecNq[c] + 2},
            AllTrue[qecLayers[c], SubsetQ[qecIdleQubits[instr, nq, #], {nq - 1, nq}] &]
        ]
    ],
    True,
    TestID -> "QEC-Idle-spectator-qubits-idle-everywhere"
]


(* ============================================================================
   Idle mechanisms
   ============================================================================ *)

(* THE REGRESSION TEST.  Idle -> 0 is the default, and there it emits nothing at
   all -- not a list of zero-probability rows.  So the mechanism list, the detector
   model, the decoder table and every rate downstream are what they were before the
   schedule existed, and the rest of the suite is unaffected by this file. *)
VerificationTest[
    With[{c = QECCode["BitFlipCode"]},
        qecIdleMech[qecInstr[c], qecNq[c], <|"Idle" -> 0|>]
    ],
    {},
    TestID -> "QEC-Idle-zero-rate-emits-nothing"
]

VerificationTest[
    With[{c = QECCode["BitFlipCode"]},
        Length[qecMechanisms[qecCodeData[c], First[QECNoiseModel["Circuit", p]], 1]]
    ],
    64,
    TestID -> "QEC-Idle-default-mechanism-count-unchanged"
]

(* With a rate, three Paulis per idle qubit per layer, and that count is exactly
   what the schedule says it should be. *)
VerificationTest[
    With[
        {c = QECCode["BitFlipCode"]},
        With[
            {instr = qecInstr[c], nq = qecNq[c]},
            Length[qecIdleMech[instr, nq, <|"Idle" -> q|>]] ===
                3 Total[Length[qecIdleQubits[instr, nq, #]] & /@ qecLayers[c]]
        ]
    ],
    True,
    TestID -> "QEC-Idle-mechanism-count-matches-the-schedule"
]

(* Each is one location per (layer, qubit), with the three Paulis as its mutually
   exclusive outcomes -- so demGroups sees groups of three, not singletons. *)
VerificationTest[
    With[
        {c = QECCode["BitFlipCode"]},
        Union[Values[Counts[Lookup[qecIdleMech[qecInstr[c], qecNq[c], <|"Idle" -> q|>], "Location"]]]]
    ],
    {3},
    TestID -> "QEC-Idle-three-paulis-per-location"
]

(* The concrete count, pinned: 6 layers with 3+3+3+2+3+4 = 18 idle qubit-steps, three
   Paulis each.  Serial extraction would have 8 steps and 28 idle qubit-steps, so
   scheduling is what makes the cost payable rather than merely visible. *)
VerificationTest[
    With[
        {c = QECCode["BitFlipCode"]},
        With[
            {instr = qecInstr[c], nq = qecNq[c]},
            {
                Length[qecIdleQubits[instr, nq, #]] & /@ qecLayers[c],
                Length[qecIdleMech[instr, nq, <|"Idle" -> q|>]]
            }
        ]
    ],
    {{3, 3, 3, 2, 3, 4}, 54},
    TestID -> "QEC-Idle-bitflip-idle-count"
]

(* The anchor is per qubit, not per layer.  Data qubit 1 of the bit-flip code is
   touched exactly once, by the CNOT at instruction 2 in layer 2, so it is idle in
   layers 1, 3, 4, 5 and 6 -- anchored before the circuit for the first, and after
   that one instruction for the other four.  This is the test that would have caught
   the wrong anchor: "the end of layer 1" is instruction 5, which in emission order
   comes after instructions 2, 3 and 4. *)
VerificationTest[
    With[
        {c = QECCode["BitFlipCode"]},
        With[
            {instr = qecInstr[c]},
            With[
                {layerOf = qecInstrLayers[instr, qecNq[c]]},
                qecAnchor[instr, layerOf, 1, #] & /@ {1, 3, 4, 5, 6}
            ]
        ]
    ],
    {0, 2, 2, 2, 2},
    TestID -> "QEC-Idle-anchor-is-per-qubit"
]

(* Every anchor is either 0 or an instruction that really touches that qubit, in a
   really earlier layer. *)
VerificationTest[
    With[
        {c = QECCode["5QubitCode"]},
        With[
            {instr = qecInstr[c], nq = qecNq[c]},
            With[
                {layers = qecLayers[c], layerOf = qecInstrLayers[instr, qecNq[c]]},
                AllTrue[
                    Catenate @ Table[
                        {L, q, qecAnchor[instr, layerOf, q, L]},
                        {L, Length[layers]}, {q, qecIdleQubits[instr, nq, layers[[L]]]}
                    ],
                    Function[t,
                        t[[3]] === 0 ||
                            (MemberQ[Rest[instr[[t[[3]]]]], t[[2]]] && layerOf[[t[[3]]]] < t[[1]])
                    ]
                ]
            ]
        ]
    ],
    True,
    TestID -> "QEC-Idle-anchor-is-well-formed"
]

(* A symbolic rate is never structurally zero, so it is always emitted: nothing can
   be concluded about it, and dropping it would silently answer a question the user
   asked symbolically. *)
VerificationTest[
    With[{c = QECCode["BitFlipCode"]},
        Length[qecIdleMech[qecInstr[c], qecNq[c], <|"Idle" -> pIdle|>]] > 0
    ],
    True,
    TestID -> "QEC-Idle-symbolic-rate-is-emitted"
]

(* Turning the rate on grows the mechanism list, and by exactly the idle count. *)
VerificationTest[
    With[
        {c = QECCode["BitFlipCode"]},
        With[
            {without = Length[qecMechanisms[qecCodeData[c], First[QECNoiseModel["Circuit", p]], 1]],
             with = Length[qecMechanisms[qecCodeData[c],
                 First[QECNoiseModel["Circuit", <|"OneQubit" -> p, "TwoQubit" -> p,
                     "Measurement" -> p, "Reset" -> p, "Idle" -> p|>]], 1]]},
            with - without ===
                Length[qecIdleMech[qecInstr[c], qecNq[c], <|"Idle" -> p|>]]
        ]
    ],
    True,
    TestID -> "QEC-Idle-wiring-adds-exactly-the-idle-mechanisms"
]


(* ============================================================================
   What an idle fault does
   ============================================================================ *)

(* Waiting *before* the check is caught: the anchor is 0, so the X is there when the
   ladder reads the qubit, and the check fires exactly as a code-capacity error would. *)
VerificationTest[
    With[{c = QECCode["BitFlipCode"]},
        qecFrame[qecInstr[c], 5, {{0, 1, {1, 0}}}]["Record"]
    ],
    {1, 0},
    TestID -> "QEC-Idle-waiting-before-the-check-is-caught"
]

(* Waiting *after* it is not caught this round -- nothing touches that qubit again --
   and yet it is not lost: the noiseless final readout sees the residual, so the
   closing detector fires.  Without that detector this fault would be invisible. *)
VerificationTest[
    With[{c = QECCode["BitFlipCode"]},
        With[
            {instr = qecInstr[c]},
            {
                qecFrame[instr, 5, {{2, 1, {1, 0}}}]["Record"],
                AnyTrue[First[qecFaultEffect[qecCodeData[c], instr, 5, 2, 1, {{2, 1, {1, 0}}}]], # === 1 &]
            }
        ]
    ],
    {{0, 0}, True},
    TestID -> "QEC-Idle-waiting-after-the-check-shows-in-the-final-readout"
]

(* And a phase idle on a data qubit of the bit-flip code is undetectable and logical
   -- the same blind spot a Z has anywhere else in this code.  Quantum distance one,
   not the three of the classical repetition code. *)
VerificationTest[
    With[{c = QECCode["BitFlipCode"]},
        With[
            {eff = qecFaultEffect[qecCodeData[QECCode["BitFlipCode"]],
                qecInstr[QECCode["BitFlipCode"]], 5, 2, 1, {{0, 1, {0, 1}}}]},
            {AllTrue[eff[[1]], # === 0 &], AnyTrue[eff[[2]], # === 1 &]}
        ]
    ],
    {True, True},
    TestID -> "QEC-Idle-phase-idle-is-undetectable-and-logical"
]


(* ============================================================================
   Rounds pipeline
   ============================================================================ *)

(* r rounds are NOT r copies of one round's schedule: an ancilla's reset for round 2
   schedules while the other ancilla of round 1 is still being measured.  This is why
   the anchors must be computed over the whole circuit and not per round, and it is
   the reason the Stim writer indexes into a whole-circuit anchor table by offset. *)
VerificationTest[
    With[
        {c = QECCode["BitFlipCode"]},
        With[
            {one = Length[qecLayers[c, 1]], two = Length[qecLayers[c, 2]]},
            {one, two, two < 2 one}
        ]
    ],
    {6, 10, True},
    TestID -> "QEC-Idle-rounds-pipeline"
]


(* ============================================================================
   The Stim bridge agrees with the model
   ============================================================================ *)

(* THE CROSS-CHECK.  Both sides read circuitIdleSlots, so the exported circuit has
   one DEPOLARIZE1 per slot and the detector model has three mechanisms per slot.
   If these ever disagree the bridge is no longer checking anything. *)
VerificationTest[
    With[
        {c = QECCode["BitFlipCode"],
         noise = QECNoiseModel["Circuit", <|"OneQubit" -> 1/1000, "TwoQubit" -> 1/1000,
             "Measurement" -> 1/1000, "Reset" -> 1/1000, "Idle" -> 1/500|>]},
        With[
            {slots = Length[qecSlots[qecInstr[c, 2], qecNq[c]]],
             dem = Length[Select[qecMechanisms[qecCodeData[c], First[noise], 2],
                 First[#["Location"]] === "Idle" &]],
             stim = StringCount[QECStimCircuit[c, noise, 2], "DEPOLARIZE1(0.002)"]},
            {dem === 3 slots, stim === slots}
        ]
    ],
    {True, True},
    TestID -> "QEC-Idle-stim-matches-the-model"
]

(* Regression: the bit-flip extraction circuit has no one-qubit gates at all (both
   checks are pure Z, so nothing needs rotating), so with Idle -> 0 -- the default --
   the exported circuit carries no DEPOLARIZE1 whatsoever.  Any appearing there would
   mean idle noise leaked into a circuit that did not ask for it. *)
VerificationTest[
    StringCount[
        QECStimCircuit[QECCode["BitFlipCode"], QECNoiseModel["Circuit", 1/1000], 2],
        "DEPOLARIZE1"],
    0,
    TestID -> "QEC-Idle-stim-clean-when-rate-is-zero"
]

(* The final round is the perfect readout and must stay perfect: no idle there, or the
   residual it is supposed to reveal would be corrupted by the act of revealing it. *)
VerificationTest[
    With[
        {src = QECStimCircuit[QECCode["BitFlipCode"],
            QECNoiseModel["Circuit", <|"OneQubit" -> 1/1000, "TwoQubit" -> 1/1000,
                "Measurement" -> 1/1000, "Reset" -> 1/1000, "Idle" -> 1/500|>], 2]},
        StringFreeQ[Last[StringSplit[src, "# --- final noiseless round ---"]], "DEPOLARIZE"]
    ],
    True,
    TestID -> "QEC-Idle-final-round-stays-noiseless"
]
