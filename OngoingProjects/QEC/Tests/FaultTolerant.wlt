(* ::Package:: *)

(* ============================================================================
   Tests/QEC/FaultTolerant.wlt

   FT(C): the ideal fault-tolerant simulation of a circuit, Got26 Def 10.6.

   The assembler itself is bookkeeping -- take each location, put its gadget in
   its place, put an error correction gadget after every preparation, gate and
   storage gadget, and none after a measurement -- so most of these tests count
   things.  Two do not.

   QEC-FT-the-gadget-for-a-logical-S-is-the-transversal-Sdg is the one the whole
   file turns on, and QEC-FT-emitting-the-gate-by-name-computes-the-conjugate is
   its physical half: assembled properly the protocol leaves the encoded qubit in
   the +1 eigenstate of Ybar, and assembled by NAME -- transversal S for a logical
   S -- it leaves it in the -1 eigenstate.  Same circuit, conjugate answer.  That
   is sec. 11.3 reaching all the way out to a state on the engine.

   The other one is QEC-FT-makes-the-logical-Bell-state: preparation and gate
   gadgets only, run with no faults, and the two encoded blocks come out
   stabilised by XbarXbar and ZbarZbar.  That is Def 10.4's requirement -- decode
   after the gadget is the same as doing the location -- checked end to end rather
   than gadget by gadget.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

QECClearCache[];

qecScope = "Wolfram`QuantumFramework`QEC`PackageScope`";

qecGateWord     = Symbol[qecScope <> "ftGateWord"];
qecReadout      = Symbol[qecScope <> "ftReadoutSupport"];
qecLocations    = Symbol[qecScope <> "ftLocations"];
qecCodeData     = Symbol[qecScope <> "codeData"];
qecFrame        = Symbol[qecScope <> "framePropagate"];
qecEngineGates  = Symbol[qecScope <> "instructionEngineGates"];
qecApplyGates   = Symbol[qecScope <> "applyGates"];
qecOneQubitOps  = Symbol[qecScope <> "$circuitOneQubitOps"];
qecTwoQubitOps  = Symbol[qecScope <> "$circuitTwoQubitOps"];

steane = QECCode["SteaneCode"];
five = QECCode["5QubitCode"];
four = QECCode["DistanceTwo", 4];

steaneData = qecCodeData[steane];

(* R, H, CNOT, M -- one of each kind of location, and one block left idling. *)
ft = QECFaultTolerant[{{"R", 1}, {"R", 2}, {"H", 1}, {"CNOT", 1, 2}, {"M", 1}, {"M", 2}}, steane];

(* The instructions of the gadgets of a given kind, in order. *)
qecGadgetInstructions[obj_, kinds_List] := obj["GadgetInstructions", kinds]

(* The state the circuit leaves, with no faults: the engine's opinion. *)
qecFinalState[instr_, nq_] :=
    qecApplyGates[PauliStabilizer[nq], qecEngineGates[instr]]

reg = QECRegister[steane, 2];
xbar = QECPauliString /@ reg["LogicalVectors"]["X"];
zbar = QECPauliString /@ reg["LogicalVectors"]["Z"];

(* Ybar = i Xbar Zbar, Hermitian, on the first block of a 14-qubit register. *)
ybar = QECPauliString @ MapAt[Mod[# + 1, 4] &, QECPauliProduct[xbar[[1]], zbar[[1]]], -1];


(* ============================================================================
   The substitution (Def 10.6)
   ============================================================================ *)

(* One gadget per location, and the error correction gadgets: after both
   preparations, after the gate, after the storage, and never after a measurement.
   Six of them for the four layers this circuit schedules into. *)
VerificationTest[
    ft["GadgetCounts"],
    <|"Preparation" -> 2, "ErrorCorrection" -> 6, "Gate" -> 2, "Storage" -> 1,
      "Measurement" -> 2|>,
    TestID -> "QEC-FT-every-location-gets-its-gadget"
]

(* The order matters as much as the count: no error correction gadget may follow a
   measurement gadget, because its output is classical. *)
VerificationTest[
    With[{kinds = #["Gadget"] & /@ Normal[ft["Gadgets"]]},
        {Last[kinds], MemberQ[Partition[kinds, 2, 1], {"Measurement", "ErrorCorrection"}]}
    ],
    {"Measurement", False},
    TestID -> "QEC-FT-no-error-correction-after-a-measurement"
]

(* A storage gadget is a wait on every qubit of the block, so it emits nothing --
   and it is still a gadget, and still earns an error correction gadget.  Charging
   for it is what keeps the overheads honest. *)
VerificationTest[
    SelectFirst[Normal[ft["Gadgets"]], #["Gadget"] === "Storage" &],
    <|"Gadget" -> "Storage", "Location" -> {"Wait", 2}, "Blocks" -> {2},
      "Instructions" -> 0|>,
    TestID -> "QEC-FT-a-wait-is-a-storage-gadget-that-emits-nothing"
]

(* The smallest FT(C) there is: one preparation, one error correction. *)
VerificationTest[
    With[{one = QECFaultTolerant[{{"R", 1}}, steane]},
        {one["GadgetCounts"], one["Qubits"], one["CircuitLocations"]}
    ],
    {<|"Preparation" -> 1, "ErrorCorrection" -> 1|>, 14, 1},
    TestID -> "QEC-FT-one-location-is-a-gadget-and-a-correction"
]

(* Nothing exotic comes out: the assembled circuit speaks the same instruction
   vocabulary every other circuit in the layer does. *)
VerificationTest[
    Complement[
        Union[First /@ ft["Instructions"]],
        Join[qecOneQubitOps, qecTwoQubitOps, {"R", "M"}]
    ],
    {},
    TestID -> "QEC-FT-the-assembled-circuit-is-an-ordinary-circuit"
]


(* ============================================================================
   Which gadget, and why it is not the one with the same name
   ============================================================================ *)

(* THE test.  The gadget for a logical S is the transversal S^dagger (sec. 11.3),
   and for a logical S^dagger it is the transversal S.  The assembler searches by
   logical action, so it finds this rather than trusting the name. *)
VerificationTest[
    {qecGateWord[steaneData, "S"], qecGateWord[steaneData, "Sdg"],
     qecGateWord[steaneData, "H"], qecGateWord[steaneData, "T"]},
    {"Sdg", "S", "H", Missing["NoGadget", "T"]},
    TestID -> "QEC-FT-the-gadget-for-a-logical-S-is-the-transversal-Sdg"
]

(* And the same claim as a state.  R, H, S on one logical qubit should leave the
   encoded qubit in the +1 eigenstate of Ybar.  Assembled by NAME -- one
   transversal S for a logical S -- it lands on -1: the conjugate circuit, which
   is exactly what sec. 11.3 warns about and what a sign-blind layer would ship. *)
VerificationTest[
    Module[{gates, named},
        gates = qecGadgetInstructions[
            QECFaultTolerant[{{"R", 1}, {"H", 1}, {"S", 1}}, steane],
            {"Preparation", "Gate"}];
        named = Join[Drop[gates, -7], Table[{"S", q}, {q, 7}]];
        {qecFinalState[gates, 14]["Expectation", ybar],
         qecFinalState[named, 14]["Expectation", ybar]}
    ],
    {1, -1},
    TestID -> "QEC-FT-emitting-the-gate-by-name-computes-the-conjugate"
]

(* Def 10.4 end to end: the gadgets for R, R, H, CNOT leave the two blocks in the
   encoded Bell state, stabilised by XbarXbar and ZbarZbar. *)
VerificationTest[
    Module[{state, prod},
        state = qecFinalState[
            qecGadgetInstructions[
                QECFaultTolerant[{{"R", 1}, {"R", 2}, {"H", 1}, {"CNOT", 1, 2}}, steane],
                {"Preparation", "Gate"}],
            28];
        prod[s1_, s2_] := QECPauliString[
            Join[QECPauliVector[QECPauliProduct[s1, s2]][[1 ;; 14]], ConstantArray[0, 14],
                 QECPauliVector[QECPauliProduct[s1, s2]][[15 ;; 28]], ConstantArray[0, 14], {0}]];
        {state["Expectation", prod[xbar[[1]], xbar[[2]]]],
         state["Expectation", prod[zbar[[1]], zbar[[2]]]]}
    ],
    {1, 1},
    TestID -> "QEC-FT-makes-the-logical-Bell-state"
]

(* With no faults nothing fires: every measurement in the assembled circuit comes
   back 0 and the frame ends where it started.  Cheap, and it is what would break
   first if a gadget were relocated onto the wrong qubits. *)
VerificationTest[
    With[{run = qecFrame[ft["Instructions"], ft["Qubits"], {}]},
        {Total[run["Record"]], Total /@ run["Frame"]}
    ],
    {0, {0, 0}},
    TestID -> "QEC-FT-with-no-faults-no-syndrome-fires"
]


(* The gadgets can be taken apart again: each kind's instructions come back on
   request, which is how anything downstream gets at "the gate gadgets only". *)
VerificationTest[
    {Length[ft["GadgetInstructions", "Gate"]],
     Union[First /@ ft["GadgetInstructions", "Gate"]],
     Union[First /@ QECFaultTolerant[{{"R", 1}, {"S", 1}}, steane]["GadgetInstructions", "Gate"]]},
    {14, {"CNOT", "H"}, {"Sdg"}},
    TestID -> "QEC-FT-the-gadgets-can-be-taken-apart-again"
]


(* ============================================================================
   The classical side of the measurement gadget
   ============================================================================ *)

(* The logical bit is the parity of the outcomes over the support of Zbar, and the
   positions are absolute in the record -- the measurement gadget's n outcomes sit
   among the error correction gadget's, so a decoder needs to be told where. *)
VerificationTest[
    {qecReadout[steaneData], ft["Readouts"], ft["Measurements"]},
    {{2, 4, 6}, <|1 -> {86, 88, 90}, 2 -> {93, 95, 97}|>, 98},
    TestID -> "QEC-FT-the-logical-bit-is-a-parity-of-named-outcomes"
]

(* Those positions really are the measurement gadget's own, not an error
   correction gadget's: they point at data qubits of block 1 and block 2. *)
VerificationTest[
    With[{ms = Cases[ft["Instructions"], {"M", q_} :> q]},
        {ms[[ft["Readouts"][1]]], ms[[ft["Readouts"][2]]]}
    ],
    {{2, 4, 6}, {9, 11, 13}},
    TestID -> "QEC-FT-the-readout-positions-point-at-the-data-block"
]


(* ============================================================================
   The three overheads (Def 10.6)
   ============================================================================ *)

(* Locations, both sides, counted the same way: instructions plus waits.  The
   circuit has six instructions and one wait -- block 2 idles while block 1 gets
   its Hadamard. *)
VerificationTest[
    {ft["CircuitLocations"], qecLocations[ft["Circuit"], 2]},
    {7, 7},
    TestID -> "QEC-FT-a-wait-counts-as-a-location-on-both-sides"
]

(* The numbers themselves, pinned: a change anywhere in the gadgets moves them, and
   that is worth noticing rather than absorbing. *)
VerificationTest[
    {ft["Qubits"], ft["Locations"], ft["Depth"], ft["Overheads"]},
    {28, 1841, 77, <|"Size" -> 263, "Qubits" -> 14, "Depth" -> 77/4|>},
    TestID -> "QEC-FT-the-three-overheads-of-Definition-10-6"
]

(* The qubit overhead is the block plus its error correction workspace, per logical
   qubit: 7 data and 7 ancilla, and the ancilla is per block so that the
   corrections of one layer can run in parallel. *)
VerificationTest[
    With[{three = QECFaultTolerant[{{"R", 1}, {"R", 2}, {"R", 3}}, steane]},
        {three["Qubits"], three["QubitOverhead"]}
    ],
    {42, 14},
    TestID -> "QEC-FT-every-block-brings-its-own-correction-workspace"
]


(* ============================================================================
   What the protocol does not have
   ============================================================================ *)

(* Three assumptions, carried rather than hidden: the encoder is not fault
   tolerant, the gate set is Clifford, and corrections live in the frame. *)
VerificationTest[
    Length[ft["OpenAssumptions"]],
    3,
    TestID -> "QEC-FT-reports-what-it-is-assuming"
]

(* Outside the Clifford group there is no transversal gadget, and chapter 13 is
   where that is fixed -- so it is refused rather than approximated. *)
VerificationTest[
    QECFaultTolerant[{{"R", 1}, {"T", 1}}, steane],
    $Failed,
    {QECFaultTolerant::gate},
    TestID -> "QEC-FT-refuses-a-gate-with-no-gadget"
]

(* Steane EC is CSS-only, so the protocol is. *)
VerificationTest[
    QECFaultTolerant[{{"R", 1}}, five],
    $Failed,
    {QECFaultTolerant::css},
    TestID -> "QEC-FT-needs-a-CSS-code"
]

(* One logical qubit per block, which is Def 10.4's own simplification and what
   makes "the gadget for logical g" a well-posed question. *)
VerificationTest[
    QECFaultTolerant[{{"R", 1}}, four],
    $Failed,
    {QECFaultTolerant::logical},
    TestID -> "QEC-FT-needs-one-logical-qubit-per-block"
]

VerificationTest[
    {QECFaultTolerant[{{"H", 1}}, steane],
     QECFaultTolerant[{{"R", 1}, {"M", 1}, {"H", 1}}, steane]},
    {$Failed, $Failed},
    {QECFaultTolerant::unprepared, QECFaultTolerant::measured},
    TestID -> "QEC-FT-a-block-has-to-be-live-to-be-used"
]

VerificationTest[
    ft["Nonsense"],
    Missing["NotFound", "Nonsense"],
    {QECFaultTolerant::noprop},
    TestID -> "QEC-FT-unknown-property-is-refused"
]
