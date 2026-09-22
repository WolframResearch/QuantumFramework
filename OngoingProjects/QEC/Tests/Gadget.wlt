(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Gadget.wlt

   The trait every circuit-bearing object answers, and the bridge from the flat
   instruction list to the framework's own QuantumCircuitOperator.

   The load-bearing test is QEC-Gadget-one-instruction-one-operator, asserted for
   all six heads at once.  A translation that drops an instruction, merges two, or
   expands one into three still draws a plausible picture -- and the picture then
   describes a different protocol from the one the noise model, the detector model
   and the Stim writer all read off the same list.  Two of those three failure
   modes were live before this trait existed: instructionEngineGates drops "R" and
   "M" (correctly, for its own purpose), and it expands "Sdg" and "Vdg" into three
   gates because the engine has no inverse for either.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

QECClearCache[];

qecScope = "Wolfram`QuantumFramework`QEC`PackageScope`";

qecGadgetGates = Symbol[qecScope <> "gadgetCircuitGates"];
qecGadgetOp    = Symbol[qecScope <> "gadgetCircuitOperator"];
qecGadgetDepth = Symbol[qecScope <> "gadgetDepth"];
qecGadgetLocs  = Symbol[qecScope <> "gadgetLocations"];
qecFtLocations = Symbol[qecScope <> "ftLocations"];

steane = QECCode["SteaneCode"];
five = QECCode["5QubitCode"];

(* One of every circuit-bearing head. *)
gadgets = {
    QECSyndromeCircuit[steane],
    QECCatState[4],
    QECPauliMeasurement[five, "XZZXI"],
    QECErrorCorrection[steane],
    QECTransversalGate[steane, "S"],
    QECFaultTolerant[{{"R", 1}, {"H", 1}, {"M", 1}}, steane]
};


(* ============================================================================
   The bridge
   ============================================================================ *)

(* THE test.  One instruction, one operator, for every head.  Everything else in
   this file is a detail of that. *)
VerificationTest[
    Table[Length[g["Instructions"]] === Length[g["QuantumCircuitOperator"]["Operators"]], {g, gadgets}],
    ConstantArray[True, 6],
    TestID -> "QEC-Gadget-one-instruction-one-operator"
]

VerificationTest[
    Table[Head[g["QuantumCircuitOperator"]], {g, gadgets}],
    ConstantArray[QuantumCircuitOperator, 6],
    TestID -> "QEC-Gadget-every-head-gives-a-circuit-operator"
]

VerificationTest[
    Table[Head[g["Diagram"]], {g, gadgets}],
    ConstantArray[Graphics, 6],
    TestID -> "QEC-Gadget-every-head-draws"
]

(* `R` and `M` are not gates, so instructionEngineGates drops them; they are
   locations, so the drawing must not.  The cat state has both, plus a herald. *)
VerificationTest[
    qecGadgetGates[{{"R", 1}, {"H", 1}, {"CNOT", 1, 2}, {"M", 2}, {"MH", 3}}][[All, 0]],
    {Rule, Rule, Rule, Rule, QuantumMeasurementOperator},
    TestID -> "QEC-Gadget-resets-and-measurements-survive-the-bridge"
]

VerificationTest[
    #["Label"] & /@ QECCatState[4]["QuantumCircuitOperator"]["Operators"],
    {"Reset"["0"], "Reset"["0"], "Reset"["0"], "Reset"["0"], "H",
     Subscript["C", "NOT"][{1}, {}], Subscript["C", "NOT"][{2}, {}],
     Subscript["C", "NOT"][{3}, {}], "Reset"["0"],
     Subscript["C", "NOT"][{2}, {}], Subscript["C", "NOT"][{3}, {}], "MH"},
    TestID -> "QEC-Gadget-the-cat-draws-its-resets-and-its-herald"
]

(* A herald is a measurement whose outcome rejects the shot, so it draws as a
   measurement that says so, rather than as an anonymous one. *)
VerificationTest[
    Last[#["Label"] & /@ QECCatState[4]["QuantumCircuitOperator"]["Operators"]],
    "MH",
    TestID -> "QEC-Gadget-a-herald-is-labelled-as-one"
]

(* Sdg is one location and draws as one box.  Sending it to the engine as three S
   is right for the engine and wrong for a picture -- and since sec. 11.3, telling
   S from S-dagger is a thing this layer does for a living. *)
VerificationTest[
    With[{ops = QECTransversalGate[steane, "Sdg"]["QuantumCircuitOperator"]["Operators"]},
        {Length[ops], DeleteDuplicates[#["Label"] & /@ ops],
         Normal[First[ops]["MatrixRepresentation"]]}
    ],
    {7, {SuperDagger["S"]}, {{1, 0}, {0, -I}}},
    TestID -> "QEC-Gadget-Sdg-is-one-box-with-the-right-matrix"
]


(* ============================================================================
   The shared counts
   ============================================================================ *)

(* The trait's bodies and the heads' properties are the same computation; this is
   the regression that keeps them that way after the six copies were removed. *)
VerificationTest[
    Table[g["Depth"] === qecGadgetDepth[g["Instructions"], g["Qubits"]],
        {g, {QECSyndromeCircuit[steane], QECPauliMeasurement[five, "XZZXI"],
             QECErrorCorrection[steane]}}],
    {True, True, True},
    TestID -> "QEC-Gadget-the-shared-depth-is-the-heads-depth"
]

(* ftLocations was the only other place that knew a location is an instruction or a
   wait; it now delegates, and this says so. *)
VerificationTest[
    With[{c = {{"R", 1}, {"R", 2}, {"H", 1}, {"CNOT", 1, 2}, {"M", 1}, {"M", 2}}},
        {qecFtLocations[c, 2], qecGadgetLocs[c, 2]}
    ],
    {7, 7},
    TestID -> "QEC-Gadget-locations-are-counted-in-one-place"
]

(* The two new properties are listed, so code[\"Properties\"] is still the answer to
   \"what can I ask this object\". *)
VerificationTest[
    Table[SubsetQ[g["Properties"], {"QuantumCircuitOperator", "Diagram"}], {g, gadgets}],
    ConstantArray[True, 6],
    TestID -> "QEC-Gadget-the-new-properties-are-listed"
]
