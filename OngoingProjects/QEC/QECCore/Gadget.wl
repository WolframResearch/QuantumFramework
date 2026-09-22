(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageScope[gadgetCircuitGates]
PackageScope[gadgetCircuitOperator]
PackageScope[gadgetDiagram]
PackageScope[gadgetDepth]
PackageScope[gadgetGateCounts]
PackageScope[gadgetMeasurements]
PackageScope[gadgetHeralds]
PackageScope[gadgetLocations]
PackageScope[$gadgetCircuitProperties]


(* ============================================================================ *)
(* What every circuit-bearing object answers.                                   *)
(*                                                                              *)
(* Six heads carry an instruction list -- the syndrome circuit, the cat state,  *)
(* the Pauli measurement, the error-correction gadget, the transversal gate and *)
(* the fault-tolerant simulation -- and each had written its own "Depth",       *)
(* "InstructionCount", "GateCounts", "Measurements" and "Heralds" with the same *)
(* body.  They are written here once and dispatched to from each head, so that  *)
(* a fix to one is a fix to all six.                                            *)
(*                                                                              *)
(* THE PART THAT IS NEW is the bridge to the framework's own circuit object.    *)
(* An instruction list is flat and fast and nobody can look at it; a            *)
(* QuantumCircuitOperator draws, composes, and runs on a QuantumState.  The     *)
(* translation is one Replace, and what is interesting about it is what it must *)
(* NOT do:                                                                      *)
(*                                                                              *)
(*   - instructionEngineGates (Circuit.wl) DROPS "R" and "M".  That is correct  *)
(*     there -- it exists to validate the unitary part against the engine, and  *)
(*     neither is a gate -- and it would be a silent lie here.  A drawn circuit *)
(*     that has lost its resets and its measurements looks like a different     *)
(*     protocol.  The framework has "Reset" and "Measurement"; the engine gate  *)
(*     list does not, which is the whole reason the two translations differ.    *)
(*                                                                              *)
(*   - "Sdg" and "Vdg" go to the engine as three "S" / three "V", because the   *)
(*     engine's vocabulary has no inverse of either.  Drawn that way one        *)
(*     location becomes three boxes and the picture stops matching the          *)
(*     instruction list it came from -- and, since sec. 11.3, "which of S and   *)
(*     S-dagger is this" is a question this layer answers for a living.         *)
(*     QuantumOperator["S", {q}]["Dagger"] is one box, labelled, with the right *)
(*     matrix.                                                                  *)
(*                                                                              *)
(*   - "MH" is a measurement whose outcome rejects the shot rather than joining *)
(*     the record.  It draws as a measurement carrying that label, so a reader  *)
(*     can see at a glance which readouts are heralds.                          *)
(*                                                                              *)
(* The invariant this rests on, and the one the tests assert: ONE INSTRUCTION,  *)
(* ONE OPERATOR.  Any translation that drops, merges or expands an instruction  *)
(* breaks the correspondence between the picture and the thing the noise model, *)
(* the detector model and the Stim writer all read.                             *)
(*                                                                              *)
(* References: the framework's QuantumCircuitOperator; Got26 sec. 11.3 (why the *)
(* dagger has to be visible); the layer's own Circuit.wl for the schedule and   *)
(* the idle slots that make a location count.                                   *)
(* ============================================================================ *)


(* The properties every circuit-bearing head answers.  Each head still declares its
   own list -- the house style is explicit dispatch -- but this is the shared part. *)
$gadgetCircuitProperties = {
    "Instructions", "Qubits", "Depth", "InstructionCount", "GateCounts",
    "QuantumCircuitOperator", "Diagram"
};


(* ---- the counts, written once ---- *)

gadgetDepth[instr_List, nq_Integer] := Length[circuitSchedule[instr, nq]]

gadgetGateCounts[instr_List] := Counts[First /@ instr]

gadgetMeasurements[instr_List] := Count[instr, {"M", _}]

gadgetHeralds[instr_List] := Count[instr, {"MH", _}]

(* A location is an instruction OR a wait (Got26 Def 10.1), and the waits come from
   the schedule rather than from the list -- which is why this is not Length. *)
gadgetLocations[instr_List, nq_Integer] := Length[instr] + Length[circuitIdleSlots[instr, nq]]


(* ---- the bridge to QuantumCircuitOperator ---- *)

gadgetCircuitGates[instr_List] := Replace[instr, {
    {"R", q_Integer} :> "Reset" -> q,
    {"M", q_Integer} :> "Measurement" -> q,
    {"MH", q_Integer} :> Wolfram`QuantumFramework`QuantumMeasurementOperator[
        "Computational", {q}, "Label" -> "MH"],
    {"Sdg", q_Integer} :> Wolfram`QuantumFramework`QuantumOperator["S", {q}]["Dagger"],
    {"Vdg", q_Integer} :> Wolfram`QuantumFramework`QuantumOperator["V", {q}]["Dagger"],
    {op_String, q_Integer} :> op -> q,
    {op_String, x_Integer, y_Integer} :> op -> {x, y}
}, {1}]

gadgetCircuitOperator[instr_List] :=
    Wolfram`QuantumFramework`QuantumCircuitOperator[gadgetCircuitGates[instr]]

(* The drawing is the circuit object's own; options pass through, so a caller can
   ask for a different size or ordering without this file knowing about it. *)
gadgetDiagram[instr_List, opts___] := gadgetCircuitOperator[instr]["Diagram", opts]
