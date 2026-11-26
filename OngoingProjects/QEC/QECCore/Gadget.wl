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
PackageScope[instructionRules]
PackageScope[instructionTuples]
PackageScope[instructionRule]
PackageScope[instructionTuple]
PackageScope[internalOp]
PackageScope[publicOp]
PackageScope[publicInstructionQ]
PackageScope[oldInstructionQ]


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


(* ---- the two spellings of an instruction ---- *)

(* Inside the layer an instruction is a flat tuple, {"CNOT", c, t}, which is what the
   frame propagator, the scheduler and every gadget builder pattern-match on.  At the
   public surface it is spelled exactly as QuantumCircuitOperator spells it, so that a
   list the layer hands out is a circuit the framework accepts as it is:

     {"R", q}        "Reset" -> q
     {"M", q}        "Measurement" -> q
     {"MH", q}       a computational QuantumMeasurementOperator on q labelled "Herald"
                     (a verification readout, kept apart from the record)
     {"Sdg", q}      SuperDagger["S"] -> q      (likewise "Vdg", "Tdg")
     {"CNOT", c, t}  "CNOT" -> {c, t}           (every other gate keeps its name)

   The two functions below are the only crossing between them, and both accept either
   spelling, so a caller never has to know which one it holds. *)
$publicOpNames = <|"R" -> "Reset", "M" -> "Measurement", "MH" -> "Herald"|>;
$internalOpNames = <|"Reset" -> "R", "Measurement" -> "M", "Herald" -> "MH"|>;
$daggerOps = <|"Sdg" -> "S", "Vdg" -> "V", "Tdg" -> "T"|>;
$undaggerOps = <|"S" -> "Sdg", "V" -> "Vdg", "T" -> "Tdg"|>;

publicOp[op_String] := Which[
    KeyExistsQ[$publicOpNames, op], $publicOpNames[op],
    KeyExistsQ[$daggerOps, op], SuperDagger[$daggerOps[op]],
    True, op
]

(* built once per qubit: a measurement operator is not free to construct *)
heraldOperator[q_Integer] := heraldOperator[q] =
    Wolfram`QuantumFramework`QuantumMeasurementOperator["Computational", {q}, "Label" -> "Herald"]

instructionRule[{"MH", q_Integer}] := heraldOperator[q]
instructionRule[{op_String, q_}] := publicOp[op] -> q
instructionRule[{op_String, qs__}] := publicOp[op] -> {qs}
instructionRule[r : ((_String | _SuperDagger) -> _)] := r
instructionRule[m_Wolfram`QuantumFramework`QuantumMeasurementOperator] := m

instructionRules[instr_List] := instructionRule /@ instr

internalOp[op_String] := Lookup[$internalOpNames, op, op]
internalOp[SuperDagger[g_String]] := Lookup[$undaggerOps, g, g <> "dg"]

instructionTuple[(op : (_String | _SuperDagger)) -> qs_List] := Prepend[qs, internalOp[op]]
instructionTuple[(op : (_String | _SuperDagger)) -> q_] := {internalOp[op], q}
instructionTuple[{q_Integer}] := {"M", q}
instructionTuple[m_Wolfram`QuantumFramework`QuantumMeasurementOperator] :=
    {If[m["Label"] === "Herald", "MH", "M"], First[m["InputOrder"]]}
instructionTuple[t : {_String, __}] := t

instructionTuples[instr_List] := instructionTuple /@ instr


(* What the public surface accepts, and the spellings it no longer does: the tuple,
   and the layer's own short names in a rule. *)
publicInstructionQ[x_] := MatchQ[x,
    ((_String | _SuperDagger) -> (_Integer | {__Integer})) | {_Integer} | _Wolfram`QuantumFramework`QuantumMeasurementOperator]
oldInstructionQ[x_] := MatchQ[x, {_String, __Integer} | (("R" | "M" | "MH" | "Sdg" | "Vdg" | "Tdg") -> _)]


(* ---- the counts, written once ---- *)

(* either spelling of an instruction list; see instructionTuples (Gadget.wl) *)
gadgetDepth[instr : {___, _Rule | _Wolfram`QuantumFramework`QuantumMeasurementOperator, ___}, rest___] := gadgetDepth[instructionTuples[instr], rest]

gadgetDepth[instr_List, nq_Integer] := Length[circuitSchedule[instr, nq]]

(* either spelling of an instruction list; see instructionTuples (Gadget.wl) *)
gadgetGateCounts[instr : {___, _Rule | _Wolfram`QuantumFramework`QuantumMeasurementOperator, ___}, rest___] := gadgetGateCounts[instructionTuples[instr], rest]

(* keyed by the public names, the ones the instruction lists themselves use *)
gadgetGateCounts[instr_List] := KeyMap[publicOp, Counts[First /@ instr]]

(* either spelling of an instruction list; see instructionTuples (Gadget.wl) *)
gadgetMeasurements[instr : {___, _Rule | _Wolfram`QuantumFramework`QuantumMeasurementOperator, ___}, rest___] := gadgetMeasurements[instructionTuples[instr], rest]

gadgetMeasurements[instr_List] := Count[instr, {"M", _}]

(* either spelling of an instruction list; see instructionTuples (Gadget.wl) *)
gadgetHeralds[instr : {___, _Rule | _Wolfram`QuantumFramework`QuantumMeasurementOperator, ___}, rest___] := gadgetHeralds[instructionTuples[instr], rest]

gadgetHeralds[instr_List] := Count[instr, {"MH", _}]

(* A location is an instruction OR a wait (Got26 Def 10.1), and the waits come from
   the schedule rather than from the list -- which is why this is not Length. *)
(* either spelling of an instruction list; see instructionTuples (Gadget.wl) *)
gadgetLocations[instr : {___, _Rule | _Wolfram`QuantumFramework`QuantumMeasurementOperator, ___}, rest___] := gadgetLocations[instructionTuples[instr], rest]

gadgetLocations[instr_List, nq_Integer] := Length[instr] + Length[circuitIdleSlots[instr, nq]]


(* ---- the bridge to QuantumCircuitOperator ---- *)

(* either spelling of an instruction list; see instructionTuples (Gadget.wl) *)
gadgetCircuitGates[instr : {___, _Rule | _Wolfram`QuantumFramework`QuantumMeasurementOperator, ___}, rest___] := gadgetCircuitGates[instructionTuples[instr], rest]

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

(* either spelling of an instruction list; see instructionTuples (Gadget.wl) *)
gadgetCircuitOperator[instr : {___, _Rule | _Wolfram`QuantumFramework`QuantumMeasurementOperator, ___}, rest___] := gadgetCircuitOperator[instructionTuples[instr], rest]

gadgetCircuitOperator[instr_List] :=
    Wolfram`QuantumFramework`QuantumCircuitOperator[gadgetCircuitGates[instr]]

(* The drawing is the circuit object's own; options pass through, so a caller can
   ask for a different size or ordering without this file knowing about it. *)
(* either spelling of an instruction list; see instructionTuples (Gadget.wl) *)
gadgetDiagram[instr : {___, _Rule | _Wolfram`QuantumFramework`QuantumMeasurementOperator, ___}, rest___] := gadgetDiagram[instructionTuples[instr], rest]

gadgetDiagram[instr_List, opts___] := gadgetCircuitOperator[instr]["Diagram", opts]
