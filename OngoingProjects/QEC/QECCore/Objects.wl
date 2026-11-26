(* ::Package:: *)

Package["Wolfram`QuantumFramework`QEC`"]

PackageExport[QECCodeQ]
PackageExport[QECSyndromeCircuitQ]
PackageExport[QECDecoderQ]
PackageExport[QECDetectorModelQ]
PackageExport[QECNoiseModelQ]
PackageExport[QECCatStateQ]
PackageExport[QECPauliMeasurementQ]
PackageExport[QECErrorCorrectionQ]
PackageExport[QECExtendedRectangleQ]
PackageExport[QECFaultTolerantCircuitQ]
PackageExport[QECRegisterQ]
PackageExport[QECThresholdQ]
PackageExport[QECPseudothresholdQ]
PackageExport[QECTransversalGateQ]

PackageScope[$QECObjectHeads]
PackageScope[qecObjectDataQ]


(* ============================================================================ *)
(* What every QEC object shares, written once, the way the framework's objects  *)
(* do it (QuantumChannel.m, QuantumState.m).                                    *)
(*                                                                              *)
(* - A predicate HQ for every head H.  An object whose data has the keys its    *)
(*   head needs is marked valid with System`Private`HoldSetValid the first time *)
(*   it evaluates, so HQ is a flag lookup from then on.                         *)
(* - obj["ValidQ"] and obj["AllProperties"].                                    *)
(* - A property the head does not know: H::undefprop, and the call is left      *)
(*   unevaluated, as QuantumState and QuantumChannel leave theirs.              *)
(*                                                                              *)
(* Loaded after every head is defined, so the catch-all is the last rule of     *)
(* each head and every named property is tried first.                           *)
(* ============================================================================ *)

(* head -> {predicate, the keys a valid object carries} *)
$qecObjects = {
    QECCode -> {QECCodeQ, {"CheckMatrix", "Phases", "Qubits"}},
    QECSyndromeCircuit -> {QECSyndromeCircuitQ, {"Code", "Instructions"}},
    QECDecoder -> {QECDecoderQ, {"Code", "Method"}},
    QECDetectorModel -> {QECDetectorModelQ, {"DetectorMatrix", "ObservableMatrix", "Probabilities"}},
    QECNoiseModel -> {QECNoiseModelQ, {"Level", "Name", "Parameters"}},
    QECCatState -> {QECCatStateQ, {"Instructions", "Qubits"}},
    QECPauliMeasurement -> {QECPauliMeasurementQ, {"Code", "Pauli", "Instructions"}},
    QECErrorCorrection -> {QECErrorCorrectionQ, {"Code", "Instructions"}},
    QECExtendedRectangle -> {QECExtendedRectangleQ, {"Code", "Location", "Instructions"}},
    QECFaultTolerantCircuit -> {QECFaultTolerantCircuitQ, {"Code", "Instructions"}},
    QECRegister -> {QECRegisterQ, {"Blocks", "Code"}},
    QECThreshold -> {QECThresholdQ, {"Method", "Value"}},
    QECPseudothreshold -> {QECPseudothresholdQ, {"Value", "Code", "Noise"}},
    QECTransversalGate -> {QECTransversalGateQ, {"Code", "Gate"}}
};

$QECObjectHeads = Append[Keys[$qecObjects], QECPauli];

SetAttributes[qecObjectDataQ, HoldFirst];
qecObjectDataQ[_[a_Association ? AssociationQ], keys_List] := AllTrue[keys, KeyExistsQ[a, #] &]
qecObjectDataQ[___] := False

$undefpropTemplate = "`1` is not a property of `2`. Use obj[\"Properties\"] for the list.";

Scan[
    Function[rule, With[{h = First[rule], q = rule[[2, 1]], keys = rule[[2, 2]]},

        q::usage = SymbolName[q] <> "[expr] gives True if expr is a valid " <> SymbolName[h] <> " object, and False otherwise.";

        q[obj_h] := System`Private`HoldValidQ[obj] || qecObjectDataQ[obj, keys];
        q[___] := False;

        obj_h /; System`Private`HoldNotValidQ[obj] && qecObjectDataQ[obj, keys] := System`Private`HoldSetValid[obj];

        h[a_Association]["ValidQ"] := q[h[a]];
        h[a_Association]["AllProperties"] := Union[h[a]["Properties"], {"Properties", "AllProperties", "ValidQ"}];

        MessageName[h, "undefprop"] = StringReplace[$undefpropTemplate, "`2`" -> SymbolName[h]];
        h[a_Association][prop_String, ___] /; (Message[MessageName[h, "undefprop"], prop]; False) := Null
    ]],
    $qecObjects
];

(* QECPauli keeps its own predicate, QECPauliQ, which also accepts a Pauli string or row. *)
QECPauli[v_List]["ValidQ"] := QECPauliQ[QECPauli[v]]
QECPauli[v_List]["AllProperties"] := Union[QECPauli[v]["Properties"], {"Properties", "AllProperties", "ValidQ"}]
QECPauli::undefprop = StringReplace[$undefpropTemplate, "`2`" -> "QECPauli"];
QECPauli[_List][prop_String, ___] /; (Message[QECPauli::undefprop, prop]; False) := Null

(* A rule whose property is a pattern is not always placed after the named ones, and the
   catch-all would then speak before a named property is tried.  Moving it to the end of
   every head's rules keeps the message for the properties that do not exist. *)
qecCatchAllQ[rule_] := ! FreeQ[rule, "undefprop"]
Scan[
    With[{sv = SubValues[#]},
        SubValues[#] = Join[Select[sv, ! qecCatchAllQ[#] &], Select[sv, qecCatchAllQ]]
    ] &,
    $QECObjectHeads
];
