Package["Wolfram`QuantumFramework`"]



(* ============================================================================ *)
(* Phase 7.1 (2026-05-06): Hybrid interop UpValues.                             *)
(*                                                                              *)
(* Cross-head dispatch so QF measurement / channel operators consume a          *)
(* PauliStabilizer or StabilizerFrame natively without forcing the caller to    *)
(* materialize via ps["State"] (which costs O(2^n)). Mirrors the existing       *)
(* PauliStabilizer integration UpValues at Stabilizer/Conversions.m:139-142.    *)
(*                                                                              *)
(* Design rationale (Phase 6 review with N. Murzin): UpValues attached to       *)
(* PauliStabilizer / StabilizerFrame -- not a Picture flag, not a QuantumBasis  *)
(* wrapper. Both alternatives would route through QuantumBasis machinery        *)
(* (KroneckerProduct[Output, Input] + MatrixInverse[ReducedMatrix]) and pay     *)
(* O(2^n)..O(8^n), defeating the formalism's O(n^2) advantage. UpValues stay   *)
(* in the tableau when they can.                                                *)
(*                                                                              *)
(* Dispatch ladder for qmo[ps]:                                                 *)
(*   1. QMO basis is a Pauli string  -> ps["M", pauli], stays in tableau.       *)
(*   2. QMO basis is non-Pauli       -> Phase 7.2: decompose in stabilizer      *)
(*                                       frame. Currently emits                  *)
(*                                       PauliStabilizer::nonpaulibasis and    *)
(*                                       falls back to                          *)
(*                                       PauliStabilizerApply[                  *)
(*                                         QuantumCircuitOperator[qmo], ps]    *)
(*                                       (the legacy generic path).             *)
(* ============================================================================ *)


(* ============================================================================ *)
(* Helper: extract a Pauli string label from a QMO, when possible.              *)
(*                                                                              *)
(* Returns the string for a Pauli-eigenbasis measurement (e.g. "XZZXI"),        *)
(* or Missing[] otherwise. Used to gate the fast path.                          *)
(* ============================================================================ *)

PackageScope[stabilizerPauliLabelFromQMO]

stabilizerPauliLabelFromQMO[qmo_QuantumMeasurementOperator] := Module[{op, label},
    op = qmo["Operator"];
    label = op["Label"];
    Which[
        StringQ[label] && StringMatchQ[label, RegularExpression["^-?[IXYZ]+$"]],
            label,
        True,
            Missing["NonPauliBasis"]
    ]
]


(* ============================================================================ *)
(* PauliStabilizer::nonpaulibasis -- fired when a non-Pauli QMO/channel acts    *)
(* on a stabilizer-form input. Phase 7.1 falls back to the generic              *)
(* PauliStabilizerApply circuit-conversion path; Phase 7.2 will instead         *)
(* return a StabilizerFrame whose components are the Pauli decomposition of     *)
(* the basis vectors, keeping the cost at O(rank * n^2) rather than O(2^n).    *)
(* ============================================================================ *)

PauliStabilizer::nonpaulibasis = "Hybrid interop: the measurement/channel basis is non-Pauli; falling back to the generic circuit path. Phase 7.2 will route this through a StabilizerFrame decomposition."


(* ============================================================================ *)
(* QuantumMeasurementOperator[qmo][ps_PauliStabilizer] (UpValue on              *)
(*   PauliStabilizer)                                                           *)
(*                                                                              *)
(* Fast path: when the QMO's Operator label parses as a Pauli string, route    *)
(* directly to ps["M", pauliString] which is the existing AG measurement       *)
(* primitive (Stabilizer/PauliMeasure.m). Stays in the tableau, O(n^2).        *)
(*                                                                              *)
(* Fallback: emit ::nonpaulibasis info message and use the legacy generic       *)
(* path PauliStabilizerApply[QuantumCircuitOperator[qmo], ps] which converts    *)
(* the QMO to a circuit and folds gates over the tableau. This was the         *)
(* unconditional behavior before Phase 7.1 (Conversions.m:140 prior to refactor)*)
(* ============================================================================ *)

qmo_QuantumMeasurementOperator[ps_PauliStabilizer ? ConcretePauliStabilizerQ] ^:= Module[{label},
    label = stabilizerPauliLabelFromQMO[qmo];
    If[StringQ[label],
        ps["M", label],
        Message[PauliStabilizer::nonpaulibasis];
        PauliStabilizerApply[QuantumCircuitOperator[qmo], ps]
    ]
]


(* StabilizerFrame: same dispatch shape but no native "M" method on a frame    *)
(* yet (Phase 7.2). For now, materialize the frame and apply the qmo on the    *)
(* materialized state.                                                          *)

qmo_QuantumMeasurementOperator[sf_StabilizerFrame] ^:= (
    Message[PauliStabilizer::nonpaulibasis];
    qmo[sf["State"]]
) /; StabilizerFrameQ[sf]


(* ============================================================================ *)
(* QuantumChannel[qc][ps_PauliStabilizer] / [sf_StabilizerFrame]                *)
(*                                                                              *)
(* Phase 7.1 fallback: materialize and dispatch with ::nonpaulibasis info.      *)
(* Phase 7.2 will identify Clifford channels (BitFlip, PhaseFlip, Pauli         *)
(* mixtures) and route them through the existing tableau-level gate updates    *)
(* without any 2^n materialization.                                             *)
(* ============================================================================ *)

qc_QuantumChannel[ps_PauliStabilizer ? ConcretePauliStabilizerQ] ^:= (
    Message[PauliStabilizer::nonpaulibasis];
    qc[ps["State"]]
)

qc_QuantumChannel[sf_StabilizerFrame] ^:= (
    Message[PauliStabilizer::nonpaulibasis];
    qc[sf["State"]]
) /; StabilizerFrameQ[sf]
