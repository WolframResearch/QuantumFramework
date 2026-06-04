Package["Wolfram`QuantumFramework`"]

PackageExport[StabilizerFrame]
PackageScope[StabilizerFrameQ]



(* ============================================================================ *)
(* StabilizerFrame: superpositions of stabilizer states.                        *)
(*                                                                              *)
(* A `StabilizerFrame` represents a state of the form                           *)
(*       Sum_i c_i |s_i\[RightAngleBracket]                                     *)
(* where each |s_i\[RightAngleBracket] is a stabilizer state and c_i are        *)
(* (possibly symbolic) complex coefficients. This is the natural object for    *)
(* T-gate-rich circuits, magic-state distillation, and stabilizer-rank          *)
(* simulation.                                                                  *)
(*                                                                              *)
(* Reference: Garc\[IAcute]a & Markov, "Simulation of Quantum Circuits via      *)
(* Stabilizer Frames" (Quipu), arxiv:1712.03554, Section 3.                     *)
(*                                                                              *)
(* Internal representation:                                                     *)
(*   StabilizerFrame[<|"Components" -> {{c_1, ps_1}, {c_2, ps_2}, ...}|>]       *)
(*                                                                              *)
(* The user-facing constructor accepts a list of {coefficient, PauliStabilizer} *)
(* pairs and stores them in the Components key.                                 *)
(* ============================================================================ *)


(* ============================================================================ *)
(* Predicate                                                                    *)
(* ============================================================================ *)

StabilizerFrameQ[StabilizerFrame[KeyValuePattern[{
    "Components" -> components_List
}]]] /; AllTrue[components, MatchQ[#, {_, _PauliStabilizer}] &] := True

StabilizerFrameQ[_] := False


(* Property contract *)
_StabilizerFrame["Properties"] = {
    "Components", "Coefficients", "Stabilizers",
    "Length", "Qubits", "GeneratorCount",
    "StateVector", "State", "InnerProduct"
}


(* Method-grade InnerProduct on a frame. *)
f_StabilizerFrame["InnerProduct", other_] := stabilizerInnerProduct[f, other]


(* ============================================================================ *)
(* Constructors                                                                 *)
(* ============================================================================ *)

(* From a list of {coefficient, PauliStabilizer} pairs *)
StabilizerFrame[components : {{_, _PauliStabilizer}..}] :=
    StabilizerFrame[<|"Components" -> components|>]

(* From a single PauliStabilizer (coefficient = 1) *)
StabilizerFrame[ps_PauliStabilizer] := StabilizerFrame[{{1, ps}}]


(* ============================================================================ *)
(* Direct property handlers                                                     *)
(* ============================================================================ *)

StabilizerFrame[assoc_Association][prop_String] /; KeyExistsQ[assoc, prop] := assoc[prop]

f_StabilizerFrame["Coefficients"] := f["Components"][[All, 1]]
f_StabilizerFrame["Stabilizers"] := f["Components"][[All, 2]]
f_StabilizerFrame["Length"] := Length[f["Components"]]
f_StabilizerFrame["Qubits"] := f["Components"][[1, 2]]["Qubits"]
f_StabilizerFrame["Qudits"] := f["Components"][[1, 2]]["Qubits"]
f_StabilizerFrame["GeneratorCount"] := f["Components"][[1, 2]]["GeneratorCount"]


(* ============================================================================ *)
(* Materialization (only practical for small n)                                 *)
(* ============================================================================ *)

f_StabilizerFrame["StateVector"] := Total[#1 * #2["State"]["StateVector"] & @@@ f["Components"]]

f_StabilizerFrame["State"] := QuantumState @ f["StateVector"]


(* ============================================================================ *)
(* Gate updates: distribute over components                                     *)
(* ============================================================================ *)

f_StabilizerFrame[gate_String, args___] := StabilizerFrame[
    {#1, #2[gate, args]} & @@@ f["Components"]
]

(* SuperDagger gates (S^\[Dagger], V^\[Dagger], T^\[Dagger]) *)
f_StabilizerFrame[SuperDagger[gate_String], args___] := StabilizerFrame[
    {#1, #2[SuperDagger[gate], args]} & @@@ f["Components"]
]


(* ============================================================================ *)
(* Non-Clifford gate: P[\[Theta]] doubles the frame size                        *)
(*                                                                              *)
(* P[\[Theta]] = diag(1, exp(I \[Theta])). Acting on a stabilizer state |s>:    *)
(*   P[\[Theta]] |s> = ((1 + e^{i\[Theta]/2})/2) |s> + ((1 - e^{i\[Theta]/2})/2) Z_q |s>     *)
(*  (up to a global phase factor).                                              *)
(* ============================================================================ *)

f_StabilizerFrame["P"[phase_], q_Integer] := With[{c = Exp[I phase / 2]},
    StabilizerFrame[Flatten[
        Function[{coeff, ps}, {{(1 + c) / 2 coeff, ps}, {(1 - c) / 2 coeff, ps["Z", q]}}] @@@ f["Components"],
        1
    ]]
]

f_StabilizerFrame["T", q_Integer] := f["P"[Pi / 2], q]
f_StabilizerFrame[SuperDagger["T"], q_Integer] := f["P"[- Pi / 2], q]


(* op -> order rewrite *)
f_StabilizerFrame[op_ -> order_] := f[op, Sequence @@ Flatten[{order}]]


(* ============================================================================ *)
(* Equality + arithmetic upvalues                                               *)
(* ============================================================================ *)

StabilizerFrame /: Equal[a_StabilizerFrame, b_StabilizerFrame] := SameQ[a["Components"], b["Components"]]


(* Addition: combine component lists *)
StabilizerFrame /: Plus[a_StabilizerFrame, b_StabilizerFrame] :=
    StabilizerFrame[Join[a["Components"], b["Components"]]]


(* Scalar multiplication *)
StabilizerFrame /: Times[c_, f_StabilizerFrame] /; FreeQ[c, _StabilizerFrame] :=
    StabilizerFrame[{c #1, #2} & @@@ f["Components"]]


(* ============================================================================ *)
(* Formatting                                                                   *)
(* ============================================================================ *)

MakeBoxes[f_StabilizerFrame ? StabilizerFrameQ, form_] ^:= With[{nComp = f["Length"]},
    BoxForm`ArrangeSummaryBox["StabilizerFrame",
        f,
        Framed["\[ScriptCapitalF]"],
        {{BoxForm`SummaryItem[{"Components: ", nComp}]},
         {BoxForm`SummaryItem[{"Qubits: ", f["Qubits"]}]}},
        {{BoxForm`SummaryItem[{"First coefficient: ", f["Coefficients"][[1]]}]}},
        form
    ]
]
