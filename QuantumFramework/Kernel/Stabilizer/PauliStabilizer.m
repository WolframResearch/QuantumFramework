Package["Wolfram`QuantumFramework`"]

PackageExport[PauliStabilizer]
PackageExport[StabilizerStateQ]
PackageScope[PauliStabilizerApply]
PackageScope[PauliStabilizerQ]
PackageScope[ConcretePauliStabilizerQ]
PackageScope[PauliTableauQ]



(* ============================================================================ *)
(* Messages                                                                     *)
(* ============================================================================ *)

PauliStabilizer::nonclifford = "Gate `1` is not a Clifford operation; PauliStabilizer cannot update the tableau. Use Method -> \"TensorNetwork\" or \"Schrodinger\" for non-Clifford circuits."



(* ============================================================================ *)
(* Predicates                                                                   *)
(* ============================================================================ *)

PauliTableauQ[t_] := ArrayQ[t, 3, MatchQ[0 | 1 | -1]] && MatchQ[Dimensions[t], {2, n_, m_}]


(* Phase 3: signs may be polynomials in F_2 over fresh measurement-outcome symbols
   (FangYing23 SymPhase). A "valid" sign is now any expression -- the implementation
   propagates symbols through Clifford gates via BitXor's canonical form. Concrete
   signs are still {-1, 1} integers; symbolic signs look like (1 - 2 s_k) etc. *)

PauliStabilizerQ[PauliStabilizer[KeyValuePattern[{
    "Signs" -> signs_List,
    "Tableau" -> tableau_ ? PauliTableauQ
}] /; Length[signs] == Dimensions[tableau][[3]]]] := True

PauliStabilizerQ[_] := False


(* Strict predicate: signs are concrete (numeric) integers in {-1, 1}.
   Used internally where symbolic signs are not yet supported (e.g. State materialization).
   The {-1, 1} typed-pattern check is moved into the Condition to avoid pattern-engine
   noise (Length/Dimensions::argx) when the structure matches but signs are symbolic. *)
ConcretePauliStabilizerQ[PauliStabilizer[KeyValuePattern[{
    "Signs" -> signs_List,
    "Tableau" -> tableau_ ? PauliTableauQ
}] /; Length[signs] == Dimensions[tableau][[3]] && MatchQ[signs, {(-1 | 1) ...}]]] := True

ConcretePauliStabilizerQ[_] := False


(* ============================================================================ *)
(* A.8 (2026-05-07): StabilizerStateQ public symbol.                            *)
(*                                                                              *)
(* True if `expr` represents an n-qubit stabilizer state. Accepts:              *)
(*   - PauliStabilizer that satisfies the structural predicate                 *)
(*   - StabilizerFrame whose components all reduce to a single PauliStabilizer  *)
(*     (rank-1 case)                                                             *)
(*                                                                              *)
(* Rejects QuantumState by default (stabilizer-state detection from a state    *)
(* vector requires 4^n tomography; use PauliStabilizer[qs] to attempt          *)
(* detection explicitly).                                                       *)
(* ============================================================================ *)

StabilizerStateQ[ps_PauliStabilizer] := PauliStabilizerQ[ps]
StabilizerStateQ[sf_StabilizerFrame] /; StabilizerFrameQ[sf] && sf["Length"] == 1 := True
StabilizerStateQ[_] := False


(* ============================================================================ *)
(* Properties contract                                                          *)
(* ============================================================================ *)

_PauliStabilizer["Properties"] = {
    "Qubits", "Qudits", "GeneratorCount",
    "Signs", "Phase",
    "X", "Z", "Tableau", "Matrix", "p", "TableauPhase",
    "Stabilizer", "StabilizerTableau", "StabilizerX", "StabilizerZ", "StabilizerSigns",
    "Destabilizer", "DestabilizerTableau", "DestabilizerX", "DestabilizerZ", "DestabilizerSigns",
    "PauliForm", "Generators", "Stabilizers", "Destabilizers",
    "PauliStrings", "PauliSymbols", "TableauForm", "StabilizerTableauForm",
    "State", "QuantumState", "Operator", "QuantumOperator", "Circuit", "QuantumCircuit",
    "SymbolicMeasure", "SubstituteOutcomes", "SampleOutcomes",
    "InnerProduct", "Expectation"
}


(* ============================================================================ *)
(* Direct property lookup: if `prop` is already a stored key, return it.        *)
(* ============================================================================ *)

PauliStabilizer[assoc_Association][prop_String] /; KeyExistsQ[assoc, prop] := assoc[prop]


(* ============================================================================ *)
(* Argument-rewrite dispatchers                                                 *)
(* ============================================================================ *)

ps_PauliStabilizer[op_ -> order_] := ps[op, Sequence @@ Flatten[{order}]]

ps_PauliStabilizer["Measure" | "M", qudits___Integer] := ps["M", {qudits}]
ps_PauliStabilizer[qudits__Integer] := ps["M", {qudits}]
ps_PauliStabilizer[qudits : {___Integer}] := ps["M", qudits]
ps_PauliStabilizer[] := ps["M", Range[ps["Qudits"]]]


(* ============================================================================ *)
(* Top-level circuit-application dispatcher                                     *)
(* (preserves Method -> "Stabilizer" wire compatibility for QuantumCircuitOperator)  *)
(* ============================================================================ *)

PauliStabilizerApply[qco_QuantumCircuitOperator, qs : Automatic | _QuantumState | _PauliStabilizer : Automatic] := Fold[
    Function[{state, gate},
        With[{
            rewrittenGate = Replace[gate, {"C", g : "NOT" | "X" | "Z" -> t_, c_, _} :> "C" <> g -> Join[c, t]]
        },
            With[{result = state[rewrittenGate]},
                Which[
                    (* normal Clifford gate update *)
                    PauliStabilizerQ[result], result,
                    (* P[\[Theta]] / T / T\[Dagger] return a Plus -- intentional non-Clifford boundary
                       (Phase 1 baseline; Phase 4 will replace with StabilizerFrame).
                       Continue silently; further gates won't reduce but no message either. *)
                    MatchQ[result, _Plus], result,
                    (* state was already a Plus from a previous P/T -- can't compose further *)
                    MatchQ[state, _Plus], state,
                    (* unknown gate: state[gate] didn't reduce to a PauliStabilizer or Plus *)
                    True, Message[PauliStabilizer::nonclifford, gate]; state
                ]
            ]
        ]
    ],
    Replace[qs, {Automatic :> PauliStabilizer[qco["Arity"]], s_QuantumState :> PauliStabilizer[s]}],
    QuantumShortcut[qco]
]
