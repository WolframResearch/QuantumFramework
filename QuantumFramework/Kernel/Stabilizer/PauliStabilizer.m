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


(* Signs may be polynomials in F_2 over fresh measurement-outcome symbols
   (FangYing23 SymPhase). A "valid" sign is any expression -- the implementation
   propagates symbols through Clifford gates via BitXor's canonical form. Concrete
   signs are still {-1, 1} integers; symbolic signs look like (1 - 2 s_k) etc. *)

(* Two storage shapes are accepted: the canonical rank-3 "Tableau" + "Signs"    *)
(* form, and the bit-packed form (Stabilizer/Packed.m) with "PackedX"/"PackedZ" *)
(* integer-row lists + "Qubits". A single Which-body inspects which keys are    *)
(* present so the structural conditions never run against the wrong shape.       *)
PauliStabilizerQ[PauliStabilizer[a_Association]] := Which[
    KeyExistsQ[a, "PackedX"],
        (* chunk-major: "PackedX"/"PackedZ" are {#chunks} lists of length-2n      *)
        (* machine-int vectors; "Signs" has the 2n entries.                       *)
        MatchQ[a["PackedX"], {__List}] && MatchQ[a["PackedZ"], {__List}] && MatchQ[a["Signs"], _List] &&
            IntegerQ[Lookup[a, "Qubits"]] &&
            Length[a["PackedX"]] == Length[a["PackedZ"]] &&
            Length[First[a["PackedX"]]] == Length[a["Signs"]],
    KeyExistsQ[a, "Tableau"] && KeyExistsQ[a, "Signs"],
        MatchQ[a["Signs"], _List] && PauliTableauQ[a["Tableau"]] &&
            Length[a["Signs"]] == Dimensions[a["Tableau"]][[3]],
    True, False
]

PauliStabilizerQ[_] := False


(* Strict predicate: signs are concrete (numeric) integers in {-1, 1}.
   Used internally where symbolic signs are not yet supported (e.g. State materialization).
   The {-1, 1} typed-pattern check is moved into the Condition to avoid pattern-engine
   noise (Length/Dimensions::argx) when the structure matches but signs are symbolic. *)
ConcretePauliStabilizerQ[ps : PauliStabilizer[a_Association]] :=
    PauliStabilizerQ[ps] && MatchQ[a["Signs"], {(-1 | 1) ...}]

ConcretePauliStabilizerQ[_] := False


(* ============================================================================ *)
(* StabilizerStateQ: True if `expr` represents an n-qubit stabilizer state.    *)
(*                                                                              *)
(* Accepts:                                                                     *)
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
        (* Short-circuit: once a non-Clifford gate has aborted the fold,        *)
        (* propagate $Failed without firing the message again for every         *)
        (* remaining gate.                                                       *)
        If[state === $Failed,
            $Failed,
            With[{
                rewrittenGate = Replace[gate, "C"[g : "NOT" | "X" | "Z" -> t_, c_, _] :> "C" <> g -> Join[c, t]]
            },
                With[{result = state[rewrittenGate]},
                    Which[
                        (* Clifford gate: PauliStabilizer in, PauliStabilizer out. *)
                        PauliStabilizerQ[result], result,
                        (* Non-Clifford gate that returns a StabilizerFrame (e.g. P[\[Theta]], T, T\[Dagger]). *)
                        StabilizerFrameQ[result], result,
                        (* Gate produced a Plus (superposition of stabilizer states): carry it through. *)
                        MatchQ[result, _Plus], result,
                        (* state was already a Plus and didn't distribute over the gate; pass it through. *)
                        MatchQ[state, _Plus], state,
                        (* Unknown / unsupported gate. *)
                        True, Message[PauliStabilizer::nonclifford, gate]; $Failed
                    ]
                ]
            ]
        ]
    ],
    Replace[qs, {Automatic :> PauliStabilizer[qco["Arity"]], s_QuantumState :> PauliStabilizer[s]}],
    QuantumShortcut[qco]
]
