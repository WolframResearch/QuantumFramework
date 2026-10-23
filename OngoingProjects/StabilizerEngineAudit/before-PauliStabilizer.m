Package["Wolfram`QuantumFramework`"]

PackageExport[PauliStabilizer]
PackageExport[StabilizerStateQ]
PackageScope[PauliStabilizerApply]
PackageScope[PauliStabilizerQ]
PackageScope[ConcretePauliStabilizerQ]
PackageScope[PauliTableauQ]
PackageScope[stabilizerGateSpecQ]
PackageScope[stabilizerGateSpecListQ]
PackageScope[stabilizerGateSpecSeqQ]
PackageScope[stabilizerDefaultTarget]
PackageScope[stabilizerNormalizeSpec]
PackageScope[stabilizerGateFold]



(* ============================================================================ *)
(* Messages                                                                     *)
(* ============================================================================ *)

PauliStabilizer::nonclifford = "Gate `1` cannot be applied to the stabilizer tableau: `2`. Use Method -> \"TensorNetwork\" or \"Schrodinger\" for operations the stabilizer engine cannot track."

PauliStabilizer::badgate = "`1` is not a recognized stabilizer gate name."

PauliStabilizer::aftermeasurement = "Gate `1` follows a measurement. ApplyCircuit returns a measurement's outcomes as an Association of tableaux, so a measurement must be the last element of the circuit."

PauliStabilizer::largegate ="Reading the `1`-qubit gate `2` from its matrix takes time and memory growing as 4^`1`: expect seconds, and minutes for an exact matrix on 12 qubits."

PauliStabilizer::nophasepoly = "`1`, so \"Compress\" -> \"PhasePolynomial\" cannot apply; falling back to the ordinary stabilizer path, which `2`."

PauliStabilizer::wires = "Circuit wires `1` lie below qubit 1, outside the tableau register."

PauliStabilizer::symbolicsigns = "Symbolic signs have no dense state to hand to Method -> `1`; use Method -> \"Stabilizer\" on this tableau."



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
(* Multi-gate application                                                       *)
(*                                                                              *)
(* ps[{spec, ...}] applies a whole circuit, routing to the compiled fast engine *)
(* ps["ApplyCircuit", ...]; ps[spec, spec, ...] is a braceless variadic sugar    *)
(* that forwards to the list form. A "spec" is an arrow rule (gate -> order), a  *)
(* bare Clifford gate name that takes a default target (1-qubit gates -> qubit 1,*)
(* two-qubit gates -> {1, 2}), or a gate call-form ("P"[theta], SuperDagger[..]).*)
(*                                                                              *)
(* These never shadow the existing forms: the list form is structurally distinct *)
(* from the single-string property accessor (ps["X"]) and from the integer-list  *)
(* measurement (ps[{1, 2}], which is the more specific {___Integer} rule above), *)
(* and the variadic form's Length >= 2 + gate-membership guard rejects single    *)
(* properties, ps[gate, target] single-gate calls, and measurement.              *)
(* ============================================================================ *)

stabilizerGateSpecQ[(lhs_ -> _)]                  := stabilizerGateSpecQ[lhs]
stabilizerGateSpecQ[g_String]                     := MemberQ[$stabilizerCliffordNames, g]
stabilizerGateSpecQ["P"[_]]                       := True
stabilizerGateSpecQ[SuperDagger["S" | "V" | "T"]] := True
stabilizerGateSpecQ[_]                            := False

stabilizerGateSpecListQ[specs_List] := Length[specs] >= 1 && AllTrue[specs, stabilizerGateSpecQ]
stabilizerGateSpecSeqQ[specs_List]  := Length[specs] >= 2 && AllTrue[specs, stabilizerGateSpecQ]

(* Default target for a bare gate token: two-qubit gates span {1, 2}, all other *)
(* (single-qubit) gates act on qubit 1.                                          *)
stabilizerDefaultTarget[g_String] := If[MemberQ[$stabilizerTwoQubitNames, g], {1, 2}, 1]
stabilizerDefaultTarget["P"[_]] := 1
stabilizerDefaultTarget[SuperDagger["S" | "V" | "T"]] := 1
stabilizerDefaultTarget[other_] := (Message[PauliStabilizer::badgate, other]; $Failed)

(* Normalize a spec to the `gate -> order` shape the engine already accepts.      *)
(* Arrow / legacy {op, q} / controlled "C"[...] specs pass through unchanged so   *)
(* the compiled encoder (Stabilizer/Compiled.m) sees only shapes it handles.      *)
stabilizerNormalizeSpec[g_String]        := g -> stabilizerDefaultTarget[g]
stabilizerNormalizeSpec["P"[phase_]]     := "P"[phase] -> 1
stabilizerNormalizeSpec[SuperDagger[g_]] := SuperDagger[g] -> 1
stabilizerNormalizeSpec[other_]          := other

ps_PauliStabilizer[specs_List] /; stabilizerGateSpecListQ[specs] :=
    ps["ApplyCircuit", stabilizerNormalizeSpec /@ specs]

ps_PauliStabilizer[specs__] /; stabilizerGateSpecSeqQ[{specs}] := ps[{specs}]


(* ============================================================================ *)
(* Top-level circuit-application dispatcher                                     *)
(* (preserves Method -> "Stabilizer" wire compatibility for QuantumCircuitOperator)  *)
(* ============================================================================ *)

Options[PauliStabilizerApply] = {"Compress" -> None};

(* A tableau's wires are its qubits, indexed from 1, so a gate wire below 1 has  *)
(* nothing to land on and the gate kernels would index past the tableau edge:   *)
(* refused up front, whichever entrance (a circuit on a tableau, Method ->       *)
(* "Stabilizer" on a state or on the default register) reached the engine. The  *)
(* dense circuit path re-seats such a circuit onto its state instead.           *)
PauliStabilizerApply[qco_QuantumCircuitOperator, qs : Automatic | _QuantumState | _PauliStabilizer : Automatic, OptionsPattern[]] :=
    With[{low = Select[qco["FullInputOrder"], # < 1 &]},
        (Message[PauliStabilizer::wires, low]; $Failed) /; low =!= {}
    ]

PauliStabilizerApply[qco_QuantumCircuitOperator, qs : Automatic | _QuantumState | _PauliStabilizer : Automatic, opts : OptionsPattern[]] :=
    (* Opt-in phase-polynomial fast build (Method -> {"Stabilizer", "Compress" ->  *)
    (* "PhasePolynomial"}): for the default |0...0> register and a diagonal circuit *)
    (* (no interior Hadamard, gates in {T,Tdg,S,Sdg,Z,CZ,CCZ,CNOT}) build the rank-1*)
    (* phase-poly backend of U|+...+> directly (the |+...+> H-layer is prepended    *)
    (* internally), never materializing the 2^t frame. sfPhasePolyFromCircuit       *)
    (* returns a Failure for anything outside the fragment, so we fall back to the  *)
    (* ordinary stabilizer trichotomy below. Because the opt-in was explicit, the    *)
    (* fallback names why first (the ordinary path evolves |0...0>, not |+...+>).     *)
    With[{compress = OptionValue[PauliStabilizerApply, {opts}, "Compress"]},
    With[{ppFrame = If[compress === "PhasePolynomial" && qs === Automatic,
        sfPhasePolyFromCircuit[qco], $Failed]},
    If[StabilizerFrameQ[ppFrame], ppFrame,
    (* Explicit "Compress" -> "PhasePolynomial" that did not produce a phase-poly     *)
    (* backend (interior Hadamard / non-diagonal gate, or a non-default input state): *)
    (* state why exactly once, then fall through to the ordinary path. The default    *)
    (* Method -> "Stabilizer" (compress =!= "PhasePolynomial") stays silent.          *)
    (If[compress === "PhasePolynomial",
        If[qs === Automatic,
            Message[PauliStabilizer::nophasepoly,
                "The circuit is not in the diagonal (no-Hadamard) fragment, or a gate's matrix is not the product of the T, S, Z, CZ and CNOT gates its label names",
                "evolves |0...0> rather than the |+...+> the phase polynomial assumes"],
            Message[PauliStabilizer::nophasepoly,
                "An explicit input state was supplied",
                "evolves that state directly, not the |+...+> the phase polynomial assumes"]]];
    With[{
        (* Register size must cover the highest wire INDEX, not the wire count:    *)
        (* a circuit like {"H" -> 2} has Arity 1 but acts on wire 2, and an        *)
        (* Arity-sized register silently drops the gate off the register edge.     *)
        init = Replace[qs, {Automatic :> PauliStabilizer[Max[qco["Arity"], qco["Max"]]], s_QuantumState :> PauliStabilizer[s]}],
        (* Each gate is read from its matrix, not its label (Stabilizer/GateMatrix.m): *)
        (* a Clifford becomes the named gates of its own tableau, a one-qubit        *)
        (* diagonal non-Clifford becomes "P"[phase], and a product of one-qubit gates *)
        (* is read factor by factor. A composite, and a symbolic gate those steps    *)
        (* cannot read, go through the gates their label names, used only when that *)
        (* reproduces the matrix. Anything else is a Missing that the fold refuses.  *)
        gateSpecs = stabilizerCircuitSpecs[qco]
    },
        (* Fast path: a pure-Clifford circuit (every gate in the compiled set)    *)
        (* over a concrete state folds in one compiled kernel call               *)
        (* (Stabilizer/Compiled.m). encodeStabilizerGates returns $Failed for any *)
        (* spec outside the set, so a frame gate or a refused gate falls through  *)
        (* to the per-gate fold.                                                  *)
        With[{gates = If[PauliStabilizerQ[init] && psConcreteFastQ[init], encodeStabilizerGates[gateSpecs], $Failed]},
            If[ ListQ[gates],
                applyCompiledFold[init, gates],
                stabilizerGateFold[init, gateSpecs]
            ]
        ]
    ]
    )
    ]]
    ]


(* ============================================================================ *)
(* The per-gate fold shared by PauliStabilizerApply and ps["ApplyCircuit"].     *)
(* A gate the engine cannot apply stops the fold with one                      *)
(* PauliStabilizer::nonclifford; later gates are not tried, so the message     *)
(* fires once. A spec refused by its matrix (a Missing from                    *)
(* stabilizerCircuitSpecs) is refused before it reaches the state, and is      *)
(* named by its QuantumShortcut token with the reason; the message shortens    *)
(* both. With branches True (ps["ApplyCircuit"]), a measurement's outcomes, an *)
(* Association of tableaux, are a result as well, and a gate after them is     *)
(* refused with PauliStabilizer::aftermeasurement. A measurement's outcomes    *)
(* elsewhere, and a channel, whose result is a mixture of tableaux, are        *)
(* refused with a reason that says where they can be applied.                  *)
(* ============================================================================ *)

stabilizerRefusal[Missing[reason_, op_QuantumOperator]] := {Replace[QuantumShortcut[op], {one_} :> one], stabilizerRefusalReason[reason]}
stabilizerRefusal[gate_] := {gate, "the stabilizer engine has no rule for it"}
stabilizerRefusal[gate_QuantumChannel, _] := {gate, "it is a channel, whose result is a mixture of tableaux; apply it to the tableau itself, as channel[ps]"}
stabilizerRefusal[gate_, _Association] := {gate, "it is a measurement with several outcomes, and only ps[\"ApplyCircuit\", specs] returns outcomes, as tableaux, for a measurement that comes last"}
stabilizerRefusal[gate_, _] := stabilizerRefusal[gate]

stabilizerGateFold[init_, specs_List, branches_ : False] := Fold[
    Function[{state, gate},
        Which[
            state === $Failed, $Failed,
            AssociationQ[state], Message[PauliStabilizer::aftermeasurement, First[stabilizerRefusal[gate]]]; $Failed,
            MissingQ[gate], Message[PauliStabilizer::nonclifford, Sequence @@ stabilizerRefusal[gate]]; $Failed,
            True, With[{result = state[Replace[gate, "C"[g : "NOT" | "X" | "Z" -> t_, c_, _] :> "C" <> g -> Join[c, t]]]},
                Which[
                    (* Clifford gate: PauliStabilizer in, PauliStabilizer out. *)
                    PauliStabilizerQ[result], result,
                    (* Non-Clifford gate that returns a StabilizerFrame (e.g. P[\[Theta]], T, T\[Dagger]). *)
                    StabilizerFrameQ[result], result,
                    (* Gate produced a Plus (superposition of stabilizer states): carry it through. *)
                    MatchQ[result, _Plus], result,
                    (* state was already a Plus and didn't distribute over the gate; pass it through. *)
                    MatchQ[state, _Plus], state,
                    (* A measurement's outcomes, where the caller accepts them. *)
                    TrueQ[branches] && AssociationQ[result] && AllTrue[result, PauliStabilizerQ], result,
                    (* Unknown / unsupported gate, or a result the fold cannot carry. *)
                    True, Message[PauliStabilizer::nonclifford, Sequence @@ stabilizerRefusal[gate, result]]; $Failed
                ]
            ]
        ]
    ],
    init,
    specs
]
