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
(*                                                                              *)
(* Gate-built frames additionally carry a "Paulis" key: one relating Pauli per  *)
(* component (see the relating-Pauli section below) that records, coherently,    *)
(* how each component is obtained from the reference component. It is what makes *)
(* the dense materialization phase-correct across components. Frames without it  *)
(* (hand-built, or combined by Plus) materialize each component independently.   *)
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
(* Relating Paulis: coherent bookkeeping for dense materialization.             *)
(*                                                                              *)
(* A relating Pauli is encoded {xbits, zbits, coeff} with coeff in {1,-1,I,-I}: *)
(* it denotes coeff * (X^{x_q} Z^{z_q}) over qubits q (so {1,1} is Y). For a    *)
(* gate-built frame each component carries one such Pauli P_i with the contract *)
(*   |component_i> = matrix[P_i] . |reference>,                                 *)
(* where the reference is component 1 (relating Pauli the identity). A Z update  *)
(* leaves the stabilizer tableau unchanged and a Clifford update transforms it   *)
(* identically for every component, so all components share the reference's      *)
(* tableau; one reference state vector plus the relating Paulis then reproduce   *)
(* each component coherently, recovering the relative phases that an independent  *)
(* +1 eigenvector per component drops. The two construction primitives are       *)
(* conjugation by a Clifford gate (P -> G.P.G^dagger) and left-multiplication by *)
(* Z_q (a frame doubling); both act only on the one or two qubits the gate       *)
(* touches, so they are precomputed below as small lookup tables.                *)
(* ============================================================================ *)

sfPauliLocal[x_, z_] := PauliMatrix[Replace[{x, z}, {{0, 0} -> 0, {1, 0} -> 1, {1, 1} -> 2, {0, 1} -> 3}]]

(* Decoders: map a conjugated/multiplied local Pauli matrix back to bits+phase. *)
(* Used once, at load, to build the lookup tables; never on the hot path.       *)
sfDecodeLocal1[mat_] := SelectFirst[
    Flatten[Table[{{x, z}, ph}, {x, 0, 1}, {z, 0, 1}, {ph, {1, -1, I, -I}}], 2],
    Chop[mat - #[[2]] sfPauliLocal @@ #[[1]]] == ConstantArray[0, {2, 2}] &
]
sfDecodeLocal2[mat_] := SelectFirst[
    Flatten[Table[{{{xa, za}, {xb, zb}}, ph},
        {xa, 0, 1}, {za, 0, 1}, {xb, 0, 1}, {zb, 0, 1}, {ph, {1, -1, I, -I}}], 4],
    Chop[mat - #[[2]] KroneckerProduct[sfPauliLocal @@ #[[1, 1]], sfPauliLocal @@ #[[1, 2]]]] == ConstantArray[0, {4, 4}] &
]

(* Clifford generators: 1-qubit, and 2-qubit with the first slot = first arg.   *)
sfGate1 = <|
    "H" -> {{1, 1}, {1, -1}} / Sqrt[2],
    "S" -> {{1, 0}, {0, I}}, "Sdg" -> {{1, 0}, {0, -I}},
    "X" -> PauliMatrix[1], "Y" -> PauliMatrix[2], "Z" -> PauliMatrix[3],
    "V" -> {{1 + I, 1 - I}, {1 - I, 1 + I}} / 2, "Vdg" -> {{1 - I, 1 + I}, {1 + I, 1 - I}} / 2
|>;
sfGate2 = <|
    "CNOT" -> {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}},
    "CZ" -> DiagonalMatrix[{1, 1, 1, -1}],
    "SWAP" -> {{1, 0, 0, 0}, {0, 0, 1, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}}
|>;

(* Conjugation P -> G.P.G^dagger and left-multiply Z.P, precomputed per gate.   *)
sfConj1 = Association @ KeyValueMap[Function[{g, m},
    g -> Association @ Flatten @ Table[
        {x, z} -> sfDecodeLocal1[m . sfPauliLocal[x, z] . ConjugateTranspose[m]],
        {x, 0, 1}, {z, 0, 1}]], sfGate1];
sfConj2 = Association @ KeyValueMap[Function[{g, m},
    g -> Association @ Flatten @ Table[
        {{xa, za}, {xb, zb}} -> sfDecodeLocal2[m . KroneckerProduct[sfPauliLocal[xa, za], sfPauliLocal[xb, zb]] . ConjugateTranspose[m]],
        {xa, 0, 1}, {za, 0, 1}, {xb, 0, 1}, {zb, 0, 1}]], sfGate2];
sfTimesZ = Association @ Flatten @ Table[
    {x, z} -> sfDecodeLocal1[PauliMatrix[3] . sfPauliLocal[x, z]], {x, 0, 1}, {z, 0, 1}];

(* Gate-name normalization to the lookup keys. *)
sfCanonicalGate[SuperDagger["S"]] := "Sdg"
sfCanonicalGate[SuperDagger["V"]] := "Vdg"
sfCanonicalGate["CX"] := "CNOT"
sfCanonicalGate[g_] := g

(* Conjugate a relating Pauli by a Clifford gate on qubits qs (1 or 2). *)
sfConjugatePauli[{xb_, zb_, coeff_}, gate_, {q_}] := With[{r = sfConj1[sfCanonicalGate[gate]][{xb[[q]], zb[[q]]}]},
    {ReplacePart[xb, q -> r[[1, 1]]], ReplacePart[zb, q -> r[[1, 2]]], coeff r[[2]]}
]
sfConjugatePauli[{xb_, zb_, coeff_}, gate_, {j_, k_}] := With[{r = sfConj2[sfCanonicalGate[gate]][{{xb[[j]], zb[[j]]}, {xb[[k]], zb[[k]]}}]},
    {ReplacePart[xb, {j -> r[[1, 1, 1]], k -> r[[1, 2, 1]]}], ReplacePart[zb, {j -> r[[1, 1, 2]], k -> r[[1, 2, 2]]}], coeff r[[2]]}
]

(* Left-multiply a relating Pauli by Z_q (used by a frame doubling). *)
sfTimesZPauli[{xb_, zb_, coeff_}, q_] := With[{r = sfTimesZ[{xb[[q]], zb[[q]]}]},
    {ReplacePart[xb, q -> r[[1, 1]]], ReplacePart[zb, q -> r[[1, 2]]], coeff r[[2]]}
]

(* Apply a relating Pauli to a dense state vector (qubit 1 = most significant). *)
sfApplyPauli[{xb_, zb_, coeff_}, v_, n_] := Module[{
    dim = 2 ^ n, xmask = FromDigits[xb, 2], iPow = Total[xb zb], targets, signs
},
    targets = 1 + BitXor[Range[0, dim - 1], xmask];
    signs = Table[(-1) ^ Mod[zb . IntegerDigits[b, 2, n], 2], {b, 0, dim - 1}];
    SparseArray[Thread[targets -> (coeff I ^ iPow) signs Normal[v]], dim]
]


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

f_StabilizerFrame["StateVector"] := With[{assoc = First[f]},
    If[ KeyExistsQ[assoc, "Paulis"],
        (* Coherent: materialize the reference once, relate the rest by Paulis. *)
        Module[{ref = assoc["Components"][[1, 2]], n, v1},
            n = ref["Qubits"];
            v1 = ref["State"]["StateVector"];
            Total @ MapThread[
                #1[[1]] * sfApplyPauli[#2, v1, n] &,
                {assoc["Components"], assoc["Paulis"]}
            ]
        ],
        (* No relating Paulis: materialize each component independently. *)
        Total[#1 * #2["State"]["StateVector"] & @@@ assoc["Components"]]
    ]
]

f_StabilizerFrame["State"] := QuantumState @ f["StateVector"]


(* ============================================================================ *)
(* Gate updates: distribute over components                                     *)
(* ============================================================================ *)

f_StabilizerFrame[gate_String, args___] := With[{assoc = First[f]},
    If[ KeyExistsQ[assoc, "Paulis"],
        StabilizerFrame[<|
            "Components" -> ({#1, #2[gate, args]} & @@@ assoc["Components"]),
            "Paulis" -> (sfConjugatePauli[#, gate, {args}] & /@ assoc["Paulis"])
        |>],
        StabilizerFrame[{#1, #2[gate, args]} & @@@ assoc["Components"]]
    ]
]

(* SuperDagger gates (S^\[Dagger], V^\[Dagger], T^\[Dagger]) *)
f_StabilizerFrame[SuperDagger[gate_String], args___] := With[{assoc = First[f]},
    If[ KeyExistsQ[assoc, "Paulis"],
        StabilizerFrame[<|
            "Components" -> ({#1, #2[SuperDagger[gate], args]} & @@@ assoc["Components"]),
            "Paulis" -> (sfConjugatePauli[#, SuperDagger[gate], {args}] & /@ assoc["Paulis"])
        |>],
        StabilizerFrame[{#1, #2[SuperDagger[gate], args]} & @@@ assoc["Components"]]
    ]
]


(* ============================================================================ *)
(* Non-Clifford gate: P[\[Theta]] doubles the frame size                        *)
(*                                                                              *)
(* P[\[Theta]] = diag(1, exp(I \[Theta])). Acting on a stabilizer state |s>:    *)
(*   P[\[Theta]] |s> = ((1 + e^{i\[Theta]/2})/2) |s> + ((1 - e^{i\[Theta]/2})/2) Z_q |s>     *)
(*  (up to a global phase factor).                                              *)
(* ============================================================================ *)

f_StabilizerFrame["P"[phase_], q_Integer] := With[{c = Exp[I phase / 2], assoc = First[f]},
    If[ KeyExistsQ[assoc, "Paulis"],
        StabilizerFrame[<|
            "Components" -> Flatten[
                Function[{coeff, ps}, {{(1 + c) / 2 coeff, ps}, {(1 - c) / 2 coeff, ps["Z", q]}}] @@@ assoc["Components"],
                1],
            "Paulis" -> Flatten[
                Function[pauli, {pauli, sfTimesZPauli[pauli, q]}] /@ assoc["Paulis"],
                1]
        |>],
        StabilizerFrame[Flatten[
            Function[{coeff, ps}, {{(1 + c) / 2 coeff, ps}, {(1 - c) / 2 coeff, ps["Z", q]}}] @@@ assoc["Components"],
            1
        ]]
    ]
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
StabilizerFrame /: Times[c_, f_StabilizerFrame] /; FreeQ[c, _StabilizerFrame] := With[{assoc = First[f]},
    StabilizerFrame[<|
        "Components" -> ({c #1, #2} & @@@ assoc["Components"]),
        If[KeyExistsQ[assoc, "Paulis"], "Paulis" -> assoc["Paulis"], Nothing]
    |>]
]


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
