Package["Wolfram`QuantumFramework`"]

PackageExport[StabilizerFrame]
PackageScope[StabilizerFrameQ]
PackageScope[frameSandwich]



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

(* Compressed phase-polynomial backend (the diagonal-fragment closed form). It    *)
(* carries no "Components" list; the master reads dispatch on which key is present.*)
StabilizerFrameQ[StabilizerFrame[KeyValuePattern[{
    "PhasePolynomial" -> _, "Variables" -> _List, "Qubits" -> _Integer,
    "LinearMap" -> _, "InverseMap" -> _, "Scale" -> _
}]]] := True

StabilizerFrameQ[_] := False


(* Property contract *)
_StabilizerFrame["Properties"] = {
    "Components", "Coefficients", "Stabilizers",
    "Length", "Qubits", "GeneratorCount",
    "StateVector", "State", "InnerProduct",
    "Expectation", "Amplitude", "Compress"
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

(* Phase-poly backend: "Qubits" is a stored key (direct lookup above); the         *)
(* conditioned rules keep the component-only "Qubits"/"Qudits" rules (which read    *)
(* f["Components"]) from ever firing on a phase-poly frame.                         *)
f_StabilizerFrame["Qubits" | "Qudits"] /; KeyExistsQ[First[f], "PhasePolynomial"] := First[f]["Qubits"]

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
    Which[
        KeyExistsQ[assoc, "PhasePolynomial"],
        (* Phase-poly backend: each amplitude is one closed-form term (flat in n).  *)
        (* The full vector is built only for cross-checking / small-n materialization.*)
        With[{n = assoc["Qubits"]},
            SparseArray[Table[phasePolyAmplitude[assoc, IntegerDigits[k, 2, n]], {k, 0, 2 ^ n - 1}], 2 ^ n]
        ],
        KeyExistsQ[assoc, "Paulis"],
        (* Coherent: materialize the reference once, relate the rest by Paulis. *)
        Module[{ref = assoc["Components"][[1, 2]], n, v1},
            n = ref["Qubits"];
            v1 = ref["State"]["StateVector"];
            Total @ MapThread[
                #1[[1]] * sfApplyPauli[#2, v1, n] &,
                {assoc["Components"], assoc["Paulis"]}
            ]
        ],
        True,
        (* No relating Paulis: materialize each component independently. *)
        Total[#1 * #2["State"]["StateVector"] & @@@ assoc["Components"]]
    ]
]

f_StabilizerFrame["State"] := QuantumState @ f["StateVector"]


(* ============================================================================ *)
(* Phase-polynomial backend: closed-form computational amplitude.               *)
(*                                                                              *)
(* <y | U | +...+> = Scale * omega^{phi(M^{-1} y)}, omega = e^{i pi/4}, where M  *)
(* is the F_2 output linear map (qubit q value = Sum_i M[q,i] x[i]) and phi the  *)
(* Z_8 phase polynomial. The preimage x* = M^{-1} y is the unique F_2 solution   *)
(* of M x = y, so the amplitude is one F_2 solve plus one substitution -- no 2^n *)
(* enumeration. For an integer phi the result is an exact eighth root of unity   *)
(* times Scale; a (would-be) symbolic phase falls back to Exp[I pi/4 phi].       *)
(* ============================================================================ *)

phasePolyAmplitude[assoc_Association, y : {(0 | 1) ..}] := With[{
    ph = assoc["PhasePolynomial"] /. Thread[assoc["Variables"] -> LinearSolve[assoc["LinearMap"], y, Modulus -> 2]],
    scale = assoc["Scale"]
},
    scale * If[IntegerQ[ph], Exp[I Pi / 4] ^ Mod[ph, 8], Exp[I (Pi / 4) ph]]
]


(* ============================================================================ *)
(* Observables: Pauli expectation and amplitudes without densifying the frame.   *)
(*                                                                              *)
(* A gate-built frame is |psi> = Sum_i c_i P_i |ref>, with P_i the relating      *)
(* Pauli of component i and |ref> the reference component (component 1). For a    *)
(* Hermitian Pauli observable P,                                                 *)
(*   <psi|P|psi> = Sum_ij conj(c_i) c_j <ref| P_i^dagger P P_j |ref>,            *)
(* and P_i^dagger P P_j is a single Pauli times a scalar in {1,i,-1,-i}; the     *)
(* surviving <ref| Pauli |ref> is exactly stabilizerExpectation on the one       *)
(* reference state (a value in {0,+1,-1}). The cost is the number of component    *)
(* pairs times a poly(n) stabilizer expectation, never a 2^n state vector.       *)
(*                                                                              *)
(* The symplectic-Pauli kernel below encodes an operator as {x, z, s}, meaning    *)
(* s X^x Z^z with s a scalar and x, z length-n bit lists. A relating Pauli       *)
(* {xb, zb, coeff} is the operator coeff I^(xb.zb) X^xb Z^zb; a Hermitian Pauli   *)
(* string (e.g. "XYZ", optional leading "-") is i^(x.z) X^x Z^z.                 *)
(* ============================================================================ *)

sfPauliMul[{x1_, z1_, s1_}, {x2_, z2_, s2_}] :=
    {BitXor[x1, x2], BitXor[z1, z2], s1 s2 (-1) ^ (z1 . x2)}

sfPauliDag[{x_, z_, s_}] := {x, z, Conjugate[s] (-1) ^ (x . z)}

sfRelatingToSymp[{xb_, zb_, coeff_}] := {xb, zb, coeff I ^ (xb . zb)}

sfStringToSymp[str_String] := With[{
    sign = If[StringStartsQ[str, "-"], -1, 1],
    chars = Characters[StringDelete[str, StartOfString ~~ "-"]]
},
    With[{
        x = Replace[chars, {"I" -> 0, "X" -> 1, "Y" -> 1, "Z" -> 0}, {1}],
        z = Replace[chars, {"I" -> 0, "X" -> 0, "Y" -> 1, "Z" -> 1}, {1}]
    },
        {x, z, sign I ^ (x . z)}
    ]
]

sfUnsignedString[x_, z_] := StringJoin @ MapThread[
    Replace[{#1, #2}, {{0, 0} -> "I", {1, 0} -> "X", {1, 1} -> "Y", {0, 1} -> "Z"}] &, {x, z}]

(* <ref| (s X^x Z^z) |ref> = s i^(-(x.z)) <ref|Phat|ref>, last factor in {0,+-1}. *)
sfSympExpectation[{x_, z_, s_}, ref_] :=
    s I ^ (- (x . z)) stabilizerExpectation[ref, sfUnsignedString[x, z]]

(* <psiA|P|psiB> for two frames that share the reference component (component 1). *)
frameSandwich[fA_StabilizerFrame, P_String, fB_StabilizerFrame] := With[{
    a = fA["Coefficients"], pa = First[fA]["Paulis"],
    b = fB["Coefficients"], pb = First[fB]["Paulis"],
    ref = fA["Components"][[1, 2]], pg = sfStringToSymp[P]
},
    Total @ Flatten @ Table[
        Conjugate[a[[i]]] b[[j]] sfSympExpectation[
            sfPauliMul[sfPauliDag[sfRelatingToSymp[pa[[i]]]], sfPauliMul[pg, sfRelatingToSymp[pb[[j]]]]],
            ref],
        {i, Length[a]}, {j, Length[b]}]
]

(* Dense Pauli-string matrix, used only by the hand-built (no relating Paulis)   *)
(* expectation fallback; the gate-built path never forms a 2^n matrix.           *)
sfPauliStringMatrix[str_String] := With[{
    sign = If[StringStartsQ[str, "-"], -1, 1],
    mats = Replace[Characters[StringDelete[str, StartOfString ~~ "-"]],
        {"I" -> PauliMatrix[0], "X" -> PauliMatrix[1], "Y" -> PauliMatrix[2], "Z" -> PauliMatrix[3]}, {1}]
},
    sign If[Length[mats] == 1, First[mats], KroneckerProduct @@ mats]
]

(* <psi|P|psi>. Gate-built frame: poly(n) Pauli sandwich. Hand-built frame: *)
(* dense fallback through the materialized state vector. A wrong-length Pauli      *)
(* string fails cleanly with the same dimension message as the bare stabilizer    *)
(* path (InnerProduct.m), rather than cascading into Thread/Dot errors. *)
f_StabilizerFrame["Expectation", P_String] /; StringMatchQ[P, RegularExpression["^-?[IXYZ]+$"]] :=
    With[{n = StringLength[StringDelete[P, StartOfString ~~ "-"]]},
        Which[
            n != f["Qubits"],
                Message[PauliStabilizer::expectationdim, n, f["Qubits"]]; $Failed,
            KeyExistsQ[First[f], "Paulis"],
                frameSandwich[f, P, f],
            True,
                With[{v = Normal @ f["StateVector"]}, Conjugate[v] . (sfPauliStringMatrix[P] . v)]
        ]
    ]

(* <y|psi> for a computational basis label y (qubit 1 = most significant). The   *)
(* reference is materialized once (O(2^n)); each amplitude is then a chi-term     *)
(* relating-Pauli sum. A label of the wrong length fails cleanly. *)
f_StabilizerFrame["Amplitude", y : {(0 | 1) ..}] /; KeyExistsQ[First[f], "Paulis"] :=
    If[ Length[y] != f["Qubits"],
        Message[PauliStabilizer::expectationdim, Length[y], f["Qubits"]]; $Failed,
        With[{a = f["Coefficients"], pa = First[f]["Paulis"], ref = f["Components"][[1, 2]]},
            With[{vref = Normal @ ref["State"]["StateVector"]},
                Total @ Table[
                    With[{xb = pa[[i, 1]], zb = pa[[i, 2]], coeff = pa[[i, 3]]},
                        a[[i]] coeff I ^ (xb . zb) (-1) ^ (zb . BitXor[y, xb]) vref[[FromDigits[BitXor[y, xb], 2] + 1]]
                    ],
                    {i, Length[a]}]
            ]
        ]
    ]

(* Phase-poly backend: <y|psi> is the closed-form term, flat in n (one F_2 solve  *)
(* + one substitution), the headline of this backend. A wrong-length label fails   *)
(* cleanly with the same dimension message as the component path. *)
f_StabilizerFrame["Amplitude", y : {(0 | 1) ..}] /; KeyExistsQ[First[f], "PhasePolynomial"] :=
    If[ Length[y] != f["Qubits"],
        Message[PauliStabilizer::expectationdim, Length[y], f["Qubits"]]; $Failed,
        phasePolyAmplitude[First[f], y]
    ]


(* "Compress": span-dimension rank cap. The component vectors u_i = P_i|ref>     *)
(* span a space of dimension <= 2^n, but a gate-built frame can carry many more   *)
(* components than that (each non-Clifford gate doubles the count). Compress       *)
(* returns an equivalent frame on a maximal linearly independent subset of        *)
(* components, with exactly-recomputed coefficients, so |psi> is unchanged. The    *)
(* Gram matrix G_ij = <u_i|u_j> is the same poly(n) sandwich (identity            *)
(* observable); its pivot columns are a maximal independent subset of the u_i     *)
(* (for a Gram G = U^dagger U, column dependencies of G equal those of U), and    *)
(* the kept coefficients solve G_BB d = G_B. c exactly. Rank is taken exactly via *)
(* RowReduce, never numeric MatrixRank. A hand-built frame is returned unchanged. *)
sfGramMatrix[f_StabilizerFrame] := With[{pa = First[f]["Paulis"], ref = f["Components"][[1, 2]]},
    Table[
        sfSympExpectation[sfPauliMul[sfPauliDag[sfRelatingToSymp[pa[[i]]]], sfRelatingToSymp[pa[[j]]]], ref],
        {i, Length[pa]}, {j, Length[pa]}]
]

(* Pivot column of each row of the reduced echelon form (the index of its leading *)
(* nonzero); zero rows yield the sentinel 0, which is then dropped. *)
sfPivotColumns[m_] := DeleteCases[
    Function[row, SelectFirst[Range[Length[row]], ! PossibleZeroQ[row[[#]]] &, 0]] /@ RowReduce[m],
    0]

f_StabilizerFrame["Compress"] := If[! KeyExistsQ[First[f], "Paulis"], f,
    Module[{g = sfGramMatrix[f], c = f["Coefficients"], pivots, d},
        pivots = sfPivotColumns[g];
        (* Simplify canonicalizes the reduced coefficients: LinearSolve returns them *)
        (* in an unreduced algebraic form that depends on how many pivots it sees     *)
        (* (e.g. (4 + 4 (-1)^(3/4))/8 on a first pass vs the equal (1 + (-1)^(3/4))/2 *)
        (* once already minimal), which would break structural idempotency            *)
        (* (Compress[Compress[f]] === Compress[f]) and the Components-based Equal.     *)
        d = Simplify @ LinearSolve[g[[pivots, pivots]], g[[pivots, All]] . c];
        StabilizerFrame[<|
            "Components" -> Thread[{d, f["Stabilizers"][[pivots]]}],
            "Paulis" -> First[f]["Paulis"][[pivots]]
        |>]
    ]
]


(* ============================================================================ *)
(* Gate updates: distribute over components                                     *)
(* ============================================================================ *)

(* args is restricted to integer targets so a trailing GATE name (e.g. f["H","S"]) *)
(* is NOT captured here -- it falls through to the multi-gate variadic rule below.  *)
(* Without this restriction sfConjugatePauli would receive a string as a qubit      *)
(* index and corrupt the relating Paulis.                                           *)
f_StabilizerFrame[gate_String, args : ___Integer] := With[{assoc = First[f]},
    If[ KeyExistsQ[assoc, "Paulis"],
        StabilizerFrame[<|
            "Components" -> ({#1, #2[gate, args]} & @@@ assoc["Components"]),
            "Paulis" -> (sfConjugatePauli[#, gate, {args}] & /@ assoc["Paulis"])
        |>],
        StabilizerFrame[{#1, #2[gate, args]} & @@@ assoc["Components"]]
    ]
]

(* SuperDagger gates (S^\[Dagger], V^\[Dagger], T^\[Dagger]) *)
f_StabilizerFrame[SuperDagger[gate_String], args : ___Integer] := With[{assoc = First[f]},
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
(*   P[\[Theta]] |s> = ((1 + e^{i\[Theta]})/2) |s> + ((1 - e^{i\[Theta]})/2) Z_q |s>.          *)
(* ============================================================================ *)

f_StabilizerFrame["P"[phase_], q_Integer] := With[{c = Exp[I phase], assoc = First[f]},
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

f_StabilizerFrame["T", q_Integer] := f["P"[Pi / 4], q]
f_StabilizerFrame[SuperDagger["T"], q_Integer] := f["P"[- Pi / 4], q]


(* op -> order rewrite *)
f_StabilizerFrame[op_ -> order_] := f[op, Sequence @@ Flatten[{order}]]


(* ============================================================================ *)
(* Multi-gate application (mirrors PauliStabilizer's list / variadic forms).    *)
(* A frame has no compiled engine, so fold per gate through the single-gate      *)
(* handlers above; each step updates the relating Paulis coherently. Normalizing *)
(* bare names routes them through the op -> order rewrite into f[gate, target].  *)
(* ============================================================================ *)

f_StabilizerFrame[specs_List] /; stabilizerGateSpecListQ[specs] :=
    Fold[#1[stabilizerNormalizeSpec[#2]] &, f, specs]

f_StabilizerFrame[specs__] /; stabilizerGateSpecSeqQ[{specs}] := f[{specs}]


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

MakeBoxes[f : StabilizerFrame[KeyValuePattern["Components" -> _]] ? StabilizerFrameQ, form_] ^:= With[{nComp = f["Length"]},
    BoxForm`ArrangeSummaryBox["StabilizerFrame",
        f,
        Framed["\[ScriptCapitalF]"],
        {{BoxForm`SummaryItem[{"Components: ", nComp}]},
         {BoxForm`SummaryItem[{"Qubits: ", f["Qubits"]}]}},
        {{BoxForm`SummaryItem[{"First coefficient: ", f["Coefficients"][[1]]}]}},
        form
    ]
]

(* Phase-polynomial backend: a rank-1 compressed frame (one F_2 linear map + one   *)
(* Z_8 phase polynomial). The summary reports the backend, the qubit count, and    *)
(* the phase polynomial rather than a component list.                              *)
MakeBoxes[f : StabilizerFrame[KeyValuePattern["PhasePolynomial" -> _]] ? StabilizerFrameQ, form_] ^:=
    BoxForm`ArrangeSummaryBox["StabilizerFrame",
        f,
        Framed["\[ScriptCapitalF]"],
        {{BoxForm`SummaryItem[{"Backend: ", "PhasePolynomial"}]},
         {BoxForm`SummaryItem[{"Qubits: ", f["Qubits"]}]}},
        {{BoxForm`SummaryItem[{"Phase polynomial: ", f["PhasePolynomial"]}]}},
        form
    ]
