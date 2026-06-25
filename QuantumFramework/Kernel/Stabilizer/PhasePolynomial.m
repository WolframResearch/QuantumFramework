Package["Wolfram`QuantumFramework`"]

PackageScope[sfPhasePolyFromCircuit]
PackageScope[sfPhasePolyFold]
PackageScope[sfPathVar]



(* ============================================================================ *)
(* Phase-polynomial carrier for the diagonal (CNOT-dihedral) circuit fragment.  *)
(*                                                                              *)
(* A circuit over {T, T^dagger, S, S^dagger, Z, CZ, CCZ, CNOT} (no Hadamard)    *)
(* acting on |+...+> = 2^{-n/2} Sum_x |x> keeps every qubit value E_q(x) a       *)
(* LINEAR F_2 form (CNOT is linear; the diagonal gates touch only the phase),    *)
(* so the output map E : F_2^n -> F_2^n is a linear bijection (CNOT generates     *)
(* GL(n, F_2)) and the accumulated phase composes into a single Z_8-valued        *)
(* pseudo-Boolean polynomial phi(x). The amplitude is then ONE closed-form term  *)
(*    <y | U | +...+> = 2^{-n/2} omega^{phi(E^{-1} y)},   omega = e^{i pi/4},     *)
(* because the constraint E(x) = y has the unique solution x* = E^{-1} y. The     *)
(* whole state is rank 1 as a phase-polynomial object (one affine support plus    *)
(* one low-degree phase polynomial), so any amplitude costs one F_2 solve and one *)
(* substitution, flat in n, where the dense state vector grows as 2^n.            *)
(*                                                                              *)
(* This file ports the diagonal subset of the standalone path-sum build          *)
(* (pathsum.wl / diagonal.wl) into the kernel as the build-time carrier for the   *)
(* StabilizerFrame "PhasePolynomial" backend. The full-angle convention matches   *)
(* the gate level: T = P[pi/4] contributes omega (phase += 1 in units of pi/4),   *)
(* S contributes 2, Z contributes 4, T^dagger contributes 7 (= -1 mod 8), S^      *)
(* dagger contributes 6, CZ adds 4 x_a x_b, CCZ adds 4 x_a x_b x_c, and CNOT      *)
(* j->k sets x_k = x_j XOR x_k (the only value update; the phase is untouched).   *)
(* ============================================================================ *)

(* Path variables sfPathVar[i], one per qubit. Idempotent (0/1-valued), so a      *)
(* multilinear pseudo-Boolean polynomial in them is reduced by x_i^k -> x_i.       *)
sfIdem[p_] := Expand[p] /. Power[sfPathVar[i_], _] :> sfPathVar[i]
sfF2[p_]   := PolynomialMod[sfIdem[p], 2]
sfZ8[p_]   := PolynomialMod[sfIdem[p], 8]

(* Boolean (XOR) F_2 polynomial -> integer-valued pseudo-Boolean form, using       *)
(* a (+) b = a + b - 2 a b on the F_2 monomials, so a linear F_2 value becomes the *)
(* integer 0/1 quantity the diagonal phase rules expect.                           *)
sfBoolToInt[p_] := sfIdem @ Fold[#1 + #2 - 2 #1 #2 &, 0, MonomialList[sfF2[p]]]

(* Fold a single diagonal gate (read in QuantumShortcut form) onto the carrier     *)
(* <|"val" -> <|q -> F_2 poly|>, "phase" -> Z_8 poly|>. A controlled-Z on any      *)
(* number of qubits applies -1 (phase 4) exactly when every involved qubit is 1,   *)
(* so CZ and CCZ share one rule (the product of the involved integer values). Any  *)
(* spec outside the diagonal+CNOT set returns $Failed, which the builder reads as  *)
(* "this circuit is not in the diagonal fragment" and the apply path falls back.   *)
sfPhasePolyFold[$Failed, _] := $Failed
sfPhasePolyFold[st_Association, spec_] := Module[{val = st["val"], phase = st["phase"]},
    Replace[spec, {
        ("T" -> {q_Integer})              :> <|st, "phase" -> sfIdem[phase + sfBoolToInt[val[q]]]|>,
        (SuperDagger["T"] -> {q_Integer}) :> <|st, "phase" -> sfIdem[phase + 7 sfBoolToInt[val[q]]]|>,
        ("S" -> {q_Integer})              :> <|st, "phase" -> sfIdem[phase + 2 sfBoolToInt[val[q]]]|>,
        (SuperDagger["S"] -> {q_Integer}) :> <|st, "phase" -> sfIdem[phase + 6 sfBoolToInt[val[q]]]|>,
        ("Z" -> {q_Integer})              :> <|st, "phase" -> sfIdem[phase + 4 sfBoolToInt[val[q]]]|>,
        ("C"["Z" -> {t_Integer}, ctrls : {___Integer}, {}]) :>
            <|st, "phase" -> sfIdem[phase + 4 (Times @@ (sfBoolToInt[val[#]] & /@ Append[ctrls, t]))]|>,
        ("C"["NOT" -> {k_Integer}, {j_Integer}, {}]) :>
            <|st, "val" -> Append[val, k -> sfF2[val[j] + val[k]]]|>,
        _ :> $Failed
    }]
]

(* Build the phase-poly backend from a diagonal circuit applied to |+...+>.        *)
(* The |+...+> H-layer is prepended implicitly by seeding val[q] = sfPathVar[q]     *)
(* (an H on |0> with a zero phase carry), so the input circuit must contain NO     *)
(* Hadamard. Returns a StabilizerFrame[<|"PhasePolynomial" -> ...|>] backend, or a  *)
(* Failure (from Enclose) when the circuit is not in the diagonal fragment (an      *)
(* interior H or any non-diagonal gate makes sfPhasePolyFold return $Failed, an     *)
(* unshortcuttable circuit fails the ConfirmAssert), so the caller falls back.      *)
sfPhasePolyFromCircuit[qco_QuantumCircuitOperator] := Enclose @ Module[
    {specs, n, st, mat},
    specs = QuantumShortcut[qco];
    ConfirmAssert[ListQ[specs] && FreeQ[specs, _Missing | _Failure]];
    n = Replace[Max[qco["Arity"], qco["Max"]], Except[_Integer ? Positive] -> $Failed];
    ConfirmAssert[IntegerQ[n] && n > 0];
    st = Fold[sfPhasePolyFold, <|"val" -> Association[Table[q -> sfPathVar[q], {q, n}]], "phase" -> 0|>, specs];
    ConfirmAssert[AssociationQ[st]];
    mat = SparseArray @ Table[Mod[Coefficient[sfF2[st["val"][q]], sfPathVar[i]], 2], {q, n}, {i, n}];
    (* The CNOT-dihedral output map is always a bijection; the check is defensive. *)
    ConfirmAssert[Mod[Det[mat], 2] === 1];
    StabilizerFrame[<|
        "PhasePolynomial" -> sfZ8[st["phase"]],
        "Variables" -> Table[sfPathVar[i], {i, n}],
        "Qubits" -> n,
        "LinearMap" -> mat,
        "InverseMap" -> Inverse[mat, Modulus -> 2],
        "Scale" -> 1 / Sqrt[2] ^ n,
        "Reference" -> PauliStabilizer[n]
    |>]
]
