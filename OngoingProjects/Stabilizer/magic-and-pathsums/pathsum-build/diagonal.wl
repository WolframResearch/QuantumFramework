(* ::Package:: *)

(* diagonal.wl : CNOT-dihedral (diagonal-fragment) phase-polynomial compressor,   *)
(* plus the isolated single-variable Gauss-sum elimination lemmas (Z_2 / Z_4 /     *)
(* odd cofactor classes + quadratic completion) developed and unit-tested here     *)
(* before they enter the Hadamard-elimination loop in Stage 3.                     *)
(*                                                                                 *)
(* THE DIAGONAL FRAGMENT. A circuit over {T, Sdg, S, Z, CZ, CCZ, CNOT} (no H)       *)
(* acting on |+...+> = 2^{-n/2} Sum_x |x>. With no internal Hadamard, each qubit    *)
(* value E_q(x) is a LINEAR F_2 form (CNOT is linear; the diagonal gates only add   *)
(* to the phase), so E : F_2^n -> F_2^n is a linear BIJECTION (CNOT generates       *)
(* GL(n, F_2)). The amplitude                                                       *)
(*     <y|U|+...+> = 2^{-n/2} Sum_{x : E(x)=y} omega^{phi(x)}                        *)
(* therefore has a UNIQUE solution x* = E^{-1}(y), and collapses to a single        *)
(* closed-form term 2^{-n/2} omega^{phi(xstar)} -- no 2^t branching, no 2^n sum. The    *)
(* whole state is rank 1 as a phase-polynomial object: one affine map + one         *)
(* degree-<=3 Z_8 polynomial psi(y) = phi(E^{-1}(y)).                                *)
(*                                                                                 *)
(* Assumes Wolfram`QuantumFramework`, pathsum.wl, harness.wl loaded.               *)

(* A qubit value E_q(x) = Sum_i M[q,i] x[i]; collect the F_2 coefficient matrix. *)
PathSum`linearMatrix[ps_, n_Integer] := Table[
    Mod[Coefficient[PathSum`f2[ps["val"][q]], x[i]], 2], {q, n}, {i, n}];

(* Compress a diagonal-fragment circuit (gates over {T,Sdg,S,Z,CZ,CCZ,CNOT}, NO H) *)
(* on |+...+>. Returns the rank-1 phase-polynomial object.                          *)
PathSum`hLayer[n_Integer] := Table[{"H", q}, {q, n}];

PathSum`diagonalCompress[n_Integer, gates_List] := Module[
    {ps, M, Minv, yvec, xforms, psi},
    ps = PathSum`build[n, Join[PathSum`hLayer[n], gates]];
    M = PathSum`linearMatrix[ps, n];
    If[Mod[Det[M], 2] === 0,
        (* not a bijection (should not happen for the CNOT-dihedral fragment) *)
        Return[<|"n" -> n, "phi" -> ps["phase"], "M" -> M, "bijective" -> False, "ps" -> ps|>]];
    Minv = Inverse[M, Modulus -> 2];
    yvec = Table[y[j], {j, n}];
    (* x[i] = (M^{-1} y)_i over F_2; substitute the integer value into phi(x). *)
    xforms = Table[PathSum`f2[Minv[[i]] . yvec], {i, n}];
    psi = PathSum`z8[ps["phase"] /. Table[x[i] -> PathSum`boolToInt[xforms[[i]]], {i, n}]];
    <|"n" -> n, "phi" -> ps["phase"], "psi" -> psi, "M" -> M, "Minv" -> Minv,
      "bijective" -> True, "ps" -> ps, "rank" -> 1|>
];

(* Closed-form amplitude <y|U|+...+>, no enumeration: solve E(x)=y over F_2,        *)
(* evaluate omega^{phi(xstar)}. (Robust per-y route via LinearSolve.)                  *)
PathSum`diagAmp[n_Integer, gates_List, y_List] := Module[{ps, M, xstar, ph},
    ps = PathSum`build[n, Join[PathSum`hLayer[n], gates]];
    M = PathSum`linearMatrix[ps, n];
    xstar = LinearSolve[M, y, Modulus -> 2];
    ph = Mod[ps["phase"] /. Thread[Table[x[i], {i, n}] -> xstar], 8];
    (1 / Sqrt[2]^n) PathSum`omega^ph
];

(* Same amplitude from the compressed psi(y) (the rank-1 object), confirming psi    *)
(* is correct: <y| = 2^{-n/2} omega^{psi(y)}.                                        *)
PathSum`diagAmpFromPoly[c_Association, yval_List] :=
    (1 / Sqrt[2]^c["n"]) PathSum`omega^Mod[c["psi"] /. Thread[Table[y[j], {j, c["n"]}] -> yval], 8];

(* Differential matchers: diagonal amplitude === QF Schrodinger of (Hlayer ++ gates),*)
(* exactly, on all 2^n outputs. *)
PathSum`diagMatchesQ[n_Integer, gates_List] := With[{full = Join[PathSum`hLayer[n], gates]},
    AllTrue[Tuples[{0, 1}, n],
        PathSum`exactEqualQ[PathSum`diagAmp[n, gates, #], PathSum`qfAmp[n, full, #]] &]];

PathSum`diagPolyMatchesQ[n_Integer, gates_List] := With[{c = PathSum`diagonalCompress[n, gates]},
    AllTrue[Tuples[{0, 1}, n],
        PathSum`exactEqualQ[PathSum`diagAmpFromPoly[c, #], PathSum`diagAmp[n, gates, #]] &]];

(* random diagonal circuit (no H): {T,Sdg,S,Z,CZ,CNOT,CCZ}. *)
PathSum`randomDiagGate[n_Integer] := Module[{kinds = Join[
    {"T", SuperDagger["T"], "S", "Z"}, If[n >= 2, {"CZ", "CNOT"}, {}], If[n >= 3, {"CCZ"}, {}]], k},
    k = RandomChoice[kinds];
    Switch[k, "CZ" | "CNOT", {k, Sequence @@ RandomSample[Range[n], 2]},
              "CCZ", {k, Sequence @@ RandomSample[Range[n], 3]},
              _, {k, RandomInteger[{1, n}]}]];
PathSum`randomDiagCircuit[n_Integer, depth_Integer] := Table[PathSum`randomDiagGate[n], {depth}];

(* ============================================================================== *)
(* Isolated single-variable Gauss-sum elimination lemmas (Stage 3 building blocks). *)
(*                                                                                 *)
(* For a phase phi linear in the idempotent variable v, write phi = phi_0 + v C,    *)
(* C = Coefficient[phi, v] the cofactor. Then                                       *)
(*     Sum_{v in {0,1}} omega^{phi} = omega^{phi_0} (1 + omega^C).                  *)
(* The cofactor CLASS (mod 8) decides what happens:                                 *)
(*   C = 4 Q (Q Boolean):  1 + (-1)^Q = 2 [Q = 0]    -- the Z_2 [HH] rule (a sign,   *)
(*                                                      eliminates to a constraint), *)
(*   C = 2 Q (Q Boolean):  1 + i^Q                   -- the Z_4 / S case: a genuine  *)
(*                                                      Gauss factor, NOT a          *)
(*                                                      constraint (the A mod 4 != 0  *)
(*                                                      branch elimQ rejects),        *)
(*   C odd (a T phase):    1 + omega^C               -- magic; the variable survives.*)
(* ============================================================================== *)

(* General single-variable reduction (closed form, exact, all cofactor classes). *)
PathSum`reduceVar1[phase_, v_] := Module[
    {phi0 = PathSum`z8[phase /. v -> 0], C = PathSum`cofactor[phase, v]},
    PathSum`omega^phi0 (1 + PathSum`omega^C)];

(* Brute-force single-variable sum, the oracle for the lemma. *)
PathSum`bruteVar1[phase_, v_] := PathSum`omega^PathSum`z8[phase /. v -> 0] + PathSum`omega^PathSum`z8[phase /. v -> 1];

(* Cofactor class label, for the unit tests. *)
PathSum`cofactorClass[phase_, v_] := With[{C = PathSum`cofactor[phase, v]},
    Which[PolynomialMod[C, 8] === 0, "free",
          PolynomialMod[C, 4] === 0, "Z2-sign",
          PolynomialMod[C, 2] === 0, "Z4-S",
          True, "odd-T"]];
