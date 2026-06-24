(* ::Package:: *)

(* pathsum.wl : symbolic path-sum simulator for Clifford+T (standalone).        *)
(*                                                                              *)
(* A Clifford+T amplitude is carried as a Gauss sum                            *)
(*     <y|C|0...0> = scale / Sqrt[2]^h  Sum_{v : E(v)=y}  omega^{phi(v)},        *)
(* with omega = e^{i pi/4} a primitive 8th root of unity, each qubit value an   *)
(* F_2 multilinear polynomial in Hadamard branch variables x[i], and the phase  *)
(* phi a degree-<=3 pseudo-Boolean polynomial valued in Z_8.                     *)
(*                                                                              *)
(* CONVENTION (resolved, Stage 0a): the full-angle phase gate                   *)
(*     "P"[theta] = diag(1, e^{i theta}),   T = "P"[pi/4] = diag(1, e^{i pi/4}). *)
(* The T rule "phase += e_q" gives omega^{e_q} = e^{i pi/4 e_q}, i.e. the        *)
(* canonical full-angle T that QuantumCircuitOperator Method->"Schrodinger"      *)
(* uses. (StabilizerFrame's internal "P" is half-angle; this simulator does NOT  *)
(* use it.) All amplitudes are exact in Z[omega, 1/Sqrt2]; comparisons use       *)
(* RootReduce[a - b] === 0, never Chop[N[...]].                                  *)
(*                                                                              *)
(* Assumes Wolfram`QuantumFramework` is already loaded by the caller.           *)

PathSum`omega = Exp[I Pi/4];

(* ---- F_2 and integer-pseudo-Boolean normal forms (idempotent variables) ---- *)
(* idempotency of 0/1 variables: x[i]^k -> x[i]. This is the native reduction    *)
(* for the Z_8 phase, where GroebnerBasis at the non-prime Modulus->8 is         *)
(* unreliable; the prime-modulus (Modulus->2) Boolean-ideal machinery enters     *)
(* only for the F_2 constraint/variety work in later stages.                     *)
PathSum`idem[p_] := Expand[p] /. Power[x[i_], _] :> x[i];
PathSum`f2[p_]   := PolynomialMod[PathSum`idem[p], 2];          (* F_2 coefficients *)
PathSum`z8[p_]   := PolynomialMod[PathSum`idem[p], 8];          (* Z_8 coefficients *)

(* Boolean (XOR) polynomial -> integer-valued pseudo-Boolean (a (+) b = a+b-2ab). *)
PathSum`boolToInt[p_] := PathSum`idem @ Fold[#1 + #2 - 2 #1 #2 &, 0, MonomialList[PathSum`f2[p]]];

(* ---- path-sum state ---- *)
(* <|"val" -> <|q -> F_2 poly|>, "phase" -> Z_8 poly, "vars" -> {x[..]},          *)
(*   "h" -> #Hadamards, "scale" -> exact 2^(#eliminations)|>                      *)
PathSum`init[n_Integer] := <|
    "val" -> Association[Table[q -> 0, {q, n}]], "phase" -> 0,
    "vars" -> {}, "h" -> 0, "scale" -> 1
|>;

PathSum`apply[ps_, g_] := Module[
    {val = ps["val"], phase = ps["phase"], vars = ps["vars"], h = ps["h"], v},
    Switch[g,
        {"H", _},          v = x[Length[vars] + 1]; AppendTo[vars, v]; h++;
                           phase += 4 PathSum`boolToInt[val[g[[2]]]] v; val[g[[2]]] = v,
        {"T", _},          phase += PathSum`boolToInt[val[g[[2]]]],          (* diag(1, omega) *)
        {SuperDagger["T"], _}, phase += 7 PathSum`boolToInt[val[g[[2]]]],    (* diag(1, omega^-1) *)
        {"S", _},          phase += 2 PathSum`boolToInt[val[g[[2]]]],        (* diag(1, i) *)
        {"Z", _},          phase += 4 PathSum`boolToInt[val[g[[2]]]],        (* diag(1, -1) *)
        {"CZ", _, _},      phase += 4 PathSum`boolToInt[val[g[[2]]]] PathSum`boolToInt[val[g[[3]]]],
        {"CCZ", _, _, _},  phase += 4 PathSum`boolToInt[val[g[[2]]]] PathSum`boolToInt[val[g[[3]]]] PathSum`boolToInt[val[g[[4]]]],
        {"CNOT", _, _},    val[g[[3]]] = PathSum`f2[val[g[[2]]] + val[g[[3]]]]
    ];
    <|ps, "val" -> val, "phase" -> PathSum`idem[phase], "vars" -> vars, "h" -> h|>
];

PathSum`build[n_Integer, gates_List] := Fold[PathSum`apply, PathSum`init[n], gates];

(* ---- exact amplitude: brute-force Gauss sum over the remaining variables ---- *)
(* This is the correctness ORACLE and the permanent fallback for the rewriter.   *)
PathSum`amp[ps_, n_Integer, y_List] := With[{vars = ps["vars"]},
    (ps["scale"] / Sqrt[2]^ps["h"]) Total[
        (PathSum`omega^Mod[ps["phase"] /. Thread[vars -> #], 8]) Boole[
            And @@ Table[Mod[ps["val"][q] /. Thread[vars -> #], 2] == y[[q]], {q, n}]] & /@
        Tuples[{0, 1}, Length[vars]]]
];

(* ============================================================================ *)
(* [HH] Hadamard-elimination (Z_2 sign-cofactor case). Native cofactor via        *)
(* Coefficient: for the multilinear phase, C_v = Coefficient[phi, x_v] and        *)
(* phi_0 = Coefficient[phi, x_v, 0]. When C_v = 4 Q with Q a linear Boolean form, *)
(*     Sum_v omega^{v C_v + phi_0} = 2 omega^{phi_0} [Q = 0],                      *)
(* so v is summed away and Q=0 imposed, with no enumeration.                      *)
(* ============================================================================ *)

PathSum`cofactor[phase_, v_] := PolynomialMod[PathSum`idem[Coefficient[PathSum`idem[phase], v]], 8];

PathSum`elimQ[phase_, v_] := Module[{A = PathSum`cofactor[phase, v], Q},
    If[PolynomialMod[A, 4] =!= 0, Return[False]];     (* cofactor must be a pure sign 4Q *)
    Q = PathSum`f2[A/4];
    Length[Variables[Q]] >= 1 && AllTrue[MonomialList[Q], Length[Variables[#]] <= 1 &]
];

PathSum`internalQ[ps_, v_] := AllTrue[Values[ps["val"]], FreeQ[#, v] &];

PathSum`elimStep[ps_] := Module[{phase = ps["phase"], v, A, Q, u, rhs, val},
    v = SelectFirst[ps["vars"], PathSum`elimQ[phase, #] && PathSum`internalQ[ps, #] &, Missing[]];
    If[MissingQ[v], Return[ps]];
    A = PathSum`cofactor[phase, v]; Q = PathSum`f2[A/4];
    u = First[Variables[Q]]; rhs = PathSum`f2[u + Q];            (* impose Q = 0: u = rhs (F_2) *)
    (* The phase treats u as an INTEGER 0/1 value, so substitute the integer value      *)
    (* boolToInt[rhs]; the F_2 form is correct only when all phase coefficients are 0    *)
    (* mod 4 (pure Z_2). val is an F_2 polynomial and takes the F_2 form directly.       *)
    phase = PathSum`z8[(phase /. v -> 0) /. (u -> PathSum`boolToInt[rhs])];
    val = PathSum`f2[#] & /@ (ps["val"] /. (u -> rhs));
    <|ps, "phase" -> phase, "val" -> val,
      "vars" -> DeleteCases[ps["vars"], v | u], "scale" -> ps["scale"] * 2|>
];

PathSum`reduce[ps_] := FixedPoint[PathSum`elimStep, ps];
