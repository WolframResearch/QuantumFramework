#!/usr/bin/env wolframscript
(* Verification demos for PathSum-Blowup-Alternatives-And-WL-Contribution.md.       *)
(* Backs Question 2: the single most original WL contribution = the exact          *)
(* symbolic-PARAMETER path sum (free gate angles -> closed-form amplitude /         *)
(* expectation / gradient FUNCTIONS over a parametrized circuit family in one       *)
(* pass), plus the bounded companion (b) Groebner/variety point count over F_2.     *)
(*                                                                                  *)
(* Every amplitude is checked EXACTLY against QF dense Method->"Schrodinger":        *)
(*   - at special angles via RootReduce in Z[omega,1/Sqrt2] (no floating point),     *)
(*   - at generic real angles numerically to < 10^-12 (a sanity floor only).         *)
(* No Quiet, no Print; each block ends with a labeled Echo of the asserted result.  *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

(* ============================================================================ *)
(* A standalone SYMBOLIC-ANGLE path-sum builder. Unlike pathsum.wl (Z_8 integer  *)
(* phases), this carries the amplitude phase as a symbolic real polynomial Phi    *)
(* with the gate factor E^{I Phi}; a free angle theta enters as E^{I theta e_q}.   *)
(*   H on q : fresh v, factor 1/Sqrt2, sign (-1)^{e_q v} = E^{I Pi e_q v}, e_q:=v.  *)
(*   P[ang] : E^{I ang e_q}  (full-angle diag(1,e^{i ang}), the QCO convention).    *)
(*   T : E^{I Pi/4 e_q}; S : E^{I Pi/2 e_q}; Z : E^{I Pi e_q}.                       *)
(*   CZ : E^{I Pi e_a e_b}; CCZ : E^{I Pi e_a e_b e_c}; CNOT : e_b := e_a (+) e_b.   *)
(* e_q is the INTEGER 0/1 value boolToInt of the F_2 wire polynomial.               *)
(* ============================================================================ *)
idem[p_] := Expand[p] /. Power[xx[i_], _] :> xx[i];
f2[p_]   := PolynomialMod[idem[p], 2];
boolToInt[p_] := idem @ Fold[#1 + #2 - 2 #1 #2 &, 0, MonomialList[f2[p]]];

initPS[n_] := <|"val" -> Association[Table[q -> 0, {q, n}]], "phase" -> 0, "vars" -> {}, "h" -> 0|>;
applyG[ps_, g_] := Module[{val = ps["val"], phase = ps["phase"], vars = ps["vars"], h = ps["h"], v},
   Switch[g,
     {"H", _},       v = xx[Length[vars] + 1]; AppendTo[vars, v]; h++;
                     phase += Pi boolToInt[val[g[[2]]]] v; val[g[[2]]] = v,
     {"P", _, _},    phase += g[[3]] boolToInt[val[g[[2]]]],     (* free or numeric angle *)
     {"T", _},       phase += (Pi/4) boolToInt[val[g[[2]]]],
     {"S", _},       phase += (Pi/2) boolToInt[val[g[[2]]]],
     {"Z", _},       phase += Pi boolToInt[val[g[[2]]]],
     {"CZ", _, _},   phase += Pi boolToInt[val[g[[2]]]] boolToInt[val[g[[3]]]],
     {"CCZ", _, _, _}, phase += Pi boolToInt[val[g[[2]]]] boolToInt[val[g[[3]]]] boolToInt[val[g[[4]]]],
     {"CNOT", _, _}, val[g[[3]]] = f2[val[g[[2]]] + val[g[[3]]]]];
   <|ps, "val" -> val, "phase" -> idem[phase], "vars" -> vars, "h" -> h|>];
buildPS[n_, gates_] := Fold[applyG, initPS[n], gates];

(* Closed-form amplitude <y|C|0...0>(theta) : a finite sum of E^{I Phi}, symbolic. *)
ampSym[ps_, n_, y_] := With[{vars = ps["vars"]},
   (1/Sqrt[2]^ps["h"]) Total[
     (Exp[I (ps["phase"] /. Thread[vars -> #])] Boole[
        And @@ Table[Mod[ps["val"][q] /. Thread[vars -> #], 2] == y[[q]], {q, n}]]) & /@
     Tuples[{0, 1}, Length[vars]]]];
ampVecSym[n_, gates_] := Table[ampSym[buildPS[n, gates], n, IntegerDigits[b, 2, n]], {b, 0, 2^n - 1}];

(* QF dense oracle, with free or numeric angles substituted in. *)
qfGateMap = {{"H", q_} :> "H" -> q, {"P", q_, a_} :> "P"[a] -> q, {"T", q_} :> "T" -> q,
   {"S", q_} :> "S" -> q, {"Z", q_} :> "Z" -> q, {"CZ", a_, b_} :> "CZ" -> {a, b},
   {"CCZ", a_, b_, c_} :> QuantumOperator[DiagonalMatrix[{1,1,1,1,1,1,1,-1}], {a, b, c}],
   {"CNOT", a_, b_} :> "CNOT" -> {a, b}};
qfVec[n_, gates_] := Normal[
   QuantumCircuitOperator[gates /. qfGateMap][QuantumState[StringRepeat["0", n]], Method -> "Schrodinger"]["StateVector"]];

exactEqualQ[a_, b_] := PossibleZeroQ[RootReduce[a - b]];
numEqual[a_, b_] := Max[Abs[N[a - b, 30]]] < 10^-12;

(* ---------------------------------------------------------------------------- *)
Echo["==== (b1) ONE free angle: H P(theta) H |0>, closed-form amplitude+expectation ====", ""];
g1[t_] := {{"H", 1}, {"P", 1, t}, {"H", 1}};
amp1 = FullSimplify[ampVecSym[1, g1[\[Theta]]]];
Echo[amp1, "(b1) <x|H P(theta) H|0> closed form (expect {(1+e^{i th})/2, (1-e^{i th})/2}):"];
Echo[FullSimplify[amp1 - {(1 + Exp[I \[Theta]])/2, (1 - Exp[I \[Theta]])/2}] === {0, 0},
   "(b1) matches the StabilizerFrame component coefficients? (expect True):"];
zExpr = FullSimplify[Conjugate[amp1] . ({{1,0},{0,-1}} . amp1), Assumptions -> \[Theta] \[Element] Reals];
Echo[zExpr, "(b1) <Z>(theta) closed form, ONE symbolic pass (expect Cos[theta]):"];
Echo[FullSimplify[zExpr - Cos[\[Theta]]] === 0, "(b1) <Z>(theta) == Cos[theta]? (expect True):"];
(* exact verification against dense at SPECIAL angles (no floating point) *)
exactAngles = {0, Pi/4, Pi/3, Pi/2, 2 Pi/3, Pi, Pi/5, Pi/7};
Echo[AllTrue[exactAngles, AllTrue[Thread[{amp1 /. \[Theta] -> #, qfVec[1, g1[#]]}],
     Apply[exactEqualQ]] &],
   "(b1) closed form === QF dense EXACTLY (RootReduce) at 8 special angles? (expect True):"];

(* ---------------------------------------------------------------------------- *)
Echo["==== (b2) THE capability: closed-form GRADIENT, free vs a re-run dense grid ====", ""];
dZ = D[zExpr, \[Theta]];   (* exact analytic derivative from the one closed form *)
Echo[dZ, "(b2) d/dtheta <Z>(theta) from the closed form (expect -Sin[theta]):"];
(* compare to a central finite difference of DENSE simulation at sampled angles *)
denseZ[t_] := Re[Conjugate[#] . ({{1,0},{0,-1}} . #)] &[qfVec[1, g1[t]]];
fd[t_, h_] := (denseZ[t + h] - denseZ[t - h]) / (2 h);
gradCheck = Table[Abs[N[dZ /. \[Theta] -> t] - fd[t, 10^-6]], {t, {0.3, 0.9, 1.7, 2.5}}];
Echo[Max[gradCheck] < 10^-5,
   "(b2) exact gradient matches dense finite-difference at 4 angles? (expect True):"];

(* ---------------------------------------------------------------------------- *)
Echo["==== (b3) TWO free angles, 2 qubits: <Z1 Z2>(a,b) in ONE pass ====", ""];
(* hardware-efficient-style layer: (H H) . P(a)_1 P(b)_2 . CZ . (H H) on |00>     *)
g2[a_, b_] := {{"H", 1}, {"H", 2}, {"P", 1, a}, {"P", 2, b}, {"CZ", 1, 2}, {"H", 1}, {"H", 2}};
amp2 = ampVecSym[2, g2[\[Alpha], \[Beta]]];
Z1Z2 = DiagonalMatrix[{1, -1, -1, 1}];
e2 = FullSimplify[Conjugate[amp2] . (Z1Z2 . amp2), Assumptions -> {\[Alpha] \[Element] Reals, \[Beta] \[Element] Reals}];
Echo[e2, "(b3) <Z1 Z2>(alpha,beta) closed form, one symbolic pass:"];
(* exact check against dense on a grid of special angles *)
grid = Tuples[{0, Pi/4, Pi/3, Pi/2, Pi}, 2];
Echo[AllTrue[grid, AllTrue[Thread[{amp2 /. {\[Alpha] -> #[[1]], \[Beta] -> #[[2]]}, qfVec[2, g2[#[[1]], #[[2]]]]}],
     Apply[exactEqualQ]] &],
   "(b3) amplitude(alpha,beta) === QF dense EXACTLY on a 25-point special-angle grid? (expect True):"];
(* a representative numeric expectation cross-check at a generic point *)
denseZZ[a_, b_] := Re[Conjugate[#] . (Z1Z2 . #)] &[qfVec[2, g2[a, b]]];
Echo[numEqual[e2 /. {\[Alpha] -> 0.7, \[Beta] -> 1.3}, denseZZ[0.7, 1.3]],
   "(b3) closed-form <Z1 Z2> == dense at a generic (0.7,1.3)? (expect True):"];

(* ---------------------------------------------------------------------------- *)
Echo["==== (b4) free angle COMPOSED WITH fixed magic (T): closed form survives ====", ""];
(* H . T . P(theta) . H |0> : a fixed T and a free P(theta) in the same pass.      *)
g3[t_] := {{"H", 1}, {"T", 1}, {"P", 1, t}, {"H", 1}};
amp3 = FullSimplify[ampVecSym[1, g3[\[Theta]]]];
Echo[amp3, "(b4) <x|H T P(theta) H|0> closed form (carries e^{i Pi/4} AND theta):"];
Echo[AllTrue[{0, Pi/5, Pi/2, Pi}, AllTrue[Thread[{amp3 /. \[Theta] -> #, qfVec[1, g3[#]]}], Apply[exactEqualQ]] &],
   "(b4) === QF dense EXACTLY at 4 special angles? (expect True):"];

(* ============================================================================ *)
(* COMPANION (a): the Groebner / variety point count over F_2 (Dawson 2004 idea, *)
(* Amy-Stinchcombe future-work technique). gap(f) = #{f=0}-#{f=1} = 2#{f=0}-2^n. *)
(* ============================================================================ *)
Echo["==== (a) Groebner/variety point count over F_2 (the bounded companion) ====", ""];
fieldEqs[vars_] := (#^2 - # &) /@ vars;
pointCountGB[polys_, vars_] := Module[{gb, mons},
   gb = GroebnerBasis[Join[polys, fieldEqs[vars]], vars, Modulus -> 2];
   If[gb === {1} || MemberQ[gb, 1], 0,
    mons = Times @@@ Subsets[vars];
    Count[mons, m_ /; PolynomialMod[Last[PolynomialReduce[m, gb, vars, Modulus -> 2]] - m, 2] === 0]]];
pointCountBrute[polys_, vars_] := Count[Tuples[{0, 1}, Length[vars]],
   a_ /; And @@ (PolynomialMod[# /. Thread[vars -> a], 2] === 0 & /@ polys)];
vars3 = {x1, x2, x3}; vars4 = {x1, x2, x3, x4};
Echo[GroebnerBasis[{x1 + x2, x1, x2 + 1}, {x1, x2}, Modulus -> 2] === {1},
   "(a) Nullstellensatz: inconsistent system GB={1} => amplitude 0? (expect True):"];
Echo[{pointCountGB[{x1 x2 + x3}, vars3], pointCountBrute[{x1 x2 + x3}, vars3],
      2 pointCountBrute[{x1 x2 + x3}, vars3] - 2^3},
   "(a) quadratic term {GB count, brute count, gap=2#-2^n}:"];
Echo[pointCountGB[{x1 x2 + x3}, vars3] === pointCountBrute[{x1 x2 + x3}, vars3],
   "(a) quadratic GB quotient-dim == brute? (expect True):"];
Echo[pointCountGB[{x1 x2 x3 + x4}, vars4] === pointCountBrute[{x1 x2 x3 + x4}, vars4],
   "(a) cubic (CCZ-type) GB quotient-dim == brute? (expect True):"];

(* ============================================================================ *)
(* BOUNDARY: the characteristic-2 wall that blocks the GF(2^k) trace route and    *)
(* forbids a Z_8 phase from enlarging the tractable class (no #P break).          *)
(* ============================================================================ *)
Echo["==== boundary: the characteristic-2 wall (no tractable-class enlargement) ====", ""];
Echo[FullSimplify[(Exp[I Pi/4] PauliMatrix[1]) . (Exp[I Pi/4] PauliMatrix[1])] === I IdentityMatrix[2],
   "(wall) (e^{i pi/4} X)^2 = i I != I (no +1 eigenspace, stabilizes nothing)? (expect True):"];
Echo[Factor[Cyclotomic[8, x], Modulus -> 2],
   "(wall) 1+x^4 over F_2 (expect (1+x)^4: omega's cyclotomic structure gone mod 2):"];

Echo["ALL DEMOS COMPLETE", "status"];
