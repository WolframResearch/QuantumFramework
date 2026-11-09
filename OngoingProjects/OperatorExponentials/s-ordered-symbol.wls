PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName, 3], "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
(* Idea 3 prototype, kept outside QF: the s-ordered symbol of an operator is exp(-(1-s)/2 Sum_i d_Ai d_Bi) applied to its normal symbol (A_i for creators, B_i for annihilators), done exactly on polynomial times Gaussian terms by adding a linear term j to the exponent (differentiate in j, then set j = 0) and the Gaussian heat formula. Run: wolframscript -file s-ordered-symbol.wls *)
out[s_, e_] := WriteString["stdout", s, " -> ", ToString[e, InputForm], "\n"];
pw[x_, k_] := GeneralizedPower[NonCommutativeMultiply, x, k];
deq[e_] := e //. {GeneralizedPower[NonCommutativeMultiply, b_, p_] :> b^p, NonCommutativeMultiply -> Times};

(* normal symbol: creators -> A_i, annihilators -> B_i, after normal ordering *)
normalSymbol[expr_, {cs_, as_}, {Asym_, Bsym_}] := Module[{ns = Thread[cs ** as], ord, numFix},
  ord = If[FreeQ[expr, E | Exp | NormalOrdered], BosonicNormalOrder[expr, {cs, as}], BosonicExpOrder[expr, {cs, as}]];
  (* a function of N standing in the middle of an ordered product is normal-ordered as :e^{(e^c - 1) N}: *)
  numFix = {
    NormalOrdered[x_, _] :> (x /. ns -> Thread[Asym Bsym] /. Thread[cs -> Asym] /. Thread[as -> Bsym]),
    Power[b_, e_] /; ! FreeQ[e, Alternatives @@ ns] && FreeQ[b, Alternatives @@ Join[cs, as]] :>
      With[{cf = Coefficient[e, ns]}, b^(e - cf . ns) Exp[(b^cf - 1) . (Asym Bsym)]],
    Exp[e_] /; ! FreeQ[e, Alternatives @@ ns] && FreeQ[e - Coefficient[e, ns] . ns, Alternatives @@ Join[cs, as]] :>
      With[{cf = Coefficient[e, ns]}, Exp[e - cf . ns] Exp[(Exp[cf] - 1) . (Asym Bsym)]]};
  (* an exponential mixing creators and annihilators was left unordered: no normal symbol *)
  If[! FreeQ[ord, Power[E, e_] /; ! FreeQ[e, Alternatives @@ cs] && ! FreeQ[e, Alternatives @@ as] && FreeQ[e, Alternatives @@ ns]], Return[$Failed, Module]];
  deq[ord /. numFix] /. Join[Thread[cs -> Asym], Thread[as -> Bsym]]
];

(* exp(lam Sum_i d_Ai d_Bi) on p(z) Exp[Q(z)], Q at most quadratic: linear term j in the exponent plus the Gaussian heat formula *)
heatTerm[t_, z_, Lm_] := Module[{p, ex, Q, cA, j0, M, src, js, k = Length[z], Id, R},
  {p, ex} = Replace[t, {Times[pp___, Power[E, q_]] :> {Times[pp], q}, Power[E, q_] :> {1, q}, other_ :> {other, 0}}];
  ex = Expand[ex];
  If[! PolynomialQ[ex, z] || Exponent[ex, z, Max] > 2 || ! PolynomialQ[p, z], Return[$Failed, Module]];
  cA = CoefficientArrays[ex, z];
  {Q, j0, M} = {cA[[1]], If[Length[cA] > 1, Normal[cA[[2]]], ConstantArray[0, k]], If[Length[cA] > 2, Normal[cA[[3]] + Transpose[cA[[3]]]], ConstantArray[0, {k, k}]]};
  js = Array[src, k]; Id = IdentityMatrix[k]; R = Inverse[Id - Lm . M];
  With[{g = Det[Id - Lm . M]^(-1/2) Exp[Q + 1/2 z . M . R . z + (j0 + js) . R . z + 1/2 (j0 + js) . R . Lm . (j0 + js)]},
    (* p(z) -> p(d/dj) acting before j -> 0 *)
    Expand[Total[KeyValueMap[#2 D[g, Sequence @@ Transpose[{js, #1}]] &, Association[CoefficientRules[p, z]]] /. Thread[js -> 0]]]]
];
sSymbol[expr_, {cs_, as_}, s_] := Module[{k = Length[cs], Asym, Bsym, z, Lm, F},
  Asym = Array[AA, k]; Bsym = Array[BB, k]; z = Join[Asym, Bsym];
  Lm = ArrayFlatten[{{0, -(1 - s)/2 IdentityMatrix[k]}, {-(1 - s)/2 IdentityMatrix[k], 0}}];
  F = normalSymbol[expr, {cs, as}, {Asym, Bsym}];
  If[F === $Failed, Return[$Failed, Module]];
  F = Expand[F, z];
  Total[heatTerm[#, z, Lm] & /@ If[Head[F] === Plus, List @@ F, {F}]] /. {AA[i_] :> Conjugate[al[i]], BB[i_] :> al[i]}
];
v = {{ad}, {a}};
(* 1. displacement: s-symbol = Exp[-s |g|^2/2 + g abar - gbar al] *)
D1 = Exp[g ad - gc a];
out["D(g), all s", Simplify[sSymbol[D1, v, s] - Exp[-s g gc/2 + g Conjugate[al[1]] - gc al[1]]]];
(* 2. s = 1 is the normal symbol, s = -1 the antinormal symbol, on a polynomial *)
P = a ** a ** ad ** ad ** a + x ad ** a ** ad;
out["s=-1 vs antinormal", Expand[sSymbol[P, v, -1] - (deq[BosonicAntinormalOrder[P, v]] /. {ad -> Conjugate[al[1]], a -> al[1]})]];
out["s=1 vs normal", Expand[sSymbol[P, v, 1] - (deq[BosonicNormalOrder[P, v]] /. {ad -> Conjugate[al[1]], a -> al[1]})]];
(* 3. C-053: Weyl symbol of e^{-tau (N + 1/2)} = sech(tau/2) e^{-2 tanh(tau/2) |al|^2} *)
out["C-053", Simplify[sSymbol[Exp[-tau (ad ** a + 1/2)], v, 0] - Sech[tau/2] Exp[-2 Tanh[tau/2] Conjugate[al[1]] al[1]], tau > 0]];
(* 4. Fock projector |n><n| = ad^n :e^{-N}: a^n / n!, Weyl symbol 2 (-1)^n e^{-2|al|^2} L_n(4|al|^2) *)
out["|n><n| Weyl, n=0..4", Table[Simplify[sSymbol[pw[ad, n] ** NormalOrdered[Exp[-ad ** a], v] ** pw[a, n]/n!, v, 0] - 2 (-1)^n Exp[-2 Conjugate[al[1]] al[1]] LaguerreL[n, 4 Conjugate[al[1]] al[1]]], {n, 0, 4}]];
(* 5. C-016: Weyl symbol of [H, B] = Poisson bracket of Weyl symbols for deg H <= 2; fails for quartic H *)
pb[f_, g_] := With[{A = Conjugate[al[1]], B = al[1]}, D[f, B] D[g, A] - D[f, A] D[g, B]];
pbcheck[H_, Bop_] := Expand[sSymbol[BosonicNormalOrder[H ** Bop - Bop ** H, v], v, 0] - (pb[sSymbol[H, v, 0], sSymbol[Bop, v, 0]] /. Conjugate[al[1]] -> AAc) /. AAc -> Conjugate[al[1]]];
out["C-016 quadratic H", pbcheck[ad ** a + k (a ** a + ad ** ad) + mu ad ** ad, ad ** ad ** a + nu a ** a ** a]];
out["C-016 quartic H must fail", pbcheck[ad ** ad ** a ** a, ad ** ad ** a + nu a ** a ** a] =!= 0];
(* 6. two modes: two-mode squeezer vacuum projector check via s=1 *)
v2 = {{ad1, ad2}, {a1, a2}};
out["two-mode D as one exponential (multimode rule open)", sSymbol[Exp[g ad1 - gc a1 + h ad2 - hc a2], v2, 0]];
out["two-mode D1 ** D2 Weyl", Simplify[sSymbol[Exp[g ad1 - gc a1] ** Exp[h ad2 - hc a2], v2, 0] - Exp[g Conjugate[al[1]] - gc al[1] + h Conjugate[al[2]] - hc al[2]]]];
