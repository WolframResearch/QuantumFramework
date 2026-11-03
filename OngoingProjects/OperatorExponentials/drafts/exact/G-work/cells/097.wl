{a, ad} = FieldVariables[{1}];
ov[x_, y_] := BosonicMatrixElement[{x, y}, 1, "Basis" -> "Coherent"];
k[x_] := FullSimplify[BosonicMatrixElement[{u, x}, Exp[Log[tT]/2 ad ** a], "Basis" -> "Coherent"]/ov[u, Sqrt[tT] x], 0 < tT < 1];
{k[al], FullSimplify[k[al] Conjugate[k[be]] Sum[(1 - tT)^l/l! (al Conjugate[be])^l, {l, 0, Infinity}] - Exp[-(1 - tT)/2 (Abs[al]^2 + Abs[be]^2 - 2 Conjugate[be] al)], 0 < tT < 1]}
