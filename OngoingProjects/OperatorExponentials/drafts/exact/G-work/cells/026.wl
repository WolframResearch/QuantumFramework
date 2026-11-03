{a, ad} = FieldVariables[{1}];
ov[x_, y_] := BosonicMatrixElement[{x, y}, 1, "Basis" -> "Coherent"];
f = ov[be, al] ov[al, ga] /. Abs[al] -> Sqrt[al Conjugate[al]] /. Conjugate[al] -> u /. al -> v;
{Simplify[u f + D[f, v] - Conjugate[be] f], Simplify[v f + D[f, u] - ga f]}
