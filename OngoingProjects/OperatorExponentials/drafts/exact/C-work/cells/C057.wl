{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
af[f_] := Conjugate[f[[1]]] a1 + Conjugate[f[[2]]] a2;
adf[f_] := f[[1]] a1d + f[[2]] a2d;
{BosonicNormalOrder[af[{f1, f2}] ** adf[{g1, g2}] - adf[{g1, g2}] ** af[{f1, f2}], {a1, a1d, a2, a2d}, "Scalars" -> {f1, f2, g1, g2}],
 Conjugate[{f1, f2}] . {g1, g2}}
