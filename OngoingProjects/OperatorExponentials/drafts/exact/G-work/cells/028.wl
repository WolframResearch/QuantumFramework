{a, ad} = FieldVariables[{1}];
ov[x_, y_] := BosonicMatrixElement[{x, y}, 1, "Basis" -> "Coherent"];
f = ov[be, al] ov[al, ga] /. Abs[al] -> Sqrt[al Conjugate[al]] /. Conjugate[al] -> u /. al -> v;
Table[Simplify[(D[Exp[u v] f, {v, m}, {u, m}] /. {u -> 0, v -> 0})/m! - Conjugate[BosonicMatrixElement[{m, 0}, DisplacementOperator[be]]] BosonicMatrixElement[{m, 0}, DisplacementOperator[ga]]], {m, 0, 6}]
