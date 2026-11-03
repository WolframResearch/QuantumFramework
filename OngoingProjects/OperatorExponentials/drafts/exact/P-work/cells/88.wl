{a, ad} = FieldVariables[{1}];
{BosonicMatrixElement[{al, be}, Exp[lam ad ** a], "Basis" -> "Coherent"], BosonicExpOrder[Exp[lam ad ** a]],
 Sum[Exp[-w k/T], {k, 0, Infinity}, Assumptions -> w/T > 0]}
