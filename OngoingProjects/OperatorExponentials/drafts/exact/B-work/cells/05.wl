{a, ad} = FieldVariables[{1}];
{BosonicMatrixElement[{m, m}, NormalOrdered[Exp[-ad ** a]]],
 Table[BosonicMatrixElement[{j, k}, NormalOrdered[Exp[-ad ** a]]], {j, 0, 3}, {k, 0, 3}]}
