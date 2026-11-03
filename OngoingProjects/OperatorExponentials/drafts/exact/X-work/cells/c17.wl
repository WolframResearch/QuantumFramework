{a, ad} = FieldVariables[{1}];
{BosonicMatrixElement[{m, n}, NormalOrdered[Exp[-ad ** a]]], BosonicMatrixElement[{z, w}, NormalOrdered[Exp[-ad ** a]], "Basis" -> "Coherent"], NormalOrdered[Exp[-ad ** a]]["Series"]}
