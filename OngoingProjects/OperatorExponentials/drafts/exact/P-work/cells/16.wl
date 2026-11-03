{a, ad} = FieldVariables[{1}];
BosonicMatrixElement[{al, al}, ad ** ad ** a, "Basis" -> "Coherent"] /. Abs[al]^2 -> al Conjugate[al]
