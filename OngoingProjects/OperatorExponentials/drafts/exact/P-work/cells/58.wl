{a, ad} = FieldVariables[{1}];
Simplify[BosonicMatrixElement[{al, be}, ad ** a, "Basis" -> "Coherent"]/BosonicMatrixElement[{al, be}, 1, "Basis" -> "Coherent"]]
