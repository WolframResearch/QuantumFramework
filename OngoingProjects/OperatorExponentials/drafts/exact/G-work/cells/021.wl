{a, ad} = FieldVariables[{1}];
Simplify[BosonicMatrixElement[{be, al}, ad ** ad ** a, "Basis" -> "Coherent"]/BosonicMatrixElement[{be, al}, 1, "Basis" -> "Coherent"]]
