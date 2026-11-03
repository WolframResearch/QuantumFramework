{a, ad} = FieldVariables[{1}];
{FullSimplify[BosonicMatrixElement[{be, al}, Exp[-s ad ** a], "Basis" -> "Coherent"] == Exp[-Abs[al]^2 (1 - Exp[-2 s])/2] BosonicMatrixElement[{be, al Exp[-s]}, 1, "Basis" -> "Coherent"], s > 0],
 FullSimplify[BosonicMatrixElement[{al, al}, Exp[-2 s ad ** a], "Basis" -> "Coherent"], s > 0]}
