{a, ad} = FieldVariables[{1}];
{FullSimplify[BosonicExpSimplify[Exp[I K t m ad ** a] ** a ** Exp[-I K t m ad ** a]], K > 0 && t > 0 && m > 0],
 FullSimplify[BosonicMatrixElement[{be, al}, Exp[-I K t m ad ** a], "Basis" -> "Coherent"] == BosonicMatrixElement[{be, al Exp[-I K t m]}, 1, "Basis" -> "Coherent"], Element[{K, t, m}, Reals]]}
