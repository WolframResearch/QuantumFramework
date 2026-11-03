{a, ad} = FieldVariables[{1}];
{BosonicNormalOrder[ad ** a ** a - a ** ad ** a, {a, ad}], BosonicNormalOrder[ad ** a ** ad - ad ** ad ** a, {a, ad}],
 FullSimplify[BosonicExpSimplify[Exp[I wt ad ** a] ** a ** Exp[-I wt ad ** a]], wt > 0],
 FullSimplify[Exp[-I wt/2] BosonicMatrixElement[{be, al}, Exp[-I wt ad ** a], "Basis" -> "Coherent"] == Exp[-I wt/2] BosonicMatrixElement[{be, al Exp[-I wt]}, 1, "Basis" -> "Coherent"], wt > 0]}
