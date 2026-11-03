{a, ad} = FieldVariables[{1}];
{SqueezeOperator[xi, Infinity, "Ordering" -> "Normal"], FullSimplify[BosonicExpOrder[Exp[(r Exp[-I th] a ** a - r Exp[I th] ad ** ad)/2]], r > 0 && Element[th, Reals]]}
