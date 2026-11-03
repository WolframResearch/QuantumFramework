{a, ad} = FieldVariables[{1}];
FullSimplify[BosonicExpOrder[Exp[r Exp[I ph] ad ** ad/2 - r Exp[-I ph] a ** a/2]], r > 0 && Element[ph, Reals]]
