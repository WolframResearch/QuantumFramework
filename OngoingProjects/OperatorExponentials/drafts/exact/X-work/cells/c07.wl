{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
FullSimplify[BosonicExpOrder[Exp[r Exp[-I th] a1 ** a2 - r Exp[I th] ad1 ** ad2]], r > 0 && Element[th, Reals]]
