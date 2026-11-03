{a, ad} = FieldVariables[{1}];
FullSimplify[{BosonicExpSimplify[Exp[-I t w ad ** a] ** a ** Exp[I t w ad ** a], {a, ad}],
  BosonicExpSimplify[Exp[-I t w ad ** a] ** ad ** Exp[I t w ad ** a], {a, ad}]}, Element[t | w, Reals]]
