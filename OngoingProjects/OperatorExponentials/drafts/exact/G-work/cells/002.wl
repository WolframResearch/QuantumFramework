{a, ad} = FieldVariables[{1}];
FullSimplify[{BosonicExpSimplify[Exp[s ad ** a] ** ad ** Exp[-s ad ** a]], BosonicExpSimplify[Exp[s a ** a] ** ad ** Exp[-s a ** a]]}, s \[Element] Reals]
