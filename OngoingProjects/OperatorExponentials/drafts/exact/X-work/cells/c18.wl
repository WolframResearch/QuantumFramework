{a, ad} = FieldVariables[{1}];
Sum[BosonicMatrixElement[{n, n}, Exp[lam ad ** a]], {n, 0, Infinity}, GenerateConditions -> True]
