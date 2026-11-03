{a, ad} = FieldVariables[{1}];
g = (r Exp[-I th] a ** a - r Exp[I th] ad ** ad)/2;
Simplify[BosonicExpSimplify[Exp[-g] ** a ** Exp[g]], r > 0]
