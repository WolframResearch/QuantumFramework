{a, ad} = FieldVariables[{1}];
g = (xb a ** a - x ad ** ad)/2;
{BosonicExpSimplify[Exp[fb a - f ad] ** a ** Exp[-(fb a - f ad)], {a, ad}],
 BosonicExpSimplify[Exp[g] ** a ** Exp[-g], {a, ad}],
 BosonicExpSimplify[Exp[g] ** ad ** Exp[-g], {a, ad}]}
