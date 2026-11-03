{a, ad} = FieldVariables[{1}];
g = al ad - Conjugate[al] a;
{BosonicExpSimplify[Exp[-g] ** a ** Exp[g]], BosonicExpSimplify[Exp[g] ** a ** Exp[-g]]}
