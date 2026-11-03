{a, ad} = FieldVariables[{1}];
dg = al ad - Conjugate[al] a;
{BosonicExpSimplify[Exp[-dg] ** a ** Exp[dg]], BosonicExpSimplify[Exp[-dg] ** ad ** Exp[dg]], BosonicExpSimplify[Exp[dg] ** a ** Exp[-dg]]}
