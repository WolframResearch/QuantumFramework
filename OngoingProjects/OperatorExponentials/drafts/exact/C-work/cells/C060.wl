{a, ad} = FieldVariables[{1}];
{BosonicExpSimplify[Exp[-al ad + Conjugate[al] a] ** a ** Exp[al ad - Conjugate[al] a]],
 BosonicExpSimplify[Exp[-al ad + Conjugate[al] a] ** (ad ** a) ** Exp[al ad - Conjugate[al] a]]}
