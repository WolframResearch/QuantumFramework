{a, ad} = FieldVariables[{1}];
{BosonicExpSimplify[Exp[I Pi ad ** a] ** DisplacementOperator[-al, Infinity] ** Exp[-I Pi ad ** a]],
 BosonicExpSimplify[DisplacementOperator[al, Infinity] ** DisplacementOperator[al, Infinity]]}
