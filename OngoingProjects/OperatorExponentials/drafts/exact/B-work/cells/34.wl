{a, ad} = FieldVariables[{1}];
{DisplacementOperator[-f, Infinity, "Ordering" -> "Normal"],
 BosonicExpSimplify[Exp[fb a - f ad] ** a ** Exp[-(fb a - f ad)], {a, ad}]}
