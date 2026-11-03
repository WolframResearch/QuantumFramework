{a, ad} = FieldVariables[{1}];
{DisplacementOperator[al, Infinity, "Ordering" -> "Normal"], DisplacementOperator[al, Infinity, "Ordering" -> "Antinormal"],
 BosonicExpSimplify[DisplacementOperator[al, Infinity] ** DisplacementOperator[-al, Infinity]]}
