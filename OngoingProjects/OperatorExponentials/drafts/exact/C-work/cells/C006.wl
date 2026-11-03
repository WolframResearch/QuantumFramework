{a, ad} = FieldVariables[{1}];
{DisplacementOperator[al, Infinity, "Ordering" -> "Normal"],
 DisplacementOperator[al, Infinity, "Ordering" -> "Antinormal"],
 BosonicExpSimplify[Exp[al ad - Conjugate[al] a] ** a ** Exp[-al ad + Conjugate[al] a]]}
