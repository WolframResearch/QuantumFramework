{a, ad} = FieldVariables[{1}];
{DisplacementOperator[al, Infinity, "Ordering" -> "Normal"], BosonicExpSimplify[Exp[-Conjugate[al] a] ** ad ** Exp[Conjugate[al] a]]}
