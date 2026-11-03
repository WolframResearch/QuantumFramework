{a, ad} = FieldVariables[{1}];
{Simplify[Table[BosonicMatrixElement[{n, 0}, DisplacementOperator[al]] - Exp[-Abs[al]^2/2] al^n/Sqrt[n!], {n, 0, 5}]],
 BosonicExpSimplify[a ** Exp[al ad]]}
