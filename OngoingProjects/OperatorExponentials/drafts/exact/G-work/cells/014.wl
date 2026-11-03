{a, ad} = FieldVariables[{1}];
{BosonicExpSimplify[Exp[-(al ad - Conjugate[al] a)] ** a ** Exp[al ad - Conjugate[al] a]],
 Table[Simplify[BosonicMatrixElement[{m, 0}, DisplacementOperator[al]] - Exp[-Abs[al]^2/2] al^m/Sqrt[m!]], {m, 0, 8}]}
