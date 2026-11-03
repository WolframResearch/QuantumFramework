{a, ad} = FieldVariables[{1}];
{Table[Simplify[BosonicMatrixElement[{m, 0}, DisplacementOperator[al]] - Exp[-Abs[al]^2/2] BosonicMatrixElement[{m, 0}, Exp[al ad]]], {m, 0, 8}],
 BosonicExpOrder[Exp[Conjugate[al] a] ** Exp[al ad]], BosonicMatrixElement[{0, 0}, BosonicExpOrder[Exp[Conjugate[al] a] ** Exp[al ad]]]}
