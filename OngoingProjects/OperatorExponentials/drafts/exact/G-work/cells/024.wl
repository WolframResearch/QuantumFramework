{a, ad} = FieldVariables[{1}];
{Table[Simplify[Conjugate[BosonicMatrixElement[{m, 0}, DisplacementOperator[al]]] - Exp[-Abs[al]^2/2] Conjugate[al]^m/Sqrt[m!]], {m, 0, 6}],
 FullSimplify[{D[z^m/Sqrt[m!], z] - Sqrt[m] z^(m - 1)/Sqrt[(m - 1)!], z z^m/Sqrt[m!] - Sqrt[m + 1] z^(m + 1)/Sqrt[(m + 1)!]}, m >= 1],
 BosonicMatrixElement[{m - 1, m}, a], BosonicMatrixElement[{m + 1, m}, ad]}
