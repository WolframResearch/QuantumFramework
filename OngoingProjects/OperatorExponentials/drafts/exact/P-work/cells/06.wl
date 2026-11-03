{a, ad} = FieldVariables[{1}];
{BosonicMatrixElement[{m, n}, a], BosonicMatrixElement[{m, n}, ad],
 Table[BosonicMatrixElement[{k, 0}, ad ** ad ** ad]/Sqrt[3!], {k, 0, 4}]}
