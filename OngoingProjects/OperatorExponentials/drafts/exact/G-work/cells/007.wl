{a, ad} = FieldVariables[{1}]; nn = ad ** a;
{BosonicNormalOrder[nn ** (nn - 1) ** (nn - 2) - ad ** ad ** ad ** a ** a ** a, {a, ad}],
 Factor[BosonicMatrixElement[{n, n}, ad ** ad ** ad ** ad ** a ** a ** a ** a]]}
