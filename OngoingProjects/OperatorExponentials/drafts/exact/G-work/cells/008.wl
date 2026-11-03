{a, ad} = FieldVariables[{1}];
X1 = (a + ad)/2; X2 = (a - ad)/(2 I);
{BosonicNormalOrder[X1 ** X2 - X2 ** X1, {a, ad}],
 Simplify[BosonicMatrixElement[{n, n}, X1 ** X1] - BosonicMatrixElement[{n, n}, X1]^2],
 Simplify[BosonicMatrixElement[{n, n}, X2 ** X2] - BosonicMatrixElement[{n, n}, X2]^2]}
