{a, ad} = FieldVariables[{1}];
d[x_] := DisplacementOperator[x, Infinity];
{BosonicExpSimplify[Exp[I Pi ad ** a] ** d[al] ** Exp[-I Pi ad ** a]], BosonicExpSimplify[d[al] ** d[al]],
 Simplify[ComplexExpand[-2 BosonicMatrixElement[{1, 1}, DisplacementOperator[2 (x + I y)]]]]}
