{a, ad} = FieldVariables[{1}];
g = (xi ad ** ad - Conjugate[xi] a ** a)/2;
{SqueezeOperator[xi, Infinity], BosonicExpSimplify[Exp[g] ** Exp[-g]] /. E^x_ :> E^Expand[x],
 Sum[x^(2 k + 1)/(2 k + 1), {k, 0, Infinity}]}
