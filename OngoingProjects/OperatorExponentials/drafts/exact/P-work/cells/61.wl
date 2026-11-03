{a, ad} = FieldVariables[{1}];
d[x_] := DisplacementOperator[x, Infinity] /. Conjugate[s_Plus] :> Map[Conjugate, s];
{Simplify[BosonicExpSimplify[d[b dt] ** d[al]]/d[al + b dt], Element[dt, Reals]],
 Exp[-g Conjugate[g]/2] BosonicMatrixElement[{n, 0}, Exp[g ad]]}
