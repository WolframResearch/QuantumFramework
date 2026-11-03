{a, ad} = FieldVariables[{1}];
d[x_] := DisplacementOperator[x, Infinity];
Simplify[BosonicExpSimplify[d[-be] ** d[al]]/(Exp[(al Conjugate[be] - Conjugate[al] be)/2] d[al - be])]
