{a, ad} = FieldVariables[{1}];
d[x_] := DisplacementOperator[x, Infinity];
BosonicExpSimplify[d[w1] ** d[w2] ** d[-w1] ** d[-w2]]
