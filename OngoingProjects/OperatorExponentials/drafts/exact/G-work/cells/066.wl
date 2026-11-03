{a, ad} = FieldVariables[{1}];
Simplify[SqueezeOperator[2 eta t, Infinity] /. Power[E, x_] :> Power[E, Expand[x]], t > 0]
