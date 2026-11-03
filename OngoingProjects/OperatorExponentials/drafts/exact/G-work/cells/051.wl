{a, ad} = FieldVariables[{1}];
s = SqueezeOperator[xi, Infinity];
{s, BosonicExpSimplify[s ** SqueezeOperator[-xi, Infinity]] /. Power[E, x_] :> Power[E, Expand[x]]}
