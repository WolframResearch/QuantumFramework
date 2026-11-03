{a, ad} = FieldVariables[{1}];
s = SqueezeOperator[r E^(I th), Infinity]; sd = SqueezeOperator[-r E^(I th), Infinity];
FullSimplify[{BosonicExpSimplify[sd ** a ** s], BosonicExpSimplify[sd ** ad ** s]}, r > 0 && th \[Element] Reals]
