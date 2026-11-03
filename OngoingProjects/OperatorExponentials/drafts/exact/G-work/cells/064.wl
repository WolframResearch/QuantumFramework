{a, ad} = FieldVariables[{1}];
s = SqueezeOperator[r E^(I th), Infinity]; sd = SqueezeOperator[-r E^(I th), Infinity];
{FullSimplify[BosonicExpSimplify[sd ** DisplacementOperator[l, Infinity] ** s], r > 0 && th \[Element] Reals], BosonicMatrixElement[{0, 0}, DisplacementOperator[z]]}
