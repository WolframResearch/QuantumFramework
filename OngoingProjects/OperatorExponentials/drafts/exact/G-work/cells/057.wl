{a, ad} = FieldVariables[{1}];
s = SqueezeOperator[r E^(I th), Infinity]; sd = SqueezeOperator[-r E^(I th), Infinity];
sa = BosonicExpSimplify[sd ** a ** s]; sad = BosonicExpSimplify[sd ** ad ** s];
x = (E^(-I v) sa + E^(I v) sad)/2;
FullSimplify[BosonicVEV[x ** x, {a, ad}] - BosonicVEV[x, {a, ad}]^2 - (Cosh[2 r] - Sinh[2 r] Cos[2 v - th])/4, r > 0 && th \[Element] Reals && v \[Element] Reals]
