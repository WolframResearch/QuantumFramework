{a, ad} = FieldVariables[{1}];
s = SqueezeOperator[r E^(I th), Infinity]; sd = SqueezeOperator[-r E^(I th), Infinity];
tr = {a -> BosonicExpSimplify[sd ** a ** s] + al, ad -> BosonicExpSimplify[sd ** ad ** s] + Conjugate[al]};
ev[o_] := BosonicVEV[o /. tr, {a, ad}];
FullSimplify[{ev[a] - al, ev[a ** a] - (al^2 - E^(I th) Sinh[r] Cosh[r]), ev[ad ** a] - (al Conjugate[al] + Sinh[r]^2)}, r > 0 && th \[Element] Reals]
