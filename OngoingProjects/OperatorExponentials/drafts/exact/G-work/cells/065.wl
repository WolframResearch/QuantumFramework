{a, ad} = FieldVariables[{1}];
s = SqueezeOperator[r E^(I th), Infinity]; sd = SqueezeOperator[-r E^(I th), Infinity];
tr = {a -> BosonicExpSimplify[sd ** a ** s] + al, ad -> BosonicExpSimplify[sd ** ad ** s] + Conjugate[al]};
ev[o_] := BosonicVEV[o /. tr, {a, ad}];
FullSimplify[ev[ad ** a ** ad ** a] - ev[ad ** a]^2 - ((Cosh[r] al - E^(I th) Sinh[r] Conjugate[al]) (Cosh[r] Conjugate[al] - E^(-I th) Sinh[r] al) + 2 Sinh[r]^2 Cosh[r]^2), r > 0 && th \[Element] Reals]
