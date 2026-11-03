{a, ad} = FieldVariables[{1}];
s = SqueezeOperator[r E^(I th), Infinity]; sd = SqueezeOperator[-r E^(I th), Infinity];
mu = Cosh[r]; nu = E^(I th) Sinh[r];
expo[e_] := Expand[e /. Power[E, x_] :> x];
lhs1 = BosonicExpSimplify[s ** DisplacementOperator[be, Infinity] ** sd];
lhs2 = BosonicExpSimplify[sd ** DisplacementOperator[al, Infinity] ** s];
FullSimplify[{expo[lhs1] - expo[DisplacementOperator[mu be - nu Conjugate[be], Infinity]], expo[lhs2] - expo[DisplacementOperator[mu al + nu Conjugate[al], Infinity]]}, r > 0 && th \[Element] Reals]
