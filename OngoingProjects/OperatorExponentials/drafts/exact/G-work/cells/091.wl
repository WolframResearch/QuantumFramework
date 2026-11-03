{a1, ad1, a2, ad2, a3, ad3} = FieldVariables[{1, 2, 3}]; v = {a1, ad1, a2, ad2, a3, ad3};
aA = l1 E^(-I w1 t) a1 + l2 E^(-I w2 t) a2 + l3 E^(-I w3 t) a3;
aAd = Conjugate[l1] E^(I w1 t) ad1 + Conjugate[l2] E^(I w2 t) ad2 + Conjugate[l3] E^(I w3 t) ad3;
x1 = (aA + aAd)/2; x2 = (aA - aAd)/(2 I);
Simplify[{BosonicNormalOrder[aA ** aAd - aAd ** aA, v], BosonicNormalOrder[x1 ** x2 - x2 ** x1, v]}]
