{a, ad, b, bd} = FieldVariables[{1, 2}]; v = {a, ad, b, bd};
aA = I ad ** ad ** b; aAd = -I bd ** a ** a;
na = ad ** a; nbb = bd ** b;
BosonicNormalOrder[aA ** aAd - aAd ** aA - (na ** (na - 1) - (4 na + 2) ** nbb), v]
