{a, ad, b, bd} = FieldVariables[{1, 2}]; v = {a, ad, b, bd};
kp = ad ** bd; km = a ** b; k0 = (ad ** a + bd ** b + 1)/2; dn = ad ** a - bd ** b;
cm[x_, y_] := BosonicNormalOrder[x ** y - y ** x, v];
{cm[k0, kp] - BosonicNormalOrder[kp, v], cm[k0, km] + BosonicNormalOrder[km, v], cm[km, kp] - BosonicNormalOrder[2 k0, v], cm[dn, kp], cm[dn, km]}
