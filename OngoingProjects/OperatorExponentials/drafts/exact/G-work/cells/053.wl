{a, ad} = FieldVariables[{1}]; v = {a, ad};
kp = ad ** ad/2; km = a ** a/2; k0 = (ad ** a + 1/2)/2;
cm[x_, y_] := BosonicNormalOrder[x ** y - y ** x, v];
{cm[k0, kp] - BosonicNormalOrder[kp, v], cm[k0, km] + BosonicNormalOrder[km, v], cm[km, kp] - BosonicNormalOrder[2 k0, v]}
