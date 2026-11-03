{a, ad} = FieldVariables[{1}];
k0 = (ad ** a + 1/2)/2; kp = ad ** ad/2; km = a ** a/2;
one = {BosonicNormalOrder[k0 ** k0 - (kp ** km + km ** kp)/2, {a, ad}],
   BosonicNormalOrder[km ** kp - kp ** km - 2 k0, {a, ad}], BosonicNormalOrder[k0 ** kp - kp ** k0 - kp, {a, ad}]};
{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
n1 = a1d ** a1; n2 = a2d ** a2;
kk0 = (n1 + n2 + 1)/2; kkp = a1d ** a2d; kkm = a1 ** a2;
two = BosonicNormalOrder[kk0 ** kk0 - (kkp ** kkm + kkm ** kkp)/2 - ((n1 - n2) ** (n1 - n2) - 1)/4, {a1, a1d, a2, a2d}];
{one, two}
