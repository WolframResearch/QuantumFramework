{a, ad} = FieldVariables[{1}];
kp = ad ** ad/2; km = a ** a/2; k0 = (ad ** a + 1/2)/2;
BosonicNormalOrder[#, {a, ad}] & /@ {k0 ** kp - kp ** k0 - kp, km ** kp - kp ** km - 2 k0, k0 ** k0 - (kp ** km + km ** kp)/2}
