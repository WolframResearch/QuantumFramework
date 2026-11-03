{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
kp = ad1 ** ad2; km = a1 ** a2; n1 = ad1 ** a1; n2 = ad2 ** a2; k0 = (n1 + n2 + 1)/2;
BosonicNormalOrder[#, {a1, ad1, a2, ad2}] & /@ {k0 ** kp - kp ** k0 - kp, km ** kp - kp ** km - 2 k0,
  k0 ** k0 - (kp ** km + km ** kp)/2 - ((n1 - n2) ** (n1 - n2) - 1)/4}
