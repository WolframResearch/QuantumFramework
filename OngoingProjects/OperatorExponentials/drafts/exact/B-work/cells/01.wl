{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
v = {a1, a1d, a2, a2d};
af = f1 a1 + f2 a2; adg = g1 a1d + g2 a2d;
BosonicNormalOrder[#, v] & /@ {af ** adg - adg ** af, a1 ** a2 - a2 ** a1, a1d ** a2d - a2d ** a1d}
