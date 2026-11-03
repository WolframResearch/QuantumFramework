{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
BosonicNormalOrder[#, {a1, a1d, a2, a2d}] & /@ {a1 ** a2 - a2 ** a1, a1d ** a2d - a2d ** a1d}
