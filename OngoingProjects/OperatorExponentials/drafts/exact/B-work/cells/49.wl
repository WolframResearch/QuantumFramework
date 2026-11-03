{a, ad} = FieldVariables[{1}];
{BosonicExpOrder[Exp[s ad ** a]] /. s -> c Log[xi], Integrate[xi^(c n - z - 1), {xi, 0, 1}]}
