{a, ad} = FieldVariables[{1}];
{BosonicExpOrder[Exp[-x ad ** a]], Simplify[{1 - E^-x, E^-x - 1} /. x -> Log[1 + 1/nb], nb > 0]}
