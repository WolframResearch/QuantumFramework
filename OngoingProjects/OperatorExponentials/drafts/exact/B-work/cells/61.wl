{a, ad} = FieldVariables[{1}];
{x1, x2} = {E^(-I t1 w) ad, E^(I t2 w) a};
BosonicNormalOrder[x2 ** x1 - x1 ** x2, {a, ad}]
