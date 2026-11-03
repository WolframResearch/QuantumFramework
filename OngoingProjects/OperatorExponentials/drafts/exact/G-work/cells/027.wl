{a, ad} = FieldVariables[{1}]; nn = ad ** a;
{BosonicNormalOrder[nn ** nn, {a, ad}], BosonicNormalOrder[nn ** nn - ad ** ad ** a ** a, {a, ad}]}
