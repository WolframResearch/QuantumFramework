{a, ad} = FieldVariables[{1}];
h = w ad ** a + g (f ad + fb a);
{BosonicNormalOrder[w (ad + g fb/w) ** (a + g f/w) - g^2 f fb/w - h, {a, ad}],
 BosonicNormalOrder[w (ad + g fb/w) ** (a + g f/w) + g^2 f fb/w - h, {a, ad}]}
