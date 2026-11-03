{a, ad} = FieldVariables[{1}];
q = Sqrt[hb/2] (a + ad); p = -I Sqrt[hb/2] (a - ad);
{BosonicNormalOrder[a ** ad - ad ** a, {a, ad}], BosonicNormalOrder[q ** p - p ** q, {a, ad}], DisplacementOperator[al, Infinity]}
