{a, ad} = FieldVariables[{1}];
q = Sqrt[hb/(2 w)] (a + ad); p = I Sqrt[hb w/2] (ad - a);
PowerExpand[{BosonicNormalOrder[a ** ad - ad ** a, {a, ad}], BosonicNormalOrder[q ** p - p ** q, {a, ad}], BosonicNormalOrder[(p ** p + w^2 q ** q)/2, {a, ad}]}]
