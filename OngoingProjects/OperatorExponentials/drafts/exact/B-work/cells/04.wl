{a, ad} = FieldVariables[{1}];
q = (a + ad)/Sqrt[2]; p = (a - ad)/(I Sqrt[2]);
BosonicNormalOrder[q ** p - p ** q, {a, ad}]
