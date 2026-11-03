{a, ad} = FieldVariables[{1}];
q = (a + ad)/Sqrt[2]; p = -I (a - ad)/Sqrt[2];
BosonicNormalOrder[p ** q - q ** p, {a, ad}]
