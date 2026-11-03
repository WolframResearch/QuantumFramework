{a, ad} = FieldVariables[{1}];
q = (a + ad)/Sqrt[2 w0]; p = -I Sqrt[w0/2] (a - ad);
BosonicNormalOrder[#, {a, ad}] & /@ {(p ** p/w0 + w0 q ** q)/4, (p ** q + q ** p)/4, (p ** p/w0 - w0 q ** q)/4}
