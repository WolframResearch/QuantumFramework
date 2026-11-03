{a, ad} = FieldVariables[{1}];
b = u a + v ad; bd = Conjugate[u] ad + Conjugate[v] a;
BosonicNormalOrder[b ** bd - bd ** b, {a, ad}]
