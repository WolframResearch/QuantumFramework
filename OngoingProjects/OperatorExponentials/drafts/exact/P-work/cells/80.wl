{a, ad} = FieldVariables[{1}];
adj[x_] := g (nu + 1) (ad ** x ** a - (ad ** a ** x + x ** ad ** a)/2) + g nu (a ** x ** ad - (a ** ad ** x + x ** a ** ad)/2);
BosonicNormalOrder[#, {a, ad}] & /@ {adj[a], adj[ad ** a]}
