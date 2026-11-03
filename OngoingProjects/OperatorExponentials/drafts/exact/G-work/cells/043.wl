{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
t = {{t11, t12}, {t21, t22}};
out = t . {a1, a2}; outd = Conjugate[t] . {ad1, ad2};
Simplify[Table[BosonicNormalOrder[out[[j]] ** outd[[k]] - outd[[k]] ** out[[j]], {a1, ad1, a2, ad2}], {j, 2}, {k, 2}] - t . ConjugateTranspose[t]]
