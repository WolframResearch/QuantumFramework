{a, ad} = FieldVariables[{1}];
s = SqueezeOperator[r E^(I th), Infinity]; sd = SqueezeOperator[-r E^(I th), Infinity];
sa = BosonicExpSimplify[sd ** a ** s]; sad = BosonicExpSimplify[sd ** ad ** s];
x1 = (sa + sad)/2; x2 = (sa - sad)/(2 I);
ev[o_] := BosonicVEV[o, {a, ad}];
cov = ev[x1 ** x2 + x2 ** x1]/2 - ev[x1] ev[x2];
vm = {{ev[x1 ** x1] - ev[x1]^2, cov}, {cov, ev[x2 ** x2] - ev[x2]^2}};
FullSimplify[{cov + Sinh[2 r] Sin[th]/4, Det[vm] - 1/16}, r > 0 && th \[Element] Reals]
