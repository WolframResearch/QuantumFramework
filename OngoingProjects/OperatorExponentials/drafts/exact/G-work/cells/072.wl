{a, ad, b, bd} = FieldVariables[{1, 2}]; v = {a, ad, b, bd};
xi = r E^(I th); g = Conjugate[xi] a ** b - xi ad ** bd;
tf[o_] := BosonicExpSimplify[Exp[-g] ** o ** Exp[g]];
tr = {a -> tf[a], ad -> tf[ad], b -> tf[b], bd -> tf[bd]};
ev[o_] := BosonicVEV[o /. tr, v];
nN = Sinh[r]^2; na = ad ** a; nbb = bd ** b;
c = (a + b)/Sqrt[2]; cd = (ad + bd)/Sqrt[2]; x1 = (c + cd)/2; x2 = (c - cd)/(2 I);
FullSimplify[{ev[na] - nN, ev[nbb] - nN, ev[na ** na] - ev[na]^2 - nN (nN + 1), ev[na ** nbb] - ev[na] ev[nbb] - nN (nN + 1),
  ev[ad ** ad ** a ** a] ev[bd ** bd ** b ** b] - 4 nN^4, ev[na ** nbb] - (2 nN^2 + nN),
  (ev[x1 ** x1] - ev[x1]^2 - E^(-2 r)/4) /. th -> 0, (ev[x2 ** x2] - ev[x2]^2 - E^(2 r)/4) /. th -> 0}, r > 0 && th \[Element] Reals]
