{a, ad} = FieldVariables[{1}];
braid = BosonicExpBraid[Exp[-Conjugate[al] a] ** Exp[-tau/2 ad ** ad]];
f = Exp[-Abs[al]^2/2]/Sqrt[Cosh[r]] Exp[al z - tau/2 (z - Conjugate[al])^2];
coefs = Table[Sqrt[m!] SeriesCoefficient[f, {z, 0, m}], {m, 0, 6}] /. tau -> E^(I th) Tanh[r];
herm = Table[Exp[-Abs[al]^2/2 - Conjugate[al]^2 E^(I th) Tanh[r]/2]/Sqrt[m! Cosh[r]] (E^(I th) Tanh[r]/2)^(m/2) HermiteH[m, (Cosh[r] al + E^(I th) Sinh[r] Conjugate[al])/Sqrt[E^(I th) Sinh[2 r]]], {m, 0, 6}];
pts = {{r -> 3/7, th -> 6/5, al -> 2/5 - 7 I/10}, {r -> 6/5, th -> -2, al -> -1 + I/3}, {r -> 1/9, th -> 3, al -> 5/4 + 2 I}};
{DisplacementOperator[al, Infinity, "Ordering" -> "Normal"], SqueezeOperator[r E^(I th), Infinity, "Ordering" -> "Normal"], braid,
 Max[Abs[N[coefs /. pts, 30] - N[herm /. pts, 30]]]}
