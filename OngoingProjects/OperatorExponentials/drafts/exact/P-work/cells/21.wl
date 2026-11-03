{a, ad} = FieldVariables[{1}];
lag[m_, n_] := Exp[-g gb/2] If[m >= n, Sqrt[n!/m!] g^(m - n) LaguerreL[n, m - n, g gb], Sqrt[m!/n!] (-gb)^(n - m) LaguerreL[m, n - m, g gb]];
nrm[m_, n_] := BosonicMatrixElement[{m, n}, DisplacementOperator[g, Infinity, "Ordering" -> "Normal"]] /. Conjugate[g] -> gb;
Union[Flatten[Table[Simplify[nrm[m, n] - lag[m, n]], {m, 0, 6}, {n, 0, 6}]]]
