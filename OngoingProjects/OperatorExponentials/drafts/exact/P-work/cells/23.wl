{a, ad} = FieldVariables[{1}];
fs[m_, n_] := Exp[-g gb/2] Sqrt[m! n!] Sum[g^(m - r) (-gb)^(n - r)/(r! (m - r)! (n - r)!), {r, 0, Min[m, n]}];
nrm[m_, n_] := BosonicMatrixElement[{m, n}, DisplacementOperator[g, Infinity, "Ordering" -> "Normal"]] /. Conjugate[g] -> gb;
Union[Flatten[Table[Simplify[nrm[m, n] - fs[m, n]], {m, 0, 6}, {n, 0, 6}]]]
