{a, ad} = FieldVariables[{1}];
der[m_, n_] := Exp[al alb/2]/Sqrt[m! n!] (-1)^m D[Exp[-al alb], {alb, m}, {al, n}];
nrm[m_, n_] := BosonicMatrixElement[{m, n}, DisplacementOperator[al, Infinity, "Ordering" -> "Normal"]] /. Conjugate[al] -> alb;
Union[Flatten[Table[Simplify[nrm[m, n] - der[m, n]], {m, 0, 5}, {n, 0, 5}]]]
