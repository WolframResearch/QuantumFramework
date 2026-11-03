{a, ad} = FieldVariables[{1}];
diag[n_] := BosonicMatrixElement[{n, n}, DisplacementOperator[g, Infinity, "Ordering" -> "Normal"]] /. Conjugate[g] -> gb;
{Union[Table[Simplify[diag[n] - Exp[-g gb/2] LaguerreL[n, g gb]], {n, 0, 8}]],
 Simplify[Exp[-y/2] Sum[(1 - x) x^n LaguerreL[n, y], {n, 0, Infinity}], 0 < x < 1]}
