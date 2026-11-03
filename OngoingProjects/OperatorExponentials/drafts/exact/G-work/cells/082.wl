{a, ad} = FieldVariables[{1}];
amp[m_, x_] := BosonicMatrixElement[{m, 0}, DisplacementOperator[x]];
Table[FullSimplify[Abs[amp[m, al] + E^(I ph) amp[m, -al]]^2/(2 (1 + Exp[-2 Abs[al]^2] Cos[ph])) - Exp[-Abs[al]^2] Abs[al]^(2 m)/m! Abs[1 + E^(I ph) (-1)^m]^2/(2 (1 + Exp[-2 Abs[al]^2] Cos[ph])), ph \[Element] Reals], {m, 0, 6}]
