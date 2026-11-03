{a, ad} = FieldVariables[{1}];
ov[x_, y_] := BosonicMatrixElement[{x, y}, 1, "Basis" -> "Coherent"];
nrm = ov[al, al] + ov[-al, -al] + E^(I ph) ov[al, -al] + E^(-I ph) ov[-al, al];
{FullSimplify[nrm - 2 (1 + Exp[-2 Abs[al]^2] Cos[ph]), ph \[Element] Reals],
 Simplify[{BosonicMatrixElement[{be, al}, a ** a, "Basis" -> "Coherent"]/ov[be, al], BosonicMatrixElement[{be, -al}, a ** a, "Basis" -> "Coherent"]/ov[be, -al]}]}
