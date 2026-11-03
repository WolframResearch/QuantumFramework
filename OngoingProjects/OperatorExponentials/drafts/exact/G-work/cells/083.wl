{a, ad} = FieldVariables[{1}];
ov[x_, y_] := BosonicMatrixElement[{x, y}, 1, "Basis" -> "Coherent"];
dp[x_, y_] := Exp[(2 be Conjugate[-y] - Conjugate[2 be] (-y))/2] ov[x, 2 be - y];
kets = {{1, al}, {E^(I ph), -al}};
w = 2/Pi Total[Flatten[Outer[Conjugate[#1[[1]]] #2[[1]] dp[#1[[2]], #2[[2]]] &, kets, kets, 1]]]/(2 (1 + Exp[-2 al^2] Cos[ph]));
target = (Exp[-2 (x - al)^2 - 2 y^2] + Exp[-2 (x + al)^2 - 2 y^2] + 2 Exp[-2 x^2 - 2 y^2] Cos[4 al y + ph])/(Pi (1 + Exp[-2 al^2] Cos[ph]));
FullSimplify[ComplexExpand[w /. be -> x + I y] - target, Element[{al, ph, x, y}, Reals]]
