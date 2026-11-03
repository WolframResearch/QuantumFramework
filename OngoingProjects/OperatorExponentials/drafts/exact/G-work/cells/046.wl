{a, ad, b, bd} = FieldVariables[{1, 2}]; v = {a, ad, b, bd};
j1 = (ad ** b + bd ** a)/2; j2 = (ad ** b - bd ** a)/(2 I); j3 = (ad ** a - bd ** b)/2; j0 = (ad ** a + bd ** b)/2;
cm[x_, y_] := BosonicNormalOrder[x ** y - y ** x, v];
{Simplify[{cm[j1, j2] - BosonicNormalOrder[I j3, v], cm[j2, j3] - BosonicNormalOrder[I j1, v], cm[j3, j1] - BosonicNormalOrder[I j2, v], cm[j0, j1], cm[j0, j2], cm[j0, j3],
  BosonicNormalOrder[j1 ** j1 + j2 ** j2 + j3 ** j3 - j0 ** j0 - j0, v]}],
 FullSimplify[{BosonicExpSimplify[Exp[-I ph j3] ** a ** Exp[I ph j3]], BosonicExpSimplify[Exp[-I ph j3] ** b ** Exp[I ph j3]]}, ph > 0]}
