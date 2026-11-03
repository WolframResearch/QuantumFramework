{a, ad} = FieldVariables[{1}]; v = {a, ad};
q1 = (ad ** ad + a ** a)/2; q2 = (ad ** ad - a ** a)/(2 I); k3 = (ad ** a + 1/2)/2;
cm[x_, y_] := BosonicNormalOrder[x ** y - y ** x, v];
Simplify[{cm[q1, q2] - BosonicNormalOrder[-4 I k3, v], cm[k3, q1] - BosonicNormalOrder[I q2, v], cm[q2, k3] - BosonicNormalOrder[I q1, v]}]
