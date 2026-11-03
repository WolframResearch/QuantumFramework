{a, ad} = FieldVariables[{1}];
q = (a + ad)/Sqrt[2]; p = -I (a - ad)/Sqrt[2];
h = (s11 q ** q + s12 (q ** p + p ** q) + s22 p ** p)/2;
heis = BosonicExpSimplify[Exp[I t h] ** a ** Exp[-I t h]];
{{fa, fb}, {fc, fd}} = MatrixExp[t {{0, 1}, {-1, 0}} . {{s11, s12}, {s12, s22}}];
{u, v} = {Coefficient[heis, a], Coefficient[heis, ad]};
Simplify[{u - (fa + fd + I (fc - fb))/2, v - (fa - fd + I (fc + fb))/2,
  ComplexExpand[u Conjugate[u] - v Conjugate[v]]}, t > 0 && s12^2 - s11 s22 > 0]
