{a, ad} = FieldVariables[{1}];
q = (a + ad)/Sqrt[2]; p = -I (a - ad)/Sqrt[2];
h = (s11 q ** q + s12 (q ** p + p ** q) + s22 p ** p)/2;
qp[x_] := Simplify[{Coefficient[x, a] + Coefficient[x, ad], I (Coefficient[x, a] - Coefficient[x, ad])}/Sqrt[2]];
gen = qp /@ {BosonicNormalOrder[I (h ** q - q ** h), {a, ad}], BosonicNormalOrder[I (h ** p - p ** h), {a, ad}]};
flow = qp /@ {BosonicExpSimplify[Exp[I t h] ** q ** Exp[-I t h]], BosonicExpSimplify[Exp[I t h] ** p ** Exp[-I t h]]};
{gen, gen - {{0, 1}, {-1, 0}} . {{s11, s12}, {s12, s22}},
 Simplify[flow - MatrixExp[t {{0, 1}, {-1, 0}} . {{s11, s12}, {s12, s22}}], t > 0 && s12^2 - s11 s22 > 0]}
