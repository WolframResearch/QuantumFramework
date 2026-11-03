{a, ad} = FieldVariables[{1}];
q = (a + ad)/Sqrt[2]; p = -I (a - ad)/Sqrt[2];
diag = Simplify[{BosonicMatrixElement[{n, n}, q ** q], BosonicMatrixElement[{n, n}, p ** p],
    BosonicMatrixElement[{n, n}, (q ** p + p ** q)/2]}, n >= 0 && Element[n, Integers]];
z[b_] := Sum[Exp[-b w (n + 1/2)], {n, 0, Infinity}];
cov = Sum[Exp[-b w (n + 1/2)] diag, {n, 0, Infinity}]/z[b];
{diag, Simplify[{z[b] - 1/(2 Sinh[b w/2]), z[2 b]/z[b]^2 - Tanh[b w/2],
   z[2 b]/z[b]^2 - 1/(2 (1/(Exp[b w] - 1)) + 1), cov - {Coth[b w/2]/2, Coth[b w/2]/2, 0}}, b > 0 && w > 0]}
