{a, ad} = FieldVariables[{1}];
toSym[e_] := e /. {ad -> zb, a -> z} /. GeneralizedPower[NonCommutativeMultiply, x_, n_] :> x^n /. NonCommutativeMultiply -> Times;
p1 = ad ** a ** a; p2 = ad ** ad ** a + a;
{s1, s2} = toSym /@ {p1, p2};
s12 = Expand[toSym[BosonicNormalOrder[p1 ** p2, {a, ad}]]];
{s12, Expand[s12 - Sum[D[s1, {z, k}] D[s2, {zb, k}]/k!, {k, 0, 4}]]}
