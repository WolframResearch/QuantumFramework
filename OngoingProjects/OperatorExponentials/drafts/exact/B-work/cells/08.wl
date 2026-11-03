{a, ad} = FieldVariables[{1}];
toSym[e_] := e /. {ad -> zb, a -> z} /. GeneralizedPower[NonCommutativeMultiply, x_, n_] :> x^n /. NonCommutativeMultiply -> Times;
q = ad ** ad ** a + 2 ad + a ** a; s = toSym[q];
{Expand[toSym[BosonicNormalOrder[a ** q, {a, ad}]] - (z s + D[s, zb])],
 Expand[toSym[BosonicNormalOrder[q ** ad, {a, ad}]] - (zb s + D[s, z])]}
