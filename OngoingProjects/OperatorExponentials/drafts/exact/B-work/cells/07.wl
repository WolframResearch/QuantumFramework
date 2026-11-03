{a, ad} = FieldVariables[{1}];
toSym[e_] := e /. {ad -> zb, a -> z} /. GeneralizedPower[NonCommutativeMultiply, x_, n_] :> x^n /. NonCommutativeMultiply -> Times;
op = a ** ad ** a ** ad;
nsym = Expand[toSym[BosonicNormalOrder[op, {a, ad}]]];
diag = BosonicMatrixElement[{al, al}, op, "Basis" -> "Coherent"] /. Abs[al]^2 -> al Conjugate[al];
{diag, Simplify[diag - (nsym /. {zb -> Conjugate[al], z -> al})]}
