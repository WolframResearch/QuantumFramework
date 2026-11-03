{a, ad} = FieldVariables[{1}];
toSym[e_] := e /. {ad -> zb, a -> z} /. GeneralizedPower[NonCommutativeMultiply, x_, n_] :> x^n /. NonCommutativeMultiply -> Times;
op = a ** ad ** a ** ad;
nsym = Expand[toSym[BosonicNormalOrder[op, {a, ad}]]];
kcoh = BosonicMatrixElement[{al, be}, op, "Basis" -> "Coherent"];
{nsym, Simplify[kcoh Exp[(Abs[al]^2 + Abs[be]^2)/2] - (nsym /. {zb -> Conjugate[al], z -> be}) Exp[Conjugate[al] be]]}
