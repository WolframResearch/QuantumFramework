{a, ad} = FieldVariables[{1}];
mm = {{m11, m12}, {m21, m22}}; id = IdentityMatrix[2];
logG = Log[Det[id - mm]/Det[id - t mm]];
tr[sg_] := Log[Det[id - mm]] + sg Sum[Tr[MatrixPower[mm, k]] t^k/k, {k, 1, 4}];
{Simplify[TrigToExp[BosonicExpSimplify[Exp[lam ad ** a] ** ad ** Exp[-lam ad ** a]]], lam > 0],
 Simplify[Normal[Series[logG - tr[+1], {t, 0, 4}]]], Simplify[Normal[Series[logG - tr[-1], {t, 0, 1}]]]}
