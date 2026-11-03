{a, ad} = FieldVariables[{1}];
toSym[e_] := e /. {ad -> zb, a -> z} /. GeneralizedPower[NonCommutativeMultiply, x_, n_] :> x^n /. NonCommutativeMultiply -> Times;
fld = {{ad, -1, 1}, {a, 1, 4}, {ad, -1, 2}, {a, 1, 5}};
tProd = toSym[BosonicNormalOrder[NonCommutativeMultiply @@ (E^(I #[[2]] #[[3]] w) #[[1]] & /@ SortBy[fld, -Last[#] &]), {a, ad}]];
{ann, cre} = Flatten[Position[fld[[All, 1]], #, {1}]] & /@ {a, ad};
dop[e_] := Sum[If[fld[[q, 3]] > fld[[p, 3]], E^(I (fld[[q, 3]] - fld[[p, 3]]) w), 0] D[e, x[q], x[p]], {q, ann}, {p, cre}];
eD = Sum[Nest[dop, Product[x[j], {j, 4}], k]/k!, {k, 0, 2}];
Expand[tProd - (eD /. x[j_] :> E^(I fld[[j, 2]] fld[[j, 3]] w) If[fld[[j, 1]] === a, z, zb])]
