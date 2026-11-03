{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
v = {a1, a1d, a2, a2d};
toSym[e_] := e /. {a1d -> zs1, a2d -> zs2, a1 -> z1, a2 -> z2} /. GeneralizedPower[NonCommutativeMultiply, x_, n_] :> x^n /. NonCommutativeMultiply -> Times;
g = th (E^(I chi) a1d ** a2 - E^(-I chi) a2d ** a1);
conj[x_] := BosonicExpSimplify[Exp[g] ** x ** Exp[-g], v];
phi = Table[Coefficient[Expand[conj[x]], y], {x, {a1, a2}}, {y, {a1, a2}}];
nA[zs_, z_] := zs[[1]]^2 z[[2]] + 3 zs[[2]] z[[1]] z[[2]];
lhs = toSym[BosonicNormalOrder[conj[a1d ** a1d ** a2] + 3 conj[a2d ** a1 ** a2], v]];
{FullSimplify[phi, Element[th, Reals]], FullSimplify[lhs - nA[Conjugate[phi] . {zs1, zs2}, phi . {z1, z2}], Element[th | chi, Reals]]}
