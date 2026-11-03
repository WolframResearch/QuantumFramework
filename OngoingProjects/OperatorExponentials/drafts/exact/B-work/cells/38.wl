{a, ad} = FieldVariables[{1}];
nA = BosonicExpOrder[Exp[l ad ** a]] /. l -> Log[1 - lam];
sdag = SqueezeOperator[-xi, Infinity, "Ordering" -> "Normal"];
moved = BosonicExpSimplify[Exp[l ad ** a] ** Exp[y ad ** ad] ** Exp[-l ad ** a], {a, ad}];
vev = BosonicMatrixElement[{0, 0}, Evaluate[BosonicExpOrder[Exp[x a ** a] ** Exp[yy ad ** ad]]]];
yy1 = moved /. Exp[q_] :> q /. GeneralizedPower[NonCommutativeMultiply, ad, 2] -> 1;
detC = Det[{{Cosh[r], (lam - 1) E^(-I phi) Sinh[r]}, {(lam - 1) E^(I phi) Sinh[r], Cosh[r]}}];
{nA, sdag, moved, vev, FullSimplify[TrigToExp[(Sech[r] vev /. {x -> E^(-I phi) Tanh[r]/2, yy -> yy1} /. {y -> E^(I phi) Tanh[r]/2, l -> Log[1 - lam]}) - detC^(-1/2)], r > 0 && 0 < lam < 1 && Element[phi, Reals]]}
