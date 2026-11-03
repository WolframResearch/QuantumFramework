{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
v = {a1, a1d, a2, a2d};
g = r (eb a1 ** a2 - e a1d ** a2d);
uN = BosonicExpOrder[Exp[g]];
b = Expand[BosonicExpSimplify[Exp[g] ** # ** Exp[-g], v, "Scalars" -> {r, e, eb}]] & /@ {a1, a2};
{phi, psi} = Table[Coefficient[bi, y], {bi, b}, {y, #}] & /@ {{a1, a2}, {a1d, a2d}} /. eb -> 1/e;
vac = uN /. {Exp[x_] /; ! FreeQ[x, a1] && FreeQ[x, a1d] :> 1, Power[s_, n_] /; s =!= E && ! FreeQ[n, a1d] :> 1} /. {a1d -> z1, a2d -> z2} /. NonCommutativeMultiply -> Times;
zz = Inverse[phi] . psi;
{uN, FullSimplify[(vac /. eb -> 1/e) - Det[phi . ConjugateTranspose[phi]]^(-1/4) Exp[-{z1, z2} . zz . {z1, z2}/2], r > 0]}
