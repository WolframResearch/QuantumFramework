{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
g = r Exp[-I th] a1 ** a2 - r Exp[I th] ad1 ** ad2;
act = Simplify[{BosonicExpSimplify[Exp[g] ** a1 ** Exp[-g]], BosonicExpSimplify[Exp[g] ** a2 ** Exp[-g]]}, r > 0];
u = Table[Coefficient[act[[i]], x], {i, 2}, {x, {a1, a2}}];
v = Table[Coefficient[act[[i]], x], {i, 2}, {x, {ad1, ad2}}];
z = Simplify[Inverse[u] . v]; w = Simplify[Conjugate[v] . Inverse[u], r > 0 && Element[th, Reals]];
{act, z, w, Inverse[u], Simplify[Det[u . ComplexExpand[ConjugateTranspose[u]]]^(-1/4), r > 0], FullSimplify[BosonicExpOrder[Exp[g]], r > 0 && Element[th, Reals]]}
