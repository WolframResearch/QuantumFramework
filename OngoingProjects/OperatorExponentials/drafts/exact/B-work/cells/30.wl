{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
v = {a1, a1d, a2, a2d};
dbl[g_] := Table[Coefficient[Expand[BosonicExpSimplify[Exp[g] ** x ** Exp[-g], v, "Scalars" -> {r, th}]], y], {x, {a1, a2, a1d, a2d}}, {y, {a1, a2, a1d, a2d}}];
s = FullSimplify[dbl[r (a1 ** a1 - a1d ** a1d)/2] . dbl[I th (a1d ** a2 + a2d ** a1)] /. {r -> Log[3], th -> Pi/4}];
{phi, psi} = {s[[1 ;; 2, 1 ;; 2]], s[[1 ;; 2, 3 ;; 4]]};
z = Inverse[phi] . psi; w = Conjugate[psi] . Inverse[phi];
Simplify[{z - Transpose[z], z . ConjugateTranspose[z] - IdentityMatrix[2] + Inverse[ConjugateTranspose[phi] . phi],
  ConjugateTranspose[w] . w - IdentityMatrix[2] + Inverse[phi . ConjugateTranspose[phi]],
  Norm[z], Norm[w], Norm[psi . Inverse[phi]]}]
