{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
v = {a1, a1d, a2, a2d};
dbl[g_] := Table[Coefficient[Expand[BosonicExpSimplify[Exp[g] ** x ** Exp[-g], v, "Scalars" -> {r, th}]], y], {x, {a1, a2, a1d, a2d}}, {y, {a1, a2, a1d, a2d}}];
s = FullSimplify[dbl[r (a1 ** a1 - a1d ** a1d)/2] . dbl[I th (a1d ** a2 + a2d ** a1)] /. {r -> Log[3], th -> Pi/4}];
{phi, psi} = {s[[1 ;; 2, 1 ;; 2]], s[[1 ;; 2, 3 ;; 4]]};
sinv = ArrayFlatten[{{ConjugateTranspose[phi], -Transpose[psi]}, {-ConjugateTranspose[psi], Transpose[phi]}}];
{s, Simplify[{s[[3 ;; 4]] - ArrayFlatten[{{Conjugate[psi], Conjugate[phi]}}],
  phi . ConjugateTranspose[phi] - psi . ConjugateTranspose[psi] - IdentityMatrix[2],
  phi . Transpose[psi] - psi . Transpose[phi], s . sinv - IdentityMatrix[4],
  ConjugateTranspose[phi] . phi - Transpose[psi] . Conjugate[psi] - IdentityMatrix[2],
  ConjugateTranspose[phi] . psi - Transpose[psi] . Conjugate[phi]}]}
