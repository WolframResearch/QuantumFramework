{a, ad} = FieldVariables[{1}];
q = (1 - E^-x)/Pi BosonicMatrixElement[{al, al}, Exp[-x ad ** a], "Basis" -> "Coherent"];
{FullSimplify[q /. x -> Log[1 + 1/nb], nb > 0],
 Table[Simplify[BosonicMatrixElement[{m, m}, Exp[l ad] ** Exp[-Conjugate[l] a]] - LaguerreL[m, l Conjugate[l]]], {m, 0, 5}],
 Simplify[Sum[(1 - t) t^m LaguerreL[m, z], {m, 0, Infinity}] /. t -> nb/(1 + nb), nb > 0]}
