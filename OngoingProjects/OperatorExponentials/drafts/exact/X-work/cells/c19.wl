{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
vars = {a1, ad1, a2, ad2};
b = Table[Sum[u[i, j] {a1, a2}[[j]] + v[i, j] {ad1, ad2}[[j]], {j, 2}], {i, 2}];
bd = Table[Sum[Conjugate[u[i, j]] {ad1, ad2}[[j]] + Conjugate[v[i, j]] {a1, a2}[[j]], {j, 2}], {i, 2}];
uu = Array[u, {2, 2}]; vv = Array[v, {2, 2}];
qi = QuantumInterferometer[{"S"[1/3 Exp[I/2]] -> 1, "S"[1/4 Exp[-I]] -> 2, "BS"[3/5, 9/10] -> {1, 2}}, 2];
{al, be, c} = First[qi]["Bogoliubov"];
{uq, vq} = ConjugateTranspose /@ {al, be};
{Simplify[Table[BosonicNormalOrder[b[[i]] ** bd[[j]] - bd[[j]] ** b[[i]], vars], {i, 2}, {j, 2}] - (uu . ConjugateTranspose[uu] - vv . ConjugateTranspose[vv])],
 Simplify[Table[BosonicNormalOrder[b[[i]] ** b[[j]] - b[[j]] ** b[[i]], vars], {i, 2}, {j, 2}] - (uu . Transpose[vv] - vv . Transpose[uu])],
 FullSimplify[uq . ConjugateTranspose[uq] - vq . ConjugateTranspose[vq] - IdentityMatrix[2]], FullSimplify[uq . Transpose[vq] - vq . Transpose[uq]]}
