{a1, a1d, a2, a2d} = FieldVariables[{1, 2}];
g = th (Exp[I ph] a2d ** a1 - Exp[-I ph] a1d ** a2);
out = {BosonicExpSimplify[Exp[g] ** a1d ** Exp[-g]], BosonicExpSimplify[Exp[g] ** a2d ** Exp[-g]]};
u = MatrixExp[{{0, -th Exp[-I ph]}, {th Exp[I ph], 0}}];
{out, Simplify[Table[Coefficient[out[[j]], {a1d, a2d}[[k]]] - u[[k, j]], {j, 2}, {k, 2}], th > 0 && Element[ph, Reals]]}
