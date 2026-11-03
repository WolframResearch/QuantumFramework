{a, ad} = FieldVariables[{1}]; nn = ad ** a;
ov[x_, y_] := BosonicMatrixElement[{x, y}, 1, "Basis" -> "Coherent"];
{FullSimplify[BosonicMatrixElement[{be, al}, Exp[I ph nn], "Basis" -> "Coherent"] == ov[be, al Exp[I ph]], ph \[Element] Reals],
 Simplify[BosonicMatrixElement[{m, n}, Exp[I Pi nn]], n \[Element] Integers],
 Simplify[BosonicMatrixElement[{m, n}, Exp[2 I Pi nn]], n \[Element] Integers],
 Simplify[BosonicMatrixElement[{be, al}, Exp[I Pi nn], "Basis" -> "Coherent"] == ov[be, -al]]}
