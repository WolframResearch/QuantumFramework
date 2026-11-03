{a, ad} = FieldVariables[{1}]; nn = ad ** a; X1 = (a + ad)/2;
ev[op_] := BosonicMatrixElement[{al, al}, op, "Basis" -> "Coherent"];
{FullSimplify[{ev[nn], ev[nn ** nn] - ev[nn]^2, ev[ad ** ad ** ad ** a ** a ** a], ev[X1 ** X1] - ev[X1]^2}],
 Table[Simplify[Abs[BosonicMatrixElement[{m, 0}, DisplacementOperator[al]]]^2 - Exp[-Abs[al]^2] Abs[al]^(2 m)/m!], {m, 0, 5}]}
