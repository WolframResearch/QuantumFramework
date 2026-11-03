{a, ad} = FieldVariables[{1}]; nn = ad ** a;
pm = nb^m/(1 + nb)^(m + 1);
th[op_] := Sum[pm BosonicMatrixElement[{m, m}, op], {m, 0, Infinity}, Assumptions -> nb > 0];
{FullSimplify[BosonicMatrixElement[{al, al}, ad ** ad ** a ** a, "Basis" -> "Coherent"]/BosonicMatrixElement[{al, al}, nn, "Basis" -> "Coherent"]^2],
 FullSimplify[th[ad ** ad ** a ** a]/th[nn]^2, nb > 0],
 Simplify[FunctionExpand[BosonicMatrixElement[{m, m}, ad ** ad ** a ** a]/BosonicMatrixElement[{m, m}, nn]^2]],
 BosonicNormalOrder[nn ** (nn - 1), {a, ad}]}
