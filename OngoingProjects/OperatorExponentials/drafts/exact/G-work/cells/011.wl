{a, ad} = FieldVariables[{1}]; nn = ad ** a;
pm = nb^m/(1 + nb)^(m + 1);
mom[op_] := FullSimplify[Sum[pm BosonicMatrixElement[{m, m}, op], {m, 0, Infinity}, Assumptions -> nb > 0], nb > 0];
fac[r_] := NonCommutativeMultiply @@ Join[ConstantArray[ad, r], ConstantArray[a, r]];
{Table[mom[fac[r]] - r! nb^r, {r, 1, 4}], mom[nn ** nn], FullSimplify[mom[nn ** nn] - mom[nn]^2]}
