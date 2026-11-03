{a, ad} = FieldVariables[{1}];
{BosonicExpOrder[Exp[lam ad ** a]], Simplify[Sum[Binomial[m, k] (E^lam - 1)^k, {k, 0, m}]],
 Table[1 + Sum[(-1)^k/k! BosonicMatrixElement[{m, m}, NonCommutativeMultiply @@ Join[ConstantArray[ad, k], ConstantArray[a, k]]], {k, 1, m}], {m, 0, 5}]}
