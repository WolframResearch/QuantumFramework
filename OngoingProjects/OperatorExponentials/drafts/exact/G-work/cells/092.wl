{a, ad, b, bd} = FieldVariables[{1, 2}];
h = mu ad ** bd + Conjugate[mu] a ** b;
term[k_] := (-I)^k (NonCommutativeMultiply @@ ConstantArray[h, k])/k!;
Table[Simplify[BosonicMatrixElement[{{j, j}, {0, 0}}, term[k]]], {k, 1, 2}, {j, 0, 2}]
