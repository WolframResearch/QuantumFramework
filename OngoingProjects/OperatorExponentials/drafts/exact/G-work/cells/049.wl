{a, ad, b, bd} = FieldVariables[{1, 2}];
g = I Pi (ad ** b + bd ** a)/4;
adO = BosonicExpSimplify[Exp[g] ** ad ** Exp[-g]]; bdO = BosonicExpSimplify[Exp[g] ** bd ** Exp[-g]];
st[nN_] := (NonCommutativeMultiply @@ Join[ConstantArray[adO, nN], ConstantArray[bdO, nN]])/nN!;
out[nN_] := Table[Simplify[BosonicMatrixElement[{{2 k, 2 nN - 2 k}, {0, 0}}, st[nN]] - I^nN Sqrt[(2 k)! (2 nN - 2 k)!]/(2^nN k! (nN - k)!)], {k, 0, nN}];
{adO, bdO, out[1], out[2], out[3], Simplify[BosonicMatrixElement[{{1, 1}, {0, 0}}, st[1]]], Simplify[BosonicMatrixElement[{{3, 1}, {0, 0}}, st[2]]]}
