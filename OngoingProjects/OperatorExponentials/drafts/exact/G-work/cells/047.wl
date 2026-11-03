{a, ad, b, bd} = FieldVariables[{1, 2}];
g = I th (ad ** b + bd ** a)/2;
bdOut = FullSimplify[BosonicExpSimplify[Exp[g] ** bd ** Exp[-g]], th > 0];
nN = 3;
{bdOut, Table[FullSimplify[BosonicMatrixElement[{{k, nN - k}, {0, 0}}, (NonCommutativeMultiply @@ ConstantArray[bdOut, nN])/Sqrt[nN!]] - Sqrt[Binomial[nN, k]] (I Sin[th/2])^k Cos[th/2]^(nN - k)], {k, 0, nN}]}
