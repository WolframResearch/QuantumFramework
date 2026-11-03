{a, ad} = FieldVariables[{1}];
x = (a + ad)/2;
nord[j_] := Total[Table[Binomial[j, i] NonCommutativeMultiply @@ Join[ConstantArray[ad, i], ConstantArray[a, j - i]], {i, 0, j}]]/2^j;
pw[j_] := NonCommutativeMultiply @@ ConstantArray[x, j];
rule[nN_] := BosonicNormalOrder[pw[nN] - Sum[nN!/((nN - 2 k)! k! 8^k) If[nN - 2 k == 0, 1, nord[nN - 2 k]], {k, 0, Floor[nN/2]}], {a, ad}];
{Table[rule[nN], {nN, 1, 6}], Table[(2 l)!/(l! 8^l) - (2 l - 1)!!/4^l, {l, 1, 5}]}
