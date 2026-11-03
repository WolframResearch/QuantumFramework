{a, ad} = FieldVariables[{1}];
{a1, ad1, a2, ad2} = FieldVariables[{1, 2}];
one = BosonicMatrixElement[{0, 0}, Exp[-Conjugate[z]/2 a ** a] ** Exp[-w/2 ad ** ad]];
zq = (zb11 a1 ** a1 + 2 zb12 a1 ** a2 + zb22 a2 ** a2)/2;
wq = (w11 ad1 ** ad1 + 2 w12 ad1 ** ad2 + w22 ad2 ** ad2)/2;
ov[k_] := BosonicMatrixElement[{{0, 0}, {0, 0}}, NonCommutativeMultiply @@ Join[Table[zq, k], Table[wq, k]]]/(k!)^2;
target = Series[Det[IdentityMatrix[2] - s {{w11, w12}, {w12, w22}} . {{zb11, zb12}, {zb12, zb22}}]^(-1/2), {s, 0, 3}];
{one, Expand[Table[ov[k] - SeriesCoefficient[target, k], {k, 1, 3}]]}
