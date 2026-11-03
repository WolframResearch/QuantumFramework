{a1, a2, b1, b2} = FieldVariables[{1, 2, 3, 4}][[{1, 3, 5, 7}]];
{ad1, ad2, bd1, bd2} = SuperDagger /@ {a1, a2, b1, b2};
z1 = {{x11, x12}, {x21, x22}}; z2 = {{y11, y12}, {y21, y22}};
cre = Total[Flatten[z2 {{ad1 ** bd1, ad1 ** bd2}, {ad2 ** bd1, ad2 ** bd2}}]];
ann = Total[Flatten[Conjugate[z1] {{a1 ** b1, a1 ** b2}, {a2 ** b1, a2 ** b2}}]];
ov[k_] := BosonicMatrixElement[{{0, 0, 0, 0}, {0, 0, 0, 0}}, NonCommutativeMultiply @@ Join[Table[ann, k], Table[cre, k]],
    "Generators" -> {a1, a2, b1, b2}]/(k!)^2;
target = Series[1/Det[IdentityMatrix[2] - s ConjugateTranspose[z1] . z2], {s, 0, 3}];
Expand[Table[ov[k] - SeriesCoefficient[target, k], {k, 1, 3}]]
