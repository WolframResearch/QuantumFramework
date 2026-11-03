{a1, a2, b1, b2} = FieldVariables[{1, 2, 3, 4}][[{1, 3, 5, 7}]];
{ad1, ad2, bd1, bd2} = SuperDagger /@ {a1, a2, b1, b2};
zm = {{z11, z12}, {z21, z22}};
cre = Total[Flatten[zm {{ad1 ** bd1, ad1 ** bd2}, {ad2 ** bd1, ad2 ** bd2}}]];
ann = Total[Flatten[Conjugate[zm] {{a1 ** b1, a1 ** b2}, {a2 ** b1, a2 ** b2}}]];
nrm[k_] := BosonicMatrixElement[{{0, 0, 0, 0}, {0, 0, 0, 0}}, NonCommutativeMultiply @@ Join[Table[ann, k], Table[cre, k]],
    "Generators" -> {a1, a2, b1, b2}]/(k!)^2;
target = Series[1/Det[IdentityMatrix[2] - s ConjugateTranspose[zm] . zm], {s, 0, 3}];
Expand[Table[nrm[k] - SeriesCoefficient[target, k], {k, 1, 3}]]
