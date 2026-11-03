{a, ad, b, bd} = FieldVariables[{1, 2}]; v = {a, ad, b, bd};
xi = r E^(I th); g = Conjugate[xi] a ** b - xi ad ** bd; mu = Cosh[r]; nu = E^(I th) Sinh[r]; tau = E^(I th) Tanh[r];
tf[o_] := BosonicExpSimplify[Exp[-g] ** o ** Exp[g]];
{FullSimplify[{BosonicNormalOrder[mu tf[a] + nu tf[bd], v], BosonicNormalOrder[mu tf[b] + nu tf[ad], v]}, r > 0 && th \[Element] Reals],
 BosonicNormalOrder[(ad ** a - bd ** b) ** g - g ** (ad ** a - bd ** b), v],
 Table[BosonicMatrixElement[{{m, m}, {0, 0}}, (NonCommutativeMultiply @@ Join[ConstantArray[ad, m], ConstantArray[bd, m]])/m!], {m, 1, 4}]}
