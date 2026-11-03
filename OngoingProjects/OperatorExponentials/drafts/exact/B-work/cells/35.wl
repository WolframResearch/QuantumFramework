{a, ad} = FieldVariables[{1}];
kern = Sech[r]^(1/2) Exp[E^(-I phi) Tanh[r] z^2/2 + Sech[r] zs z - E^(I phi) Tanh[r] zs^2/2];
{SqueezeOperator[xi, Infinity, "Ordering" -> "Normal"],
 BosonicExpSimplify[Exp[(xb a ** a - x ad ** ad)/2] ** a ** Exp[-(xb a ** a - x ad ** ad)/2], {a, ad}],
 FullSimplify[Table[BosonicMatrixElement[{m, n}, SqueezeOperator[r E^(I phi)]] - Sqrt[m! n!] SeriesCoefficient[kern, {zs, 0, m}, {z, 0, n}], {m, 0, 4}, {n, 0, 4}], r > 0 && -Pi < phi < Pi]}
