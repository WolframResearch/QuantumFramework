{a, ad} = FieldVariables[{1}]; nn = ad ** a;
{BosonicNormalOrder[nn ** (nn - 1) - ad ** ad ** a ** a, {a, ad}], FullSimplify[Exp[-I Pi m (m - 1)], m \[Element] Integers],
 Simplify[Exp[-I K t m (m - 1)] - Exp[-I K t m^2] Exp[I K t m]]}
