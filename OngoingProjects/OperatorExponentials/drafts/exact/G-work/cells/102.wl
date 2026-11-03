{a, ad} = FieldVariables[{1}];
d2 = FullSimplify[BosonicMatrixElement[{al, al}, Exp[-x ad ** a], "Basis" -> "Coherent"] /. al -> rho E^(I phi), rho > 0 && phi \[Element] Reals];
{d2, Integrate[2 rho d2, {rho, 0, Infinity}, Assumptions -> x > 0], Sum[E^(-x m), {m, 0, Infinity}]}
