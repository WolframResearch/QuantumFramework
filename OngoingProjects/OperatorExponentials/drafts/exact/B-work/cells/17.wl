{a, ad} = FieldVariables[{1}];
tr = Integrate[Exp[(E^-b - 1) (x^2 + y^2)]/Pi, {x, -Infinity, Infinity}, {y, -Infinity, Infinity}, Assumptions -> b > 0];
{BosonicExpOrder[Exp[-b ad ** a]], tr, Simplify[tr - Sum[E^(-b n), {n, 0, Infinity}], b > 0]}
