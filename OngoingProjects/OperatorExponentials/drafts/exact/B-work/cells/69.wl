{a, ad} = FieldVariables[{1}];
{BosonicExpOrder[Exp[-x ad ** a]] /. x -> beta (w + mu),
 Integrate[Exp[(E^(-beta (w + mu)) - 1) (p^2 + q^2)]/Pi, {p, -Infinity, Infinity}, {q, -Infinity, Infinity}, Assumptions -> beta (w + mu) > 0]}
