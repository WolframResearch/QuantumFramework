{a, ad} = FieldVariables[{1}];
pm = (1 - E^-x) BosonicMatrixElement[{m, m}, Exp[-x ad ** a]];
{pm, FullSimplify[pm == nb^m/(1 + nb)^(m + 1) /. nb -> 1/(E^x - 1), x > 0 && m >= 0],
 Simplify[Sum[pm, {m, 0, Infinity}, Assumptions -> x > 0]], Simplify[Sum[Exp[-x (m + 1/2)], {m, 0, Infinity}, Assumptions -> x > 0] - E^(-x/2)/(1 - E^-x)]}
