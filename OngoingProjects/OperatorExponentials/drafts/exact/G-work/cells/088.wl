{a, ad, b, bd} = FieldVariables[{1, 2}];
el = BosonicMatrixElement[{{m + q, m}, {m + q + 1, m + 1}}, a ** b];
{el, FullSimplify[eta^(m + 1)/Sqrt[(m + 1)! (m + q + 1)!] el - eta eta^m/Sqrt[m! (m + q)!], m >= 0 && q >= 0],
 BosonicMatrixElement[{{m + q, m}, {m + q, m}}, ad ** a - bd ** b],
 FullSimplify[Sum[x^(2 m)/(m! (m + q)!), {m, 0, Infinity}] - x^(-q) BesselI[q, 2 x], x > 0 && q \[Element] Integers && q >= 0]}
