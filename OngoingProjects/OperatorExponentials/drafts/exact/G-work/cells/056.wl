{a, ad} = FieldVariables[{1}];
s = SqueezeOperator[r E^(I th), Infinity]; sd = SqueezeOperator[-r E^(I th), Infinity];
mu = Cosh[r]; nu = E^(I th) Sinh[r];
null = FullSimplify[BosonicNormalOrder[mu BosonicExpSimplify[sd ** a ** s] + nu BosonicExpSimplify[sd ** ad ** s], {a, ad}], r > 0 && th \[Element] Reals];
coef[m_] := FullSimplify[BosonicMatrixElement[{2 m, 0}, SqueezeOperator[r E^(I th)]] - Sqrt[(2 m)!]/(2^m m! Sqrt[Cosh[r]]) (-E^(I th) Tanh[r])^m, r > 0 && th \[Element] Reals];
{null, Table[coef[m], {m, 0, 5}], Table[BosonicMatrixElement[{2 m + 1, 0}, SqueezeOperator[r E^(I th)]], {m, 0, 2}],
 FullSimplify[Sum[(2 m)! Tanh[r]^(2 m)/(2^(2 m) (m!)^2 Cosh[r]), {m, 0, Infinity}], r > 0]}
