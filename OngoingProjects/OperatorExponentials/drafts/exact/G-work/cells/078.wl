{a, ad} = FieldVariables[{1}];
amp = BosonicMatrixElement[{n - 1, n}, a];
lhs = Exp[I t f[n - 1]] amp Exp[-I t f[n]];
rhs = amp Exp[-I t (f[n] - f[n - 1])];
kerr = f -> Function[x, K x (x - 1)];
{amp, Simplify[lhs - rhs], Simplify[(lhs /. kerr) - amp Exp[-2 I K t (n - 1)]],
 Simplify[(lhs /. kerr) - BosonicMatrixElement[{n - 1, n - 1}, Exp[-2 I K t ad ** a]] amp], Simplify[(lhs /. kerr)/(amp Exp[-2 I K t n])]}
