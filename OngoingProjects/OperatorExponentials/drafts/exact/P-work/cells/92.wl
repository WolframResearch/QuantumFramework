{a, ad} = FieldVariables[{1}];
dm[m_, n_] := ComplexExpand[BosonicMatrixElement[{m, n}, DisplacementOperator[gx + I gy, Infinity, "Ordering" -> "Normal"]]];
f = Sum[u^m v^n dm[m, n] ComplexExpand[Conjugate[dm[m, n]]], {m, 0, 5}, {n, 0, 5}];
rhs = Normal[Series[Exp[-(gx^2 + gy^2) (1 - u) (1 - v)/(1 - u v)]/(1 - u v), {u, 0, 5}, {v, 0, 5}]];
{Simplify[Expand[f - rhs]], f /. {u -> 0, v -> 0}}
