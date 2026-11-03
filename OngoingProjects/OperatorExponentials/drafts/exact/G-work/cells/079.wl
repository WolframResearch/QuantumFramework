{a, ad} = FieldVariables[{1}];
f[x_] := K x (x - 1);
cn[n_] := Exp[-rho^2/2] rho^n E^(I n p)/Sqrt[n!]; cnc[n_] := Exp[-rho^2/2] rho^n E^(-I n p)/Sqrt[n!];
term[m_] := FullSimplify[FunctionExpand[(cnc[n - m] cn[n] BosonicMatrixElement[{n - m, n}, NonCommutativeMultiply @@ ConstantArray[a, m]] Exp[I t (f[n - m] - f[n])]) /. n -> k + m], k >= 0 && rho > 0];
mom[m_] := Sum[Evaluate[term[m]], {k, 0, Infinity}];
xv = (a E^(-I v) + ad E^(I v))/2;
{Table[FullSimplify[mom[m] - (rho E^(I p))^m Exp[-I K t m (m - 1)] Exp[rho^2 (Exp[-2 I m K t] - 1)], rho > 0 && Element[{p, K, t}, Reals]], {m, 1, 3}],
 BosonicNormalOrder[xv ** xv, {a, ad}]}
