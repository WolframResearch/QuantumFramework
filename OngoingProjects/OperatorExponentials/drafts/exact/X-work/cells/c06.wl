{a, ad} = FieldVariables[{1}];
tau = Exp[I th] Tanh[r];
formula[m_] := If[EvenQ[m], Cosh[r]^(-1/2) (-tau/2)^(m/2) Sqrt[m!]/(m/2)!, 0];
{Table[FullSimplify[BosonicMatrixElement[{m, 0}, SqueezeOperator[r Exp[I th], Infinity]] - formula[m], r > 0 && Element[th, Reals]], {m, 0, 10}], FullSimplify[Sum[Abs[formula[2 k]]^2, {k, 0, Infinity}], r > 0 && Element[th, Reals]]}
