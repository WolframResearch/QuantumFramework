{a, ad, b, bd} = FieldVariables[{1, 2}];
tau = E^(I th) Tanh[r];
cm = Table[BosonicMatrixElement[{{m, m}, {0, 0}}, (NonCommutativeMultiply @@ Join[ConstantArray[ad, m], ConstantArray[bd, m]])/m!] (-tau)^m/Cosh[r], {m, 1, 4}];
ent = -Sum[(1 - t) t^m (Log[1 - t] + m Log[t]), {m, 0, Infinity}, Assumptions -> 0 < t < 1];
{FullSimplify[Abs[cm]^2 - Table[Sech[r]^2 Tanh[r]^(2 m), {m, 1, 4}], r > 0 && th \[Element] Reals],
 FullSimplify[(ent /. t -> nN/(nN + 1)) - ((nN + 1) Log[nN + 1] - nN Log[nN]), nN > 0], FullSimplify[Tanh[r]^2 - Sinh[r]^2/(Sinh[r]^2 + 1)]}
