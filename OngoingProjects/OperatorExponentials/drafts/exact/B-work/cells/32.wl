{a, ad} = FieldVariables[{1}];
{sN, fN} = {SqueezeOperator[xi, Infinity, "Ordering" -> "Normal"], DisplacementOperator[-f, Infinity, "Ordering" -> "Normal"]};
mid = BosonicExpOrder[BosonicExpBraid[Exp[E^(-I phi) Tanh[r] a ** a/2] ** Exp[-f ad]]];
c = Sech[r]^(1/2) E^(-f fb/2) BosonicMatrixElement[{0, 0}, Evaluate[mid]];
{ph, ps} = {Cosh[r], E^(I phi) Sinh[r]};
formula = (ph ph)^(-1/4) Exp[(f^2 Conjugate[ps]/ph - 2 f fb + fb^2 ps/ph)/4];
cplx = {f -> fx + I fy, fb -> fx - I fy};
{sN, fN, mid, c, FullSimplify[ComplexExpand[Abs[c /. cplx]] - ComplexExpand[formula /. cplx], r > 0]}
