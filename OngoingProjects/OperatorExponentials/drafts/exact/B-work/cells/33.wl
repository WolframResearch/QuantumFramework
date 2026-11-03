{a, ad} = FieldVariables[{1}];
tau = E^(I phi) Tanh[r]; taub = E^(-I phi) Tanh[r];
mid = BosonicExpOrder[BosonicExpBraid[Exp[taub a ** a/2] ** Exp[-f ad]]];
mv = FullSimplify[BosonicExpSimplify[Exp[l ad ** a] ** Exp[-f ad] ** Exp[-l ad ** a], {a, ad}] /. l -> -Log[Cosh[r]], r > 0];
u = Sech[r]^(1/2) E^(-f fb/2) E^(f^2 taub/2) Exp[-tau ad ** ad/2] ** mv ** Sech[r]^(ad ** a) ** Exp[taub (a ** a - 2 f a)/2] ** Exp[fb a];
kern = u /. Sech[r]^(ad ** a) -> Exp[Sech[r] zs z] /. {ad -> zs, a -> z} /. NonCommutativeMultiply -> Times;
{ph, ps} = {Cosh[r], E^(I phi) Sinh[r]};
expo = (Conjugate[ps]/ph z^2 + 2 z zs/ph - ps/ph zs^2)/2 + z (fb - Conjugate[ps] f/ph) - zs f/ph;
{mid, mv, FullSimplify[PowerExpand[Log[kern/(kern /. {z -> 0, zs -> 0})]] - expo, r > 0 && Element[phi, Reals]]}
