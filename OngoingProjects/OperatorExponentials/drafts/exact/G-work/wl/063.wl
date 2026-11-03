ov0 = Exp[-Abs[d]^2/2 - E^(I th) Conjugate[d]^2 Tanh[r]/2]/Sqrt[Cosh[r]];
c2m = Sqrt[(2 m)!]/(2^m m! Sqrt[Cosh[r]]) (-E^(I th) Tanh[r])^m;
sum = Exp[-Abs[d]^2/2] Sum[Conjugate[d]^(2 m)/Sqrt[(2 m)!] c2m, {m, 0, Infinity}];
{FullSimplify[sum - ov0, r > 0 && th \[Element] Reals],
 FullSimplify[ComplexExpand[ov0 Conjugate[ov0]/Pi - Exp[-Abs[d]^2 - Tanh[r] Re[E^(-I th) d^2]]/(Pi Cosh[r]) /. d -> x + I y], r > 0 && Element[{th, x, y}, Reals]]}
