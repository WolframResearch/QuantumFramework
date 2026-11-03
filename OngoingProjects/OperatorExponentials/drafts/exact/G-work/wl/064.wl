quad = ComplexExpand[(Cosh[r] l + E^(I th) Sinh[r] Conjugate[l]) (Cosh[r] Conjugate[l] + E^(-I th) Sinh[r] l)/2 /. l -> x + I y];
mm = Simplify[D[quad, {{x, y}, 2}]];
lin = ComplexExpand[(Conjugate[l] dd - l Conjugate[dd]) /. {l -> x + I y, dd -> p + I q}];
kk = Simplify[D[lin, {{x, y}}]/I];
w = Simplify[(2 Pi/Sqrt[Det[mm]]) Exp[-kk . Inverse[mm] . kk/2]/Pi^2, r > 0 && th \[Element] Reals];
xy = ComplexExpand[{Re[#], Im[#]} &[E^(-I th/2) (p + I q)]];
target = 2/Pi Exp[-2 E^(2 r) xy[[1]]^2 - 2 E^(-2 r) xy[[2]]^2];
{Simplify[Det[mm]], Simplify[lin - I kk . {x, y}], FullSimplify[w - target, r > 0 && Element[{th, p, q}, Reals]]}
