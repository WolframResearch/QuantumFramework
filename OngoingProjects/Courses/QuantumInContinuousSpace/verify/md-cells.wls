Off[Integrate::idiv];
WriteString[$Output,"=== 1.1 ===\n"];
\[Psi][x_] := Sqrt[\[Kappa]] Exp[-\[Kappa] Abs[x]]
dens = Simplify[Abs[\[Psi][x]]^2, \[Kappa] > 0 && x \[Element] Reals]
WriteString[$Output,"=== 1.2 ===\n"];
\[Phi][x_] := Sech[x]
norm2 = Integrate[Abs[\[Phi][x]]^2, {x, -Infinity, Infinity}]
\[Psi]N = \[Phi][x]/Sqrt[norm2];
Integrate[Abs[\[Psi]N]^2, {x, -Infinity, Infinity}]
WriteString[$Output,"=== 1.2 ===\n"];
Integrate[Abs[Exp[I k x]]^2, {x, -Infinity, Infinity}, Assumptions -> Element[{k, x}, Reals]]
WriteString[$Output,"=== 1.3 ===\n"];
\[Psi][x_] := Sech[x]/Sqrt[2];
Integrate[Abs[\[Psi][x]]^2, {x, -a, a}, Assumptions -> a > 0]
WriteString[$Output,"=== 1.4 ===\n"];
u[x_] := 2 \[Kappa]^(3/2) x Exp[-\[Kappa] x];
\[Phi][x_] := Exp[I q x] u[x];
ovl = Integrate[Conjugate[\[Phi][x]] u[x], {x, 0, Infinity}, Assumptions -> \[Kappa] > 0 && q \[Element] Reals]
Simplify[{Abs[ovl], Abs[ovl]^2}, \[Kappa] > 0 && q \[Element] Reals]
WriteString[$Output,"=== 1.5 ===\n"];
\[Psi][x_] := Sech[x - x0]/Sqrt[2];
mx = Integrate[x Abs[\[Psi][x]]^2, {x, -Infinity, Infinity}, Assumptions -> x0 \[Element] Reals];
mx2 = Integrate[x^2 Abs[\[Psi][x]]^2, {x, -Infinity, Infinity}, Assumptions -> x0 \[Element] Reals];
{mx, Simplify[Sqrt[mx2 - mx^2], x0 \[Element] Reals]}
WriteString[$Output,"=== 1.6 ===\n"];
\[Rho] = g[x, t] f[x, t];
j = (g[x, t] D[f[x, t], x] - f[x, t] D[g[x, t], x])/(2 I);
cont = D[\[Rho], t] + D[j, x] /.
   {Derivative[0, 1][f][x, t] -> (I/2) D[f[x, t], {x, 2}],
    Derivative[0, 1][g][x, t] -> -(I/2) D[g[x, t], {x, 2}]};
Simplify[cont]
WriteString[$Output,"=== 1.7 ===\n"];
\[Psi][x_] := Sqrt[\[Kappa]] Exp[-\[Kappa] Abs[x]];
grid = Subdivide[-12., 12., 6000]; dx = grid[[2]] - grid[[1]];
Total[Abs[(\[Psi][#] /. \[Kappa] -> 1.) & /@ grid]^2] dx
WriteString[$Output,"=== 1.8 ===\n"];
\[Psi][x_] := Sech[x] Exp[I \[Beta] x^2]/Sqrt[2];
\[Psi]c[x_] := Sech[x] Exp[-I \[Beta] x^2]/Sqrt[2];
Simplify[Abs[Exp[I \[Alpha]] \[Psi][x]]^2 - Abs[\[Psi][x]]^2, Element[{\[Alpha], \[Beta], x}, Reals]]
WriteString[$Output,"=== 1.8 ===\n"];
j = Simplify[(\[Psi]c[x] D[\[Psi][x], x] - \[Psi][x] D[\[Psi]c[x], x])/(2 I), Element[{\[Beta], x}, Reals]];
Simplify[j/(\[Psi]c[x] \[Psi][x]), Element[{\[Beta], x}, Reals]]
WriteString[$Output,"=== 1.9 ===\n"];
boost = Exp[I (v x - v^2 t/2)];
Simplify[Abs[boost]^2, Element[{v, x, t}, Reals]]
WriteString[$Output,"=== 1.9 ===\n"];
\[Psi]sol[x_] := Sech[x]/Sqrt[2];
Simplify[Abs[boost \[Psi]sol[x - v t]]^2 - Abs[\[Psi]sol[x - v t]]^2, Element[{v, x, t}, Reals]]
WriteString[$Output,"=== 2.1 ===\n"];
xop[f_] := x f;
pop[f_] := -I D[f, x];
{xop[Exp[I k x]], pop[Exp[I k x]]}
WriteString[$Output,"=== 2.2 ===\n"];
Simplify[xop[pop[f[x]]] - pop[xop[f[x]]]]
WriteString[$Output,"=== 2.3 ===\n"];
\[Psi][x_] := Sqrt[\[Kappa]] Exp[-\[Kappa] Abs[x]] Exp[I k x];
mp = Integrate[Conjugate[\[Psi][x]] (-I D[\[Psi][x], x]), {x, -Infinity, Infinity}, Assumptions -> \[Kappa] > 0 && k \[Element] Reals];
p2 = Integrate[Abs[D[\[Psi][x], x]]^2, {x, -Infinity, Infinity}, Assumptions -> \[Kappa] > 0 && k \[Element] Reals];
{Simplify[mp, \[Kappa] > 0 && k \[Element] Reals], Simplify[Sqrt[p2 - mp^2], \[Kappa] > 0 && k \[Element] Reals]}
WriteString[$Output,"=== 2.4 ===\n"];
\[Phi][x_] := Sech[x]/Sqrt[2];
\[Chi][x_] := Sech[x - 1] Exp[I x]/Sqrt[2];
Limit[-I Conjugate[\[Phi][x]] \[Chi][x], x -> Infinity] -
  Limit[-I Conjugate[\[Phi][x]] \[Chi][x], x -> -Infinity]
WriteString[$Output,"=== 2.4 ===\n"];
Chop[NIntegrate[Conjugate[\[Phi][x]] (-I D[\[Chi][x], x]), {x, -Infinity, Infinity}, WorkingPrecision -> 30] -
     NIntegrate[Conjugate[-I D[\[Phi][x], x]] \[Chi][x], {x, -Infinity, Infinity}, WorkingPrecision -> 30]]
WriteString[$Output,"=== 2.5 ===\n"];
FourierTransform[1, x, q]
WriteString[$Output,"=== 2.5 ===\n"];
FourierTransform[1/Sqrt[2 Pi] Exp[I p x], x, pp, FourierParameters -> {0, -1}]
WriteString[$Output,"=== 2.6 ===\n"];
\[Psi][x_] := Sech[x]/Sqrt[2];
\[Psi]T = FullSimplify[FourierTransform[\[Psi][x], x, p, FourierParameters -> {0, -1}]]
WriteString[$Output,"=== 2.6 ===\n"];
{Integrate[Abs[\[Psi]T]^2, {p, -Infinity, Infinity}],
 Integrate[p^2 Abs[\[Psi]T]^2, {p, -Infinity, Infinity}]}
WriteString[$Output,"=== 2.7 ===\n"];
\[Psi]sol[x_] := Sech[x]/Sqrt[2];
dxSol = Sqrt[Integrate[x^2 Abs[\[Psi]sol[x]]^2, {x, -Infinity, Infinity}]];
dpSol = Sqrt[Integrate[Abs[D[\[Psi]sol[x], x]]^2, {x, -Infinity, Infinity}]];
dxSol dpSol
WriteString[$Output,"=== 2.7 ===\n"];
\[Psi]G[x_] := Pi^(-1/4) \[Sigma]^(-1/2) Exp[-x^2/(2 \[Sigma]^2)];
dxG = Sqrt[Integrate[x^2 Abs[\[Psi]G[x]]^2, {x, -Infinity, Infinity}, Assumptions -> \[Sigma] > 0]];
dpG = Sqrt[Integrate[Abs[D[\[Psi]G[x], x]]^2, {x, -Infinity, Infinity}, Assumptions -> \[Sigma] > 0]];
Simplify[dxG dpG, \[Sigma] > 0]
WriteString[$Output,"=== 2.8 ===\n"];
sym[h_] := (xop[pop[h]] + pop[xop[h]])/2;
Simplify[sym[h[x]]]
WriteString[$Output,"=== 2.8 ===\n"];
u[x_] := Pi^(-1/4) Exp[-(x - 1)^2/2];
w[x_] := Pi^(-1/4) Exp[-(x + 1)^2/2] Exp[I x];
{Integrate[Conjugate[u[x]] sym[w[x]], {x, -Infinity, Infinity}] -
   Integrate[Conjugate[sym[u[x]]] w[x], {x, -Infinity, Infinity}],
 Integrate[Conjugate[u[x]] xop[pop[w[x]]], {x, -Infinity, Infinity}] -
   Integrate[Conjugate[xop[pop[u[x]]]] w[x], {x, -Infinity, Infinity}]}
WriteString[$Output,"=== 2.9 ===\n"];
(* box [0,L] with twisted boundary condition \[Psi](L) = e^{i\[Theta]} \[Psi](0) *)
bt = -I (Conjugate[pL] sL - Conjugate[p0] s0) /. {sL -> Exp[I \[Theta]] s0, pL -> Exp[I \[Theta]] p0};
Simplify[bt, \[Theta] \[Element] Reals]
WriteString[$Output,"=== 2.9 ===\n"];
(* half-line [0,\[Infinity]): \[Psi](\[Infinity]) = 0 but \[Psi](0) is free *)
-I (0 - Conjugate[p0] s0)
WriteString[$Output,"MD-CELLS-RAN-OK\n"];
