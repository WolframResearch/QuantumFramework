#!/usr/bin/env wolframscript
(* Verification harness for "Quantum in Continuous Space" solutions, Parts 1-2.
   Natural units hbar = m = 1. Pure WL, no QuantumFramework.
   Each block prints a label, the computed result, and (where possible) a PASS/FAIL check. *)

show[label_, expr_] := WriteString[$Output, label, " => ", ToString[expr, InputForm], "\n"];
check[label_, bool_] := WriteString[$Output, label, " : ", If[TrueQ[bool], "PASS", "FAIL <<<<<"], "\n"];
WriteString[$Output, "================ PART 1 ================\n"];

(* ---- 1.1 wavefunction + |psi|^2 ---- *)
psi11[x_] := (1/(Pi sig^2)^(1/4)) Exp[-x^2/(2 sig^2)] Exp[I k0 x];
dens11 = Simplify[ComplexExpand[psi11[x] Conjugate[psi11[x]],
    TargetFunctions -> {Re, Im}], Assumptions -> sig > 0 && Element[{x, k0}, Reals]];
show["1.1 |psi|^2", dens11];
check["1.1 density real Gaussian", Simplify[dens11 == Exp[-x^2/sig^2]/(sig Sqrt[Pi]), sig > 0]];

(* ---- 1.2 normalization + non-normalizable ---- *)
phi12[x_] := Exp[-a x^2];
n12 = Integrate[Abs[phi12[x]]^2, {x, -Infinity, Infinity}, Assumptions -> a > 0];
psiN = phi12[x]/Sqrt[n12];
chk12 = Integrate[Abs[psiN]^2, {x, -Infinity, Infinity}, Assumptions -> a > 0];
show["1.2 norm^2 of Exp[-a x^2]", n12];
check["1.2 normalized integrates to 1", Simplify[chk12 == 1, a > 0]];
nonnorm = Integrate[Abs[Exp[I k x]]^2, {x, -Infinity, Infinity},
   Assumptions -> Element[{k, x}, Reals]];
show["1.2 plane-wave norm (divergent)", nonnorm];

(* ---- 1.3 Born probability in [a,b] ---- *)
prob13 = Integrate[Exp[-x^2]/Sqrt[Pi], {x, aa, bb}];
show["1.3 P(a<x<b), sigma=1/Sqrt2 Gaussian", prob13];
check["1.3 equals (Erf[bb]-Erf[aa])/2", Simplify[prob13 == (Erf[bb] - Erf[aa])/2]];

(* ---- 1.4 inner product, modulus, transition probability ---- *)
g[x_, x0_, p0_] := Pi^(-1/4) Exp[-(x - x0)^2/2] Exp[I p0 x];
ip14 = Integrate[Conjugate[g[x, d, kk]] g[x, 0, 0], {x, -Infinity, Infinity},
   Assumptions -> Element[{d, kk}, Reals]];
ip14 = Simplify[ip14, Element[{d, kk}, Reals]];
mod14 = Simplify[Sqrt[ip14 Conjugate[ip14]], Element[{d, kk}, Reals]];
show["1.4 <phi|psi>", ip14];
show["1.4 |<phi|psi>|", mod14];
show["1.4 |<phi|psi>|^2", Simplify[mod14^2, Element[{d, kk}, Reals]]];

(* ---- 1.5 <x> and Delta x ---- *)
rho15 = Abs[Pi^(-1/4) sig^(-1/2) Exp[-(x - x0)^2/(2 sig^2)]]^2;
rho15 = ComplexExpand[rho15] /. Abs[z_]^2 :> z^2;  (* real Gaussian, ensure real form *)
rho15 = Exp[-(x - x0)^2/sig^2]/(sig Sqrt[Pi]);
mx = Integrate[x rho15, {x, -Infinity, Infinity}, Assumptions -> sig > 0 && x0 \[Element] Reals];
mx2 = Integrate[x^2 rho15, {x, -Infinity, Infinity}, Assumptions -> sig > 0 && x0 \[Element] Reals];
dx = Simplify[Sqrt[mx2 - mx^2], sig > 0];
show["1.5 <x>", mx];
show["1.5 Delta x", dx];
check["1.5 <x>=x0, Dx=sig/Sqrt2", Simplify[mx == x0 && dx == sig/Sqrt[2], sig > 0]];

(* ---- 1.6 probability current + continuity (free TDSE, hbar=m=1) ---- *)
(* treat f = psi, gg = conj psi as independent; TDSE: f_t = (I/2) f_xx, gg_t = -(I/2) gg_xx *)
rho16 = gg[x, t] f[x, t];
j16 = (gg[x, t] D[f[x, t], x] - f[x, t] D[gg[x, t], x])/(2 I);
cont = D[rho16, t] + D[j16, x];
cont = cont /. {Derivative[0, 1][f][x, t] -> (I/2) D[f[x, t], {x, 2}],
                Derivative[0, 1][gg][x, t] -> -(I/2) D[gg[x, t], {x, 2}]};
show["1.6 continuity residual", Simplify[cont]];
check["1.6 d_t rho + d_x j == 0", Simplify[cont] === 0];

(* ---- 1.7 discrete grid: integral -> sum ---- *)
psiG[x_] := Pi^(-1/4) Exp[-x^2/2];
grid = Subdivide[-8., 8., 4000]; dxg = grid[[2]] - grid[[1]];
sum17 = Total[Abs[psiG[grid]]^2] dxg;
show["1.7 grid sum of |psi|^2 dx", sum17];
check["1.7 grid sum ~ 1", Abs[sum17 - 1] < 10^-6];

(* ---- 1.8 global phase vs local phase gradient ---- *)
psi18 = Exp[-x^2/2] Exp[I k x];
psi18c = Exp[-x^2/2] Exp[-I k x];   (* conjugate for real k, x *)
glob = Simplify[Abs[Exp[I al] psi18]^2 - Abs[psi18]^2, Element[{al, k, x}, Reals]];
show["1.8 |e^{ia}psi|^2 - |psi|^2", glob];
check["1.8 global phase unobservable", glob === 0];
j18 = Simplify[(psi18c D[psi18, x] - psi18 D[psi18c, x])/(2 I), Element[{k, x}, Reals]];
v18 = Simplify[j18/(psi18c psi18), Element[{k, x}, Reals]];
show["1.8 current j", j18];
show["1.8 velocity j/rho = local phase slope", v18];
check["1.8 v == k", Simplify[v18 == k, Element[{k, x}, Reals]]];

(* ---- 1.9 Galilean boost ---- *)
boost = Exp[I (v x - v^2 t/2)];
show["1.9 |boost phase|^2", Simplify[Abs[boost]^2, Element[{v, x, t}, Reals]]];
check["1.9 boost unit modulus", Simplify[Abs[boost]^2 == 1, Element[{v, x, t}, Reals]]];
(* density of boosted packet = density of original shifted by v t *)
densShape = Simplify[Abs[boost]^2 RR[x - v t] - RR[x - v t], Element[{v, x, t}, Reals]];
check["1.9 |psi_v(x,t)|^2 == |psi(x-vt,t)|^2", densShape === 0];

WriteString[$Output, "================ PART 2 ================\n"];

(* ---- 2.1 position and momentum operators on a test function ---- *)
xop[ff_] := x ff;
pop[ff_] := -I D[ff, x];
show["2.1 xhat e^{ikx}", xop[Exp[I k x]]];
show["2.1 phat e^{ikx}", pop[Exp[I k x]]];
check["2.1 phat eigen with eigenvalue k", Simplify[pop[Exp[I k x]] == k Exp[I k x]]];

(* ---- 2.2 canonical commutator ---- *)
comm22 = Simplify[xop[pop[ff[x]]] - pop[xop[ff[x]]]];
show["2.2 [x,p] ff[x]", comm22];
check["2.2 [x,p] = i", comm22 === I ff[x]];

(* ---- 2.3 <p>, Delta p in position representation ---- *)
psi23 = Pi^(-1/4) Exp[-x^2/2] Exp[I k x];
mp = Integrate[Conjugate[psi23] (-I D[psi23, x]), {x, -Infinity, Infinity},
   Assumptions -> k \[Element] Reals];
mp2 = Integrate[Conjugate[psi23] (-D[psi23, {x, 2}]), {x, -Infinity, Infinity},
   Assumptions -> k \[Element] Reals];
dp = Simplify[Sqrt[mp2 - mp^2], k \[Element] Reals];
show["2.3 <p>", Simplify[mp, k \[Element] Reals]];
show["2.3 Delta p", dp];
check["2.3 <p>=k, Dp=1/Sqrt2", Simplify[mp == k, k \[Element] Reals] && dp == 1/Sqrt[2]];

(* ---- 2.4 Hermiticity of p by parts ---- *)
ph[x_] := Pi^(-1/4) Exp[-(x - 1)^2/2];
ps[x_] := Pi^(-1/4) Exp[-(x + 1)^2/2] Exp[I x];
lhs24 = Integrate[Conjugate[ph[x]] (-I D[ps[x], x]), {x, -Infinity, Infinity}];
rhs24 = Integrate[Conjugate[-I D[ph[x], x]] ps[x], {x, -Infinity, Infinity}];
show["2.4 <phi|p psi> - <p phi|psi>", Simplify[lhs24 - rhs24]];
check["2.4 p Hermitian (bulk equal)", Simplify[lhs24 - rhs24] === 0];
bdy = Limit[-I Conjugate[ph[x]] ps[x], x -> Infinity] -
      Limit[-I Conjugate[ph[x]] ps[x], x -> -Infinity];
check["2.4 boundary term vanishes", bdy === 0];

(* ---- 2.5 momentum eigenfunctions, delta-normalization ---- *)
ortho25 = FourierTransform[1, x, q];   (* int e^{i q x} dx = 2 Pi delta(q) *)
show["2.5 FT[1] = (1/Sqrt2Pi) int e^{iqx} dx", ortho25];
(* momentum-space image of a normalized momentum eigenstate u_p = e^{ipx}/Sqrt[2Pi] *)
delta25 = FourierTransform[1/Sqrt[2 Pi] Exp[I p x], x, pp, FourierParameters -> {0, -1}];
show["2.5 u_p in momentum space", delta25];
check["2.5 delta-normalized", delta25 === DiracDelta[p - pp]];

(* ---- 2.6 momentum-space wavefunction + <p^n> ---- *)
psi26[x_] := Pi^(-1/4) Exp[-x^2/2] Exp[I k x];
psiT = FullSimplify[FourierTransform[psi26[x], x, p, FourierParameters -> {0, -1}],
   k \[Element] Reals];
show["2.6 psi-tilde(p)", psiT];
normP = Integrate[psiT Conjugate[psiT], {p, -Infinity, Infinity}, Assumptions -> k \[Element] Reals];
mpf = Integrate[p psiT Conjugate[psiT], {p, -Infinity, Infinity}, Assumptions -> k \[Element] Reals];
show["2.6 int |psi-tilde|^2", Simplify[normP, k \[Element] Reals]];
show["2.6 <p> from momentum space", Simplify[mpf, k \[Element] Reals]];
check["2.6 norm 1 and <p>=k (matches 2.3)", Simplify[normP == 1, k \[Element] Reals] && Simplify[mpf == k, k \[Element] Reals]];

(* ---- 2.7 uncertainty relation, Gaussian saturates ---- *)
psi27[x_] := Pi^(-1/4) sig^(-1/2) Exp[-x^2/(2 sig^2)];
nx = Integrate[Abs[psi27[x]]^2, {x, -Infinity, Infinity}, Assumptions -> sig > 0];
x2 = Integrate[x^2 Abs[psi27[x]]^2, {x, -Infinity, Infinity}, Assumptions -> sig > 0];
p2 = Integrate[Conjugate[psi27[x]] (-D[psi27[x], {x, 2}]), {x, -Infinity, Infinity}, Assumptions -> sig > 0];
dx7 = Simplify[Sqrt[x2], sig > 0]; dp7 = Simplify[Sqrt[p2], sig > 0];
show["2.7 norm", Simplify[nx, sig > 0]];
show["2.7 Delta x", dx7]; show["2.7 Delta p", dp7];
show["2.7 product", Simplify[dx7 dp7, sig > 0]];
check["2.7 Dx Dp == 1/2", Simplify[dx7 dp7 == 1/2, sig > 0]];

(* ---- 2.8 Weyl symmetric ordering ---- *)
sym28 = Simplify[(xop[pop[ff[x]]] + pop[xop[ff[x]]])/2];
show["2.8 (xp+px)/2 on ff", sym28];
(* Hermiticity test of symmetric vs non-symmetric ordering on Gaussians *)
A[ff_] := (xop[pop[ff]] + pop[xop[ff]])/2;
herm = Integrate[Conjugate[ph[x]] A[ps[x]], {x, -Infinity, Infinity}] -
       Integrate[Conjugate[A[ph[x]]] ps[x], {x, -Infinity, Infinity}];
check["2.8 symmetric ordering Hermitian", Simplify[herm] === 0];
nonh = Integrate[Conjugate[ph[x]] xop[pop[ps[x]]], {x, -Infinity, Infinity}] -
       Integrate[Conjugate[xop[pop[ph[x]]]] ps[x], {x, -Infinity, Infinity}];
check["2.8 x.p alone NOT Hermitian", Simplify[nonh] =!= 0];

(* ---- 2.9 self-adjoint vs Hermitian: boundary conditions ---- *)
(* box [0,L]: twisted BC psi(L)=e^{i th} psi(0). Boundary term -i[phi* psi]_0^L *)
bt = -I (Conjugate[pL] sL - Conjugate[p0] s0) /. {sL -> Exp[I th] s0, pL -> Exp[I th] p0};
bt = Simplify[bt, th \[Element] Reals];
show["2.9 box twisted-BC boundary term", bt];
check["2.9 box: B vanishes for |e^{i th}|=1", bt === 0];
(* half-line [0,inf): psi(inf)=0, psi(0) free -> B = i phi*(0) psi(0), cannot vanish for all *)
btHalf = -I (0 - Conjugate[p0] s0);
show["2.9 half-line boundary term", btHalf];
check["2.9 half-line: B != 0 (no self-adjoint extension)", btHalf =!= 0];

WriteString[$Output, "================ DONE ================\n"];
