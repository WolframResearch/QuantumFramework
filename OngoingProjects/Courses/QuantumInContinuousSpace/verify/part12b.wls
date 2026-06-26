#!/usr/bin/env wolframscript
(* Verify the NON-Gaussian states for the revised Parts 1-2. Natural units hbar=m=1. *)
show[l_, e_] := WriteString[$Output, l, " => ", ToString[e, InputForm], "\n"];
check[l_, b_] := WriteString[$Output, l, " : ", If[TrueQ[b], "PASS", "FAIL <<<<<"], "\n"];

WriteString[$Output, "======== CUSP STATE  psi = Sqrt[k] Exp[-k|x|]  (delta-potential bound state) ========\n"];
psiC[x_] := Sqrt[k] Exp[-k Abs[x]];
normC = Integrate[Abs[psiC[x]]^2, {x, -Infinity, Infinity}, Assumptions -> k > 0];
check["cusp normalized", Simplify[normC == 1, k > 0]];
densC = Simplify[Abs[psiC[x]]^2, k > 0 && x \[Element] Reals];
show["cusp density", densC];
(* 1.3-style localization probability P(|x|<a) -- use cusp here too as alt; main 1.3 uses sech below *)
PC = Integrate[densC, {x, -a, a}, Assumptions -> k > 0 && a > 0];
show["cusp P(|x|<a)", Simplify[PC, k > 0 && a > 0]];
(* moving cusp: <p>, <p^2> via Integral of |psi'|^2 (boundary-term-free), Delta p *)
psiCk[x_] := Sqrt[k] Exp[-k Abs[x]] Exp[I kk x];
mpC = Integrate[Conjugate[psiCk[x]] (-I D[psiCk[x], x]), {x, -Infinity, Infinity}, Assumptions -> k > 0 && kk \[Element] Reals];
p2C = Integrate[Abs[D[psiCk[x], x]]^2, {x, -Infinity, Infinity}, Assumptions -> k > 0 && kk \[Element] Reals];
show["moving cusp <p>", Simplify[mpC, k > 0 && kk \[Element] Reals]];
show["moving cusp <p^2>", Simplify[p2C, k > 0 && kk \[Element] Reals]];
show["moving cusp Delta p", Simplify[Sqrt[p2C - mpC^2], k > 0 && kk \[Element] Reals]];
(* numeric grid (k=1): resolve the cusp *)
gC = Subdivide[-12., 12., 6000]; dxC = gC[[2]] - gC[[1]];
sumC = Total[Abs[(psiC[#] /. k -> 1.) & /@ gC]^2] dxC;
show["cusp grid sum (k=1)", sumC];
check["cusp grid ~1", Abs[sumC - 1] < 10^-5];
(* cusp momentum-space wavefunction = Lorentzian (the 1D-hydrogen momentum distribution) *)
psiCt = FullSimplify[FourierTransform[psiC[x], x, p, FourierParameters -> {0, -1}], k > 0];
show["cusp psi-tilde(p) (Lorentzian)", psiCt];
normCp = Integrate[Abs[psiCt]^2, {p, -Infinity, Infinity}, Assumptions -> k > 0];
p2Cp = Integrate[p^2 Abs[psiCt]^2, {p, -Infinity, Infinity}, Assumptions -> k > 0];
check["cusp momentum norm 1", Simplify[normCp == 1, k > 0]];
show["cusp <p^2> in momentum space", Simplify[p2Cp, k > 0]];
check["cusp <p^2>=k^2 both ways", Simplify[p2Cp == k^2, k > 0]];

WriteString[$Output, "======== SECH SOLITON  psi = sech(x)/Sqrt[2]  (NLS soliton / Poschl-Teller ground) ========\n"];
psiS[x_] := Sech[x]/Sqrt[2];
normS = Integrate[Abs[psiS[x]]^2, {x, -Infinity, Infinity}];
check["sech normalized (int sech^2 = 2)", normS == 1];
show["int sech^2", Integrate[Sech[x]^2, {x, -Infinity, Infinity}]];
PS = Integrate[Abs[psiS[x]]^2, {x, -a, a}, Assumptions -> a > 0];
show["sech P(|x|<a)", Simplify[PS, a > 0]];
check["sech P=Tanh[a]", Simplify[PS == Tanh[a], a > 0]];
(* displaced sech: <x>, Delta x *)
psiSd[x_] := Sech[x - x0]/Sqrt[2];
mxS = Integrate[x Abs[psiSd[x]]^2, {x, -Infinity, Infinity}, Assumptions -> x0 \[Element] Reals];
mx2S = Integrate[x^2 Abs[psiSd[x]]^2, {x, -Infinity, Infinity}, Assumptions -> x0 \[Element] Reals];
show["sech <x>", Simplify[mxS, x0 \[Element] Reals]];
show["sech Delta x", Simplify[Sqrt[mx2S - mxS^2], x0 \[Element] Reals]];
check["sech <x>=x0, Dx=Pi/(2Sqrt3)", Simplify[mxS == x0, x0 \[Element] Reals] && Simplify[Sqrt[mx2S - mxS^2] == Pi/(2 Sqrt[3]), x0 \[Element] Reals]];
(* sech momentum space: FT of sech is sech (self-similar) *)
psiSt = FullSimplify[FourierTransform[psiS[x], x, p, FourierParameters -> {0, -1}]];
show["sech psi-tilde(p) (a sech!)", psiSt];
normSp = Integrate[Abs[psiSt]^2, {p, -Infinity, Infinity}];
p2Sp = Integrate[p^2 Abs[psiSt]^2, {p, -Infinity, Infinity}];
check["sech momentum norm 1", normSp == 1];
show["sech <p^2> (momentum space)", p2Sp];
(* sech <p^2> via |psi'|^2 *)
p2Sx = Integrate[Abs[D[psiS[x], x]]^2, {x, -Infinity, Infinity}];
show["sech <p^2> (position space, int|psi'|^2)", p2Sx];
check["sech <p^2>=1/3 both ways", p2Sp == 1/3 && p2Sx == 1/3];
(* uncertainty product for sech: Pi/6, contrast with Gaussian 1/2 *)
dxS = Pi/(2 Sqrt[3]); dpS = Sqrt[1/3];
show["sech Dx Dp", Simplify[dxS dpS]];
check["sech Dx Dp = Pi/6 > 1/2", Simplify[dxS dpS == Pi/6] && (N[Pi/6] > 0.5)];

WriteString[$Output, "======== HALF-LINE HYDROGEN 1s RADIAL  u = 2 k^{3/2} x Exp[-k x], x>0 ========\n"];
u[x_] := 2 k^(3/2) x Exp[-k x];
normU = Integrate[Abs[u[x]]^2, {x, 0, Infinity}, Assumptions -> k > 0];
check["radial u normalized", Simplify[normU == 1, k > 0]];
(* recoil overlap: <u | e^{i q x} u> = sudden-momentum-kick form factor (COMPLEX) *)
ovl = Integrate[Conjugate[u[x]] Exp[I q x] u[x], {x, 0, Infinity}, Assumptions -> k > 0 && q \[Element] Reals];
ovl = Simplify[ovl, k > 0 && q \[Element] Reals];
show["radial recoil overlap <u|e^{iqx}|u>", ovl];
show["  |overlap|", Simplify[Sqrt[ovl Conjugate[ovl]], k > 0 && q \[Element] Reals]];
show["  |overlap|^2", Simplify[ovl Conjugate[ovl], k > 0 && q \[Element] Reals]];
check["overlap = 8k^3/(2k - i q)^3", Simplify[ovl == 8 k^3/(2 k - I q)^3, k > 0 && q \[Element] Reals]];

WriteString[$Output, "======== GAUSSIAN (2.7 saturating reference only) ========\n"];
psiG[x_] := Pi^(-1/4) sig^(-1/2) Exp[-x^2/(2 sig^2)];
dxG = Sqrt[Integrate[x^2 Abs[psiG[x]]^2, {x, -Infinity, Infinity}, Assumptions -> sig > 0]];
dpG = Sqrt[Integrate[Abs[D[psiG[x], x]]^2, {x, -Infinity, Infinity}, Assumptions -> sig > 0]];
check["Gaussian Dx Dp = 1/2 (saturates)", Simplify[dxG dpG == 1/2, sig > 0]];

WriteString[$Output, "======== CHIRPED SOLITON  psi = sech(x) Exp[i beta x^2]/Sqrt2  (position-dependent velocity) ========\n"];
psiCh[x_] := Sech[x] Exp[I beta x^2]/Sqrt[2];
psiChc[x_] := Sech[x] Exp[-I beta x^2]/Sqrt[2];
jCh = Simplify[(psiChc[x] D[psiCh[x], x] - psiCh[x] D[psiChc[x], x])/(2 I), Element[{beta, x}, Reals]];
vCh = Simplify[jCh/(psiChc[x] psiCh[x]), Element[{beta, x}, Reals]];
show["chirped current j", jCh];
show["chirped velocity v(x) = local phase slope", vCh];
check["chirped v(x) = 2 beta x (position-dependent!)", Simplify[vCh == 2 beta x, Element[{beta, x}, Reals]]];

WriteString[$Output, "======== DONE ========\n"];
