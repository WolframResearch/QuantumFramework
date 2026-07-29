#!/usr/bin/env wolframscript
(* Verification harness for Solutions-Part-03.md. hbar = m = 1.
   Every stated result asserted. NO Quiet: expected messages are DECLARED in the
   VerificationTest, so the test checks that they fire. *)

report[vt_] := WriteString[$Output,
   If[vt["Outcome"] === "Success", "PASS  ", "FAIL  "], vt["TestID"], "\n"];

(* ================= 3.1 Poschl-Teller FAMILY: bound tower vs scattering ================= *)
Vpt[x_, \[Lambda]_] := -\[Lambda] (\[Lambda] + 1) Sech[x]^2/2;
asl = Assumptions -> \[Lambda] > 0 && Element[x, Reals];
(* the GENERAL SYMBOLIC statement: arbitrary lambda, not one well *)
report@VerificationTest[
  FullSimplify[-(1/2) D[Sech[x]^\[Lambda], {x, 2}] + Vpt[x, \[Lambda]] Sech[x]^\[Lambda]
    - (-\[Lambda]^2/2) Sech[x]^\[Lambda], asl], 0,
  TestID -> "3.1 Sech^lam is the exact ground state at E0=-lam^2/2 for ARBITRARY lam"];
report@VerificationTest[
  FullSimplify[Integrate[Sech[x]^(2 \[Lambda]), {x, -Infinity, Infinity}, asl], asl],
  Sqrt[Pi] Gamma[\[Lambda]]/Gamma[1/2 + \[Lambda]],
  TestID -> "3.1 squared norm is a closed form in lam"];
(* and it must reduce to Part 1's hand computation at lam = 1 *)
report@VerificationTest[
  FullSimplify[Sqrt[Pi] Gamma[\[Lambda]]/Gamma[1/2 + \[Lambda]] /. \[Lambda] -> 1], 2,
  TestID -> "3.1 lam=1 recovers Part 1's Int Sech^2 = 2"];
(* THE TOWER as a theorem: symbolic in lam AND n at once, not a list of integer instances *)
H[\[Lambda]_][f_] := -(1/2) D[f, {x, 2}] + Vpt[x, \[Lambda]] f;
report@VerificationTest[
  FullSimplify[H[\[Lambda]][LegendreP[\[Lambda], \[Lambda] - n, Tanh[x]]]
    - (-(\[Lambda] - n)^2/2) LegendreP[\[Lambda], \[Lambda] - n, Tanh[x]], asl], 0,
  TestID -> "3.1 FULL tower: E_n = -(lam-n)^2/2 for ARBITRARY lam AND ARBITRARY n"];
report@VerificationTest[
  FullSimplify[H[\[Lambda]][Sech[x]^(\[Lambda] - 1) Tanh[x]]
    - (-(\[Lambda] - 1)^2/2) Sech[x]^(\[Lambda] - 1) Tanh[x],
   Assumptions -> \[Lambda] > 1 && Element[x, Reals]], 0,
  TestID -> "3.1 first excited state Sech^(lam-1) Tanh, symbolic in lam"];
(* the elementary instances: lam=2 is where a SECOND state first appears *)
report@VerificationTest[
  {FullSimplify[H[2][Sech[x]^2] - (-2) Sech[x]^2],
   FullSimplify[H[2][Sech[x] Tanh[x]] - (-1/2) Sech[x] Tanh[x]]}, {0, 0},
  TestID -> "3.1 lam=2 tower: E0=-2 AND E1=-1/2"];
(* nothing is special about integer lam for the TOWER (only for transparency) *)
report@VerificationTest[
  FullSimplify[H[5/2][Sech[x]^(3/2) Tanh[x]] - (-9/8) Sech[x]^(3/2) Tanh[x]], 0,
  TestID -> "3.1 NON-integer lam=5/2 carries its tower too"];
(* SHAPE INVARIANCE as the real theorem, not the tautology E_1(2) == E_0(1).
   A[lam] and Ad[lam] factorize H, and reversing the product gives the SAME family at lam-1,
   which FORCES the ladder rather than noticing it. *)
A[\[Lambda]_][f_] := (1/Sqrt[2]) (D[f, x] + \[Lambda] Tanh[x] f);
Ad[\[Lambda]_][f_] := (1/Sqrt[2]) (-D[f, x] + \[Lambda] Tanh[x] f);
report@VerificationTest[
  FullSimplify[Ad[\[Lambda]][A[\[Lambda]][f[x]]] - (H[\[Lambda]][f[x]] + \[Lambda]^2/2 f[x])], 0,
  TestID -> "3.1 factorization H_lam = Ad A - lam^2/2 (so no state lies below -lam^2/2)"];
report@VerificationTest[
  FullSimplify[A[\[Lambda]][Ad[\[Lambda]][f[x]]] - (H[\[Lambda] - 1][f[x]] + \[Lambda]^2/2 f[x])], 0,
  TestID -> "3.1 SHAPE INVARIANCE: the partner of V_lam IS V_{lam-1}"];
report@VerificationTest[
  FullSimplify[A[\[Lambda]][Sech[x]^\[Lambda]], Assumptions -> \[Lambda] > 0], 0,
  TestID -> "3.1 A annihilates the ground state (so -lam^2/2 is attained, not just bounded)"];

\[Psi]k[x_] := Exp[I k x] (Tanh[x] - I k);
report@VerificationTest[
  FullSimplify[-(1/2) D[\[Psi]k[x], {x, 2}] + Vpt[x, 1] \[Psi]k[x] - (k^2/2) \[Psi]k[x],
   Assumptions -> Element[{k, x}, Reals]], 0,
  TestID -> "3.1 scattering state solves same eqn at E=k^2/2, arbitrary k"];
dens = FullSimplify[Abs[\[Psi]k[x]]^2, Assumptions -> Element[{k, x}, Reals]];
report@VerificationTest[dens, k^2 + Tanh[x]^2, TestID -> "3.1 scattering density"];
(* REFLECTIONLESS with content: r = 0, not merely equal asymptotic density (necessary,
   not sufficient). Stripping the carrier leaves a CONSTANT at each end: no e^{-2ikx} piece. *)
report@VerificationTest[
  {Limit[\[Psi]k[x] Exp[-I k x], x -> -Infinity, Assumptions -> Element[k, Reals]],
   Limit[\[Psi]k[x] Exp[-I k x], x -> Infinity, Assumptions -> Element[k, Reals]]},
  {-1 - I k, 1 - I k},
  TestID -> "3.1 r=0: psi_k e^{-ikx} tends to a constant at BOTH ends"];
report@VerificationTest[
  FullSimplify[Abs[(1 - I k)/(-1 - I k)]^2, Element[k, Reals]], 1,
  TestID -> "3.1 |t|^2 = 1 for every k: perfectly transparent"];
report@VerificationTest[
  {Limit[dens, x -> Infinity], Limit[dens, x -> -Infinity]}, {1 + k^2, 1 + k^2},
  TestID -> "3.1 equal asymptotic density (necessary but NOT sufficient for r=0)"];
(* CONTROL for the r=0 argument: the limit having a value is the whole content. Give the state
   a reflected piece and the same limit does NOT exist. Without this the reader has only the
   author's word that "tends to a constant" means anything. *)
report@VerificationTest[Limit[Exp[-2 I k x], x -> -Infinity, Assumptions -> k > 0], Indeterminate,
  TestID -> "3.1 CONTROL: a reflected e^{-2ikx} piece makes the limit Indeterminate"];
(* LEVINSON: the phase shift counts the bound states, tying 3.1's two halves together *)
tk = FullSimplify[ComplexExpand[Arg[(1 - I k)/(-1 - I k)], TargetFunctions -> {Re, Im}], k > 0];
d0 = Limit[tk, k -> 0, Direction -> "FromAbove"]; dInf = Limit[tk, k -> Infinity];
report@VerificationTest[{d0, dInf, (d0 - dInf)/Pi}, {Pi, 0, 1},
  TestID -> "3.1 LEVINSON: (delta(0)-delta(inf))/Pi = 1 = the bound-state count at lam=1"];
(* the direction is physics, not pedantry: k>0 is a wavenumber, and two-sided the threshold
   sits on a branch cut. Assert that the two-sided limit genuinely fails. *)
report@VerificationTest[Limit[tk, k -> 0], Indeterminate,
  TestID -> "3.1 two-sided limit at threshold is Indeterminate (branch cut at k=0)"];
(* WHY transparency is quantized: R vanishes identically in k exactly at integer lam.
   R MUST be tied to the equation, not asserted next to Sin[Pi lam]: a mutant with
   Sinh -> Cosh is physically wrong yet passes every test that only inspects Sin[Pi lam].
   Unitarity R = 1 - |t|^2 against the family's transmission amplitude is the link. *)
tAmp[\[Lambda]_, k_] := Gamma[1 + \[Lambda] - I k] Gamma[-\[Lambda] - I k]/(Gamma[1 - I k] Gamma[-I k]);
refl[\[Lambda]_, k_] := Sin[Pi \[Lambda]]^2/(Sin[Pi \[Lambda]]^2 + Sinh[Pi k]^2);
report@VerificationTest[
  FullSimplify[tAmp[1, k] - (1 - I k)/(-1 - I k), Assumptions -> k > 0], 0,
  TestID -> "3.1 the family's t reduces to the lam=1 amplitude derived from psi_k"];
report@VerificationTest[
  FullSimplify[1 - Abs[tAmp[\[Lambda], k]]^2 - refl[\[Lambda], k],
   Assumptions -> k > 0 && \[Lambda] > 0], 0,
  TestID -> "3.1 UNITARITY ties R to the amplitude: 1 - |t|^2 - R = 0 for arbitrary lam, k"];
(* the mutation that used to pass must now fail *)
report@VerificationTest[
  FullSimplify[1 - Abs[tAmp[\[Lambda], k]]^2 - Sin[Pi \[Lambda]]^2/(Sin[Pi \[Lambda]]^2 + Cosh[Pi k]^2),
   Assumptions -> k > 0 && \[Lambda] > 0] === 0, False,
  TestID -> "3.1 MUTANT Sinh->Cosh is REJECTED (the old test could not tell)"];
report@VerificationTest[FullSimplify[refl[\[Lambda], k] /. \[Lambda] -> 2, Assumptions -> k > 0], 0,
  TestID -> "3.1 R = 0 identically in k at integer lam"];
report@VerificationTest[
  {N@Table[refl[\[Lambda], 1], {\[Lambda], {1, 2, 3}}],
   AllTrue[N@Table[refl[\[Lambda], 1], {\[Lambda], {1/2, 3/2, 13/10}}], # > 10.^-3 &]},
  {{0., 0., 0.}, True},
  TestID -> "3.1 transparency is QUANTIZED: zero at integer lam, nonzero at half-integer"];
(* LEVINSON's half-unit: the lam=1 integer count is the EXCEPTION, not the rule *)
report@VerificationTest[
  {Round[N[Arg[tAmp[13/10, 10.^-8]]/Pi], 0.001], Ceiling[13/10]}, {-0.5, 2},
  TestID -> "3.1 at generic lam, delta(0)/Pi = n_b - 1/2 (mod 2), NOT n_b: the 1D half-unit"];
report@VerificationTest[
  {Abs[tAmp[13/10, 10.^-8]]^2 < 10.^-10, Abs[Abs[tAmp[1, 10.^-8]]^2 - 1] < 10.^-10}, {True, True},
  TestID -> "3.1 generic lam reflects TOTALLY at threshold; integer lam is a zero-energy resonance"];
(* the half-bound state: the n=lam member of the tower, at E=0, bounded, not normalizable *)
report@VerificationTest[
  {FullSimplify[LegendreP[1, 0, Tanh[x]]],
   FullSimplify[-(1/2) D[Tanh[x], {x, 2}] + Vpt[x, 1] Tanh[x]]}, {Tanh[x], 0},
  TestID -> "3.1 the half-bound state IS Tanh[x] = the n=lam tower member, solving H psi = 0"];
pn = Integrate[dens, {x, -a, a}, Assumptions -> a > 0 && Element[k, Reals]];
report@VerificationTest[FullSimplify[pn - 2 (a + a k^2 - Tanh[a])], 0,
  TestID -> "3.1 partial norm = 2(a + a k^2 - Tanh a)"];
report@VerificationTest[Limit[pn, a -> Infinity, Assumptions -> Element[k, Reals]], Infinity,
  TestID -> "3.1 scattering norm diverges (not in L^2)"];

(* ================= 3.2 infinite well ================= *)
{es, fs} = DEigensystem[{-(1/2) D[u[x], {x, 2}], DirichletCondition[u[x] == 0, True]},
   u[x], {x, 0, L}, 4, Assumptions -> L > 0, Method -> "Normalize"];
report@VerificationTest[FullSimplify[es - Table[n^2 Pi^2/(2 L^2), {n, 1, 4}], L > 0],
  {0, 0, 0, 0}, TestID -> "3.2 DEigensystem spectrum = n^2 Pi^2/(2L^2)"];
report@VerificationTest[
  FullSimplify[Integrate[#^2, {x, 0, L}, Assumptions -> L > 0], L > 0] & /@ fs, {1, 1, 1, 1},
  TestID -> "3.2 DEigensystem eigenfunctions normalized"];
(* THE TRAP: no boundary condition returns the NEUMANN spectrum, incl. a spurious E=0 *)
report@VerificationTest[
  FullSimplify[DEigenvalues[-(1/2) D[u[x], {x, 2}], u[x], {x, 0, L}, 4, Assumptions -> L > 0]
    - Table[n^2 Pi^2/(2 L^2), {n, 0, 3}], L > 0], {0, 0, 0, 0},
  TestID -> "3.2 TRAP: no BC gives Neumann spectrum with spurious E=0"];

\[Psi]n[n_, x_] := Sqrt[2/L] Sin[n Pi x/L];
asn = Assumptions -> L > 0 && Element[n, Integers] && n >= 1;
report@VerificationTest[
  FullSimplify[-(1/2) D[\[Psi]n[n, x], {x, 2}] - (n^2 Pi^2/(2 L^2)) \[Psi]n[n, x], asn], 0,
  TestID -> "3.2 eigenvalue equation for ARBITRARY n"];
report@VerificationTest[FullSimplify[Integrate[\[Psi]n[n, x]^2, {x, 0, L}, asn], asn], 1,
  TestID -> "3.2 prefactor normalized for ARBITRARY n"];
(* why the arbitrary-n identity is needed at all: a symbolic COUNT is refused *)
report@VerificationTest[
  Head@DEigenvalues[{-(1/2) D[u[x], {x, 2}], DirichletCondition[u[x] == 0, True]},
    u[x], {x, 0, L}, n, Assumptions -> L > 0 && n > 0], DEigenvalues,
  TestID -> "3.2 DEigenvalues refuses a symbolic count, so 4 levels are not a theorem"];
(* the Neumann trap is SILENT: that is the whole danger, so assert the silence.
   Check[expr, "fired"] returns "fired" iff ANY message was generated. *)
report@VerificationTest[
  Check[DEigenvalues[-(1/2) D[u[x], {x, 2}], u[x], {x, 0, L}, 4, Assumptions -> L > 0],
   "MESSAGE-FIRED"] =!= "MESSAGE-FIRED", True,
  TestID -> "3.2 TRAP is SILENT: the Neumann substitution emits no message"];

(* ================= 3.3 orthonormality and completeness ================= *)
report@VerificationTest[
  FullSimplify[Integrate[\[Psi]n[m, x] \[Psi]n[n, x], {x, 0, L},
     Assumptions -> L > 0 && Element[{m, n}, Integers] && m >= 1 && n >= 1 && m != n],
   Assumptions -> Element[{m, n}, Integers] && m != n], 0,
  TestID -> "3.3 <m|n> = 0 for arbitrary m != n"];
(* f MUST be asymmetric. With the symmetric x(L-x) every even c_n vanishes, so Parseval
   cannot detect a basis missing every even mode: assert that failure, then avoid it. *)
cSym = FullSimplify[Integrate[\[Psi]n[n, x] x (L - x), {x, 0, L}, asn], asn];
report@VerificationTest[
  {FullSimplify[cSym /. n -> 2, L > 0],
   FullSimplify[Sum[cSym^2, {n, 1, Infinity}], Assumptions -> L > 0],
   FullSimplify[Sum[(cSym /. n -> 2 j - 1)^2, {j, 1, Infinity}], Assumptions -> L > 0]},
  {0, L^5/30, L^5/30},
  TestID -> "3.3 TRAP: with symmetric f, Parseval CANNOT see every even mode deleted"];
(* fC, NOT f: 3.3's completeness target must not shadow the generic test function f that
   3.1/3.7/3.9 use. With f[x_]:=x^2(L-x) leaking, 3.9's Integrate[delta f]=f[0]=0 and 3.7's
   commutator go concrete, and both tests then pass VACUOUSLY (0===0). *)
fC[x_] := x^2 (L - x);
cn = FullSimplify[Integrate[\[Psi]n[n, x] fC[x], {x, 0, L}, asn], asn];
report@VerificationTest[FreeQ[FullSimplify[cn /. n -> 2, L > 0], 0] || (FullSimplify[cn /. n -> 2, L > 0] =!= 0),
  True, TestID -> "3.3 the asymmetric fC has NONZERO even coefficients"];
(* Parseval closes EXACTLY and symbolically, AND now discriminates *)
report@VerificationTest[
  {FullSimplify[Sum[cn^2, {n, 1, Infinity}], Assumptions -> L > 0],
   Integrate[fC[x]^2, {x, 0, L}, Assumptions -> L > 0]}, {L^7/105, L^7/105},
  TestID -> "3.3 Parseval: sum |c_n|^2 = int f^2 = L^7/105"];
report@VerificationTest[
  FullSimplify[Sum[(cn /. n -> 2 j - 1)^2, {j, 1, Infinity}], Assumptions -> L > 0], L^7/120,
  TestID -> "3.3 deleting the even modes now FALLS SHORT (L^7/120 < L^7/105): the test discriminates"];
(* the prose claims dropping m != n returns a BARE 0, falsely asserting <n|n> = 0.
   Assert the trap, or the reader has only the author's word for it. *)
report@VerificationTest[
  FullSimplify[Integrate[\[Psi]n[m, x] \[Psi]n[n, x], {x, 0, L},
     Assumptions -> L > 0 && Element[{m, n}, Integers] && m >= 1 && n >= 1],
   Assumptions -> Element[{m, n}, Integers] && m >= 1 && n >= 1], 0,
  TestID -> "3.3 TRAP: without m != n the kernel returns a bare 0, contradicting the norm"];
(* the resummation does NOT close symbolically -> documented rung 4, several random draws *)
res = Sum[cn \[Psi]n[n, x], {n, 1, Infinity}];
report@VerificationTest[
  FullSimplify[res - fC[x], Assumptions -> L > 0 && 0 < x < L] === 0,
  False, TestID -> "3.3 the resummation does NOT close under FullSimplify (PolyLog form)"];
report@VerificationTest[
  Union@Table[With[{lv = RandomReal[{0.5, 4}]},
     With[{xv = RandomReal[{0.01, 0.99}] lv},
      Chop[N[(res - fC[x]) /. {L -> lv, x -> xv}], 10^-10]]], {8}], {0},
  TestID -> "3.3 resummation reproduces x^2(L-x), 8 random (L,x) draws"];

(* ================= 3.4 finite well ================= *)
fEven[u0_][u_?NumericQ] := u Sin[u] - Sqrt[u0^2 - u^2] Cos[u];
fOdd[u0_][u_?NumericQ] := u Cos[u] + Sqrt[u0^2 - u^2] Sin[u];
(* NOTE: `fs` here shadows the global `fs` bound by 3.2 above. That is deliberate and it is
   the shipped bytes: With is HoldAll and lexical, so the inner binding wins and the global is
   untouched. Renaming it in the harness would verify a transcription, not the document. *)
bracket[f_, u0_] := With[{us = Subdivide[10.^-9, u0 (1 - 10.^-9), 4000]},
   With[{fs = f /@ us}, Pick[Partition[us, 2, 1], Sign[Most[fs] Rest[fs]], -1]]];
rootsOf[f_, u0_] := (u /. FindRoot[f[u], {u, #[[1]], #[[2]]}, Method -> "Brent"]) & /@ bracket[f, u0];
uAll[u0_] := Sort[Join[rootsOf[fEven[u0], u0], rootsOf[fOdd[u0], u0]]];
nBound[u0_] := Ceiling[2 u0/Pi];
report@VerificationTest[Length[uAll[5.]], nBound[5.],
  TestID -> "3.4 root count matches closed form at u0=5"];
report@VerificationTest[
  Table[Length[uAll[u0]] == nBound[u0], {u0, {0.05, 0.5, 1.5, 2.0, 3.2, 7.9, 12.8, 20.1, 50.3}}],
  ConstantArray[True, 9], TestID -> "3.4 count formula holds across u0 = 0.05 .. 50.3"];
Esym[u0_, u_, a_] := -(u0^2 - u^2)/(2 a^2);
eFinite = Esym[5., uAll[5.], 1.];
report@VerificationTest[And @@ Thread[-12.5 < eFinite < 0.], True,
  TestID -> "3.4 every bound level lies strictly between -V0 and 0"];
report@VerificationTest[
  With[{us = uAll[200.]}, Max@Abs[(us[[1 ;; 4]]^2/us[[1]]^2) - {1., 4., 9., 16.}] < 10.^-4], True,
  TestID -> "3.4 deep-well limit recovers the infinite-well ratios {1,4,9,16}"];
report@VerificationTest[Length[uAll[0.01]], 1,
  TestID -> "3.4 EDGE: however shallow, a 1D well binds exactly one state"];
(* the deep-well ratios are claimed to 5 significant figures, so assert THAT, not a
   loose absolute band that would pass for a wrong claim too *)
report@VerificationTest[
  With[{us = uAll[200.]}, Max@Abs[(us[[1 ;; 4]]^2/us[[1]]^2)/{1., 4., 9., 16.} - 1.] < 10.^-5],
  True, TestID -> "3.4 deep-well ratios agree to 5 significant figures (relative)"];
(* GATE (r): DEigensystem declines INSTANTLY on the piecewise well. It does not hang.
   The Piecewise is written out in full rather than via well[x]: `well` is not defined until
   the 3.5 block below, so referring to it here would hand DEigenvalues an UNDEFINED symbol,
   which it also declines. The test would then pass while testing nothing, which is exactly
   the failure the prose it guards once committed. Positive control below. *)
report@VerificationTest[
  Head@DEigenvalues[{-(1/2) D[u[x], {x, 2}] + Piecewise[{{-12.5, Abs[x] < 1.}}, 0.] u[x],
     DirichletCondition[u[x] == 0, True]}, u[x], {x, -8., 8.}, 4], DEigenvalues,
  TestID -> "3.4 GATE: DEigenvalues returns UNEVALUATED on the piecewise finite well"];
report@VerificationTest[
  First@AbsoluteTiming[
    DEigenvalues[{-(1/2) D[u[x], {x, 2}] + Piecewise[{{-12.5, Abs[x] < 1.}}, 0.] u[x],
      DirichletCondition[u[x] == 0, True]}, u[x], {x, -8., 8.}, 4]] < 30., True,
  TestID -> "3.4 GATE: it declines fast (loud failure), it does not grind"];
(* The gate above CANNOT distinguish the real well from any unknown symbol, so on its own it is
   weak evidence: it is consistent with "DEigenvalues refuses anything it does not recognise". *)
report@VerificationTest[
  Head@DEigenvalues[{-(1/2) D[u[x], {x, 2}] + qqqUndefined[x] u[x],
     DirichletCondition[u[x] == 0, True]}, u[x], {x, -8., 8.}, 4], DEigenvalues,
  TestID -> "3.4 GATE control: an undefined potential is refused identically, so the gate alone proves little"];
(* The POSITIVE control is what gives the gate teeth, and it must carry a NONZERO potential:
   3.2's infinite well has no potential term at all, so it can only exclude "the solver is
   useless", not the live hypothesis "the solver refuses ANY potential". A constant potential
   closes exactly, so refusal is specific to the PIECEWISE form and not to potentials as such. *)
report@VerificationTest[
  DEigenvalues[{-(1/2) D[u[x], {x, 2}] + 1 u[x], DirichletCondition[u[x] == 0, True]},
   u[x], {x, 0, 1}, 2], {(2 + Pi^2)/2, (2 + 4 Pi^2)/2},
  TestID -> "3.4 GATE positive control: a NONZERO (constant) potential closes exactly"];
(* GATE (B): with depth and half-width SYMBOLIC it does not decline, it GRINDS.
   Bounded with TimeConstrained and asserted as $Aborted: the prose preaches bounding the
   probe, so the harness must bound it too rather than take the author's word. *)
report@VerificationTest[
  TimeConstrained[
   DEigenvalues[{-(1/2) D[u[x], {x, 2}] + Piecewise[{{-V0, Abs[x] < a}}, 0] u[x],
     DirichletCondition[u[x] == 0, True]}, u[x], {x, -8, 8}, 4,
    Assumptions -> V0 > 0 && 0 < a < 8], 40], $Aborted,
  TestID -> "3.4 GATE: symbolic depth/half-width GRINDS (aborted at 40 s), it does not decline"];

(* ================= 3.5 NDEigensystem vs the analytic roots ================= *)
well[x_] := Piecewise[{{-12.5, Abs[x] < 1.}}, 0.];
ndAll[box_, mcm_] := Sort@NDEigenvalues[{-(1/2) D[u[x], {x, 2}] + well[x] u[x],
    DirichletCondition[u[x] == 0., True]}, u[x], {x, -box, box}, 4,
   Method -> {"SpatialDiscretization" -> {"FiniteElement",
      {"MeshOptions" -> {MaxCellMeasure -> mcm}}},
     "Eigensystem" -> {"Arnoldi", "Shift" -> -12.5}}];
(* assert the TRAP is real: the default misses the bound spectrum entirely *)
report@VerificationTest[
  Max@NDEigenvalues[{-(1/2) D[u[x], {x, 2}] + well[x] u[x],
      DirichletCondition[u[x] == 0., True]}, u[x], {x, -8., 8.}, 4] > -1.,
  True, TestID -> "3.5 TRAP: unshifted NDEigenvalues misses every bound state"];
(* the trap's danger is the SILENCE, not the wrong numbers: assert it separately *)
report@VerificationTest[
  Check[NDEigenvalues[{-(1/2) D[u[x], {x, 2}] + well[x] u[x],
     DirichletCondition[u[x] == 0., True]}, u[x], {x, -8., 8.}, 4],
   "MESSAGE-FIRED"] =!= "MESSAGE-FIRED", True,
  TestID -> "3.5 TRAP is SILENT: no message accompanies the wrong spectrum"];
report@VerificationTest[
  Max@Abs[(ndAll[12., 0.005] - eFinite)/eFinite] < 10.^-6, True,
  TestID -> "3.5 shifted NDEigensystem agrees with the 3.4 transcendental roots"];
(* refinement sweep reported as a SEQUENCE: the error must fall monotonically.
   chnpdef fires on the coarse mesh; DECLARE it so the test asserts that it fires
   rather than merely tolerating it. *)
report@VerificationTest[
  With[{errs = (Max@Abs[(ndAll[8., #] - eFinite)/eFinite]) & /@ {0.5, 0.1, 0.02}},
   And @@ Thread[Rest[errs] < Most[errs]]], True, {Eigensystem::chnpdef},
  TestID -> "3.5 mesh refinement sweep: error falls monotonically"];
(* the weakest-bound state is the one that needs the biggest box *)
report@VerificationTest[
  With[{e4 = (Last@ndAll[#, 0.005]) & /@ {4., 8., 12.}},
   And @@ Thread[Abs[Rest[e4] - Last[eFinite]] < Abs[Most[e4] - Last[eFinite]]]], True,
  TestID -> "3.5 box sweep: shallowest state converges as the box grows"];

(* ================= 3.6 Numerov and shooting ================= *)
(* the shipped bytes: FoldList over a sliding window of {k^2, f} pairs. No index counter (the
   Python accent), and k^2 is kept in the window so the middle coefficient is computed as
   2(1 - 5 h^2 k^2/12) directly. Eliminating k^2 and writing it as 12 - 10 f is algebraically
   identical but subtracts two numbers near 12 and 10 to get a small one: that cancellation
   costs ~8x accuracy and makes the h^4 sweep DEGRADE at the fine end (measured). *)
numerov[En_?NumericQ, v_, {x0_, x1_}, npts_Integer] := With[{h = (x1 - x0)/(npts - 1)},
   With[{ks = 2 (En - v[x0 + # h]) & /@ Range[0, npts - 1]},
    With[{w = Partition[Transpose[{ks, 1 + h^2 ks/12}], 3, 1]},
     Prepend[Last /@ FoldList[
        With[{fm = #2[[1, 2]], k0 = #2[[2, 1]], fp = #2[[3, 2]]},
          {#1[[2]], (2 (1 - 5 h^2 k0/12) #1[[2]] - fm #1[[1]])/fp}] &,
        {0., 1.*^-6}, w], 0.]]]];
shoot[En_?NumericQ] := Last[numerov[En, 0. &, {0., 1.}, 2001]];
report@VerificationTest[
  Max@Abs[Table[(u /. FindRoot[shoot[u], {u, n^2 Pi^2/2. - 0.4, n^2 Pi^2/2. + 0.4},
        Method -> "Brent"]) - n^2 Pi^2/2., {n, 1, 4}]] < 10.^-8, True,
  TestID -> "3.6 Numerov reproduces the exact infinite-well spectrum"];
shootND[En_?NumericQ] := NDSolveValue[{u''[x] == -2 En u[x], u[0] == 0, u'[0] == 1}, u[1], {x, 0, 1}];
(* worst level is n=4 at ~2.2e-6, ~6100x (nearly four orders) WORSE than Numerov: that
   ordering is the point of running both, so the threshold records the measured value. *)
report@VerificationTest[
  Max@Abs[Table[(v /. FindRoot[shootND[v], {v, n^2 Pi^2/2. - 0.4, n^2 Pi^2/2. + 0.4},
        Method -> "Brent"]) - n^2 Pi^2/2., {n, 1, 4}]] < 10.^-5, True,
  TestID -> "3.6 NDSolve shooting agrees (independent machinery)"];
(* the prose says "four orders"; the measured gap is ~6100x = 10^3.8, so assert THAT band,
   not a loose >100x that would pass for a much weaker claim *)
report@VerificationTest[
  With[{eNum = Max@Abs[Table[(u /. FindRoot[shoot[u], {u, n^2 Pi^2/2. - 0.4, n^2 Pi^2/2. + 0.4},
         Method -> "Brent"]) - n^2 Pi^2/2., {n, 1, 4}]],
    eND = Max@Abs[Table[(v /. FindRoot[shootND[v], {v, n^2 Pi^2/2. - 0.4, n^2 Pi^2/2. + 0.4},
         Method -> "Brent"]) - n^2 Pi^2/2., {n, 1, 4}]]},
   10.^3.5 < eND/eNum < 10.^4.], True,
  TestID -> "3.6 Numerov beats general-purpose NDSolve by ~4 orders (h^4 tuned to this ODE)"];
(* h^4: the FIRST two halvings must cut the error by ~16. The prose claims the sequence
   16.0, 16.0, 14.0, 2.5, so assert the DECAY of the ratio too: the 3rd is already short of 16
   and the 4th has left the pattern entirely. A test that only bounded 12 < r < 20 would pass
   while the prose's roundoff-floor story was false. *)
shootN[En_?NumericQ, np_Integer] := Last[numerov[En, 0. &, {0., 1.}, np]];
report@VerificationTest[
  With[{rs = With[{errs = Table[Abs[(u /. FindRoot[shootN[u, np], {u, Pi^2/2. - 0.4, Pi^2/2. + 0.4},
           Method -> "Brent"]) - Pi^2/2.], {np, {51, 101, 201, 401, 801}}]},
      Most[errs]/Rest[errs]]},
   {15.5 < rs[[1]] < 16.5, 15.5 < rs[[2]] < 16.5, 13.5 < rs[[3]] < 14.5, 2. < rs[[4]] < 3.}],
  {True, True, True, True},
  TestID -> "3.6 Numerov h^4 then roundoff floor: ratios 16.0, 16.0, 14.0, 2.5"];
(* THE TRAP: Last of an unevaluated NumericQ-guarded call returns the last ARGUMENT *)
report@VerificationTest[Last[numerov[En, 0. &, {0., 1.}, 2001]], 2001,
  TestID -> "3.6 TRAP: symbolic En leaves Last returning npts, so FindRoot sees a constant"];
(* the returned list must pair with the grid: Prepend restores psi_0 = 0, the seed the fold
   consumes rather than emits. Without it every plotted point is off by one step. *)
report@VerificationTest[Length@numerov[4.9348, 0. &, {0., 1.}, 201], 201,
  TestID -> "3.6 numerov returns npts values (psi_0 restored), so it pairs with the grid"];
(* V = 0 makes k^2 CONSTANT, so the infinite well never exercises the varying-k^2 weights that
   are Numerov's entire reason to exist. Assert the degeneracy, then test on a real well. *)
report@VerificationTest[
  SameQ @@ (2 (4.9348 - 0.) & /@ Range[0, 5]), True,
  TestID -> "3.6 V=0 makes k^2 constant: the benchmark does NOT exercise the scheme"];
report@VerificationTest[
  SameQ @@ (2 (0. - (-3 Sech[#]^2 &)[-8. + # 0.004]) & /@ Range[0, 5]), False,
  TestID -> "3.6 the PT well makes k^2 vary: THIS exercises the scheme"];
shootPT[En_?NumericQ] := Last[numerov[En, -3 Sech[#]^2 &, {-8., 8.}, 4001]];
(* brmp fires on BOTH PT roots: Brent brackets to machine precision while the shot endpoint,
   which is an exponentially amplified tail value, still exceeds the absolute tolerance.
   Declared rather than suppressed: it is the scheme telling the truth about its own limits. *)
report@VerificationTest[
  Abs[(u /. FindRoot[shootPT[u], {u, -2.3, -1.7}, Method -> "Brent"]) - (-2)] < 10.^-9, True,
  {FindRoot::brmp},
  TestID -> "3.6 Numerov on the PT lam=2 well recovers the exact E0 = -2 (varying potential)"];
(* the weakly bound state is the hard one: kappa=1 has not decayed by |x|=8, and shooting
   through a long forbidden region amplifies the growing solution. brmp is EXPECTED. *)
(* no message declared here: brmp fires on the E0 bracket, not this one. The E1 root is
   nonetheless 5 orders LESS accurate (1.3e-6 vs 6e-12), which is the physics: kappa = 1 has
   not decayed by |x| = 8, so the shot amplifies the growing solution across a long forbidden
   region. Accuracy and loudness are independent; the quiet root is the inaccurate one. *)
report@VerificationTest[
  With[{e1 = u /. FindRoot[shootPT[u], {u, -0.8, -0.3}, Method -> "Brent"]},
   {Abs[e1 + 1/2] < 10.^-5, Abs[e1 + 1/2] > 10.^-8}], {True, True},
  TestID -> "3.6 the weakly bound PT state is only 1e-6-hard: shooting's real limitation"];

(* ================= 3.7 parity ================= *)
(* W and Vpt are the SHIPPED symbols: W is the generic (unassigned) potential and Vpt is the
   Poschl-Teller FAMILY bound by 3.1 above, exactly as in the document. Renaming either would
   test a transcription instead of the document. The generic claim is only meaningful because
   W is unassigned, so assert that first.
   NOT Head[W] === Symbol: that is TRUE of every symbol, assigned or not, so it cannot fail.
   Emptiness of the value lists is the property actually claimed. *)
report@VerificationTest[{OwnValues[W], DownValues[W]}, {{}, {}},
  TestID -> "3.7 W is genuinely unassigned"];
(* f is ALSO a generic test function here and in 3.9, and 3.3's completeness target must not
   have shadowed it. If it had, comm and 3.9's Integrate[delta f] go concrete and pass vacuously. *)
report@VerificationTest[{OwnValues[f], DownValues[f]}, {{}, {}},
  TestID -> "3.7/3.9 f is genuinely unassigned (3.3 used fC, not f)"];
report@VerificationTest[Vpt[x, 1], -Sech[x]^2,
  TestID -> "3.7 Vpt[x,1] is the PT well of 3.1, in scope and even"];
comm = Simplify[(-(1/2) D[f[-x], {x, 2}] + W[x] f[-x])
    - ((-(1/2) D[f[x], {x, 2}] + W[x] f[x]) /. x -> -x)];
report@VerificationTest[Simplify[comm - f[-x] (W[x] - W[-x])], 0,
  TestID -> "3.7 [H,P] = (W(x)-W(-x)) f(-x) for GENERIC W, i.e. the check CAN fail"];
report@VerificationTest[comm === 0, False,
  TestID -> "3.7 the generic commutator is NOT zero (the claim is not vacuous)"];
report@VerificationTest[Simplify[comm /. W[-x] -> W[x]], 0,
  TestID -> "3.7 [H,P] = 0 when W is even"];
report@VerificationTest[
  Simplify[(-(1/2) D[f[-x], {x, 2}] + Vpt[x, 1] f[-x])
    - ((-(1/2) D[f[x], {x, 2}] + Vpt[x, 1] f[x]) /. x -> -x)], 0,
  TestID -> "3.7 [H,P] = 0 for the concrete even V = -Sech^2"];
(* the CONSEQUENCE, tested on states rather than left as an inference: 3.4 SPLIT the well into
   even/odd sectors by assumption, so the split must be verified on the well's own eigenstates. *)
{fvals, ffuns} = NDEigensystem[{-(1/2) D[u[x], {x, 2}] + well[x] u[x],
    DirichletCondition[u[x] == 0., True]}, u[x], {x, -12., 12.}, 4,
   Method -> {"SpatialDiscretization" -> {"FiniteElement",
      {"MeshOptions" -> {MaxCellMeasure -> 0.005}}},
     "Eigensystem" -> {"Arnoldi", "Shift" -> -12.5}}];
pExp[g_] := NIntegrate[g (g /. x -> -x), {x, -12., 12.}]/NIntegrate[g^2, {x, -12., 12.}];
report@VerificationTest[
  Round[pExp[#], 0.001] & /@ ffuns[[Ordering[fvals]]], {1., -1., 1., -1.},
  TestID -> "3.7 <P> = +1,-1,+1,-1 on the finite well's OWN eigenstates (3.4's split verified)"];
(* the band must be tight enough to mean something: the raw values are within 2e-6 of +-1 *)
report@VerificationTest[
  Max@Abs[(Abs[pExp[#]] & /@ ffuns[[Ordering[fvals]]]) - 1.] < 10.^-5, True,
  TestID -> "3.7 |<P>| = 1 to 1e-5, not merely to the Round band"];
(* CONTROL: offset the well so V is NOT even, and <P> must stop being +-1. Without this the
   test could be measuring the solver's symmetry rather than the potential's. *)
{ovals, ofuns} = NDEigensystem[{-(1/2) D[u[x], {x, 2}]
     + Piecewise[{{-12.5, Abs[x - 0.3] < 1.}}, 0.] u[x],
    DirichletCondition[u[x] == 0., True]}, u[x], {x, -12., 12.}, 4,
   Method -> {"SpatialDiscretization" -> {"FiniteElement",
      {"MeshOptions" -> {MaxCellMeasure -> 0.005}}},
     "Eigensystem" -> {"Arnoldi", "Shift" -> -12.5}}];
(* ncvb: the offset state is NOT even, so g(x) g(-x) has a kink the quadrature struggles with.
   That is the control working, not failing; declare it. *)
report@VerificationTest[
  AnyTrue[Abs[pExp[#]] & /@ ofuns[[Ordering[ovals]]], # < 0.9 &], True,
  {NIntegrate::ncvb},
  TestID -> "3.7 CONTROL: shift the well off centre and <P> stops being +-1"];

(* ================= 3.8 node theorem ================= *)
(* SYMBOLIC in n and L: Solve returns the whole node family, not 5 pinned values *)
report@VerificationTest[
  FullSimplify[Sin[n Pi (kk L/n)/L], Assumptions -> Element[{n, kk}, Integers]], 0,
  TestID -> "3.8 psi_n vanishes at x = k L/n for ARBITRARY integer k, n, symbolic L"];
report@VerificationTest[
  Length@Solve[Sin[n Pi x/L] == 0 && 0 < x < L, x, Reals,
    Assumptions -> Element[n, Integers] && n >= 1 && L > 0], 2,
  TestID -> "3.8 Solve returns the node family symbolically (even + odd branches)"];
report@VerificationTest[Table[CountRoots[HermiteH[n, y], y], {n, 0, 6}], {0, 1, 2, 3, 4, 5, 6},
  TestID -> "3.8 oscillator state n has n nodes (different potential, same law)"];
(* CountRoots is not vacuously the degree: it counts REAL roots and so can fail *)
report@VerificationTest[CountRoots[y^6 + 1, y], 0,
  TestID -> "3.8 CountRoots counts REAL roots (degree 6, zero real): the check can fail"];
(* THE REAL TEST: a potential with NO closed-form eigenstates. Sin and HermiteH wear their
   nodes on their sleeve; the finite well does not. *)
nodes[f_] := With[{v = Table[f, {x, -11.9, 11.9, 0.005}]},
   With[{keep = Pick[v, UnitStep[Abs[v] - 1.*^-6 Max[Abs[v]]], 1]},
    Count[Sign[Most[keep] Rest[keep]], -1]]];
report@VerificationTest[nodes /@ ffuns[[Ordering[fvals]]], {0, 1, 2, 3},
  TestID -> "3.8 node theorem holds on the finite well, which has NO closed-form eigenstates"];
(* mutation controls: the counter must be able to be wrong *)
report@VerificationTest[{nodes[Sin[6 Pi x/24]], nodes[1. + 0 x], nodes[Sin[Pi x/24]]}, {5, 0, 1},
  TestID -> "3.8 node counter controls: 5 zeros -> 5, constant -> 0, 1 zero -> 1"];
(* the MASK is the load-bearing part and the controls above never touch it (they have large
   amplitude everywhere). Unmasked, the FEM tails give 153 nodes for the nodeless ground state.
   So show the answer is threshold-INDEPENDENT across orders of magnitude: a plateau, not a
   fudge factor tuned until it agreed. *)
nodesTol[f_, tol_] := With[{v = Table[f, {x, -11.9, 11.9, 0.005}]},
   With[{keep = Pick[v, UnitStep[Abs[v] - tol Max[Abs[v]]], 1]},
    Count[Sign[Most[keep] Rest[keep]], -1]]];
report@VerificationTest[
  nodesTol[First@ffuns[[Ordering[fvals]]], 0] > 100, True,
  TestID -> "3.8 unmasked, roundoff in the tails fakes >100 nodes on the NODELESS ground state"];
report@VerificationTest[
  Table[nodesTol[#, tol] & /@ ffuns[[Ordering[fvals]]], {tol, {10.^-14, 10.^-10, 10.^-6, 10.^-2, 0.5}}],
  ConstantArray[{0, 1, 2, 3}, 5],
  TestID -> "3.8 the mask is a PLATEAU: same count for tol spanning 1e-14 to 0.5"];

(* ================= 3.9 delta well ================= *)
report@VerificationTest[Integrate[DiracDelta[x] f[x], {x, -e, e}, Assumptions -> e > 0], f[0],
  TestID -> "3.9 Integrate is exact on a point measure"];
(* G7 declared, not suppressed: NIntegrate MISSES it and says so *)
report@VerificationTest[NIntegrate[DiracDelta[x] Exp[-x^2], {x, -1, 1}], 0.,
  {NIntegrate::izero}, TestID -> "3.9 G7: NIntegrate returns 0. on a delta (message declared)"];
\[Psi]R[x_] := Sqrt[k] Exp[-k x];
\[Psi]L[x_] := Sqrt[k] Exp[+k x];
report@VerificationTest[
  Solve[{\[Psi]R'[0] - \[Psi]L'[0] == -2 g \[Psi]R[0], k > 0}, k, Reals],
  {{k -> ConditionalExpression[g, g > 0]}}, TestID -> "3.9 jump condition gives kappa = g"];
report@VerificationTest[
  -k^2/2 /. First@Solve[{\[Psi]R'[0] - \[Psi]L'[0] == -2 g \[Psi]R[0], k > 0}, k, Reals],
  ConditionalExpression[-g^2/2, g > 0], TestID -> "3.9 E = -g^2/2"];
report@VerificationTest[
  Integrate[k Exp[-2 k Abs[x]], {x, -Infinity, Infinity}, Assumptions -> k > 0], 1,
  TestID -> "3.9 cusp state normalized (= the Part-1 state at kappa=g)"];
(* the delta well IS the zero-width limit of 3.4's finite well at fixed g = 2 a V0.
   Prose asserted this; compute it. Two routes sharing no machinery must meet at -g^2/2. *)
report@VerificationTest[
  With[{es = Table[With[{u0 = a Sqrt[2 (1/(2 a))]},
        -(u0^2 - First[uAll[u0]]^2)/(2 a^2)], {a, {0.5, 0.1, 0.02, 0.004}}]},
   {And @@ Thread[Abs[Rest[es] + 0.5] < Abs[Most[es] + 0.5]], Abs[Last[es] + 0.5] < 0.005}],
  {True, True},
  TestID -> "3.9 finite well -> delta well: E -> -g^2/2 as a -> 0 at fixed g = 2 a V0"];
(* the k>0 constraint is load-bearing: without it Solve offers the spurious kappa=0 first *)
report@VerificationTest[
  MemberQ[k /. Solve[\[Psi]R'[0] - \[Psi]L'[0] == -2 g \[Psi]R[0], k], 0], True,
  TestID -> "3.9 TRAP: unconstrained Solve admits the spurious root kappa=0"];
(* and it is returned FIRST, which is what makes First@Solve silently wrong *)
report@VerificationTest[
  First[k /. Solve[\[Psi]R'[0] - \[Psi]L'[0] == -2 g \[Psi]R[0], k]], 0,
  TestID -> "3.9 TRAP: the spurious root is returned FIRST"];
(* GATE (hh): the symbolic route really is closed, which is why the jump condition is the route.
   InverseFunction::ifun fires while it declines; declare it rather than suppress it. *)
report@VerificationTest[
  Head@DSolveValue[-(1/2) u''[x] - g DiracDelta[x] u[x] == En u[x], u[x], x], DSolveValue,
  {InverseFunction::ifun},
  TestID -> "3.9 GATE: DSolve returns unevaluated on the delta potential"];

(* ================= 3.10 quantum bouncer ================= *)
report@VerificationTest[
  FullSimplify[-(1/2) D[AiryAi[(2 F)^(1/3) (x - En/F)], {x, 2}]
     + F x AiryAi[(2 F)^(1/3) (x - En/F)] - En AiryAi[(2 F)^(1/3) (x - En/F)],
   Assumptions -> F > 0], 0,
  TestID -> "3.10 Airy solves linear-V TISE for ARBITRARY E and F (the ODE quantizes nothing)"];
Ebounce[n_, F_] := -(F^2/2)^(1/3) AiryAiZero[n];
report@VerificationTest[
  Max@Abs[(Sort@NDEigenvalues[{-(1/2) D[u[x], {x, 2}] + x u[x],
        DirichletCondition[u[x] == 0, True]}, u[x], {x, 0, 25.}, 4,
       Method -> {"SpatialDiscretization" -> {"FiniteElement",
          {"MeshOptions" -> {MaxCellMeasure -> 0.002}}}}]
      - Table[N@Ebounce[n, 1], {n, 1, 4}])/Table[N@Ebounce[n, 1], {n, 1, 4}]] < 10.^-5,
  True, TestID -> "3.10 NDEigensystem reproduces the Airy-zero spectrum"];
(* the F^(2/3) scaling is the experimentally measured signature; assert it symbolically *)
report@VerificationTest[
  FullSimplify[Ebounce[n, F]/Ebounce[n, 1] - F^(2/3), Assumptions -> F > 0], 0,
  TestID -> "3.10 levels scale as F^(2/3), symbolically in F"];
(* INVARIANT: unit norm on the bound state, not just its energy. Closes in closed form. *)
report@VerificationTest[
  Abs[NIntegrate[AiryAi[2^(1/3) x + AiryAiZero[1]]^2, {x, 0, Infinity}]
    - N[AiryAiPrime[AiryAiZero[1]]^2/2^(1/3)]] < 10.^-8, True,
  TestID -> "3.10 bouncer norm: Int Ai^2 = Ai'(a_n)^2/(2F)^(1/3), closed form"];
(* the wall condition that PRODUCED the spectrum is satisfied identically, every n *)
report@VerificationTest[Table[AiryAi[AiryAiZero[n]], {n, 1, 3}], {0, 0, 0},
  TestID -> "3.10 wall condition psi_n(0) = Ai(a_n) = 0 exactly, for every n and every F"];
(* CLASSICAL LIMIT with content: n -> Infinity against WKB. hbar -> 0 would be vacuous. *)
report@VerificationTest[
  With[{errs = Table[N@Abs[(3 Pi (n - 1/4)/2)^(2/3) (1/2)^(1/3)/Ebounce[n, 1] - 1],
      {n, {1, 5, 20, 100}}]},
   And @@ Thread[Rest[errs] < Most[errs]] && Last[errs] < 10.^-6], True,
  TestID -> "3.10 correspondence: WKB -> exact as n -> Infinity (err 7.6e-3 down to 4.7e-7)"];
(* and the WKB/exact ratio is F-INDEPENDENT, so this is a statement in n alone *)
report@VerificationTest[
  FreeQ[FullSimplify[((3 Pi (n - 1/4)/2)^(2/3) (F^2/2)^(1/3))/(-(F^2/2)^(1/3) AiryAiZero[n]),
    Assumptions -> F > 0], F], True,
  TestID -> "3.10 WKB/exact ratio is independent of F"];
(* the Maslov -1/4 is load-bearing: drop it and the ground state is an order worse *)
report@VerificationTest[
  N@Abs[(3 Pi (1 - 1/4)/2)^(2/3) (1/2)^(1/3)/Ebounce[1, 1] - 1]
    < N@Abs[(3 Pi 1/2)^(2/3) (1/2)^(1/3)/Ebounce[1, 1] - 1]/5, True,
  TestID -> "3.10 the Maslov index -1/4 matters: dropping it is >5x worse at n=1"];
(* the prose claims the DEFAULT mesh is NOT converged (0.2% at n=1, 1.8% at n=4).
   That is a claim about being wrong, so it needs a test as much as the right answer does. *)
report@VerificationTest[
  With[{rel = (Sort@NDEigenvalues[{-(1/2) D[u[x], {x, 2}] + x u[x],
         DirichletCondition[u[x] == 0, True]}, u[x], {x, 0, 25.}, 4]
       - Table[N@Ebounce[n, 1], {n, 1, 4}])/Table[N@Ebounce[n, 1], {n, 1, 4}]},
   {10.^-3 < Abs[First[rel]] < 10.^-2, 10.^-2 < Abs[Last[rel]] < 10.^-1}], {True, True},
  TestID -> "3.10 default mesh is NOT converged: ~0.2% at n=1 degrading to ~1.8% at n=4"];

(* the 3.10 closing note: DEigensystem DOES return on the bouncer, as Root objects agreeing with
   the Airy zeros, but emits internal ReplaceAll::reps. Declared, so the test asserts it fires.
   Three traps in one test, each of which made it FAIL while the physics was right:
   (1) Both sides are numericized to 20 digits FIRST: comparing an exact Root against an exact
       AiryAiZero symbolically exhausts $MaxExtraPrecision (Max::meprec) and the inequality
       never evaluates to True.
   (2) reps fires repeatedly, which itself triggers General::stop, a SECOND message.
   (3) VerificationTest matches declared messages by MULTIPLICITY, not set membership: reps
       fires 3 times (once per eigenvalue), so declaring it once still fails. Hence 3 copies.
       Brittle by construction, which is the honest cost of not reaching for Quiet. *)
report@VerificationTest[
  Max@Abs[N[DEigenvalues[{-(1/2) D[u[x], {x, 2}] + x u[x], DirichletCondition[u[x] == 0, True]},
       u[x], {x, 0, 25}, 3], 20] - N[Table[Ebounce[n, 1], {n, 1, 3}], 20]] < 10.^-12, True,
  {ReplaceAll::reps, ReplaceAll::reps, ReplaceAll::reps, General::stop},
  TestID -> "3.10 DEigensystem returns Root objects agreeing with the Airy zeros (reps declared)"];

WriteString[$Output, "PART03-HARNESS-DONE\n"];
