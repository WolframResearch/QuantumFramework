#!/usr/bin/env wolframscript
(* Verification demos for the WL-symbolic-novelty re-examination.            *)
(* Each block prints a labeled result; correctness is asserted with ===.     *)
(* No Quiet, no Print. Core WL only (no QF paclet needed for these claims).  *)

omega = Exp[I Pi/4];

(* ======================================================================== *)
(* CANDIDATE (a): the variety / Groebner point-count over F_2.              *)
(*   amplitude of a Toffoli-Hadamard circuit  = gap(f) = #{f=0} - #{f=1}.    *)
(* ======================================================================== *)

(* --- (a1) Nullstellensatz: GB = {1}  <=>  no F_2 solution  <=>  amp 0. --- *)
inconsistent = GroebnerBasis[{x + y, x, y + 1}, {x, y}, Modulus -> 2];
Echo[inconsistent, "a1 GB of inconsistent F2 system (expect {1})"];
Echo[inconsistent === {1}, "a1 amplitude-zero detected? (expect True)"];

consistent = GroebnerBasis[{x + y + z}, {x, y, z}, Modulus -> 2];
Echo[consistent, "a1 GB of one linear F2 constraint (= Gaussian elim)"];

(* --- (a2) LINEAR (Clifford-sign / [HH]) fragment: point count via rank,    *)
(*          genuinely polynomial, NO 2^n enumeration. {x1+x2+x3=0} in F_2^3. *)
M = {{1, 1, 1}};                            (* coefficient matrix over F_2 *)
nLin = 3; rk = MatrixRank[M, Modulus -> 2];
countLinRank = 2^(nLin - rk);               (* affine solution count *)
countLinBrute = Count[Tuples[{0, 1}, 3], a_ /; Mod[Total[a], 2] == 0];
Echo[{rk, countLinRank, countLinBrute}, "a2 {rank, count via rank, brute} (expect {1,4,4})"];
Echo[countLinRank === countLinBrute, "a2 rank count == brute? (expect True)"];

(* --- (a3) NONLINEAR (magic) point count via GB quotient dimension.        *)
(*   #{x in F_2^n : f(x)=0} = # standard monomials of (f, x_i^2-x_i).       *)
(*   A squarefree monomial m is standard iff its normal form equals m.      *)
fieldEqs[vars_] := (#^2 - # &) /@ vars;
pointCountGB[polys_, vars_] := Module[{gb, mons},
   gb = GroebnerBasis[Join[polys, fieldEqs[vars]], vars, Modulus -> 2];
   If[gb === {1} || MemberQ[gb, 1], 0,
    mons = Times @@@ Subsets[vars];          (* the 2^n multilinear monomials *)
    Count[mons, m_ /; PolynomialMod[
        Last[PolynomialReduce[m, gb, vars, Modulus -> 2]] - m, 2] === 0]]];
pointCountBrute[polys_, vars_] := Count[Tuples[{0, 1}, Length[vars]],
    a_ /; And @@ (PolynomialMod[# /. Thread[vars -> a], 2] === 0 & /@ polys)];

vars3 = {x1, x2, x3};
fNonlin = x1 x2 + x3;                        (* a genuine degree-2 magic term *)
cGB = pointCountGB[{fNonlin}, vars3];
cBr = pointCountBrute[{fNonlin}, vars3];
gapVal = 2 cBr - 2^3;                         (* gap(f) = #{f=0} - #{f=1} *)
Echo[{cGB, cBr, gapVal}, "a3 {GB-quotient count, brute count, gap(f)}"];
Echo[cGB === cBr, "a3 GB quotient-dim == brute point count? (expect True)"];

(* cross-check on a cubic (CCZ-type) term too *)
vars4 = {x1, x2, x3, x4};
fCubic = x1 x2 x3 + x4;
Echo[pointCountGB[{fCubic}, vars4] === pointCountBrute[{fCubic}, vars4],
    "a3 cubic-term quotient==brute? (expect True)"];

(* ======================================================================== *)
(* CANDIDATE (b): SYMBOLIC-PARAMETER amplitude (free angle theta).          *)
(*   One symbolic pass -> closed-form amplitude & expectation FUNCTION.     *)
(* ======================================================================== *)
(* Path-sum form of H P(theta) H |0>, P(theta)=diag(1,e^{i theta}):         *)
(*   <x|HPH|0> = (1/2) sum_{a in F_2} e^{i theta a} (-1)^{a x}.             *)
ampTheta[x_] := (1/2) Sum[Exp[I th a] (-1)^(a x), {a, 0, 1}];
amp0 = FullSimplify[ampTheta[0]];
amp1 = FullSimplify[ampTheta[1]];
Echo[{amp0, amp1}, "b amplitudes <0|..|0>,<1|..|0> as closed forms in theta"];

(* expectation <Z>(theta) as a closed-form trig function, one symbolic pass *)
H = {{1, 1}, {1, -1}}/Sqrt[2];
Pth = {{1, 0}, {0, Exp[I th]}};
Zm = {{1, 0}, {0, -1}};
psi = H . Pth . H . {1, 0};
zExp = FullSimplify[Conjugate[psi] . (Zm . psi), Assumptions -> th \[Element] Reals];
Echo[zExp, "b <Z>(theta) closed form (expect Cos[theta])"];
Echo[FullSimplify[zExp - Cos[th]] === 0, "b <Z>(theta) == Cos[theta]? (expect True)"];

(* the amplitudes equal the StabilizerFrame coefficients (1 +/- e^{i th})/2 *)
Echo[FullSimplify[{amp0, amp1} - {(1 + Exp[I th])/2, (1 - Exp[I th])/2}],
    "b match StabilizerFrame coeffs? (expect {0,0})"];

(* ======================================================================== *)
(* CANDIDATE (c): GF(2^k) trace character + the characteristic-2 wall.      *)
(* ======================================================================== *)
FF = FiniteField[2, 4];
elts = FF[#] & /@ Range[0, 15];
traceSum = Total[(-1)^(FiniteFieldElementTrace /@ elts)];
Echo[traceSum, "c additive character sum sum_{x in GF(2^4)} (-1)^Tr(x) (expect 0)"];
Echo[traceSum === 0, "c trace balanced? (expect True)"];

(* char-2 obstruction: a primitive 8th-root phase cannot sit in a sign.     *)
g = omega PauliMatrix[1];                     (* e^{i pi/4} X *)
gsq = FullSimplify[g . g];
Echo[gsq, "c (e^{i pi/4} X)^2 (expect i*Identity, NOT Identity)"];
Echo[gsq === I IdentityMatrix[2], "c equals i*I (not I)? (expect True)"];
Echo[FullSimplify[Eigenvalues[g]], "c eigenvalues of e^{i pi/4} X (no +1 eigenvalue)"];

(* the cyclotomic structure of omega collapses completely mod 2 *)
Echo[Cyclotomic[8, x], "c minimal polynomial of omega (expect 1+x^4)"];
Echo[Factor[1 + x^4, Modulus -> 2], "c 1+x^4 over F_2 (expect (1+x)^4: structure gone)"];

(* ======================================================================== *)
(* CANDIDATE (d) honest limit: RootSum gives LINEAR root-of-unity sums but   *)
(* NOT the quadratic Gauss sum sum_k e^{i pi k^2 / n}.                       *)
(* ======================================================================== *)
linSum = RootSum[1 + #^4 &, # &];             (* sum of 8th-roots' partner: a Vieta sum *)
Echo[linSum, "d RootSum[1+#^4&, #&] (a linear/Vieta power sum, closed form)"];
quadGauss = Table[{n, Simplify[Sum[Exp[I Pi k^2/n], {k, 0, n - 1}]]}, {n, 2, 6}];
Echo[quadGauss, "d quadratic Gauss sums sum_k e^{i pi k^2/n} (closed forms, NOT RootSum-shaped)"];

Echo["ALL DEMOS COMPLETE", "status"];
