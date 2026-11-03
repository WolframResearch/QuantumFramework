(* Exact (non-QF) checks for section X. No Fock truncation: finite group matrices, Jordan-Wigner fermions, Bargmann functions, exact integrals. *)
emit[k_, v_] := WriteString["stdout", k, " -> ", ToString[InputForm[v]], "\n"];
dag = ConjugateTranspose;
assum = r > 0 && Element[ph, Reals];
xi = r Exp[I ph];

(* X-004 compact: spin-1/2 and spin-1, symbolic r, ph *)
su2[jp_, jm_, j0_, z_] := MatrixExp[xi jp - Conjugate[xi] jm] - MatrixExp[z jp] . MatrixExp[Log[1 + Abs[z]^2] j0] . MatrixExp[-Conjugate[z] jm];
sp = {{0, 1}, {0, 0}}; sm = Transpose[sp]; sz = {{1, 0}, {0, -1}}/2;
s1p = Sqrt[2] {{0, 1, 0}, {0, 0, 1}, {0, 0, 0}}; s1m = Transpose[s1p]; s1z = DiagonalMatrix[{1, 0, -1}];
alg[p_, m_, z_] := {z . p - p . z - p, z . m - m . z + m, p . m - m . p - 2 z};
emit["su2 algebra spin1/2, spin1", {alg[sp, sm, sz], alg[s1p, s1m, s1z]}];
emit["su2 spin1/2 0<r<pi/2", FullSimplify[su2[sp, sm, sz, Exp[I ph] Tan[r]], 0 < r < Pi/2 && Element[ph, Reals]]];
emit["su2 spin1 0<r<pi/2", FullSimplify[su2[s1p, s1m, s1z, Exp[I ph] Tan[r]], 0 < r < Pi/2 && Element[ph, Reals]]];
emit["su2 spin1/2 pi/2<r<pi: LHS and RHS", FullSimplify[{MatrixExp[xi sp - Conjugate[xi] sm], MatrixExp[xi sp - Conjugate[xi] sm] - su2[sp, sm, sz, Exp[I ph] Tan[r]]}, Pi/2 < r < Pi && Element[ph, Reals]]];
emit["su2 spin1 pi/2<r<pi", FullSimplify[su2[s1p, s1m, s1z, Exp[I ph] Tan[r]], Pi/2 < r < Pi && Element[ph, Reals]]];

(* X-004 noncompact: 2x2 nonunitary rep K0 = sz, K+ = sp, K- = -sm *)
emit["su11 algebra", {sz . sp - sp . sz - sp, sz . (-sm) - (-sm) . sz - sm, (-sm) . sp - sp . (-sm) - 2 sz}];
emit["su11 identity", FullSimplify[With[{z = Exp[I ph] Tanh[r]}, MatrixExp[xi sp + Conjugate[xi] sm] - MatrixExp[z sp] . MatrixExp[Log[1 - Abs[z]^2] sz] . MatrixExp[Conjugate[z] sm]], assum]];

(* Jordan-Wigner *)
jw[n_] := Table[KroneckerProduct @@ Join[ConstantArray[{{1, 0}, {0, -1}}, j - 1], {{{0, 1}, {0, 0}}}, ConstantArray[IdentityMatrix[2], n - j]], {j, n}];
{c1, c2} = jw[2]; i4 = IdentityMatrix[4]; vac4 = UnitVector[4, 1];
Jp = dag[c1] . dag[c2]; Jm = c2 . c1; J0 = (dag[c1] . c1 + dag[c2] . c2 - i4)/2;
emit["fermion pair algebra; Jm-Jp^dag; algebra with c1.c2", {alg[Jp, Jm, J0], Jm - dag[Jp], alg[Jp, c1 . c2, J0][[3]]}];
emit["fermion pair state (all r)", FullSimplify[MatrixExp[xi Jp - Conjugate[xi] Jm] . vac4 - (Cos[r] vac4 + Exp[I ph] Sin[r] dag[c1] . dag[c2] . vac4), assum]];

(* X-003 two-mode squeezed vacuum in Bargmann form: a_i -> d/dz_i, ad_i -> z_i. psi(s) = S2(s xi)|00>, d psi/ds = G psi *)
psi = Sech[s r] Exp[-Exp[I ph] Tanh[s r] z1 z2];
emit["two-mode vacuum ODE residual, initial value", {FullSimplify[D[psi, s] - (r Exp[-I ph] D[psi, z1, z2] - r Exp[I ph] z1 z2 psi), assum && s > 0], psi /. s -> 0}];
emit["(z1 z2)^n/n! coefficient is |n,n> amplitude/sqrt(n!^2) check", Simplify[SeriesCoefficient[Exp[-tau z1 z2], {z1, 0, 3}, {z2, 0, 3}] - (-tau)^3/3!]];

(* X-006 boson normalization: Bargmann norm^2 = 1/Sqrt[Det A], A real 2N x 2N; compare Det A with Det(I - Z Zbar) *)
gaussDet[n_] := Module[{X = Array[x, {n, n}], Y = Array[y, {n, n}], A},
  X = (X + Transpose[X])/2; Y = (Y + Transpose[Y])/2;
  A = ArrayFlatten[{{IdentityMatrix[n] + X, -Y}, {-Y, IdentityMatrix[n] - X}}];
  Expand[Det[A] - ComplexExpand[Det[IdentityMatrix[n] - (X + I Y) . (X - I Y)]]]];
emit["boson Gaussian det identity N=1,2,3", gaussDet /@ {1, 2, 3}];
(* exponent check: -|z|^2 - Re(z^T Z z) = -v^T A v *)
emit["boson exponent form N=2", Module[{n = 2, X, Y, A, zz, vv},
   X = Array[x, {n, n}]; X = (X + Transpose[X])/2; Y = Array[y, {n, n}]; Y = (Y + Transpose[Y])/2;
   A = ArrayFlatten[{{IdentityMatrix[n] + X, -Y}, {-Y, IdentityMatrix[n] - X}}];
   zz = Array[p, n] + I Array[q, n]; vv = Join[Array[p, n], Array[q, n]];
   Expand[ComplexExpand[-Conjugate[zz] . zz - Re[zz . (X + I Y) . zz]] + vv . A . vv]]];

(* X-006 fermion Gaussian vacuum, 4 modes, exact rational complex antisymmetric Z *)
fs = jw[4]; vac16 = UnitVector[16, 1];
Zf = {{0, 1/2 + I/3, -1/5, 2/7 - I/4}, {0, 0, 3/4 I, -1/3}, {0, 0, 0, 1/6 + I/2}, {0, 0, 0, 0}}; Zf = Zf - Transpose[Zf];
psiU = MatrixExp[-1/2 Sum[Zf[[i, j]] dag[fs[[i]]] . dag[fs[[j]]], {i, 4}, {j, 4}]] . vac16;
emit["fermion vacuum: annihilation (exact), (norm^2)^2 - det(I+ZZ^dag)", {Union[Flatten[Simplify[Table[(fs[[i]] + Sum[Zf[[i, j]] dag[fs[[j]]], {j, 4}]) . psiU, {i, 4}]]]], Simplify[(Conjugate[psiU] . psiU)^2 - Det[IdentityMatrix[4] + Zf . dag[Zf]]]}];

(* X-007 number-conserving, bosons, Bargmann: normal-ordered :exp(sum z_i B_ij d_j): f = f(z + B^T z) (Taylor, exact on a cubic), flow of e^{t sum z_i h_ij d_j} is f(e^{t h^T} z) *)
hB = Array[b, {2, 2}]; zv = {z1, z2};
fpoly = Sum[cf[i, j] z1^i z2^j, {i, 0, 3}, {j, 0, 3 - i}];
noOp[f_] := Sum[1/k! Sum[(Times @@ (zv[[#[[1]]]] & /@ t)) (Times @@ (hB[[#[[1]], #[[2]]]] & /@ t)) D[f, Sequence @@ ({zv[[#[[2]]]], 1} & /@ t)], {t, Tuples[Tuples[{1, 2}, 2], k]}], {k, 0, 3}];
emit["normal-ordered Taylor identity on generic cubic", Expand[noOp[fpoly] - (fpoly /. Thread[zv -> zv + Transpose[hB] . zv])]];
hh = Array[hm, {2, 2}];
flowRes = Module[{E1 = MatrixExp[t Transpose[hh]], F},
   F = ff @@ (E1 . zv);
   D[F, t] - Sum[zv[[i]] hh[[i, j]] D[F, zv[[j]]], {i, 2}, {j, 2}]];
emit["flow ODE residual at a random numeric point", Chop[(flowRes /. ff -> (Exp[#1^2 #2] + #1 #2^3 &)) /. MapAt[N, {hm[1, 1] -> 1/3, hm[1, 2] -> -2/5 + I/7, hm[2, 1] -> 3/4, hm[2, 2] -> -1/2 + I/3, t -> 7/10, z1 -> 2/3 - I/5, z2 -> -1/4 + I}, {All, 2}], 10^-12]];

(* X-007 fermions (3 modes, full 8x8 JW): e^{c^dag h c} vs CAR normal-ordered exponential; projectors; trace *)
gs = jw[3]; i8 = IdentityMatrix[8]; proj0 = Outer[Times, UnitVector[8, 1], UnitVector[8, 1]];
ferNO[B_] := i8 + Sum[1/kk! Sum[With[{is = t[[;; kk]], js = t[[kk + 1 ;;]]},
      (Times @@ MapThread[B[[#1, #2]] &, {is, js}]) (Dot @@ (dag[gs[[#]]] & /@ is)) . (Dot @@ (gs[[#]] & /@ Reverse[js]))],
     {t, Tuples[Range[3], 2 kk]}], {kk, 1, 3}];
emit["fermion projector :e^{-sum N}: and prod(1-N_i) vs |0><0| (exact)", {ferNO[-IdentityMatrix[3]] - proj0 // Flatten // Union, (Dot @@ ((i8 - dag[#] . #) & /@ gs)) - proj0 // Flatten // Union}];
(* exact h with rational e^h: h = P.diag(Log 2, Log 3, Log(1/5)).P^-1 *)
Pm = {{1, I, 0}, {1/2, 1, 1}, {0, -1, 2}}; hf = Pm . DiagonalMatrix[{Log[2], Log[3], Log[1/5]}] . Inverse[Pm]; Mf = Pm . DiagonalMatrix[{2, 3, 1/5}] . Inverse[Pm];
eFer = MatrixExp[Sum[hf[[i, j]] dag[gs[[i]]] . gs[[j]], {i, 3}, {j, 3}]];
emit["fermion e^{c^dag h c} - :exp[c^dag(e^h-1)c]: (exact); Tr - det(1+e^h)", {Union[Flatten[Simplify[eFer - ferNO[Mf - IdentityMatrix[3]]]]], Simplify[Tr[eFer] - Det[IdentityMatrix[3] + Mf]]}];

(* X-007 boson trace, two modes: Tr Gamma(M) = Integral <z|Gamma(M)|z> e^{-|z|^2} d^4/pi^2 with <z|Gamma(M)|w> = e^{zbar.M.w} (unnormalized). Exact rational non-Hermitian M. *)
Mb = {{1/3, 1/5 + I/7}, {-1/4, I/6}};
zc = {p1 + I q1, p2 + I q2};
expo = ComplexExpand[-Conjugate[zc] . (IdentityMatrix[2] - Mb) . zc];
emit["boson two-mode trace integral vs 1/det(I-M)", Simplify[Integrate[Exp[expo], {p1, -Infinity, Infinity}, {q1, -Infinity, Infinity}, {p2, -Infinity, Infinity}, {q2, -Infinity, Infinity}]/Pi^2 - 1/Det[IdentityMatrix[2] - Mb]]];
emit["spectral radius and norm of Mb", N[{Max[Abs[Eigenvalues[Mb]]], Norm[Mb]}]];

(* X-009 fermions: symbolic anticommutators for b = u c + v c^dag, 2 modes *)
uS = Array[uu, {2, 2}]; vS = Array[vv, {2, 2}];
cl = {c1, c2};
bF = Table[Sum[uS[[i, j]] cl[[j]] + vS[[i, j]] dag[cl[[j]]], {j, 2}], {i, 2}];
bFd = Table[Sum[Conjugate[uS[[i, j]]] dag[{c1, c2}[[j]]] + Conjugate[vS[[i, j]]] {c1, c2}[[j]], {j, 2}], {i, 2}];
emit["fermion {b_i,b_j^dag} - (uu^dag+vv^dag)_ij, {b_i,b_j} - (uv^T+vu^T)_ij", {Union[Flatten[Simplify[Table[bF[[i]] . bFd[[j]] + bFd[[j]] . bF[[i]] - (uS . dag[uS] + vS . dag[vS])[[i, j]] i4, {i, 2}, {j, 2}]]]], Union[Flatten[Simplify[Table[bF[[i]] . bF[[j]] + bF[[j]] . bF[[i]] - (uS . Transpose[vS] + vS . Transpose[uS])[[i, j]] i4, {i, 2}, {j, 2}]]]]}];

(* X-009 fermions: full JW matrices, random complex even quadratic H, 3 and 4 modes (complete operators) *)
mx[m_] := Max[Abs[Flatten[N[m]]]];
ferX009[n_, seed_] := Module[{cs = jw[n], Af, Bf, H, U, uu1, vv1, Zz, Ww, hh1, c0, left, mid},
  SeedRandom[seed];
  Af = RandomComplex[{-1 - I, 1 + I}, {n, n}]; Af = (Af + dag[Af])/2;
  Bf = RandomComplex[{-1 - I, 1 + I}, {n, n}]; Bf = (Bf - Transpose[Bf])/2;
  H = Sum[Af[[i, j]] dag[cs[[i]]] . cs[[j]] + 1/2 (Bf[[i, j]] dag[cs[[i]]] . dag[cs[[j]]] + Conjugate[Bf[[i, j]]] cs[[j]] . cs[[i]]), {i, n}, {j, n}];
  U = MatrixExp[-I H];
  uu1 = Table[Tr[dag[cs[[k]]] . (U . cs[[i]] . dag[U])]/2^(n - 1), {i, n}, {k, n}];
  vv1 = Table[Tr[cs[[k]] . (U . cs[[i]] . dag[U])]/2^(n - 1), {i, n}, {k, n}];
  Zz = Inverse[uu1] . vv1; Ww = Conjugate[vv1] . Inverse[uu1]; hh1 = MatrixLog[Inverse[uu1]]; c0 = U[[1, 1]];
  left = MatrixExp[-1/2 Sum[Zz[[i, j]] dag[cs[[i]]] . dag[cs[[j]]], {i, n}, {j, n}]]; mid = MatrixExp[Sum[hh1[[i, j]] dag[cs[[i]]] . cs[[j]], {i, n}, {j, n}]];
  <|"canonical action" -> Max[Table[mx[U . cs[[i]] . dag[U] - Sum[uu1[[i, j]] cs[[j]] + vv1[[i, j]] dag[cs[[j]]], {j, n}]], {i, n}]],
    "uu^dag+vv^dag-I" -> mx[uu1 . dag[uu1] + vv1 . dag[vv1] - IdentityMatrix[n]], "uv^T+vu^T" -> mx[uu1 . Transpose[vv1] + vv1 . Transpose[uu1]],
    "Z,W antisym" -> {mx[Zz + Transpose[Zz]], mx[Ww + Transpose[Ww]]}, "|c0|-det(uu^dag)^(1/4)" -> Abs[Abs[c0] - Det[uu1 . dag[uu1]]^(1/4)],
    "factorization (eta/2 = -1/2)" -> mx[c0 left . mid . MatrixExp[-1/2 Sum[Ww[[i, j]] cs[[i]] . cs[[j]], {i, n}, {j, n}]] - U],
    "with +1/2 instead" -> mx[c0 left . mid . MatrixExp[1/2 Sum[Ww[[i, j]] cs[[i]] . cs[[j]], {i, n}, {j, n}]] - U]|>];
emit["fermion X-009", {ferX009[3, 21], ferX009[4, 22]}];

