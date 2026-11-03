(* Exact (non-truncated) checks for the C section that do not go through a QF call.
   Finite-dimensional group matrices, spin matrices, Jordan-Wigner fermion matrices,
   exterior algebra, Bargmann differential realization, exact series and sums. *)
res = <||>;
SetAttributes[put, HoldRest];
put[k_, v_] := (res[k] = TimeConstrained[v, 300, $TimedOut]);
mx[x_] := Max[Abs[Flatten[{N[x]}]]];
(* C-020: Siegel action, Im part and cocycle, two degrees of freedom *)
j4 = ArrayFlatten[{{0, IdentityMatrix[2]}, {-IdentityMatrix[2], 0}}];
sym4a = {{0.5, 0.2, -0.1, 0.3}, {0.2, -0.4, 0.6, 0.1}, {-0.1, 0.6, 0.3, -0.2}, {0.3, 0.1, -0.2, 0.7}};
sym4b = {{-0.3, 0.1, 0.4, 0.}, {0.1, 0.2, -0.5, 0.3}, {0.4, -0.5, 0.6, 0.1}, {0., 0.3, 0.1, -0.2}};
f1 = MatrixExp[j4 . sym4a]; f2 = MatrixExp[0.7 j4 . sym4b];
blocks[f_] := {f[[1 ;; 2, 1 ;; 2]], f[[1 ;; 2, 3 ;; 4]], f[[3 ;; 4, 1 ;; 2]], f[[3 ;; 4, 3 ;; 4]]};
act[f_, g_] := With[{b = blocks[f]}, (b[[3]] + b[[4]] . g) . Inverse[b[[1]] + b[[2]] . g]];
dl[f_, g_] := With[{b = blocks[f]}, Det[b[[1]] + b[[2]] . g]];
g0 = I IdentityMatrix[2] + {{0.3 + 0.1 I, 0.1}, {0.1, -0.2 + 0.2 I}};
With[{b = blocks[f1]},
  put["C020-ImGamma", mx[Im[act[f1, g0]] - Inverse[ConjugateTranspose[b[[1]] + b[[2]] . g0]] . Im[g0] . Inverse[b[[1]] + b[[2]] . g0]]]];
put["C020-action-composition", mx[act[f1 . f2, g0] - act[f1, act[f2, g0]]]];
put["C020-cocycle", mx[dl[f1 . f2, g0] - dl[f1, act[f2, g0]] dl[f2, g0]]];
put["C020-symmetric-posIm", {mx[act[f1, g0] - Transpose[act[f1, g0]]], Eigenvalues[Im[act[f1, g0]]]}];
(* C-021 Riccati, finite difference, S = sym4a = {{K, L},{L^T, G}} *)
With[{kk = sym4a[[1 ;; 2, 1 ;; 2]], ll = sym4a[[1 ;; 2, 3 ;; 4]], gg = sym4a[[3 ;; 4, 3 ;; 4]], h = 10.^-5, t0 = 0.4},
  With[{gm0 = act[MatrixExp[t0 j4 . sym4a], g0]},
   put["C021-Riccati", mx[(act[MatrixExp[(t0 + h) j4 . sym4a], g0] - act[MatrixExp[(t0 - h) j4 . sym4a], g0])/(2 h) -
       (-kk - ll . gm0 - gm0 . Transpose[ll] - gm0 . gg . gm0)]]]];

(* C-037 clock and shift *)
With[{nn = 5, om = Exp[2 Pi I/5]},
  With[{zc = DiagonalMatrix[om^Range[0, nn - 1]], xs = RotateRight[IdentityMatrix[nn]]},
   put["C037-ZX=wXZ", mx[zc . xs - om xs . zc]];
   put["C037-powers", {mx[MatrixPower[zc, nn] - IdentityMatrix[nn]], mx[MatrixPower[xs, nn] - IdentityMatrix[nn]]}]]];

(* C-038 .. C-042 SU(2), spin j = 3/2 *)
jj = 3/2; ms = Range[-jj, jj]; dd = Length[ms];
jz = N@DiagonalMatrix[ms];
jp = N@Normal@SparseArray[Table[{i + 1, i} -> Sqrt[(jj - ms[[i]]) (jj + ms[[i]] + 1)], {i, dd - 1}], {dd, dd}];
jm = Transpose[jp];
put["C038-commutators", {mx[jz . jp - jp . jz - jp], mx[jp . jm - jm . jp - 2 jz],
   mx[jz . jz + (jp . jm + jm . jp)/2 - jj (jj + 1) IdentityMatrix[dd]]}];
With[{nv = Normalize[{0.3, -0.5, 0.8}], t = 1.3},
  With[{aa = I nv . PauliMatrix /@ {1, 2, 3}},
   put["C039-su2exp", {Det[aa], mx[MatrixExp[t aa] - (Cos[t] IdentityMatrix[2] + aa Sin[t])]}]]];
put["C040-Gauss", Simplify[{{1, \[Beta]/\[Delta]}, {0, 1}} . {{1/\[Delta], 0}, {0, \[Delta]}} . {{1, 0}, {\[Gamma]/\[Delta], 1}} /.
    \[Gamma] -> (\[Alpha] \[Delta] - 1)/\[Beta]]];
With[{xi = 0.7 Exp[0.4 I]},
  With[{r = Abs[xi], ze = xi/Abs[xi] Tan[Abs[xi]]},
   put["C041-normal-form", mx[MatrixExp[xi jp - Conjugate[xi] jm] -
       MatrixExp[ze jp] . MatrixExp[-2 Log[Cos[r]] jz] . MatrixExp[-Conjugate[ze] jm]]];
   put["C041-lowest", mx[MatrixExp[xi jp - Conjugate[xi] jm] . UnitVector[dd, 1] -
       (1 + Abs[ze]^2)^(-jj) MatrixExp[ze jp] . UnitVector[dd, 1]]]]];
scs[z_] := (1 + Abs[z]^2)^(-jj) MatrixExp[z jp] . UnitVector[dd, 1];
nvec[z_] := {2 Re[z], 2 Im[z], Abs[z]^2 - 1}/(1 + Abs[z]^2);
With[{et = 0.4 - 0.9 I, ze = -0.3 + 0.5 I},
  put["C042-overlap", {Conjugate[scs[et]] . scs[ze] - (1 + Conjugate[et] ze)^(2 jj)/((1 + Abs[et]^2)^jj (1 + Abs[ze]^2)^jj),
     Abs[Conjugate[scs[et]] . scs[ze]]^2 - ((1 + nvec[et] . nvec[ze])/2)^(2 jj)}]];

(* C-045 Bergman realization *)
km[g_] := D[g, z]; kp[g_] := 2 k z g + z^2 D[g, z]; k0[g_] := k g + z D[g, z];
put["C045-algebra", Simplify[{km[kp[ff[z]]] - kp[km[ff[z]]] - 2 k0[ff[z]],
    k0[kp[ff[z]]] - kp[k0[ff[z]]] - kp[ff[z]], k0[km[ff[z]]] - km[k0[ff[z]]] + km[ff[z]],
    k0[k0[ff[z]]] - (kp[km[ff[z]]] + km[kp[ff[z]]])/2 - k (k - 1) ff[z]}]];

(* C-048 noncompact Gauss factorization in the 2x2 representation of su(1,1) *)
kp2 = {{0, 1}, {0, 0}}; km2 = {{0, 0}, {-1, 0}}; k02 = {{1/2, 0}, {0, -1/2}};
put["C048-2x2-algebra", {km2 . kp2 - kp2 . km2 - 2 k02, k02 . kp2 - kp2 . k02 - kp2, k02 . km2 - km2 . k02 + km2}];
put["C048-2x2-symbolic", Simplify[MatrixExp[r E^(I th) kp2 - r E^(-I th) km2] -
    MatrixExp[E^(I th) Tanh[r] kp2] . MatrixExp[-2 Log[Cosh[r]] k02] . MatrixExp[-E^(-I th) Tanh[r] km2], r > 0 && Element[th, Reals]]];
(* C-050 time-dependent SU(1,1), constant coefficients, 2x2 rep *)
With[{lam = 0.3 + 0.2 I, mu = 0.7, tf = 1.5},
  Module[{sol, om, ga, beta, uu},
   sol = NDSolve[{I om'[t] == Conjugate[lam] om[t]^2 + mu om[t] + lam, om[0] == 0,
       ga'[t] == -lam Conjugate[om[t]] - Conjugate[lam] om[t] - mu, ga[0] == 0}, {om, ga}, {t, 0, tf},
      PrecisionGoal -> 12, AccuracyGoal -> 12];
   beta[t_] := With[{w = om[t] /. First[sol]}, w/Abs[w] ArcTanh[Abs[w]]];
   uu[t_] := MatrixExp[-I t (lam kp2 + Conjugate[lam] km2 + mu k02)];
   put["C050-oneTime", Table[mx[uu[t] - MatrixExp[beta[t] kp2 - Conjugate[beta[t]] km2] . MatrixExp[I (ga[t] /. First[sol]) k02]], {t, {0.5, 1.0, 1.5}}]];
   put["C050-twoTime-corrected", mx[uu[1.5] . Inverse[uu[0.5]] -
       MatrixExp[beta[1.5] kp2 - Conjugate[beta[1.5]] km2] . MatrixExp[I ((ga[1.5] - ga[0.5]) /. First[sol]) k02] .
        MatrixExp[-(beta[0.5] kp2 - Conjugate[beta[0.5]] km2)]]];
   put["C050-twoTime-asPrinted", mx[uu[1.5] . Inverse[uu[0.5]] -
       MatrixExp[beta[1.5] kp2 - Conjugate[beta[1.5]] km2] . MatrixExp[I ((ga[1.5] - ga[0.5]) /. First[sol]) k02] .
        MatrixExp[-(beta[1.5] kp2 - Conjugate[beta[1.5]] km2)]]]]];

(* C-051 Duhamel *)
With[{l0 = {{0.2, 0.5, -0.1}, {0.3, -0.4, 0.2}, {0.1, 0.6, 0.3}}, l1 = {{-0.3, 0.1, 0.4}, {0.2, 0.5, -0.6}, {0.7, -0.2, 0.1}}, t0 = 0.3, h = 10.^-5},
  put["C051-Duhamel", mx[(MatrixExp[l0 + (t0 + h) l1] - MatrixExp[l0 + (t0 - h) l1])/(2 h) -
      Table[NIntegrate[(MatrixExp[s (l0 + t0 l1)] . l1 . MatrixExp[(1 - s) (l0 + t0 l1)])[[i, j]], {s, 0, 1}], {i, 3}, {j, 3}]]]];

(* C-052 singular oscillator lowest weight *)
With[{s = 1/2 + Sqrt[1/4 + g^2]},
  put["C052-lowestWeight", Simplify[(-D[x^s Exp[-x^2/2], {x, 2}]/2 + (x^2 + g^2/x^2) x^s Exp[-x^2/2]/2)/(x^s Exp[-x^2/2]) - (2 s + 1)/2, x > 0 && g > 0]]];

(* ---------- fermions: Jordan-Wigner ---------- *)
sm = {{0, 1}, {0, 0}}; zz = {{1, 0}, {0, -1}}; i2 = IdentityMatrix[2];
cOp[j_, n_] := KroneckerProduct @@ Join[ConstantArray[zz, j - 1], {sm}, ConstantArray[i2, n - j]];
dag = ConjugateTranspose;
acomm[x_, y_] := x . y + y . x;

With[{cs = Table[cOp[j, 3], {j, 3}]},
  put["C064-CAR", {Max[Table[mx[acomm[cs[[i]], dag[cs[[j]]]] - KroneckerDelta[i, j] IdentityMatrix[8]], {i, 3}, {j, 3}]],
     Max[Table[mx[acomm[cs[[i]], cs[[j]]]], {i, 3}, {j, 3}]]}]];

With[{cs = Table[cOp[j, 2], {j, 2}], u = {0.3 + 0.2 I, -0.5 + 0.1 I}, v = {0.4 - 0.3 I, 0.2 + 0.6 I}},
  With[{ll = Sum[u[[j]] dag[cs[[j]]] + v[[j]] cs[[j]], {j, 2}], uv = u . v},
   put["C065-Lsquared", mx[ll . ll - uv IdentityMatrix[4]]];
   put["C065-exp", mx[MatrixExp[ll] - (Cosh[Sqrt[uv]] IdentityMatrix[4] + Sinh[Sqrt[uv]]/Sqrt[uv] ll)]]]];

With[{cs = Table[cOp[j, 3], {j, 3}], h = {{0.3 + 0.1 I, 0.2, -0.4 I}, {0.5, -0.2, 0.1 + 0.3 I}, {0.1 I, 0.6, 0.4}}},
  With[{hq = Sum[h[[i, j]] dag[cs[[i]]] . cs[[j]], {i, 3}, {j, 3}], par = KroneckerProduct[zz, zz, zz]},
   put["C079-Tr", Tr[MatrixExp[hq]] - Det[IdentityMatrix[3] + MatrixExp[h]]];
   put["C079-Str", Tr[par . MatrixExp[hq]] - Det[IdentityMatrix[3] - MatrixExp[h]]]]];

(* C-080, C-082, C-083 BdG, two modes *)
With[{cs = Table[cOp[j, 2], {j, 2}],
      mm = {{0.7, 0.2 - 0.3 I}, {0.2 + 0.3 I, -0.4}}, lL = {{0, 0.5 + 0.2 I}, {-0.5 - 0.2 I, 0}}},
  Module[{hh, sig, cvecF, uu, fF2, aa, bb, gam, omega, hd = 10.^-5},
   hh = (Sum[mm[[i, j]] dag[cs[[i]]] . cs[[j]] + lL[[i, j]] dag[cs[[i]]] . dag[cs[[j]]] - Conjugate[lL[[i, j]]] cs[[i]] . cs[[j]], {i, 2}, {j, 2}])/2 -
     Tr[mm]/4 IdentityMatrix[4];
   sig = ArrayFlatten[{{mm/2, lL}, {-Conjugate[lL], -Conjugate[mm]/2}}];
   cvecF = Join[cs, dag /@ cs];
   uu[t_] := MatrixExp[-I t hh];
   fF2[t_] := MatrixExp[I t sig];
   put["C080-hermitian", mx[hh - dag[hh]]];
   put["C080-heisenberg", Max[Table[mx[uu[0.9] . cvecF[[k]] . dag[uu[0.9]] - Sum[fF2[0.9][[k, j]] cvecF[[j]], {j, 4}]], {k, 4}]]];
   aa[t_] := fF2[t][[1 ;; 2, 1 ;; 2]]; bb[t_] := fF2[t][[1 ;; 2, 3 ;; 4]];
   gam[t_] := Transpose[bb[t]] . Inverse[Transpose[aa[t]]];
   omega = UnitVector[4, 1];
   put["C082-antisym", mx[gam[0.6] + Transpose[gam[0.6]]]];
   put["C082-pairVacuum", mx[uu[0.6] . omega - Sqrt[Det[aa[0.6]]] MatrixExp[Sum[gam[0.6][[i, j]] dag[cs[[i]]] . dag[cs[[j]]], {i, 2}, {j, 2}]/2] . omega]];
   put["C083-Riccati", mx[I (gam[0.6 + hd] - gam[0.6 - hd])/(2 hd) -
       ((mm . gam[0.6] + gam[0.6] . Conjugate[mm])/2 + lL + gam[0.6] . Conjugate[lL] . gam[0.6])]];
   put["C083-amplitude", {omega . uu[0.6] . omega - Sqrt[Det[aa[0.6]]],
      I (Sqrt[Det[aa[0.6 + hd]]] - Sqrt[Det[aa[0.6 - hd]]])/(2 hd) + Tr[Conjugate[lL] . gam[0.6] + mm/2]/2 Sqrt[Det[aa[0.6]]]}]]];

(* C-087 fermion pair rotation *)
With[{cs = Table[cOp[j, 2], {j, 2}], xi = 0.9 Exp[0.3 I]},
  With[{jp2 = dag[cs[[1]]] . dag[cs[[2]]], jm2 = cs[[2]] . cs[[1]],
        j32 = (dag[cs[[1]]] . cs[[1]] + dag[cs[[2]]] . cs[[2]] - IdentityMatrix[4])/2},
   put["C087-su2", {mx[jp2 . jm2 - jm2 . jp2 - 2 j32], mx[j32 . jp2 - jp2 . j32 - jp2]}];
   put["C087-vacuum", mx[MatrixExp[xi jp2 - Conjugate[xi] jm2] . UnitVector[4, 1] -
       (Cos[Abs[xi]] UnitVector[4, 1] + xi/Abs[xi] Sin[Abs[xi]] jp2 . UnitVector[4, 1])]]]];

(* C-071 / C-072 / C-075: Grassmann labels realized as extra anticommuting creators *)
With[{c = cOp[1, 5], g = dag[cOp[2, 5]], gs = dag[cOp[3, 5]], al = dag[cOp[4, 5]], als = dag[cOp[5, 5]]},
  Block[{tt, cd = dag[c]}, tt[x_, xs_] := MatrixExp[cd . x - xs . c];
   put["C071-BCH", mx[tt[g, gs] - MatrixExp[-gs . g/2] . MatrixExp[cd . g] . MatrixExp[-gs . c]]];
   put["C071-product", mx[tt[al, als] . tt[g, gs] - tt[al + g, als + gs] . MatrixExp[-(als . g + al . gs)/2]]];
   put["C071-shift", {mx[Inverse[tt[g, gs]] . c . tt[g, gs] - (c + g)], mx[Inverse[tt[g, gs]] . cd . tt[g, gs] - (cd + gs)],
      mx[Inverse[tt[g, gs]] - tt[-g, -gs]]}];
   put["C075-normal-antinormal", {mx[MatrixExp[cd . g] . MatrixExp[-gs . c] - tt[g, gs] . MatrixExp[gs . g/2]],
      mx[MatrixExp[-gs . c] . MatrixExp[cd . g] - tt[g, gs] . MatrixExp[-gs . g/2]]}];
   put["C072-eigen", mx[c . MatrixExp[cd . g] . UnitVector[32, 1] - g . MatrixExp[cd . g] . UnitVector[32, 1]]]]];

(* ---------- Grassmann algebra ---------- *)
gclean[x_] := Select[x, ! (NumericQ[#] && Abs[#] < 10^-12) &];
gm[x_, y_] := gclean@Merge[Flatten[KeyValueMap[Function[{k1, v1}, KeyValueMap[Function[{k2, v2},
        If[IntersectingQ[k1, k2], Nothing, Sort[Join[k1, k2]] -> Signature[Join[k1, k2]] v1 v2]], y]], x]], Total];
gadd[xs___] := gclean@Merge[{xs}, Total];
one = <|{} -> 1|>;
gen[i_] := <|{i} -> 1|>;
gexp[x_] := With[{b = Lookup[x, Key[{}], 0], nil = KeyDrop[x, {{}}]},
   Exp[b] (gadd @@ MapIndexed[#1/((First[#2] - 1)!) &, NestList[gm[#, nil] &, one, 8]])];
deriv[i_, x_] := gclean@Merge[KeyValueMap[Function[{kk2, v}, With[{p = FirstPosition[kk2, i]},
       If[MissingQ[p], Nothing, DeleteCases[kk2, i] -> (-1)^(First[p] - 1) v]]], x], Total];
rderiv[i_, x_] := gclean@Merge[KeyValueMap[Function[{kk2, v}, With[{p = FirstPosition[kk2, i]},
       If[MissingQ[p], Nothing, DeleteCases[kk2, i] -> (-1)^(Length[kk2] - First[p]) v]]], x], Total];
pairedInt[pairs_, x_] := Fold[deriv[#2, #1] &, x, Reverse[Flatten[Reverse /@ pairs]]];
gdiff[x_, y_] := Max[Append[Abs[N[Values[gadd[x, -1 y]]]], 0]];
glin[coef_, idx_] := gadd @@ MapThread[#1 gen[#2] &, {coef, idx}];
gbil[l_, m_, r_] := gadd @@ Flatten[Table[m[[i, j]] gm[gen[l[[i]]], gen[r[[j]]]], {i, Length[l]}, {j, Length[r]}]];

(* C-066 Leibniz and left/right derivative *)
With[{f = gadd[0.3 gen[1], -0.7 gm[gm[gen[1], gen[2]], gen[3]], 0.2 gen[4]], g = gadd[one, 0.5 gm[gen[1], gen[3]], 0.4 gen[2]],
      fe = gadd[0.6 gm[gen[1], gen[2]], -0.3 gm[gen[1], gen[4]], 0.9 one]},
  put["C066-Leibniz-oddf", gdiff[deriv[1, gm[f, g]], gadd[gm[deriv[1, f], g], -1 gm[f, deriv[1, g]]]]];
  put["C066-right-left", {gdiff[rderiv[1, f], deriv[1, f]], gdiff[rderiv[1, fe], -1 deriv[1, fe]]}]];

(* C-067 Berezin determinant rule and paired normalization *)
With[{am = {{0.4, -0.2, 0.7}, {0.1, 0.5, -0.3}, {0.6, 0.2, 0.8}}},
  With[{at = Table[glin[am[[i]], {1, 2, 3}], {i, 3}]},
   With[{fA = gadd[gm[gm[at[[1]], at[[2]]], at[[3]]], 0.3 gm[at[[1]], at[[2]]], 0.5 one]},
    put["C067-detRule", Lookup[Fold[deriv[#2, #1] &, fA, {3, 2, 1}], Key[{}], 0] - Det[am] *
       Lookup[Fold[deriv[#2, #1] &, gadd[gm[gm[gen[1], gen[2]], gen[3]], 0.3 gm[gen[1], gen[2]], 0.5 one], {3, 2, 1}], Key[{}], 0]]]]];
put["C067-pairedGaussian", pairedInt[{{1, 3}, {2, 4}}, gexp[-1 gbil[{3, 4}, IdentityMatrix[2], {1, 2}]]]];

(* C-068 fermionic Gaussian with sources, n = 2: theta = 1,2; theta* = 3,4; gamma = 5,6; gamma* = 7,8 *)
With[{bm = {{0.8 + 0.2 I, -0.3}, {0.5 I, 0.6 - 0.4 I}}},
  With[{lhs = pairedInt[{{1, 3}, {2, 4}}, gexp[gadd[-1 gbil[{3, 4}, bm, {1, 2}], gbil[{7, 8}, IdentityMatrix[2], {1, 2}],
         gbil[{3, 4}, IdentityMatrix[2], {5, 6}]]]],
        rhs = Det[bm] gexp[gbil[{7, 8}, Inverse[bm], {5, 6}]]},
   put["C068-sources", gdiff[lhs, rhs]];
   put["C068-noSource", Lookup[pairedInt[{{1, 3}, {2, 4}}, gexp[-1 gbil[{3, 4}, bm, {1, 2}]]], Key[{}], 0] - Det[bm]]]];

(* C-069 Pfaffian Gaussian, n = 2: gamma = 1,2; gamma* = 3,4; xi = 5,6; xi* = 7,8 *)
pf4[m_] := m[[1, 2]] m[[3, 4]] - m[[1, 3]] m[[2, 4]] + m[[1, 4]] m[[2, 3]];
With[{kK = {{0, 0.4 + 0.1 I}, {-0.4 - 0.1 I, 0}}, lL = {{0, -0.3 + 0.2 I}, {0.3 - 0.2 I, 0}}, mM = {{0.9, 0.2 I}, {-0.1, 0.7 + 0.3 I}}},
  With[{lam = ArrayFlatten[{{kK, -Transpose[mM]}, {mM, lL}}]},
   With[{phi = gbil[{1, 2, 3, 4}, lam/2, {1, 2, 3, 4}],
         src = gadd[gbil[{3, 4}, IdentityMatrix[2], {5, 6}], -1 gbil[{7, 8}, IdentityMatrix[2], {1, 2}]]},
    With[{lhs = pairedInt[{{1, 3}, {2, 4}}, gexp[gadd[phi, src]]],
          rhs = (-1)^(2 (2 - 1)/2) pf4[lam] gexp[gbil[{7, 8, 5, 6}, Inverse[lam]/2, {7, 8, 5, 6}]]},
     put["C069-Pfaffian", gdiff[lhs, rhs]];
     put["C069-PfSquared", pf4[lam]^2 - Det[lam]]]]]];
With[{lam0 = ArrayFlatten[{{0, IdentityMatrix[2]}, {-IdentityMatrix[2], 0}}]},
  put["C069-example", {pf4[lam0], pairedInt[{{1, 3}, {2, 4}}, gexp[-1 gbil[{3, 4}, IdentityMatrix[2], {1, 2}]]]}]];

(* C-098 Grassmann Fourier transform *)
relabel[x_, rules_] := Association[KeyValueMap[Function[{kk2, v}, With[{nk = kk2 /. rules}, Sort[nk] -> Signature[nk] v]], x]];
four1[f_] := pairedInt[{{1, 2}}, gm[gexp[gadd[gm[gen[2], gen[3]], -1 gm[gen[4], gen[1]]]], f]];
inv1[ff_] := pairedInt[{{3, 4}}, gm[gexp[gadd[gm[gen[4], gen[1]], -1 gm[gen[2], gen[3]]]], ff]];
mon1 = {one, gen[1], gen[2], gm[gen[1], gen[2]]};
put["C098-table-n1", four1 /@ mon1];
put["C098-inverse-n1", gdiff[inv1[four1[#]], #] & /@ mon1];
put["C098-square-n1", gdiff[relabel[four1[relabel[four1[#], {3 -> 1, 4 -> 2}]], {3 -> 1, 4 -> 2}], #] & /@ mon1];
four2[f_] := pairedInt[{{1, 3}, {2, 4}}, gm[gexp[gadd[gbil[{3, 4}, IdentityMatrix[2], {5, 6}], -1 gbil[{7, 8}, IdentityMatrix[2], {1, 2}]]], f]];
mon2 = Fold[gm, one, gen /@ #] & /@ Subsets[{1, 2, 3, 4}];
put["C098-F1-n2", four2[one]];
put["C098-square-n2", Max[gdiff[relabel[four2[relabel[four2[#], {5 -> 1, 6 -> 2, 7 -> 3, 8 -> 4}]], {5 -> 1, 6 -> 2, 7 -> 3, 8 -> 4}], #] & /@ mon2]];
put["C098-inverse-n2", Max[gdiff[pairedInt[{{5, 7}, {6, 8}}, gm[gexp[gadd[gbil[{7, 8}, IdentityMatrix[2], {1, 2}], -1 gbil[{3, 4}, IdentityMatrix[2], {5, 6}]]], four2[#]]], #] & /@ mon2]];

(* C-099 holomorphic reproducing kernel, n = 2: theta = 1,2; gamma = 3,4; gamma* = 5,6 *)
eK = gm[gadd[one, -1 gm[gadd[gen[1], -1 gen[3]], gen[5]]], gadd[one, -1 gm[gadd[gen[2], -1 gen[4]], gen[6]]]];
put["C099-reproduce", Table[gdiff[pairedInt[{{3, 5}, {4, 6}}, gm[Fold[gm, one, gen /@ s], eK]], Fold[gm, one, gen /@ (s /. {3 -> 1, 4 -> 2})]],
   {s, {{}, {3}, {4}, {3, 4}}}]];


(* ---------- new exact checks replacing truncated-matrix checks ---------- *)

(* C-010: coherent resolution of identity, diagonal Fock element; off-diagonal vanish by the angular integral *)
put["C010-radial", FullSimplify[Integrate[2 r^(2 n + 1) Exp[-r^2], {r, 0, Infinity}, Assumptions -> n >= 0]/n!, n >= 0]];

(* C-012: Mehler kernel, exact Taylor coefficients in w through order 10 *)
hf[k_, x_] := (2^k k! Sqrt[Pi])^(-1/2) HermiteH[k, x] Exp[-x^2/2];
put["C012-Mehler-series", Simplify[Normal[Series[Pi^(-1/2) (1 - w^2)^(-1/2) Exp[-(1 + w^2)/(2 (1 - w^2)) (x^2 + y^2) + 2 w x y/(1 - w^2)], {w, 0, 10}]] -
    Sum[w^k hf[k, x] hf[k, y], {k, 0, 10}]]];
put["C012-tau-form", FullSimplify[(Exp[-tau/2] Pi^(-1/2) (1 - w^2)^(-1/2) Exp[-(1 + w^2)/(2 (1 - w^2)) (x^2 + y^2) + 2 w x y/(1 - w^2)] /. w -> Exp[-tau]) /
    ((2 Pi Sinh[tau])^(-1/2) Exp[-(x^2 + y^2)/(2 Tanh[tau]) + x y/Sinh[tau]]), tau > 0]];

(* C-022 / C-054: one-mode quadratic propagator U = exp(-i t H), H = (s11 Q^2 + s12 (QP+PQ) + s22 P^2)/2.
   H in su(1,1): checked by QF normal ordering; U disentangled exactly in the 2x2 representation,
   then coherent matrix elements of e^{c K+} e^{d K0} e^{f K-} are closed forms. *)
{a, ad} = FieldVariables[{1}];
qO = (a + ad)/Sqrt[2]; pO = -I (a - ad)/Sqrt[2];
put["C022-H-in-K", BosonicNormalOrder[(s11 qO ** qO + s12 (qO ** pO + pO ** qO) + s22 pO ** pO)/2 -
     ((s11 + 2 I s12 - s22)/2 ad ** ad/2 + (s11 - 2 I s12 - s22)/2 a ** a/2 + (s11 + s22) (ad ** a + 1/2)/2), {a, ad}]];
ss = {{9/10, 3/10}, {3/10, 3/5}}; tt = 4/5; jJ = {{0, 1}, {-1, 0}};
{s11v, s12v, s22v} = {ss[[1, 1]], ss[[1, 2]], ss[[2, 2]]};
gpath[t_] := MatrixExp[-I t ((s11v + 2 I s12v - s22v)/2 kp2 + (s11v - 2 I s12v - s22v)/2 km2 + (s11v + s22v) k02)];
put["C022-g22-path-arg", Max[Abs[Arg[N[gpath[#][[2, 2]]] & /@ Range[0, tt, tt/40]]]]];
gg = N[gpath[tt], 40];
cK = gg[[1, 2]]/gg[[2, 2]]; fK = -gg[[2, 1]]/gg[[2, 2]];
meGauss[al_, be_] := Exp[cK Conjugate[al]^2/2 + fK be^2/2] gg[[2, 2]]^(-1/2) Exp[-(Abs[al]^2 + Abs[be]^2)/2 + Conjugate[al] be/gg[[2, 2]]];
fF = N[MatrixExp[tt jJ . ss], 40];
{fa, fb, fc, fd} = {fF[[1, 1]], fF[[1, 2]], fF[[2, 1]], fF[[2, 2]]};
zal[z_] := (z[[1]] + I z[[2]])/Sqrt[2];
sig[z_, zp_] := z[[1]] zp[[2]] - z[[2]] zp[[1]];
atlas022[z_, zp_] := Module[{yY = fa + fd + I (fb - fc), zZ = fa - fd + I (fb + fc), wW, xx, xi},
   wW = zZ/yY; {xx, xi} = z - fF . zp;
   Sqrt[2] yY^(-1/2) Exp[-I sig[fF . zp, z]/2] Exp[-(xx^2 + xi^2)/4] Exp[-(xi + I xx) wW (xi + I xx)/4]];
pts = {{{3/10, -1/5}, {1/10, 3/5}}, {{-1/2, 2/5}, {2/5, 0}}, {{1, 1/2}, {-7/10, 3/10}}, {{0, 0}, {0, 0}}};
put["C022-coherentME", mx[(meGauss[zal[#[[1]]], zal[#[[2]]]] - atlas022[#[[1]], #[[2]]]) & /@ pts]];
put["C022-sample", N[{meGauss[zal[pts[[1, 1]]], zal[pts[[1, 2]]]], atlas022[pts[[1, 1]], pts[[1, 2]]]}, 12]];
put["C022-sensitivity-phase-flipped", mx[(meGauss[zal[#[[1]]], zal[#[[2]]]] - atlas022[#[[1]], #[[2]]] Exp[I sig[fF . #[[2]], #[[1]]]]) & /@ pts]];
(* Herman-Kluk: (C_t / 2 pi) Integral <alpha|phi_{F z}> <phi_z|beta> dq dp, Gaussian integral in closed form *)
ovl[x_, xc_, y_, yc_] := Exp[-(x xc + y yc)/2 + xc y];
hkME[al_, be_] := Module[{z = {q, p}, gam, gamc, zt, ztc, ex, mM, jv, c0, ev},
   gam = (fF . z)[[1]]/Sqrt[2] + I (fF . z)[[2]]/Sqrt[2]; gamc = (fF . z)[[1]]/Sqrt[2] - I (fF . z)[[2]]/Sqrt[2];
   zt = (q + I p)/Sqrt[2]; ztc = (q - I p)/Sqrt[2];
   ex = Expand[-(al Conjugate[al] + gam gamc)/2 + Conjugate[al] gam - (zt ztc + be Conjugate[be])/2 + ztc be];
   mM = -D[ex, {z, 2}]; jv = D[ex, {z}] /. Thread[z -> 0]; c0 = ex /. Thread[z -> 0];
   ev = Eigenvalues[mM];
   {Min[Eigenvalues[Re[mM]]], Sqrt[(fa + fd + I (fc - fb))/2]/(2 Pi) 2 Pi Times @@ (ev^(-1/2)) Exp[jv . Inverse[mM] . jv/2 + c0]}];
put["C054-HK-ReM-min", Min[hkME[zal[#[[1]]], zal[#[[2]]]][[1]] & /@ pts]];
put["C054-HK-sample", N[hkME[zal[pts[[2, 1]]], zal[pts[[2, 2]]]][[2]], 12]];
put["C054-HK-vs-exact", mx[(hkME[zal[#[[1]]], zal[#[[2]]]][[2]] - meGauss[zal[#[[1]]], zal[#[[2]]]]) & /@ pts]];

(* C-032: driven oscillator H = N + 1/2 + f Q in the Bargmann realization a -> d/dz, a^dag -> z.
   psi(t,z) = e^{i delta_t} [D(alpha_t) R_t h](z) for arbitrary analytic h must solve i d_t psi = H psi. *)
qt[t_] := -f + (q0 + f) Cos[t] + p0 Sin[t];
pt[t_] := -(q0 + f) Sin[t] + p0 Cos[t];
put["C032-classical-orbit", Simplify[{D[qt[t], t] - pt[t], D[pt[t], t] + qt[t] + f}]];
del[t_] = Integrate[pt[s] D[qt[s], s] - ((pt[s]^2 + qt[s]^2)/2 + f qt[s]), {s, 0, t}] - (qt[t] pt[t] - q0 p0)/2;
alt[t_] := (qt[t] + I pt[t])/Sqrt[2]; altc[t_] := (qt[t] - I pt[t])/Sqrt[2];
psi[t_, z_] := Exp[I del[t]] Exp[-alt[t] altc[t]/2] Exp[alt[t] z] Exp[-I t/2] hh[Exp[-I t] (z - altc[t])];
hamB[ps_] := z D[ps, z] + ps/2 + f (z ps + D[ps, z])/Sqrt[2];
put["C032-Schroedinger", FullSimplify[(I D[psi[t, z], t] - hamB[psi[t, z]])/(Exp[I del[t]] Exp[-alt[t] altc[t]/2] Exp[alt[t] z] Exp[-I t/2])]];
put["C032-initial", Simplify[{del[0], qt[0] - q0, pt[0] - p0}]];

WriteString["stdout", ToString[InputForm[res]], "\n"];
