(* Exact Wolfram Language checks for the P-section entries that have no exact QF route.
   No truncated Fock matrix appears anywhere: finite matrices are spin, 2x2 group, 2N x 2N
   coefficient, or Jordan-Wigner fermion matrices, all with exact entries. *)
res = <||>;
put[k_, e_] := (res[k] = TimeConstrained[e, 300, "TIMEOUT"]);
SetAttributes[put, HoldRest];
cj = Conjugate; id[n_] := IdentityMatrix[n];

(* 13: coherent resolution of identity, radial integral for every n *)
put["E13radial", Simplify[Integrate[2 rr Exp[-rr^2] rr^(2 n)/n!, {rr, 0, Infinity}, Assumptions -> n >= 0]]];
(* 14: e^{z^2/2} has Bargmann coefficients c_{2k} = Sqrt[(2k)!]/(2^k k!); the norm series *)
put["E14normSeriesConverges", SumConvergence[(2 k)!/(4^k (k!)^2), k]];
(* 15: coordinate to Bargmann kernel, summed exactly *)
put["E15", {Sum[tg^nn HermiteH[nn, xh]/nn!, {nn, 0, Infinity}],
   PowerExpand[(Exp[-xh^2/2] Sum[tg^nn HermiteH[nn, xh]/nn!, {nn, 0, Infinity}] /. {tg -> z/Sqrt[2], xh -> q/Sqrt[hb]}) -
     Exp[-z^2/2 + Sqrt[2/hb] z q - q^2/(2 hb)]]}];
(* 19: smoothing the stated P symbol of D(gamma) gives the stated Q symbol *)
put["E19", Simplify[Integrate[Exp[-((x0 - x)^2 + (y0 - y)^2)] Exp[(g1^2 + g2^2)/2 + 2 I (g2 x - g1 y)],
       {x, -Infinity, Infinity}, {y, -Infinity, Infinity}, Assumptions -> Element[{x0, y0, g1, g2}, Reals]]/Pi /
     Exp[-(g1^2 + g2^2)/2 + 2 I (g2 x0 - g1 y0)]]];
(* 25: thermal P -> W with kernel 2 e^{-2|a-b|^2}, compared with W = Tr[T rho] summed exactly *)
put["E25PtoW", Simplify[Integrate[2 Exp[-2 ((x0 - x)^2 + y^2)] Exp[-(x^2 + y^2)/nb]/nb,
      {x, -Infinity, Infinity}, {y, -Infinity, Infinity}, Assumptions -> nb > 0 && Element[x0, Reals]]/Pi]];
put["E25WfromTrace", Simplify[2 Exp[-2 x0^2] Sum[(1 - xx) (-xx)^m LaguerreL[m, 4 x0^2], {m, 0, Infinity}] /.
     xx -> nb/(1 + nb), nb > 0]];

(* 29: su(2) holomorphic realization and exact spin matrices *)
jp2[f_] := -z^2 D[f, z] + 2 j z f; jm2[f_] := D[f, z]; j02[f_] := z D[f, z] - j f;
put["E29diff", Simplify[{j02[jp2[f[z]]] - jp2[j02[f[z]]] - jp2[f[z]], j02[jm2[f[z]]] - jm2[j02[f[z]]] + jm2[f[z]],
    jp2[jm2[f[z]]] - jm2[jp2[f[z]]] - 2 j02[f[z]]}]];
spin[jj_] := Module[{ms = Range[-jj, jj], jpm},
   jpm = Table[If[mr == mc + 1, Sqrt[(jj - mc) (jj + mc + 1)], 0], {mr, ms}, {mc, ms}];
   {jpm, Transpose[jpm], DiagonalMatrix[ms]}];
put["E29spin", Table[With[{s = spin[jj]}, Union[Flatten[{s[[3]] . s[[1]] - s[[1]] . s[[3]] - s[[1]],
      s[[3]] . s[[2]] - s[[2]] . s[[3]] + s[[2]], s[[1]] . s[[2]] - s[[2]] . s[[1]] - 2 s[[3]]}]]], {jj, {1/2, 1, 3/2, 2}}]];

(* 30: SU(2) disentangling in the faithful 2x2 matrices, symbolic r, phi *)
{p2, m2, z2} = {{{0, 1}, {0, 0}}, {{0, 0}, {1, 0}}, {{1/2, 0}, {0, -1/2}}};
lhs2 = MatrixExp[r Exp[I ph] p2 - r Exp[-I ph] m2];
zt = Exp[I ph] Tan[r]; ztb = Exp[-I ph] Tan[r]; lg = Log[1 + Tan[r]^2];
mexp[c_, m_] := MatrixExp[c m];
nrm2 = mexp[zt, p2] . mexp[lg, z2] . mexp[-ztb, m2];
ant2 = mexp[-ztb, m2] . mexp[-lg, z2] . mexp[zt, p2];
put["E30chart", FullSimplify[{nrm2 - lhs2, ant2 - lhs2}, 0 < r < Pi/2 && Element[ph, Reals]]];
put["E30beyond", FullSimplify[{nrm2 + lhs2, ant2 + lhs2}, Pi/2 < r < Pi && Element[ph, Reals]]];

(* 31, 32: SU(2) coherent state, overlap, annihilator for j = 3/2 and 2, symbolic zeta, eta *)
cs2[jj_, x_] := (1 + x xb)^(-jj) MatrixExp[x spin[jj][[1]]] . UnitVector[2 jj + 1, 1];
put["E31", Table[Simplify[cs2[jj, ze] - (1 + ze xb)^(-jj) Table[Sqrt[Binomial[2 jj, nn]] ze^nn, {nn, 0, 2 jj}]],
    {jj, {3/2, 2}}]];
put["E32overlap", Table[Simplify[(MatrixExp[ze spin[jj][[1]]] . UnitVector[2 jj + 1, 1] /. ze -> zb) .
       (MatrixExp[et spin[jj][[1]]] . UnitVector[2 jj + 1, 1]) - (1 + zb et)^(2 jj)], {jj, {3/2, 2}}]];
put["E32annihilator", Table[With[{s = spin[jj]}, Simplify[(s[[2]] + 2 ze s[[3]] - ze^2 s[[1]]) .
       MatrixExp[ze s[[1]]] . UnitVector[2 jj + 1, 1]]], {jj, {3/2, 2}}]];
(* 33: SU(2) resolution of identity, diagonal element for all n, j *)
put["E33", FullSimplify[(2 jj + 1) Binomial[2 jj, nn] Integrate[tt^nn (1 + tt)^(-2 jj - 2), {tt, 0, Infinity},
      Assumptions -> nn >= 0 && 2 jj - nn >= 0], Element[nn, Integers] && 0 <= nn <= 2 jj]];

(* 34: su(1,1): holomorphic realization and the ladder coefficients, symbolic k and n *)
kp1[f_] := z^2 D[f, z] + 2 k z f; km1[f_] := D[f, z]; k01[f_] := z D[f, z] + k f;
put["E34diff", Simplify[{k01[kp1[f[z]]] - kp1[k01[f[z]]] - kp1[f[z]], km1[kp1[f[z]]] - kp1[km1[f[z]]] - 2 k01[f[z]],
    k01[k01[f[z]]] - (kp1[km1[f[z]]] + km1[kp1[f[z]]])/2 - k (k - 1) f[z]}]];
up[nn_] := Sqrt[(nn + 1) (2 k + nn)]; dn[nn_] := Sqrt[nn (2 k + nn - 1)];
put["E34ladder", FullSimplify[{up[nn]^2 - dn[nn]^2 - 2 (k + nn) (* [K-,K+] on |n> *),
    (k + nn)^2 - (dn[nn]^2 + up[nn]^2)/2 - k (k - 1) (* Casimir on |n> *)}, k > 0 && nn >= 0]];
(* 35: SU(1,1) disentangling in the faithful 2x2 matrices, corrected and printed middle signs *)
{p3, m3, z3} = {{{0, 1}, {0, 0}}, {{0, 0}, {-1, 0}}, {{1/2, 0}, {0, -1/2}}};
put["E35comm", {z3 . p3 - p3 . z3 - p3, m3 . p3 - p3 . m3 - 2 z3}];
lhs3 = MatrixExp[r Exp[I ph] p3 - r Exp[-I ph] m3];
zh = Exp[I ph] Tanh[r]; zhb = Exp[-I ph] Tanh[r]; lh = Log[1 - Tanh[r]^2];
put["E35corrected", FullSimplify[{mexp[zh, p3] . mexp[lh, z3] . mexp[-zhb, m3] - lhs3,
     mexp[-zhb, m3] . mexp[-lh, z3] . mexp[zh, p3] - lhs3}, r > 0 && Element[ph, Reals]]];
put["E35printed", FullSimplify[{mexp[zh, p3] . mexp[-lh, z3] . mexp[-zhb, m3] - lhs3,
     mexp[-zhb, m3] . mexp[lh, z3] . mexp[zh, p3] - lhs3}, r > 0 && Element[ph, Reals]]];
(* 36-38: SU(1,1) coherent states, symbolic k *)
put["E36coeff", FullSimplify[Product[up[jx], {jx, 0, nn - 1}]/nn! - Sqrt[Gamma[nn + 2 k]/(nn! Gamma[2 k])],
    k > 0 && Element[nn, Integers] && nn >= 0]];
put["E36norm", Simplify[Sum[Gamma[nn + 2 k]/(nn! Gamma[2 k]) x^nn, {nn, 0, Infinity}], 0 < x < 1 && k > 0]];
cc[nn_] := ze^nn Sqrt[Gamma[nn + 2 k]/(nn! Gamma[2 k])];
(* coefficient of |k,n> in (K- - 2 ze K0 + ze^2 K+)|ze;k>, divided by cc[n]; ratios of neighbours from Gamma[x+1] = x Gamma[x] *)
put["E37annihilator", {FullSimplify[FunctionExpand[(cc[nn + 1]/cc[nn])^2] - ze^2 (nn + 2 k)/(nn + 1), k > 0 && nn >= 1],
   FullSimplify[FunctionExpand[(cc[nn - 1]/cc[nn])^2] - nn/((nn + 2 k - 1) ze^2), k > 0 && nn >= 1],
   PowerExpand[Simplify[ze Sqrt[(nn + 2 k)/(nn + 1)] up[nn] - 2 ze (k + nn) + ze^2 Sqrt[nn/(nn + 2 k - 1)]/ze dn[nn], k > 0 && nn >= 1]]}];
put["E37overlap", Simplify[Sum[Gamma[nn + 2 k]/(nn! Gamma[2 k]) (zb et)^nn, {nn, 0, Infinity}], Abs[zb et] < 1 && k > 0]];
put["E38", FullSimplify[(2 k - 1) Gamma[nn + 2 k]/(nn! Gamma[2 k]) Integrate[tt^nn (1 - tt)^(2 k - 2), {tt, 0, 1},
      Assumptions -> nn >= 0 && Element[k, Reals]]]];

(* exact group elements by the Cayley transform of a Lie-algebra element *)
cay[l_] := (id[Length[l]] + l/2) . Inverse[id[Length[l]] - l/2];
(* 43: bosonic Bogoliubov, N = 2: L = {{C, A}, {Conj A, Conj C}}, C anti-Hermitian, A symmetric *)
cB = {{I/3, 1/2 + I/5}, {-1/2 + I/5, -I/7}}; aB = {{1/4 + I/3, -1/5 + I/2}, {-1/5 + I/2, 2/7 - I/6}};
mB = cay[ArrayFlatten[{{cB, aB}, {cj[aB], cj[cB]}}]];
kB = DiagonalMatrix[{1, 1, -1, -1}];
{uB, vB} = {mB[[1 ;; 2, 1 ;; 2]], mB[[1 ;; 2, 3 ;; 4]]};
zB = Inverse[uB] . vB;
put["E43", Union[Flatten[Simplify[{mB[[3 ;; 4, 1 ;; 2]] - cj[vB], mB[[3 ;; 4, 3 ;; 4]] - cj[uB],
      uB . ConjugateTranspose[uB] - vB . ConjugateTranspose[vB] - id[2], uB . Transpose[vB] - vB . Transpose[uB],
      mB . kB . ConjugateTranspose[mB] - kB,
      Inverse[mB] - ArrayFlatten[{{ConjugateTranspose[uB], -Transpose[vB]}, {-ConjugateTranspose[vB], Transpose[uB]}}],
      ConjugateTranspose[uB] . uB - Transpose[vB] . cj[vB] - id[2], ConjugateTranspose[uB] . vB - Transpose[vB] . cj[uB],
      zB - Transpose[zB], id[2] - zB . ConjugateTranspose[zB] - Inverse[ConjugateTranspose[uB] . uB]}]]]];
(* 44: L K + K L^dagger for a general block L *)
put["E44", Module[{c = Array[cc0, {2, 2}], aa = Array[aa0, {2, 2}], l},
   l = ArrayFlatten[{{c, aa}, {cj[aa], cj[c]}}];
   Simplify[(l . kB + kB . ConjugateTranspose[l]) - ArrayFlatten[{{c + ConjugateTranspose[c], -aa + Transpose[aa]},
        {cj[aa] - ConjugateTranspose[aa], -cj[c] - Transpose[c]}}]]]];
(* 48: fractional-linear action, composition and image *)
act[m_, zz_] := With[{u = m[[1 ;; 2, 1 ;; 2]], v = m[[1 ;; 2, 3 ;; 4]]}, Inverse[u + zz . cj[v]] . (v + zz . cj[u])];
mB2 = cay[ArrayFlatten[{{{{-I/5, 1/3}, {-1/3, I/2}}, {{1/6, I/4}, {I/4, -1/3}}}, {cj[{{1/6, I/4}, {I/4, -1/3}}], cj[{{-I/5, 1/3}, {-1/3, I/2}}]}}]];
z0 = {{1/5 + I/7, -1/6}, {-1/6, I/4}};
put["E48", With[{img = act[mB, z0]}, {Simplify[img - Transpose[img]],
    Simplify[act[mB, act[mB2, z0]] - act[mB2 . mB, z0]],
    Eigenvalues[id[2] - img . ConjugateTranspose[img]] // N[#, 20] &}]];
(* 49: bosonic bilinear symbols from the annihilator (a + Z a^dagger)|Z> = 0, exact rational Z *)
z49 = {{1/4 + I/5, 1/6 - I/7}, {1/6 - I/7, -1/5 + I/9}};
nn49 = Array[nv, {2, 2}];
sol49 = First[Solve[Flatten[Table[nv[kk, ii] == Sum[cj[z49[[kk, l]]] z49[[ii, m]] (KroneckerDelta[l, m] + nv[m, l]), {l, 2}, {m, 2}],
       {kk, 2}, {ii, 2}]], Flatten[nn49]]];
nEx = nn49 /. sol49; (* nEx[[k, i]] = <a_k^dagger a_i> *)
aEx = Table[-Sum[z49[[jj, kk]] (KroneckerDelta[ii, kk] + nEx[[kk, ii]]), {kk, 2}], {ii, 2}, {jj, 2}]; (* <a_i a_j> *)
aF = -z49 . Inverse[id[2] - ConjugateTranspose[z49] . z49];
bF = id[2]/2 + z49 . Inverse[id[2] - ConjugateTranspose[z49] . z49] . ConjugateTranspose[z49];
put["E49", {Simplify[aEx - aF], Simplify[Transpose[nEx] + id[2]/2 - bF], Simplify[nEx + id[2]/2 - bF],
   Simplify[{bF . bF - aF . cj[aF] - id[2]/4, bF . aF - aF . cj[bF], ConjugateTranspose[bF] - bF, Transpose[aF] - aF}]}];

(* Jordan-Wigner fermions, exact integer matrices *)
sm = {{0, 1}, {0, 0}}; sz = {{1, 0}, {0, -1}};
jw[nm_] := Table[KroneckerProduct @@ Join[Table[sz, i - 1], {sm}, Table[id[2], nm - i]], {i, nm}];
anti[x_, y_] := x . y + y . x; comm[x_, y_] := x . y - y . x;
vac[nm_] := UnitVector[2^nm, 1];
(* 50: fermionic Bogoliubov, N = 3, Cayley element of L = {{C, A}, {Conj A, Conj C}}, C anti-Hermitian, A skew *)
cF = {{I/3, 1/2, I/5}, {-1/2, -I/4, 1/3 + I/6}, {I/5, -1/3 + I/6, I/2}};
aFm = {{0, 1/3 + I/4, -1/5}, {-1/3 - I/4, 0, I/6}, {1/5, -I/6, 0}};
mF = cay[ArrayFlatten[{{cF, aFm}, {cj[aFm], cj[cF]}}]];
{uF, vF} = {mF[[1 ;; 3, 1 ;; 3]], mF[[1 ;; 3, 4 ;; 6]]};
c3 = jw[3]; cd3 = ConjugateTranspose /@ c3;
b3 = Table[Sum[uF[[i, j]] c3[[j]] + vF[[i, j]] cd3[[j]], {j, 3}], {i, 3}]; bd3 = ConjugateTranspose /@ b3;
put["E50", {Union[Flatten[Simplify[Table[{anti[b3[[i]], bd3[[j]]] - KroneckerDelta[i, j] id[8], anti[b3[[i]], b3[[j]]]}, {i, 3}, {j, 3}]]]],
   Union[Flatten[Simplify[{mF[[4 ;; 6, 1 ;; 3]] - cj[vF], uF . ConjugateTranspose[uF] + vF . ConjugateTranspose[vF] - id[3],
       uF . Transpose[vF] + vF . Transpose[uF], ConjugateTranspose[mF] . mF - id[6],
       Inverse[uF] . vF + Transpose[Inverse[uF] . vF]}]]]}];
(* 51: fermionic quadratic commutator table, two modes *)
c2 = jw[2]; cd2 = ConjugateTranspose /@ c2;
eF[i_, j_] := cd2[[i]] . c2[[j]] - KroneckerDelta[i, j] id[4]/2; fF[i_, j_] := c2[[i]] . c2[[j]]; gF[i_, j_] := cd2[[i]] . cd2[[j]];
dd = KroneckerDelta;
put["E51", Union[Flatten[Table[{comm[eF[i, j], eF[k, l]] - (dd[j, k] eF[i, l] - dd[i, l] eF[k, j]),
      comm[eF[i, j], fF[k, l]] - (-dd[i, k] fF[j, l] - dd[i, l] fF[k, j]),
      comm[eF[i, j], gF[k, l]] - (dd[j, k] gF[i, l] + dd[j, l] gF[k, i]),
      comm[fF[i, j], gF[k, l]] - (-dd[j, k] eF[l, i] + dd[i, k] eF[l, j] + dd[j, l] eF[k, i] - dd[i, l] eF[k, j]),
      comm[fF[i, j], fF[k, l]], comm[gF[i, j], gF[k, l]]}, {i, 2}, {j, 2}, {k, 2}, {l, 2}]]]];
(* 52, 54: fermionic pair-exponential state on four modes, exact rational skew Z, W *)
c4 = jw[4]; cd4 = ConjugateTranspose /@ c4;
skew[v_] := Module[{m = ConstantArray[0, {4, 4}], pr = Subsets[Range[4], {2}]},
   Scan[(m[[pr[[#, 1]], pr[[#, 2]]]] = v[[#]]; m[[pr[[#, 2]], pr[[#, 1]]]] = -v[[#]]) &, Range[6]]; m];
zS = skew[{1/3 + I/4, -1/5, I/6, 1/2 - I/3, 1/7, -1/4 + I/5}];
wS = skew[{-1/4, 1/6 + I/3, 1/5 - I/7, -I/3, 1/3 + 1/8 I, 1/6}];
ket4[zz_] := MatrixExp[-Sum[zz[[i, j]] cd4[[i]] . cd4[[j]], {i, 4}, {j, 4}]/2] . vac[4];
put["E52", {Simplify[(cj[ket4[zS]] . ket4[zS])^2 - Det[id[4] + zS . ConjugateTranspose[zS]]],
   Union[Flatten[Simplify[Table[(c4[[i]] + Sum[zS[[i, j]] cd4[[j]], {j, 4}]) . ket4[zS], {i, 4}]]]],
   Simplify[(cj[ket4[zS]] . ket4[wS])^2 - Det[id[4] + wS . ConjugateTranspose[zS]]],
   Simplify[cj[ket4[zS]] . ket4[wS] - Sqrt[Det[id[4] + wS . ConjugateTranspose[zS]]]]}];
put["E54", With[{nk = cj[ket4[zS]] . ket4[zS]}, With[{ev = Function[op, cj[ket4[zS]] . op . ket4[zS]/nk]},
    With[{aF4 = -zS . Inverse[id[4] + ConjugateTranspose[zS] . zS],
       bF4 = zS . Inverse[id[4] + ConjugateTranspose[zS] . zS] . ConjugateTranspose[zS] - id[4]/2},
     {Simplify[aF4 - Table[ev[c4[[j]] . c4[[i]]], {i, 4}, {j, 4}]],
      Simplify[aF4 - Table[ev[c4[[i]] . c4[[j]]], {i, 4}, {j, 4}]] // Max[Abs[N[#]]] &,
      Simplify[bF4 - Table[ev[cd4[[j]] . c4[[i]]] - dd[i, j]/2, {i, 4}, {j, 4}]],
      Union[Flatten[Simplify[{bF4 . bF4 - aF4 . cj[aF4] - id[4]/4, bF4 . aF4 - aF4 . cj[bF4], ConjugateTranspose[bF4] - bF4,
          Transpose[aF4] + aF4}]]]}]]]];
(* 56: log-determinant series, symbolic 3x3 to order t^5 *)
x3 = Array[xx3, {3, 3}];
put["E56", Simplify[Normal[Series[-Log[Det[id[3] - t x3]], {t, 0, 5}]] - Sum[t^m Tr[MatrixPower[x3, m]]/m, {m, 1, 5}]]];
(* 57: real symplectic 4x4 from the Cayley transform of J S, S real symmetric *)
jS = ArrayFlatten[{{0, id[2]}, {-id[2], 0}}];
mS = cay[jS . {{1/2, 1/3, 0, 1/5}, {1/3, -1/4, 1/6, 0}, {0, 1/6, 1/3, -1/7}, {1/5, 0, -1/7, 1/2}}];
{aS, bS, cS, dS} = {mS[[1 ;; 2, 1 ;; 2]], mS[[1 ;; 2, 3 ;; 4]], mS[[3 ;; 4, 1 ;; 2]], mS[[3 ;; 4, 3 ;; 4]]};
zSg = {{1/3 + I, 1/4}, {1/4, -1/5 + 2 I}};
put["E57", With[{img = (aS . zSg + bS) . Inverse[cS . zSg + dS]}, {Transpose[mS] . jS . mS - jS,
    Union[Flatten[{Transpose[aS] . cS - Transpose[cS] . aS, Transpose[bS] . dS - Transpose[dS] . bS,
       Transpose[aS] . dS - Transpose[cS] . bS - id[2]}]], Simplify[img - Transpose[img]], Eigenvalues[Im[img]] // N[#, 20] &}]];
(* 63: Riccati dynamics in the holomorphic realization; zeta, zetabar, phi functions of t *)
psi = Exp[-I ph[t]] (1 - ze[t] zb[t])^k (1 - ze[t] zz)^(-2 k);
kpH[f_] := zz^2 D[f, zz] + 2 k zz f; kmH[f_] := D[f, zz]; k0H[f_] := zz D[f, zz] + k f;
ham[f_] := I (be kpH[f] - bb kmH[f] - I gm k0H[f]);
flow = {ze'[t] -> be - bb ze[t]^2 - I gm ze[t], zb'[t] -> bb - be zb[t]^2 + I gm zb[t]};
put["E63", {Simplify[(I D[psi, t] - ham[psi])/psi /. flow /. ph'[t] -> -I k (bb ze[t] - be zb[t] + I gm)],
   Simplify[(I D[psi, t] - ham[psi])/psi /. flow /. ph'[t] -> -I k (-bb ze[t] + be zb[t] + I gm)]}];
(* 65: singular oscillator ground state *)
put["E65", Simplify[(-D[x^aa Exp[-x^2/2], {x, 2}]/2 + (x^2/2 + aa (aa - 1)/(2 x^2)) x^aa Exp[-x^2/2])/(x^aa Exp[-x^2/2])]];
(* 67: affine Gaussian norm by the Bargmann Gaussian integral, zeta = zr + I zi, alpha = ar + I ai *)
put["E67", Module[{w = x + I y, f, n2},
   f = Exp[(zr + I zi) w^2/2 + (ar + I ai) w];
   n2 = Integrate[Exp[-(x^2 + y^2)] ComplexExpand[f Conjugate[f]], {x, -Infinity, Infinity}, {y, -Infinity, Infinity},
      Assumptions -> Element[{zr, zi, ar, ai}, Reals] && zr^2 + zi^2 < 1]/Pi;
   {Simplify[n2 ((1 - zr^2 - zi^2)^(1/4) Exp[-(ar^2 + ai^2 + ComplexExpand[Re[(zr - I zi) (ar + I ai)^2]])/(2 (1 - zr^2 - zi^2))])^2,
      Element[{zr, zi, ar, ai}, Reals] && zr^2 + zi^2 < 1],
    Simplify[n2 ((1 - zr^2 - zi^2)^(-1/2) Exp[-(ar^2 + ai^2)/2])^2 /. {zi -> 0, ai -> 0}, -1 < zr < 1]}]];
(* 74, 76-78: fermionic pairs on modes a1, a2, b1, b2 (Jordan-Wigner order a1, a2, b1, b2) *)
{fa1, fa2, fb1, fb2} = c4; {fad1, fad2, fbd1, fbd2} = cd4;
put["E74", With[{jpp = fad1 . fbd1, jmm = fb1 . fa1, j0p = (fad1 . fa1 + fbd1 . fb1 - id[16])/2},
    {jpp . jmm - jmm . jpp - 2 j0p} . vac[4] // Simplify]];
pairKet[zz_] := MatrixExp[Sum[zz[[i, j]] {fad1, fad2}[[i]] . {fbd1, fbd2}[[j]], {i, 2}, {j, 2}]] . vac[4];
zP = {{1/2 + I/3, -1/4}, {I/5, 1/3 - I/6}}; wP = {{-1/3, 1/5 + I/2}, {1/7, -I/4}};
put["E76", Simplify[cj[pairKet[zP]] . pairKet[zP] - Det[id[2] + ConjugateTranspose[zP] . zP]]];
put["E77", Simplify[cj[pairKet[zP]] . pairKet[wP] - Det[id[2] + ConjugateTranspose[zP] . wP]]];
put["E78", Simplify[cj[pairKet[zP]] . MatrixPower[(id[16] + (tt - 1) fad1 . fa1) . (id[16] + (tt - 1) fad2 . fa2), 1] . pairKet[zP] /
      cj[pairKet[zP]] . pairKet[zP] - Det[id[2] + tt ConjugateTranspose[zP] . zP]/Det[id[2] + ConjugateTranspose[zP] . zP]]];
(* 81: attenuation solution solves the characteristic-function PDE; x, y stand for eta, etabar *)
chi = Exp[-(1 - Exp[-gm t]) (nu + sg) x y] f0[Sqrt[Exp[-gm t]] x, Sqrt[Exp[-gm t]] y];
put["E81", Simplify[D[chi, t] + gm ((x D[chi, x] + y D[chi, y])/2 + (nu + sg) x y chi)]];
(* 83: the corrected kernel is a normalized Gaussian and its Fourier transform is the decay factor *)
put["E83", {Integrate[Exp[-((x - Sqrt[p] u)^2 + (y - Sqrt[p] v)^2)/s]/(Pi s), {x, -Infinity, Infinity}, {y, -Infinity, Infinity},
     Assumptions -> s > 0 && p > 0 && Element[{u, v}, Reals]],
   Simplify[Integrate[Exp[I (k1 x + k2 y)] Exp[-(x^2 + y^2)/s]/(Pi s), {x, -Infinity, Infinity}, {y, -Infinity, Infinity},
      Assumptions -> s > 0 && Element[{k1, k2}, Reals]]]}];
(* 86: populations read off the displaced thermal kernel: rho_nn = n! [z1^n z2b^n] R *)
rK = Exp[(mu u1 u2 + bt u1 + btb u2 - bt btb)/(1 + mu)]/(1 + mu);
put["E86", Table[Simplify[nn! SeriesCoefficient[rK, {u1, 0, nn}, {u2, 0, nn}] -
      mu^nn/(1 + mu)^(nn + 1) Exp[-bt btb/(1 + mu)] LaguerreL[nn, -bt btb/(mu (1 + mu))], mu > 0], {nn, 0, 6}]];
(* 89: proper-time oscillator trace *)
put["E89", Simplify[Sum[Exp[-s (2 m + 1)], {m, 0, Infinity}, Assumptions -> s > 0] - 1/(2 Sinh[s]), s > 0]];
(* 90, 91: Laguerre generating functions *)
put["E90", Simplify[Normal[Series[(1 - t)^(-kk - 1) Exp[-t x/(1 - t)], {t, 0, 6}]] - Sum[t^nn LaguerreL[nn, kk, x], {nn, 0, 6}]]];
put["E91", Table[Simplify[Normal[Series[(1 + t)^mm Exp[-t x], {t, 0, 10}]] - Sum[t^nn LaguerreL[nn, mm - nn, x], {nn, 0, 10}]], {mm, 0, 6}]];
WriteString["stdout", StringRiffle[KeyValueMap[#1 <> " :: " <> ToString[InputForm[#2]] &, res], "\n"], "\n"];
