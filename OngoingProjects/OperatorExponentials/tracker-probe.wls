PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName, 3], "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
SetAttributes[chk, HoldRest];
chk[id_, e_] := Print[StringPadRight[id, 22],Replace[TimeConstrained[Quiet[e], 120, "TIMEOUT"], {True -> "PASS", False -> "FAIL", x_ :> "OTHER " <> ToString[Short[InputForm[x], 1]]}]];
v = {{bd}, {b}}; v2 = {{bd1, bd2}, {b1, b2}}; v3 = {{bd1, bd2, bd3}, {b1, b2, b3}};
zeroQ[x_, vv_] := Simplify[BosonicNormalOrder[x, vv]] === 0;
expo[x_] := Replace[x, {Exp[y_] :> y, c_ Exp[y_] :> y + Log[c]}];
q = (b + bd)/Sqrt[2]; p = (b - bd)/(I Sqrt[2]);
nn = bd ** b;

chk["P-014", With[{r12 = BosonicExpSimplify[Exp[u bd - uc b] ** Exp[w bd - wc b], v], r21 = BosonicExpSimplify[Exp[w bd - wc b] ** Exp[u bd - uc b], v]},
    Simplify[r12/r21 - E^(u wc - uc w)] === 0]];
chk["P-006", Simplify[BosonicNormalOrder[expo[BosonicExpSimplify[Exp[I PP q] ** Exp[-I QQ p], v]] - (I PP QQ/2 + I (PP q - QQ p)), v]] === 0];
chk["C-004", zeroQ[BosonicExpSimplify[Exp[-I QQ p] ** q ** Exp[I QQ p], v] - (q - QQ), v]];
chk["C-061", With[{x = expo[BosonicExpSimplify[Exp[-(f bd - fc b)] ** Exp[(f + eps g) bd - (fc + eps gc) b], v]]},
    zeroQ[Coefficient[Expand[x], eps] - (g bd - gc b + (fc g - f gc)/2), v]]];
chk["G-001", With[{qq = Sqrt[h/(2 w)] (b + bd), pq = I Sqrt[h w/2] (bd - b)}, Simplify[BosonicNormalOrder[qq ** pq - pq ** qq, v] - I h, h > 0 && w > 0] === 0]];
chk["X-003 null", MatchQ[BosonicExpSimplify[(b1 + z bd2) ** Exp[-z bd1 ** bd2], v2], _NonCommutativeMultiply ? (Last[#] === b1 &)]];
chk["C-046", With[{K0 = (bd1 ** b1 + bd2 ** b2 + 1)/2, Kp = bd1 ** bd2, Km = b1 ** b2},
    zeroQ[K0 ** K0 - (Kp ** Km + Km ** Kp)/2 - ((bd1 ** b1 - bd2 ** b2) ** (bd1 ** b1 - bd2 ** b2) - 1)/4, v2]]];
chk["G-068", Simplify[BosonicExpSimplify[Exp[-(r (b1 ** b2 - bd1 ** bd2))] ** b1 ** Exp[r (b1 ** b2 - bd1 ** bd2)], v2] - (b1 Cosh[r] - bd2 Sinh[r]), r > 0] === 0];
chk["G-071", With[{c = (b1 + b2)/Sqrt[2], cd = (bd1 + bd2)/Sqrt[2], d = (b1 - b2)/Sqrt[2], dd = (bd1 - bd2)/Sqrt[2]},
    zeroQ[expo[BosonicExpSimplify[Exp[(x c ** c - y cd ** cd)/2] ** Exp[-(x d ** d - y dd ** dd)/2], v2]] - (x b1 ** b2 - y bd1 ** bd2), v2]]];
chk["G-095", BosonicMatrixElement[{{1, 1}, {0, 0}}, 1 - I t eta bd1 ** bd2 - (t eta)^2 bd1 ** bd2 ** bd1 ** bd2/2, v2] === -I eta t];
chk["B-053", With[{bs = th (bd1 ** b2 - bd2 ** b1)}, zeroQ[BosonicExpSimplify[Exp[bs] ** (bd1 ** b1 + b2) ** Exp[-bs], v2] -
    (BosonicExpSimplify[Exp[bs] ** (bd1 ** b1) ** Exp[-bs], v2] + BosonicExpSimplify[Exp[bs] ** b2 ** Exp[-bs], v2]), v2]]];
chk["G-037", With[{sp = bd2 ** b3, sm = bd3 ** b2, s3 = bd2 ** b2 - bd3 ** b3, H = gg (b1 ** bd2 ** b3 + bd1 ** bd3 ** b2), Ntot = bd1 ** b1 + bd2 ** b2},
    And[zeroQ[sp ** sm - sm ** sp - s3, v3], zeroQ[s3 ** sp - sp ** s3 - 2 sp, v3], zeroQ[H ** Ntot - Ntot ** H, v3]]]];
chk["G-046", With[{Jp = bd1 ** b2, Jm = bd2 ** b1, J3 = (bd1 ** b1 - bd2 ** b2)/2},
    And[zeroQ[Jp ** Jm - Jm ** Jp - 2 J3, v2], zeroQ[J3 ** Jp - Jp ** J3 - Jp, v2],
        Simplify[BosonicExpSimplify[Exp[I ph J3] ** b1 ** Exp[-I ph J3], v2] - Exp[-I ph/2] b1, ph \[Element] Reals] === 0]]];
chk["C-024", Simplify[Coefficient[Expand[((AA q + BB p) + I (CC q + DD p))/Sqrt[2]], b] - (AA + DD + I (CC - BB))/2] === 0];
chk["C-018", With[{G = x b ** b + y bd ** bd + w bd ** b, W = Sqrt[w^2 - 4 x y]},
    Simplify[BosonicExpSimplify[Exp[G] ** b ** Exp[-G], v] - (b Cosh[W] - (w b + 2 y bd) Sinh[W]/W)] === 0]];
chk["P2-010", With[{r = expo[BosonicExpSimplify[Exp[be bd - bec b] ** Exp[al bd - alc b], v]]},
    Simplify[r - ((al + be) bd - (alc + bec) b + (bec al - be alc)/2)] === 0 || Simplify[r - ((al + be) bd - (alc + bec) b + (be alc - bec al)/2)] === 0]];
chk["P2-012", And[zeroQ[(q ** q + p ** p)/4 - (nn + 1/2)/2, v], zeroQ[(q ** p + p ** q)/4 - I (bd ** bd - b ** b)/4, v]]];
chk["P2-042", BosonicNormalOrder[b ** (ep bd + epc b) - (ep bd + epc b) ** b, v] === ep];
chk["B-075", And[Simplify[BosonicMatrixElement[{al, be}, Exp[c Log[xi] nn], v, "Basis" -> "Coherent"] - Exp[xi^c Conjugate[al] be - (Conjugate[al] al + Conjugate[be] be)/2]] === 0,
    Simplify[BosonicMatrixElement[{n, n}, (c nn - z)^-1, v] - 1/(c n - z)] === 0]];
chk["G-025", Simplify[BosonicMatrixElement[{be, al}, bd, v, "Basis" -> "Coherent"] - Conjugate[be] BosonicMatrixElement[{be, al}, 1, v, "Basis" -> "Coherent"]] === 0];
chk["P2-034", With[{L = Function[X, gam (nu + 1) (bd ** X ** b - (nn ** X + X ** nn)/2) + gam nu (b ** X ** bd - (b ** bd ** X + X ** b ** bd)/2)]},
    And[zeroQ[L[b] + gam b/2, v], zeroQ[L[nn] - (gam nu - gam nn), v]]]];
chk["P2-048", Simplify[Sum[BosonicMatrixElement[{k, k}, Exp[-s (2 nn + 1)], v], {k, 0, Infinity}] - 1/(2 Sinh[s])] === 0];
chk["C-012", Simplify[BosonicMatrixElement[{al, be}, Exp[Log[w] nn], v, "Basis" -> "Coherent"] - Exp[w Conjugate[al] be - (Conjugate[al] al + Conjugate[be] be)/2]] === 0];
chk["B-097", Simplify[BosonicNormalOrder[(E^(-I om t2) b) ** (E^(I om t1) bd), v] - E^(-I om (t2 - t1)) (nn + 1)] === 0];
chk["C-001", BosonicNormalOrder[p ** q - q ** p, v] === -I];
chk["C-057", BosonicNormalOrder[(fc1 b1 + fc2 b2 + fc3 b3) ** (g1 bd1 + g2 bd2 + g3 bd3) - (g1 bd1 + g2 bd2 + g3 bd3) ** (fc1 b1 + fc2 b2 + fc3 b3), v3] === fc1 g1 + fc2 g2 + fc3 g3];
chk["X-004", ! FreeQ[BosonicExpOrder[Exp[x bd1 ** b2 - y bd2 ** b1], v2], Cosh]];
chk["P-024", Simplify[BosonicMatrixElement[{{j + m + 1, j - m - 1}, {j + m, j - m}}, bd1 ** b2, v2] - Sqrt[j - m] Sqrt[j + m + 1]] === 0];
chk["P-026", zeroQ[BosonicExpSimplify[Exp[z bd1 ** b2] ** bd2 ** Exp[-z bd1 ** b2], v2] - (bd2 - z bd1), v2] ||
    zeroQ[BosonicExpSimplify[Exp[z bd1 ** b2] ** bd2 ** Exp[-z bd1 ** b2], v2] - (bd2 + z bd1), v2]];
chk["C-042", zeroQ[BosonicExpSimplify[Exp[z bd1 ** b2] ** (bd2 ** b1) ** Exp[-z bd1 ** b2], v2] - (z bd1 + bd2) ** (b1 - z b2), v2]];
chk["P-030 k=1/4", zeroQ[BosonicExpSimplify[Exp[z bd ** bd/2] ** ((b ** b)/2) ** Exp[-z bd ** bd/2], v] - ((b - z bd) ** (b - z bd))/2, v]];
chk["Idea2 <a|D|b>", With[{ov = Function[{x, y}, Exp[Conjugate[x] y - (x Conjugate[x] + y Conjugate[y])/2]]},
    Simplify[(BosonicMatrixElement[{al, be}, Exp[g bd - Conjugate[g] b], v, "Basis" -> "Coherent"] /. Abs[q_]^2 :> q Conjugate[q]) -
        Exp[(g Conjugate[be] - Conjugate[g] be)/2] ov[al, be + g] /. Conjugate[s_Plus] :> Conjugate /@ s] === 0]];
chk["G-012 E+E", And[Table[BosonicMatrixElement[{k, k}, bd ** ((nn + 1)^-1) ** b, v], {k, 0, 4}] === {0, 1, 1, 1, 1},
    BosonicMatrixElement[{1, 0}, bd ** ((nn + 1)^-1) ** b, v] === 0]];
chk["B-020 |1><1|", And[Table[BosonicMatrixElement[{i, k}, bd ** NormalOrdered[Exp[-nn], v] ** b, v], {i, 0, 2}, {k, 0, 2}] === {{0, 0, 0}, {0, 1, 0}, {0, 0, 0}},
    Simplify[BosonicMatrixElement[{al, be}, bd ** NormalOrdered[Exp[-nn], v] ** b, v, "Basis" -> "Coherent"] -
        Conjugate[al] be Exp[-(Conjugate[al] al + Conjugate[be] be)/2]] === 0]];
chk["shorthand {b}", (BosonicNormalOrder[b ** SuperDagger[b], {b}] /. SuperDagger[b] -> bd) === BosonicNormalOrder[b ** bd, v]];
chk["Idea5 shift", And[BosonicExpSimplify[b ** (nn + 1)^(-1/2), v] === (nn + 2)^(-1/2) ** b,
    BosonicExpSimplify[(nn + 1)^(-1/2) ** bd, v] === bd ** (nn + 2)^(-1/2),
    Simplify[BosonicExpSimplify[GeneralizedPower[NonCommutativeMultiply, b, 2] ** Exp[c nn], v] - Exp[c (nn + 2)] ** GeneralizedPower[NonCommutativeMultiply, b, 2]] === 0]];
chk["Idea5 shift matrix", With[{d = 40, am = N @ SparseArray[Band[{1, 2}] -> Sqrt[Range[39]], {40, 40}], fN = N @ DiagonalMatrix[Table[#[k], {k, 0, 39}]] &},
    Max[Abs[Normal[am . fN[1/Sqrt[# + 1] &] - fN[1/Sqrt[# + 2] &] . am][[1 ;; 10, 1 ;; 10]]]] < 10^-12]];
chk["Idea5 <1|a f(N)|2>", Simplify[BosonicMatrixElement[{1, 2}, BosonicExpSimplify[b ** (nn + 1)^(-1/2), v], v] - Sqrt[2/3]] === 0];
chk["G-012 EE+", And[BosonicExpSimplify[(nn + 1)^(-1/2) ** b ** bd ** (nn + 1)^(-1/2), v] === 1,
    BosonicExpSimplify[(nn + 1)^(-1/2) ** b ** (nn + 1)^(-1/2) ** b ** bd ** (nn + 1)^(-1/2) ** bd ** (nn + 1)^(-1/2), v] === 1]];
chk["G-012 E+E kept", BosonicExpSimplify[bd ** (nn + 1)^(-1/2) ** (nn + 1)^(-1/2) ** b, v] === bd ** (nn + 1)^(-1) ** b];
chk["f-deformed", And[Simplify[BosonicExpSimplify[b ** g[nn] ** g[nn] ** bd, v] - (nn + 1) g[nn + 1]^2] === 0,
    Simplify[BosonicExpSimplify[g[nn] ** bd ** b ** g[nn], v] - nn g[nn]^2] === 0]];
chk["Holstein-Primakoff", Simplify[BosonicExpSimplify[Sqrt[2 j - nn] ** b ** bd ** Sqrt[2 j - nn], v] - (2 j - nn) (nn + 1)] === 0];
chk["b bd untouched", BosonicExpSimplify[b ** bd, v] === b ** bd];
nc[0, _] = 0; nc[_, 0] = 0; nc[x_, y_] := x ** y;
entry[{l_, f_}] := Which[f === 0, 0, l === 1, f /. n -> nn, True, l ** (f /. n -> nn)];
residualZeroQ[h_, u_] := With[{U = Map[entry, u, {2}], dU = Map[entry[{#[[1]], D[#[[2]], t]}] &, u, {2}]},
    AllTrue[Flatten[Map[Table[Simplify[BosonicMatrixElement[{i, k}, BosonicExpSimplify[#, v], v]], {i, 0, 3}, {k, 0, 3}] &,
        lam Table[Sum[nc[h[[i, k]], U[[k, j]]], {k, 2}], {i, 2}, {j, 2}] - I dU, {2}]], # === 0 &]];
chk["G-038 JC  HU = i dU", residualZeroQ[{{0, b}, {bd, 0}},
    {{{1, Cos[lam t Sqrt[n + 1]]}, {b, -I Sin[lam t Sqrt[n]]/Sqrt[n]}}, {{bd, -I Sin[lam t Sqrt[n + 1]]/Sqrt[n + 1]}, {1, Cos[lam t Sqrt[n]]}}}]];
chk["G-041 Raman", residualZeroQ[{{0, nn}, {nn, 0}}, {{{1, Cos[lam t n]}, {1, -I Sin[lam t n]}}, {{1, -I Sin[lam t n]}, {1, Cos[lam t n]}}}]];
toAbs = Abs[z_]^k_ :> (z Conjugate[z])^(k/2);
coh[{x_, y_}, e_] := BosonicMatrixElement[{x, y}, e, v, "Basis" -> "Coherent"] /. toAbs;
kernel[e_, {z_, w_}] := coh[{z, w}, e] Exp[(z Conjugate[z] + w Conjugate[w])/2];
gaussTrace[f_] := With[{r = Association @ CoefficientRules[Expand[PowerExpand[Log[f]]] /. Conjugate[al] -> u /. al -> w, {w, u}]},
    With[{aa = -Lookup[r, Key[{1, 1}], 0], bb = Lookup[r, Key[{1, 0}], 0], cc = Lookup[r, Key[{0, 1}], 0]}, Exp[Lookup[r, Key[{0, 0}], 0]] Exp[bb cc/aa]/aa]];
disp[z_] := Exp[z bd - Conjugate[z] b];
moment[chi_, m_, k_] := (-1)^k D[chi, {g, m}, {gb, k}] /. {g -> 0, gb -> 0};
chk["B-026 kernel on |1>", Simplify[(D[gaussTrace[kernel[disp[h], {z, al}] Exp[s Conjugate[al]] Exp[-al Conjugate[al]]], s] /. s -> 0) -
    (Conjugate[z] - Conjugate[h]) Exp[-h Conjugate[h]/2] Exp[Conjugate[z] h]] === 0];
chk["B-027 kernel product", Simplify[gaussTrace[kernel[disp[h], {z, al}] kernel[disp[k], {al, w0}] Exp[-al Conjugate[al]]] -
    kernel[disp[h] ** disp[k], {z, w0}]] === 0];
chk["P-020 moments |be>", With[{chi = coh[{be, be}, Exp[g bd] ** Exp[-gb b]]},
    Simplify[Table[moment[chi, m, k] - Conjugate[be]^m be^k, {m, 0, 2}, {k, 0, 2}]] === ConstantArray[0, {3, 3}]]];
chk["P-020 moments |2>", With[{chi = BosonicMatrixElement[{2, 2}, Exp[g bd] ** Exp[-gb b], v]},
    Table[moment[chi, m, k] - BosonicMatrixElement[{2, 2}, NonCommutativeMultiply @@ Join[ConstantArray[bd, m], ConstantArray[b, k]], v], {m, 0, 2}, {k, 0, 2}] === ConstantArray[0, {3, 3}]]];
chk["G-024 Fock from diagonal", And @@ (With[{s = Exp[al Conjugate[al]] coh[{al, al}, #[[1]]] /. Conjugate[al] -> u /. al -> w},
    Simplify[Table[(D[s, {u, m}, {w, k}] /. {u -> 0, w -> 0})/Sqrt[m! k!] - (#[[2]][m, k] /. toAbs), {m, 0, 3}, {k, 0, 3}]] === ConstantArray[0, {4, 4}]] & /@
    {{disp[gg], BosonicMatrixElement[{#1, #2}, DisplacementOperator[gg]] &}, {Exp[c nn], BosonicMatrixElement[{#1, #2}, Exp[c nn], v] &}})];
chk["C-031 translation coeff", Simplify[ExpToTrig[gaussTrace[coh[{al, al}, Exp[-I t/2] disp[-a0] ** Exp[-I t nn]]] -
    1/(2 I Sin[t/2]) Exp[I a0 Conjugate[a0] Cot[t/2]/2]], 0 < t < Pi] === 0];
chk["ExpOrder D(-a) e^{-itN}", MatchQ[BosonicExpOrder[disp[-a0] ** Exp[-I t nn], v], _. Power[E, _ ? (FreeQ[#, b] &)] ** _NormalOrdered ** Power[E, _ ? (FreeQ[#, bd] &)]]];
chk["G-040 dispersive", With[{u = Exp[-I chi t (bd1 ** b1 + bd3 ** b3 ** (bd1 ** b1 - bd2 ** b2))]},
    {BosonicMatrixElement[{{1, 0, n}, {1, 0, n}}, u, v3], BosonicMatrixElement[{{0, 1, n}, {0, 1, n}}, u, v3], BosonicMatrixElement[{{1, 0, n}, {0, 1, n}}, u, v3]} ===
        {E^(-I chi (1 + n) t), E^(I chi n t), 0}]];
chk["B-106 free Xi", Simplify[Sum[BosonicMatrixElement[{{k1, k2}, {k1, k2}}, Exp[-be0 ((e1 - mu) bd1 ** b1 + (e2 - mu) bd2 ** b2)], v2], {k1, 0, Infinity}, {k2, 0, Infinity}] -
    1/((1 - Exp[-be0 (e1 - mu)]) (1 - Exp[-be0 (e2 - mu)])), be0 (e1 - mu) > 0 && be0 (e2 - mu) > 0] === 0];
chk["X-007 diagonal M", Simplify[Sum[BosonicMatrixElement[{{k1, k2}, {k1, k2}}, Exp[Log[m1] bd1 ** b1 + Log[m2] bd2 ** b2], v2], {k1, 0, Infinity}, {k2, 0, Infinity}] -
    1/Det[IdentityMatrix[2] - DiagonalMatrix[{m1, m2}]], 0 < m1 < 1 && 0 < m2 < 1] === 0];
dtrace[e_] := FourierTransform[ComplexExpand[coh[{be, be}, e] /. be -> x + I y], {x, y}, {kx, ky}, FourierParameters -> {1, 1}]/Pi /. {kx -> 0, ky -> 0};
recon[e_, A_] := Simplify[Integrate[Simplify[dtrace[e]] (coh[{z, w0}, disp[-(gx + I gy)]] /. {Conjugate[gx] -> gx, Conjugate[gy] -> gy}),
    {gx, -Infinity, Infinity}, {gy, -Infinity, Infinity}]/Pi - coh[{z, w0}, A]] === 0;
chk["P-017 Tr D", Simplify[dtrace[disp[gx + I gy]] - Pi DiracDelta[gx] DiracDelta[gy]] === 0];
chk["P-017 Q, P symbols", With[{sym = Exp[g Conjugate[al] - Conjugate[g] al]},
    And[Simplify[coh[{al, al}, disp[g]] - Exp[-g Conjugate[g]/2] sym] === 0,
        Simplify[(DisplacementOperator[g, Infinity, v, "Ordering" -> "Antinormal"] /. {b -> al, bd -> Conjugate[al]} /. NonCommutativeMultiply -> Times /. toAbs) - Exp[g Conjugate[g]/2] sym] === 0,
        Simplify[gaussTrace[Exp[g Conjugate[g]/2] sym Exp[Conjugate[z] al + Conjugate[al] w0 - al Conjugate[al]]] Exp[-(z Conjugate[z] + w0 Conjugate[w0])/2] - coh[{z, w0}, disp[g]]] === 0]]];
chk["P-017 Tr D(a) D(b)+", With[{t = Simplify[dtrace[disp[gx + I gy] ** disp[-(hx + I hy)]]]},
    Sort[Cases[t, _DiracDelta, Infinity]] === Sort[{DiracDelta[gx - hx], DiracDelta[gy - hy]}] && Simplify[(t /. _DiracDelta -> 1) /. {gx -> hx, gy -> hy}] === Pi]];
chk["P-017 reconstruction", And[recon[disp[gx + I gy] ** b, b], recon[disp[gx + I gy] ** bd, bd], recon[disp[gx + I gy] ** b ** bd, b ** bd]]];
chk["coherent D ** a+", Simplify[coh[{be, be}, disp[g] ** bd] - (Conjugate[be] - Conjugate[g]) coh[{be, be}, disp[g]]] === 0];
realG = {Conjugate[gx] -> gx, Conjugate[gy] -> gy, Conjugate[bx] -> bx, Conjugate[by] -> by};
wave = Exp[ComplexExpand[(gx + I gy) (ax - I ay) - (gx - I gy) (ax + I ay)]];
chk["G-027 P of |be>", With[{chi = coh[{bx + I by, bx + I by}, Exp[(gx + I gy) bd] ** Exp[-(gx - I gy) b]] /. realG},
    Simplify[FourierTransform[ComplexExpand[chi/wave], {gx, gy}, {kx, ky}, FourierParameters -> {1, 1}]/Pi^2 /. {kx -> 0, ky -> 0}] === DiracDelta[ax - bx] DiracDelta[ay - by]]];
chk["G-027 P of |m> weak", Table[Simplify[Integrate[Exp[ax^2 + ay^2]/m! Nest[(D[#, {ax, 2}] + D[#, {ay, 2}])/4 &, DiracDelta[ax] DiracDelta[ay], m] wave,
    {ax, -Infinity, Infinity}, {ay, -Infinity, Infinity}] - BosonicMatrixElement[{m, m}, Exp[(gx + I gy) bd] ** Exp[-(gx - I gy) b], v]], {m, 0, 3}] === {0, 0, 0, 0}];
kernG = coh[{z, w0}, disp[-(gx + I gy)]] /. realG;
ipG = Integrate[#, {gx, -Infinity, Infinity}, {gy, -Infinity, Infinity}]/Pi &;
chk["P-020 rho |be>", Simplify[ipG[(coh[{bx + I by, bx + I by}, disp[gx + I gy]] /. realG) kernG] - coh[{z, bx + I by}, 1] coh[{bx + I by, w0}, 1], Element[{bx, by}, Reals]] === 0];
polyGauss[poly_, gauss_] := Total[(#2 D[gaussTrace[gauss Exp[s1 al + t1 Conjugate[al]]], {s1, #1[[1]]}, {t1, #1[[2]]}] /. {s1 -> 0, t1 -> 0}) & @@@
    CoefficientRules[Expand[poly /. Conjugate[al] -> cu /. al -> cw], {cw, cu}]];
chk["P-020 rho |m>", Table[With[{chi = BosonicMatrixElement[{m, m}, DisplacementOperator[al]] /. toAbs},
    Simplify[polyGauss[Simplify[chi Exp[al Conjugate[al]/2]], Exp[-al Conjugate[al]/2] coh[{z, w0}, disp[-al]]] - coh[{z, 0}, 1] coh[{0, w0}, 1] (Conjugate[z] w0)^m/m!]], {m, 0, 3}] === {0, 0, 0, 0}];
chk["P-030 measure k=3/4", With[{st = BosonicMatrixElement[{2 # + 1, 1}, Exp[ze bd ** bd/2], v] /. ze -> r Exp[I ph] &},
    Table[Simplify[1/(2 Pi) Integrate[r (1 - r^2)^(-1/2) ComplexExpand[st[m] Conjugate[st[n]]], {r, 0, 1}, {ph, 0, 2 Pi}]], {m, 0, 2}, {n, 0, 2}] === IdentityMatrix[3]]];
chk["X-003 TMSV Fock series", With[{s2 = Exp[r E^(I th) bd1 ** bd2 - r E^(-I th) b1 ** b2], as = r > 0 && th \[Element] Reals},
    And[Table[Simplify[BosonicMatrixElement[{{k, k}, {0, 0}}, s2, v2, Assumptions -> as] - E^(I k th) Tanh[r]^k Sech[r], as], {k, 0, 3}] === {0, 0, 0, 0},
        BosonicMatrixElement[{{1, 2}, {0, 0}}, s2, v2, Assumptions -> as] === 0]]];
chk["VEV exponentials", And[Simplify[BosonicVEV[disp[g], v] - Exp[-g Conjugate[g]/2]] === 0,
    BosonicVEV[Exp[r (bd1 ** bd2 - b1 ** b2)], v2, Assumptions -> r > 0] === Sech[r]]];
chk["Idea1 quad3", FreeQ[BosonicExpOrder[Exp[x b ** b + y bd ** bd + w bd ** b], v], Power[E, x b ** b + y bd ** bd + w bd ** b]]];
chk["Idea1 cn + linear", Simplify[TrigToExp[Table[BosonicMatrixElement[{m, n}, BosonicExpOrder[Exp[w nn + mu bd + nu b], v], v] -
    BosonicMatrixElement[{m, n}, E^(mu nu (E^w - 1 - w)/w^2) Exp[mu (E^w - 1)/w bd] ** Exp[w nn] ** Exp[nu (E^w - 1)/w b], v], {m, 0, 2}, {n, 0, 2}]], w > 0] === ConstantArray[0, {3, 3}]];
chk["C-029 vacuum phase", Module[{al, be, ga, t, u, v0, F},
    v0 = BosonicMatrixElement[{0, 0}, BosonicExpOrder[Exp[-I t (al q ** q + be p ** p + ga (q ** p + p ** q))/2], v], v];
    F = MatrixExp[t {{ga, be}, {-al, -ga}}];
    u = (F[[1, 1]] + F[[2, 2]] - I (F[[2, 1]] - F[[1, 2]]))/2;
    And @@ (Simplify[TrigToExp[v0^-2 - u], t > 0 && #] === 0 & /@ Flatten[Outer[And, {al be - ga^2 > 0, al be - ga^2 < 0}, {al + be > 0, al + be < 0}]])]];
chk["root continuous", Module[{t, v0},
    v0 = BosonicMatrixElement[{0, 0}, BosonicExpOrder[Exp[-I t (nn + (1/5 + I/10) bd ** bd + (1/5 - I/10) b ** b)], v], v];
    Simplify[Limit[v0, t -> Pi Sqrt[5]/2, Direction -> "FromBelow"] - Limit[v0, t -> Pi Sqrt[5]/2, Direction -> "FromAbove"]] === 0]];
chk["C-032 driven osc", Module[{f, t, u},
    u = BosonicExpOrder[Exp[-I t (nn + 1/2 + f q)], v];
    Simplify[TrigToExp[Table[BosonicMatrixElement[{m, n}, u, v] - BosonicMatrixElement[{m, n},
        Exp[-I t/2 + I t f^2/2] Exp[-f/Sqrt[2] bd + f/Sqrt[2] b] ** Exp[-I t nn] ** Exp[f/Sqrt[2] bd - f/Sqrt[2] b], v], {m, 0, 2}, {n, 0, 2}]], t \[Element] Reals] === ConstantArray[0, {3, 3}]]];
chk["C-050 Riccati", Module[{l, lb, mu, t, u, wing, om, rr, x},
    u = BosonicExpOrder[Exp[-I t (l bd ** bd/2 + lb b ** b/2 + mu (nn + 1/2)/2)], v];
    wing[o_] := Coefficient[Expand[First @ Cases[u, Exp[z_] /; ! FreeQ[z, o] && FreeQ[z, nn] :> z, {0, Infinity}]], GeneralizedPower[NonCommutativeMultiply, o, 2]];
    om = 2 wing[bd]; rr = wing[b];
    x = Exp[2 Coefficient[Expand[First @ Cases[u, g : Power[_, z_] /; ! FreeQ[z, nn] :> Replace[g, {Power[E, y_] :> y, Power[s_, y_] :> y Log[s]}], {0, Infinity}]], nn]] - 2 rr om;
    Limit[om, t -> 0] === 0 && Simplify[TrigToExp[{I D[om, t] - (lb om^2 + mu om + l), D[-I Log[x], t] - 2 l rr/x + lb om + mu}]] === {0, 0}]];
{chk022, chk054} = Module[{al, be, t, F, fa, fb, fc, fd, al2, al2b, qf, labs = {z1, z2, w1, w2}, x, xi, y, w, same, hk},
    F = MatrixExp[t {{0, be}, {-al, 0}}];
    {fa, fb, fc, fd} = Flatten[F];
    al2[{r_, s_}] := (r + I s)/Sqrt[2]; al2b[{r_, s_}] := (r - I s)/Sqrt[2];
    qf = (BosonicMatrixElement[{aw, aw2}, BosonicExpOrder[Exp[-I t (al q ** q + be p ** p)/2], v], v, "Basis" -> "Coherent"] /. Abs[r_]^2 :> r Conjugate[r]) /.
        {Conjugate[aw] -> al2b[{z1, z2}], Conjugate[aw2] -> al2b[{w1, w2}]} /. {aw -> al2[{z1, z2}], aw2 -> al2[{w1, w2}]};
    same[o_] := And @@ (Simplify[TrigToExp[D[Log[qf], {labs}] - D[Log[o], {labs}]], Element[labs, Reals] && t > 0 && #] === {0, 0, 0, 0} &&
        Simplify[TrigToExp[(qf/o /. Thread[labs -> 0])^2], t > 0 && #] === 1 & /@ {al > 0 && be > 0, al < 0 && be < 0, al > 0 && be < 0 && al + be > 0, al > 0 && be < 0 && al + be < 0}) && Limit[qf/o /. Thread[labs -> 0], t -> 0] === 1;
    y = fa + fd + I (fb - fc); w = (fa - fd + I (fb + fc))/y;
    {x, xi} = {z1, z2} - F . {w1, w2};
    hk = Module[{zq, zp, ex, mm, bb},
        ex = -(F . {zq, zp}) . (F . {zq, zp})/4 - (zq^2 + zp^2)/4 + al2b[{z1, z2}] al2[F . {zq, zp}] + al2b[{zq, zp}] al2[{w1, w2}] - (z1^2 + z2^2 + w1^2 + w2^2)/4;
        mm = -D[ex, {{zq, zp}, 2}]; bb = D[ex, {{zq, zp}}] /. {zq -> 0, zp -> 0};
        ((fa + fd + I (fc - fb))/2)^(1/2)/Sqrt[Det[mm]] Exp[bb . Inverse[mm] . bb/2 + (ex /. {zq -> 0, zp -> 0})]];
    {same[Sqrt[2] y^(-1/2) Exp[-I ((F . {w1, w2})[[1]] z2 - (F . {w1, w2})[[2]] z1)/2] Exp[-(x^2 + xi^2)/4] Exp[-(xi + I x)^2 w/4]], same[hk]}];
chk["C-022 coherent el", chk022];
chk["C-054 Herman-Kluk", chk054];
Print["--- ideas still open (expect FAIL / unevaluated) ---"];
chk["Idea3 symbol", FreeQ[NormalOrdered[Exp[c nn], v]["Symbol", 0], "Symbol"]];
chk["BME a f(N) direct", FreeQ[BosonicMatrixElement[{1, 2}, b ** (nn + 1)^(-1/2), v], BosonicMatrixElement]];
chk["mode between",BosonicExpSimplify[(nn + 1)^(-1/2) ** b ** bd2 ** bd ** (nn + 1)^(-1/2), {{bd, bd2}, {b, b2}}] === bd2];
