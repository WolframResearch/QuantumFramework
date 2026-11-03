(* Exact WL checks outside QF for section B. No Fock truncation anywhere. *)
res = <||>;
SetAttributes[put, HoldRest];
put[k_String, v_] := (res[k] = TimeConstrained[v, 400, $TimedOut]);
ri = {x, -Infinity, Infinity}; rj = {y, -Infinity, Infinity};

(* B-022: real Gaussian, normalized measure pi^(-n/2) d^n f *)
With[{am = {{2, 1/2}, {1/2, 1}}},
  put["B022", Simplify[Integrate[Exp[-{x, y} . am . {x, y}]/Pi, ri, rj] - Det[am]^(-1/2)]]];

(* B-023: one mode, complex A11, A22, A12 = A21 = 1, complex sources, measure d^2f/pi *)
With[{am = {{1/5 + I/10, 1}, {1, 1/10 - I/5}}, p = {1/3 + I/4, -1/2 + I/5}},
  With[{lhs = Integrate[Exp[-{x + I y, x - I y} . am . {x + I y, x - I y}/2 + p . {x + I y, x - I y}]/Pi, ri, rj],
        rhs = Det[{{am[[2, 1]], am[[2, 2]]}, {am[[1, 1]], am[[1, 2]]}}]^(-1/2) Exp[p . Inverse[am] . p/2]},
    put["B023", {N[lhs, 20], N[rhs, 20], N[lhs - rhs, 50]}]]];

(* B-027: normal-symbol product integral at z0 = 3/10 + i/5, exact *)
With[{z0 = 3/10 + I/5, s1 = zb z^2, s2 = zb^2 z + z},
  With[{s12 = 2 z zb + z^3 zb + 4 z^2 zb^2 + z^3 zb^3},
    put["B027int", Simplify[Integrate[(s1 /. {zb -> Conjugate[z0], z -> x + I y}) (s2 /. {zb -> x - I y, z -> z0}) Exp[-((x - Re[z0])^2 + (y - Im[z0])^2)]/Pi, ri, rj] - (s12 /. {zb -> Conjugate[z0], z -> z0})]]]];

(* exact two-mode Bogoliubov blocks: one-mode squeeze, passive unitary, two-mode squeeze *)
dbl[ph_, ps_] := ArrayFlatten[{{ph, ps}, {Conjugate[ps], Conjugate[ph]}}];
sq1 = dbl[DiagonalMatrix[{5/4, 1}], DiagonalMatrix[{3/4 (3 + 4 I)/5, 0}]];
pv = dbl[{{3/5, -4/5 (3 - 4 I)/5}, {4/5 (3 + 4 I)/5, 3/5}}, ConstantArray[0, {2, 2}]];
tms = dbl[5/3 IdentityMatrix[2], 4/3 I {{0, 1}, {1, 0}}];
sB = sq1 . pv . tms;
{phB, psB} = {sB[[1 ;; 2, 1 ;; 2]], sB[[1 ;; 2, 3 ;; 4]]};
fB = {1/2 + I/3, -1/4 + I};
ct = ConjugateTranspose; tp = Transpose;
put["relations", Simplify[{phB . ct[phB] - psB . ct[psB] - IdentityMatrix[2], phB . tp[psB] - psB . tp[phB]}]];

(* B-048: transformed vacuum with shift, two modes: annihilated by b, and its Gaussian norm fixes |c| as in (4.25) *)
With[{zz = Inverse[phB] . psB, zv = {z1, z2}},
  With[{f0 = Exp[-zv . zz . zv/2 - zv . Inverse[phB] . fB]},
    put["B048pde", Simplify[phB . (D[f0, #] & /@ zv) + psB . zv f0 + fB f0]];
    With[{q = ComplexExpand[Re[-zv . zz . zv/2 - zv . Inverse[phB] . fB] 2 - (x1^2 + y1^2 + x2^2 + y2^2) /. {z1 -> x1 + I y1, z2 -> x2 + I y2}]},
      With[{xv = {x1, y1, x2, y2}},
        With[{qq = -D[q, {xv, 2}]/2, bb = D[q, {xv}] /. Thread[xv -> 0]},
          With[{norm2 = Pi^2/Sqrt[Det[qq]] Exp[bb . Inverse[qq] . bb/4]/Pi^2,
                cmod = Det[phB . ct[phB]]^(-1/4) Exp[Join[fB, Conjugate[fB]] . ArrayFlatten[{{Inverse[tp[phB]] . ct[psB], -IdentityMatrix[2]}, {-IdentityMatrix[2], Inverse[ct[phB]] . tp[psB]}}] . Join[fB, Conjugate[fB]]/4]},
            put["B048norm", {Simplify[Expand[q - (-xv . qq . xv + bb . xv)]], Simplify[Im[cmod]], N[{1/Sqrt[norm2], cmod}, 30], FullSimplify[1/Sqrt[norm2] - cmod]}]]]]]]];

(* B-049: the Gaussian kernel intertwines, b U = U a and b* U = U a*, two modes, exact *)
With[{av = {p1, p2}, sv = {s1, s2}, phi = phB, psi = psB, f = fB},
  With[{k = Exp[(av . Conjugate[psi] . Inverse[phi] . av + av . Inverse[tp[phi]] . sv + sv . Inverse[phi] . av - sv . Inverse[phi] . psi . sv)/2 + av . (Conjugate[f] - Conjugate[psi] . Inverse[phi] . f) - sv . Inverse[phi] . f]},
    put["B049pde", Simplify[{phi . (D[k, #] & /@ sv) + psi . sv k + f k - av k,
      Conjugate[psi] . (D[k, #] & /@ sv) + Conjugate[phi] . sv k + Conjugate[f] k - (D[k, #] & /@ av)}/k]]]];

(* B-071: vacuum column of e^{itH} in the Bargmann realization a -> d/dz, a* -> z *)
vacResidual[am_, cm_, f_, t_, half_] := Module[{n = Length[f], l, m, md, ph, ps, phd, psd, v, g, gd, r, rd, ell, elld, cdot, zv, ex, hpsi},
  l = ArrayFlatten[{{-cm, -am}, {Conjugate[am], Conjugate[cm]}}];
  m = MatrixExp[I t l]; md = I l . m;
  {ph, ps, phd, psd} = {m[[1 ;; n, 1 ;; n]], m[[1 ;; n, n + 1 ;;]], md[[1 ;; n, 1 ;; n]], md[[1 ;; n, n + 1 ;;]]};
  v = Join[-f, Conjugate[f]];
  g = (Inverse[l] . (m - IdentityMatrix[2 n]) . v)[[1 ;; n]]; gd = (I m . v)[[1 ;; n]];
  r = -Inverse[ph] . ps; rd = Inverse[ph] . phd . Inverse[ph] . ps - Inverse[ph] . psd;
  ell = -Inverse[ph] . g; elld = Inverse[ph] . phd . Inverse[ph] . g - Inverse[ph] . gd;
  cdot = -Tr[Inverse[ph . MatrixExp[I t cm]] . (phd . MatrixExp[I t cm] + ph . (I cm) . MatrixExp[I t cm])]/2 + I (half ell . Conjugate[am] . ell + Conjugate[f] . ell);
  zv = Array[zz, n];
  ex = zv . r . zv/2 + ell . zv;
  hpsi = zv . am . zv/2 + Sum[Conjugate[am][[j, k]] (D[ex, zv[[j]], zv[[k]]] + D[ex, zv[[j]]] D[ex, zv[[k]]]), {j, n}, {k, n}]/2 + Sum[zv[[i]] cm[[i, j]] D[ex, zv[[j]]], {i, n}, {j, n}] + f . zv + Conjugate[f] . (D[ex, #] & /@ zv);
  Expand[cdot + zv . rd . zv/2 + elld . zv - I hpsi]
];
With[{res1 = vacResidual[{{3/10 + 2 I/5}}, {{1}}, {1/2 + I/3}, t, 1/2]},
  put["B071oneModeExact", Simplify[res1, Element[t, Reals]]]];
With[{res0 = vacResidual[{{3/10 + 2 I/5}}, {{1}}, {1/2 + I/3}, 63/100, 1]},
  put["B071oneModePrinted", N[res0, 20]]];
With[{am = {{1/2, 1/5 + I/10}, {1/5 + I/10, -1/3 + I/4}}, cm = {{1, 1/4 - I/6}, {1/4 + I/6, 3/2}}, f = {1/3 - I/2, 1/5 + I/7}},
  put["B071twoMode", {Max[Abs[CoefficientList[vacResidual[N[am, 60], N[cm, 60], N[f, 60], N[7/10, 60], 1/2], Array[zz, 2]]]],
    Max[Abs[CoefficientList[vacResidual[N[am, 60], N[cm, 60], N[f, 60], N[7/10, 60], 1], Array[zz, 2]]]]}]];

(* B-085: Riccati fixed point, exact positive matrices *)
With[{be = {{2, 1/2}, {1/2, 1}}, ga = {{1, 1/3}, {1/3, 3/2}}},
  With[{gh = MatrixPower[ga, 1/2], gi = MatrixPower[ga, -1/2]},
    With[{w = gi . MatrixPower[gh . be . gh, 1/2] . gi},
      put["B085", {RootReduce[Flatten[w . ga . w - be]], RootReduce[Flatten[w - tp[w]]], RootReduce[Eigenvalues[w]] // N}]]]];

(* B-107: e^{a* h a} maps exp(z.w) to exp(z.e^h w), so its kernel is exp(a* e^h a); two modes, non-normal complex h *)
With[{h = {{-1, 1/2 + I/3}, {1/5, -2}}},
  With[{gfun = Exp[{z1, z2} . MatrixExp[s h] . {w1, w2}]},
    put["B107flow", Simplify[D[gfun, s] - Sum[{z1, z2}[[i]] h[[i, j]] D[gfun, {z1, z2}[[j]]], {i, 2}, {j, 2}]]]]];
With[{m = {{1/2, 1/8}, {1/8, 1/4}}},
  put["B107trace", Simplify[Integrate[Exp[-({x1, x2} . (IdentityMatrix[2] - m) . {x1, x2} + {y1, y2} . (IdentityMatrix[2] - m) . {y1, y2})]/Pi^2,
     {x1, -Infinity, Infinity}, {x2, -Infinity, Infinity}, {y1, -Infinity, Infinity}, {y2, -Infinity, Infinity}] - 1/Det[IdentityMatrix[2] - m]]]];

WriteString["stdout", StringRiffle[KeyValueMap[#1 <> " -> " <> ToString[InputForm[#2]] &, res], "\n\n"], "\n"];
